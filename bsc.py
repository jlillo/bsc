import sys
import os.path
import argparse
import gzip
import shutil

import math
import time
from scipy import integrate
from math import *
import scipy.optimize as op
from scipy.special import erfc
from scipy.optimize import curve_fit
from scipy import stats
import scipy.integrate as integrate
from scipy.optimize import curve_fit

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec # GRIDSPEC !
from matplotlib.pyplot import *

from astropy.io import ascii
from astropy.table import Table, Column
from astroML import stats
from astroML.stats import sigmaG
from astropy import constants as c
from astropy.io.votable import parse_single_table
from astropy.io import fits
from astropy import coordinates
from astroquery.simbad import Simbad
from astrobase.services import trilegal
from astroquery.mast import Catalogs

import jlillo_pypref


def cli():
    """command line inputs

    Get parameters from command line

    Returns
    -------
    Arguments passed by command line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("object", help="Host name resolvable with Simbad")
    parser.add_argument("planet", help="Planet ID: b, c, d, etc.")
    parser.add_argument("depth", help="Transit depth of the planet in ppm")
    parser.add_argument("inst", help="Instrument ID: Instrument_Filter")
    parser.add_argument("file", help="Contrast file with sep,contrast columns")
    parser.add_argument("-C", "--COORD", help="Use coordinates", default=False)
    parser.add_argument("-M", "--MAG", help="Reference magnitude (Tmag/Gmag)")
    args = parser.parse_args()
    return args


def get_coord(tic):
    """
    Get TIC corrdinates

    Returns
    -------
    TIC number
    """
    try:
        catalog_data = Catalogs.query_object(objectname="TIC"+tic, catalog="TIC")
        ra   = catalog_data[0]["ra"]
        dec  = catalog_data[0]["dec"]
        Tmag = catalog_data[0]["Tmag"]
        return ra, dec, Tmag
    except:
    	print("ERROR: No gaia ID found for this TIC")

def simbad_query_radec(obname):
    res = Simbad.query_object(obname)
    SimRA 	= tuple(str(res["RA"][0]).split(' '))
    SimDEC 	= tuple(str(res["DEC"][0]).split(' '))


    if len(SimRA) == 2:
    	SimRA = (SimRA[0], str(int(np.floor(float(SimRA[1])))).zfill(2), str(int(float(SimRA[1]) % np.floor(float(SimRA[1])) )).zfill(2)   )

    if len(SimDEC) == 2:
    	SimDEC = (SimDEC[0], str(int(np.floor(float(SimDEC[1])))).zfill(2), str(int(float(SimDEC[1]) % np.floor(float(SimDEC[1])) )).zfill(2)   )

    c = coordinates.SkyCoord(ra='%sh%sm%ss' % SimRA,
    						 dec='%sd%sm%ss' % SimDEC,
    						 frame='icrs')

    return c.ra.deg, c.dec.deg


def interception(x1,y1,y2):
    """
    Get the interception of a line and a function
    """

    # --- Substract one to the other
    zeros = y1-y2

    # --- Find zeros
    intercept = np.array([]).astype('int')
    for ii in range(len(x1)):
       if zeros[ii-1]*zeros[ii] < 0.: intercept = np.append(intercept,ii)

    if len(intercept) > 1:
        intercept = intercept[1:]
    else:
        intercept = []

    if len(intercept) != 0:
        nintercepts = len(intercept)
        int = np.zeros(nintercepts)
        for jj in range(len(intercept)):
           suby1 = np.array([y1[intercept[jj]-1],y1[intercept[jj]] ])
           subx1 = np.array([x1[intercept[jj]-1],x1[intercept[jj]] ])
           plt.plot(subx1,suby1,'o',c='red')
           srt = np.argsort(suby1)
           int_tmp = np.interp(y2,suby1[srt],subx1[srt])
           int[jj] = int_tmp

    # print(int,intercept)
    # plt.plot(x1,y1)
    # plt.axhline(y2)
    # plt.axvline(int)
    # plt.savefig('tmp2.pdf')
    # plt.close()

    if len(intercept) == 0: int = -1

    return int


def uncompress(filename,outname='trilegal_result_tmp.txt'):
    """
    Uncompress the gz files downloaded from TRILEGAL and rename it
    """
    with gzip.open(filename, 'rb') as f_in:
        with open(outname, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

def trilegal_query(ra,dec,maglim,outname,filtersystem='stroemgren',fieldsize=1.0):
    """
    Get the TRILEGAL catalog
    """

    triquery = trilegal.query_radecl(ra,dec,
                                     magnitude_limit = maglim,
                                     field_deg2      = fieldsize,
                                     filtersystem    = filtersystem,
                                     verbose         = True)
    uncompress(triquery['tablefile'],outname=outname)


def do_hist(catname,filter):
    """
    Perform histogram on TRILEGAL query
    """
    tritab = np.genfromtxt(catname,names=True)
    mags   = tritab[filter]
    hist, loc = np.histogram(mags,bins=2500,range=(4,29))
    return loc,hist


def get_density5(m1,m2,catname,filter,fieldsize=1.):
    """
    2D map of the density of stars in the field depending on mag contrast and set_position
    """
    if m1 >= m2:
        print,'BSC_ERROR: get_density5 --> m1 should be smaller than m2'
        sys.exit()

    loc, hist = do_hist(catname,filter)
    elem = np.where((loc[0:len(hist)] >= m1) & (loc[0:len(hist)] < m2))[0]
    dens = np.sum(hist[elem])/ (fieldsize*3600.**2) # stars per arcsec^2 in the magnitude range
    return dens


def get_filter(inst):
    """
    Get the filter and filtersystem to serarch for in TRILEGAL

    Output:
        filtersystem    Name of the filter system in trilegal_query
        filter          Name of filter
        conv            Conversion factor between TESS and filter. Obtained
                        from linear regression of Tmag and filter magnitudes
                        from the TIC catalog.
    """
    if inst == 'ZORRO_562':
        filtersystem, filter, conv = 'sloan', 'r', 0.954
    elif inst == 'ZORRO_832':
        filtersystem, filter,conv = 'sloan', 'z', 0.920
    elif inst == 'NIRC2_K':
        filtersystem, filter,conv = '2mass', 'Ks', 0.920
    elif inst == 'NIRI_K':
        filtersystem, filter,conv = '2mass', 'Ks', 0.920
    elif inst == 'AstraLux_SDSSz':
        filtersystem, filter,conv = 'sloan', 'z', 0.920
    else:
        print('BSC_ERROR --> Instrument '+inst+' not defined. Please define...')
        sys.exit()
    return filtersystem, filter,conv


def get_sensitivity(t):
    detection = t['detection']
    dist_arr  = t['dist_arr']
    dmag_arr  = t['dmag_arr']
    sens = dist_arr*0.0

    for i,dd in enumerate(dist_arr):
    	sens[i] = np.interp( 0.7, np.cumsum(detection[i,:]), dmag_arr)

    if np.abs(sens[-2]-sens[-1]) > 1:
    	dist_arr = dist_arr[:-1]
    	sens = sens[:-1]

    # for i in range(len(dist_arr)):
    #     print(dist_arr[i],-1.*sens[i])

    return dist_arr,-1.*sens


def BSC(host,planet,ra,dec,depth,inst,kepmag,file_contrast):

    detectability=40.
    completeness=40.

    filtersystem,filter,conv = get_filter(inst)

    delta_max = conv * 2.5*np.log10(depth)  #* 0.947   ; eq 7 in Morton 2011 corrected to be in ∆m_SDSSi band.
    mi        = kepmag # (0.947*kepmag+0.510)
    fieldsize = 1. # square degree

    ycomplete = mi-completeness
    ydetect   = mi-detectability
    magstep   = 0.05
    m_comp    = ycomplete

    if file_contrast.split('.')[-1] == 'npz':
        table_cont = np.load(file_contrast)
        sep,s3 = get_sensitivity(table_cont)
    else:
        table_cont = np.genfromtxt(file_contrast,dtype=None,encoding='utf-8')
        sep,s3 = table_cont[:,0], -1.*table_cont[:,1]

    max_sep = np.max(sep) # arcsec
    maxcontrast = 10 if np.min(s3) > -10 else np.abs(np.floor(np.min(s3))-1)

    # plt.plot(sep,-s3)
    # plt.axhline(delta_max)
    # plt.savefig('tmp2.pdf')
    # plt.close()
    # sys.exit()

    # - Interpolated values
    xbins = 10
    xx    = np.linspace(np.min(sep), max_sep, xbins)
    dalpha =  (max_sep - np.min(sep))/xbins
    Lsens = np.interp(xx, sep, s3)

    # Get TRILEGAL catalog:
    maglim = mi+maxcontrast#-delta_max
    outname = 'ra'+str(round(ra,6))+'_dec'+str(round(dec,6))+\
              '_dm'+str(round(maglim,2))+'_'+filtersystem+\
              '_fov'+str(fieldsize)+'.txt'
    if os.path.isfile(outname) == False:
        print('Running TRILEGAL...')
        trilegal_query(ra,dec,maglim,outname,
                       filtersystem = filtersystem,
                       fieldsize = 1.0)

    # - Get the rho general parameter
    rho = get_density5(mi,mi-delta_max,outname,filter,fieldsize=1.)

    ybins = 10
    ytmp = np.linspace(0,maxcontrast,ybins)
    prob = np.zeros((len(xx),len(ytmp)))
    for ii in range(len(xx)):
      for jj in range(len(ytmp)):
         density = get_density5(kepmag,kepmag+ytmp[jj]+0.5,outname,filter,fieldsize=1.)
         prob[ii,jj]=np.pi*xx[ii]**2*density
    plt.imshow(prob)
    plt.savefig('tmp.png')
    plt.close()

    rho_max2 = interception(sep,s3,delta_max)

    # - Find the maximum separation corresponding to ∆m_max (=delta_max) and its region (subxx,subyy)
    if delta_max < np.min(s3):
        print('\n \t --> delta_max is BELOW the contrast curve...')
        subxx, subyy = sep,s3 #xx, Lsens
        rho_max2 = max_sep
        rho_max  = max_sep
        flag = 0
    else:
        print('\n \t --> delta_max is ABOVE the contrast curve...')
        rho_max2 = interception(sep,s3,delta_max)
        print(rho_max2,Lsens,delta_max)
        if ((len(rho_max2) > 0) & (np.max(Lsens) > delta_max)):
            flag = 0
            rho_max = rho_max2[0]
            subxx = sep[np.where(sep <= rho_max2[0])[0]]
            subyy = s3[np.where(sep <= rho_max2[0])[0]]
        else:
            flag = 1
            rho_max = rho_max2

    # ==== INITIAL probability ---> P_{BS,0}
    PrTotal = np.pi*max_sep**2 *rho
    if flag == 0:
        rho_tmp   = np.zeros(len(xx))
        rho_tmp2  = np.zeros(len(xx))
        xcenters  = np.zeros(len(xx))
        dmag      = np.zeros(len(xx))
        incomp    = np.zeros(len(xx))

        # --- Loop for each bin in ANGULAR SEPARATION
        for ii in range(len(xx)-1):
            # --- Center of the bin
            xcenters[ii] = (xx[ii+1]+xx[ii])/2.
            # --- Magnitudes to be considered
            mags  = np.array([Lsens[ii], m_comp, delta_max])
            # --- Select the UPPER magnitude (above it, all stars are detected)
            maxim = np.where(mags == np.max(mags))[0]
            mymag = mags[maxim]
            # --- Calculate the density of stars from mi to mi-mymag
            if mi < mi-mymag:
                rho_tmp[ii] = get_density5(mi,mi-mymag,outname,filter,fieldsize=1.)
            else:
                rho_tmp[ii] = 0.0
            # --- If there exists an INCOMPLETENESS region then get the incompleteness
            incomp[ii] = 0.0


        PrAstraLux = PrTotal - (2.*np.pi*dalpha * (np.sum(xcenters*rho_tmp) + np.sum(xcenters*incomp) )  )




    else:
        PrAstraLux = 0.0
        rho_max = np.min(sep)

    print('\n\t**Results**')
    print('\t','----------------------------------------')
    print('\t','{0: <30}'.format('Raw Prob: '),'{0: <10}'.format(str(round(PrTotal*100,3))+'%'))
    print('\t','{0: <30}'.format('Prob of undetected BS: '),'{0: <10}'.format(str(round(PrAstraLux*100,3))+'%'))
    print('\t','{0: <30}'.format('Improvement: '),'{0: <10}'.format(str(round((PrTotal-PrAstraLux)/PrTotal * 100,1))+'%'))
    print('\t','----------------------------------------')


    # ==============================================================================
    rho_max2 = np.atleast_1d(rho_max2)
    sep_range = np.max(sep)-np.min(sep)

    fig = plt.figure(figsize=(6.93, 5))
    gs = gridspec.GridSpec(1, 3 , height_ratios=[1], width_ratios=[1,0.1,0.1])
    gs.update(left=0.15, right=0.95, bottom=0.15, top=0.93, wspace=0.2, hspace=0.2)

    ax0 = plt.subplot(gs[0,0])


    splot = plt.imshow(np.transpose(prob*1.e3),extent=(0.0,max_sep,-maxcontrast,0.0),
                       aspect='auto',cmap='RdPu')

    plt.plot(sep,s3,c='white',lw=3,zorder=0)
    plt.plot(sep,s3,c='k',lw=2,zorder=1)

    plt.axhline(delta_max,ls=':',c='k')
    plt.text(0.6*sep_range,delta_max+0.1,
             r'$\Delta m_{\rm max} = $ '+str(round(delta_max,2))+' mag',
             color='k',fontsize=14)
    plt.axvline(np.min(sep),ls=':',c='k')


    # -- Line of interception between Deltam_max and the sensitivity line.
    if delta_max >= np.min(s3):
        plt.axvline(rho_max,ls='--',lw=1)
    # -- If delta_max is too large, simply set it to the minimum plotted forplotting purposes
    if delta_max < -maxcontrast:
        plotdelta_max = -maxcontrast
    else:
        plotdelta_max = delta_max
    # -- If all is correct, plot the shaded region non-detected by
    #    AstraLux: between the sensitivity line and the delta_max.
    if flag == 0:
      if len(rho_max2) == 1:
        plt.fill_between(subxx,subyy*0.+plotdelta_max,subyy,
                         facecolor='lightgreen',alpha=0.3,hatch="/")
    else:
        if (len(rho_max2) % 2) == 0:
            if max_sep-rho_max2[len(rho_max2)] < 0.03:
                rho_max3=np.array([0.0,rho_max])
            else:
                rho_max3=np.array([0.0,rho_max,max_sep])
        elif len(rho_max2) % 2 == 1:
            rho_max3=np.array([0.0,rho_max])

        for kk in range(len(rho_max3)-1):
            tmp_plot = np.where((xx > rho_max3[kk]) & (xx < rho_max3[kk+1]))[0]
            plt.fill_between(xx[tmp_plot],Lsens[tmp_plot]*0.+plotdelta_max,
                             Lsens[tmp_plot],color='white',alpha=0.3)


    plt.text(0.1*sep_range, -0.5,
            r'P(blended source) = '+str(round(PrAstraLux*100,3))+'%',
            color='k', fontsize=14)

    plt.ylim(-9,0)
    plt.xlabel('Angular Separation [arcsec]')
    plt.ylabel('Contrast [mag]')
    plt.title('BSC for '+host+' '+planet)

    # colorbar
    cbax = plt.subplot(gs[0,1])
    plt.colorbar(cax=cbax, orientation='vertical',label=r'P$_{\rm aligned} (\times 10^{-3})$') # , ticks=[0.1,0.5,1]
    #plt.colorbar(cax=cbaxes, orientation='horizontal', ticks=[0.1,0.5,1],label=r'Age (Gyr)')
    # cbax.xaxis.set_label_position('top')
    # cbax.xaxis.set_ticks_position('top')
    # cbax.xaxis.set_major_formatter(formatter)


    path = os.path.dirname(file_contrast)
    if len(path) == 0: path='./'
    plt.savefig(os.path.join(path,host+planet+'_'+inst+'_EBlimits.pdf'))
    plt.close()

    return PrAstraLux

#,/bar,btit='P!daligned!n('+textoidl('\Deltam = \Deltam_i')+','+textoidl('\alpha = \alpha_i')+')')
#,chars = 1.3,position = [0.05,0.05,0.95,0.95]


if __name__ == '__main__':

    print('=====================================')
    print('     Blended Source Confidence       ')
    print('              (BSC)                  ')
    print('          by J.Lillo-Box             ')
    print('=====================================')


    args = cli()

    # ===== Host star
    host = args.object
    if args.COORD:
        coords = args.COORD
        ra, dec = coords.split(',')[0], coords.split(',')[1]
        refmag = float(args.MAG)
    elif host[0:3] == 'TIC':
        ra,dec,_ = get_coord(host[3:])
        refmag   = float(args.MAG)
    else:
        try:
            ra,dec = simbad_query_radec(host)
        except:
            print('--> BSC_ERROR: Object '+host+' could not be resolved.')
            print('               Try a Simbad-readable name')
            sys.exit()

        refmag = float(args.MAG)

    if refmag == None:
        print('--> BSC_ERROR: Refeference magnitude is not defined.')
        print('               Either use a TIC ID or add manually with --MAG option ')
        sys.exit()

    ra  = float(ra)
    dec = float(dec)

    # ===== Planet
    _planets = args.planet
    planets  = _planets.split(',')

    _depths  = args.depth
    depths   = np.array(_depths.split(',')).astype('float')  * 1e-6

    # ===== Instrument & Contrast
    inst    = args.inst   # Instrument ID
    file_contrast = args.file

    # planet  = 'TOI-220b'    # planet name
    # depth   = 905.84e-6     # transit depth
    # ra      = 91.799699     # 83.26925
    # dec     = -61.99721     # -26.72387
    # inst    = 'ZORRO_562'   # Instrument ID
    # kepmag    = 9.6878      # TESS mag
    # file_contrast = 'TOI220_20200109_562.dat'



    count = 1

    for planet,depth in zip(planets,depths):
        print('\nRunning planet '+str(count)+'/'+str(len(planets))+'...')

        print('\nParameters:')
        print('\t','{0: <20}'.format('Planet'.ljust(20, '-')),'{0: <30}'.format(host+planet))
        print('\t','{0: <20}'.format('Depth [ppm]'.ljust(20, '-')),'{0: <30}'.format(depth*1e6))
        print('\t','{0: <20}'.format('RA'.ljust(20, '-')),'{0: <30}'.format(ra))
        print('\t','{0: <20}'.format('DEC'.ljust(20, '-')),'{0: <30}'.format(dec))
        print('\t','{0: <20}'.format('Ref.Mag.'.ljust(20, '-')),'{0: <30}'.format(refmag))
        print('\t','{0: <20}'.format('Inst_ID'.ljust(20, '-')),'{0: <30}'.format(inst))
        print('\t','{0: <20}'.format('FileContrast'.ljust(20, '-')),'{0: <30}'.format(file_contrast))

        BSC(host,planet,ra,dec,depth,inst,refmag,file_contrast)
        count+=1







# trilegal.query_radecl(83.26925, -26.72387,magnitude_limit=14.0,field_deg2=0.1,verbose=True)
