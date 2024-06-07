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
import glob
from termcolor import colored

import bsc


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

if '__main__':

    print('=====================================')
    print('             BSC wraper       ')
    print('              (BSC)                  ')
    print('          by J.Lillo-Box             ')
    print('=====================================')


    # =====
    # PATHS
    # =====
    root = '/Users/lillo_box/00_projects/11__HighResTESS/data'
    night = 'XXXXXX'
    pathred = '../data/11_REDUCED/XXXXXX/' 
    files = glob.glob(root+'/22_ANALYSIS/XXXXXX/Sensitivity/*.npz')
    files = np.sort(files)

    # =====
    # INFO
    # =====
    astralux_file = 'astralux_tess.csv'
    data = np.genfromtxt(astralux_file,delimiter=',',encoding='utf-8',dtype=None,names=True)
    targets = data['TOI']

    utargets = np.unique(targets)

    pids = np.array(['b','c','d','e','f','g','k'])
    inst = 'AstraLux_SDSSz'

    # ====================
    # LOOP for each image
    # ====================
    for f,_file in enumerate(files):

        print(colored("=======================================================================",'blue'))
        print(colored(_file,'blue'))
        print(colored("=======================================================================",'blue'))

        # ===== Extracting info from filename
        file = os.path.basename(_file)
        TOIname = file.split('_')[0]
        if '-' not in TOIname:
            TOIname = 'TOI-'+TOIname[3:]

        if len(file.split('_')) == 7:
            epoch = '_'+file.split('_')[1]
            date = file.split('_')[4]
        else:
            epoch = ''
            date = file.split('_')[3]

        target = int(TOIname[4:])

        these = np.atleast_1d(np.where(targets.astype(int) == target)[0])
        print(these,target,targets)
        Nplanets = len(these)
        
        planets = pids[:Nplanets]
        depths  = data["Transit_Depth_Value"][these]* 1e-6
        host    = 'TIC'+str(data['TIC'][these[0]])
        refmag  = data['TMag_Value'][these[0]]
        ra,dec,_ = get_coord(host[3:])

        file_contrast = _file
        


        if refmag == None:
            print('--> BSC_ERROR: Refeference magnitude is not defined.')
            print('               Either use a TIC ID or add manually with --MAG option ')
            sys.exit()

        ra  = float(ra)
        dec = float(dec)

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

            PrAstraLux = bsc.BSC(host,planet,ra,dec,depth,inst,refmag,file_contrast)

            path = os.path.dirname(file_contrast)
            print("Saving...")
            np.savez(path+'/'+TOIname+'_'+planet+'_'+inst+'_EBlimits',PrAstraLux=PrAstraLux)
            count+=1



