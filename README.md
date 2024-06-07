# bsc - Blended Source Confidence
 Estimate the probability that a given planet host candidate has a 
 blended stellar companion up to a contrast magnitude that can mimic 
 your planetary transit. This code provides paper-ready figures and 
 the probability that your source is isolated from threatening sources.

 ![alt text](https://github.com/jlillo/bsc/blob/f9a6791c17c88f844a1977fc9130205cea68cc4f/logo_bsc.png.001.png)


## Installation & Requirenments
 Clone this folder or download it to your computer. That's it!

The following modules are also required in your python path. 
```
astrobase --> > 0.5.3
astroquery --> > 0.4.6
astroML --> > 1.0.2
```

***bsc*** is written in both Python3.7 but it has also been tested with Python3.11

## Usage

***bsc*** requires a sensitivity curve (aka contrast curve) from a high-spatial 
resolution image that can constrain the probability of your star having a blended source.

The basic input are:

| Input  | Comment |
| ------------- | ------------- |
| `host_ID`  | Host name resolvable with Simbad  |
| `planet_ID`  | Planet ID: b, c, d, etc.  |
| `depth`  | Depth(s) of the planet(s) transit [ppm]  |
| `inst`  | Instrument ID in the format Instrument_Filter. The available options are: AstraLux_SDSSi, NIRI_K, NIRC2_K, ZORRO_832, ZORRO_562  |
| `file`  | Contrast file with sep,contrast columns. This can be either an nzp file with keys "dist_arr" and "dmag_arr" or a 2-column text file with a first column being the separation array and the 2nd column being the contrast.  |

Optional inputs:

| Input  | Comment |
| ------------- | ------------- |
| `--COORD`  | Use coordinates. Default=False  |
| `--MAG`  | Reference magnitude (Tmag/Gmag)  |

To run the code:

```
python bsc.py TOI-2128 b 372 AstraLux_SDSSz TOI-2128_1_SDSSz__000_0100__Sensitivity.npz --MAG 6.68
```

This will return a probability of a blended source capable of mimicking the transit of TOI-2128.01 of 0.034%. 
The code also returns a paper-ready figure showing the contrast curve (black solid line) and the color map
indicating the probability of a chance-aligned source in the separation-contrast parameter space.  The BSC 
probability is claculated as the integral of the space beween the horizontal line (indicating the maximum
contrast a blended source could have to be able to mimic the planetary transity if it is an eclipsing binary)
and the contrast curve. 

![alt text](https://github.com/jlillo/bsc/blob/1da025b550f353ec6611f968894dc86d6ec126d6/TOI-2128b_AstraLux_SDSSz_EBlimits.jpg)

