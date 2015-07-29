# RV

The BDNYC RV fitting code (Vivienne Baldassare 2011-2012, Adric Riedel 2012-2015)
it is:
* Wavelength-independent (but the wavelength ranges for the object, and the RV standard, must match or at least overlap - it has been used on Optical and NIR data)
* Resolution-independent (it has been successfully used on datasets with R=5000 (SALT Red Side Spectrograph) and R=50,000 (Magellan MIKE))
* Unit independent (seriously... it computes the kilometers-per-second-per-pixel based on the input data, so it can handle data in Angstroms and Microns)
* Fits on both bands and lines
 * And therefore all objects should be paired up with RV standards of a similar spectral type, so the lines and bands in both spectra are close to the same.
It uses cross-correlation to determine the RV shift between the spectrum of an object, and the spectrum of a standard star.
It takes out a third-order polynomial fit and regularizes the data, to reduce the effect of scaling and calibration-related differences between the two spectra. 
It uses 1000 Monte Carlo tries with different additive scaled (to the uncertainty array) noise to remove the effects of noise on the data.

##################################
How to use:
##################################

import find_rv.py, and then run find_rv.radial_velocity with the following arguments:

find_rv.radial_velocity(wavelength_object,flux_object,uncertainty_object, # spectral data on your target (uncertainty is NOT SNR)
              wavelength_standard,flux_standard,uncertainty_standard, # spectral data on a star with known RV, in the same (or at least overlapping) wavelength region
              name_obj,name_standard, # their names, for the plots (see point 2 below for details)
              rv_standard,rv_uncertainty_std, # the radial velocity of the standard, for the plots
              order, # for the plots. Should be the same for both.
              cross_correlation_width, # see point 3 below for what this means
              rv_cut_flag,rv_cut_start,rv_cut_end) # see point 4 below for what this means

Generally, I have used this code with a wrapper that reads in and cleans up input data and runs it through radial_velocity().
Two such wrappers are included in this repository as examples:
NIRSPEC.py, which is designed to read in text files (with names and orders extracted from the filename, so it assumes the files have come from the Keck II NIRSPEC pipeline)
FIRE.py, which uses astropy.io.fits to read in fits files from the Magellan FIREhose v2 pipeline. - Oddities include that the FIRE data is in two .fits files (which this program assumes have well-defined names, and the wavelengths are specified on a log scale

Both NIRSPEC.py and FIRE.py expect the same command-line arguments:
> XXXXX.py standard_spectrum standard_rv standard_rv_uncertainty object_spectrum crop_flag [cross_correlation_width rv_cut_start rv_cut_end]
The items in [] are optional. (see point 5 below for what 'crop flag' means)

######################################
How it works:
######################################
Point 1:

Main operation of find_rv.radial_velocity involves:
1. Constructing a new wavelength array based on the overlap between the two spectra, at 10x input resolution
2. Interpolating both input datasets onto that wavelength array
3. In a loop for N monte carlo tries:
  4. Regularize data (subtract mean, divide by standard deviation)
  5. Remove cubic fit to data, to remove the slope of the data
  6. Cross-correlate the RV standard and object spectrum
  7. Chop out just the region (of generally size 200 units) around the peak of the cross-correlation function.
  8. Fit a gaussian+linear curve to the cross correlation peak. The location of the gaussian peak should be the optimal rv shift (in pixels)
  9. Store that peak location
10. Fit a gaussian to the set of N gaussian peaks. The center of the gaussian should be the correct RV shift; the width of the gaussian the uncertainty based on noise.
11. Convert the pixel shifts to actual kilometer-per-second RVs.
12. Make some plots demonstrating the fit
13. Return the rv and rv uncertainty to the calling thread.

Point 2:
The plots are optional, in fact many of the plots are currently commented out of find_rv.radial_velocity right now. If they are removed (and they are useful, so I don't recommend it), the names of the stars, the order, and the RV standard's RV will become unnecessary

Point 3:
The size of the region chopped out of the cross-correlation function and fit with the gaussian+linear peak can be changed. 
