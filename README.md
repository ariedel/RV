# RV

The BDNYC RV fitting code (Vivienne Baldassare 2011-2012, Adric Riedel 2012-2015)
It uses cross-correlation to determine the RV shift between the spectrum of an object, and the spectrum of a standard star.
It takes out a third-order polynomial fit and regularizes the data, to reduce the effect of scaling and calibration-related differences between the two spectra. 
It uses 1000 Monte Carlo tries with different additive scaled (to the uncertainty array) noise to remove the effects of noise on the data.
It assumes all input data has been corrected for heliocentric/barycentric radial velocity. This must be done externally, either in the data reduction, or as part of a wrapper that prepares the data for this routine.
It is:
* Wavelength-independent (but the wavelength ranges for the object, and the RV standard, must match or at least overlap - it has been used on Optical and NIR data)
* Resolution-independent (it has been successfully used on datasets with R=5000 (SALT Red Side Spectrograph) and R=50,000 (Magellan MIKE))
* Unit independent (it computes the kilometers-per-second-per-pixel based on the input data, so it can handle data in Angstroms and Microns)
* Fits on both bands and lines. In general, it is most sensitive to spectral lines, so make sure all cosmic rays are removed.
 * And therefore all objects should be paired up with RV standards of a similar spectral type, so the lines and bands in both spectra are close to the same.
FINAL NOTE: The radial velocity uncertainties output from this code only relate to the noise in the spectra. There are other sources of error (wavelength calibration errors, for instance) that affect the spectra. It is highly recommended to take a weighted mean of multiple different runs (different standards, different spectra of the target star, different orders if the spectra are from an echelle) to get the TRUE radial velocity uncertainty.

#############
## Update
#############
find_rv_outliers.py will now run until 1000 iterations have been generated within set limits, for use with noisy spectra where some of the crossmatches are clearly unphysical outliers. This is an update to the older version that would simply remove everything outside those limits and calculate based on a smaller number of points.

##################################
How to use:
##################################

import find_rv_outliers.py, and then run find_rv_outliers.radial_velocity with the following arguments:

find_rv_outliers.radial_velocity(wavelength_object,flux_object,uncertainty_object, # spectral data on your target (uncertainty is NOT SNR)
              wavelength_standard,flux_standard,uncertainty_standard, # spectral data on a star with known RV, in the same (or at least overlapping) wavelength region
              name_obj,name_standard, # their names, for the plots (see point 2 below for details)
              rv_standard,rv_uncertainty_std, # the radial velocity of the standard, for the plots
              order, # for the plots. Should be the same for both.
              cross_correlation_width, # set to 200 for 'default', see point 3 below for what this means
              pixel_cut_flag,pixel_cut_start,pixel_cut_end) # set all to zero by default. See point 4 below for what this means

Generally, I have used this code with a wrapper that reads in and cleans up input data and runs it through radial_velocity().
Two such wrappers are included in this repository as examples:
NIRSPEC.py, which is designed to read in text files (with names and orders extracted from the filename, so it assumes the files have come from the Keck II NIRSPEC pipeline)
FIRE.py, which uses astropy.io.fits to read in fits files from the Magellan FIREhose v2 pipeline. - Unique features include that the FIRE data is split across two .fits files (which this program assumes have well-defined names relative to each other), and the wavelengths are specified on a log scale.

Both NIRSPEC.py and FIRE.py expect the same command-line arguments:
> XXXXX.py standard_spectrum standard_rv standard_rv_uncertainty object_spectrum crop_flag [cross_correlation_width rv_cut_start rv_cut_end]
The items in [] are optional. (see point 5 below for what 'crop flag' means)

######################################
How it works:
######################################
Point 1:

Main operation of find_rv_outliers.radial_velocity involves:
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
The plot (step 12) can be considered optional. It shows three graphs: the top shows the spectra, shifted by the fitted RV and overlayed. The middle plot shows the last of the N cross-correlations with gaussian+linear fit and residuals. The bottom graph shows a histogram of the the distribution of pixel shifts and THEIR gaussian fit - this is useful for checking the quality of the final RV result. If they are removed (and they are useful, so I don't recommend it), the names of the stars, the order, and the RV standard's RV will become unnecessary. In particular, "order" is only in there so that the output plots don't overwrite each other. There are several other plots commented out of the data here, anyway.

Point 3:
The size of the region chopped out of the cross-correlation function (Step 7 above) and fit with the gaussian+linear peak can be changed. It must be set if you run find_rv.radial_velocity directly (we recommend 200), and defaults to 200 if not specified when you run either of the wrappers (NIRSPEC.py and FIRE.py). Tweaking this could be useful if the cross-correlation function has a complicated shape and the gaussian+linear fit (step 8) is a visibly poor fit seen on the plot (step 12).

Point 4:
Occasionally, when the SNR of a spectrum is low, the cross-correlation function will occasionally favor a bad match - matching up some similar-but-incorrect region of the spectrum, or a spike in the noise. The most common symptom is that the output RV is in the hundreds or thousands of km/s, the output RV uncertainty is similarly enormous, and the bottom panel of the output plot (from step 12) show two spikes in the pixel shift histogram, rather than a smooth bell curve- one near zero, and one with a pixel shift of +/- thousands of pixels. The assumption is that the spike near 0 pixel shift is the actual RV, and the spike at some enormous pixel shift is some sort of off-lock, where the process was fooled into lining up the spectra in some obvious way (the top plot of the two shifted spectra should also look obviously wrong).
If you set the last three arguments in find_rv.radial_velocity to 1,minimum_pixel_shift,maximum_pixel_shift, you can force the final gaussian fit (step 10) to ignore all values outside those boundaries, isolating just the good RVs. 
USE THIS OPTION SPARINGLY AND WITH CARE. If misused, you could theoretically force the output RV to be whatever you wanted.

find_rv_outliers.py will force the program to keep running until 1000 points are generated within the range [minimum_pixel_shift,maximum_pixel_shift] and should produce more robust results if such a restriction is necessary.

Point 5:
Both NIRSPEC.py and FIRE.py use the 'crop flag' to trim outliers from a spectrum. If you set 'crop flag' to 1, any points more than 3 standard deviations from the mean of the spectrum will be removed, prior to the arrays being sent to find_rv.radial_velocity. 
USE THIS OPTION WITH CARE: This could remove actual useful data from a spectrum that has deep features or lines. It's probably better to write your own more clever routine to remove outliers. Or do it by hand.


