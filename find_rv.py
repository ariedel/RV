#!/usr/bin/env python
# encoding: utf-8
"""
find_rv.py

Created by Vivienne Baldassare on 2012-03-30.
Copyright (c) 2012 __MyCompanyName__. All rights reserved.

Description: 
	This code finds the radial velocity of a target when supplied with data for the target and data for a standard object
	whose radial velocity is known.

Usage:
	Note: Data is not corrected for heliocentric velocities.
	
	Inputs:
		wv_obj, fx_obj, and sig_obj are arrays containing data for the the wavelength, flux, and flux uncertainty of the target.
		wv_std, fx_std, and sig_std are arrays containing data for the the wavelength, flux, and flux uncertainty of the standard.

	Example:
		>>> import find_rv
		>>> find_rv.radial_velocity(wv_obj,fx_obj,sig_obj,wv_std,fx_std,sig_std,rv_std,rv_std_err,obj_name,std_name)

"""

from array import array
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import scipy
import scipy.optimize as op
import scipy.ndimage

def quadratic(x,a,b,c):
    return a + b*x + c*x**2

def cubic(x,a,b,c,d):
	return a + b*x + c*x**2 + d*x**3

def quartic(x,a,b,c,d,e):
	return a + b*x + c*x**2 + d*x**3 + e*x**4


def lsf_rotate(deltav,vsini,epsilon=None,velgrid=None):
    # based on the IDL routine LSF_ROTATE.PRO
    if epsilon == None:
        epsilon = 0.6

    e1 = 2.0*(1.0 - epsilon)
    e2 = np.pi*epsilon/2.0
    e3 = np.pi*(1.0 - epsilon/3.0)

    npts = np.ceil(2*vsini/deltav)
    if npts % 2 == 0:
        npts += 1
    nwid = np.floor(npts/2)
    x = np.arange(npts) - nwid
    x = x*deltav/vsini  
    x1 = np.abs(1.0 - x**2)
    if velgrid != None:
        velgrid = x*vsini
        return (e1*np.sqrt(x1) + e2*x1)/e3,velgrid

    return (e1*np.sqrt(x1) + e2*x1)/e3

def vsini_show(wv_obj,fx_obj,sig_obj,wv_std,fx_std,sig_std,name):
    # Find where standard and object overlap ---------------
	wv_min = max([min(wv_std),min(wv_obj)])
	wv_max = min([max(wv_std),max(wv_obj)])

	n_pix_std = len(wv_std)

	acoef_std = (n_pix_std*10 -1)/(math.log(wv_max) - math.log(wv_min))
	bcoef_std = (n_pix_std*10) - (acoef_std * math.log(wv_max))

	arr = np.arange(n_pix_std*10)+1
	wv_ln_std = np.exp((arr - bcoef_std)/acoef_std)
	
# Interpolate data onto same ln wavelength scale -------------------------------

	fx_interp_std = np.interp(wv_ln_std, wv_std, fx_std) 
	fx_interp_obj = np.interp(wv_ln_std, wv_obj, fx_obj)
	sig_interp_std = np.interp(wv_ln_std, wv_std, sig_std) # AR 2012.1018 Also need to rebin sig
	sig_interp_obj = np.interp(wv_ln_std, wv_obj, sig_obj) # AR 2012.1018 Also need to rebin sig

# Rebin Data ----------------------------

	wv_arr_std=np.asarray(wv_ln_std,dtype=float)
	fx_arr_obj=np.asarray(fx_interp_obj,dtype=float)
	fx_arr_std=np.asarray(fx_interp_std,dtype=float)
	sig_arr_obj=np.asarray(sig_interp_obj,dtype=float)
	sig_arr_std=np.asarray(sig_interp_std,dtype=float)
	
        datalen = len(fx_arr_obj)
        pix_scale = (2.99792458*10**5)/acoef_std

        fig = plt.figure(2)
        ax1 = fig.add_subplot(111)
        fx_arr_obj = fx_arr_obj / np.mean(fx_arr_obj)
        ax1.plot(wv_arr_std,fx_arr_obj,'k')
        ax1.plot(wv_arr_std,fx_arr_std+2,'k')
        vsini = [10,20,30,40,50]
        for i in vsini:
            print('alive')
            kernel = lsf_rotate(pix_scale,i)
            convolve = np.correlate(fx_arr_obj,kernel,mode='same')
            convolve = convolve / np.mean(convolve)
            resid = fx_arr_obj - convolve
            ax1.plot(wv_arr_std,convolve+i/30.0,label='{0:} km/s'.format(i))
            #ax1.text(1.225,0.9+i/40.0,'{0:}: {1:5.2f}'.format(i,np.std(resid,ddof=1)),color='r')
        plt.legend()
        fig.savefig("vsini_{0:}.png".format(name))
        fig.clf()
        plt.close()


def radial_velocity(wv_obj,fx_obj,sig_obj,wv_std,fx_std,sig_std,obj_name,std_name,rv_std,rv_std_err,order,xcorr_width,cut,cutstart,cutend):

# Step 1: Fix the spectra:
# * Select only the region in which they overlap
# * Make a new stretched wavelength array (for sub-pixel precision work)
# * Interpolate the data onto the new wavelength array
# * Remove large scale slopes so we only compare line and band features

# Find where standard and object overlap ---------------
        wv_min = max([min(wv_std),min(wv_obj)])
        wv_max = min([max(wv_std),max(wv_obj)])

        n_pix_std = len(wv_std)

# Creates ln standard wavelength array ---------------------------------
        # AR 2013.0423 The wavelength array only covers the overlap region.  Also, I'm folding the rebinning by 10 into this statement.
	acoef_std = (n_pix_std*10 -1)/(math.log(wv_max) - math.log(wv_min))
	bcoef_std = (n_pix_std*10) - (acoef_std * math.log(wv_max))

	arr = np.arange(n_pix_std*10)+1
	wv_ln_std = np.exp((arr - bcoef_std)/acoef_std)
	
# AR 2012.1018: Find the conversion between pixels and velocity.  This will vary from instrument
#  to instrument and spectral order to spectral order, so we should preferentially calculate this
#  based on the actual input spectrum.
# AR 2013.0422: Change the calculation to happen AFTER the corrected wavelength scale has been made
        # Find the average pixel/spectrum offset
	#  Note: even though it's called micron_per_pix, it will still work if the wavelengths are
	#  angstroms instead (it really converts <wavelength unit> to km/s)

# Interpolate data onto same ln wavelength scale -------------------------------

	fx_interp_std = np.interp(wv_ln_std, wv_std, fx_std) 
	fx_interp_obj = np.interp(wv_ln_std, wv_obj, fx_obj)
	sig_interp_std = np.interp(wv_ln_std, wv_std, sig_std) # AR 2012.1018 Also need to rebin sig
	sig_interp_obj = np.interp(wv_ln_std, wv_obj, sig_obj) # AR 2012.1018 Also need to rebin sig

# Rebin Data ----------------------------

	wv_arr_std=np.asarray(wv_ln_std,dtype=float)
	fx_arr_obj=np.asarray(fx_interp_obj,dtype=float)
	fx_arr_std=np.asarray(fx_interp_std,dtype=float)
	sig_arr_obj=np.asarray(sig_interp_obj,dtype=float)
	sig_arr_std=np.asarray(sig_interp_std,dtype=float)
	
        datalen = len(fx_arr_obj)

# Step 2: Measure vsini:
#  Note that as of 2015.0605, this doesn't actually work.

# AR 2014.0922: For vsini: 
# In a loop: 
#   Take the standard spectrum
#   broaden it to width X
#   autocorrelate,
#   measure width of gaussian Y (this is supposed to give you a means of translating between width-of-cross-correlation and vsini)
# Fit function solving Y for X.
# For each cross correlation of object and standard:
# Determine vsini

        pix_scale = (2.99792458*10**5)/acoef_std

        vsinirange = [1,2,5,10,20,30,40,50,60,80,100,100]
        widthrange = []
        for v in vsinirange:
            # Make convolution kernel for v km/s
            kernel = lsf_rotate(pix_scale,v)
            # Broaden the standard spectrum
            fx_obj_wide = np.correlate(fx_arr_obj, kernel, mode='same')
            # Rectify the spectrum
            fx_obj_orig = (fx_arr_obj - np.mean(fx_arr_obj))/np.std(fx_arr_obj,ddof=1)
            fx_obj_wide = (fx_obj_wide - np.mean(fx_obj_wide))/np.std(fx_obj_wide,ddof=1)

            # Remove a cubic (flatten the spectrum)
            coeff,pcov = op.curve_fit(cubic,wv_arr_std,fx_obj_wide)
            fx_obj_wide = fx_obj_wide - (coeff[0] + coeff[1]*wv_arr_std + coeff[2]*wv_arr_std**2 + coeff[3]*wv_arr_std**3)
            coeff,pcov = op.curve_fit(cubic,wv_arr_std,fx_obj_orig)
            fx_obj_orig = fx_obj_orig - (coeff[0] + coeff[1]*wv_arr_std + coeff[2]*wv_arr_std**2 + coeff[3]*wv_arr_std**3)

            # Cross-correlate the spectrum with its broadened self
            ycorr = np.correlate(fx_obj_orig, fx_obj_wide, mode='full')
            # Now determine where the peak is (should be near 0)
            length = len(ycorr)
            xcorr = np.arange(length) - length//2
            xmid = np.argmax(ycorr)
            ymax = np.max(ycorr)
            # Chop out just the portion of the array near the peak
            xcorr_min=xmid-xcorr_width
            xcorr_max=xmid+xcorr_width
            ycorr1=ycorr[xcorr_min:xcorr_max]	#isolate section of array with gaussian
            xcorr1=xcorr[xcorr_min:xcorr_max]       #isolate the same section of the pixel range

            # set up initial values for gaussian fitting via chi2
            sig = 10
            sky = np.min(ycorr1)/1.2
            #                print ycorr1[-1],ycorr1[0],xcorr1[-1],xcorr1[0]
            sky2 = (ycorr1[-1]-ycorr1[0])/(xcorr1[-1]-xcorr1[0])
            lnamp = np.log(ymax/1.2-sky)	# guess some values
            mean = xcorr[xmid]
            
            amp = np.exp(lnamp)
            sig2 = sig**2
            # suggestion from D. Hogg 12/15/12: Add extra linear feature to fit.
            # suggestion from D. Hogg 12/15/12: operate on ln(amp) so that the amplitude CANNOT be negative.
            def chi2(p):	#define gaussian function for fitting
                sig2=p[2] ** 2
                m = (np.exp(p[0]) * np.exp(-0.5 * (xcorr1 - p[1]) ** 2 / sig2)) + p[3] + p[4]*xcorr1
                return (ycorr1 - m)

            # Fit the gaussian.
            popt, ier = op.leastsq(chi2, [lnamp, mean, sig, sky, sky2])
            lnamp, mean, sig, sky, sky2 = popt
            
            amp = np.exp(lnamp)
            # record the width
            widthrange.append(sig)

        # Plot all the widths to get a width-vsini curve
        vsinicoeff,popt = op.curve_fit(quartic,np.asarray(widthrange),np.asarray(vsinirange))

        relationx = np.arange(50,200,1)
        relationy = vsinicoeff[0]+vsinicoeff[1]*relationx+vsinicoeff[2]*relationx**2+vsinicoeff[3]*relationx**3+vsinicoeff[4]*relationx**4
	figv = plt.figure(1)
        axv = figv.add_subplot(211)
        axv.scatter(widthrange,vsinirange)
        axv.plot(relationx,relationy)
        #ax.text(70,100,"{0:} {1:} {2:} {3:} {4:}".format(vsinicoeff))

# 3. Cross-correlate the data, using 5000 trials:
#  * Generate two random gaussian noises scaled to the uncertainty on the fluxes
#  * Apply those gaussian noises to the standard and target stars
#  * Cross-correlate the standard and target stars
#  * Find and then cut out just the part of the cross-correlation curve near the maximum
#  * Set up gaussian
#  * Fit gaussian to that center part
#  * Save fitted parameters (pixel shift aka mean of gaussian, width aka stddev of gaussian)
#  * Repeat 5000 times

# Cross correlation loop -------------------------------- 
	pix_shift=np.zeros(5000)	#initialize array for pixel shift values
	pix_width=np.zeros(5000)	#initialize array for pixel width values
	l = 0


        # using the xrange generator rather than making a full list saves memory
        for l in xrange(5000):
            # prepare the randomized data
            # GETTING ARRAYS READY FOR CROSS CORRELATION

    
            # Randomize noise:
            # create gaussian distribution of random numbers b/t 1 and -1, multiply err by numbers, add numbers to flux
            # I have drastically simplified the arrays here AR 2013.0319
            # AR 2013.0318: There was a problem, previously: noise was a fixed value, not linked to the known error values

            # AR 2013.0321: Speed fix.  Rather than step through the array and generate one
            #  normally-distributed error value scaled to the SNR at that point, I will generate an
            #  array of normally-distributed error values scaled to 1, and then multiply by the SNR:
            #  One array generation, one array multiplication.

            rand_dist = np.random.normal(loc=0.0,scale=1.0,size=datalen)
            rand_dist2 = np.random.normal(loc=0.0,scale=1.0,size=datalen)
        
            fx_temp_obj = np.asarray(fx_arr_obj + rand_dist * sig_arr_obj)
            fx_temp_std = np.asarray(fx_arr_std + rand_dist2 * sig_arr_std)
            mean_obj=np.mean(fx_temp_obj)
            mean_std=np.mean(fx_temp_std)
            stddev_obj=np.std(fx_temp_obj,ddof=1)
            stddev_std=np.std(fx_temp_std,ddof=1)

            # Regularize data (subtract mean, divide by std dev) (Should definitely be done AFTER noise was added)
            fx_reg_temp_obj = fx_temp_obj-mean_obj
            fx_reg_temp_obj = fx_reg_temp_obj/stddev_obj
            fx_reg_temp_std = fx_temp_std-mean_std
            fx_reg_temp_std = fx_reg_temp_std/stddev_std

    # curve fit - remove a cubic AR 2012.1113
            coeff,pcov = op.curve_fit(cubic,wv_arr_std,fx_reg_temp_obj)
            fx_reg_temp_obj = fx_reg_temp_obj - (coeff[0] + coeff[1]*wv_arr_std + coeff[2]*wv_arr_std**2 + coeff[3]*wv_arr_std**3)
            coeff,pcov = op.curve_fit(cubic,wv_arr_std,fx_reg_temp_std)
            fx_reg_temp_std = fx_reg_temp_std - (coeff[0] + coeff[1]*wv_arr_std + coeff[2]*wv_arr_std**2 + coeff[3]*wv_arr_std**3)
            
            # CROSS CORRELATION 
            
            # compute the cross-correlation between the two spectra

            ycorr = np.correlate(fx_reg_temp_obj, fx_reg_temp_std, mode='full')
            # time required: 0.045 seconds average

            #http://stackoverflow.com/questions/12323959/fast-cross-correlation-method-in-python
            #conv1 = np.zeros(datalen * 2)
            #conv1[datalen/2:datalen/2+datalen] = fx_reg_temp_obj
            #conv2 = fx_reg_temp_std[::-1]
            #ycorr = signal.fftconvolve(conv1,conv2, mode='valid')
            # time required: 0.006 seconds average, but it segfaults by the third try.

            ## slight smoothing AR 2013.0315
            #ycorr = scipy.ndimage.filters.gaussian_filter1d(ycorr,11)
                
            # create the x offset axis (same length as ycorr, with 0 in the MIDDLE)
            length = len(ycorr)
            xcorr = np.arange(length) - length//2
            # AR 2012.1126 Select a tiny piece around the maximum to fit with a gaussian.
            xmid = np.argmax(ycorr)
            ymax = np.max(ycorr)
            # now take just the portion of the array that matters
            xcorr_min=xmid-xcorr_width
            xcorr_max=xmid+xcorr_width
            ycorr1=ycorr[xcorr_min:xcorr_max]	#isolate section of array with gaussian
            xcorr1=xcorr[xcorr_min:xcorr_max]       #isolate the same section of the pixel range
            ycorr2=ycorr[xcorr_min-50:xcorr_max+50]
            xcorr2=xcorr[xcorr_min-50:xcorr_max+50]

            # suggestion from D. Hogg 12/15/12: Add extra linear feature to fit.
            # suggestion from D. Hogg 12/15/12: operate on ln(amp) so that the amplitude CANNOT be negative.
            def chi2(p):	#define gaussian function for fitting
                sig2=p[2] ** 2
                m = (np.exp(p[0]) * np.exp(-0.5 * (xcorr1 - p[1]) ** 2 / sig2)) + p[3] + p[4]*xcorr1
                return (ycorr1 - m)
            
            # set up initial values for chi2
            sig = 10
            sky = np.min(ycorr1)/1.2
            #                print ycorr1[-1],ycorr1[0],xcorr1[-1],xcorr1[0]
            sky2 = (ycorr1[-1]-ycorr1[0])/(xcorr1[-1]-xcorr1[0])
            lnamp = np.log(ymax/1.2-sky)	# guess some values
            mean = xcorr[xmid]
            
            amp = np.exp(lnamp)
            sig2 = sig**2

            popt, ier = op.leastsq(chi2, [lnamp, mean, sig, sky, sky2])
            lnamp, mean, sig, sky, sky2 = popt
            
            amp = np.exp(lnamp)
            
            print_num=l%500		#prints data every 200 fits
            if print_num == 0:
                #fig = plt.figure(l)
                #ax = fig.add_subplot(111)
                #my_gauss = (amp * (np.exp(-0.5 * ((xcorr1 - mean) ** 2) / sig**2))) + sky + sky2 * xcorr1
                #ax.plot(xcorr1,my_gauss,'r--')
                #ax.plot(xcorr2,ycorr2,'#000000')
                #ax.plot(xcorr1,ycorr1-my_gauss,'#00CC00')
                ##if abs(mean - xcorr[xmid]) > 5:
                ##    print "Mean is off",mean,xcorr[xmid]
                #figname='rv_{0:}_{1:}_{2:}_{3:}.png'.format(std_name,obj_name,order,l)
                #ax.set_xlim(xcorr[xcorr_min-50],xcorr[xcorr_max+50])
                #fig.savefig(figname)
                #fig.clf()
                #plt.close()
                print "amp={0: 12.4f}  mu={1: 10.4f}  sig={2: 9.4f}  sky={3: 11.4f}  sky2={4: 8.4f}".format(amp,mean,sig,sky,sky2)
            
            # if ier < 5:
            pix_shift[l] = mean
            # I'm calculating the vsini now because I need errors, and the vsini calculation is not linear.
            pix_width[l] = vsinicoeff[0] + vsinicoeff[1] * sig + vsinicoeff[2] * sig**2 + vsinicoeff[3] * sig**3 + vsinicoeff[4] * sig**4

# End cross correlation loop --------------------------------- 

        fig = plt.figure(2)
        ax = fig.add_subplot(111)
        ax.plot(xcorr,ycorr,'k')
        figname='rv_{0:}_{1:}_{2:}_xcorr.png'.format(std_name,obj_name,order)
        fig.savefig(figname)
        fig.clf()
        plt.close()
	#print len(ycorr)	
	#print pix_shift

	#pix_shift=np.array(pix_shift)
	#(mu,sigma)=norm.fit(pix_shift)	# get mean and std dev of array of pixel shift values
        # Running time = 0.286 seconds

        pix_shift2 = np.copy(pix_shift)

        # cut the outliers from the pixel shift
        if cut == 1:
            print cutstart,cutend
            pix_shift2 = pix_shift2[np.where(pix_shift2 > np.float(cutstart))]
            pix_shift2 = pix_shift2[np.where(pix_shift2 < np.float(cutend))]

        #print len(pix_shift),len(pix_shift2),pix_shift

        mu = np.mean(pix_shift2)
        sigma = np.std(pix_shift2,ddof=1)

        vsini = np.mean(pix_width)
        vsini_err = np.std(pix_width,ddof=1)

        axh = figv.add_subplot(212)
	n, bins, patches=axh.hist(pix_width,bins=30,normed=1.0,facecolor='green',align='mid')
        figv.savefig('vsiniplot.png')
        plt.clf()
        plt.close()
        
        # Running time =

# Transform pixel shift to shift in radial velocity -------------------------------- 

        # AR 2013.0423: The actually appropriate method requires a speed-of-light correction
	#rv_meas=pix_to_kms*mu # AR 2012.1018: no longer NIRSPEC specific
	#rv_meas_err=pix_to_kms*sigma # AR 2012.1018: no longer NIRSPEC specific
	#print "vshift=",rv_meas
        rv_meas = (2.99792458*10**5 * mu)/acoef_std
        rv_meas_err = (2.99792458*10**5 * sigma)/acoef_std

# Apply shift to arrays --------------------------------
# AR 2013.0423: This is a terrible way to do it.  Instead, we will shift THE WAVELENGTH ARRAY in the same fashion as a heliocentric radial velocity correction.
#	
#	fx_rebin_list_obj=fx_reg_temp_obj.tolist()
#	fx_rebin_list_std=fx_reg_temp_std.tolist()
#	#	print 'mu=',mu
#	if mu >= 0:
#		val= abs(mu)	# so we can shift properly
#		i=0
#		while i < val:
#			del fx_rebin_list_obj[0]
#			fx_rebin_list_obj.append(1)
#			i=i+1
#		print 'mu is negative'		
#	elif mu < 0:
#		val=mu
#		i=0
#		while i < val:
#			del fx_rebin_list_std[0]	
#			fx_rebin_list_std.append(1)
#			i=i+1
#			#	print 'mu=',mu		
        wv_rvcorr_obj=wv_arr_std * (1 - rv_meas/(2.99792458*10**5))

## Create plots --------------------------------- 
# Plot object and standard so you can clearly see that shift exists --------------------------------
	fig = plt.figure(1)
	
        # AR 2013.0703 Regularize the spectra for display purposes in the final graph
        # I'm using the mean and stddev of the last random-added attempt so it won't be perfect...
        fx_reg_obj = fx_arr_obj-mean_obj
        fx_reg_obj = fx_reg_obj/stddev_obj
        fx_reg_std = fx_arr_std-mean_std
        fx_reg_std = fx_arr_std/stddev_std

	#Plots target and standard with shift applied
	ax1 = fig.add_subplot(311)
	ax1.plot(wv_rvcorr_obj, fx_reg_obj, 'red')
	ax1.plot(wv_arr_std, fx_reg_std, 'blue')
	ax1.set_xlabel('wavelength (microns)')
	ax1.set_ylabel('normalized flux')
	target = 'Target: %s' %(obj_name)
	standard = 'Standard: %s' %(std_name)
	ax1.annotate(target,xy=(.7,.9),xycoords='axes fraction',xytext=(.6,.9),textcoords='axes fraction',color='red') 
	ax1.annotate(standard,xy=(.7,.8),xycoords='axes fraction',xytext=(.6,.8),textcoords='axes fraction',color='blue') 
	
        sig2=sig ** 2	
        my_gauss = (amp * (np.exp(-0.5 * ((xcorr1 - mu) ** 2) / sig2))) + sky + sky2 * xcorr1

	#Plots example of gaussian fit to cross correlation function
	ax2 = fig.add_subplot(312)
	ax2.plot(xcorr1,  ycorr1, 'k.')
	ax2.plot(xcorr1, my_gauss, 'r--', linewidth=2)
        ax2.plot(xcorr1,ycorr1-my_gauss,'#00CC00')
	ax2.set_xlabel('example of fit to cross correlation function')
        ax2.set_xlim(xcorr[xcorr_min-50],xcorr[xcorr_max+50])	
	#print pix_shift


## Plot histogram of pixel shift values -------------------------------- 
	ax3 = fig.add_subplot(313)
	n, bins, patches=plt.hist(pix_shift,bins=30,normed=1.0,facecolor='green',align='mid') 
	#Plot best fit gaussian over histogram
	y=mlab.normpdf(bins,mu,sigma)
	ax3.plot(bins,y,'r--',linewidth=2)
	ax3.set_xlabel('radial velocity of target (pixels)')
	ax3.set_ylabel('frequency (normalized)')
	rad='RV = %.3f +/- %.3f' %(rv_meas,rv_meas_err)
        corr = 'RV (corr) = %.3f +/- %.3f' %(rv_std + rv_meas, (rv_std_err**2 + rv_meas_err**2)**(0.5))
        vsinistr = 'VsinI = %.3f +/- %.3f' % (vsini,vsini_err)
	ax3.annotate(rad,xy=(.66,.9),xycoords='axes fraction',xytext=(.66,.9),textcoords='axes fraction',color='black')
	ax3.annotate(corr,xy=(.6,.8),xycoords='axes fraction',xytext=(.60,.8),textcoords='axes fraction',color='black')
        ax3.annotate(vsinistr,xy=(.6,.6),xycoords='axes fraction',xytext=(.60,.6),textcoords='axes fraction',color='black')
	ax3.annotate('{0:+5.2f} {1: 5.2f}'.format(mu,sigma),xy=(.05,.9),xycoords='axes fraction',xytext=(.05,.9),textcoords='axes fraction',color='black')
	ax3.annotate('{0:5.3f} km/s/pix'.format((2.99792458*10**5)/acoef_std),xy=(.05,.8),xycoords='axes fraction',xytext=(.05,.8),textcoords='axes fraction',color='black')
	fig.subplots_adjust(hspace=.3)

	figname='rv_%s_%s_%d.png' %(std_name,obj_name,order)
	fig.savefig(figname)
	fig.clf()
        plt.close()
	
	#plt.figure(l+1)
	#plt.hist(pix_shift)
	
#END RADIAL VELOCITY FUNCTION -----------------------------------------
	return rv_meas,rv_meas_err,vsini,vsini_err
	
