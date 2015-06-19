import math
import numpy
#import BDNYC
import scipy
import random
import find_rv
#from astrolibpy import astrolib
from astrolibpy.astrolib import baryvel
from astropy.io import ascii
from astropy.io import fits
from astropy import coordinates
from astropy import time
from astropy import units
#from astrolib import baryvel
import sys
import matplotlib.pyplot as plt

# Weighted Standard Devation taken from http://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy
def wstddev(x,w):
    average = numpy.average(x,weights=w)
    variance = numpy.dot(w,(x-average)**2)/w.sum()
    return average,math.sqrt(variance)

# AR 2012.1203: Remove telluric features entirely.
# AR 2012.1227: Now a slightly more generic one-d clipping routine.  Assumes a[0] is the index being tested
def clip(a,start,end):
    # Returns the linked arrays 'a' with the entries from 'start' to 'end' clipped out.
    
    #change 3 arrays:
    #    (wave[0], wave[1], wave[2]...; flux[0], flux[1], flux[2]...; er[0], err[1], err[2]...)
    #into 1 array of paired tuples:
    #    ([wave[0],flux[0], err[0]], [wave[1],flux[1],err[1]], [wave[2],flux[2],err[2]], ...).
    b=zip(*a)
    indices = []
    c = []
    # make a new array of just the good elements
    for i in range(len(b)):
        wave = b[i][0]
        #print(wave)
        if not ( wave > start and wave < end):
            c.append(b[i])
            indices.append(i)

            out = zip(*c)
    # return an un-zipped array
    return numpy.asarray(out)

def main(argv=None):
    if argv is None:
        argv = sys.argv
        
    std_path = argv[1]
    
    rv_std = float(argv[2])
    rv_std_err = float(argv[3])
    obj_path = argv[4]
    crop = float(argv[5])
    try:
        xcorr_width = float(argv[6])
    except IndexError:
        xcorr_width = float(200)
    try:
        trimstart = argv[7]
        trimend = argv[8]
        trim = 1
    except IndexError:
        trimstart = 0
        trimend = 0
        trim = 0

    cut = 0
    cutstart = 0
    cutend = 0

    spec_obj= fits.open(obj_path)
    spec_obj_E = fits.open("{0:}_E.fits".format(obj_path.split('_')[0]))
    spec_std= fits.open(std_path)
    spec_std_E = fits.open("{0:}_E.fits".format(std_path.split('_')[0]))

    # put both the name and the UNUM in the title of the file.
    starname_obj = spec_obj[0].header['OBJECT']
    starname_std = spec_std[0].header['OBJECT']
    #order_obj = int(argv[4].split('.')[0][-1])
    #order_std = int(argv[1].split('.')[0][-1])

    #if order_obj != order_std:
    #    raise ValueError

    #spec_obj.close()
    #spec_std.close()
    obj_flux = spec_obj[0].data
    obj_unc = spec_obj_E[0].data
    wstart = spec_obj[0].header['CRVAL1']
    wdel = spec_obj[0].header['CDELT1']
    obj_wave = 10**(wstart + numpy.arange(len(obj_flux))*wdel)

    std_flux = spec_std[0].data
    std_unc = spec_std_E[0].data
    wstart = spec_std[0].header['CRVAL1']
    wdel = spec_std[0].header['CDELT1']
    std_wave = 10**(wstart + numpy.arange(len(std_flux))*wdel)


    if crop == 1:
        # find spikes.  Select all that are NOT.
        mask = numpy.where(abs(obj_flux-numpy.mean(obj_flux)) < 3.0 * numpy.std(obj_flux,ddof=1))
        obj_wave = obj_wave[mask]
        obj_flux = obj_flux[mask]
        obj_unc = obj_unc[mask]
        #print(wave_obj)
    stardata_obj=[]
    stardata_std=[]

    if trim == 1:
        #print(obj_wave)
        keep = numpy.where( (obj_wave > float(trimstart)) & (obj_wave < float(trimend)))
        #print(keep)
        obj_wave = obj_wave[keep]
        obj_flux = obj_flux[keep]
        obj_unc = obj_unc[keep]
        std_wave= std_wave[keep]
        std_flux = std_flux[keep]
        std_unc = std_unc[keep]

    # First: Cross correlate to zero
    #telluric = fits.open('adric.fits')
    #tell = telluric[0].data
    #tel_wave = numpy.asarray(zip(*tell)[0])
    #tel_flux = numpy.asarray(zip(*tell)[1])
    #tel_unc = numpy.zeros_like(tel_flux)*tel_flux

    outputname='std_{0:}_obj_{1:}.txt'.format(starname_std,starname_obj)
    g=open(outputname,'w')

    #orders = [[1,8310,13500],[2,14100,17900],[3,19300,24800]]
    orders = [[1,8310,8800],[2,8800,9200],[3,9200,9600],[4,9600,10000],[5,10000,10400],[6,10400,10800],[7,10800,11200],[8,11200,11600],[9,11600,12000],[10,12400,12800],[11,12800,13150],[12,13150,13500],[13,14100,14450],[14,14450,14800],[15,14800,15200],[16,15200,15600],[17,15600,16000],[18,16000,16400],[19,16400,16800],[20,16800,17200],[21,17200,17550],[22,17550,17900],[23,19300,19650],[24,19650,20000],[25,20000,20400],[26,20400,20800],[27,20800,21200],[28,21200,21600],[29,21600,22000],[30,22000,22400],[31,22400,22800],[32,22800,23200],[33,23200,23600],[34,23600,24000],[35,24000,24400],[36,244000,24800]]

    for i in range(len(orders)):
        obj_order = orders[i][0]
        
        split = numpy.where((std_wave > orders[i][1]) & (std_wave < orders[i][2]))[0]
        tstd_wave = std_wave[split]
        tstd_flux = std_flux[split]
        tstd_unc = std_unc[split]
        split = numpy.where((obj_wave > orders[i][1]) & (obj_wave < orders[i][2]))[0]
        tobj_wave = obj_wave[split]
        tobj_flux = obj_flux[split]
        tobj_unc = obj_unc[split]


        rv_meas,rv_meas_err=find_rv.radial_velocity(numpy.asarray(tobj_wave),numpy.asarray(tobj_flux),numpy.asarray(tobj_unc),numpy.asarray(tstd_wave),numpy.asarray(tstd_flux),numpy.asarray(tstd_unc),starname_obj,starname_std,rv_std,rv_std_err,obj_order,xcorr_width,cut,cutstart,cutend)

        g.write('rv_meas = {0:>+.2f} +/- {1:>+.2f}\n'.format(rv_meas, rv_meas_err))
        # correct radial velocity to standard
        # we can calculate the correct measured RV here.
        rv_obj = rv_std + rv_meas
        # we only want to add together the per-measurement errors and weight them appropriately.
        #  Leave the fixed error on the RV standard out of this. (add in at the end)
        rv_obj_err = numpy.sqrt(rv_meas_err**2 + rv_std_err**2)

        outstring1 = 'std:{1:} obj:{0:}  Order {2:>2d} rv={3:>6.3f} +/- {4:>6.3f},\n'.format(starname_obj,starname_std,obj_order,rv_obj,rv_obj_err)
        print outstring1 
        g.write(outstring1)

    g.close()

if __name__ == "__main__":
    main()
