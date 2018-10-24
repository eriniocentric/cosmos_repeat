#!/usr/bin/python

import sys
import os
import os.path
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd

from astropy.io import fits
from astropy.table import Table, Column, MaskedColumn, hstack
from astropy.io import ascii
from astropy.nddata import NDData

rasep=6
decsep=6
wavesep=6

cat=ascii.read('/work/05178/cxliu/maverick/detect/test/catalog/cat.cx3',names=['ID','RA','DEC','wave','sn','chi2','amplitude','sigma','continuum','ncat0','nalines','ngood','fwhm'])

#sort by decreasing S/N to start on highest S/N obs

cat[::-1].sort('sn')

ncat=np.size(cat)
cat['multi_idx']=-1*np.ones(ncat)
cat['nobs']=-1*np.ones(ncat)

output=list()
detectlist=list()
#for loop to group objects according to separations defined above

groupcount=1

for count in range(0,ncat):
    if cat['multi_idx'][count] < 0:
        sel=np.where((abs(cat['RA'][count]-cat['RA'])<(rasep/3600.)) & (abs(cat['DEC'][count]-cat['DEC'])<(decsep/3600.)) & (abs(cat['wave'][count]-cat['wave'])<wavesep)) 
        cat['multi_idx'][sel]=groupcount
        cat['nobs'][count]=np.size(sel)

        
        output.append([groupcount,np.median(cat['RA'][sel]),np.median(cat['DEC'][sel]),np.median(cat['wave'][sel]),np.median(cat['sn'][sel]),np.std(cat['sn'][sel]),np.max(cat['nalines'][sel]),np.max(cat['ncat0'][sel]),np.max(cat['ngood'][sel])])
#        if np.size(sel) > 3: detectlist.append([groupcount,cat['ID'][sel]])
        groupcount+=1

#save output table 


ascii.write(cat,'cosdeep_group_info.dat',overwrite=True)

sel2=np.where( (cat['nobs']>1) )
nuniq=np.size(np.unique(cat['multi_idx'][sel2]))

print "number of unique objects found in 2 or more shots is", nuniq


outputarr = Table(rows=output, names=('multi_idx','ra','dec','wave','sn','snrms','nalines','nobs','ngood'),dtype=('i4','f8','f8','f8','f8','f8','i4','i4','i4') )

#first plot completeness

nbin=20
snarray=np.arange(nbin)*2
compl=np.zeros(nbin)

compl_2frame=np.zeros(nbin)
compl_3frame=np.zeros(nbin)
compl_5frame=np.zeros(nbin)

compl_alines=np.zeros(nbin)

for count in xrange(nbin-1):
    sel=np.where( (outputarr['sn']>snarray[count]) & (outputarr['sn']<snarray[count+1]) & (outputarr['ngood']>0)) 
    if np.size(sel) > 0: 
        compl[count]=float(np.sum(outputarr['nobs'][sel]))/float(np.sum(outputarr['ngood'][sel]))


    sel2=np.where( (outputarr['sn']>snarray[count]) & (outputarr['sn']<snarray[count+1]) & (outputarr['ngood']>0) & (outputarr['nobs']>1))
    if np.size(sel2) > 0:
        compl_2frame[count]=float(np.sum(outputarr['nobs'][sel2]))/float(np.sum(outputarr['ngood'][sel2]))

    sel3=np.where( (outputarr['sn']>snarray[count]) & (outputarr['sn']<snarray[count+1]) & (outputarr['ngood']>0) & (outputarr['nobs']>2))
    if np.size(sel3) > 0:
        compl_3frame[count]=float(np.sum(outputarr['nobs'][sel3]))/float(np.sum(outputarr['ngood'][sel3]))
        
    sel5=np.where( (outputarr['sn']>snarray[count]) & (outputarr['sn']<snarray[count+1]) & (outputarr['ngood']>0) & (outputarr['nobs']>4))
    if np.size(sel5) > 0:
        compl_5frame[count]=float(np.sum(outputarr['nobs'][sel5]))/float(np.sum(outputarr['ngood'][sel5]))

    



plt.figure()
sel=np.where(compl>0.1)
plt.plot(snarray[sel],compl[sel],label='all cat0 objects')
sel2=np.where(compl_2frame>0.1)
plt.plot(snarray[sel2],compl_2frame[sel2],label='Objects detected in 2 or more frames')
sel3=np.where(compl_3frame>0.1)
plt.plot(snarray[sel3],compl_3frame[sel3],label='Objects detected in 3 or more frames')
sel5=np.where(compl_5frame>0.1)
plt.plot(snarray[sel5],compl_5frame[sel5],label='Objects detected in 5 or more frames')
plt.xlabel('Median S/N')
plt.ylabel('N_detected/N_observed')
plt.legend()
plt.savefig("completeness_detects.png")
plt.close()

plt.figure()
plt.hist(cat['sn'],label='All COSDEEP detects catalogs')
sel=np.where( (outputarr['nobs']>1) )
plt.hist(cat['sn'][sel],label='Objects in 2 or more frames')
sel=np.where( (outputarr['nobs']>4) )
plt.hist(cat['sn'][sel],label='Ojbects in 5 or more frames')
plt.yscale('log')
plt.ylabel('N Detections')
plt.xlabel('S/N')
plt.legend()
plt.savefig("hist_SN.png")
plt.close()

