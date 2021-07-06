import numpy as np
import os
import matplotlib.pyplot as plt
import sys  



import sncosmo # we'll use this to get SN distances
import snana # we'll use this to read supernova light curve files
from astropy.io import fits
import astropy.units as u
from astropy.table import Table
import gzip
import warnings
import get_vpec
class LC:

    def __init__(self):
        self.lightcurve = None
        self.meta = None
        self.zpoff = None
        self.snsed = None
        self.filtertrans = None
        self.primarysed = None
        self.salt3_model = None
        self.offset_lightcurve = None
    
    # there was a method for this above, but this will read in
    # the lightcurve data
    def get_lc_from_snana_file(self,lc_filename, row_entry_label='OBS'):
        # reading in gzipped files is a little more annoying
        if lc_filename.endswith('.gz'):
            with gzip.open(lc_filename,'rt') as fin:
                meta, tables = sncosmo.read_snana_ascii(fin,default_tablename=row_entry_label)
        else:
            meta, tables = sncosmo.read_snana_ascii(lc_filename,default_tablename=row_entry_label)
        self.lightcurve = tables['OBS']
        self.meta = meta
    
    # filter information - mainly the transmission functions
    def get_filter_info(self,kcor_filename):
        warnings.simplefilter('ignore')
        hdu = fits.open(kcor_filename)
        self.zpoff_tab = Table.read(kcor_filename,format='fits')
        self.zpoff = {row['Filter Name'].strip():{col:row[col] for col in [x for x in self.zpoff_tab.colnames if x!='Filter Name']}\
                     for row in self.zpoff_tab}
        self.snsed = hdu[2].data
        temp_filtertrans = hdu[5].data
        self.filtertrans = {temp_filtertrans.names[i]:np.array([x[i] for x in temp_filtertrans]) for i in range(len(temp_filtertrans.names))}
        self.primarysed = hdu[6].data
        hdu.close()
    
    # now, let's initialize the SALT3 model
    def setup_salt3_model(self,modeldir=None):
        # Milky Way dust
        dust = sncosmo.F99Dust()
        if modeldir is None:
            self.salt3_model = sncosmo.Model('salt3',effects=[dust],effect_names=['mw'],effect_frames=['obs'])
        else:
            self.salt3_model = sncosmo.Model(sncosmo.SALT3Source(modeldir=modeldir),
                                            effects=[dust],effect_names=['mw'],effect_frames=['obs'])
        # we need to know the redshift and the MW dst
        if self.meta is not None:
            if 'REDSHIFT_HELIO' in self.meta.keys():
                self.salt3_model.set(z=self.meta['REDSHIFT_HELIO'])
            if 'MWEBV' in self.meta.keys():
                self.salt3_model.set(mwebv=self.meta['MWEBV'])
            
            
    
    def register_sncosmo_bandpass(self,filtwave,filttrans,filtname,unit=u.angstrom):
        band = sncosmo.Bandpass(filtwave,
                                filttrans,
                                wave_unit=u.angstrom,name=filtname)
        sncosmo.register(band, force=True)
    def apply_zpoffs(self):
        if self.lightcurve is None or self.zpoff is None:
            print("Must read lc and kcor before applying offsets.")
            return
        
        if 'SEARCH_PEAKMJD' in self.meta.keys():
            t0 = self.meta['SEARCH_PEAKMJD']
        else:
            t0 = self.meta['PEAKMJD']
        
        self.offset_lightcurve = Table(rows=None,names=['mjd','band','flux','fluxerr','zp','zpsys'],
                     dtype=('f8','S1','f8','f8','f8','U5'),
                     meta={'t0':t0,'z':self.meta['REDSHIFT_HELIO']})

        for i in range(len(self.lightcurve)):
            m = self.lightcurve['MJD'][i]
            flt = self.lightcurve['FLT'][i]
            if flt in self.zpoff.keys():
                sys = self.zpoff[flt]['Primary Name'].strip()
                fullfilt = flt[:]
            elif flt in [k[-1] for k in self.zpoff.keys()]:
                keys_abbrev = np.array([k[-1] for k in self.zpoff.keys()])
                key = np.array([*self.zpoff.keys()])[keys_abbrev == flt][0]
                sys = self.zpoff[key]['Primary Name'].strip()
                fullfilt = key
            else: continue
            #if sys == 'BD17': sys = 'Vega' # SNANA uses mag of Vega in the BD17 system
            flx = self.lightcurve['FLUXCAL'][i]
            flxe = self.lightcurve['FLUXCALERR'][i]
            self.offset_lightcurve.add_row((m,flt,flx*10**(0.4*self.zpoff[fullfilt]['Primary Mag']),
                              flxe*10**(0.4*self.zpoff[fullfilt]['Primary Mag']),
                              27.5,sys))
        
        filters = np.unique(self.offset_lightcurve['band'])
        for band in filters:
            if band not in self.filtertrans and flt in [k[-1] for k in self.filtertrans.keys()]:
                keys_abbrev = np.array([k[-1] for k in self.filtertrans.keys()])
                fullband = np.array([*self.filtertrans.keys()])[keys_abbrev == band][0]
            else: fullband = band[:]
            self.register_sncosmo_bandpass(self.filtertrans['wavelength (A)'],
                                      self.filtertrans[fullband],
                                      band)
def salt3_lc_fit(flc,lc_filename=None,kcor_filename=None,salt3_directory=None,fit_method='fit_lc'):
        if flc.lightcurve is None:
            if lc_filename is None:
                print('Must read in light first or provide filename.')
                return
            flc.get_lc_from_snana_file(lc_filename)
        if flc.zpoff is None:
            if kcor_filename is None:
                print('Must read in kcor info first or provide filename.')
                return
            flc.get_filter_info(kcor_filename)
        flc.apply_zpoffs()
        if flc.offset_lightcurve is None:
            print('Applying offsets failed.')
            return
        if flc.salt3_model is None:
            flc.setup_salt3_model(salt3_directory)

            
        if fit_method == 'fit_lc':
            res, fit = sncosmo.fit_lc(flc.offset_lightcurve,flc.salt3_model,['t0','x0','x1','c'],
                                      modelcov=True)

        sncosmo.plot_lc(flc.offset_lightcurve,fit,errors=res)
        plt.savefig('fit.png',format='png')
        plt.show()
        
        return res