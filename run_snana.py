#!/usr/bin/env python
# D. Jones - example for running snlc_fit.exe for a given VERSION
# and plotting the output.  This is from the SuperNova ANAlysis
# code (SNANA), which is an annoying but extremely powerful tool
# for SN Ia distance measurements
#
# steps are
# 1. take input version and kcor file
# 2. write NML input file for that version and kcor file
# 3. write peculiar velocities and add to NML file - ignore for now
# 4. run the nml file
# 5. plot the output as a PDF file
# 6. run SALT2mu to go from light curve parameters to distances
# 7. parse the outputs and create a list of distances
#
# annoying SNANA things
# 1. each version refers to a directory with a .LIST file that contains all
#    the SN light curve files and a .README that contains some documentation
# 2. SNANA will run each directory individually
# 3. if you get error "Could not find IDSURVEY for SURVEY=CfA3K", open $SNDATA_ROOT/SURVEY.DEF and add your SURVEY with a random integer after it
# 4. zHEL is heliocentric redshift, zCMB is CMB-frame redshift, zHD is Hubble diagram redshift which is just zCMB but adding any peculiar velocity
#    corrections that you might have provided.
# 4. more annoying clarifications to come soon
# 
# example run
# python run_snana.py /Users/David/Dropbox/research/hubbleconst/CfA3_NATSYS_KEPLERCAM kcor_Hicken2009.fits hjkl -d /Users/David/Dropbox/research/hubbleconst

import numpy as np
import argparse
import os
import astropy.table as at

_nml_template = """

&SNLCINP
  PRIVATE_DATA_PATH = '<data_dir>'
  USE_MINOS = False
  VERSION_PHOTOMETRY = '<version>'
  SNTABLE_LIST = 'FITRES LCPLOT(text:key)'
  TEXTFILE_PREFIX = '<output_dir>/<version>'
  KCOR_FILE = '<kcorfile>'
  NFIT_ITERATION = 3
  INTERP_OPT = 1
  ZTOL_HELIO2CMB = 0.1
  OPT_SETPKMJD = 5
  OPT_MWEBV = 3
  OPT_MWCOLORLAW = 99
  ABORT_ON_NOEPOCHS = False
  ABORT_ON_TRESTCUT = False
  ABORT_ON_DUPLCID = False
  CUTWIN_NEPOCH = 5
  CUTWIN_TREST = -20.0,60.0
  CUTWIN_TRESTMIN = -200.0,10.0
  CUTWIN_TRESTMAX = 5,99
  CUTWIN_TRESTRANGE = 10,999
  CUTWIN_MWEBV = 0.0,0.25
  CUTWIN_SNRMAX = 5.0,100000000.0
  CUTWIN_NFILT_SNRMAX = 2.0,99.0
&END

&FITINP
  FITMODEL_NAME = 'SALT3.K21'
  FILTLIST_FIT = '<filters>'
  FUDGEALL_ITER1_MAXFRAC = 0.01
  PRIOR_MJDSIG = 10
  PRIOR_SHAPE_RANGE = -4.0,4.0
  FITWIN_SHAPE = -3.0,3.0
  FITWIN_COLOR = -0.3,0.3
  FITWIN_TREST = -15.0,45.0
&END

"""

class run_snana:
    def __init__(self):
        self.debug = False
        
    def add_args(self,parser=None, usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage, conflict_handler="resolve")

        parser.add_argument('lcdir', type=str, default=None,
                            help='directory containing the lightcurve files')
        parser.add_argument('kcorfile', type=str, default=None,
                            help='kcor filename')
        parser.add_argument('filters', type=str, default=None,
                            help='filters')
        parser.add_argument('-o','--outputdir', type=str, default='output',
                            help='output directory for the lightcurve fit parameters and distances')
        parser.add_argument('-d','--datadir', type=str, default='whatever',
                            help='parent directory where your data lives')

        return parser

    def make_nml_file(self):

        lcversion = self.options.lcdir.split('/')[-1]
        
        nmltext = _nml_template.replace('<filters>',self.options.filters).\
            replace('<version>',lcversion).\
            replace('<output_dir>',self.options.outputdir).\
            replace('<kcorfile>',self.options.kcorfile).\
            replace('<data_dir>',self.options.datadir)

        outputnml = f'{self.options.outputdir}/{lcversion}.NML'
        with open(outputnml,'w') as fout:
            print(nmltext,file=fout)

        return outputnml
            
    def main(self):
        # the steps that run lightcurve fitting

        # make the LC fitting input file
        nmlfilename = self.make_nml_file()

        # do the LC fitting
        print(f"running: snlc_fit.exe {nmlfilename}")
        #os.system(f"snlc_fit.exe {nmlfilename}")

        # plot the lightcurves
        lcversion = self.options.lcdir.split('/')[-1]
        print(f"running: python plot_snana.py -v {lcversion} -a -20 -f {nmlfilename} --plotAll")
        #os.system(f"python plot_snana.py -v {lcversion} -a -20 -f {nmlfilename} --plotAll")

        # get the distances
        print(f"running: SALT2mu.exe salt2mu.default file={self.options.outputdir}/{lcversion}.FITRES.TEXT")
        os.system(f"SALT2mu.exe salt2mu.default file={self.options.outputdir}/{lcversion}.FITRES.TEXT")

        fitresfile = f"{self.options.outputdir}/{lcversion}.FITRES.TEXT"
        data = at.Table.read(fitresfile,format='ascii')
        # and then you can plot the data with, e.g.
        #plt.plot(data['zHD'],data['mB'] + 0.14*data['x1'] - 3.1*data['c'] + 19.36,'o')
        
        
if __name__ == "__main__":

    usagestring = "python run_snana.py <lightcurve_directory> <kcorfile>"

    rs = run_snana()
    parser = rs.add_args()
    args = parser.parse_args()
    rs.options = args
    
    rs.main()
