# preliminaries
import numpy as np
import os
import matplotlib.pyplot as plt
import sys  
sys.path.insert(0, r'C:\Users\payto\hubbleconst')

import sncosmo # we'll use this to get SN distances
import snana # we'll use this to read supernova light curve files
from astropy.io import fits
import astropy.units as u
from astropy.table import Table
import gzip
import warnings
import get_vpec
import LC_CLASS
# let's get ready to measure distances relative to
# a cosmological model with H0=70, cosmic matter = 0.3
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(70,0.3)

from os import listdir
from os.path import isfile, join
path = 'SNDATA_ROOT/research_2nd_sample/'
fileslist = [f for f in listdir(path) if isfile(join(path, f))]

sns = []
for path in fileslist:
    file = os.path.expandvars('SNDATA_ROOT/research_2nd_sample/'+ path)
    print(file)
    sn_ = snana.SuperNova(file)
    sns.append(sn_)
CSP_1_dict = ['2007jd','2007ai','2006hx','2006D','2004gc','2005hj','2009al','2007sr','2008hj','2009aa','2007ol','2007on','2007hu','2007st','2007as','2006is','2008ff','2005W','2005bg','2006ef','2008ia','2005al','2007le','2006hb','2005ku','2008bc','2008hv','2007nq','2005be','2007ba','2005lu','2008O','2008fl','2008fw','2008bi','2005na','2009Y','2005am','2008cc','2004gs','2007hj','2008fu','2007jg','2005mc','2008gg','2007jh','2008go','2006gt','2009ab','2008R','2006ot','2009ag','2009dc','2006eq','2009cz','2007al','2008hu']
CSP_2_dict = ['2008gp','2006bh','2005M','2006lu','2006et','2006ev','2005bo','2004ef','2004ey','2006ej','2008bq','2006py','2004gu']

def procedure(sn):
    kcordict = {'CFA3':os.path.expandvars('$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_Hicken2009.fits'),
                'CSP_1':os.path.expandvars('$SNDATA_ROOT/kcor/CSP/CSPDR3/kcor_CSPDR3_BD17.fits.gz'),
                'CSP_2':os.path.expandvars('$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_CSPDR2.fits'),
                'CFA4p1':os.path.expandvars('$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_CFA4p1.fits'),
                'CfA4p1':os.path.expandvars('$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_CFA4p1.fits'),
                'LOWZ':os.path.expandvars('$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_OTHER_LOWZ.fits'),
                'SDSS':os.path.expandvars('$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_SDSS.fits'),
                'CFA3S':os.path.expandvars('$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_Hicken2009.fits'),
                'CfA3S':os.path.expandvars('$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_Hicken2009.fits'),
                'CFA3K':os.path.expandvars('$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_Hicken2009.fits'),
                'CfA3K':os.path.expandvars('$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_Hicken2009.fits'),
                'CfA4p2':os.path.expandvars('$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_CfA4p2.fits'),
                'CFA4p2':os.path.expandvars('$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_CfA4p2.fits'),
                'CfA1':os.path.expandvars('$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_Riess1999.fits'),
                'CFA1':os.path.expandvars('$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_Riess1999.fits'),
                'CFA2':os.path.expandvars('$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_Jha2006.fits'),
                'CfA2':os.path.expandvars('$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_Jha2006.fits'),
                'PS1MD':os.path.expandvars('$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_PS1MD.fits'),
                'FOUNDATION_DR1':os.path.expandvars('$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_Foundation_DR1.fits'),
                'FOUNDATION':os.path.expandvars('$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_Foundation_DR1.fits'),
                'CSPDR2':os.path.expandvars('$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_CSPDR2.fits')}
    SURVEY_KEY = sn.SURVEY.split('(')[0]
    if SURVEY_KEY == 'CSP':
        if sn.SNID in CSP_1_dict:
            SURVEY_KEY = 'CSP_1'
        elif sn.SNID in CSP_2_dict:
            SURVEY_KEY = 'CSP_2'
        else:
            raise 
    if SURVEY_KEY not in kcordict:
        print(f'COULD NOT FIND {sn.SURVEY} IN kcordict USING KEY {SURVEY_KEY}')
        return    
    kcorfile = kcordict[SURVEY_KEY]
    kcor = fits.open(kcorfile)

    #for i,k in enumerate(kcor):
        #if i==0:
            #continue
        #print(k.data.names)
    print(sn.SURVEY,kcor[5].data.names)
        # cosmetic things
    plt.rcParams['figure.figsize'] = (4,12)
    plt.subplots_adjust(hspace=0.1)

    # plot the SN "light curve" for each band
    for i,filt in enumerate(np.unique(sn.FLT)):
        # subplot for each filter
        ax = plt.subplot(len(np.unique(sn.FLT)),1,i+1)

        # fluxes and errors
        print(filt)
        ax.errorbar(sn.MJD[sn.FLT == filt],sn.FLUXCAL[sn.FLT == filt],yerr=sn.FLUXCALERR[sn.FLT == filt],fmt='o')

        # compute the effective wavelength so we know the "color" of each filter
        if filt not in kcor[5].data.names and len(filt) == 1:
            for name in kcor[5].data.names:
                if name[-1] == filt:
                    kcorfilt = name
        else:
            kcorfilt = filt
            
        lameff = np.sum(kcor[5].data[kcorfilt]*kcor[5].data['wavelength (A)'])/kcor[5].data[kcorfilt].sum()
        # using ax.text instead of title makes the spacing better
        ax.text(0.05,0.25,f"{filt} ($\lambda_{{eff}} = {lameff:.0f} \AA$)",ha='left',va='center',transform=ax.transAxes,bbox={'facecolor':'1.0','edgecolor':'1.0'})

        # cosmetic stuff
        ax.set_xlabel('MJD')
        ax.set_ylabel('flux')
        ax.set_yticks([]) # absolute flux doesn't matter for now
    plt.show()
    flc=LC_CLASS.LC()
    res = LC_CLASS.salt3_lc_fit(flc, sn.datfile, kcordict[SURVEY_KEY])
    return res
reses = []
for sn in sns:
    print(sn.SNID)
    reses.append(procedure(sn))

x0s = []
x1s= []
cs = []
i = 0
for sn in sns:
    x0 = reses[i]['parameters'][np.array(reses[i]['param_names']) == 'x0']
    x1 = reses[i]['parameters'][np.array(reses[i]['param_names']) == 'x1']
    c = reses[i]['parameters'][np.array(reses[i]['param_names']) == 'c']
    x0s.append(x0)
    x1s.append(x1)
    cs.append(c)
    print(sn.SNID)
    print(x0,x1,c)
    i+=1

mus = []
salt3alpha = 0.14
salt3beta = 3.1
i = 0
#for res in reses:
for sn in sns:
    mu = -2.5*np.log10(x0s[i]) + 10.635 + salt3alpha*x1s[i] - salt3beta*cs[i] + 19.635
    mus.append(mu[0])
    print(sn.SNID)
    #print(mu)
    print(f"mu = {mu[0]:.3f} mag")
    i +=1

muerrs = []
for i,sn in enumerate(sns):
    muerr = np.sqrt((2.5/np.log(10)*reses[i]['errors']['x0']/x0)**2.+ 
                    salt3alpha**2. * reses[i]['errors']['x1']**2. + 
                    salt3beta**2. * reses[i]['errors']['c']**2.)
    muerrs.append(muerr[0])
    print(sn.SNID)
    print(f"mu_err = {muerr[0]:.3f} mag")

sns_RA=[]
sns_DEC = []
for i,sn in enumerate(sns):
    if isinstance(sn.RA,str):
        sn.RA = float(sn.RA.split()[0])
    sns_RA.append(sn.RA)
    i+=1
print(sns_RA)

for i,sn in enumerate(sns):
    if 'DECL' in sn.__dict__.keys():
        sn.__dict__["DEC"] = sn.__dict__['DECL']
        del sn.__dict__['DECL']

print(sn.__dict__.keys())
print(sn.__dict__)

for i,sn in enumerate(sns):
    if isinstance(sn.DEC,str):
        sn.DEC = float(sn.DEC.split()[0])
    sns_DEC.append(sn.DEC)
    i+=1
    #if i>n: 
        #break
        
print(sns_DEC)

zcmbs = []
zs = []
for i,sn in enumerate(sns):
    z = reses[i]['parameters'][np.array(reses[i]['param_names']) == 'z'][0]
    zs.append(z)
    #if isinstance(z,str):
        #z = float(z.split()[0])
    #zs.append(z)
zs

type(sns_RA)
sns_RA = np.array(sns_RA,dtype= np.float64)
sns_DEC = np.array(sns_DEC,dtype = np.float64)
zs = np.array(zs)

for z,sn in zip(zcmbs,sns):
    print(sn.SNID,z)

# convert to CMB redshift and apply peculiar velocities
vpec_mapfile_default = os.path.expandvars('$SNDATA_ROOT/models/VPEC/twomass++_velocity_LH11.npy')
vc = get_vpec.VelocityCorrection(vpec_mapfile_default)
for i,sn in enumerate(sns):
    zcmb = vc.convert_helio_to_cmb(sns_RA[i],sns_DEC[i],zs[i])
    zcmbs.append(zcmb)
#vpec,vsys = get_vpec.main(sns_RA,sns_DEC,zcmb,vpec_mapfile=vpec_mapfile_default)
# we'll use vpec later, omitted for now

    print(f"cosmological distance is {cosmo.distmod(zcmbs[i]).value:.3f} mag")

#fig, axs = plt.subplots(2)
fig, axs = plt.subplots(2, 1,gridspec_kw={'height_ratios':[4,1]})
fig.set_figwidth(10)
z = np.linspace(0,0.1,100)
zerr = 250/3e5*5.0/np.log(10)*(1.0+z)/(z*(1.0+z/2.0))
#plt.plot(z,cosmo.distmod(z).value, c = 'blue')
avg = np.median(mus - cosmo.distmod(zcmbs).value,)
residuals = (mus- cosmo.distmod(zcmbs).value) - avg
#axs[1].rcParams\n",
axs[0].plot(z,cosmo.distmod(z).value+avg, c = 'orange')
axs[0].plot(z,cosmo.distmod(z).value+avg+zerr, '--' , c = 'purple')
axs[0].plot(z,cosmo.distmod(z).value+avg-zerr, '--', c = 'purple')
axs[0].scatter(zcmbs,mus,marker = '.', s=.1)
#axs[1].scatter(residuals, mus)
axs[1].scatter(zcmbs,residuals)
axs[1].plot(z,zerr, '--' , c = 'purple')
axs[1].plot(z,-zerr, '--', c = 'purple')
plt.rcParams['figure.figsize'] = (4,12)
axs[0].errorbar(zcmbs,mus,yerr=muerrs,fmt = '.')
axs[1].errorbar(zcmbs,residuals,yerr=muerrs,fmt = '.')
axs[0].set_xlim([0.00,0.11])
#axs[1].set_xlim([0.005,0.02])
axs[0].set_ylim([32,39])
#axs[1].set_ylim([-0.5,0.5])
axs[0].set_xlabel('Redshift')
axs[0].set_ylabel('Distance (Mag)')
axs[1].set_xlabel('Redshift')
axs[1].set_ylabel('Hubble Residual')
axs[1].axhline(0, c = 'orange')
plt.show()
