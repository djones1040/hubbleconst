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
path1 = 'SNDATA_ROOT/research_payto/'
fileslist1 = [f for f in listdir(path1) if isfile(join(path1, f))]

from os import listdir
from os.path import isfile, join
path2 = 'SNDATA_ROOT/research_2nd_sample/'
fileslist2 = [f for f in listdir(path2) if isfile(join(path2, f))]

sns1 = []
for path1 in fileslist1:
    file1 = os.path.expandvars('SNDATA_ROOT/research_payto/'+ path1)
    print(file1)
    sn_1 = snana.SuperNova(file1)
    sns1.append(sn_1)

sns2 = []
for path2 in fileslist2:
    file2 = os.path.expandvars('SNDATA_ROOT/research_2nd_sample/'+ path2)
    print(file2)
    sn_2 = snana.SuperNova(file2)
    sns2.append(sn_2)

CSP_1_dict = ['2007jd','2007ai','2006hx','2006D','2004gc','2005hj','2009al','2007sr','2008hj','2009aa','2007ol','2007on','2007hu','2007st','2007as','2006is','2008ff','2005W','2005bg','2006ef','2008ia','2005al','2007le','2006hb','2005ku','2008bc','2008hv','2007nq','2005be','2007ba','2005lu','2008O','2008fl','2008fw','2008bi','2005na','2009Y','2005am','2008cc','2004gs','2007hj','2008fu','2007jg','2005mc','2008gg','2007jh','2008go','2006gt','2009ab','2008R','2006ot','2009ag','2009dc','2006eq','2009cz','2007al','2008hu','2006bh']
CSP_2_dict = ['2008gp','2005M','2006lu','2006et','2006ev','2005bo','2004ef','2004ey','2006ej','2008bq','2006py','2004gu',]



def procedure1(sn):
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
                'PS1_LOWZ_COMBINED':os.path.expandvars('$SNDATA_ROOT/kcor/PS1/Pantheon/kcor_PS1_LOWZ_COMBINED.fits'),
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
reses1 = []
for sn in sns1:
    print(sn.SNID)
    reses1.append(procedure1(sn))
    

def procedure2(sn):
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
                'PS1_LOWZ_COMBINED':os.path.expandvars('$SNDATA_ROOT/kcor/PS1/Pantheon/kcor_PS1_LOWZ_COMBINED.fits'),
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
reses2 = []
for sn in sns2:
    print(sn.SNID)
    reses2.append(procedure2(sn))

x0s_1 = []
x1s_1= []
cs_1 = []
for i,sn in enumerate(sns1):
    x0_1 = reses1[i]['parameters'][np.array(reses1[i]['param_names']) == 'x0']
    x1_1 = reses1[i]['parameters'][np.array(reses1[i]['param_names']) == 'x1']
    c_1 = reses1[i]['parameters'][np.array(reses1[i]['param_names']) == 'c']
    x0s_1.append(x0_1)
    x1s_1.append(x1_1)
    cs_1.append(c_1)
    print(sn.SNID)
    print(x0_1,x1_1,c_1)

x0s_2 = []
x1s_2= []
cs_2 = []
for i,sn in enumerate(sns2):
    x0_2 = reses2[i]['parameters'][np.array(reses2[i]['param_names']) == 'x0']
    x1_2 = reses2[i]['parameters'][np.array(reses2[i]['param_names']) == 'x1']
    c_2 = reses2[i]['parameters'][np.array(reses2[i]['param_names']) == 'c']
    x0s_2.append(x0_2)
    x1s_2.append(x1_2)
    cs_2.append(c_2)
    print(sn.SNID)
    print(x0_2,x1_2,c_2)

mus_1 = []
salt3alpha = 0.14
salt3beta = 3.1
i = 0
#for res in reses:
for sn in sns1:
    mu_1 = -2.5*np.log10(x0s_1[i]) + 10.635 + salt3alpha*x1s_1[i] - salt3beta*cs_1[i] + 19.635
    mus_1.append(mu_1[0])
    print(sn.SNID)
    #print(mu)
    print(f"mu = {mu_1[0]:.3f} mag")
    i +=1

mus_2 = []
salt3alpha = 0.14
salt3beta = 3.1
i = 0
#for res in reses:
for sn in sns2:
    mu_2 = -2.5*np.log10(x0s_2[i]) + 10.635 + salt3alpha*x1s_2[i] - salt3beta*cs_2[i] + 19.635
    mus_2.append(mu_2[0])
    print(sn.SNID)
    #print(mu)
    print(f"mu = {mu_2[0]:.3f} mag")
    i +=1

muerrs_1 = []
for i,sn in enumerate(sns1):
    muerr_1 = np.sqrt((2.5/np.log(10)*reses1[i]['errors']['x0']/x0_1)**2.+ 
                    salt3alpha**2. * reses1[i]['errors']['x1']**2. + 
                    salt3beta**2. * reses1[i]['errors']['c']**2.)
    muerrs_1.append(muerr_1[0])
    print(sn.SNID)
    print(f"mu_err = {muerr_1[0]:.3f} mag")

muerrs_2 = []
for i,sn in enumerate(sns2):
    muerr_2 = np.sqrt((2.5/np.log(10)*reses2[i]['errors']['x0']/x0_2)**2.+ 
                    salt3alpha**2. * reses2[i]['errors']['x1']**2. + 
                    salt3beta**2. * reses2[i]['errors']['c']**2.)
    muerrs_2.append(muerr_2[0])
    print(sn.SNID)
    print(f"mu_err = {muerr_2[0]:.3f} mag")

sns_RA_1=[]
sns_DEC_1 = []
for i,sn in enumerate(sns1):
    if isinstance(sn.RA,str):
        sn.RA = float(sn.RA.split()[0])
    sns_RA_1.append(sn.RA)
    i+=1
print(sns_RA_1)

sns_RA_2=[]
sns_DEC_2 = []
for i,sn in enumerate(sns2):
    if isinstance(sn.RA,str):
        sn.RA = float(sn.RA.split()[0])
    sns_RA_2.append(sn.RA)
    i+=1
print(sns_RA_2)

for i,sn in enumerate(sns1):
    if 'DECL' in sn.__dict__.keys():
        sn.__dict__["DEC"] = sn.__dict__['DECL']
        del sn.__dict__['DECL']

for i,sn in enumerate(sns2):
    if 'DECL' in sn.__dict__.keys():
        sn.__dict__["DEC"] = sn.__dict__['DECL']
        del sn.__dict__['DECL']

for i,sn in enumerate(sns1):
    if isinstance(sn.DEC,str):
        sn.DEC = float(sn.DEC.split()[0])
    sns_DEC_1.append(sn.DEC)
    i+=1
    #if i>n: 
        #break
        
print(sns_DEC_1)

for i,sn in enumerate(sns2):
    if isinstance(sn.DEC,str):
        sn.DEC = float(sn.DEC.split()[0])
    sns_DEC_2.append(sn.DEC)
    i+=1
    #if i>n: 
        #break
        
print(sns_DEC_2)

zcmbs_1 = []
zs_1 = []
for i,sn in enumerate(sns1):
    z_1 = reses1[i]['parameters'][np.array(reses1[i]['param_names']) == 'z'][0]
    zs_1.append(z_1)
    #if isinstance(z,str):
        #z = float(z.split()[0])
    #zs.append(z)
zs_1

zcmbs_2 = []
zs_2 = []
for i,sn in enumerate(sns2):
    z_2 = reses2[i]['parameters'][np.array(reses2[i]['param_names']) == 'z'][0]
    zs_2.append(z_2)
    #if isinstance(z,str):
        #z = float(z.split()[0])
    #zs.append(z)
zs_2

type(sns_RA_1)
sns_RA_1 = np.array(sns_RA_1,dtype= np.float64)
sns_DEC_1 = np.array(sns_DEC_1,dtype = np.float64)
zs_1 = np.array(zs_1)

type(sns_RA_2)
sns_RA_2 = np.array(sns_RA_2,dtype= np.float64)
sns_DEC_2 = np.array(sns_DEC_2,dtype = np.float64)
zs_2 = np.array(zs_2)

# convert to CMB redshift and apply peculiar velocities
vpec_mapfile_default = os.path.expandvars('$SNDATA_ROOT/models/VPEC/twomass++_velocity_LH11.npy')
vc = get_vpec.VelocityCorrection(vpec_mapfile_default)
for i,sn in enumerate(sns1):
    zcmb_1 = vc.convert_helio_to_cmb(sns_RA_1[i],sns_DEC_1[i],zs_1[i])
    zcmbs_1.append(zcmb_1)
#vpec,vsys = get_vpec.main(sns_RA,sns_DEC,zcmb,vpec_mapfile=vpec_mapfile_default)
# we'll use vpec later, omitted for now

    print(f"cosmological distance is {cosmo.distmod(zcmbs_1[i]).value:.3f} mag")

# convert to CMB redshift and apply peculiar velocities
vpec_mapfile_default = os.path.expandvars('$SNDATA_ROOT/models/VPEC/twomass++_velocity_LH11.npy')
vc = get_vpec.VelocityCorrection(vpec_mapfile_default)
for i,sn in enumerate(sns2):
    zcmb_2 = vc.convert_helio_to_cmb(sns_RA_2[i],sns_DEC_2[i],zs_2[i])
    zcmbs_2.append(zcmb_2)
#vpec,vsys = get_vpec.main(sns_RA,sns_DEC,zcmb,vpec_mapfile=vpec_mapfile_default)
# we'll use vpec later, omitted for now

    print(f"cosmological distance is {cosmo.distmod(zcmbs_2[i]).value:.3f} mag")

#fig, axs = plt.subplots(2)
fig, axs = plt.subplots(2, 1,gridspec_kw={'height_ratios':[4,1]})
fig.set_figwidth(10)
z = np.linspace(0,0.1,100)
zerr = 250/3e5*5.0/np.log(10)*(1.0+z)/(z*(1.0+z/2.0))
#plt.plot(z,cosmo.distmod(z).value, c = 'blue')
avg1 = np.median(mus_1 - cosmo.distmod(zcmbs_1).value,)
avg2 = np.median(mus_2 - cosmo.distmod(zcmbs_2).value,)
residuals_1 = (mus_1- cosmo.distmod(zcmbs_1).value) - avg1
residuals_2 = (mus_2- cosmo.distmod(zcmbs_2).value) - avg2
#axs[1].rcParams\n",
axs[0].plot(z,cosmo.distmod(z).value+avg1, c = 'purple')
axs[0].plot(z,cosmo.distmod(z).value+avg1+zerr, '--' , c = 'purple')
axs[0].plot(z,cosmo.distmod(z).value+avg1-zerr, '--', c = 'purple')
axs[0].plot(z,cosmo.distmod(z).value+avg2, c = 'orange')
axs[0].plot(z,cosmo.distmod(z).value+avg2+zerr, '--' , c = 'orange')
axs[0].plot(z,cosmo.distmod(z).value+avg2-zerr, '--', c = 'orange')
axs[0].scatter(zcmbs_1,mus_1,marker = 'o', s=20, c = 'blue')
axs[0].scatter(zcmbs_2,mus_2,marker = '.', s=10,c = 'green')
#axs[1].scatter(residuals, mus)
axs[1].scatter(zcmbs_1,residuals_1,c = 'blue' , s=20)
axs[1].scatter(zcmbs_2,residuals_2, s = 10, marker = '.', c = 'green')
axs[1].plot(z,zerr, '--' , c = 'purple')
axs[1].plot(z,-zerr, '--', c = 'purple')
plt.rcParams['figure.figsize'] = (4,12)
axs[0].errorbar(zcmbs_1,mus_1,yerr=muerrs_1,fmt = 'o',c='blue')
axs[0].errorbar(zcmbs_2,mus_2,yerr=muerrs_2,fmt = '.',c='green')
axs[1].errorbar(zcmbs_1,residuals_1,yerr=muerrs_1,fmt = 'o',c='red',zorder = 10)
axs[1].errorbar(zcmbs_2,residuals_2,yerr=muerrs_2,fmt = '.',c='green')
axs[0].set_xlim([0.005,0.1])
axs[1].set_xlim([0.005,0.1])
axs[0].set_ylim([32,39])
axs[1].set_ylim([-0.5,0.5])
axs[0].set_xlabel('Redshift')
axs[0].set_ylabel('Distance (Mag)')
axs[1].set_xlabel('Redshift')
axs[1].set_ylabel('Hubble Residual')
axs[1].axhline(0, c = 'orange')
plt.show()