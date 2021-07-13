#!/usr/bin/env python
# Get SN list, files, surveys, kcors from latest PantheonPlus FITRES file

import numpy as np
from txtobj import txtobj
import os
import getmu
import glob
import snana
import astropy.units as u
import astropy.table as at

from dustmaps.sfd import SFDQuery
sfd = SFDQuery()
from astropy.coordinates import SkyCoord

# let's get the survey names
surveynamedict = {}
with open(os.path.expandvars('$SNDATA_ROOT/SURVEY.DEF')) as fin:
    for line in fin:
        line = line.replace('\n','')
        if line.startswith('SURVEY:'):
            surveynamedict[line.split()[2]] = line.split()[1]

def cfa3_naturalsys():

    # keplercam, 4shooter
    filterdict = {'1':('CFA3K-U/f','CFA3S-U/a'),
                  '2':('CFA3K-B/h','CFA3S-B/b'),
                  '3':('CFA3K-V/j','CFA3S-V/c'),
                  '4':('','CFA3S-R/d'),
                  '5':('','CFA3S-I/e'),
                  '13':('CFA3K-r/k',''),
                  '14':('CFA3K-i/l','')}
    
    # loop through the file
    snid = ''
    with open('cfa3lightcurves.naturalsystem.txt') as fin:
        for line in fin:
            line = line.replace('\n','')
            if snid and not line.startswith('sn'):
                filt += [line.split()[0]]
                mjd += [float(line.split()[1])]
                mag += [float(line.split()[2])]
                magerr += [float(line.split()[3])]
            
            elif line.startswith('sn'):
                newsnid = f'20{line[2:]}'
                if snid == '':
                    snid = newsnid[:]
                    filt,mjd,mag,magerr = [],[],[],[]
                    continue
                
                # get the pantheon version
                panfile = glob.glob(os.path.expandvars(f'$SNDATA_ROOT/lcmerge/Pantheon_LOWZ_TEXT/CFA3*{snid}*'))[0]
                sn = snana.SuperNova(panfile)

                # make sure something hasn't gone horribly wrong
                if sn.SNID != snid: raise RuntimeError(f'SN file mismatch for {snid}!')
                if sn.SURVEY == 'PS1_LOWZ_COMBINED(CFA3K)':
                    dirname = 'CfA3_NATSYS_KEPLERCAM'
                elif sn.SURVEY == 'PS1_LOWZ_COMBINED(CFA3S)':
                    dirname = 'CfA3_NATSYS_4SHOOTER2'
                else:
                    print(sn.SNID,sn.SURVEY)
                    snid = newsnid[:]
                    filt,mjd,mag,magerr = [],[],[],[]
                    continue
                    
                    #raise RuntimeError(f'Survey {sn.SURVEY} isn\'t CfA3!')


                        
                fout = open(f'{dirname}/{sn.SNID}.dat','w')
                # use the pantheon header

                with open(panfile) as fin2:
                    for line2 in fin2:
                        line2 = line2.replace('\n','')
                        if line2.startswith('SNID:') or \
                           line2.startswith('IAUC:') or \
                           line2.startswith('RA:') or \
                           line2.startswith('DECL:') or \
                           line2.startswith('REDSHIFT_HELIO:') or \
                           line2.startswith('REDSHIFT_CMB:') or \
                           line2.startswith('REDSHIFT_FINAL:') or \
                           line2.startswith('HOSTGAL_LOGMASS:') or \
                           line2.startswith('PRIVATE(COSMO):') or \
                           line2.startswith('SEARCH_PEAKMJD:'):
                            print(line2,file=fout)
                        elif line2.startswith('SURVEY:'):
                            if 'KEPLERCAM' in dirname:
                                print('SURVEY: CfA3K',file=fout)
                            else:
                                print('SURVEY: CfA3S',file=fout)
                        elif line2.startswith('MWEBV:'):
                            sc = SkyCoord(sn.RA,sn.DECL,unit=u.deg)
                            ebv = sfd(sc)*0.86
                            print(f'MWEBV: {ebv:.3f} +- {ebv*0.05:.3f}',file=fout)
                        elif line2.startswith('FILTERS:'):
                            if 'KEPLERCAM' in dirname:
                                print('FILTERS: fhjkl',file=fout)
                            else:
                                print('FILTERS: abcde',file=fout)
                print(f"""NOBS: {len(mjd)}
NVAR: 7
VARLIST:  MJD  FLT FIELD   FLUXCAL   FLUXCALERR  MAG    MAGERR
#=============================================================""",file=fout)

                for f,m,mg,mge in zip(filt,mjd,mag,magerr):
                    flux = 10**(-0.4*(mg-27.5))
                    fluxerr = flux*0.4*np.log(10)*mge

                    if 'KEPLERCAM' in dirname: myfilt = filterdict[f][0]
                    else: myfilt = filterdict[f][1]

                    print(f'OBS: {m:.5f} {myfilt} VOID {flux:.3f} {fluxerr:.3f} {mg:.3f} {mge:.3f}',file=fout)
                print('END:',file=fout)
                filt,mjd,mag,magerr = [],[],[],[]
                
                snid = newsnid[:]


    # check MW E(B-V)
    
    # write the LC, use SALT3TRAIN filter codes
    # send to keplercam or 4shooter dirs depending on pantheon directory path

def cfa4_naturalsys():

    # keplercam, 4shooter
    filterdict = {'U':('CFA41-U','CFA42-U'),
                  'B':('CFA41-B','CFA42-B'),
                  'V':('CFA41-V','CFA42-V'),
                  "r'":('CFA41-r','CFA42-r'),
                  "i'":('CFA41-i','CFA42-i')}
    
    # loop through the file
    data = at.Table.read('cfa4.lc.natsystem.sort.ascii',format='ascii')
    for snid in np.unique(data['SN']):
        iSN = data['SN'] == snid

        nobs = len(data[iSN])

        # get the pantheon version
        panfile = glob.glob(os.path.expandvars(f'$SNDATA_ROOT/lcmerge/Pantheon_LOWZ_TEXT/CFA4*{snid}*'))
        if not len(panfile): continue
        panfile = panfile[0]
        sn = snana.SuperNova(panfile)

        # make sure something hasn't gone horribly wrong
        if sn.SNID != snid: raise RuntimeError(f'SN file mismatch for {snid}!')
        if sn.SURVEY == 'PS1_LOWZ_COMBINED(CFA4p1)':
            dirname = 'CfA4p1'
        elif sn.SURVEY == 'PS1_LOWZ_COMBINED(CFA4p2)':
            dirname = 'CfA4p2'
        
        filt,mjd,mag,magerr = [],[],[],[]
        for d in data[iSN]:
            if d['Filter'] == "u'": continue
            filt += [d['Filter']]
            mjd += [d['MJD']]
            mag += [d['Mag']]
            magerr += [d['sigma']]

                        
        with open(f'{dirname}/{sn.SNID}.dat','w') as fout:
        # use the pantheon header

            with open(panfile) as fin2:
                for line2 in fin2:
                    line2 = line2.replace('\n','')
                    if line2.startswith('SNID:') or \
                       line2.startswith('IAUC:') or \
                       line2.startswith('RA:') or \
                       line2.startswith('DECL:') or \
                       line2.startswith('REDSHIFT_HELIO:') or \
                       line2.startswith('REDSHIFT_CMB:') or \
                       line2.startswith('REDSHIFT_FINAL:') or \
                       line2.startswith('HOSTGAL_LOGMASS:') or \
                       line2.startswith('PRIVATE(COSMO):') or \
                       line2.startswith('SEARCH_PEAKMJD:'):
                        print(line2,file=fout)
                    elif line2.startswith('SURVEY:'):
                        if 'p2' in dirname:
                            print('SURVEY: CfA4p2',file=fout)
                        else:
                            print('SURVEY: CfA4p1',file=fout)
                    elif line2.startswith('MWEBV:'):
                        sc = SkyCoord(sn.RA,sn.DECL,unit=u.deg)
                        ebv = sfd(sc)*0.86
                        print(f'MWEBV: {ebv:.3f} +- {ebv*0.05:.3f}',file=fout)
                    elif line2.startswith('FILTERS:'):
                        if 'p1' in dirname:
                            print('FILTERS: UBVRI',file=fout)
                        else:
                            print('FILTERS: UBVRI',file=fout)
            print(f"""NOBS: {len(mjd)}
NVAR: 7
VARLIST:  MJD  FLT FIELD   FLUXCAL   FLUXCALERR  MAG    MAGERR
#=============================================================""",file=fout)

            for f,m,mg,mge in zip(filt,mjd,mag,magerr):
                flux = 10**(-0.4*(mg-27.5))
                fluxerr = flux*0.4*np.log(10)*mge

                if 'p1' in dirname: myfilt = filterdict[f][0]
                else: myfilt = filterdict[f][1]

                print(f'OBS: {m:.5f} {myfilt} VOID {flux:.4f} {fluxerr:.4f} {mg:.4f} {mge:.4f}',file=fout)
            print('END:',file=fout)
            #import pdb; pdb.set_trace()

    # check MW E(B-V)
    
    # write the LC, use SALT3TRAIN filter codes
    # send to keplercam or 4shooter dirs depending on pantheon directory path

def CfA12_conv():

    filtdict = {'B':'Bessell-U/U',
                'C':'Bessell-B/B',
                'D':'Bessell-V/V',
                'E':'Bessell-R/R',
                'F':'Bessell-I/I',
                'G':'Bessell-U/U',
                'H':'Bessell-B/B',
                'I':'Bessell-V/V',
                'J':'Bessell-R/R',
                'K':'Bessell-I/I'}
    
    files = glob.glob(os.path.expandvars('$SNDATA_ROOT/lcmerge/Pantheon_LOWZ_TEXT/LOWZ_JRK07*'))
    for f in files:
        sn = snana.SuperNova(f)

        if sn.SURVEY == 'PS1_LOWZ_COMBINED(CFA1)': dirname = 'CfA1'
        elif sn.SURVEY == 'PS1_LOWZ_COMBINED(CFA2)': dirname = 'CfA2'
        else: continue
        
        with open(f) as fin, open(f"{dirname}/{sn.SNID}.dat",'w') as fout:
            for line in fin:
                line = line.replace('\n','')
                if line.startswith('OBS'):
                    lineparts = line.split()
                    lineparts[2] = filtdict[lineparts[2]]
                    print(' '.join(lineparts),file=fout)
                elif line.startswith('SURVEY:'):
                    print(f'SURVEY: {dirname}',file=fout)
                elif line.startswith('FILTERS:'):
                    print('FILTERS: UBVRI',file=fout)
                elif line.startswith('MWEBV:'):
                    sc = SkyCoord(sn.RA,sn.DECL,unit=u.deg)
                    ebv = sfd(sc)*0.86
                    print(f'MWEBV: {ebv:.3f} +- {ebv*0.05:.3f}',file=fout)
                else:
                    print(line,file=fout)
            
def main():

    fr = txtobj('PanPlus_Latest.FITRES',fitresheader=True)
    fr = getmu.mkcuts(fr,fitprobmin=0)
    iz = fr.zHD < 0.1

    count = 0
    snid_outlist = []
    filelist = []
    kcorlist = []
    surveylist = []
    with open('sn_HF_list.txt','w') as fout:
        print('# SNID survey kcorfile snfile',file=fout)
        for j,i in enumerate(fr.CID):
            if fr.zHD[j] > 0.1: continue
            if fr.IDSURVEY[j] not in [50,51,52,56,57,58,59]:
                # let's stay away from KAIT, Swift, SupernovaFactory, and generic low-z

                # first look in the SALT3 directories:
                if 'SDSS' in surveynamedict[str(int(fr.IDSURVEY[j]))]:
                    salt3file = glob.glob(os.path.expandvars(f'$SNDATA_ROOT/lcmerge/SALT3TRAIN_K21/*SDSS*/*{fr.CID[j]}*'))
                else:
                    salt3file = glob.glob(os.path.expandvars(f'$SNDATA_ROOT/lcmerge/SALT3TRAIN_K21/*/*{fr.CID[j]}*'))
                if len(salt3file) == 1:  
                    filelist += [salt3file[0]]
                    snid_outlist += [fr.CID[j]]
                    sn = snana.SuperNova(salt3file[0])
                    surveylist += [sn.SURVEY]
                    kcorlist += [f'$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_{sn.SURVEY}.fits']
                    continue
                if len(salt3file) > 1: import pdb; pdb.set_trace()
                
                # next, look in CSPDR3
                cspfile = glob.glob(os.path.expandvars(f'$SNDATA_ROOT/lcmerge/CSPDR3/*{fr.CID[j]}*'))
                if len(cspfile) == 1:
                    filelist += [cspfile[0]]
                    snid_outlist += [fr.CID[j]]
                    sn = snana.SuperNova(cspfile[0])
                    surveylist += [sn.SURVEY]

                    kcorlist += ['$SNDATA_ROOT/kcor/CSP/CSPDR3/kcor_CSPDR3_BD17.fits.gz']
                    continue

                # any stragglers in Foundation?
                foundfile = glob.glob(os.path.expandvars(f'$SNDATA_ROOT/lcmerge/Foundation_DR1/*{fr.CID[j]}*'))
                if len(foundfile) == 1:
                    filelist += [foundfile[0]]
                    sn = snana.SuperNova(foundfile[0])
                    surveylist += [sn.SURVEY]

                    snid_outlist += [fr.CID[j]]
                    kcorlist += ['$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_Foundation_DR1.fits']
                    continue

                # PS1
                ps1file = glob.glob(os.path.expandvars(f'$SNDATA_ROOT/lcmerge/Pantheon_PS1MD_TEXT/*{fr.CID[j]}*'))
                if len(ps1file) == 1:
                    filelist += [ps1file[0]]
                    sn = snana.SuperNova(ps1file[0])
                    surveylist += [sn.SURVEY]
                    
                    snid_outlist += [fr.CID[j]]
                    kcorlist += ['$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_PS1MD.fits']
                    continue

                # SDSS
                #sdssfile = glob.glob(os.path.expandvars(f'$SNDATA_ROOT/lcmerge/SDSS_HOLTZ08/*{fr.CID[j]}*'))
                #if len(sdssfile) == 1:
                #    filelist += [sdssfile[0]]
                #    sn = snana.SuperNova(sdssfile[0])
                #    surveylist += [sn.SURVEY]

                #    snid_outlist += [fr.CID[j]]
                #    kcorlist += ['$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_SDSS.fits']
                #    continue

                
                # now let's check CfA, but make sure we have the natural system LCs
                # good directories should be...
                cfafile = glob.glob(os.path.expandvars(f'CfA3_NATSYS*/*{fr.CID[j]}*'))
                if len(cfafile) == 1:
                    #print(fr.CID[j])
                    filelist += [cfafile[0]]
                    sn = snana.SuperNova(cfafile[0])
                    surveylist += [sn.SURVEY]
                    
                    snid_outlist += [fr.CID[j]]
                    kcorlist += ['$SNDATA_ROOT/kcor/SALT3TRAIN_K21/kcor_Hicken2009.fits']
                    continue

                # finally, other CfA
                cfafile = glob.glob(os.path.expandvars(f'CfA*/*{fr.CID[j]}*'))
                if len(cfafile) == 1:
                    #print(fr.CID[j])
                    filelist += [cfafile[0]]
                    sn = snana.SuperNova(cfafile[0])
                    surveylist += [sn.SURVEY]

                    snid_outlist += [fr.CID[j]]
                    if sn.SURVEY == 'CfA1':
                        kcorlist += ['kcor_Riess1999.fits']
                    elif sn.SURVEY == 'CfA2':
                        kcorlist += ['kcor_Jha2006.fits']
                    elif sn.SURVEY.startswith('CfA4'):
                        kcorlist += [f'kcor_{sn.SURVEY}.fits']
                    else: raise RuntimeError(f'survey {sn.SURVEY} unknown!')
                    
                    continue

                print(f'SN {fr.CID[j]} not found!')
                #import pdb; pdb.set_trace()
                continue
    
        print(len(np.unique(snid_outlist)))
        print(len(snid_outlist),len(surveylist),len(kcorlist),len(filelist))
        for snid,survey,kcor,filen in zip(snid_outlist,surveylist,kcorlist,filelist):
            print(f'{snid} {survey} {kcor} {filen}',file=fout)
    
    import pdb; pdb.set_trace()
    
if __name__ == "__main__":
    #cfa3_naturalsys()
    #cfa4_naturalsys()
    #CfA12_conv()
    main()
