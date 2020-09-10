import matplotlib.text as text
import pandas as pd
import numpy as np
import pylab as pl
import scipy as sp
import re

from rpy2.robjects import pandas2ri
from rpy2.robjects.vectors import ListVector, FloatVector
import scipy.optimize as optz
import seaborn as sns
from scipy import stats
import rpy2.robjects.packages as rpackages
import seaborn as sns
from statsmodels import robust
from matplotlib.backends.backend_pdf import PdfPages

from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import rpy2.robjects as R

pandas2ri.activate()
httk = rpackages.importr('httk')

def solvePBTK(casrn,days=90,dose=None,doses_per_day=1,
              species='Rat',out_unit='uM',rest_clear=False):
    kwargs = {'chem.cas':casrn,
            'days':float(days),'plots':False,
            'suppress.messages':True,
            'species':species,
            'default.to.human':True,
            'out.units':out_unit,
            'daily.dose':float(dose),
            'restrictive.clearance':rest_clear,
            'doses.per.day':float(doses_per_day)}
    
    R0 = httk.solve_pbtk(**kwargs)
    return pd.DataFrame(R0,columns=['time', 'Agutlumen', 'Cart', 'Cven', 'Clung', 'Cgut', 'Cliver',
       'Ckidney', 'Crest', 'Ametabolized', 'Atubules', 'Cplasma', 'AUC_v'])


def calcOralEquiv(casrn,conc,q=0.5,species='Rat',units_in='uM',
                  units_out='mg',rest_clear=False):
    if type(q)==list:
        q=FloatVector(q)
    if type(conc)==list or type(conc)==pd.Series:
        conc=FloatVector(conc)
        
    kwargs={'conc':conc,
            'chem.cas':casrn,
            'which.quantile':q,
            'species':species,
            'input.units':units_in,
            'output.units':units_out,
            'restrictive.clearance':rest_clear,
            'suppress.messages':True
           }   
    X = httk.calc_mc_oral_equiv(**kwargs)
    #return pandas2ri.ri2py_listvector(X)
    return X

def getChemInfo(species='Rat',def_human=True):
    kwargs = {'info':"CAS",'species':species,
              'default.to.human':def_human}
    return list(httk.get_cheminfo(**kwargs))

def getQ(x,q0=0.5):
    oed = x.pop('oed')
    x['oed']=next(i['oral_dose'] for i in oed if i['q']==q0)
    return x



def solveOED(casrn,days=90,CR=None,InVitro=[],species='Rat',rest_clear=True):
    # CR is a dataframe containing columns 
    #  casrn
    #  timeh 
    #  in vitro uM critical concentrations given in InVitro
    
    Res=[]
    InVt0= InVitro
    CT = ['cmax','cave','auc']
    LD = list(CR.timeh.unique()/24)
    for d in [1e-2,1,5,10,50,1e2,1e3]:
        CXt = solvePBTK(casrn,days=days,dose=d,rest_clear=rest_clear,species=species)
        #CXt.loc[CXt.index[1:],'AUC_v']=sp.integrate.cumtrapz(CXt.Cven,CXt.time)

        for t in LD:
            d0 = days-t
            CXi = CXt[CXt.time>=d0]
            auc_t = CXi.AUC_v.max()
            cmax_t= CXi.Cven.max()
            cave_t= CXi.Cven.mean()

            Res.append(dict(oral_dose=d,last_days=t,auc=auc_t,
                            cmax=cmax_t,cave=cave_t))

    IVIV0 = pd.DataFrame(Res)
    return IVIV0

    OED = pd.DataFrame(np.zeros((3*len(LD),len(InVt0))),columns=InVt0)
    ind = [[i,j] for i in CT for j in LD]
    OED.index=pd.MultiIndex.from_arrays(list(zip(*ind)),names=['conc_type','last_days'])

    T = []
    #return OED,IVIV0
    for ld,(Iviv) in IVIV0.groupby(['last_days']):
        CRi = CR[(CR['timeh']==(ld*24))]
        #print t,Iviv.shape
        #T.append(t*24)
        for ct in CT:
            C = CRi[InVt0]
            if C.shape[0]==0:continue
            if ct=='auc': C *= ld
            OED.loc[(ct,ld)]=sp.interp(C,Iviv[ct],Iviv['oral_dose'])[0]
    
    OED = OED.reset_index()
    
    OED.loc[:,'timeh']=OED.last_days*24
    OED.loc[:,'casrn']=casrn
    OED.loc[:,'sim_days']=days
    
    return OED

def simIVIVEdata(casrn,species='Rat',
                 rest_clear=True,
                 doses_per_day=1,ivt_units='uM',
                OD=[1e-2,1,10,1e2,1e3],
                IVt_days=[],IVv_days=[],
                replace=True,
                dbc_sav=None):
    """
    Generate the data for IVIVE 
    
    days_ivv, dose_ivv, days_ivt, dose_ivt 
    
    In vitro conc of hepatocytes is venous conc in rat at time-matched 
    time points (time-matched)

    * Dose_ivt at ti = Cmax at ti 
    * Dose_ivt at ti = Cave at ti 
    * Dose_ivt * ti  = AUCv at ti 

    In vitro conc of hepatocytes is venous conc in rats time-averaged across 
    T (time-averaged)

    * Dose_ivt at ti = Cmax from 0 to Ti 
    * Dose_ivt at ti = Cave from 0 to Ti 
    * Dose_ivt * ti  = AUCv at Ti 

    In vitro conc of hepatocytes is venous conc in rats (final-time-matched)

    * Dose_ivt at ti = Cmax from Ti-ti:Ti 
    * Dose_ivt at ti = Cave from Ti-ti:Ti 
    * Dose_ivt * ti  = mean AUCv Ti-ti:Ti     
    """
    def mkTimeMatched():
        d0=days_ivt-1                
        CXi = CXt[(CXt.time>d0) & (CXt.time<=days_ivt)]
        auc_t = CXi.AUC_v.max()
        cmax_t= CXi.Cven.max()
        cave_t= CXi.Cven.mean()

        return dict(dose_ivv=dose_ivv,days_ivv=days_ivv,
                    days_ivt=days_ivt,
                    method='time-matched',
                    auc=auc_t,
                    cmax=cmax_t,
                    cave=cave_t)
        
    def mkTimeAveraged():
        CXi = CXt[(CXt.time<=days_ivv)]
        auc_t = CXi.AUC_v.max()
        cmax_t= CXi.Cven.max()
        cave_t= CXi.Cven.mean()

        return dict(dose_ivv=dose_ivv,days_ivv=days_ivv,
                    days_ivt=days_ivt,
                    method='time-averaged',
                    auc=auc_t,
                    cmax=cmax_t,
                    cave=cave_t)
        
    def mkTimeFinalMatched():
        d0=days_ivv-days_ivt               
        CXi = CXt[(CXt.time>d0) & (CXt.time<=days_ivv)]
        auc_t = CXi.AUC_v.max()
        cmax_t= CXi.Cven.max()
        cave_t= CXi.Cven.mean()

        return dict(dose_ivv=dose_ivv,days_ivv=days_ivv,
                    days_ivt=days_ivt,
                    method='time-final-matched',
                    auc=auc_t,
                    cmax=cmax_t,
                    cave=cave_t)
        
        
    Res=[]
    CT = ['cmax','cave','auc']
    
    for dose_ivv in OD:
        CXt = solvePBTK(casrn,days=np.max(IVv_days),dose=dose_ivv,
                        rest_clear=rest_clear,
                        out_unit=ivt_units,
                        doses_per_day=doses_per_day,
                        species=species)
        for days_ivv in IVv_days:
            for days_ivt in IVt_days:
                if days_ivt + 1 > days_ivv: continue
                Res.append(mkTimeMatched())
                Res.append(mkTimeAveraged())
                Res.append(mkTimeFinalMatched())

    IV0 = pd.DataFrame(Res)
    IV1 = IV0.melt(id_vars=['days_ivt','days_ivv','dose_ivv','method'],
                     var_name='dose_metric',
                     value_name='dose_ivt')
    IV1.insert(0,'casrn',casrn)
    IV1.insert(1,'species',species)
    
    if dbc_sav:
        old = dbc_sav.find_one(dict(casrn=casrn,species=species))
        if replace and old:
            dbc_sav.delete_one(old)
            
        dbc_sav.insert_one(dict(casrn=casrn,species=species,
                                rest_clear=rest_clear,
                                doses_per_day=doses_per_day,
                                ivt_units=ivt_units,
                                oral_doses=OD,
                                invitro_days=IVt_days,
                                invivo_days=IVv_days,
                                ivivd=IV1.to_dict('records')
                               ))
    else:
        return dict(casrn=casrn,species=species,
                                rest_clear=rest_clear,
                                doses_per_day=doses_per_day,
                                ivt_units=ivt_units,
                                oral_doses=OD,
                                invitro_days=IVt_days,
                                invivo_days=IVv_days,
                                ivivd=IV1.to_dict('records')
                               )

def getIVIVEdata(chm_id,species,
                 days_ivt,days_ivv,
                 dose_metric=None,
                 method=None,
                 dbc_iviv=None):
    """
    Get the IVIV data stored by simIVIVEdata. Currently, 
    assumes that only one record per casrn 
    
    One of the issues with matching by days_ivt is that the input and 
    the database stored time has to match at 3 dp
    """
    Q = {'$or':[{k:chm_id,'species':species} for k in ['casrn','dsstox_sid','name']]}
    D = dbc_iviv.find_one(Q)
    if not D:
        raise ValueError("No record for: chem_id="+str(chm_id)+
                         ",species="+str(species))
    IV = pd.DataFrame(D['ivivd'])
    IV.loc[:,'days_ivt'] = IV.days_ivt.round(3)
    q = "days_ivt==%f and days_ivv==%f" % (np.round(days_ivt,3),days_ivv)
    if dose_metric:
        q += "and dose_metric=='%s'" % dose_metric
    if method:
        q += "and method=='%s'" % method        
    
    IVi= IV.query(q)
    #for i in ['casrn','dsstox_sid','name']:
    #    if i not in IVi: IVi.insert(0,i,D[i])
    
    if IVi.shape[0]==0:
        raise ValueError("No data for:" + q)
    else:
        return IVi
        
def calcAED(chm_id,species,days_ivt,days_ivv,dose_ivt,
            dose_metric='auc',
            method='time-matched',
            dbc_iviv=None
           ):

    """
    Generate the AED using IVIVD:-
    
    Use IVIVD to make days_ivt, days_ivv, dose_metric 
    
    AED = F_chem(dose_invitro,day_invitro, day_invivo,dose_metric)
    
    Methods:
    day_invitro = ti
    day_invivo  = Ti
    
    In vitro conc of hepatocytes is venous conc in rat at time-matched 
    time points (time-matched)

        Dose_ivt at ti = Cmax at ti 
        Dose_ivt at ti = Cave at ti 
        Dose_ivt * ti  = AUCv at ti 

    In vitro conc of hepatocytes is venous conc in rats time-averaged across 
    T (time-averaged)

        Dose_ivt at ti = Cmax from 0 to Ti 
        Dose_ivt at ti = Cave from 0 to Ti 
        Dose_ivt * ti  = AUCv at Ti 

    In vitro conc of hepatocytes is venous conc in rats (final-time-matched)
        Dose_ivt at ti = Cmax from Ti-ti:Ti 
        Dose_ivt at ti = Cave from Ti-ti:Ti 
        Dose_ivt * ti  = mean AUCv Ti-ti:Ti 

 
    Note: 
    * For dose_metric == auc: dose_invitro *= day_invitro 
    
    """
    IVIVi = getIVIVEdata(chm_id,species,days_ivt,days_ivv,dose_metric,
                         method,dbc_iviv)

    if dose_metric=='auc':
        dose_ivt = dose_ivt * days_ivt

    AED = sp.interp(dose_ivt,IVIVi.dose_ivt,IVIVi.dose_ivv)
    
    return AED


def calcAllAED(chm_id,species,days_ivt,days_ivv,dose_ivt,
               css=False,
               dbc_iviv=None
           ):

    """
    Generate the AED using IVIVD:-
    
    Use IVIVD to make days_ivt, days_ivv, dose_metric 
    
    AED = F_chem(dose_invitro,day_invitro, day_invivo,dose_metric)
    
    Methods:
    day_invitro = ti
    day_invivo  = Ti
    
    In vitro conc of hepatocytes is venous conc in rat at time-matched 
    time points (time-matched)

        Dose_ivt at ti = Cmax at ti 
        Dose_ivt at ti = Cave at ti 
        Dose_ivt * ti  = AUCv at ti 

    In vitro conc of hepatocytes is venous conc in rats time-averaged across 
    T (time-averaged)

        Dose_ivt at ti = Cmax from 0 to Ti 
        Dose_ivt at ti = Cave from 0 to Ti 
        Dose_ivt * ti  = AUCv at Ti 

    In vitro conc of hepatocytes is venous conc in rats (final-time-matched)
        Dose_ivt at ti = Cmax from Ti-ti:Ti 
        Dose_ivt at ti = Cave from Ti-ti:Ti 
        Dose_ivt * ti  = mean AUCv Ti-ti:Ti 

 
    Note: 
    * For dose_metric == auc: dose_invitro *= day_invitro 
    
    """
    IVIV = getIVIVEdata(chm_id,species,days_ivt,days_ivv,dbc_iviv=dbc_iviv)
    AED = []
    
    for (name,dtxsid,casrn,method,dose_metric),IVIVi in \
            IVIV.groupby(['name','dsstox_sid','casrn','method','dose_metric']):
        if dose_metric=='auc':
            d_ivt = dose_ivt * days_ivt
        else: 
            d_ivt = dose_ivt            
        
        Aed = sp.interp(d_ivt,IVIVi.dose_ivt,IVIVi.dose_ivv)
        AED.append(pd.DataFrame(dict(name=name,casrn=casrn,
                                     dsstox_sid=dtxsid,
                                     days_ivt=days_ivt,days_ivv=days_ivv,
                                     ivive_method=method,dose_metric=dose_metric,
                                     dose_ivt=dose_ivt,aed=list(Aed))))
        
    if css:
        for q in [0.05,0.5,0.95]:
            Aed = calcOralEquiv(casrn,list(dose_ivt),q=q)
            AED.append(pd.DataFrame(dict(name=name,casrn=casrn,
                                     dsstox_sid=dtxsid,
                                     days_ivt=days_ivt,days_ivv=days_ivv,
                                     dose_metric='q%.2f'%q,ivive_method='css',
                                     dose_ivt=dose_ivt,aed=list(Aed))))
    AED1 = pd.concat(AED)
    return AED1



def calcIntD(chm_id,species,days_ivt,days_ivv,dose_ivv,
             rest_clear=True,
             dbc_iviv=None
             ):

    """
    Generate the internal dose using IVIVD:-
    
    Use IVIVD to make days_ivt, days_ivv, dose_metric 
    
    
    Note: 
    * For dose_metric == auc: dose_invitro *= day_invitro 
    
    """
    Y = simIVIVEdata(chm_id,species=species,OD=[min(dose_ivv),max(dose_ivv)],
                     rest_clear=rest_clear,IVt_days=[days_ivt],IVv_days=[days_ivv])
    IVIVi = pd.DataFrame(Y['ivivd'])
    
    Res = []
    for (dose_metric,method,days_ivt,days_ivv),X_i in \
        IVIVi.groupby(['dose_metric','method','days_ivt','days_ivv']):
        InD = sp.interp(dose_ivv,X_i.dose_ivv,X_i.dose_ivt)
        if dose_metric == 'auc':
            InD = InD/days_ivt
        Res.append(pd.DataFrame(dict(dose_metric=dose_metric,method=method,dose_ivv=dose_ivv,
                                     days_ivv=days_ivv,days_ivt=days_ivt,dose_internal=InD)))
        
    return pd.concat(Res)

def runCalcAllAED(name=None,ft=None,days_ivt=None,days_ivv=None,
                  potency_val=None,assay=None,platform=None,
                  potency_name=None,dbg=False,
                  potency_meth=None,
                  dbc_iviv=None,dbc_sav=None,dbc_err=None,
                  **kwds):
    if dbg: print(name,ft,days_ivt,days_ivv,potency_name,potency_val)
    R = []
    for di in days_ivv:
        try:
            Yi = calcAllAED(name,species='Rat',days_ivt=days_ivt,
                        days_ivv=di,
                        dose_ivt=pd.Series(potency_val),
                        css=True,
                        dbc_iviv=dbc_iviv)
        except:
            dbc_err.insert_one(dict(name=name,species='Rat',days_ivt=days_ivt,
                                    days_ivv=di,
                                    dose_ivt=potency_val,
                                    ft=ft,assay=assay,platform=platform,
                                    css=True))
        else:
            R.append(Yi)
    if len(R)==0: return
    
    Y = pd.concat(R)
    if potency_name: Y.insert(Y.shape[1],'potency_name',potency_name)
    if ft: Y.insert(Y.shape[1],'ft',ft)
    if assay: Y.insert(Y.shape[1],'assay',assay)
    if platform: Y.insert(Y.shape[1],'platform',platform)
    if potency_meth: Y.insert(Y.shape[1],'potency_meth',potency_meth)
    #Y.rename(columns=dict(method='ivive_method'),inplace=True)
    
    if dbc_sav:
        dbc_sav.insert_many(Y.to_dict('records'))
    else:
        return Y
