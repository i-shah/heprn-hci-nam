
import matplotlib.text as text
import pandas as pd
import numpy as np
import pylab as pl
import scipy as sp
import re

import rpy2
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri,FloatVector,ListVector
#from rpy2.rinterface import RNULLType
import scipy.optimize as optz
import seaborn as sns
from scipy import stats
import rpy2.robjects.packages as rpackages
import seaborn as sns
from statsmodels import robust
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl
import sys

from box import Box

def hillF(x,tp,ga,gw):
    return tp/(1+10**((ga-x)*gw))

def gnlsF(x,tp,ga,gw,lw,la):
    gn = 1/(1+10**((ga-x)*gw))
    ls = 1/(1+10**((x-la)*lw))
    return tp*gn*ls

def calc_Resp(Tcpl):
    BF = Tcpl['best_fit']
    b0 = Tcpl['cr_info']['bmad']
    ch = BF['ch']
    
    model=None
    if BF['model']=='hill':
        kw = {o:v for o,v in iter(BF.items()) if o in ['tp','ga','gw'] }
        model=hillF
    elif BF['model']=='gnls':
        kw = {o:v for o,v in iter(BF.items()) if o in ['tp','ga','gw','la','lw'] }
        model=gnlsF
    elif BF['model']=='cnst':
        kw = {}
        def model(x): 
            y=np.median(Tcpl['cr_data']['resp'])
            return [y]*len(x)
    
    # CR Data
    cr_c = np.array(Tcpl['cr_data']['conc'])
    cr_r = np.array(Tcpl['cr_data']['resp'])    
    c0 = cr_c.min()-0.8

    # Response
    ci,cf = cr_c.min()-1,cr_c.max()+1
    C = np.linspace(ci,cf,50)
    R = model(C,**kw)

    return cr_c,cr_r,model,kw

def calc_BMDs(Tcpl,BMR=dict(E=[10,20,30],Z=[1,2,3]),
              ret='dict',add_info=False,
              dbg=False):
    """
    Calculate benchmark doses corresponding to bmrs:-
    E: fractional efficacy (top) 
    Z: number of standard deviations (assumes response is in units of Z)
    """
    
    BF = Tcpl['best_fit']
    if BF['model']=='cnst':
        return pd.DataFrame()
        
    C,R,model,kw = calc_Resp(Tcpl)
    ci,cf = C.min(),C.max()

    BMD=[]
    if 'E' in BMR:
        for e0 in BMR.get('E'):
            def bmdf(c): return model(c,**kw)-e0*kw['tp']/100.0
            try:
                soln = optz.fsolve(bmdf,[ci,cf])
                soln = soln[np.logical_and(soln>ci,soln<cf)]
                bmd0 = np.min(soln)
            except:
                bmd0 = None
                if dbg: print("Failed E %0.2f" % e0)
            else:
                BMD.append(dict(bmr_type='E',bmr=e0,bmd=bmd0,bmd_um=10**bmd0))
    
    if 'Z' in BMR:
        for z0 in BMR.get('Z'):
            def bmdf(c): return model(c,**kw)-z0
            try:
                soln = optz.fsolve(bmdf,[ci,cf])
                soln = soln[np.logical_and(soln>ci,soln<cf)]
                bmd0 = np.min(soln)
            except:
                bmd0 = None
                if dbg: print("Failed Z %0.2f" % z0)
            else:
                BMD.append(dict(bmr_type='Z',bmr=z0,bmd=bmd0,bmd_um=10**bmd0))
    
    if len(BMD)==0:
        return pd.DataFrame()
    
    if add_info: ret = 'df'
    if ret=='df':
        
        DF = pd.DataFrame(BMD)
        if add_info:
            for i in ['name','timeh','ft','ft_ch','ch','ft_type']:
                DF.loc[:,i]=Tcpl[i]
        return DF
    else:
        return BMD
    
def plotTcpl(Tcpl,ax=None,show_data=True,fnsz=8,r_max=None,
             show_bmds=dict(E=[10,20,50],Z=[1,2,3]),show_legend=True,
             cols=dict(data='green',fit='seagreen',ac50='red',bmad='skyblue'),
             title_fmt="%(name)s %(ft)s %(timeh)ih"
            ):
    BF = Tcpl['best_fit']
    b0 = Tcpl['cr_info']['bmad']
    ch = BF['ch']
    
    if not ax:
        ax = pl.subplot(111)
    
    # CR Data
    cr_c,cr_r,model,kw = calc_Resp(Tcpl)

    c0 = cr_c.min()-0.8

    if show_data:
        ax.scatter(cr_c,cr_r*ch,marker='+',s=30,c=cols['data'],alpha=0.7)

    # Plot model results
    ci,cf = cr_c.min()-1,cr_c.max()+1
    C = np.linspace(ci,cf,50)
    R = ch * model(C,**kw)
    ax.plot(C,R,color=cols['fit'],label='fit',alpha=0.7)

    if not r_max:
        r_max = ch * np.max([max(R),cr_r.max(),5])
    else:
        r_max = ch * r_max
        
    # SHow the bmads
    ax.hlines(b0*np.array([-2,-2,-1,0,1,2,3]),ci,cf,
              lw=0.3,color=cols['bmad'],linestyle='-.')
    
    # Show AC50
    #if BF.has_key('ga'):
    #    ax.vlines(BF['ga'],0,r_max,color='red',lw=1)
    #    ax.text(c0,r_max*0.9,r"AC50=%5.2f$\mu$M" % 10**BF['ga'],
    #           fontdict=dict(size=fnsz,color=cols['ac50'],alpha=0.6))

    # Find the BMD10
    if show_bmds and BF['model']!='cnst':
        colE = mpl.cm.ScalarMappable(mpl.colors.Normalize(vmin=0,vmax=100),
                                 cmap=mpl.cm.Reds) 
        colZ = mpl.cm.ScalarMappable(mpl.colors.Normalize(vmin=0,vmax=5),
                                 cmap=mpl.cm.Blues) 
        
        BMDS = calc_BMDs(Tcpl,BMR=show_bmds)
        for bmd in BMDS:
            colF = colE if bmd['bmr_type']=='E' else colZ
            ax.vlines(bmd['bmd'],-1*ch,r_max*1.5,lw=2,
                      color=colF.to_rgba(bmd['bmr']), 
                      label=r'$%(bmr_type)s_{%(bmr)i}$ %(bmd_um).1f$\mu$M' % bmd)
            #ax.text(c0,r_max*0.7,r"BMD10=%5.2f$\mu$M" % 10**bmd10,
            #        fontdict=dict(size=fnsz,color=cols['bmd']))

    ax.set_xlim(ci,cf)
    
    #Conc=np.linspace(cr_c.min()-0.5,cr_c.max()+0.5,9)
    Conc=np.array([0.1,0.5,1,5,10,25,50,100,250])
    Clab=['0.1','0.5','1','5','10','25','50','100','250']
    #Clab=['%5.2f'%c for c in Conc]
    
    ax.set_xticks(np.log10(Conc), minor=False)
    xl=ax.set_xticklabels(Clab,rotation=90)
    #for tick in ax.get_xticklines(): 
    #    tick.set_visible(False)
    for tick in ax.get_xticklabels(): 
        tick.set_fontsize(6)

    if show_legend: ax.legend(fontsize='xx-small')

    if ch>0:
        ax.set_ylim(-0.2*r_max,r_max*1.5)    
    else:
        ax.set_ylim(r_max*1.5,-0.2*r_max)    
        
    ax.set_title(title_fmt % Tcpl,   
                 fontdict=dict(size=fnsz,family='monospace'))
    


def plotAllFits(chem,timeh=[6,24,48,72],
                dbc_cr=None,
                dbc_qual=None,
                title_fnsz=24,
                r_max=None,dbg=False,xs=3,ys=2,
                Models=['hill','gnls','cnst'],
                show_bmds=dict(E=50),
                **kwargs):
    if not dbc_cr:
        raise ValueError
    
    if dbc_cr.count_documents(dict(name=chem))==0: 
        raise LookupError
            
    
    
    Assays = [ 'ES',  'MF', 'LM','St','DT', 'NS', 'Ap', 'CN']
    nr,nc = len(timeh)+1,len(Assays)
    fig = pl.figure(figsize=(nc*xs,nr*ys))
    Q={}
    
    i=nc
    ax = pl.subplot(nr,nc,1)
    ax.text(5,5,chem,fontsize=title_fnsz)
    ax.set_xlim(0,60)
    ax.set_ylim(0,10)
    ax.set_axis_off()
    HQ=None
    for th in timeh:
        for ft in Assays:
            Q=dict(name=chem,timeh=th,ft=ft)
            Q['best_fit.model']={'$in':['hill','gnls']}
            Fit = None
            nh = dbc_cr.count_documents(Q)

            if dbg: print(("%(name)s %(ft)s %(timeh)i " % Q + str(nh)))
            if nh==1:
                Fit = dbc_cr.find_one(Q)
            elif nh==2:
                Q['best_fit.model']='hill'
                Fit = dbc_cr.find_one(Q)                
            elif nh==0:
                Q.pop('best_fit.model')
                Q['ch']=1
                Fit = dbc_cr.find_one(Q)
            if Fit['best_fit']['model']!='cnst':
                HQ = dbc_qual.find_one(dict(name=chem,timeh=th,ft=ft))
                
            i+=1
            ax = pl.subplot(nr,nc,i)
            plotTcpl(Fit,ax=ax,r_max=r_max,
                     show_bmds=show_bmds if HQ and HQ['hit']==1 else None,
                     **kwargs)
    
    pl.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9,
                wspace=0.5, hspace=0.5)
    
    return fig    