import matplotlib.text as text
import pandas as pd
import numpy as np
import pylab as pl
import scipy as sp
import re

#from IPython.core.debugger import set_trace

import rpy2
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri,FloatVector,ListVector
from rpy2.rinterface import NULL
import scipy.optimize as optz
import seaborn as sns
from scipy import stats
import rpy2.robjects.packages as rpackages
import seaborn as sns
from statsmodels import robust
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl
from functools import reduce

rtcpl = rpackages.importr('tcpl')

def fix_type(x):
    t = type(x)
    if t in [rpy2.rinterface.NARealType,
             rpy2.rinterface.NAIntegerType
            ]:
        return None
    else:
        return x
             
def as_dict(vector):
    """Convert an RPy2 ListVector to a Python dict"""
    result = {}
    for i, name in enumerate(vector.names):
        if isinstance(vector[i], ListVector) :
            result[name] = as_dict(vector[i])
        elif type(vector[i])!=rpy2.rinterface.NULLType:
            if len(vector[i]) == 1:
                x = vector[i][0]
                if not np.isnan(x):
                    result[name]=x
            else:
                result[name] = list(vector[i])
    return result

def calc_TcplFits(X,FT=None,SB0=None,
                  chem_name='name',lconc='lconc',
                  model_exclude='cnst',
                  conc_limit=True,
                  direction={'_up':1,'_dn':-1},
                  info=None):
    
    def _p05(X):  return np.percentile(X,5)
    def _p95(X): return np.percentile(X,95)
    def _p50(X): return np.percentile(X,50)
    def _n(X): return len(X)

    if not FT: 
        raise ValueError("Provide a list of feature(response) columns to calculate fits")
    if not len(X.columns.intersection(FT)) == len(FT):
        raise ValueError("The response data do not contain all features")
        
    if not SB0:
        raise ValueError("Provide a dict of features(response):1 bmad or std of background  values")
    # Test the SB0 to ensure all values are present
    for ft in FT:
        if not SB0.get(ft)>0:
            raise ValueError("SB0 does not have value for %s" % ft)

        
    HITS=[]
    RES=[]
    GA_ci=[]

    for (chem,timeh),X_ct in X.groupby([chem_name,'timeh']):
        LC = X_ct[lconc]
        ci,cf=LC.min(),LC.max()
        Res = []

        for ft in FT:
            Y = X_ct[ft]
            x0= SB0[ft]
            for d,c in direction.items():
                Hit=tcplFit(LC,Y*c,x0*c)
                # Is it resp_max > 3*bmad?
                CRinfo=Hit['cr_info']
                CRinfo['resp_max_gt_3bmad']=CRinfo['resp_max']>3*CRinfo['bmad']
                CR = Hit['best_fit']
                CR.update(dict(ft=ft,ft_ch=ft+d,name=chem,timeh=timeh,ch=c,bmad=x0*c))
                Res.append(CR)
                #Hit.pop('fits')
                Hit.update(dict(ft=ft,ft_ch=ft+d,name=chem,timeh=timeh,ch=c))
                if info: 
                    Hit.update(info)
                    CR.update(info)
                    
                HITS.append(Hit)
        q = ""
        if model_exclude: q+="model!='%s'" % model_exclude
        if conc_limit: 
            if len(q)>0: q+= " and " 
            q+= "ga>=%f and ga<=%f" % (ci,cf)
            
        GA = pd.DataFrame(Res).query(q)\
                    .groupby(['name','timeh'])\
                    .aggregate(dict(ga=[_p05,_p50,_p95,_n]))
        
        GA.columns = [''.join(i) for i in list(GA.columns)]
        GA['n_ft'] = len(FT)
            
        GA_ci.append(GA)
        RES.append(pd.DataFrame(Res).query(q))
    
    return HITS,pd.concat(RES),pd.concat(GA_ci)
                
def tcplFit(conc,resp,bmad=False):
    kwargs={'logc':FloatVector(conc),
            'resp':FloatVector(resp),
           }
    if bmad:
        kwargs['bmad']=1.0*bmad
    else:
        kwargs['force.fit']=True
        
    Y0 = rtcpl.tcplFit(**kwargs)
    Y1= pd.Series(as_dict(Y0))
    F0 = []
    for m in ['cnst','hill','gnls']:
        fit=Y1.pop(m)
        K = [i for i in Y1.index if i.startswith(m)]
        Y2= Y1[K]
        Y2.index=[re.sub('%s_?'%m,'',i) for i in K]
        Y2['model']=m
        if fit==fit:
            F0.append(Y2.to_dict())
        Y1= Y1.drop(K)
        
    Y1['bmad']=bmad
    R0 = {}
    R0['fits']=F0
    #set_trace()
    R0['cr_info']=Y1.to_dict()
    R0['cr_data']=dict(conc=list(conc),resp=list(resp))

    # FIgure out best fit
    F1 = pd.DataFrame(F0)
    R0['best_fit']=F1.sort_values('aic').iloc[0].to_dict()
    
    return R0

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
    
def tcplPlot(Tcpl,ax=None,show_data=True,fnsz=8,r_max=None,
             show_bmds=dict(E=[10,20,50],Z=[1,2,3]),show_legend=True,
             cols=dict(data='green',fit='seagreen',ac50='red',bmad='skyblue')
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
        ax.scatter(cr_c,cr_r,marker='+',s=30,c=cols['data'],alpha=0.7)

    # Plot model results
    ci,cf = cr_c.min()-1,cr_c.max()+1
    C = np.linspace(ci,cf,50)
    R = model(C,**kw)
    ax.plot(C,R,color=cols['fit'],label='fit',alpha=0.7)

    if not r_max:
        r_max = np.max([max(R),cr_r.max()])

    # SHow the bmads
    ax.hlines(b0*np.array([-2,-2,-1,1,2,3]),ci,cf,
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
            ax.vlines(bmd['bmd'],0,r_max,lw=2,
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
        tick.set_fontsize(8)

    if show_legend: ax.legend(fontsize='xx-small')
    
    ax.set_ylim(-0.2,r_max*1.1)    
    ax.set_title("%(name)s %(ft_ch)s %(timeh)ih" % Tcpl,   
                 fontdict=dict(size=fnsz,family='monospace'))
    
    
def plotFits(chem,timeh=None,assay=None,assay_comp=None,dbc=None,
             nc=5,r_max=None,
             Models=['hill','gnls','cnst'],
             **kwargs):
    if not dbc:
        raise ValueError
    Q=dict(name=chem)
    if timeh: 
        if type(timeh)==list:
            Q['timeh']={'$in':timeh}
        else:
            Q['timeh']=timeh
    if assay: Q['ft']=assay
    if assay_comp: 
        if type(assay_comp)==list:
            Q['ft_ch']={'$in':assay_comp}
        else:
            Q['ft_ch']=assay_comp
    Q['best_fit.model']={'$in':Models}
    
    N = dbc.find(Q).count()
    if N==0: 
        raise LookupError
     
    nc = nc if N>nc else N
    nr = 1+int(1.0*N/nc)
    fig = pl.figure(figsize=(nc*3,nr*2))
    
    Fits = [i for i in dbc.find(Q)]
    Fits_df = pd.DataFrame([i['best_fit'] for i in Fits])
    # Find out the max top
    if not r_max:
        r_max = Fits_df.tp.max() * 1.2
    
    i=0

    for TcplFit in Fits:
        i+=1
        ax = pl.subplot(nr,nc,i)
        tcplPlot(TcplFit,ax=ax,r_max=r_max,**kwargs)
    
    pl.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9,
                wspace=0.5, hspace=0.5)
    
    return fig

def getAssayCompOrder():
    Asy=['MF', 'St', 'ES', 'LM', 'DT', 'NS', 'Ap', 'CN', 'X']
    AC =reduce(lambda a,b: a+b, [[i+j for j in ['_up','_dn']] for i in Asy])
    return dict(list(zip(AC,list(range(len(AC))))))
    
def plotChemFits(chem,figfile=None,r_max=None,ft_type='single',
                 dbc=None,dbg=False,nc=5,
                 Models=['hill','gnls','cnst'],
                 **kwargs):
    if not dbc:
        raise ValueError
    
    Q=dict(name=chem)
    if ft_type: Q['ft_type']=ft_type
    #Q['best_fit.model']={'$in':['hill','gnls']}
    ASCMP = dbc.find(Q).distinct('ft_ch')
    #if len(ASCMP)<=1:
    #    ASCMP += ['MF_dn','St_up','ES_up','LM_dn','Ap_up','CN_dn','X_up']
    #    ASCMP = list(set(ASCMP))
        
    ACO = getAssayCompOrder() 
    ASCMP = sorted(ASCMP,cmp=lambda a,b: 1 if ACO[a]>ACO[b] else -1)
    
    if dbg: print("FT",len(ASCMP))
    #pdf = PdfPages(figfile)
    
    nr = int(np.ceil(1.0*len(ASCMP)/nc))
    
    pl.close('all')
    fig,AX=pl.subplots(nr,nc,sharey='row',figsize=(nc*3,nr*2))
    
    colTh = mpl.cm.ScalarMappable(mpl.colors.Normalize(vmin=0,vmax=72),
                                  cmap=mpl.cm.Reds) 

    for i,AC in enumerate(ASCMP):
        if dbg: print(i,AC)
        Q['ft_ch']=AC
        if not r_max: 
            r_max=getRmax(Q,dbc)
        j = int(np.floor(1.0*i/nc))
        ax = AX[j,i % nc]
        for TcplFit in dbc.find(Q).sort('timeh'):
            if dbg: print(TcplFit['name'])
            tcol = colTh.to_rgba(TcplFit['timeh'])    
            tcplPlot(TcplFit,ax=ax,r_max=r_max,
                     cols=dict(data=tcol,fit=tcol,bmad='grey',bmd='red'),
                     **kwargs)
    
        
    pl.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9,
                wspace=0.5, hspace=0.5)
    
    return fig


def getAssayCompOrder():
    Asy=['MF', 'St', 'ES', 'LM', 'DT', 'NS', 'Ap', 'CN', 'X']
    AC =reduce(lambda a,b: a+b, [[i+j for j in ['_up','_dn']] for i in Asy])
    return dict(list(zip(AC,list(range(len(AC))))))
    
def plotChemFits(chem,figfile=None,r_max=None,ft_type='single',
                 dbc=None,dbg=False,nc=5,
                 Models=['hill','gnls','cnst'],
                 **kwargs):
    if not dbc:
        raise ValueError
    
    Q=dict(name=chem)
    if ft_type: Q['ft_type']=ft_type
    #Q['best_fit.model']={'$in':['hill','gnls']}
    ASCMP = dbc.find(Q).distinct('ft_ch')
    #if len(ASCMP)<=1:
    #    ASCMP += ['MF_dn','St_up','ES_up','LM_dn','Ap_up','CN_dn','X_up']
    #    ASCMP = list(set(ASCMP))
        
    ACO = getAssayCompOrder() 
    ASCMP = sorted(ASCMP,cmp=lambda a,b: 1 if ACO[a]>ACO[b] else -1)
    
    if dbg: print("FT",len(ASCMP))
    #pdf = PdfPages(figfile)
    
    nr = int(np.ceil(1.0*len(ASCMP)/nc))
    
    pl.close('all')
    fig,AX=pl.subplots(nr,nc,sharey='row',figsize=(nc*3,nr*2))
    
    colTh = mpl.cm.ScalarMappable(mpl.colors.Normalize(vmin=0,vmax=72),
                                  cmap=mpl.cm.Reds) 

    for i,AC in enumerate(ASCMP):
        if dbg: print(i,AC)
        Q['ft_ch']=AC
        if not r_max: 
            r_max=getRmax(Q,dbc)
        j = int(np.floor(1.0*i/nc))
        ax = AX[j,i % nc]
        for TcplFit in dbc.find(Q).sort('timeh'):
            if dbg: print(TcplFit['name'])
            tcol = colTh.to_rgba(TcplFit['timeh'])    
            tcplPlot(TcplFit,ax=ax,r_max=r_max,
                     cols=dict(data=tcol,fit=tcol,bmad='grey',bmd='red'),
                     **kwargs)
    
        
    pl.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9,
                wspace=0.5, hspace=0.5)
    
    return fig

def plotChemFits0(chem,figfile=None,r_max=None,ft_type='single',
                 dbc=None,dbg=False,
                 Models=['hill','gnls','cnst'],
                 **kwargs):
    if not dbc:
        raise ValueError
    
    Q=dict(name=chem)
    if ft_type: Q['ft_type']=ft_type
    #Q['best_fit.model']={'$in':['hill','gnls']}
    ACO = getAssayCompOrder() 
    ASCMP = dbc.find(Q).distinct('ft_ch')
    if len(ASCMP)<=1:
        ASCMP += ['MF_dn','St_up','ES_up','LM_dn','Ap_up','CN_dn','X_up']
        ASCMP = list(set(ASCMP))
        
    ASCMP = sorted(ASCMP,cmp=lambda a,b: 1 if ACO[a]>ACO[b] else -1)
    
    if dbg: print("FT",len(ASCMP))
    #pdf = PdfPages(figfile)
    
    nr,nc=len(ASCMP),6
    pl.close('all')
    fig,AX=pl.subplots(nr,nc,sharey='row',figsize=(nc*3,nr*2))
    
    for i,AC in enumerate(ASCMP):
        if dbg: print(i,AC)
        Q['ft_ch']=AC
        r_max=getRmax(Q,dbc)
        for j,TcplFit in enumerate(dbc.find(Q).sort('timeh')):
            if dbg: print(TcplFit['name'],i,j)
            ax = AX[i,j]
            tcplPlot(TcplFit,ax=ax,r_max=r_max,**kwargs)
    
    
    pl.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9,
                wspace=0.5, hspace=0.5)
    
    return fig

def getRmax(Qi,dbc=None):
    r_max = max([i['best_fit']['tp'] for i in 
                   dbc.find(Qi,{'_id':0,'best_fit.tp':1})
                 if 'tp' in i['best_fit']
                ] or [1])
    r_max = 1.2 if r_max!=r_max else r_max * 1.2

    return r_max
