import pandas as pd
import numpy as np
import scipy as sp
import pylab as pl
from collections import *
from copy import copy
import os
import subprocess
import numpy.linalg as LA
from scipy.stats import iqr

from statsmodels import robust 
import seaborn as sns
from matplotlib import cm


def calcWellChanges(plate_id,dbc_ch=None,dbc_raw=None,
                   ctrl_neg='ctrl-',ctrl_pos='ctrl+',
                   test_chem='chem_test',ft_name='FA0',
                   normal='all'):
    """
    Use plate level data to generate L2FC, Pct values at a replicate level
    
    KWArgs:
    normal = median plate median / 
             ctrl-  negative control
             
    """
    C1 = ['dsstox_sid','sample_id','timeh','plate_id','stype',
      'conc','FA0','FA1','FN0','FN1','raw_value','row','col',
      'l2fc_dmso','l2fc_plt','z1_dmso','z1_plt'
     ]
    def stp5(X): np.percentile(X,5)
    def stp95(X): np.percentile(X,95)
    def p02(X): return np.percentile(X,2)
    def p98(X): return np.percentile(X,98)
    def p50(X): return np.percentile(X,50)
    def cv(X): return np.std(X)/np.mean(X)
    def outHi(X): return np.percentile(X,75)+1.5*iqr(X)
    def outLo(X): return np.percentile(X,25)-1.5*iqr(X)

    Plate = dbc_raw.find_one(dict(plate_id=plate_id),dict(_id=0))
    Wells=pd.DataFrame(Plate.pop('wells'))
    Wells.loc[:,'lconc']=np.log10(Wells.conc)
    P_neg=Wells[Wells.stype==ctrl_neg]
    P_pos =Wells[Wells.stype==ctrl_pos]
    P_chem=Wells[Wells.stype.isin([test_chem,ctrl_pos])]

    ST_neg=P_neg.groupby([ft_name])\
                .aggregate(dict(raw_value=[np.mean,np.median,np.std,robust.mad,p02,p98,p50]))
    ST_neg.columns=['ctrl_n_'+i[1] for i in ST_neg.columns]
    ST_neg = ST_neg.reset_index()
    
    ST_pos=P_pos.groupby([ft_name])\
                .aggregate(dict(raw_value=[np.mean,np.median,np.std,robust.mad,p02,p98,p50]))
    ST_pos.columns=['ctrl_p_'+i[1] for i in ST_pos.columns]
    ST_pos = ST_pos.reset_index()

    ST_plt =Wells.groupby([ft_name])\
                 .aggregate(dict(raw_value=[np.mean,np.median,np.std,robust.mad,p02,p98,p50]))
    ST_plt.columns=['plt_'+i[1] for i in ST_plt.columns]
    ST_plt = ST_plt.reset_index()

    CH = Wells[['dsstox_sid','sample_id','name','timeh',
                        'stype',ft_name,'lconc','conc','raw_value']]
    CH=CH.merge(ST_neg[[ft_name,'ctrl_n_mean','ctrl_n_std',
                        'ctrl_n_mad']],
                on=ft_name)\
            .merge(ST_pos[[ft_name,'ctrl_p_mean','ctrl_p_std',
                                'ctrl_p_mad']],
                        on=ft_name)\
            .merge(ST_plt[[ft_name,'plt_mean','plt_std','plt_mad']],
                       on=ft_name)

    CH.loc[:,'l2fc_ctrl']=np.log2(CH.raw_value/CH.ctrl_n_mean)
    CH.loc[:,'z_ctrl']   =(CH.raw_value- CH.ctrl_n_mean)/CH.ctrl_n_std
    CH.loc[:,'l2fc_plt'] =np.log2(CH.raw_value/CH.plt_mean)
    CH.loc[:,'z_plt']   =(CH.raw_value-CH.plt_mean)/CH.plt_std
    CH.loc[:,'pv_n_ctrl']  =sp.stats.norm.cdf(CH.raw_value,CH.ctrl_n_mean,CH.ctrl_n_std)
    CH.loc[:,'pv_plt']   =sp.stats.norm.cdf(CH.raw_value,CH.plt_mean,CH.plt_std)
    
    CH['pct_ctrl']= (CH.raw_value-CH.ctrl_n_mean)/(CH.ctrl_p_mean-CH.ctrl_n_mean)
    CH['pct_plt'] = (CH.raw_value-CH.plt_mean)/CH.plt_mean
    CH['pct_ctrl1']= CH.raw_value/CH.ctrl_n_mean
    CH['pct_plt1'] = CH.raw_value/CH.plt_mean

    
    Plate['stats_ctrl_n']  = ST_neg.reset_index().to_dict('records')
    Plate['stats_ctrl_p']  = ST_pos.reset_index().to_dict('records')
    Plate['stats_plt']  = ST_plt.reset_index().to_dict('records')    
    Plate['chem_ch'] = CH.reset_index().to_dict('records')
    
    if dbc_ch:
        dbc_ch.insert_one(Plate)
    else:
        return Plate


def calcChanges(plate_id,dbc_ch=None,dbc_raw=None,
               ctrl_neg='ctrl-',ctrl_pos='ctrl+',
               test_chem='chem_test',ft_name='FA0',
               normal='all'):
    """
    Use plate level data to generate L2FC, Pct values
    
    KWArgs:
    normal = median plate median / 
             ctrl-  negative control
             
    """
    C1 = ['dsstox_sid','sample_id','timeh','plate_id','stype',
      'conc','FA0','FA1','FN0','FN1','raw_value','row','col',
      'l2fc_dmso','l2fc_plt','z1_dmso','z1_plt'
     ]
    def stp5(X): np.percentile(X,5)
    def stp95(X): np.percentile(X,95)
    def p02(X): return np.percentile(X,2)
    def p98(X): return np.percentile(X,98)
    def p50(X): return np.percentile(X,50)
    def cv(X): return np.std(X)/np.mean(X)
    def outHi(X): return np.percentile(X,75)+1.5*iqr(X)
    def outLo(X): return np.percentile(X,25)-1.5*iqr(X)

    Plate = dbc_raw.find_one(dict(plate_id=plate_id),dict(_id=0))
    Wells=pd.DataFrame(Plate.pop('wells'))
    Wells.loc[:,'lconc']=np.log10(Wells.conc)
    P_neg=Wells[Wells.stype==ctrl_neg]
    P_pos =Wells[Wells.stype==ctrl_pos]
    P_chem=Wells[Wells.stype.isin([test_chem,ctrl_pos])]

    ST_neg=P_neg.groupby([ft_name])\
                .aggregate(dict(raw_value=[np.mean,np.median,np.std,robust.mad,p02,p98,p50]))
    ST_neg.columns=['ctrl_n_'+i[1] for i in ST_neg.columns]
    ST_neg = ST_neg.reset_index()
    
    ST_pos=P_pos.groupby([ft_name])\
                .aggregate(dict(raw_value=[np.mean,np.median,np.std,robust.mad,p02,p98,p50]))
    ST_pos.columns=['ctrl_p_'+i[1] for i in ST_pos.columns]
    ST_pos = ST_pos.reset_index()

    ST_plt =Wells.groupby([ft_name])\
                 .aggregate(dict(raw_value=[np.mean,np.median,np.std,robust.mad,p02,p98,p50]))
    ST_plt.columns=['plt_'+i[1] for i in ST_plt.columns]
    ST_plt = ST_plt.reset_index()

    CH = Wells.groupby(['dsstox_sid','sample_id','name','timeh',
                        'stype',ft_name,'lconc','conc'])\
            .aggregate(dict(raw_value=[np.mean,cv]))
    CH.columns=['raw_'+i[1] for i in CH.columns]
    CH = CH.reset_index()
    CH=CH.merge(ST_neg[[ft_name,'ctrl_n_mean','ctrl_n_std',
                        'ctrl_n_mad']],
                on=ft_name)\
    .merge(ST_pos[[ft_name,'ctrl_p_mean','ctrl_p_std',
                        'ctrl_p_mad']],
                on=ft_name)\
    .merge(ST_plt[[ft_name,'plt_mean','plt_std','plt_mad']],
               on=ft_name)

    CH.loc[:,'l2fc_ctrl']=np.log2(CH.raw_mean/CH.ctrl_n_mean)
    CH.loc[:,'z_ctrl']   =(CH.raw_mean- CH.ctrl_n_mean)/CH.ctrl_n_std
    CH.loc[:,'l2fc_plt'] =np.log2(CH.raw_mean/CH.plt_mean)
    CH.loc[:,'z_plt']   =(CH.raw_mean-CH.plt_mean)/CH.plt_std
    CH.loc[:,'pv_n_ctrl']  =sp.stats.norm.cdf(CH.raw_mean,CH.ctrl_n_mean,CH.ctrl_n_std)
    CH.loc[:,'pv_plt']   =sp.stats.norm.cdf(CH.raw_mean,CH.plt_mean,CH.plt_std)
    
    CH['pct_ctrl']= (CH.raw_mean-CH.ctrl_n_mean)/(CH.ctrl_p_mean-CH.ctrl_n_mean)
    CH['pct_plt'] = (CH.raw_mean-CH.plt_mean)/CH.plt_mean
    CH['pct_ctrl1']= CH.raw_mean/CH.ctrl_n_mean
    CH['pct_plt1'] = CH.raw_mean/CH.plt_mean

    
    Plate['stats_ctrl_n']  = ST_neg.reset_index().to_dict('records')
    Plate['stats_ctrl_p']  = ST_pos.reset_index().to_dict('records')
    Plate['stats_plt']  = ST_plt.reset_index().to_dict('records')    
    Plate['chem_ch'] = CH.reset_index().to_dict('records')
    
    if dbc_ch:
        dbc_ch.insert_one(Plate)
    else:
        return Plate

def calcFC(plate_id,dbc_fc=None,dbc_raw=None,
           ctrl_neg='ctrl-',ctrl_pos='ctrl+',
           test_chem='chem_test',ft_name='FA0',
           normal='all'):
    """
    Use plate level data to generate L2FC value
    
    KWArgs:
    normal = median plate median / 
             ctrl-  negative control
             
    """
    C1 = ['dsstox_sid','sample_id','timeh','plate_id','stype',
      'conc','FA0','FA1','FN0','FN1','raw_value','row','col',
      'l2fc_dmso','l2fc_plt','z1_dmso','z1_plt'
     ]
    def stp5(X): np.percentile(X,5)
    def stp95(X): np.percentile(X,95)
    def cv(X): return np.std(X)/np.mean(X)
    

    Plate = dbc_raw.find_one(dict(plate_id=plate_id),dict(_id=0))
    Wells=pd.DataFrame(Plate.pop('wells'))
    Wells.loc[:,'lconc']=np.log10(Wells.conc)
    P_neg=Wells[Wells.stype==ctrl_neg]
    P_pos =Wells[Wells.stype==ctrl_pos]
    P_chem=Wells[Wells.stype==test_chem]

    ST_neg=P_neg.groupby([ft_name])\
                .aggregate(dict(raw_value=[np.mean,np.median,np.std,robust.mad]))
    ST_neg.columns=['ctrl_'+i[1] for i in ST_neg.columns]
    ST_neg = ST_neg.reset_index()
    
    ST_plt =Wells.groupby([ft_name])\
                 .aggregate(dict(raw_value=[np.mean,np.median,np.std,robust.mad]))
    ST_plt.columns=['plt_'+i[1] for i in ST_plt.columns]
    ST_plt = ST_plt.reset_index()

    FC = Wells.groupby(['dsstox_sid','sample_id','name','timeh',
                        'stype',ft_name,'lconc','conc'])\
            .aggregate(dict(raw_value=[np.mean,cv]))
    FC.columns=['raw_'+i[1] for i in FC.columns]
    FC = FC.reset_index()
    FC=FC.merge(ST_neg[[ft_name,'ctrl_mean','ctrl_std','ctrl_mad']],
                left_on=ft_name,right_on=ft_name)\
        .merge(ST_plt[[ft_name,'plt_mean','plt_std','plt_mad']],
               left_on=ft_name,right_on=ft_name)

    FC['l2fc_ctrl']=np.log2(FC.raw_mean/FC.ctrl_mean)
    FC['z_ctrl']   =(FC.raw_mean-FC.ctrl_mean)/FC.ctrl_std
    FC['l2fc_plt'] =np.log2(FC.raw_mean/FC.plt_mean)
    FC['z_plt']   =(FC.raw_mean-FC.plt_mean)/FC.plt_std
    FC['pv_ctrl']  =sp.stats.norm.cdf(FC.raw_mean,FC.ctrl_mean,FC.ctrl_std)
    FC['pv_plt']   =sp.stats.norm.cdf(FC.raw_mean,FC.plt_mean,FC.plt_std)
    
    
    Plate['dmso_stats']  = ST_neg.reset_index().to_dict('records')
    Plate['plate_stats']  = ST_plt.reset_index().to_dict('records')    
    Plate['chem_fc'] = FC.reset_index().to_dict('records')
    
    if dbc_fc:
        dbc_fc.insert_one(Plate)
    else:
        return Plate
    
def calcFCSm(plate_id,dbc_fc=None,dbc_sm=None,
           ctrl_neg='ctrl-',ctrl_pos='ctrl+',
           test_chem='chem_test',ft_name='FA0',
             smooth='smooth_w7',
           normal='all'):
    """
    Use plate level data to generate L2FC value
    
    KWArgs:
    normal = median plate median / 
             ctrl-  negative control
             
    """
    C1 = ['dsstox_sid','sample_id','timeh','plate_id','stype',
      'conc','FA0','FA1','FN0','FN1','raw_value','row','col',
      'l2fc_dmso','l2fc_plt','z1_dmso','z1_plt'
     ]
    def stp5(X): np.percentile(X,5)
    def stp95(X): np.percentile(X,95)
    def cv(X): return np.std(X)/np.mean(X)
    

    Plate = dbc_sm.find_one(dict(plate_id=plate_id),dict(_id=0))
    Wells=pd.DataFrame(Plate.pop('wells'))
    Smooth=pd.DataFrame(Plate.get('smooth'))

    Wells.loc[:,'lconc']=np.log10(Wells.conc)
    P_neg=Wells[Wells.stype==ctrl_neg]
    P_pos =Wells[Wells.stype==ctrl_pos]
    P_chem=Wells[Wells.stype==test_chem]

    ST_neg=P_neg.groupby([ft_name])\
                .aggregate(dict(raw_value=[np.mean,np.median,np.std,robust.mad]))
    ST_neg.columns=['ctrl_'+i[1] for i in ST_neg.columns]
    ST_neg = ST_neg.reset_index()
    
    ST_plt =Smooth.groupby([ft_name])\
                     .aggregate({smooth:[np.mean,np.median,
                                         np.std,robust.mad]})
    
    ST_plt.columns=['plt_'+i[1] for i in ST_plt.columns]
    ST_plt = ST_plt.reset_index()

    FC=Smooth.merge(ST_neg[[ft_name,'ctrl_mean','ctrl_std',
                          'ctrl_mad']],
                left_on=ft_name,right_on=ft_name)\
        .merge(ST_plt[[ft_name,'plt_mean','plt_std','plt_mad']],
               left_on=ft_name,right_on=ft_name)

    FC['l2fc_ctrl']=np.log2(FC[smooth]/FC.ctrl_mean)
    FC['z_ctrl']   =(FC[smooth]-FC.ctrl_mean)/FC.ctrl_std
    FC['l2fc_plt'] =np.log2(FC[smooth]/FC.plt_mean)
    FC['z_plt']   =(FC[smooth]-FC.plt_mean)/FC.plt_std
    FC['pv_ctrl']  =sp.stats.norm.cdf(FC[smooth],FC.ctrl_mean,FC.ctrl_std)
    FC['pv_plt']   =sp.stats.norm.cdf(FC[smooth],FC.plt_mean,FC.plt_std)
    
    
    Plate['dmso_stats']  = ST_neg.reset_index().to_dict('records')
    Plate['plate_stats']  = ST_plt.reset_index().to_dict('records')    
    Plate['chem_fc'] = FC.reset_index().to_dict('records')
    
    if dbc_fc:
        dbc_fc.insert_one(Plate)
    else:
        return Plate
    
def getChemFC(chem,feat='FA0',val='l2fc_ctrl',pivot=False,
              ind=None,fc_key='chem_ch',
              what='traj',add_t0=False,use_times=None,
              dbc_fc=None,dbc_ft=None):
    
    Q = dict()
    
    Features=pd.DataFrame(list(dbc_ft.find({},dict(_id=0))))
    X=pd.concat([pd.DataFrame(i.get(fc_key)) for i in dbc_fc.find(dict(chemicals=chem),{'_id':0})])
    X= X[X.name==chem]
    FC1=X.merge(Features,on='FA0')

    if use_times:
            FC1 = FC1[(FC1.timeh.isin(use_times))]

    I1 = ind if ind else ['name','conc','timeh']

    if what=='traj' and not ind:
        I1 = ['name','conc','timeh']
    else:
        I1 = ['name','timeh','conc']
        
    if not pivot: 
        return FC1
    
    if feat=='FA0':
        FC2=FC1.pivot_table(index=I1,
                           columns='FA0',
                           values=val,fill_value=0)
    elif feat=='FA1':
        FC2=FC1.pivot_table(index=I1,
                           columns='FA1',
                           values=val,fill_value=0,
                           aggfunc=np.mean)
        
    if add_t0 and what=='traj':
        Ind1=list(set([(i[0],i[1]) for i in FC2.index]))
    
        for (ch,co) in Ind1:
            FC2.loc[(ch,co,0)]=0

        FC2.sort_index(inplace=True)

    
    return FC2

def vizChemTraj(chem,pv=False,vm=2,show_concs=None,
                nr=2,nc=5,dx=2,AX=[],title=1,
                **kwds):
    FCp = getChemFC(chem,pivot=True,what='traj',**kwds)
    feat=kwds.get('feat')
    val =kwds.get('val')

    if pv:
        suf= val.split('_')[1]
        pv_va = 'pv_'+ suf if suf=='plt' else 'pv_n_ctrl'
        kw = {k:v for k,v in iter(kwds.items()) if k.startswith('dbc')}
        PV = getChemFC(chem,feat=feat,val=pv_va,
                       pivot=True,**kw)

    if feat=='FA0':
        C1 = ['MF','St','ES','LM','DT','NS1','NS2','Ap','CN1','CN2']
    elif feat=='FA1':
        C1 = ['MF','St','ES','LM','DT','NS','Ap','CN']
        
    if show_concs:
        FCp = FCp.select(lambda I: I[1] in show_concs)
        if pv:
            PV = PV.select(lambda I: I[1] in show_concs)
    
    Ind=[j for j in 
         sorted(set([(i[0],i[1]) for i in list(FCp.index)]),
           key=lambda x: x[1])
         if j[1]!=0]

    
    if len(AX)==0:
        fig,AX = pl.subplots(nrows=nr,ncols=nc,figsize=(nc*dx,nr*dx),
                             sharex=True,sharey=True)
    else:
        fig = None
        
    G = []
    for ir in range(nr):
        for ic in range(nc):
            ii = ir*nc + ic
            ind= Ind[ii]
            X = FCp.loc[ind]
            X = X[C1].loc[X.index[::-1]]
            ax = AX[ir][ic] if nr>1 else AX[ic]
            if pv:
                Y = PV.loc[ind]
                Y = Y[C1].loc[Y.index[::-1]]
                Y1 = pd.DataFrame(np.array([['']*Y.shape[1]]*Y.shape[0]),
                                  columns=Y.columns)
                Y1.index=Y.index
                Y1[Y<pv]='*'
                g = sns.heatmap(X,vmin=-vm,vmax=vm,cmap=cm.RdYlBu_r,
                                annot=Y1,fmt='s',
                                lw=1,ax=ax,cbar=False)
            else:
                g = sns.heatmap(X,vmin=-vm,vmax=vm,cmap=cm.RdYlBu_r,
                                lw=1,ax=ax,cbar=False)
            if ir==nr-1: 
                ax.set_xlabel('')
            else: 
                ax.set_xlabel('')

            if ic==0 : 
                ax.set_ylabel('timeh [h]')
            if title==1: 
                ax.set_title(r'%s %6.2f$\mu$M' % (chem,ind[1]))
            elif title==2:
                ax.set_title(r'%6.2f$\mu$M' % ind[1])
     
    if fig:
        pl.subplots_adjust(wspace=0.3,hspace=0.3,top=0.9,bottom=0.1)
    
    return fig


def vizChemCR(chem,pv=False,vm=2,**kwds):
    FCp = getChemFC(chem,pivot=True,what='cr',**kwds)
    feat=kwds.get('feat')
    val =kwds.get('val')

    if pv:
        suf= val.split('_')[1]
        pv_va = 'pv_'+ suf if suf=='plt' else 'pv_n_ctrl'
        kw = {k:v for k,v in iter(kwds.items()) if k.startswith('dbc')}
        PV = getChemFC(chem,feat=feat,val=pv_va,pivot=True,                         what='cr',**kw)

    if feat=='FA0':
        C1 = ['MF','St','ES','LM','DT','NS1','NS2','Ap','CN1','CN2']
    elif feat=='FA1':
        C1 = ['MF','St','ES','LM','DT','NS','Ap','CN']
    
    Ind=sorted(set([(i[0],i[1]) for i in list(FCp.index)]),
           key=lambda x: x[1])

    nr,nc=1,6
    #nc=5 if feat=='FA0' else 4
    
    fig,AX = pl.subplots(nrows=nr,ncols=nc,figsize=(nc*3,nr*3),
                         sharex=True,sharey=True)

    for ir in range(nr):
        for ic in range(nc):
            ii = ir*nc + ic
            
            ind= Ind[ii]
            X = FCp.loc[ind]
            X = X[C1].loc[X.index[::-1]]
            ax = AX[ic]
            
            if pv:
                Y = PV.loc[ind]
                Y = Y[C1].loc[Y.index[::-1]]
                Y1 = pd.DataFrame(np.array([['']*Y.shape[1]]*Y.shape[0]),
                                  columns=Y.columns)
                Y1.index=Y.index
                Y1[Y<pv]='*'
                g = sns.heatmap(X,vmin=-vm,vmax=vm,cmap=cm.RdYlBu_r,
                                annot=Y1,fmt='s',
                                lw=1,ax=ax,cbar=False)
            else:
                g = sns.heatmap(X,vmin=-vm,vmax=vm,cmap=cm.RdYlBu_r,
                                lw=1,ax=ax,cbar=False)
                
            if ir==nr-1: 
                ax.set_xlabel('endpoints')
            else: 
                ax.set_xlabel('')
            if ic==0 : ax.set_ylabel('$\mu$M')
            ax.set_title(r'%s %ih' % (chem[:10],ind[1]))

    pl.subplots_adjust(wspace=0.3,hspace=0.3,top=0.9,bottom=0.1)
    
    return fig