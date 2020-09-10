import pandas as pd
import numpy as np
import scipy as sp
import pylab as pl
from collections import *
from copy import copy
import os
import subprocess
import numpy.linalg as LA

from statsmodels import robust 
import seaborn as sns
from matplotlib import cm
from scipy.stats import iqr

def calcPct(plate_id,dbc_pc=None,dbc_sm=None,
           ctrl_neg='ctrl-',ctrl_pos='ctrl+',
           test_chem='chem_test',ft_name='FA0',
           smooth='smooth_w7',
           normal='all'):
    """
    Use plate level data to generate percentage change value
    
    KWArgs:
    normal = median plate median / 
             ctrl-  negative control
             
    """
    def p02(X): return np.percentile(X,2)
    def p98(X): return np.percentile(X,98)
    def p50(X): return np.percentile(X,50)

    Plate = dbc_sm.find_one(dict(plate_id=plate_id),dict(_id=0))
    Wells=pd.DataFrame(Plate.pop('wells'))
    Smooth=pd.DataFrame(Plate.get('smooth'))

    Wells.loc[:,'lconc']=np.log10(Wells.conc)
    P_neg=Wells[Wells.stype==ctrl_neg]
    P_pos =Wells[Wells.stype==ctrl_pos]
    P_chem=Wells[Wells.stype==test_chem]

    ST_neg=P_neg.groupby([ft_name])\
                .aggregate(dict(raw_value=[np.mean,np.median,np.std,
                                           p02,p98,p50]))
    ST_neg.columns=['ctrl_'+i[1] for i in ST_neg.columns]
    ST_neg = ST_neg.reset_index()
    
    ST_plt =Smooth.groupby([ft_name])\
                     .aggregate({smooth:[np.mean,np.median,
                                         np.std,p02,p98,p50]})
    
    ST_plt.columns=['plt_'+i[1] for i in ST_plt.columns]
    ST_plt = ST_plt.reset_index()

    PC=Smooth.merge(ST_neg[[ft_name,'ctrl_mean','ctrl_std','ctrl_median',
                            'ctrl_p02','ctrl_p98']],
                    on=ft_name)\
        .merge(ST_plt[[ft_name,'plt_mean','plt_std','plt_median',
                       'plt_p02','plt_p98']],
               on=ft_name)

    PC['ctrl_pct']= (PC[smooth]-PC.ctrl_mean)/PC.ctrl_mean
    PC['plt_pct'] = (PC[smooth]-PC.plt_mean)/PC.plt_mean
    PC['ctrl_pct1']= PC[smooth]/PC.ctrl_mean
    PC['plt_pct1'] = PC[smooth]/PC.plt_mean
    
    PC['ctrl_pv']  =sp.stats.norm.cdf(PC[smooth],PC.ctrl_mean,PC.ctrl_std)
    PC['plt_pv']   =sp.stats.norm.cdf(PC[smooth],PC.plt_mean,PC.plt_std)
    
    
    Plate['dmso_stats']  = ST_neg.reset_index().to_dict('records')
    Plate['plate_stats']  = ST_plt.reset_index().to_dict('records')    
    Plate['chem_pc'] = PC.reset_index().to_dict('records')
    
    if dbc_pc:
        dbc_pc.insert_one(Plate)
    else:
        return Plate

def getChemPC(chem,feat='FA0',val='ctrl_pct',pivot=False,
              ind=None,
              what='traj',add_t0=False,use_times=None,
              dbc_pc=None,dbc_ft=None):
    
    Q = dict()
    
    Features=pd.DataFrame(list(dbc_ft.find({},dict(_id=0))))
    PC1=[[j for j in i['chem_pc'] if j['name']==chem] 
             for i in 
                dbc_pc.find({'chem_pc':{'$elemMatch':dict(name=chem)}})]
    PC1=pd.concat(list(map(pd.DataFrame,PC1)))\
        .merge(Features,left_on='FA0',right_on='FA0')

    if use_times:
            PC1 = PC1[(PC1.timeh.isin(use_times))]

    I1 = ind if ind else ['name','conc','timeh']

    if what=='traj' and not ind:
        I1 = ['name','conc','timeh']
    else:
        I1 = ['name','timeh','conc']
        
    if not pivot: 
        return PC1
    
    if feat=='FA0':
        PC2=PC1.pivot_table(index=I1,
                           columns='FA0',
                           values=val,fill_value=1.0)
    elif feat=='FA1':
        PC2=PC1.pivot_table(index=I1,
                           columns='FA1',
                           values=val,fill_value=1.0,
                           aggfunc=np.mean)
        
    if add_t0 and what=='traj':
        Ind1=list(set([(i[0],i[1]) for i in PC2.index]))
    
        for (ch,co) in Ind1:
            PC2.loc[(ch,co,0)]=0

        PC2.sort_index(inplace=True)

    
    return PC2
