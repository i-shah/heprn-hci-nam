import pylab as pl
import matplotlib
from matplotlib.collections import PatchCollection
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.text as text
import matplotlib.font_manager as fm 
from matplotlib.patches import Ellipse, Circle
import matplotlib.cm as cm
from matplotlib import gridspec
from scipy.interpolate import splev, splrep
import seaborn as sns
import numpy as np
import pandas as pd
from box import Box as ddict

pl.style.use('seaborn-white')

    
def drawIVIVE(LPR,POD_nmtx,POD_effect,POD_hcift,
              fig_text=ddict(a=ddict(title=r'(a) $POD_{nam}$ vs $POD_{tox}$',
                                     xlab='Dose [mg/kg/day]'),
                             b=ddict(title=r'(b) Hepatic LOAELs',
                                     xlab='Hepatic effects'),
                             c=ddict(title=r'(c) HCI AEDs',
                                     xlab='HCI Endpoints')
                            ),
              Conci=ddict(heprn_max_conc=ddict(color='mediumblue',
                                     marker='<',
                                     dy=0,z=25,size=30,lw=0.8,
                                     label='$Max_{HepRn}$')
                        ),
              Podi=ddict(heprn_ac50=ddict(color='mediumblue',
                                     marker='D',
                                     dy=0,z=30,size=40,lw=0.8,
                                     label='$POD_{HepRn}$'),
                         heprnatg_ac50=ddict(color='blueviolet',
                                       marker='d',
                                       dy=0.2,z=8,size=30,lw=0.5,
                                       label='$POD_{HepRnAtg}$'),
                        tx_ac50=ddict(color='forestgreen',
                                 marker='8',
                                 dy=-0.2,z=5,size=30,lw=0.4,
                                 label='$POD_{ToxCast}$'),
                        loael=ddict(color='crimson',
                                 marker='s',
                                 dy=0,z=15,size=50,lw=0.5,
                                 label='$POD_{Tox}$'),
                       ), 
              lpr_info=ddict(LPR_heprn_ac50=ddict(color='mediumblue',
                                     label='$LPR_{HepRn}$'),
                    LPR_heprnatg_ac50=ddict(color='blueviolet',
                                   label='$LPR_{HepRnAtg}$'),
                    LPR_tx_ac50=ddict(color='forestgreen',
                             label='$LPR_{ToxCast}$')),
              cmap=ddict(tox=cm.hot,hts=cm.Blues_r),
              fig_sz=(17,15),fig_file=None,fig_dpi=600,
              pod_range=ddict(xmin=1e-2,xmax=1e4),
              fill_hm=False,
              bb_legend=[0.75,0.1,0.2,0.1]
             ):
    """
    Three panels with
    a) IVIVE: ranges for POD_nam and POD_tox for chemicals
    b) HTS: AED for potency values for significant effects
    c) TOX: Effect types with Dose values
    
    2x5 grid with 
    1,0: Chemical names
    1,1: IVIVE
    0,1: POD_nam:kde, POD_tox:kde
    0,2: LPR: kde
    1,2: LPR
    0,3: Tox:hist
    1,3: Tox
    0,4: HTS:hist
    1,4: HTS
    """

    fig=pl.figure(figsize=fig_sz)
    GS = gridspec.GridSpec(2,6,
                           width_ratios=[0.2,0.6,0,0.5,0.2,0.1],
                           height_ratios=[0.1,1])
    
    # Panel(a)   
    
    drawHistPOD(POD_nmtx,pl.subplot(GS[0,1]),pod_info=Podi,pod_range=pod_range,
             title=fig_text.a.title)
#     drawBoxPOD(POD_nmtx,pl.subplot(GS[0,1]),pod_info=Podi,pod_range=pod_range,
#              title=fig_text.a.title)    
    drawChems(LPR,pl.subplot(GS[1,0]))    
    drawPOD(POD_nmtx,pl.subplot(GS[1,1]),pod_info=Podi,
            conc_info=Conci,pod_range=pod_range,
            lg_bb=fig_text.a.lg_bb,
            xlab=fig_text.a.xlab)

    
    # Panel (b)    
#     drawHistLPR(LPR,ax=pl.subplot(GS[0,2]),
#                 lpr_info=lpr_info,
#                 title=fig_text.b.title)
#     drawLPR(LPR,ax=pl.subplot(GS[1,2]),lpr_info=lpr_info,cmap=cmap.lpr,
#             xlab=fig_text.b.xlab)
    
    # Panel (c)
    drawBar(POD_effect,pl.subplot(GS[0,3]),cmap=cmap.tox,
            title=fig_text.c.title)
    hm_tox=drawHM(POD_effect,pl.subplot(GS[1,3]),cmap=cmap.tox,
                  xlab=fig_text.c.xlab,
                  pod_range=pod_range,
                  fill=fill_hm)

    # Panel (d)
    drawBar(POD_hcift,pl.subplot(GS[0,4]),cmap=cmap.hts,
            title=fig_text.d.title)
    
    hm_hts=drawHM(POD_hcift,pl.subplot(GS[1,4]),cmap=cmap.hts,
                  xlab=fig_text.d.xlab,
                  pod_range=pod_range,
                  fill=fill_hm)
    
    # Add the legend for heatmaps
    ax = pl.subplot(GS[1,5])
    ax.set_axis_off()
    
    # HCI
    ax1 = fig.add_axes(bb_legend)
    ax1.set_axis_off()
    cb = pl.colorbar(hm_hts,ax=ax1)
    cb.ax.set_xlabel('Dose [mg/kg/day]')
    #cb.ax.set_yticklabels([0.1,1,10,100,1000,10000])
    
    pl.subplots_adjust(wspace=0.05,hspace=0.03,top=0.8,bottom=0.1)
    if fig_file:
        fig.savefig(fig_file,bbox_inches='tight',dpi=fig_dpi)

    
def drawChems(X0,ax,x0=20):
    X=X0.reset_index()
    for i in range(X.shape[0]):
        ax.text(x0,i+0.5,X.name[i],ha='right',va='center',color='black',
                fontproperties=fm.FontProperties(size=11))
    ax.set_ylim(0,X.shape[0])
    ax.set_xlim(0,x0)
    ax.set_axis_off()
    

def drawPOD(X,ax,pod_info=None,conc_info=None,
            pod_range=None,xlab=r"Dose [mg/kg/day]",title=None,lg_bb=None):
    yoff=0.5
    VL = np.logspace(np.log10(pod_range.xmin),np.log10(pod_range.xmax),num=pod_range.num)
    #XT = ['0.001','0.01','0.1','1','10','100','1000']
    #XT = 10**VL
    ax.hlines(np.arange(X.shape[0])+0.5,
              pod_range.xmin,pod_range.xmax,
              colors='grey',lw=0.5)
    ax.vlines(VL,-1,X.shape[0],colors='grey',lw=0.5)
    
    for k,info in pod_info.items():
        ax.scatter(X[k],np.arange(X.shape[0])+yoff+info.dy,
                   s=info.size,c=info.color,marker=info.marker,
                   label=info.label,lw=0,
                   zorder=info.z,alpha=1)        
#         ax.hlines(np.arange(X.shape[0])+yoff+info.dy,
#                   X[p05c],X[p95c],colors=info.color,
#                   zorder=info.z,
#                   lw=info.lw,alpha=0.6)

    for k,info in conc_info.items():
        ax.scatter(X[k],np.arange(X.shape[0])+yoff+info.dy,
                   s=info.size,c=info.color,marker=info.marker,
                   label=info.label,lw=0,
                   zorder=info.z,alpha=1)        


    ax.set_xscale('log')
    ax.set_ylim(0,X.shape[0])
    ax.set_xlim(pod_range.xmin,pod_range.xmax)
    ax.set_xlabel(xlab)
    if title: ax.set_title(title)
    
    #ax.set_xticklabels(VL,rotation=90)
    for tick in ax.get_yticklabels(): 
        tick.set_visible(False)
    ax.xaxis.tick_bottom()
    leg = ax.legend(loc=3,fontsize='medium',fancybox=False,
                    bbox_to_anchor=lg_bb,
                    framealpha=2,facecolor='grey')
    #leg.get_frame().set_facecolor('white')
    

def drawLPR(X,ax,lpr_info=None,cmap=None,title=None,xlab=None,fn_sz=10):
    K1 = list(lpr_info.keys())
    X1 = X[K1].reset_index(drop=True)
    #X1.index=[range(1,X1.shape[0]+1)]
    C1 = [lpr_info[k].color for k in K1]
    #X1.plot.barh(color=C1,ax=ax,legend=False,width=1.0,lw=0,alpha=0.7)    
    X1.plot.barh(color=C1,ax=ax,legend=False,width=0.9,lw=0,alpha=0.9,grid=True)    
    
    #ax.hlines(np.arange(X.shape[0])+0.5,X1.min().min(),X1.max().max(),colors='grey',lw=0.5)
    #ax.vlines([-1,0,1,2,3],0,len(X),colors='grey',lw=0.5)
    ax.set_xticks([-1,0,1,2,3])
    ax.set_xticklabels([-1,0,1,2,3],rotation=90)
    for tick in ax.get_xticklabels(): 
        tick.set_visible(True)
        tick.set_fontsize(10)
        
    for tick in ax.get_yticklines(): 
        tick.set_visible(False)
    for tick in ax.get_yticklabels(): 
        tick.set_visible(False)
    
    if xlab: ax.set_xlabel(xlab)
    if title: ax.set_title(title,fontsize=16)  
    ax.set_ylim(-0.5,X.shape[0]-0.5)
    #ax.set_xlim(xmin,-xmin)    
    

def drawHM(X,ax,cmap=None,xlab='Hepatic effects',fn_sz=10,fill=False,
           pod_range=None,title=None):
    ax.xaxis.tick_bottom()
    ax.set_ylim(0,X.shape[0])
    Nrm = matplotlib.colors.LogNorm(vmin=pod_range.xmin,vmax=pod_range.xmax)
    myCol = cm.ScalarMappable(Nrm,cmap=cmap)
    if fill: X = X.fillna(fill)
    hm=ax.pcolor(X,norm=Nrm,cmap=cmap,lw=1,edgecolors='#bcbcbc')
    
    #ax.set_axis_off()
    ax.set_xticks(np.arange(X.shape[1])+0.5, minor=False)
    xl=ax.set_xticklabels(X.columns,rotation=90)
    for tick in ax.get_xticklines(): 
        tick.set_visible(False)
    for tick in ax.get_xticklabels(): 
        tick.set_fontsize(fn_sz)

    for tick in ax.get_yticklines(): 
        tick.set_visible(False)
    for tick in ax.get_yticklabels(): 
        tick.set_visible(False)
    ax.set_xlabel(xlab)
    if title: ax.set_title(title)    
        
    return hm
    
def drawBar(X,ax,cmap=None,title=None):
    N = pd.DataFrame(X.notnull().sum()).reset_index()
    N.columns=['x','y']

    col=cmap(128)
    sns.barplot(x='x',y='y',data=N,ax=ax,color=col)
    for tick in ax.get_xticklines(): 
        tick.set_visible(False)
    for tick in ax.get_xticklabels(): 
        tick.set_visible(False)
    ax.set_axis_off()
    if title: ax.set_title(title,fontsize=16)    
#     ax.vlines(XT,0,len(N)+0.5,
#               linewidth=0.3,linestyle='-',color='#676767')
#     ax.hlines(np.arange(0,len(N)),0,15,
#               linewidth=0.3,linestyle='-',color='#676767')
#     ax.set_xticks(XT)
#     ax.set_xlabel('# chemicals')
#     ax.set_ylabel('')

def drawHistPOD(X,ax,pod_info=None,pod_range=None,title=None):    
    Bins = np.logspace(-3,4,num=50)
    
    D=[]
    C=[]
    for k,info in pod_info.items():
        Xi = X[k]
        D.append(Xi[Xi.notnull()])
        C.append(info.color)
    ax.hist(D,bins=Bins,histtype='bar',stacked=True,
            color=C,alpha=0.9,rwidth=0.9)
    ax.set_xscale('log')
    ax.set_xlim(pod_range.xmin,pod_range.xmax)
    ax.set_axis_off()    
    if title: ax.set_title(title,fontsize=16)        

def drawBarPOD(X,ax,pod_info=None,pod_range=None,title=None):    
    Bins = np.logspace(-3,4,num=50)
    
    X.plot.box()
    ax.hist(D,bins=Bins,histtype='bar',stacked=True,
            color=C,alpha=0.9,rwidth=0.9)
    ax.set_xscale('log')
    ax.set_xlim(pod_range.xmin,pod_range.xmax)
    ax.set_axis_off()    
    if title: ax.set_title(title,fontsize=16)        
        
def drawHistLPR(X,ax,lpr_info=None,title=None):  
    K1 = list(lpr_info.keys())
    X1 = X[K1]

    xmin,xmax=X1.min().min(),X1.max().max()
    Bins = np.linspace(xmin,xmax,num=10)
    
    D=[]
    C=[]
    for k,info in lpr_info.items():
        Xi = X[k]
        D.append(Xi[Xi.notnull()])
        C.append(info.color)
    ax.hist(D,bins=Bins,histtype='bar',stacked=True,
            color=C,alpha=0.9,rwidth=0.9)
    
    ax.set_xlim(xmin,xmax)
    ax.set_axis_off()    
    if title: ax.set_title(title,fontsize=16)        
        
def drawHistLPR1(X,ax,bins=10,title=None,cmap=None):    
    ax.hist(X[X.notnull()],bins=bins,histtype='bar',
            color=cmap(20),alpha=0.9,rwidth=0.9)
    #ax.set_xlim(pod_range.xmin,pod_range.xmax)
    ax.set_axis_off()    
    if title: ax.set_title(title,fontsize=16)                