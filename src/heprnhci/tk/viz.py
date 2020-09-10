import pylab as pl
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

pl.style.use('seaborn-white')


def drawInvitroInvivo(X_cr,fig_sz=(15,8),fig_file=None,
                      cols=['grey','orange','firebrick','green'],
                      adverse=True,
                      ax1=None,ax2=None
                     ):

    xmin = 1e-2
    xmax = 1e4

    fig=pl.figure(figsize=fig_sz)
    gs1 = gridspec.GridSpec(1,4,width_ratios=[0.2,1,0.2,0.1])

    # Chemical names
    ax = pl.subplot(gs1[0,0])
    for i in range(X_cr.shape[0]):
        ax.text(20,i,X_cr.name[i],ha='right',va='center',color='black',
                fontproperties=fm.FontProperties(size=11))
    ax.set_ylim(-1,X_cr.shape[0])
    ax.set_xlim(0,20)
    ax.set_axis_off()

    # doses
    col_tx,col_heprn,col_loael,col_hepatg = cols
    ax = pl.subplot(gs1[0,1])
    ax.hlines(list(range(X_cr.shape[0])),xmin,xmax,colors='grey',lw=0.5)
    ax.vlines(np.logspace(-3,4,num=8),-1,X_cr.shape[0],colors='grey',lw=0.5)

    # TX
    ax.scatter(X_cr.tx_ac50_p50,np.arange(X_cr.shape[0])+0.25,
               s=25,c=col_tx,marker='D',
               label=r'$POD_{ToxCast}$',lw=0,
               zorder=5,alpha=1)

    ax.hlines(np.arange(X_cr.shape[0])+0.25,
              X_cr.tx_ac50_p05,X_cr.tx_ac50_p95,colors=col_tx,
              lw=1,alpha=1)


    # HepRn
    ax.scatter(X_cr.heprn_ac50_p50,np.arange(X_cr.shape[0])-0.25,s=30,
               c=col_heprn,marker='o',
               label=r'$POD_{HepRn}$',lw=0,
               zorder=5,alpha=1)
    ax.hlines(np.arange(X_cr.shape[0])-0.25,X_cr.heprn_ac50_p05,X_cr.heprn_ac50_p95,
              colors=col_heprn,lw=1,alpha=1)

    # HepRnAtg
    ax.scatter(X_cr.heprnatg_ac50_p50,
               np.arange(X_cr.shape[0])-0.25,s=30,
               c=col_hepatg,marker='o',
               label=r'$POD_{HepRnAtg}$',lw=0,
               zorder=5,alpha=1)
    ax.hlines(np.arange(X_cr.shape[0])-0.25,
              X_cr.heprnatg_ac50_p05,X_cr.heprnatg_ac50_p95,
              colors=col_hepatg,lw=1,alpha=1)
    
    # LOAEL
    ax.scatter(X_cr.loael_p50,list(range(X_cr.shape[0])),s=80,
               c=col_loael,marker='<',lw=0,
               zorder=5,alpha=0.8,label=r'$POD_{sub}$')
    ax.hlines(np.arange(X_cr.shape[0]),X_cr.loael_p05,X_cr.loael_p95,
              colors=col_loael,lw=1.5,alpha=1)

    ax.set_xscale('log')
    ax.set_ylim(-1,X_cr.shape[0])
    ax.set_xlim(xmin,xmax)
    ax.set_xlabel(r"Dose [mg/kg/day]")
    ax.set_xticklabels(['0.001','0.01','0.1','1','10','100','1000','10000'])

    #for tick in ax.get_yticklines(): 
    #    tick.set_visible(True)
    for tick in ax.get_yticklabels(): 
        tick.set_visible(False)
    ax.legend(loc=3,fontsize='small',fancybox=True,
              framealpha=1.0,facecolor='white')
                  
    pl.subplots_adjust(wspace=0.05,hspace=0.03,top=0.8,bottom=0.1)
    if fig_file:
        fig.savefig(fig_file,bbox_inches='tight',dpi=600)
        
def drawInvitroInvivo4(X_cr,fig_sz=(15,8),fig_file=None,
                      adverse=True):

    xmin = 1e-2
    xmax = 1e4

    fig=pl.figure(figsize=fig_sz)
    gs1 = gridspec.GridSpec(1,4,width_ratios=[0.2,1,0.2,0.1])

    # Chemical names
    ax = pl.subplot(gs1[0,0])
    for i in range(X_cr.shape[0]):
        ax.text(20,i,X_cr.name[i],ha='right',va='center',color='black',
                fontproperties=fm.FontProperties(size=11))
    ax.set_ylim(-1,X_cr.shape[0])
    ax.set_xlim(0,20)
    ax.set_axis_off()

    # doses
    col_tx,col_heprn,col_loael = 'green','blue','firebrick'
    ax = pl.subplot(gs1[0,1])
    ax.hlines(list(range(X_cr.shape[0])),xmin,xmax,colors='grey',lw=0.5)
    ax.vlines(np.logspace(-3,4,num=8),-1,X_cr.shape[0],colors='grey',lw=0.5)

    # TX
    ax.scatter(X_cr.tx_ac50_p50,np.arange(X_cr.shape[0])+0.25,
               s=25,c=col_tx,marker='D',
               label=r'$AED_{toxcast}$',lw=0,
               zorder=5,alpha=1)

    ax.hlines(np.arange(X_cr.shape[0])+0.25,
              X_cr.tx_ac50_p05,X_cr.tx_ac50_p95,colors=col_tx,
              lw=1,alpha=1)


    # HepRn
    ax.scatter(X_cr.heprn_ac50_p50,np.arange(X_cr.shape[0])-0.25,s=30,
               c=col_heprn,marker='o',
               label=r'$AED_{heprn}$',lw=0,
               zorder=5,alpha=1)
    ax.hlines(np.arange(X_cr.shape[0])-0.25,X_cr.heprn_ac50_p05,X_cr.heprn_ac50_p95,
              colors=col_heprn,lw=1,alpha=1)
    
    # LOAEL
    ax.scatter(X_cr.loael_p50,list(range(X_cr.shape[0])),s=80,
               c=col_loael,marker='<',lw=0,
               zorder=5,alpha=0.8,label=r'LOAEL')
    ax.hlines(np.arange(X_cr.shape[0]),X_cr.loael_p05,X_cr.loael_p95,
              colors=col_loael,lw=1.5,alpha=1)

    ax.set_xscale('log')
    ax.set_ylim(-1,X_cr.shape[0])
    ax.set_xlim(xmin,xmax)
    ax.set_xlabel(r"Dose [mg/kg/day]")
    ax.set_xticklabels(['0.001','0.01','0.1','1','10','100','1000','10000'])

    #for tick in ax.get_yticklines(): 
    #    tick.set_visible(True)
    for tick in ax.get_yticklabels(): 
        tick.set_visible(False)
    ax.legend(loc=3,fontsize='small',fancybox=True,
              framealpha=1.0,facecolor='white')
                  
    pl.subplots_adjust(wspace=0.05,hspace=0.03,top=0.8,bottom=0.1)
    if fig_file:
        fig.savefig(fig_file,bbox_inches='tight',dpi=600)
        

def drawInvitroInvivo2(X_cr,fig_sz=(15,8),fig_file=None,
                      adverse=True):

    xmin = -3
    xmax = 4

    fig=pl.figure(figsize=fig_sz)
    gs1 = gridspec.GridSpec(1,4,width_ratios=[0.2,1,0.2,0.1])

    # Chemical names
    ax = pl.subplot(gs1[0,0])
    for i in range(X_cr.shape[0]):
        ax.text(20,i,X_cr.preferred_name[i],ha='right',va='center',color='black',
                fontproperties=fm.FontProperties(size=10))
    ax.set_ylim(-1,X_cr.shape[0])
    ax.set_xlim(0,20)
    ax.set_axis_off()

    # Concentrations
    col_ccr,col_lec,col_loael,col_noael = 'goldenrod','dodgerblue','firebrick','forestgreen'
    ax = pl.subplot(gs1[0,1])
    ax.hlines(list(range(X_cr.shape[0])),xmin,xmax,colors='grey',lw=0.5)
    ax.vlines(list(range(xmin,xmax)),0,X_cr.shape[0],colors='grey',lw=0.5)

    # Ccr
    ax.scatter(X_cr.c_cr_p50,list(range(X_cr.shape[0])),s=X_cr.cr_conf*50,c=col_ccr,marker='D',
               label=r'$OED_{cr}$',lw=0,
               zorder=5,alpha=1)
    #X_cr['C_cr_lo'] =X_cr['C_cr_lo'].apply(lambda x: x if x>0.1 else 0.1) 
    ax.hlines(list(range(X_cr.shape[0])),X_cr.c_cr_p05,X_cr.c_cr_p95,colors=col_ccr,
              lw=1.5,alpha=1)

    # LEC
    ax.scatter(X_cr.lec_min,list(range(X_cr.shape[0])),s=30,c=col_lec,marker='s',lw=0,
               zorder=5,alpha=0.9,label=r'$OED_{LEC}$')
    # LOAEL
    if adverse:
        ax.scatter(X_cr.loael,list(range(X_cr.shape[0])),s=50,
                   c=col_loael,marker='<',lw=0,
                   zorder=5,alpha=0.8,label=r'LOAEL')

        # LOAEL
        ax.scatter(X_cr.noael,list(range(X_cr.shape[0])),s=50,
                   c=col_noael,marker='>',lw=0,
                   zorder=5,alpha=0.9,label=r'NOAEL')
    else:
        ax.scatter(X_cr.lel,list(range(X_cr.shape[0])),s=50,
                   c=col_loael,marker='<',lw=0,
                   zorder=5,alpha=0.8,label=r'LEL')

        # LOAEL
        ax.scatter(X_cr.nel,list(range(X_cr.shape[0])),s=50,
                   c=col_noael,marker='>',lw=0,
                   zorder=5,alpha=0.9,label=r'NEL')

    #ax.set_xscale('log')
    ax.set_ylim(-1,X_cr.shape[0])
    ax.set_xlim(xmin,xmax)
    ax.set_xlabel(r"Oral dose [mg/kg/day]")
    ax.set_xticklabels(['0.001','0.01','0.1','1','10','100','1000'])

    #for tick in ax.get_yticklines(): 
    #    tick.set_visible(True)
    for tick in ax.get_yticklabels(): 
        tick.set_visible(False)
    ax.legend(loc=3,fontsize='small',fancybox=True,
              framealpha=1.0,facecolor='white')
                  
    pl.subplots_adjust(wspace=0.05,hspace=0.03,top=0.8,bottom=0.1)
    if fig_file:
        fig.savefig(fig_file,bbox_inches='tight',dpi=600)
        
        
def drawInvitroInvivo1(X_cr,fig_sz=(15,8),fig_file=None,
                      adverse=True):

    xmin = -3
    xmax = 4

    fig=pl.figure(figsize=fig_sz)
    gs1 = gridspec.GridSpec(1,4,width_ratios=[0.2,1,0.2,0.1])

    # Chemical names
    ax = pl.subplot(gs1[0,0])
    for i in range(X_cr.shape[0]):
        ax.text(20,i,X_cr.preferred_name[i],ha='right',va='center',color='black',
                fontproperties=fm.FontProperties(size=10))
    ax.set_ylim(-1,X_cr.shape[0])
    ax.set_xlim(0,20)
    ax.set_axis_off()

    # Concentrations
    col_ccr,col_lec,col_loael,col_noael = 'goldenrod','dodgerblue','firebrick','forestgreen'
    ax = pl.subplot(gs1[0,1])
    ax.hlines(list(range(X_cr.shape[0])),xmin,xmax,colors='grey',lw=0.5)
    ax.vlines(list(range(xmin,xmax)),0,X_cr.shape[0],colors='grey',lw=0.5)

    # Ccr
    ax.scatter(X_cr.c_cr_p50,list(range(X_cr.shape[0])),s=X_cr.cr_conf*50,c=col_ccr,marker='D',
               label=r'$OED_{cr}$',lw=0,
               zorder=5,alpha=1)
    #X_cr['C_cr_lo'] =X_cr['C_cr_lo'].apply(lambda x: x if x>0.1 else 0.1) 
    ax.hlines(list(range(X_cr.shape[0])),X_cr.c_cr_p05,X_cr.c_cr_p95,colors=col_ccr,
              lw=1.5,alpha=1)

    # LEC
    ax.scatter(X_cr.lec_min,list(range(X_cr.shape[0])),s=30,c=col_lec,marker='s',lw=0,
               zorder=5,alpha=0.9,label=r'$OED_{LEC}$')
    # LOAEL
    if adverse:
        ax.scatter(X_cr.loael,list(range(X_cr.shape[0])),s=50,
                   c=col_loael,marker='<',lw=0,
                   zorder=5,alpha=0.8,label=r'LOAEL')

        # LOAEL
        ax.scatter(X_cr.noael,list(range(X_cr.shape[0])),s=50,
                   c=col_noael,marker='>',lw=0,
                   zorder=5,alpha=0.9,label=r'NOAEL')
    else:
        ax.scatter(X_cr.lel,list(range(X_cr.shape[0])),s=50,
                   c=col_loael,marker='<',lw=0,
                   zorder=5,alpha=0.8,label=r'LEL')

        # LOAEL
        ax.scatter(X_cr.nel,list(range(X_cr.shape[0])),s=50,
                   c=col_noael,marker='>',lw=0,
                   zorder=5,alpha=0.9,label=r'NEL')

    #ax.set_xscale('log')
    ax.set_ylim(-1,X_cr.shape[0])
    ax.set_xlim(xmin,xmax)
    ax.set_xlabel(r"Oral dose [mg/kg/day]")
    ax.set_xticklabels(['0.001','0.01','0.1','1','10','100','1000'])

    #for tick in ax.get_yticklines(): 
    #    tick.set_visible(True)
    for tick in ax.get_yticklabels(): 
        tick.set_visible(False)
    ax.legend(loc=3,fontsize='small',fancybox=True,
              framealpha=1.0,facecolor='white')
                  
    pl.subplots_adjust(wspace=0.05,hspace=0.03,top=0.8,bottom=0.1)
    if fig_file:
        fig.savefig(fig_file,bbox_inches='tight',dpi=600)