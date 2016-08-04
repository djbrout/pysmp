import os
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
np.set_printoptions(threshold=np.nan)
import sys



def lightcurve(mjd,fitmag,fitmagerr,fakemag,fakemagerr,fitflux,fitfluxerr,fakeflux,fakefluxerr,filter,filename,title=''):
    plt.clf()
    fig,ax = plt.subplots(3,1,sharex=True,figsize=(7,10))
    
    ax[0].set_title(title)
    ax[0].set_xlim(min(mjd[mjd>0])-50,min(mjd[mjd>0])+250.)

    ax[0].scatter(fakemag*0.,fakemag,color='black',marker='+',label='Fake Mag')
    ax[0].scatter(fitmag*0.,fitmag,color='black',label='Fit Mag')
    ax[0].invert_yaxis()
    ax[0].set_ylabel('Mag')
    ax[0].legend(fontsize=10)
    ff = fakemag[(fakemag<50) & (fakemag>0)]
    ax[0].set_ylim(max(ff[ff<50])+1.,min(ff[ff<50]-1.))

    for f in np.unique(filter):
        if f == 'g':
            color = 'green'
        if f == 'r':
            color = 'red'
            print 'redddddd'
            raw_input()
        if f == 'i':
            color = 'blue'
        if f == 'z':
            color = 'pink'
        ww = (filter == f) & (abs(fitmag - fakemag) < 10.)
        ax[0].errorbar(mjd[ww],fakemag[ww],fakemagerr[ww],color=color,marker='+')
        ax[0].errorbar(mjd[ww],fitmag[ww],fitmagerr[ww],color=color,fmt='o',alpha=.3)
        ff = fakemag[ww]
        ax[1].errorbar(mjd[ww],fakemag[ww]-fitmag[ww],fitmagerr[ww],color=color,fmt='o')
        ax[2].scatter(mjd[ww],(fakeflux[ww]-fitflux[ww])/fitfluxerr[ww],color=color)
 
    s = (fakeflux-fitflux)/fitfluxerr
    r = fakemag-fitmag
    try:
        ax[1].set_ylim(min(r[r>-3.])-.5,max(r[r<3])+.5)
    except:
        ax[1].set_ylim(-1.,1.)
    try:
        ax[2].set_ylim(min(s[s>-3.])+.25,max(s[s<3])-.25)
    except:
        ax[2].set_ylim(-2.,2.)

    ax[1].set_xlim(min(mjd[mjd>0])-50,min(mjd[mjd>0])+250.)
    ax[1].set_ylabel('Fake - Fit Mag')
    ax[1].plot([min(mjd),max(mjd)],[0,0],color='black')
    ax[2].plot([min(mjd),max(mjd)],[0,0],color='black')
    ax[2].set_ylabel('Fake - Fit Flux / Err')
    ax[2].set_xlim(min(mjd[mjd>0])-50,min(mjd[mjd>0])+250.)
    ax[2].set_xlabel('MJD')


    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)
    plt.setp([a.get_xticklabels() for a in fig.axes[0:1]], visible=False)

    plt.savefig(filename)
    print filename


def pslightcurve(mjd, fitmag, fitmagerr, fakemag, fakemagerr, fitflux, fitfluxerr, fakeflux, fakefluxerr, filter,
               filename, title=''):
    plt.clf()
    fig, ax = plt.subplots(3, 1, sharex=True, figsize=(7, 10))

    ax[0].set_title(title)
    ax[0].set_xlim(min(mjd[mjd > 0]) - 50, min(mjd[mjd > 0]) + 250.)

    ax[0].scatter(fakemag * 0., fakemag, color='black', marker='+', label='Fake Mag')
    ax[0].scatter(fitmag * 0., fitmag, color='black', label='Fit Mag')
    ax[0].invert_yaxis()
    ax[0].set_ylabel('Mag')
    ax[0].legend(fontsize=10)
    ff = fakemag[(fakemag < 50) & (fakemag > 0)]
    ax[0].set_ylim(max(ff[ff < 50]) + 1., min(ff[ff < 50] - 1.))

    for f in np.unique(filter):
        if f == 'g':
            color = 'green'
        if f == 'r':
            color = 'red'
            print 'redddddd'
            raw_input()
        if f == 'i':
            color = 'blue'
        if f == 'z':
            color = 'pink'
        #ww = (filter == f) & (abs(fitmag - fakemag) < 10.)
        aa=np.ones(len(mjd))
        ww = aa == 1
        ax[0].scatter(mjd,fakemag)
        ax[0].errorbar(mjd[ww], fakemag[ww], fakemagerr[ww], color=color, marker='+')
        ax[0].errorbar(mjd[ww], fitmag[ww], fitmagerr[ww], color=color, fmt='o', alpha=.3)
        ff = fakemag[ww]
        ax[1].errorbar(mjd[ww], fakemag[ww] - fitmag[ww], fitmagerr[ww], color=color, fmt='o')
        ax[2].scatter(mjd[ww], (fakeflux[ww] - fitflux[ww]) / fitfluxerr[ww], color=color)

    s = (fakeflux - fitflux) / fitfluxerr
    r = fakemag - fitmag
    try:
        ax[1].set_ylim(min(r[r > -3.]) - .5, max(r[r < 3]) + .5)
    except:
        ax[1].set_ylim(-1., 1.)
    try:
        ax[2].set_ylim(min(s[s > -3.]) + .25, max(s[s < 3]) - .25)
    except:
        ax[2].set_ylim(-2., 2.)

    ax[1].set_xlim(min(mjd[mjd > 0]) - 50, min(mjd[mjd > 0]) + 250.)
    ax[1].set_ylabel('Fake - Fit Mag')
    ax[1].plot([min(mjd), max(mjd)], [0, 0], color='black')
    ax[2].plot([min(mjd), max(mjd)], [0, 0], color='black')
    ax[2].set_ylabel('Fake - Fit Flux / Err')
    ax[2].set_xlim(min(mjd[mjd > 0]) - 50, min(mjd[mjd > 0]) + 250.)
    ax[2].set_xlabel('MJD')

    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)
    plt.setp([a.get_xticklabels() for a in fig.axes[0:1]], visible=False)

    plt.savefig(filename)
    print filename

def lightcurveoverlay(mjd,fitmag,fitmagerr,fakemag,fitflux,fitfluxerr,fakeflux,filter,filename,compare,title=''):
    plt.clf()
    fig,ax = plt.subplots(3,1,sharex=True,figsize=(7,10))

    ax[0].set_title(title)
    ax[0].set_xlim(min(mjd[mjd>0])-50,min(mjd[mjd>0])+250.)

    ax[0].scatter(fakemag*0.,fakemag,color='black',marker='+',label='Fake Mag')
    ax[0].scatter(fitmag*0.,fitmag,color='black',label='Fit Mag')
    ax[0].invert_yaxis()
    ax[0].set_ylabel('Mag')
    ax[0].legend(fontsize=10)
    ff = fakemag[(fakemag<50) & (fakemag>0)]
    ax[0].set_ylim(max(ff[ff<50])+1.,min(ff[ff<50]-1.))
    ax[0].axhline(20)
    colors = ['green','red','blue','orange']
    
    for i,c in enumerate(np.unique(compare)):
        color = colors[i]
        print c
        ww = (abs(fitmag - fakemag) < 10.) & (compare == c) & (fitmagerr < 1.)
        ax[0].scatter(mjd[ww],fakemag[ww],color='black',marker='+')
        ax[0].errorbar(mjd[ww],fitmag[ww],fitmagerr[ww],color=color,fmt='o',alpha=.4,label=c)
        ff = fakemag[ww]
        ax[1].errorbar(mjd[ww],fakemag[ww]-fitmag[ww],fitmagerr[ww],color=color,alpha=.4,fmt='o')
        ax[2].scatter(mjd[ww],(fakeflux[ww]-fitflux[ww])/fitfluxerr[ww],color=color,alpha=.4)

    s = (fakeflux-fitflux)/fitfluxerr
    r = fakemag-fitmag
    ax[0].legend(fontsize=10)
    try:
        ax[1].set_ylim(min(r[r>-3.])-.5,max(r[r<3])+.5)
    except:
        ax[1].set_ylim(-1.,1.)
    try:
        ax[2].set_ylim(min(s[s>-3.])+.25,max(s[s<3])-.25)
    except:
        ax[2].set_ylim(-2.,2.)

    ax[1].set_xlim(min(mjd[mjd>0])-50,min(mjd[mjd>0])+250.)
    ax[1].set_ylabel('Fake - Fit Mag')
    ax[1].plot([min(mjd),max(mjd)],[0,0],color='black')
    ax[2].plot([min(mjd),max(mjd)],[0,0],color='black')
    ax[2].set_ylabel('Fake - Fit Flux / Err')
    ax[2].set_xlim(min(mjd[mjd>0])-50,min(mjd[mjd>0])+250.)
    ax[2].set_xlabel('MJD')


    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)
    plt.setp([a.get_xticklabels() for a in fig.axes[0:1]], visible=False)

    plt.savefig(filename)
    print filename
