# Import libraries
import sys
import os
import string
import linecache
from math import *
from mpl_toolkits import mplot3d
from mpl_toolkits.axes_grid1.inset_locator import (InsetPosition, mark_inset)
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt

def parsefname(fnm):
    "Extracts file name from path"
    parsed = fnm.split('/')
    return parsed[-1]

def file2buffer(fnm):
    "Reads file into buffer"
    f = open(fnm,'r')
    buf = []
    for line in f.readlines():
        buf.append(line)
    f.close()
    return buf

def extractdatacols(buf):
    "Gets data columns from buffer"
    coli=[]
    coef=[]
    colx=[]
    colf=[]
    for ii in range(len(buf)):
        buftmp=buf[ii].split()
        coli.append(int(buftmp[0]))
        coef.append(float(buftmp[1]))
        colx.append(float(buftmp[2]))
        colf.append(float(buftmp[3]))
    return coli,coef,colx,colf

def extractchebdata(buf):
    "Gets bounds and Chebyshev coefs"
    # Bounds: Chebyshev coordinate mappings
    btmp=buf[1].split()
    bound = [float(btmp[ii]) for ii in range(len(btmp))]
    # Filter parameters: either homo/lumo or mu 
    ptmp=buf[3].split()
    param = [float(ptmp[ii]) for ii in range(len(ptmp))]
    # Filter Chebyshev coefficients
    coef=[]
    for ii in range(5,len(buf)):
        buftmp=buf[ii].split()
        coef.append(float(buftmp[0]))
    return bound,param,coef

def evalcheb(x,coefs,trunc):
    "Evaluate Chebyshev series via Clenshaw recursion"
    tx=2.0*x
    ov=0.0
    nv=0.0
    for ii in range(trunc-1,0,-1):
        sv=nv
        nv=tx*nv-ov+coefs[ii]
        ov=sv
    fx=x*nv-ov+coefs[0]
    return fx

def genxvals(N):
    "Generates evenly spaced points over [-1,1]"
    inc = 2.0/(N-1)
    xvals = [-1+ii*inc for ii in range(N)]
    return xvals

def getinsetranges(f):
    "Gets parameters for plotting the inset"
    st=0
    fi=len(f)-1
    tol=0.005
    # Start and end of filter gap
    for ii in range(len(f)):
        if abs(f[ii]-1.0) > tol: break
        st=ii
    for ii in range(len(f)-1,-1,-1):
        if abs(f[ii]) > tol: break
        fi=ii
    stb=max(0,st - (fi+1-st)//2)
    fib=min(len(f),fi + (fi+1-st)//2 + 1)
    return [stb,fib]

def plotthedata(i,c,x,f,t,p,girg):
    "Create MatPlotLib plot of data"
    fig, axs = plt.subplots(2,1)
    # Chebyshev coefficient magnitudes
    axs[0].plot(i,c)
    axs[0].vlines([t],0,1,transform=axs[0].get_xaxis_transform(), colors='r')
    axs[0].set_xlim(i[0], i[-1])
    axs[0].set_xlabel('i_coef')
    axs[0].set_ylabel('log_10 |coef|')
    axs[0].grid(True)
    # Filter reconstructed from coefficients
    axs[1].plot(x,f)
    axs[1].set_xlim(min(x),max(x))
    axs[1].set_ylim(min(-0.1,min(f)-0.1),max(1.1,max(f)+0.1))
    axs[1].set_xlabel('E (hartree)')
    axs[1].set_ylabel('filter')
    # Inset (filter cutoff region)
    if girg[1]-girg[0] > 0:
       ax1ins = plt.axes([0,0,1,1])
       ip = InsetPosition(axs[1], [0.45,0.35,0.5,0.5])
       ax1ins.set_axes_locator(ip)
       ax1ins.minorticks_on()
       ax1ins.grid(True)
       xr=[x[ii] for ii in range(girg[0],girg[1])]
       fr=[f[ii] for ii in range(girg[0],girg[1])]
       ax1ins.plot(xr,fr, color='tab:blue', linewidth=0.1, alpha=0.7)
       ax1ins.plot(xr,fr, color='tab:blue', marker='o', markersize=2)
       ax1ins.vlines(p,0,1,transform=ax1ins.get_xaxis_transform(), colors='magenta')
       mark_inset(axs[1], ax1ins, loc1=2, loc2=4, fc="none", ec='0.5')
    plt.show()
    return

def makethegraph(fnm):
    "Top routine for making graph"
    data = file2buffer(fnm)
    b,p,c=extractchebdata(data)
    # N is oversampled energy axis for filter
    N=len(c)*4
    # t is number of coefficients to use (can be truncated at t<len(c))
    t=len(c)
    z=genxvals(N)
    x=[z[ii]*b[1]+b[0] for ii in range(N)]
    f=[evalcheb(z[ii],c,t) for ii in range(N)]
    i=[ii for ii in range(len(c))]
    girg = getinsetranges(f)
    cl = [log10(abs(c[ii])) for ii in range(len(c))]
    plotthedata(i,cl,x,f,t,p,girg)
    return

if len(sys.argv) == 2:
   fnm = sys.argv[1]
   makethegraph(fnm)
else:
   print("syntax: '$> %s <filter data file name>'" %(sys.argv[0]))

