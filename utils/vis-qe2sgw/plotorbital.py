# Import libraries
import sys
import os
import string
import linecache
from math import *
from mpl_toolkits import mplot3d
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt

# Table of elements for plotting
# [symbol, size, alpha, color]
elemtable = [
['H' , 40, 0.51, 'teal'],
['He', 40, 0.51, 'hotpink'],
['Li', 45, 0.51, 'silver'],
['Be', 45, 0.51, 'slategray'],
['B' , 45, 0.51, 'saddlebrown'],
['C' , 45, 0.51, 'black'],
['N' , 45, 0.51, 'blue'],
['O' , 45, 0.51, 'red'],
['F' , 45, 0.51, 'khaki'],
['Ne', 45, 0.51, 'orangered'],
['Na', 50, 0.51, 'silver'],
['Mg', 50, 0.51, 'slategray'],
['Al', 50, 0.51, 'gray'],
['Si', 50, 0.51, 'darkgray'],
['P' , 50, 0.51, 'firebrick'],
['S' , 50, 0.51, 'yellow'],
['Cl', 50, 0.51, 'yellowgreen'],
['Ar', 50, 0.51, 'cyan'],
['K' , 55, 0.51, 'silver'],
['Ca', 55, 0.51, 'slategray'],
['Sc', 55, 0.51, 'gray'],
['Ti', 55, 0.51, 'gray'],
['V' , 55, 0.51, 'gray'],
['Cr', 55, 0.51, 'gray'],
['Mn', 55, 0.51, 'gray'],
['Fe', 55, 0.51, 'gray'],
['Co', 55, 0.51, 'gray'],
['Ni', 55, 0.51, 'gray'],
['Cu', 55, 0.51, 'darkred'],
['Zn', 55, 0.51, 'gray'],
['Ga', 55, 0.51, 'gray'],
['Ge', 55, 0.51, 'gray'],
['As', 55, 0.51, 'rosybrown'],
['Se', 55, 0.51, 'dimgray'],
['Br', 55, 0.51, 'maroon'],
['Kr', 55, 0.51, 'steelblue'],
['Rb', 60, 0.51, 'silver'],
['Sr', 60, 0.51, 'slategray'],
['Y' , 60, 0.51, 'gray'],
['Zr', 60, 0.51, 'gray'],
['Nb', 60, 0.51, 'gray'],
['Mo', 60, 0.51, 'gray'],
['Tc', 60, 0.51, 'gray'],
['Ru', 60, 0.51, 'gray'],
['Rh', 60, 0.51, 'gray'],
['Pd', 60, 0.51, 'gray'],
['Ag', 60, 0.51, 'silver'],
['Cd', 60, 0.51, 'gray'],
['In', 60, 0.51, 'gray'],
['Sn', 60, 0.51, 'gray'],
['Sb', 60, 0.51, 'lightlategray'],
['Te', 60, 0.51, 'darkslategray'],
['I' , 60, 0.51, 'purple'],
['Xe', 60, 0.51, 'dodgerblue'],
['Cs', 65, 0.51, 'silver'],
['Ba', 65, 0.51, 'slategray'],
['La', 65, 0.51, 'gray'],
['Ce', 65, 0.51, 'gray'],
['Pr', 65, 0.51, 'gray'],
['Nd', 65, 0.51, 'gray'],
['Pm', 65, 0.51, 'gray'],
['Sm', 65, 0.51, 'gray'],
['Eu', 65, 0.51, 'gray'],
['Gd', 65, 0.51, 'gray'],
['Tb', 65, 0.51, 'gray'],
['Dy', 65, 0.51, 'gray'],
['Ho', 65, 0.51, 'gray'],
['Er', 65, 0.51, 'gray'],
['Tm', 65, 0.51, 'gray'],
['Yb', 65, 0.51, 'gray'],
['Lu', 65, 0.51, 'gray'],
['Hf', 65, 0.51, 'gray'],
['Ta', 65, 0.51, 'gray'],
['W' , 65, 0.51, 'gray'],
['Re', 65, 0.51, 'gray'],
['Os', 65, 0.51, 'gray'],
['Ir', 65, 0.51, 'gray'],
['Pt', 65, 0.51, 'gray'],
['Au', 65, 0.51, 'gold'],
['Hg', 65, 0.51, 'gray'],
['Tl', 65, 0.51, 'gray'],
['Pb', 65, 0.51, 'gray'],
['Bi', 65, 0.51, 'gray'],
['Po', 65, 0.51, 'gray'],
['At', 65, 0.51, 'gray'],
['Rn', 65, 0.51, 'aqua'],
['Fr', 70, 0.51, 'silver'],
['Ra', 70, 0.51, 'slategray'],
['Ac', 70, 0.51, 'gray'],
['Th', 70, 0.51, 'gray'],
['Pa', 70, 0.51, 'gray'],
['U' , 70, 0.51, 'gray'],
['Np', 70, 0.51, 'gray'],
['Pu', 70, 0.51, 'gray'],
['Am', 70, 0.51, 'gray'],
['Cm', 70, 0.51, 'gray'],
['Bk', 70, 0.51, 'gray'],
['Cf', 70, 0.51, 'gray'],
['Es', 70, 0.51, 'gray'],
['Fm', 70, 0.51, 'gray'],
['Md', 70, 0.51, 'gray'],
['No', 70, 0.51, 'gray'],
['Lr', 70, 0.51, 'gray'],
['Rf', 70, 0.51, 'gray'],
['Db', 70, 0.51, 'gray'],
['Sg', 70, 0.51, 'gray'],
['Bh', 70, 0.51, 'gray'],
['Hs', 70, 0.51, 'gray'],
['Mt', 70, 0.51, 'gray'],
['Ds', 70, 0.51, 'gray'],
['Rg', 70, 0.51, 'gray'],
['Cn', 70, 0.51, 'gray'],
['Nh', 70, 0.51, 'gray'],
['Fl', 70, 0.51, 'gray'],
['Mc', 70, 0.51, 'gray'],
['Lv', 70, 0.51, 'gray'],
['Ts', 70, 0.51, 'gray'],
['Og', 70, 0.51, 'gray'],
]

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

def procsgwinpsections(fnm,orbindx):
    "Parse 'sgwinp.txt' for sections"
    f = open(fnm,'r')
    sections=['$GEOMETRY',
              '$GRID',
              '$ORBITALS'
             ]
    geombuf=[]
    gridparbuf=[]
    ampbuf=[]
    kk=0
    foundorb=0
    lastsectread=None
    for line in f.readlines():
        kk=kk+1
        spl=line.split()
        foundsection=False
        for ii in range(len(sections)):
            if spl[0] == sections[ii]:
               lastsectread=ii
               foundsection=True
               break
        if foundsection: continue
        elif lastsectread == 0: geombuf.append(line)
        elif lastsectread == 1: gridparbuf.append(line)
        elif lastsectread == 2:
           if spl[0] == 'ORB' or spl[0] == 'DENS': # check for queried orbital
              if int(spl[1]) == orbindx: foundorb=kk
           elif kk-1 == foundorb: ampbuf.append(line)
    f.close()
    if foundorb == 0: sys.exit("ERROR: orbital %d not found in 'sgwinp.txt'" %(orbindx))
    w = extractgridamplitudes(ampbuf,0)
    npts,incs = extractgridparameters(gridparbuf)
    atoms,coords = extractatoms(geombuf)
    return w,npts,incs,atoms,coords

def extractgridparameters(buf):
    "Extracts grid parameters from buffer"
    shft=0  # line where grid data begins
    cshft=1 # data column with n-points, increments
    npts=[]
    incs=[]
    for ii in range(3):
        pt=buf[shft+ii].split()
        sp=buf[shft+3+ii].split()
        npts.append(int(pt[cshft]))
        incs.append(float(sp[cshft]))
    return npts,incs

def extractgridparametersALT(fnm):
    "Extracts grid parameters from file, assuming they come first"
    f = open(fnm,'r')
    indx=0
    buf = []
    for line in f.readlines():
       buf.append(line)
       indx += 1
       if indx == 6: break
    f.close()
    npts=[]
    incs=[]
    cshft=1 # data column with n-points, increments
    for ii in range(3):
        pt=buf[ii].split()
        sp=buf[3+ii].split()
        npts.append(int(pt[cshft]))
        incs.append(float(sp[cshft]))
    return npts,incs

def extractgridamplitudes(buf,ampline):
    "Extracts amplitudes from buffer"
    amptmp = buf[ampline].split()
    amps = [float(amptmp[ii]) for ii in range(len(amptmp))]
    return amps

def extractgridamplitudesALT(fnm,ampline):
    "Extracts grid amplitudes from file (linecache method)"
    amptmp = linecache.getline(fnm,ampline).split()
    amps = [float(amptmp[ii]) for ii in range(len(amptmp))]
    return amps

def generategrid(npts,incs):
    "Generate xyz values for grid"
    xs=[]
    ys=[]
    zs=[]
    xshft = npts[0]*incs[0]/2
    yshft = npts[1]*incs[1]/2
    zshft = npts[2]*incs[2]/2
    for kk in range(npts[2]):
        zval=incs[2]*kk-zshft
        for jj in range(npts[1]):
            yval=incs[1]*jj-yshft
            for ii in range(npts[0]):
                xval=incs[0]*ii-xshft
                xs.append(xval)
                ys.append(yval)
                zs.append(zval)
    return xs,ys,zs

def getgridlimits(npts,incs):
    "Computes limits of grid"
    lims=[[-npts[ii]*incs[ii]/2,npts[ii]*incs[ii]/2] for ii in range(3)]
    return lims

def filtergriddata(x,y,z,w,iso):
    "Filter out small values"
    xf=[]
    yf=[]
    zf=[]
    wf=[]
    for ii in range(len(w)):
        if abs(w[ii]) > iso:
           xf.append(x[ii])
           yf.append(y[ii])
           zf.append(z[ii])
           wf.append(w[ii])
    return xf,yf,zf,wf

def extractatoms(buf):
    "Gets atomic positions from buffer"
    atoms=[]
    coords=[]
    for ii in range(len(buf)):
        ltmp = buf[ii].split()
        xyz = [float(ltmp[jj+1]) for jj in range(3)]
        atoms.append(ltmp[0])
        coords.append(xyz)
    return atoms,coords

def drawheatmap(x,y,z,w,lims,atoms,coords,fnm):
   "Generates the plot via MatPlotLib"
   xlims=lims[0]
   ylims=lims[1]
   zlims=lims[2]
   # Plot points
   fig = plt.figure(figsize = (10, 7))
   ax = plt.axes(projection ="3d")
   ax.set_xlim(xlims)
   ax.set_ylim(ylims)
   ax.set_zlim(zlims)
   ax.set_xlabel("x")
   ax.set_ylabel("y")
   ax.set_zlabel("z")
   if len(w) > 0:
      vmax=max(w)
      vmin=min(w)
      vvmax=max(abs(vmax),abs(vmin))
      vvmin=-vvmax
      vdif=vvmax-vvmin
      vscl = [(w[ii] - vvmin)/vdif for ii in range(len(w))]
      ax.scatter3D(x, y, z,
                   s=10,
                   alpha=0.1,
                   c=cm.seismic(vscl))
   if coords != None:
      coordelem = makeelementarray(atoms,coords)
      for jj in range(len(elemtable)):
          if len(coordelem[jj]) > 0:
             xc = [coordelem[jj][ii][0] for ii in range(len(coordelem[jj]))]
             yc = [coordelem[jj][ii][1] for ii in range(len(coordelem[jj]))]
             zc = [coordelem[jj][ii][2] for ii in range(len(coordelem[jj]))]
             ax.scatter3D(xc, yc, zc, 
                          s=elemtable[jj][1],
                          alpha=elemtable[jj][2],
                          color=elemtable[jj][3])
   plt.title(fnm)
   plt.show()
   return

def makeelementarray(atoms,coords):
    "Makes arrays for each element for MatPlotLib"
    coordelem=[[] for jj in range(len(elemtable))]
    for ii in range(len(atoms)):
        for jj in range(len(elemtable)):
            if atoms[ii] == elemtable[jj][0]:
               coordelem[jj].append(coords[ii])
               break
    return coordelem

def makethegraph(wfnm,wfid,cfnm):
    "Main routine for making the graph from input file"
    iso=0.06
    if parsefname(wfnm) == 'sgwinp.txt':
       w,npts,incs,atoms,coords=procsgwinpsections(wfnm,wfid)
    else:
       if parsefname(wfnm) == 'wf.txt': w = extractgridamplitudesALT(wfnm,11+2*wfid)
       elif parsefname(wfnm) == 'orbj.txt': w = extractgridamplitudesALT(wfnm,8+wfid)
       else: w = extractgridamplitudesALT(wfnm,9)
       npts,incs = extractgridparametersALT(wfnm)
       if cfnm != None:
          cbuf = file2buffer(cfnm)
          atoms,coords = extractatoms(cbuf)
       else:
           atoms= None
           coords = None
    x,y,z = generategrid(npts,incs)
    xf,yf,zf,wf = filtergriddata(x,y,z,w,iso)
    lims = getgridlimits(npts,incs)
    drawheatmap(xf,yf,zf,wf,lims,atoms,coords,wfnm)
    return

if len(sys.argv) == 2: # single orbital or density
   wfnm = sys.argv[1]
   wfid = None
   cfnm = None
   makethegraph(wfnm,wfid,cfnm)
elif len(sys.argv) == 3:
   wfnm = sys.argv[1]
   if parsefname(wfnm) == 'wf.txt' or parsefname(wfnm) == 'orbj.txt' or parsefname(wfnm) ==  'sgwinp.txt': # selected orbital
      if sys.argv[2].isdigit():
         wfid = int(sys.argv[2])
         cfnm = None
         makethegraph(wfnm,wfid,cfnm)
      else: print("syntax: '$> %s %s <wf index>'" %(sys.argv[0],sys.argv[1]))
   else: # single orbital or density and coordinates
      if not sys.argv[2].isdigit():
         wfid = None
         cfnm = sys.argv[2]
         makethegraph(wfnm,wfid,cfnm)
      else: print("syntax: '$> %s %s <coordinate file name>'" %(sys.argv[0],sys.argv[1]))
elif len(sys.argv) == 4: # selected orbital and coordinates
   wfnm = sys.argv[1]
   if parsefname(wfnm) == 'wf.txt' or parsefname(wfnm) == 'orbj.txt':
      wfid = int(sys.argv[2])
      cfnm = sys.argv[3]
      makethegraph(wfnm,wfid,cfnm)
   else: print("syntax: '$> %s %s <coordinate file name>'" %(sys.argv[0],sys.argv[1]))
else:
    print("syntax: '$> %s <wf file name> <(optional) wf index)> <(optional) coordinate file name>'" %(sys.argv[0]))

