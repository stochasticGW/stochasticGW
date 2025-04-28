# Python script for shifting xyz coordinates in .xyz file
# (useful for preparing shifted coordinates for preliminary DFT calculation)
import sys

def file2buffer(fnm):
    "Reads file into buffer"
    f = open(fnm,'r')
    buf = []
    for line in f.readlines():
        buf.append(line)
    f.close()
    return buf

def getcoordvals(buf):
    "Extracts coordinate values from file"
    buftmp=buf[0].split()
    nvals = int(buf[0])
#    print('nvals = %d' %(nvals))
    elem=[]
    xyz_raw=[]
    for ii in range(nvals):
        buftmp=buf[2+ii].split()
        if (len(buftmp) == 4):
           elem.append(buftmp[0])
           xyz = [float(buftmp[jj+1]) for jj in range(3)]
           xyz_raw.append(xyz)
        else: print("bad input line: %d" %(2+ii))
    return elem,xyz_raw

def shiftcoordvals(xyz_raw,dxyz):
    "Shifts coordinate values"
    xyz_shifted=[]
    for ii in range(len(xyz_raw)):
        xyzs = [xyz_raw[ii][jj]+dxyz[jj] for jj in range(len(dxyz))]
        xyz_shifted.append(xyzs)
    return xyz_shifted

def printshifteddata(buf,elem,xyz_shifted):
    "Print shifted .xyz data"
    print("%s" %(buf[0].rstrip('\n')))
    print("%s" %(buf[1].rstrip('\n')))
    for ii in range(len(elem)):
        print("%s  %16.10f  %16.10f %16.10f" %(elem[ii],xyz_shifted[ii][0],xyz_shifted[ii][1],xyz_shifted[ii][2]))
    return

def shiftcoordsmain(fnm,dxyz):
   "Main .xyz shifting script"
   buf=file2buffer(fnm)
   elem,xyz_raw = getcoordvals(buf)
   xyz_shifted = shiftcoordvals(xyz_raw,dxyz)
   printshifteddata(buf,elem,xyz_shifted)
   return

if len(sys.argv) == 5:
    fnm=sys.argv[1]
    dxyz=[float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4])]
    shiftcoordsmain(fnm,dxyz)
else:
    print("syntax: '$> %s <.xyz file name> <x-shift> <y-shift> <z-shift>'" %(sys.argv[0]))

