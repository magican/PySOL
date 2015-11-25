import os
import sys
import numpy
import netCDF4
from cerbere.mapper.ncfile import NCFile
from matplotlib import pyplot as plt
if __name__ =='__main__':
    filou = "/home/cerdata/provider/neodc/l4/esacci_sst/2010/03/02/20100302120000-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_LT-v02.0-fv01.0.nc"
    print filou
    nc = NCFile(filou)
    filout = '/tmp/ostia_sst.png'
    filout2 = '/tmp/ostia_sst2.png'
    sst = nc.read_values('analysed_sst')
    sst = numpy.squeeze(sst)
    print sst.shape
    nc.close()
    
    nc = netCDF4.Dataset(filou)
    sst2 = nc.variables['analysed_sst'][:]
    sst2 = numpy.squeeze(sst2)
    print sst2.shape
    nc.close()
    
    #plot
    plt.figure()
    plt.pcolor(sst[200:800,100:600])
    plt.colorbar()
    plt.savefig(filout)
    #plot
    plt.figure()
    plt.pcolor(sst2[200:800,100:600])
    plt.colorbar()
    plt.savefig(filout2)
    print filout