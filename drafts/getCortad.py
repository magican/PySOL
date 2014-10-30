import os, sys, datetime, string
import numpy as np
from netCDF4 import Dataset
import numpy.ma as ma
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
from pylab import *
#import laplaceFilter
import mpl_util

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2011, 6, 25)
__modified__ = datetime.datetime(2012, 4, 10)
__version__  = "1.0"
__status__   = "Development, 25.6.2011, 10.04.2012"

def makeMap(latStart,latEnd,lonStart,lonEnd,longitude,latitude,filledSST,name,dateString):

    plt.figure(figsize=(12,10),frameon=False)

    if lonStart< 0 and lonEnd < 0:
        lon_0= - (abs(lonEnd)+abs(lonStart))/2.0
    elif lonStart > 0 and lonEnd > 0:
        lon_0=(abs(lonEnd)+abs(lonStart))/2.0
    else:
        lon_0=((lonEnd)+(lonStart))/2.0

    map = Basemap(llcrnrlat=latStart,urcrnrlat=latEnd,\
            llcrnrlon=lonStart,urcrnrlon=lonEnd,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=latStart,lon_0=lon_0)

    lon,lat=meshgrid(longitude,latitude)

    x, y = map(lon,lat)
    map.drawcoastlines()
    map.fillcontinents(color='grey')
    map.drawcountries()
    #map.bluemarble()

    levels = np.arange(-2.0, 20.0, 0.5)

    CS1 = map.contourf(x,y,filledSST,levels,
                        cmap=cm.get_cmap('jet',len(levels)-1),
                        extend='both',
                        alpha=1.0,
                        origin='lower',
                        rasterized=True)


    CS1.axis='tight'
    map.drawmeridians(np.arange(lon.min(),lon.max(),10),labels=[0,0,0,1])
    map.drawparallels(np.arange(lat.min(),lat.max(),4),labels=[1,0,0,0])


    plt.colorbar(CS1,orientation='vertical',extend='both', shrink=0.5)

    plt.title('CoRTAD FilledSST - %s'%(dateString))
    plotfile='/home/mag/SST_'+str(dateString)+'_'+str(name)+'.png'
    plt.savefig(plotfile)
    plt.show()



def extratCORTAD(name,latStart,latEnd,lonStart,lonEnd):

    """
    Info on the different tiles used to identoify a region is found here:
    http://www.nodc.noaa.gov/SatelliteData/Cortad/TileMap.jpg
    """
    base="http://data.nodc.noaa.gov/opendap/cortad/Version4/"
    base="http://data.nodc.noaa.gov/opendap/cortad/Version3/"

    refDate = datetime.datetime(1982,1,1,0,0,0)

    if name=="North Sea":
        file1="cortadv3_row01_col07.h5"
        file2="cortadv3_row01_col08.h5"
        tiles=2
        direction=1 # add the tiles
    if name=="Lion":
#        file1="cortadv4_row02_col07.nc"
#        file2="cortadv4_row02_col08.nc"
        file1="cortadv3_row02_col07.h5"
        file2="cortadv3_row02_col08.h5"
        tiles=2
        direction=1 # add the tiles

    filename1=base+file1
    filename2=base+file2
    print "Extrating data from openDAP: %s"%(filename1)
    print "Extrating data from openDAP: %s"%(filename2)

    cdf1=Dataset(filename1)
    cdf2=Dataset(filename2)

    longitude1=np.squeeze(cdf1.variables["Longitude"][:])
    latitude1=np.squeeze(cdf1.variables["Latitude"][:])

    longitude2=np.squeeze(cdf2.variables["Longitude"][:])
    latitude2=np.squeeze(cdf2.variables["Latitude"][:])
    time=np.squeeze(cdf2.variables["Time"][:])
    print "Finished extracting the grid"

    wantedYears=[2010]
#    wantedMonths=[1,2,3,4,5,6,7,8,9,10,11,12]
    wantedMonths=[12]


    for t in range(len(time)):
        currentTime=refDate + datetime.timedelta(days=time[t])
        #print currentTime
        if currentTime.year in wantedYears and currentTime.month in wantedMonths:
            print "Current time inspecting: %s"%(currentTime)
            myIndex=t; dateString="%s"%(currentTime)

            filledSST1=cdf1.variables["FilledSST"][:,:,myIndex]
            filledSST2=cdf2.variables["FilledSST"][:,:,myIndex]


            offset1=cdf1.variables["FilledSST"].__getattribute__('add_offset')
            offset2=cdf2.variables["FilledSST"].__getattribute__('add_offset')

            filledSST1=filledSST1 + abs(offset1)
            filledSST1=np.rot90(filledSST1,3) # rotate 90 degrees 3 times
            filledSST1=np.fliplr(filledSST1) # then flip left right

            filledSST2=filledSST2 + abs(offset2)
            filledSST2=np.rot90(filledSST2,3) # rotate 90 degrees 3 times
            filledSST2=np.fliplr(filledSST2) # then flip left right

            filledSST=concatenate((filledSST1,filledSST2),axis=direction)

            if direction==1:
                longitude=concatenate((longitude1,longitude2),axis=direction)
                latitude=latitude1
            else:
                longitude=longitude1
                latitude=concatenate((latitude1,latitude2),axis=direction)

            print "Created matrix for SST for area: %s"%(name)
            print "Matrix size: ",filledSST.shape
            print "Longitude size: ",longitude.shape
            print "Latitude size: ",latitude.shape
            print "Long min: %s Long max: %s"%(longitude.min(),longitude.max())
            print "Lat min: %s Lat max: %s"%(latitude.min(),latitude.max())
            print "------------------------------\n"

            makeMap(latStart,latEnd,lonStart,lonEnd,longitude,latitude,filledSST,name,dateString)


def main():

    latStart = 41
    latEnd   = 43
    lonStart = 2
    lonEnd   = 4

    name='Lion'
    name='North Sea'

    extratCORTAD(name,latStart,latEnd,lonStart,lonEnd)




if __name__ == "__main__":

    main()
