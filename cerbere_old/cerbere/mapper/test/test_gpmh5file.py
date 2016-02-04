'''
Created on 4 mai 2015

@author: jfpiolle
'''

from cerbere.mapper import gpmhdf5file
from matplotlib import pyplot as plt
import pdb
import netCDF4
import numpy
input = '/home/cerdata/provider/jaxa/satellite/l2/gpm/kupr/2014/05/GPMCOR_KUR_1405201429_1602_001275_L2S_DU2_03B.h5'
# input = '/home/cerdata/provider/jaxa/satellite/l2/gpm/kupr/2015/03/GPMCOR_KUR_1503111732_1905_005867_L2S_DU2_03B.h5'
f = gpmhdf5file.GPMHDF5File(input)
# f = gpmhdf5file.GPMHDF5File('/home/cerdata/provider/jaxa/satellite/l2/gpm/kapr/2015/03/GPMCOR_KAR_1503192259_0032_005995_L2S_DA2_03B.h5', subproduct="HS")
testfield = '/NS/SLV/sigmaZeroCorrected'
# testfield = "/MS/SLV/precipRateNearSurface"
# testfield = "/HS/SLV/precipRateNearSurface"
# testfield = '/NS/SLV/precipRate'

f.open()

basic = True
youwantplot = False
fieldnames = f.get_fieldnames()
if basic==True:
    
    for fname in fieldnames:
        print fname
    
    print "\nGlobal Attributes"
    print "-------------------"
    attrs = f.read_global_attributes()
    for attr in attrs:
        print attr, f.read_global_attribute(attr)
    
    print "\Dimensions of field : ", testfield
    print "---------------------------------------------------"
    dims = f.get_dimensions(testfield)
    for dim in dims:
        print dim, f.get_dimsize(dim)
    
    dims = f.get_dimensions('lat')
    print dims
    dims = f.get_dimensions('lon')
    print dims
    
    print "\Attributes of field : ", testfield
    print "---------------------------------------------------"
    attrs = f.read_field_attributes(testfield)
    print attrs
    for att in attrs:
        print att, attrs[att]
    
    
    
    print "\nField ", testfield
    print "------------------------------------"
    field = f.read_field(testfield)
    print field
    
    print "\nRead 'lat'"
    print "------------------------------------"
    field = f.read_field('lat')
    print field
    print f.read_values('lat', slices=(slice(0, 10), slice(0, 10)))
    print f.read_values('lon', slices=(slice(0, 10), slice(0, 10)))
    
    print "\nTime coverage"
    print "---------------"
    print f.get_start_time()
    print f.get_end_time()
    
    print "\nRead Time"
    print "---------------"
    field_time = f.read_field('time')
    print field_time
    print f.read_values('time', slices=(slice(0, 10), slice(0, 10)))
    times = f.read_values('time')
    mimi = netCDF4.num2date(numpy.amin(times),field_time.units)
    mama = netCDF4.num2date(numpy.amax(times),field_time.units)
    print 'max time',mama,'min time',mimi
    print input
    print field_time.units
    print 'start attr',f.get_start_time()
    print 'stop attr',f.get_end_time()
    
    print "\nRead ", testfield
    print "---------------"
    print f.read_field(testfield)
    print f.read_values(testfield, slices=(slice(0, 10), slice(0, 10)))
    
if youwantplot==True:

    from cerbere.datamodel.swath import Swath
    swath = Swath()
    swath.load(f)
    
    
    # print "\PLOT"
    # from cerbereutils.plot.mapping import CerMap, SCATTER
    # from mpl_toolkits.basemap import Basemap
    # m = CerMap(swath, testfield, contouring=SCATTER)
    # filout = '/tmp/toto.png'
    # m.save(filout)
    # print 'file output:',filout
    
    print '3d plot'
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.basemap import Basemap
    
    limS=-80
    limN=80
    limW=-180
    limE=180
    projection = Basemap(projection='merc',llcrnrlat=limS,urcrnrlat=limN,llcrnrlon=limW, urcrnrlon=limE,resolution='i')
    
    fig = plt.figure(figsize=(15,15))
    ax = Axes3D(fig)
    
    
    ax.azim = 240
    ax.elev = 70
    ax.dist = 4
    
    
    ax.add_collection3d(projection.drawcoastlines(linewidth=0.25))
    ax.add_collection3d(projection.drawcountries(linewidth=0.35))
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    
    # for z in range(f.get_dimsize('nbin')):
    #     pdb.set_trace()
    elev = f.read_values('/NS/PRE/elevation')
    elev = numpy.ravel(elev)
    lon = numpy.ravel(swath.get_lon())
    lat = numpy.ravel(swath.get_lat())
    precip = numpy.ravel(swath.get_values('/NS/SLV/precipRate')[:,:,0])
    # for lo in swath.get_lon()[0]
    x, y = projection(lon, lat)
    ax.scatter(x,y, elev, c=elev,s=2,lw=0.01)
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    filout2 = '/tmp/scatter3d.png'
    plt.savefig(filout2,dpi=500)
    print filout2
    # ax.plot(x, y, z, label='parametric curve')
    # import mpld3
    # html_output = '/tmp/3dscatter.html'
    # mpld3.save_html(fig, html_output)
    # print html_output