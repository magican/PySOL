'''test lecture of bufr file'''

#setenv PYTHONPATH ~/git/cerbere/:/opt/lib/python2.6/site-packages/
#setenv BUFR_TABLES /opt/lib/emos/bufrtables/
import time
import logging
from cerbere.mapper.knmiscatbufrfile import KNMIScatBUFRFile
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    start = time.time()
    #from cerbereutils.plot.mapping import Map
    # f = '/home/cerdata/provider/knmi/satellite/iss/rapidscat/l2b/bufr_knmi/25km/2015/101/rapid_20150411_151310_iss____03111_o_250_ovw_l2.bufr'
    f = '/home/cerdata/provider/knmi/satellite/iss/rapidscat/l2b/bufr_knmi/25km/2015/110/rapid_20150420_211249_iss____03255_o_250_ovw_l2.bufr'
    # f = '/home/cerdata/provider/knmi/satellite/iss/rapidscat/l2b/bufr_knmi/25km/2015/110/rapid_20150420_224518_iss____03256_o_250_ovw_l2.bufr'
    # f = '/tmp/ascat_20150531_013300_metopa_44686_eps_o_coa_ovw.l2_bufr'
    # f = '/tmp/ascat_20150601_062400_metopa_44703_eps_o_coa_ovw.l2_bufr'
    # f = '/tmp/ascat_20150601_135839_metopb_14021_ear_o_125_sva_ovw.l2_bufr'
    bufrfile = KNMIScatBUFRFile(url=f)
    bufrfile.open()
    
    print "\nFIELD NAMES"
    fieldnames = bufrfile.get_fieldnames()
    for fieldname in fieldnames:
        print 'var:',fieldname
    
    print "\nATTRIBUTES"
    attrs = bufrfile.read_global_attributes()
    # for att in attrs:
    #     print att, bufrfile.read_global_attribute(att)
    
    print "\nFIELDS"
    # for fieldname in fieldnames:
    #     print bufrfile.read_field(fieldname)
    
    print "\nSUBSAMPLING"
    #print bufrfile.read_values('lat', slices=(slice(0,10), slice(0,10)))
    #print bufrfile.read_values('time', slices=(slice(0,10), slice(0,10)))
    
    print "\nTIME COVERAGE"
    print bufrfile.get_start_time()
    print bufrfile.get_end_time()
    
    # LOAD THROUGH MODEL
    print "\n\nLOAD SWATH MODEL"
    from cerbere.datamodel.swath import Swath
    swath = Swath()
    swath.load(bufrfile)
    
    
    
    data = swath.get_lat()
    print "#####", data.min(), data.max(), data.count(), data
    data = swath.get_lon()
    print "#####", data.min(), data.max(), data.count(), data
    print swath.get_bbox()
    print 'lecture time elapsed',time.time()-start
    special_plot = True
    if special_plot:
        import numpy
        import pdb
        import matplotlib.pyplot as plt
    #     from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.basemap import Basemap
        
        limS=-89
        limN=89
        limW=-179
        limE=179
    #     projection = Basemap(projection='ortho',lon_0=-110,lat_0=-30,resolution='l')
    #     projection = Basemap(projection='merc',llcrnrlat=limS,urcrnrlat=limN,llcrnrlon=limW, urcrnrlon=limE,resolution='i')
        projection = Basemap(projection='npstere',boundinglat=10,lon_0=270,resolution='l')
        projection.drawcoastlines(linewidth=0.25)
    #     projection.drawcountries(linewidth=0.25)
    #     projection.drawmapboundary(fill_color='aqua')
    #     projection.fillcontinents(color='yellow',lake_color='aqua')
    #     fig = plt.figure()
    #     pdb.set_trace()
        lon0 = swath.get_lon()
        lat0 = swath.get_lat()
        lon = numpy.ma.ravel(lon0)
        lat = numpy.ma.ravel(lat0)
        
        # for lo in swath.get_lon()[0]
        
    #     print lon[0,0],lat[0,0]
        print 'lon max',numpy.ma.amax(lon),'lon min',numpy.ma.amin(lon),lon0.shape
        print 'lat max',numpy.ma.amax(lat),'lat min',numpy.ma.amin(lat),lat0.shape
        x,y = projection(lon,lat)
    #     for ww in range(4):
    #         ww_str = str(ww)
    #         winds0 = swath.get_values('wind_speed_at_10_m_'+ww_str)[:,:]
    #         winds = numpy.ma.ravel(winds0)
    #         winds.mask = False
    #         print 'wind_speed',ww_str,winds.count()
    #     winds0 = swath.get_values('wind_speed_at_10_m')
        winds0 = swath.get_values('model_wind_direction_at_10m')
        winds = numpy.ma.ravel(winds0)
        print winds
        print 'winds max',numpy.ma.amax(winds),'winds min',numpy.ma.amin(winds),winds0.shape,winds.count()
        
    #     pdb.set_trace()
    #     projection.plot(x,y,'r.',alpha=0.9)
        
        projection.scatter(x,y, s=0.1, c=winds,alpha=0.7,lw=0.001)
        plt.colorbar()
    #     ax = plt.gca()
    #     ax.set_xlabel('X Label')
    #     ax.set_ylabel('Y Label')
    #     ax.set_zlabel('Z Label')
        filout2 = '/tmp/knmi_bufr.png'
        plt.savefig(filout2,dpi=500)
        print filout2
    else:
        print "\PLOT"
        from cerbereutils.plot.mapping import CerMap, SCATTER
        m = CerMap(swath, 'model_wind_speed_at_10m', contouring=SCATTER)
        filout = '/tmp/toto.png'
        m.save(filout)
        print 'file output:',filout
        m = 0
    
    end = time.time()
    print 'time elapsed',end-start