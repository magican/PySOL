from cerbere.mapper.gribfile import GribFile
from cerbere.datamodel.grid import Grid
import numpy
import pygrib
from matplotlib import pyplot as plt
import pdb
def test_attr(hh):
    attrs = hh.read_field_attributes('shtfl')
    for aa in attrs.keys():
        print aa,'=',attrs[aa]
    return

def test_lonlat(hh):
    print hh.get_fieldnames()
    lats = hh.read_values('lat')
    print 'lat:',lats.shape,lats
    lons = hh.read_values('lon')
    print 'lon:',lons
    lon_raw = hh._messages['shtfl'][0].longitudes
    # lon_raw = numpy.reshape(lon_raw,lats.shape)
    pcolor_lon = False
    if pcolor_lon:
        fig = plt.figure()
        fig.add_subplot(121)
        plt.pcolor(lon_raw,vmin=-180, vmax=360)
        plt.colorbar()
        fig.add_subplot(122)
        plt.pcolor(lons,vmin=-180, vmax=360)
        filout = '/tmp/lon_raw_cfsr.png'
        plt.savefig(filout)
        print filout
    pcolor_lat = False
    if pcolor_lat:
        lat_raw = hh._messages['shtfl'][0].latitudes
        lat_raw = numpy.reshape(lat_raw,lats.shape)
        fig = plt.figure()
        fig.add_subplot(121)
        plt.pcolor(lat_raw,vmin=-90, vmax=90)
        fig.add_subplot(122)
        plt.pcolor(lats,vmin=-90, vmax=90)
        plt.colorbar()
        filout = '/tmp/lat_raw_cfsr.png'
        plt.savefig(filout)
        print filout
    return
# pdb.set_trace()

# print 'lon1-lon0',lons[1] - lons[0]
# diff = numpy.diff(lons)
# zarb = list(diff==diff[50])
# fig = plt.figure(figsize=(24,24))
# plt.plot(diff[1:],'b.',markersize=2)
# new_diff = diff
# print zarb
# new_diff = numpy.ma.masked_where(zarb,diff,copy=False)
# plt.plot(new_diff,'ro',markersize=2)
# # filout = '/tmp/toto.png'
# plt.savefig(filout,dpi=800)
# print filout
# print 'time:',hh.read_values('time').shape
# sensinle_heat = hh.read_values('shtfl')

def recup_var_names_with_pygrib(ff):
    grbs=pygrib.open(ff)
    print dir(grbs)
    grbs.message(6).data()[1].shape
    for dodo in range(grbs.messages):
        print grbs.message(dodo+1).name
#     pdb.set_trace()
    grb = grbs.select(name='Momentum flux, u component')[0]
    data=grb.values
    lat,lon = grb.latlons()
    grbs.close()
    print data.shape
    return grbs.messages

def recup_var_name_with_mapper(hh):
    fields = hh.get_fieldnames()
    print fields
#     for fifi in fields:
#         fofo = hh.read_field(fifi)
#         attrs = hh.read_field_attributes(fifi)
#         for aat in attrs.keys():
# #             if 'abb' in aat or 'ariabl' in aat:
#             print fifi,aat,'=',attrs[aat] #shortNameECMF,parameterName,nameECMF,unitsECMF,units,parameterUnits
#         pdb.set_trace()
#         print fofo
    return len(fields)

def testgrid(hh):
    att_show = False
    if att_show==True:
        attribs = hh._messages['shtfl'][0].keys()
        # print 'attrib field',attribs
        for ii,k in enumerate(attribs):
            if not k in [
                    'values',
                    'distinctLongitudes',
                    'distinctLatitudes',
                    'latitudes',
                    'longitudes',
                    'latLonValues',
                    'codedValues'
                    ]:
        #         attrs[k] = hh._messages['shtfl'][0][k]
                try:
                    val = hh._messages['shtfl'][0][k]
                    print '       att',ii+1,'/',len(attribs),' ',k,':',val
                except:
                    print 'no values for key',k
                
    # lon_v1 = hh._messages['shtfl'][0]['longitudes']
    # lon_v2 = hh._messages['shtfl'][0]['distinctLongitudes']
    # lon_v3 = hh._messages['shtfl'][0]['latLonValues']
    # print '1-2',lon_v1==lon_v2
    # print '1-3',lon_v1==lon_v3
    # print '2-3',lon_v2==lon_v3
    # print '1',lon_v1
    # print '2',lon_v2
    # print '3',lon_v3
    # pdb.set_trace()
    gg = Grid()
    gg.load(hh)
    print "\PLOT"
    # from cerbereutils.plot.mapping import CerMap, SCATTER
    # m = CerMap(gg, 'shtfl', contouring=SCATTER)
    from mpl_toolkits.basemap import Basemap
    fig = plt.figure()
    limS=-80
    limN=85
    limW=-179
    limE=179
    #     projection = Basemap(projection='ortho',lon_0=-110,lat_0=-30,resolution='l')
    projection = Basemap(projection='merc',llcrnrlat=limS,urcrnrlat=limN,llcrnrlon=limW, urcrnrlon=limE,resolution='i')
    
    # projection = Basemap(projection='npstere',boundinglat=10,lon_0=270,resolution='l')
    projection.drawcoastlines(linewidth=0.25)
    scatter = False
    if scatter:
        lons = numpy.ravel(lons)
        lats = numpy.ravel(lats)
    else:
        lons,lats = numpy.meshgrid(lons,lats)
    x,y = projection(lons,lats)
    pdb.set_trace()
    if scatter:
        param = numpy.ravel(gg.get_values('shtfl'))
        projection.scatter(x,y, s=0.1, c=param,alpha=0.7,lw=0.001)
    else:
        param = gg.get_values('u')
        projection.contourf(x,y,param,50)
    plt.colorbar()
    filout = '/tmp/toto.png'
    plt.savefig(filout,dpi=800)
    print 'file output:',filout
    return
if __name__ == '__main__':
#     ff = '/home/cercache/project/oceanheatflux/data/sources/fluxes/cfsr/1990/flxf06.gdas.1990123118.grb2'
    ff = "/home/cercache/users/mirror/home/oo7/oo/ncep/cfsr/grb/2002/flxf06.gdas.2002081618.grb2"
    hh = GribFile(url=ff,mode='r')
    hh.open()
    taille_vraie = recup_var_names_with_pygrib(ff)
    taille_proxy = recup_var_name_with_mapper(hh)
    print 'taille_vraie',taille_vraie
    print 'taille_proxy',taille_proxy
    test_attr(hh)
#     
    hh.close()