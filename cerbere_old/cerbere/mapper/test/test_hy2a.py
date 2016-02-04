"""
@author: Antoine Grouazel
@date: 18 may 2015
test: ok
"""
from cerbere.mapper.hy2ansoashdffile import HY2ANSOASHDFFile
from cerbere.datamodel.swath import Swath
from matplotlib import pyplot as plt
from numpy import amax,amin,arange,ravel
from mpl_toolkits.basemap import Basemap
ff = '/home/cerdata/provider/nsoas/satellite/l2/hy2/scatterometer/data/2015/134/H2A_SM2B20150514_18221.h5'
ff = '/home/cerdata/provider/nsoas/satellite/l2/hy2/scatterometer/data/2015/134/H2A_SM2B20150514_18218.h5'
# ff = '/home/cerdata/provider/nsoas/satellite/l2/hy2/scatterometer/data/2015/175/H2A_SM2B20150624_18776.h5'
handler = HY2ANSOASHDFFile(url=ff,mode='r')
handler.open()
my_swath = Swath()
my_swath.load(handler)

print 'mapper fieldnames',handler.get_fieldnames()
ws_selec = my_swath.get_values('wind_speed_selection')
wd_selec = my_swath.get_values('wind_dir_selection')
print 'min',amin(ws_selec),'max',amax(ws_selec),'shape',ws_selec.shape,'type',type(ws_selec)
lon = my_swath.get_lon()
lat = my_swath.get_lat()
print 'lon shape',lon.shape
print my_swath.get_bbox()
print 'feature',my_swath
print 'selction',my_swath.get_values('wind_dir_selection')
print 'solution1',my_swath.get_values('wind_dir_solution_1')
sol1 = handler.read_values('wind_dir_solution_1')
sol2 = handler.read_values('wind_dir_solution_2')
print 'mapper solution1',handler.read_values('wind_dir_solution_1')
print 'is equal sol1 and sol2',(sol1==sol2).all()
print 'time',handler.read_values('time')
print 'start_time',handler.get_start_time()
print 'stop_time',handler.get_end_time()

handler.close()

#plot
fig = plt.figure()
delta=5
lonmi=amin(lon)
lonma=amax(lon)
latmi=amin(lat)
latma=amax(lat)
mW=lonmi-delta
mE=lonma+delta
mN=latma+delta
mS=latmi-delta
mW=-179
mE=179
mN=80
mS=-80
print 'west=',mW, 'east=',mE,'north=',mN,'south=',mS
plt.hold('true')
m = Basemap(projection='merc',llcrnrlat=mS,urcrnrlat=mN,\
                llcrnrlon=mW, urcrnrlon=mE,resolution='c')
m.drawmapboundary(fill_color='#85A6D9')
m.drawcoastlines(color='#6D5F47', linewidth=.4)
m.fillcontinents(color='white',lake_color='#85A6D9')
m.drawcountries()
#m.shadedrelief()
meridians = arange(mW,mE,30)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10,fmt='%i')
parallels = arange(mS,mN,40)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10,fmt='%i')
filout = '/tmp/hy2a_test.png'
# plt.pcolor(ws_selec)
x,y = m(ravel(lon),ravel(lat))
m.scatter(x,y,c=wd_selec,s=0.1,lw=0.001,alpha=0.5)
plt.colorbar()
plt.savefig(filout,dpi=600)

print filout