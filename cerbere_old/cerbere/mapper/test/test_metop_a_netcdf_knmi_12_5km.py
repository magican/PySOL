from cerbere.mapper.ncfile import NCFile
from cerbere.mapper.knmil2ifrncfile import KNMIL2IFRNCFile#pour lire fichier de denis
from cerbere.mapper.knmil2ncfile import KNMIL2NCFile#pour lire fichier netcdf du knmi
from cerbere.datamodel.swath import Swath
from cerbereutils.plot.mapping import CerMap, SCATTER
"""validated 02/06/2015 agrouaze"""
ff = '/home/cerdata/provider/knmi/satellite/metop-a/ascat/l2b/12.5km/2015/154/ascat_20150603_001500_metopa_44728_eps_o_coa_2300_ovw.l2.nc'
hh = KNMIL2NCFile(ff,collection='ASCAT-A-L2-12_5km')
matching = {'cell':'NUMCELLS','row':'NUMROWS'}
# hh = NCFile(ff,geodim_matching=matching)
print hh.get_fieldnames()
print 'collection_id',hh.get_collection_id()
ss = Swath()
ss.load(hh)
m = CerMap(ss, 'wind_speed', contouring=SCATTER)
filout = '/tmp/toto.png'
m.save(filout)
print 'file output:',filout