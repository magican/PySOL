from cerbere.mapper.rapidscatpodaacncfile import RapidScatPODAACNCFile
import netCDF4

ff='/home/cerdata/provider/podaac/satellite/l2b/iss/rapidscat/netcdf/v1.1/2015/053/rs_l2b_v1.1_02369_201504161704.nc'
print ff
nc=RapidScatPODAACNCFile(ff)
tf=nc.read_field('time')
print tf
toto=nc.read_values('time')
tata = netCDF4.num2date(toto,tf.units)
print tata
print nc.get_start_time()
print nc.get_end_time()
