from cerbere.mapper.gribfile import GribFile
from cerbere.datamodel.grid import Grid

g = Grid()
gf = GribFile(url='/home1/taveeg/data/private/mpc-sentinel1/ancillary_data/ECMWF_0125/D1D07290000080800001')
g.load(gf)
