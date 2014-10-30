#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 18:05:58 2012

@author: mag
"""

import numpy as np
import pyresample as pr

filename = "/home/mag/Documents/repos/solab/PySOL/drafts/ssmis_swath.npz"
data = np.load(filename)['data']
lons = data[:, 0].astype(np.float64)
lats = data[:, 1].astype(np.float64)
tb37v = data[:, 2].astype(np.float64)
    
area_def = pr.utils.parse_area_file('/home/mag/Documents/repos/solab/PySOL/drafts/areas.cfg', 'ease_sh')[0]

swath_def = pr.geometry.SwathDefinition(lons, lats)

result = pr.kd_tree.resample_nearest(swath_def, tb37v, area_def, \
    radius_of_influence=20000, fill_value=None)

bmap = pr.plot.area_def2basemap(area_def)
col = bmap.imshow(result, origin='upper')
bmap.drawcoastlines()
bmap.drawmeridians(arange(-180,180,15),labels=[0,0,1,0],fontsize=12, color='k')
bmap.drawparallels(arange(-90,90,15),labels=[0,0,0,0],fontsize=12, color='k')

bmng = bmap.bluemarble()

pr.plot.save_quicklook('tb37v_quick.png', area_def, result, label='Tb 37v (K)')




import numpy as np
from pyresample import image, geometry
area_def = geometry.AreaDefinition('areaD', 'Europe (3km, HRV, VTC)', 'areaD',
                               {'a': '6378144.0', 'b': '6356759.0',
                                'lat_0': '50.00', 'lat_ts': '50.00',
                                'lon_0': '8.00', 'proj': 'stere'},
                               5924, 7930,
                               [-1370912.72, -909968.64,
                                1029087.28, 1490031.36])
msg_area = geometry.AreaDefinition('msg_full', 'Full globe MSG image 0 degrees',
                               'msg_full',
                               {'a': '6378169.0', 'b': '6356584.0',
                                'h': '35785831.0', 'lon_0': '0',
                                'proj': 'geos'},
                               5924, 7930,
                               [-5568742.4, -5568742.4,
                                5568742.4, 5568742.4])
                                
data = S_VV_ABS

msg_con_quick = image.ImageContainerQuick(data, msg_area)

area_con_quick = msg_con_quick.resample(area_def)

result_data_quick = area_con_quick.image_data

msg_con_nn = image.ImageContainerNearest(data, msg_area, radius_of_influence=50000)

area_con_nn = msg_con_nn.resample(area_def)

result_data_nn = area_con_nn.image_data

result = pr.kd_tree.resample_nearest(msg_area, data.ravel(), \
    area_def, radius_of_influence=50000, epsilon=100, nprocs=4)