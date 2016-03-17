import os
import gc

import numpy as np
import numpy.ma as ma
import math
from netCDF4 import Dataset as ncDataset

__author__ = 'Denis Spiridonov'


def array_size_normalize(input_array, masked=False):
    avail_resolution_list = [int(256*math.pow(2, i)) for i in range(15)]
    # print avail_resolution_list
    shape = input_array.shape
    # print 'shape: ', shape
    current_pixels = max(shape)
    # print 'current pixels: ', current_pixels
    necessary_pixels = min(filter(lambda x: current_pixels < x,
                                  avail_resolution_list))
    # print 'necessary_pixels: ', necessary_pixels

    source_array = input_array[:]

    # src =
    #
    # xxx
    # xxx
    #
    # top_array =
    #
    # xxx
    # xxx
    #
    # res =
    #
    # xxx
    # xxx
    # xxx
    # xxx
    #
    # right_array =
    #
    # x
    # x
    # x
    # x
    #
    # res 4x4 =
    #
    # xxxx
    # xxxx
    # xxxx
    # xxxx
    if masked:
        top_array = np.zeros((necessary_pixels-shape[0], shape[1]))
        top_array = np.where(top_array == 0, 99999, top_array)
        right_array = np.zeros((necessary_pixels, necessary_pixels-shape[1]))
        right_array = np.where(right_array == 0, 99999, right_array)

        source_array = ma.concatenate((top_array, source_array), axis=0)
        # print source_array.shape
        source_array = ma.concatenate((source_array, right_array), axis=1)
        source_array = np.ma.masked_equal(source_array, 99999)
        # print source_array.shape
    if not masked:
        top_array = np.zeros((necessary_pixels-shape[0], shape[1]))
        top_array = np.where(top_array == 0, 99999, top_array)
        right_array = np.zeros((necessary_pixels, necessary_pixels-shape[1]))
        right_array = np.where(right_array == 0, 99999, right_array)

        source_array = np.concatenate((top_array, source_array), axis=0)
        # print source_array.shape
        source_array = np.concatenate((source_array, right_array), axis=1)
        # print source_array.shape
        source_array = np.ma.masked_equal(source_array, 99999)

    del top_array, right_array
    return source_array


def create_base_tiles(input_array):
    number_tiles = input_array.shape[0]/256

    lines = []
    all_list = []

    for i in range(number_tiles):
        for j in range(number_tiles):
            # print '[%d:%d' % (j*256, (j+1)*256), ',%d:%d]'%(i*256, (i+1)*256)
            lines.append(input_array[j*256:(j+1)*256, i*256:(i+1)*256])

        all_list.append(lines)
        lines = []
    return all_list


def decrease_zoom(source_array):
    rows, cols = source_array.shape
    rows_2 = rows/2
    cols_2 = cols/2
    sh = rows_2, rows//rows_2, cols_2, cols//cols_2
    return source_array.reshape(sh).mean(-1).mean(1)


def concatenate_and_decrease(tile00, tile01, tile10, tile11):
    tiles = np.concatenate((np.concatenate((tile00, tile01), axis=0),
                            np.concatenate((tile10, tile11), axis=0)),
                           axis=1)
    return decrease_zoom(tiles)


def write_tile_to_nc(nc_variable, tiles_array, variable, zoom, _min=0, _max=4,
                     polarization=None):
    # print len(tiles_array)
    for row_number in range(len(tiles_array)):
        tile_row = tiles_array[row_number]
        for col_number in range(len(tile_row)):
            m = tile_row[col_number]
            m = ma.where(m <= _min, _min + 0.0001, m)
            m = ma.where(m >= _max, _max, m)

            # 0-254 , 255 for mask
            m = (m-_min)/float(_max-_min) * (2**8-2)

            m = np.where(m.mask, (2**8-1), m)
            m = np.uint8(m)
            # m = np.ma.masked_where(m == 255, m)
            if polarization is None:
                nc_variable[variable, zoom, row_number, col_number, :] = m[:]
            else:
                nc_variable[variable, polarization, zoom,
                            row_number, col_number, :] = m[:]


def create_dataset(output_path, max_zoom, variables_list,
                   polarizations=None):
    if os.path.isfile(output_path):
        return ncDataset(output_path, 'a', format='NETCDF4')

    max_x_tiles = 2**max_zoom

    dataset = ncDataset(output_path, 'w', format='NETCDF4')
    dataset.createDimension('vars', len(variables_list))
    if polarizations is not None:
        dataset.createDimension('polarizations', 4)
    dataset.createDimension('zoom', max_zoom+1)
    dataset.createDimension('x', max_x_tiles)
    dataset.createDimension('y', max_x_tiles)
    dataset.createDimension('shape0', 256)
    dataset.createDimension('shape1', 256)

    if polarizations is not None:
        dims = ('vars', 'polarizations', 'zoom', 'x', 'y', 'shape0', 'shape1')
    else:
        dims = ('vars', 'zoom', 'x', 'y', 'shape0', 'shape1')

    vars_var = dataset.createVariable('Variables', 'S1', ('vars',),
                                      zlib=True, complevel=6)
    vars_var[:] = np.array(variables_list, dtype='string')[:]

    if polarizations is not None:
        polar_var = dataset.createVariable('Polarizations', 'S1',
                                           ('polarizations',),
                                           zlib=True, complevel=6)
        polar_var[:] = np.array(polarizations, dtype='string')[:]

    # u1 = NC_UBYTE 0-255
    dataset.createVariable('Data', 'u1', dims, zlib=True, complevel=6)
    return dataset


def create_nc_tiles(reprojected_array, variable, nc_path, polarization=None):
    var_list = ['roughness', 'wind_speed', 'sigma0']
    minmax = {  'roughness': [-1, 1],
                'wind_speed': [0, 35],
                'sigma0': [-35, 5]}
    p_mask = {'roughness': False,
              'wind_speed': True,
              'sigma0': False}

    polarization_list = ['hh', 'vv', 'hv', 'vh']

    if variable in var_list:
        variable_num = var_list.index(variable)
    else:
        print 'Variable not found in var_list'
        return False

    if polarization is not None:
        polarization_num = polarization_list.index(polarization)
    else:
        print 'Polarization not found in polarization_list'
        return False

    minmax_value = minmax[variable]
    # outDimVar[link_num, :] = varin[:]

    source_array = array_size_normalize(reprojected_array, p_mask[variable])


    number_tiles = source_array.shape[0]/256
    zoom_level = int(math.log(number_tiles, 2))
#    if not False in reprojected_array.mask:
#        return zoom_level
    print 'zoom_level: ', zoom_level

    dataset = create_dataset(nc_path, zoom_level, var_list, polarization_list)
    datasetVar = dataset.variables['Data']

    for zoom in range(zoom_level, -1, -1):
        print 'Start zoom level %d' % zoom
        # print source_array.shape

        tiles_array = create_base_tiles(source_array)

        write_tile_to_nc(datasetVar, tiles_array, variable_num, zoom,
                         _min=minmax_value[0], _max=minmax_value[1],
                         polarization=polarization_num)
        del tiles_array

        if zoom:
            source_array = decrease_zoom(source_array)
        print
        gc.collect()

    dataset.close()
    return zoom_level
