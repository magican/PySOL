# -*- coding: utf-8 -*-
"""
"""


import re
from collections import OrderedDict
from dateutil import parser
import numpy as np
from netCDF4 import Dataset
from ast import literal_eval

from .ncfile import NCFile
from .. import READ_ONLY
from ..datamodel.field import Field
from ..datamodel.variable import Variable


# TODO:
# - read_field_attributes() : change long name when action=merge_var (beam->incidence)
# - not implemented functions to check (end of file)
# (get_bbox() : only for incidence ?)
# - add docstrings


def get_macrocycle_angles(url):
    """
    """
    dataset = Dataset(url)
    if 'macrocycle_angle' in dataset.ncattrs():
        angles = dataset.macrocycle_angle
        if isinstance(angles, np.ndarray):
            angles = angles.tolist()
        elif isinstance(angles, basestring) and ';' in angles:
            angles = angles.replace(';', ',')
            angles = literal_eval('[{}]'.format(angles))
        else:
            dataset.close()
            raise Exception('Unexpected macrocycle_angle attribute : {}'.format(angles))
    elif len(dataset.dimensions.get('n_beam', [])) == 6:
        angles = [0, 2, 4, 6, 8, 10]
    else:
        dataset.close()
        raise Exception('Could not find macrocycle angles.')
    dataset.close()
    return angles


class IncidenceError(Exception):
    pass


# The mapper virtually modifies the format in order to give
# dimensions/variables/attributes for a given incidence angle instead of
# a given beam (which is different for non nominal macrocycles).
# To be used with a Trajectory datamodel.
class SWIMCWWICL1NCFile(NCFile):
    """
    """
    def __init__(self, url, incidence, mode=READ_ONLY):
        """
        """
        if mode != READ_ONLY:
            raise Exception('Only mode {} is accepted.'.format(READ_ONLY))
        super(SWIMCWWICL1NCFile, self).__init__(url=url, mode=mode)
        self._incidence = incidence
        self._l1type = None
        self._mcycle = OrderedDict()
        self._vars = OrderedDict()
        self._vars2attr = OrderedDict()
        self._dims = OrderedDict()
        self._parse()
        return

    def _parse(self):
        """
        """
        # Infer L1A or L1B
        title = super(SWIMCWWICL1NCFile, self).read_global_attribute('title')
        if 'L1a' in title:
            self._l1type = 'L1A'
        elif 'L1b' in title:
            self._l1type = 'L1B'
        elif 'AUX_METEO' in title:
            self._l1type = 'AUXMETEO'
        else:
            raise Exception('Could not infer L1A or L1B.')
        # Get macrocycle
        mc_angle = super(SWIMCWWICL1NCFile, self).read_global_attribute('macrocycle_angle')
        if isinstance(mc_angle, np.ndarray):
            pass
        elif isinstance(mc_angle, basestring) and ';' in mc_angle:
            mc_angle = mc_angle.replace(';', ',')
            mc_angle = np.array(literal_eval('[{}]'.format(mc_angle)))
        elif mc_angle is None and self._l1type == 'AUXMETEO':
            if len(self.get_handler().dimensions['n_beam']) == 6:
                mc_angle = np.array([0, 2, 4, 6, 8, 10])
            else:
                raise Exception('Could not guess macrocycle angles.')
        else:
            raise Exception('Unexpected macrocycle_angle attribute.')
        self._mcycle['mcycle'] = {'angle': mc_angle}
        # Parse dimensions in order to finc n_beam* dims
        nat_dims = self.get_handler().dimensions
        beam_dnames = []
        for nat_dname in nat_dims.keys():
            if re.match(r'^n_beam.*$', nat_dname) is None:
                continue
            if nat_dname == 'n_beam':
                b_angle = mc_angle.copy()
            elif nat_dname == 'n_beamM1':
                b_angle = mc_angle[np.where(mc_angle > 0)]
            elif nat_dname == 'n_beam_l1b':
                b_angle = mc_angle[np.where(mc_angle >= 6)]
            else:
                txt = '{} is a new n_beam* dim ?'
                raise Exception(txt.format(nat_dname))
            beam_dnames.append(nat_dname)
            self._mcycle[nat_dname] = {'angle': b_angle}
        # Get index/count for mcycle and n_beam*
        for d in self._mcycle.keys():
            index = np.where(self._mcycle[d]['angle'] == self._incidence)[0].tolist()
            if len(index) != 0:
                if any([i + index[0] != index[i] for i in range(len(index))]):
                    txt = 'Not adjacent beams for incidence={} ?'
                    raise Exception(txt.format(self._incidence))
            self._mcycle[d]['index'] = index
            self._mcycle[d]['count'] = len(index)
        no_inc = (self._l1type == 'L1A' and \
                  self._mcycle['n_beam']['count'] == 0) or \
                  (self._l1type == 'L1B' and \
                   self._mcycle['n_beam_l1b']['count'] == 0)
        if no_inc:
            txt = 'incidence {} not found in {}'
            raise IncidenceError(txt.format(self._incidence, self.get_url()))
        # Parse variables
        nat_vars = self.get_handler().variables
        for nat_vname in nat_vars.keys():
            nat_dnames = list(nat_vars[nat_vname].dimensions)
            nat_dvals = [len(nat_dims[nat_dname]) for nat_dname in nat_dnames]
            is_numbered = re.match(r'^.*_[0-5]$', nat_vname) is not None
            has_nmcycles = 'n_mcycles' in nat_dnames
            beam_dname = [d for d in beam_dnames if d in nat_dnames]
            if len(beam_dname) > 1:
                txt = '{} with more than one n_beam* dim ?'
                raise Exception(txt.format(nat_vname))
            has_beam = len(beam_dname) == 1
            if is_numbered:
                if not has_nmcycles:
                    txt = '{} without n_mcycles dim ?'
                    raise Exception(txt.format(nat_vname))
                if has_beam:
                    txt = '{} with n_beam* dim ?'
                    raise Exception(txt.format(nat_vname))
                split_nat_vname = nat_vname.split('_')
                vname = '_'.join(split_nat_vname[0:-1])
                vbeam = split_nat_vname[-1]
                if int(vbeam) >= self._mcycle['mcycle']['index'][0] and \
                   int(vbeam) <= self._mcycle['mcycle']['index'][-1]:
                    dnames = []
                    dvals = []
                    for nat_dname, nat_dval in zip(nat_dnames, nat_dvals):
                        if nat_dname == 'n_mcycles':
                            dnames.append('time')
                        elif re.match(r'^.*_[0-5]$', nat_dname) is not None:
                            split_nat_dname = nat_dname.split('_')
                            dname = '_'.join(split_nat_dname[0:-1])
                            dbeam = split_nat_dname[-1]
                            if vbeam != dbeam:
                                txt = '{} with {} dim ?'
                                raise Exception(txt.format(nat_vname, nat_dname))
                            dnames.append(dname)
                        else:
                            dnames.append(nat_dname)
                        dvals.append(nat_dval)
                    if vname not in self._vars.keys():
                        self._vars[vname] = {'nat_vnames': [nat_vname],
                                             'nat_dnames': [nat_dnames],
                                             'nat_dvals': [nat_dvals],
                                             'dnames': dnames,
                                             'dvals': dvals,
                                             'action': 'merge_var'}
                    else:
                        chk_dnames = self._vars[vname]['dnames']
                        if len(dnames) != len(chk_dnames) or \
                           any([d1 != d2 for d1, d2 in zip(dnames, chk_dnames)]):
                            txt = '{} not consistent with {}_*'
                            raise Exception(txt.format(nat_vname, vname))
                        chk_dvals = self._vars[vname]['nat_dvals'][-1]
                        if len(dvals) != len(chk_dvals) or \
                           any([d1 != d2 for d1, d2 in zip(dvals, chk_dvals)]):
                            txt = '{} not consistent with {}_*'
                            raise Exception(txt.format(nat_vname, vname))
                        tmp_nat_vnames = list(self._vars[vname]['nat_vnames'])
                        tmp_nat_vnames.append(nat_vname)
                        tmp_nat_vnames.sort()
                        vindex = tmp_nat_vnames.index(nat_vname)
                        self._vars[vname]['nat_vnames'].insert(vindex, nat_vname)
                        self._vars[vname]['nat_dnames'].insert(vindex, nat_dnames)
                        self._vars[vname]['nat_dvals'].insert(vindex, nat_dvals)
                        time_pos = dnames.index('time')
                        self._vars[vname]['dvals'][time_pos] += dvals[time_pos]
            elif has_nmcycles:
                if not has_beam:
                    txt = '{} without n_beam* dim ?'
                    raise Exception(txt.format(nat_vname))
                beam_dname = beam_dname[0]
                if self._mcycle[beam_dname]['count'] == 0:
                    continue
                dnames = []
                dvals = []
                for nat_dname, nat_dval in zip(nat_dnames, nat_dvals):
                    if nat_dname == 'n_mcycles':
                        dnames.append('time')
                        dvals.append(nat_dval * self._mcycle[beam_dname]['count'])
                    elif nat_dname == beam_dname:
                        pass
                    else:
                        dnames.append(nat_dname)
                        dvals.append(nat_dval)
                self._vars[nat_vname] = {'nat_vnames': [nat_vname],
                                         'nat_dnames': [nat_dnames],
                                         'nat_dvals': [nat_dvals],
                                         'dnames': dnames,
                                         'dvals': dvals,
                                         'beam_dname': beam_dname,
                                         'action': 'merge_dim'}
            elif has_beam:
                if len(nat_dnames) != 1:
                    txt = '{} with n_beam* dim not alone ?'
                    raise Exception(txt.format(nat_vname))
                beam_dname = beam_dname[0]
                if self._mcycle[beam_dname]['count'] != 0:
                    self._vars2attr[nat_vname] = {'nat_vnames': [nat_vname],
                                                  'nat_dnames': [nat_dnames],
                                                  'nat_dvals': [nat_dvals],
                                                  'beam_dname': beam_dname}
            else:
                self._vars[nat_vname] = {'nat_vnames': [nat_vname],
                                         'nat_dnames': [nat_dnames],
                                         'nat_dvals': [nat_dvals],
                                         'dnames': nat_dnames,
                                         'dvals': nat_dvals,
                                         'action': 'none'}
        for vname in self._vars.keys():
            for dname, dval in zip(self._vars[vname]['dnames'], \
                                   self._vars[vname]['dvals']):
                if dname not in self._dims:
                    self._dims[dname] = dval
                elif self._dims[dname] != dval:
                    txt = 'Inconsistent dimension in {}'
                    raise Exception(txt.format(vname))
        return

    def get_geolocation_field(self, fieldname):
        """
        """
        if fieldname in ['time', 'lon', 'lat']:
            return fieldname
        else:
            return None

    def get_matching_dimname(self, geodimname):
        """
        """
        if geodimname in ['time']:
            return geodimname
        else:
            return None

    def get_standard_dimname(self, geodimname):
        """
        """
        return geodimname

    def open(self):
        """
        """
        return super(SWIMCWWICL1NCFile, self).open()

    def read_values(self, fieldname, slices=None):
        """
        """
        assert isinstance(slices, (type(None), list))
        if fieldname in self._vars.keys():
            _var = self._vars[fieldname]
            action = _var['action']
            if action == 'none':
                nat_vname = _var['nat_vnames'][0]
                values = super(SWIMCWWICL1NCFile, self).read_values(nat_vname, slices)
            elif action == 'merge_dim' or action == 'merge_var':
                # Make filled slices for fieldname
                dvals = _var['dvals']
                if slices is None:
                    tmp_slices = [slice(None)] * len(dvals)
                else:
                    tmp_slices = list(slices)
                filled_slices = super(SWIMCWWICL1NCFile, self)._fill_slices(tmp_slices, dvals)
                # Adapt slices and read native fieldname(s)
                # 'merge_dim' : incidence is split in a n_beam* dimension
                # 'merge_var' : incidence is split between *_{beam} variables
                # For both cases, the "beam" dimension is put after n_mcycles dimension
                # in order to facilitate the merge (simple reshape)
                nat_slices = list(filled_slices)
                nat_vnames = _var['nat_vnames']
                nat_dnames = _var['nat_dnames']
                if action == 'merge_dim':
                    # add a slice for n_beam* dimensions
                    beam_dname = _var['beam_dname']
                    beam_dpos = nat_dnames[0].index(beam_dname)
                    index = self._mcycle[beam_dname]['index']
                    nat_slices.insert(beam_dpos, slice(index[0], index[-1] + 1, 1))
                    # cancel n_mcycles slice (is input time slice)
                    nmcycle_dpos = nat_dnames[0].index('n_mcycles')
                    nat_slices[nmcycle_dpos] = slice(None)
                    # read variable and put n_beam* dim after n_mcycles dim
                    values = super(SWIMCWWICL1NCFile, self).read_values(nat_vnames[0], nat_slices)
                    values = np.rollaxis(values, beam_dpos, nmcycle_dpos + 2)
                elif action == 'merge_var':
                    # cancel n_mcycles slice (is input time slice)
                    nmcycle_dpos = nat_dnames[0].index('n_mcycles')
                    nat_slices[nmcycle_dpos] = slice(None)
                    # read variables, create a virtual "beam" dim and
                    # put it after n_mcycles dim
                    beam_values = []
                    for nat_vname in nat_vnames:
                        b_values = super(SWIMCWWICL1NCFile, self).read_values(nat_vname, nat_slices)
                        beam_values.append(b_values[np.newaxis, :])
                    values = np.ma.concatenate(beam_values, axis=0)
                    values = np.rollaxis(values, 0, nmcycle_dpos + 2)
                # Merge n_mcycles / "beam" dims into a time dim
                old_shape = values.shape
                new_shape = []
                for i in range(len(old_shape)):
                    if i == nmcycle_dpos:
                        new_shape.append(old_shape[i] * old_shape[i + 1])
                    elif i == nmcycle_dpos + 1:
                        pass
                    else:
                        new_shape.append(old_shape[i])
                values = values.reshape(new_shape, order='C')
                # Apply time slice
                if slices is not None:
                    time_slices = [slice(None)] * len(dvals)
                    dnames = _var['dnames']
                    time_dpos = dnames.index('time')
                    time_slices[time_dpos] = filled_slices[time_dpos]
                    values = values[time_slices]
                # Check all dimensions (to be removed in future)
                shape = values.shape
                if len(shape) != len(filled_slices):
                    txt = '{} with unexpected dims after merge.'
                    raise Exception(txt.format(fieldname))
                for shp, sli in zip(shape, filled_slices):
                    start, stop, step = sli.start, sli.stop, sli.step
                    if stop is None:
                        stop = -1
                    dimsize = 1 + (abs(stop - start) - 1) / abs(step)
                    if shp != dimsize:
                        txt = '{} with unexpected dims after merge.'
                        raise Exception(txt.format(fieldname))
            else:
                raise Exception('Unknown action for {}'.format(fieldname))
        elif fieldname in ['time', 'lon', 'lat'] and self._l1type == 'L1A':
            if slices is None:
                ntime = self.get_dimsize('time')
                slices = [slice(ntime)]
            nswath = self.get_dimsize('n_swath')
            slices.append(slice(nswath / 2, nswath / 2 + 1))
            if fieldname == 'time':
                time_field = self.read_field('time')
                # Get conversions factors
                # (ref time is the same -> no offsets)
                units = time_field.units
                nr_units = self.read_field('time_nr').units
                swath_units = self.read_field('time_swath').units
                if units.split(' ')[0] == 'milliseconds':
                    if nr_units.split(' ')[0] == 'seconds':
                        nr_fac = 1000.
                    elif nr_units.split(' ')[0] == 'microseconds':
                        nr_fac = 0.001
                    else:
                        raise Exception('time_nr units ?')
                    if swath_units.split(' ')[0] == 'microseconds':
                        swath_fac = 0.001
                    else:
                        raise Exception('time_swath units ?')
                else:
                    raise Exception('time units ?')
                # Read
                time_nr = self.read_values('time_nr', slices=[slices[0]])
                time_swath = self.read_values('time_swath', slices=slices)
                time_swath = time_swath.reshape(time_swath.shape[0])
                dtype = time_field.datatype
                values = time_nr.astype(dtype) * nr_fac + \
                         time_swath.astype(dtype) * swath_fac
            elif fieldname == 'lon':
                values = self.read_values('lon_l1a', slices=slices)
                values = values.reshape(values.shape[0])
            elif fieldname == 'lat':
                values = self.read_values('lat_l1a', slices=slices)
                values = values.reshape(values.shape[0])
        elif fieldname == 'time' and self._l1type == 'L1B':
            time_field = self.read_field('time')
            # Get conversions factors
            # (ref time is the same -> no offsets)
            units = time_field.units
            l1b_units = self.read_field('time_l1b').units
            if units.split(' ')[0] == 'milliseconds':
                if l1b_units.split(' ')[0] == 'seconds':
                    l1b_fac = 1000.
                else:
                    raise Exception('time_l1b units ?')
            else:
                raise Exception('time units ?')
            # Read
            time_l1b = self.read_values('time_l1b', slices=slices)
            dtype = time_field.datatype
            values = time_l1b.astype(dtype) * l1b_fac
        elif fieldname == 'lon' and self._l1type == 'L1B':
            values = self.read_values('lon_l1b', slices=slices)
        elif fieldname == 'lat' and self._l1type == 'L1B':
            values = self.read_values('lat_l1b', slices=slices)
        else:
            txt = 'Field {} not existing in mapper or NetCDF file.'
            raise Exception(txt.format(fieldname))
        return values

    def read_fillvalue(self, fieldname):
        """
        """
        # Is this function sometimes used ?
        attrs = self.read_field_attributes(fieldname)
        if '_FillValue' in attrs:
            fillvalue = attrs['_FillValue']
        elif 'fill_value' in attrs:
            fillvalue = attrs['fill_value']
        elif 'missing_value' in attrs:
            fillvalue = attrs['missing_value']
        else:
            fillvalue = None
        return fillvalue

    def read_field(self, fieldname):
        """
        """
        # Get NetCDF Variable
        naming_auth = self.get_naming_authority()
        attrs = self.read_field_attributes(fieldname)
        if 'long_name' in attrs:
            descr = attrs['long_name']
        else:
            descr = None
        if 'standard_name' in attrs:
            stdname = attrs['standard_name']
        else:
            stdname = None
        var = Variable(
            shortname=fieldname,
            description=descr,
            authority=naming_auth,
            standardname=stdname
            )
        dims = self.get_full_dimensions(fieldname)
        if 'scale_factor' in attrs:
            datatype = np.dtype(attrs['scale_factor'])
        elif 'add_offset' in attrs:
            datatype = np.dtype(attrs['add_offset'])
        else:
            if fieldname == 'time':
                datatype = np.dtype('float64')
            else:
                if fieldname in ['lon', 'lat']:
                    vname = fieldname + '_' + self._l1type.lower()
                else:
                    vname = fieldname
                _var = self._vars[vname]
                nat_vname = _var['nat_vnames'][0]
                datatype = self.get_handler().variables[nat_vname].dtype
        if '_FillValue' in attrs:
            fillvalue = attrs['_FillValue']
        elif 'fill_value' in attrs:
            fillvalue = attrs['fill_value']
        elif 'missing_value' in attrs:
            fillvalue = attrs['missing_value']
        else:
            fillvalue = None
        field = Field(var, dims, datatype=datatype, fillvalue=fillvalue)
        field.attach_storage(self.get_field_handler(fieldname))
        # MetaData
        field.attributes = {}
        field.units = None
        if 'units' in attrs:
            field.units = attrs['units']
        field.valid_min = None
        field.valid_max = None
        if 'valid_min' in attrs and 'valid_max' in attrs:
            # sometimes expressed as a string : need type cast
            try:
                field.valid_min = np.array(attrs['valid_min']).astype(
                    field.datatype)
                field.valid_max = np.array(attrs['valid_max']).astype(
                    field.datatype)
                if 'scale_factor' in attrs:
                    field.valid_min = field.valid_min * attrs['scale_factor']
                    field.valid_max = field.valid_max * attrs['scale_factor']
                if 'add_offset' in attrs:
                    field.valid_min = field.valid_min + attrs['add_offset']
                    field.valid_max = field.valid_max + attrs['add_offset']
            except:
                field.valid_min = attrs['valid_min']
                field.valid_max = attrs['valid_max']
                # logging.error("invalid valid_min or valid_max : %s %s",
                #               field.valid_min,
                #               field.valid_max)
        for att in attrs:
            if att not in ['units', 'scale_factor', 'add_offset',
                           '_FillValue', 'valid_min', 'valid_max']:
                field.attributes[att] = attrs[att]
        return field

    def get_fieldnames(self):
        """
        """
        return self._vars.keys()

    def get_dimensions(self, fieldname=None):
        """
        """
        if fieldname is None:
            dims = self._dims.keys()
        elif fieldname in ['time', 'lon', 'lat']:
            dims = ['time']
        elif fieldname in self.get_fieldnames():
            dims = self._vars[fieldname]['dnames']
        else:
            txt = 'Field {} not existing in NetCDF file or mapper.'
            raise Exception(txt.format(fieldname))
        return tuple(dims)

    # def get_full_dimensions(self, fieldname=None):
    #     Parent's method should work

    def get_dimsize(self, dimname):
        """
        """
        if dimname not in self._dims.keys():
            txt = 'Dimension {} not existing in NetCDF file or mapper.'
            raise Exception(txt.format(dimname))
        return self._dims[dimname]

    def read_field_attributes(self, fieldname):
        """
        """
        if fieldname in self._vars.keys():
            _var = self._vars[fieldname]
            nat_vname = _var['nat_vnames'][0]
            attrs = super(SWIMCWWICL1NCFile, self).read_field_attributes(nat_vname)
            if fieldname == 'time_nr':
                units = attrs['units']
                if 'since' not in units:
                    ref_time = self.read_global_attribute('reference_time')
                    units += ' since ' + ref_time
                split_units = units.split(' ')
                if split_units[0] == 's':
                    units = 'seconds ' + ' '.join(split_units[1:])
                elif split_units[0] == 'us':
                    units = 'microseconds ' + ' '.join(split_units[1:])
                else:
                    raise Exception('Which unit for {} ?'.format(fieldname))
                attrs['units'] = units
            elif fieldname == 'time_swath':
                units = attrs['units']
                if units == u'\xb5s':
                    units = 'microseconds'
                elif units == 'us':
                    units = 'microseconds'
                else:
                    raise Exception('Which unit for {} ?'.format(fieldname))
                attrs['units'] = units
            elif fieldname == 'time_l1b':
                units = attrs['units']
                if 'since' not in units:
                    if 'reference_time' in self.read_global_attributes():
                        ref_time = self.read_global_attribute('reference_time')
                    else: # assume same than L1A
                        ref_time = '2009-01-01T00:00:00'
                    units += ' since ' + ref_time
                split_units = units.split(' ')
                if split_units[0] == 's':
                    units = 'seconds ' + ' '.join(split_units[1:])
                else:
                    raise Exception('Which unit for {} ?'.format(fieldname))
                attrs['units'] = units
        elif fieldname == 'time':
            if self._l1type == 'L1A':
                temp = self.read_field_attributes('time_nr')
                long_name = 'Swath middle gate time'
            elif self._l1type == 'L1B':
                temp = self.read_field_attributes('time_l1b')
                long_name = 'Footprint center time'
            attrs = OrderedDict()
            attrs['long_name'] = long_name
            attrs['_FillValue'] = temp['_FillValue']
            split_units = temp['units'].split(' ')
            attrs['units'] = 'milliseconds ' + ' '.join(split_units[1:])
        elif fieldname == 'lon':
            if self._l1type == 'L1A':
                temp = self.read_field_attributes('lon_l1a')
                long_name = 'Swath middle gate longitude'
            elif self._l1type == 'L1B':
                temp = self.read_field_attributes('lon_l1b')
                long_name = 'Footprint center longitude'
            attrs = OrderedDict()
            attrs['long_name'] = long_name
            attrs['_FillValue'] = temp['_FillValue']
            attrs['units'] = temp['units']
        elif fieldname == 'lat':
            if self._l1type == 'L1A':
                temp = self.read_field_attributes('lat_l1a')
                long_name = 'Swath middle gate latitude'
            elif self._l1type == 'L1B':
                temp = self.read_field_attributes('lat_l1b')
                long_name = 'Footprint center latitude'
            attrs = OrderedDict()
            attrs['long_name'] = long_name
            attrs['_FillValue'] = temp['_FillValue']
            attrs['units'] = temp['units']
        else:
            txt = 'Field {} not existing in mapper or NetCDF file.'
            raise Exception(txt.format(fieldname))
        return attrs

    def read_global_attributes(self):
        """
        """
        attrs = self._vars2attr.keys() + \
                super(SWIMCWWICL1NCFile, self).read_global_attributes()
        return attrs

    def read_global_attribute(self, name):
        """
        """
        if name in self._vars2attr.keys():
            _var = self._vars2attr[name]
            nat_vname = _var['nat_vnames'][0]
            atts = super(SWIMCWWICL1NCFile, self).read_values(nat_vname)
            atts = atts[self._mcycle[_var['beam_dname']]['index']]
            if atts.min() != atts.max():
                txt = '{} with different values for one incidence ?'
                raise Exception(txt.format(name))
            att = atts[0]
        elif name in super(SWIMCWWICL1NCFile, self).read_global_attributes():
            att = super(SWIMCWWICL1NCFile, self).read_global_attribute(name)
            if name in ['reference_time', 'creation_date',
                        'start_time', 'stop_time',
                        'start_date', 'stop_date']:
                att = att.replace(';', '')
                att = att.replace('UTC=', '')
            elif isinstance(att, basestring) and ';' in att:
                try:
                    att = att.replace(';', ',')
                    att = np.array(literal_eval('[{}]'.format(att)))
                except:
                    pass
            if isinstance(att, np.ndarray):
                if name in ['macrocycle_angle', 'macrocycle_beam']:
                    pass
                else:
                    nvalues = att.size
                    nbeams = [v['angle'].size for v in self._mcycle.values()]
                    if nvalues in nbeams:
                        dim = self._mcycle.keys()[nbeams.index(nvalues)]
                        if self._mcycle[dim]['count'] == 0:
                            att = None
                        else:
                            att = att[self._mcycle[dim]['index']]
                            if att.min() != att.max():
                                txt = '{} with different values for one incidence ?'
                                raise Exception(txt.format(name))
                            att = att[0]
        else:
            txt = 'Attribute {} not existing in NetCDF file.'
            raise Exception(txt.format(name))
        return att

    def get_start_time(self):
        """
        """
        if 'time_coverage_start' in self.read_global_attributes():
            attrdate = self.read_global_attribute('time_coverage_start')
        else:
            start_date = self.read_global_attribute('start_date')
            start_time = self.read_global_attribute('start_time')
            attrdate = start_date + 'T' + start_time
        return parser.parse(attrdate)

    def get_end_time(self):
        """
        """
        if 'time_coverage_end' in self.read_global_attributes():
            attrdate = self.read_global_attribute('time_coverage_end')
        else:
            stop_date = self.read_global_attribute('stop_date')
            stop_time = self.read_global_attribute('stop_time')
            attrdate = stop_date + 'T' + stop_time
        return parser.parse(attrdate)

    # def get_spatial_resolution_in_deg(self):
    #     """
    #     """
    #     raise NotImplementedError

    # def get_product_version(self):
    #     """
    #     """
    #     raise NotImplementedError

    # def get_bbox(self):
    #     """
    #     """
    #     raise NotImplementedError

    # def get_cycle_number(self):
    #     """
    #     """
    #     raise NotImplementedError

    # def get_orbit_number(self):
    #     """
    #     """
    #     raise NotImplementedError
