# encoding: utf-8
"""
cerbere.mapper.hy2ansoashdffile
============================

Mapper class for the HY2-A scatterometer HDF files from NSOAS

:copyright: Copyright 2014 Ifremer / Cersat.
:license: Released under GPL v3 license, see :ref:`license`.

.. sectionauthor:: agrouaze <antoine.grouazel@ifremer.fr>
.. codeauthor:: agrouaze
"""
from datetime import datetime
from collections import OrderedDict

#from numpy import dtype, float32, resize
import numpy
from netCDF4 import num2date,date2num
import pdb
from .. import READ_ONLY, DEFAULT_TIME_UNITS
from .hdf5file import HDF5File
from .abstractmapper import AbstractMapper
from ..datamodel.field import Field
from ..datamodel.variable import Variable
import netCDF4
# list_variable_artifical_name = {
#                                 '/wind_speed_selection':'wind_speed_selection',
#                                 '/wvc_lat':'wvc_lat',
#                                 '/wvc_lon':'wvc_lon',
#                                 /row_time
#                                 /wind_dir_selection
#                                 /max_likelihood_est
#                                 }
time_units = "days since 1950-01-01T00:00:00Z"
class HY2ANSOASHDFFile(HDF5File):
    def __init__(self, url=None, mode=READ_ONLY, **kwargs):
        """
        """
        HDF5File.__init__(self, url=url, mode=mode, **kwargs)
        return
        
    def get_fieldnames(self):
        '''
        Returns the list of geophysical fields stored for the feature
        note: in the hdf files there is only wind_dir(y,x,4), wind_speed(y,x,4) and wind_*_selection(y,x) 
        '''
        fields = super(HY2ANSOASHDFFile,self).get_fieldnames()
        for xx,fifi in enumerate(fields):
            if fifi[0]=='/':
                fields[xx] = fifi[1:]
        #fields = self.get_handler().datasets().keys()
        # remove here time/space information to keep only geophysical fields
        if 'row_time' in fields:
            fields.remove('row_time')
        if 'wvc_lat' in fields:
            fields.remove('wvc_lat')
        if 'wvc_lon' in fields:
            fields.remove('wvc_lon')
        if 'wind_dir' in fields:
            fields.remove('wind_dir')
            fields.remove('wind_speed') #fields splitted into solution1 2 3 4
        #fields.append['wind_speed_solution_1']
        #fields.append['wind_speed_solution_2']
        #fields.append['wind_speed_solution_3']
        #fields.append['wind_speed_solution_4']
        if 'wind_speed_solution_1' not in fields:
            fields.extend(['wind_speed_solution_1', 'wind_speed_solution_2','wind_speed_solution_3','wind_speed_solution_4'])
            fields.extend(['wind_dir_solution_1', 'wind_dir_solution_2','wind_dir_solution_3','wind_dir_solution_4'])
        return fields


    def get_geolocation_field(self, fieldname):
        matching = {'lat': 'wvc_lat',
                    'lon': 'wvc_lon',
                    'time': 'row_time'
                    }
        if fieldname in matching:
            return matching[fieldname]
        else:
            return fieldname
            

    def get_matching_dimname(self, geodimname):
        """
        Return the equivalent name in the current data mapper for a standard
        feature dimension, or None if not found.
        """
        matching = {
                'row': 'l2b_expected_wvc_rows',
                'cell': 'Number of Pixel Control Points',
                'time': 'row_time',
                'vec_1': 'vec',
                'vec_2': 'vec',
                }
        if geodimname in matching:
            return matching[geodimname]
        else:
            return None

    def get_standard_dimname(self, geodimname):
        """Returns the equivalent standard dimension for a geolocation
        dimension in the current data mapper.
        """
        matching = {            
                'l2b_expected_wvc_rows': 'row',
                'Pixels per Scan Line': 'cell',
                'Number of Pixel Control Points': 'cell',
                'row_time': 'time'
                }
        if geodimname in matching:
            return matching[geodimname]
        else:
            return None
            
    def get_dimsize(self,dimname):
        #special case for hy2ansoashdffile only row and cell dimensions
        if dimname == 'time':
            return self.get_dimsize('row')
        elif dimname == 'cell':
            return self.get_handler().get('/wind_speed').shape[1]
        elif dimname == 'row':
            return self.get_handler().get('/wind_speed').shape[0]


    def get_start_time(self):
        date = self.read_global_attribute('rangeBeginningTime')
        date = date[0]
        try:
            res = datetime.strptime(date,
                    '%Y%m%dT%H:%M:%S.%f')
        except:
            subset = self.read_values('time')>0
            tmp = self.read_values('time')[subset]
            min = numpy.amin(tmp )
#             tmp = self.read_values('time')[0,0]
            res = netCDF4.num2date(min,time_units)
        return res

    def get_end_time(self):
        date = self.read_global_attribute('rangeEndingTime')
        date = date[0]
        try:
            res = datetime.strptime(date,'%Y%m%dT%H:%M:%S.%f')
        except:
            subset = self.read_values('time')>0
            max = numpy.amax(self.read_values('time')[subset])
            res = netCDF4.num2date(max,time_units)
        return res
        
    def get_full_dimensions(self,fieldname):
        """
        return the dimension names and sizes of a field

        :param fieldname: name of the field
        :type fieldname: str

        :rtype: OrderedDict
        :return: Ordered dictionary where keys are the dimension names
            and values are their respective size
        """
        res = self.get_dimensions(fieldname=fieldname,full=True)
        return res
    
    def read_field(self, fieldname): #add this to filled inexistant attr for rowtime
        native = self.get_geolocation_field(fieldname)
        if native==None:
            native = fieldname
        if native=='row_time':
            namingauth = None
            descr='date of the measure'
            variable= Variable(
                        shortname=fieldname,
                        description=descr,
                        authority=namingauth,
                        standardname=None
                        )
            dims = self.get_dimensions(fieldname)
            #gotchas trick to avoid inheritance from last call to Field
            toto = Field(variable=None,datatype=None,dimensions=None,attributes=None, units=None, valid_min=None, valid_max=None)
            toto.units='rara'
            rec = Field(
                variable,
                dims,
                datatype=numpy.dtype('float64')#float#'float64'
                )
            rec.attach_storage(self.get_field_handler(fieldname))
            # MetaData
            rec.units = time_units
            rec.valid_min = datetime(2010,1,1).toordinal()
            rec.valid_max = datetime(2040,1,1).toordinal()
            return rec
        elif 'wind_speed_solution_' in native:
            namingauth = None
            descr = 'wind speed solution '+fieldname[-1]
            variable = Variable(
                        shortname=fieldname,
                        description=descr,
                        authority=namingauth,
                        standardname=None
                        )
            dims = self.get_dimensions('/wind_speed_selection')
            #gotchas trick to avoid inheritance from last call to Field
            toto = Field(variable=None,datatype=None,dimensions=None,attributes=None, units=None, valid_min=None, valid_max=None)
            roc = Field(
                variable,
                dims,
                datatype=numpy.dtype('float64')#float#'float64'
                )

            aatrs = self.read_field_attributes('/wind_speed_selection')
            #roc.attach_storage(self.get_field_handler(native))
            for tt in aatrs:
                roc.attributes[tt]=aatrs[tt]
            return roc
        elif 'wind_dir_solution_' in native:
            namingauth = None
            descr='wind direction solution '+fieldname[-1]
            variable= Variable(
                        shortname=fieldname,
                        description=descr,
                        authority=namingauth,
                        standardname=None
                        )
            dims=self.get_dimensions('/wind_dir_selection')
            #gotchas trick to avoid inheritance from last call to Field
            toto=Field(variable=None,datatype=None,dimensions=None,attributes=None, units=None, valid_min=None, valid_max=None)
            roc = Field(
                variable,
                dims,
                datatype=numpy.dtype('float64')#float#'float64'
                )

            aatrs = self.read_field_attributes('/wind_dir_selection')
            roc.attach_storage(self.get_field_handler(native))
            for tt in aatrs:
                roc.attributes[tt]=aatrs[tt]
            return roc
        else:
            res=super(HY2ANSOASHDFFile,self).read_field(fieldname=fieldname)
            return res
            
    def get_dimensions(self,fieldname,full=False):
        res = OrderedDict([
            ('row',self.get_dimsize('row')),
            ('cell',self.get_dimsize('cell')),
            ])
        return res
        
    def read_field_attributes(self, fieldname): 
        if fieldname in ['cell','row']:
            return None
        elif fieldname in ['row_time','time']:
            res=self.get_handler().get('row_time').attributes
            return res
        elif '_solution' in fieldname:
            fifi=self.read_field(fieldname)
            res=fifi.attributes
            return res
        else:
            res=super(HY2ANSOASHDFFile,self).read_field_attributes(fieldname=fieldname)
            return res
            
    def read_values(self, fieldname, slices=None, indices=None):
        native = self.get_geolocation_field(fieldname)
        if native==None:
            native=fieldname
        if native == 'row_time':
            rowtime = self.get_handler().get('row_time')
            rowtime = rowtime[:]
            res = numpy.ma.array([])
            for i,gg in enumerate(rowtime):
                if len(gg)!=21: #case where the time is not defined blank line
                    totp = datetime(1950,1,1,1)#.toordinal()
                    res = numpy.append(res,totp)
                else:
                    tn = datetime.strptime(gg,'%Y%m%dT%H:%M:%S.%f')
                    res = numpy.append(res,tn)
            res = date2num(res,'days since 1950-01-01T00:00:00Z','proleptic_gregorian')
            res = numpy.ma.array(res,dtype=float)
            azimuth_dim,range_dim = self.read_values('/wvc_lat').shape
            res = numpy.tile(res,(range_dim,1)).transpose()
            return res
        elif native in  ['wvc_lat','wvc_lon']:
            res = self.get_handler().get(native)
            res = numpy.ma.masked_array(res)
            res[res==0] = numpy.ma.masked
            if native=='wvc_lon':
                res[res>180] = res[res>180]-360
            return res
        elif 'wind_speed_solution' in native:
            solution_number = int(native[-1])
            
            tmp = self.get_handler().get('/wind_speed')
            tmp = numpy.ma.masked_array(tmp)
            res = tmp[:,:,solution_number-1]
            return res
        elif 'wind_dir_solution' in native:
            solution_number = int(native[-1])
            
            tmp = self.get_handler().get('/wind_dir')
            tmp = numpy.ma.masked_array(tmp)
            res = tmp[:,:,solution_number-1]
            return res
        else:
            if native[0]!='/':
                native = '/'+native
            print 'native trasnformed',native
            res = super(HY2ANSOASHDFFile, self).read_values(native)
            res = numpy.ma.masked_array(res)
            if '_selection' in native:
                res = res/100.0
            return res
            
    def get_product_version(self):
        """return the product version"""
        attrs = self.read_global_attributes()
        return attrs['longName']
        
    def get_orbit_number(self):
        attrs = self.read_global_attributes()
        tmp = attrs['inputPointer'][0]
        res = tmp.split('_')[2]
        return res[0:5]
        
    def get_cycle_number(self):
        return None
    
