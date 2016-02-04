#This routine reads version-4 RSS MW-only, Optimal Interpolated SST fusion daily files 
#	You must UNZIP FILES before reading them
#
#	INPUT
# 	file_name  with path in form mw.yyyy.doy.v04.0
#	where yyyy = year
#		  doy  = day of year
#
#	OUTPUT  
#	sst_data 	 (a 1440 x 720 real*4 array of Sea Surface Temperature)
#	error_data	 (a 1440 x 720 real*4 array of the interpolation error estimate)
#	mask_data  	 (a 1440 x 720 integer*1 array of data masking information)
#			mask_data 
#				bit 0 = 1 for land			
#				bit 1 = 1 for ice
#				bit 2 = 1 for IR data used
#				bit 3 = 1 for MW data used	
#				bit 4 = 1 for bad data	
#				bits 5 - 7 not used

#        'land' : 'Is Land'
#         'ice' : 'Is Ice'
#     'uses_ir' : 'IR Data Used'
#     'uses_mw' : 'MW Data Used'
#         'bad' : 'Bad Data'
#   'longitude' : 'Grid Cell Center Longitude'
#    'latitude' : 'Grid Cell Center Latitude'
#      'nodata' : 'Is there no data?'
#
# 	  xcell=grid cell values between 1 and 1440
#     ycell=grid cell values between 1 and  720
#     dx=360./1440.
#     dy=180./720.
#	  Center of grid cell Longitude  = dx*xcell-dx/2.    (degrees east)
#	  Center of grid cell Latitude   = dy*ycell-(90+dy/2.)  (-90 to 90)
#
#	Please read the data description on www.remss.com/measurements/sea-surface-temperature
#	To contact RSS support: http://www.remss.com/support

from bytemaps import sys
from bytemaps import Dataset
from bytemaps import Verify

from bytemaps import btest
from bytemaps import get_data


class OISST(Dataset):
    """ Read daily MW -only OI SST bytemaps. """
    """
    Public data:
        filename = name of data file
        missing  = fill value used for missing data;
                  if None, then fill with byte codes (251-255)
        dimensions = dictionary of dimensions for each coordinate
        variables = dictionary of data for each variable
    """

    def __init__(self, filename, missing=-999.):
        """
        Required arguments:
            filename = name of data file to be read (string)
                
        Optional arguments:
            missing = fill value for missing data,
                      default is the value used in verify file
        """
        self.filename = filename
        self.missing = missing
        Dataset.__init__(self)

    # Dataset:

    def _attributes(self):
        return ['coordinates','long_name','units','valid_min','valid_max']

    def _coordinates(self):
        return ('variable','latitude','longitude')

    def _shape(self):
        return (3,720,1440)

    def _variables(self):
        return ['sst','error','mask','land','ice','uses_ir','uses_mw','bad',
                'longitude','latitude','nodata']                

    # _default_get():
    
    def _get_index(self,var):
        return {'sst' : 0,
                'error' : 1,
                'mask' : 2,
                }[var]

    def _get_scale(self,var):
        return {'sst' : 0.15,
                'error' : 0.005,
                }[var]

    def _get_offset(self,var):
        return {'sst':-3.0,
                }[var]

    # _get_ attributes:
    
    def _get_long_name(self,var):
        return {'sst' : 'Sea Surface Temperature',
                'error' : 'Interpolation Error Estimate',
                'mask' : 'Data Mask',
                'land' : 'Is Land',
                'ice' : 'Is Ice',
                'uses_ir' : 'IR Data Used',
                'uses_mw' : 'MW Data Used',
                'bad' : 'Bad Data',
                'longitude' : 'Grid Cell Center Longitude',
                'latitude' : 'Grid Cell Center Latitude',
                'nodata' : 'Is there no data?',
                }[var]

    def _get_units(self,var):
        return {'sst': 'deg C',
                'error' : 'deg C',
                'mask' : 'no dimensions',
                'land' : 'True or False',
                'ice' : 'True or False',
                'uses_ir' : 'True or False',
                'uses_mw' : 'True or False',
                'bad' : 'True or False',
                'longitude' : 'degrees east',
                'latitude' : 'degrees north',
                'nodata' : 'True or False',
                }[var]

    def _get_valid_min(self,var):
        return {'sst' : -3.0,
                'error' : 0.0,
                'longitude' : 0.0,
                'latitude' : -90.0,
                'land' : False,
                'ice' : False,
                'nodata' : False,
                'uses_ir' : False,
                'uses_mw' : False,
                'bad' : False,
                'mask' : None,
                }[var]

    def _get_valid_max(self,var):
        return {'sst' : 34.5,
                'error' : 1.25,
                'longitude' : 360.0,
                'latitude' : 90.0,
                'land' : True,
                'ice' : True,
                'nodata' : True,
                'uses_ir' : True,
                'uses_mw' : True,
                'bad' : True,
                'mask' : None,
                }[var]

    # _get_ variables:

    def _get_land(self,var,bmap):
        indx = self._get_index('mask')
        return get_data(btest(bmap,ipos=0),indx=indx)
        
    def _get_ice(self,var,bmap):
        indx = self._get_index('mask')
        return get_data(btest(bmap,ipos=1),indx=indx)

    def _get_uses_ir(self,var,bmap):
        indx = self._get_index('mask')
        return get_data(btest(bmap,ipos=2),indx=indx)
    
    def _get_uses_mw(self,var,bmap):
        indx = self._get_index('mask')
        return get_data(btest(bmap,ipos=3),indx=indx)

    def _get_bad(self,var,bmap):
        indx = self._get_index('mask')
        return get_data(btest(bmap,ipos=4),indx=indx)
    

class VerifyDataset(Verify):
    """ Contains info for verification. """
    
    def __init__(self,dataset):
        self.filename = 'verify_sst_v04.txt'
        self.ilon1 = 770
        self.ilon2 = 775
        self.ilat1 = 474
        self.ilat2 = 478
        self.variables = ['sst']
        self.startline = { 'sst' : 9 }
        Verify.__init__(self,dataset)
        

if __name__ == '__main__':
    """ Automated testing. """

    # read:
    sst = OISST('mw.fusion.2004.140.v04.0.gz')
    if not sst.variables: sys.exit('problem reading file')
    
    # verify:
    verify = VerifyDataset(sst)
    if verify.success: print 'successful verification for daily'
    else: sys.exit('verification failed for daily')
    print

    print 'all tests completed successfully'
    
