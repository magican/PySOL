#	this example program calls the python code for reading SST 
#   files from RSS, Version-4 files.  It is currently set to write out 
#   the data in the verification file,  verify_sst_v04.txt.
#	To test, place the data files in the same directory as 
#   this code along with the verify_sst_v04.txt file
#   if problems with this code, contact support@remss.com

from mw_sst_v4 import OISST

def read_data(filename='mw.fusion.2004.140.v04.0.gz'):
    dataset = OISST(filename, missing=missing)
    if not dataset.variables: sys.exit('problem reading file')
    return dataset

ilon = (769,774)
ilat = (473,477)
avar = 'sst'
missing = -999.

#----------------------------------------------------------------------------

def show_dimensions(ds):
    print
    print 'Dimensions'
    for dim in ds.dimensions:
        print ' '*4, dim, ':', ds.dimensions[dim]

def show_variables(ds):
    print
    print 'Variables:'
    for var in ds.variables:
        print ' '*4, var, ':', ds.variables[var].long_name

def show_validrange(ds):
    print
    print 'Valid min and max and units:'
    for var in ds.variables:
        print ' '*4, var, ':', \
              ds.variables[var].valid_min, 'to', \
              ds.variables[var].valid_max,\
              '(',ds.variables[var].units,')'

def show_somedata(ds):
    print
    print 'Show some data for:',avar
    print 'index range: (' + \
          str(ilat[0]) + ':' + str(ilat[1]) + ' ,' + \
          str(ilon[0]) + ':' + str(ilon[1]) + ')'
    print ds.variables[avar][ilat[0]:ilat[1]+1, ilon[0]:ilon[1]+1]

#----------------------------------------------------------------------------

if __name__ == '__main__':    
    import sys    
    dataset = read_data()
    show_dimensions(dataset)
    show_variables(dataset)
    show_validrange(dataset)
    show_somedata(dataset)
    print
    print 'done'
