==============================
Working with ODYSSEA data
==============================

**ODYSSEA** SST data can be considered from two feature perspectives :
  * GridTimeSeries
  * Grid

Using ODYSSEA data as GridTimeSeries
====================================

We can ignore the fact that each grid at each time step is stored in an individual file, by using a storage class aggregating all files altogether. Only the file pattern is to be provided to this class.

.. doctest::

 >>> import data.gridtimeteries
 >>> import mapper.urlseries
 >>> 
 >>> l4_pattern = '/home2/taveeg/data/operation/project/myocean/sst-tac/odyssea/v2/med/analysed_sst_002/%Y/%j/%Y%2m%2d-IFR-L4_GHRSST-SSTfnd-ODYSSEA-MED_002-v2.0-fv1.0.nc'
 >>> odysseaFiles = mapper.urlseries.URLSeries( urlpattern=l4_pattern )
 >>> odysseaTimeSeries = data.gridtimeseries.GridTimeSeries()
 >>> odysseaTimeSeries.load( odysseaFiles )



Using ODYSSEA data as grid
==========================

If no temporal aspect is considered, and one wants to work on one or a few independant time steps, ODYSSEA data grids can be accessed directly.

.. doctest::

 >>> import datamodel.grid
 >>> import mapper.ncfile
 >>> ncf = mapper.ncfile.NCFile( URL = "/home2/taveeg/data/operation/project/myocean/sst-tac/odyssea/v2/med/analysed_sst_002/2011/260/20110917-IFR-L4_GHRSST-SSTfnd-ODYSSEA-MED_002-v2.0-fv1.0.nc" )
 >>> g = datamodel.grid.Grid()
 >>> g.load( ncf )

Displaying data
---------------

.. doctest::

 >>> g.display_map( 'analysed_sst', output='toto.png' )

