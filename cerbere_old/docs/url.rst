Using URL and URL series
========================

URL series
----------
URL series are used for featured spanning over multiple files : this can generally be the case for instance of time series (GridTimeSeries and PointTimeSeries)
 or trajectories. In such case information must be provided to the mapper class in order to select the proper file(s) when reading data.

* simple url pattern

a simple url pattern where time is expressed in a pythonic way (using python time directive as defined by datetime package). 
By convention, the specified time expresses the reference time of the first feature record in the file, which will be :
  * the time of the grid for a file containing a single grid
  * the time of the first grid for a file containing multiple grids
  * the time of the first record for a point time series
 

* url pattern with start and end time
 
This can be used for files containing multiple records spanning over a period of time given in the file name. The time directives have 
to be appended with the #S (start) and #E (end) that they represent.

 
