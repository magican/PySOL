=====================================
Examples for some CERSAT gridded data
=====================================

CO2 Exchange coefficient
------------------------

.. doctest::

 >>> import datamodel.grid
 >>> import mapper.ncfile
 >>> ncf = mapper.ncfile.NCFile( URL = "/home3/datacer/htdocs/gridded/cersat/k_co2/quikscat/weekly/2009/200910050000-200910120000.nc" )
 >>> g = datamodel.grid.Grid()
 >>> g.load( ncf )
 >>> val = g.get_values('k_liss_merlivat')
 >>> g.display_map(data=val, palette=plt.cm.rainbow, pretty=True, range=[0,0.15,0.005], output='toto.png')


Mean wind fields
----------------

.. doctest::
 
 >>> import datamodel.grid
 >>> import mapper.ncfile
 >>> ncf = mapper.ncfile.NCFile( URL = "/home3/datacer/htdocs/gridded/cersat/wind/quikscat/daily/2009/200908040000-200908050000.nc" )
 >>> g = datamodel.grid.Grid()
 >>> g.load( ncf )
 >>> val = g.get_values('wind_speed')
 >>> g.display_map(data=val, palette='medspiration', pretty=True, range=[0,15,0.5], output='toto.png')


