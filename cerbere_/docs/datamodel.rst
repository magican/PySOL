==========
Data Model
==========

The data model correspond to the main sampling patterns usually available for Earth Observation data. For each typical sampling pattern, cerbere provides a unique representation for them, so that common data handling or display operations are implemented uniquely and through a common interface.

Structure
=========

The following table describes the dimensions and geolocation variables associated with each data model :

+-------------------------+--------------------+---------------------------------+--------------------------------+
| Pattern                 | Dimensions (size)  | Geolocation fields (dimensions) | Geophysical field (dimensions) |
+=========================+====================+=================================+================================+
| Swath                   |                    | - time (row,cell)               | V (row, cell)                  |
|                         | - row (y)          | - lat (row,cell)                |                                |
|                         | - cell (x)         | - lon (row,cell)                |                                |
+-------------------------+--------------------+---------------------------------+--------------------------------+
| Image                   | - time (1)         | - time (time)                   | V (row, cell)                  |
|                         | - row (y)          | - lat (row,cell)                |                                |
|                         | - cell (x)         | - lon (row,cell)                |                                |
+-------------------------+--------------------+---------------------------------+--------------------------------+
| grid                    | - time (1)         | - time (time)                   | V (y, x)                       |
|                         | - y (y)            | - lat (y,x)                     |                                |
|                         | - x (x)            | - lon (y,x)                     |                                |
+-------------------------+--------------------+---------------------------------+--------------------------------+
| grid time series        | - time (t)         | - time (time)                   | V (time, y, x)                 |
|                         | - y (y)            | - lat (y,x)                     |                                |
|                         | - x (x)            | - lon (y,x)                     |                                |
+-------------------------+--------------------+---------------------------------+--------------------------------+
| point collection        | - station (x)      | - time (station)                |                                |
|                         |                    | - lat (station)                 |                                |
|                         |                    | - lon (station)                 |                                |
+-------------------------+--------------------+---------------------------------+--------------------------------+
| time series at a point  | - station (1)      | - time (t)                      |                                |
|                         | - time (t)	       | - lat (station)                 |                                |
|                         |                    | - lon (station)                 |                                |
+-------------------------+--------------------+---------------------------------+--------------------------------+
| trajectory              | - time (t)         | - time (t)                      |                                |
|                         |                    | - lat (t)                       |                                |
|                         |                    | - lon (t)                       |                                |
+-------------------------+--------------------+---------------------------------+--------------------------------+
| spatial section         |                    |                                 |                                |
+-------------------------+--------------------+---------------------------------+--------------------------------+
| section time series     |                    |                                 |                                |
+-------------------------+--------------------+---------------------------------+--------------------------------+


