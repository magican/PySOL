Tested datasets
===============

**cerbere** has been tested with the following list of products:

NCFile
------

+------------------------------------------------------+---------------+
| Dataset                                              + Datamodel     |
+======================================================+===============+
| SMOS SSS by CEC-OS (N. Reul)                         | Grid          |
+------------------------------------------------------+---------------+
| MLD climatology by Ifremer (Boyer)                   | Grid          |
+------------------------------------------------------+---------------+
| WaveWatch3/Hindcast by Ifremer                       | Grid          |
+------------------------------------------------------+---------------+

Example::

    >>> from cerbere.datamodel.grid import Grid
    >>> from cerbere.mapper.ncfile import NCFile
    >>> 
    >>> ncf = NCFile(url='/home/salinity1/public/data/cecos/sss_smos_l3_V02/2012/Daily_composite/Half_degree/SSS_SMOS_L3_Daily_0.5deg_CATDS_CECOS_2012.12.30_V02.nc')
    >>> g = Grid()
    >>> g.load(ncf)

GHRSSTNCFile
------------

+------------------------------------------------------+---------------+
| Dataset                                              + Datamodel     |
+======================================================+===============+
| VIIRS L2P by NAVO                                    | Swath         |
+------------------------------------------------------+---------------+
| AMSR2 L2P by JAXA                                    | Swath         |
+------------------------------------------------------+---------------+
| (A)ATSR L2P by ESA SST CCI/ARC                       | Swath         |
+------------------------------------------------------+---------------+

SAFESLFile
------------

Mapper for Sentinel-3 SLSTR files in SAFE format.

This case is quite complex because of the different sub-products
(with different grids and views) existing within a single SLSTR SAFE container. Because these
sub-products don't all overlap (e.g. they don't have the same array dimensions nor same
pixel locations), we can not extract a single miniprod from a SAFE container.

The data in felyx are read with `**cerbere**<http://cerbere.readthedocs.org/>`_ python package.
The `**cerbere** strategy for SLSTR products is to provide several mappers adapted to
each product/subproduct combination. A mapper will merge the nadir and oblique view fields as
if they were in the same product file. 

There are 4 different possible mappers for a S3A_SL_1_RBT product which means you
will have to configure 4 datasets in felyx if you want to extract miniprod for all possible
subproducts. This is the more tricky case. Fortunately there is only one possible mapper
for any other S3 SLSTR product!

+------------------------------------------------+------------------------+---------------+
| Dataset                                        | Mapper class           | Datamodel     |
+================================================+========================+===============+
| S3A_SL_2_WCT Nadir (in)                        | SAFESLIRFile           | Swath         |
| S3A_SL_1_RBT 1km Nadir (in)                    |                        |               |
| S3A_SL_2_WCT Oblique (io)                      |                        |               |
| S3A_SL_1_RBT 1km Oblique (io)                  |                        |               |
+------------------------------------------------+------------------------+---------------+
| S3A_SL_1_RBT 500m & SWIR A Stripe Nadir (an)   | SAFESL500AFile         | Swath         |
| S3A_SL_1_RBT 500m & SWIR A Stripe Oblique (ao) |                        |               |
+------------------------------------------------+------------------------+---------------+
| S3A_SL_1_RBT 500m & SWIR B Stripe Nadir (bn)   | SAFESL500BFile         | Swath         |
| S3A_SL_1_RBT 500m & SWIR B Stripe Oblique (bo) |                        |               |
+------------------------------------------------+------------------------+---------------+
| S3A_SL_1_RBT TDI Nadir (cn)                    | SAFESL500TDIFile       | Swath         |
| S3A_SL_1_RBT TDI Oblique (co)                  |                        |               |
+------------------------------------------------+------------------------+---------------+


