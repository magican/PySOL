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

+------------------------------------------------+------------------------+---------------+
| Dataset                                        | Mapper class           | Datamodel     |
+================================================+========================+===============+
| S3A_SL_2_WCT Nadir (in)                        | SAFESLIRNadirFile      | Swath         |
| S3A_SL_1_RBT 1km Nadir (in)                    |                        |               |
+------------------------------------------------+------------------------+---------------+
| S3A_SL_2_WCT Oblique (io)                      | SAFESLIRObliqueFile    | Swath         |
| S3A_SL_1_RBT 1km Oblique (io)                  |                        |               |
+------------------------------------------------+------------------------+---------------+
| S3A_SL_1_RBT 500m & SWIR A Stripe Nadir (an)   | SAFESL500ANadirFile    | Swath         |
+------------------------------------------------+------------------------+---------------+
| S3A_SL_1_RBT 500m & SWIR A Stripe Oblique (ao) | SAFESL500AObliqueFile  | Swath         |
+------------------------------------------------+------------------------+---------------+
| S3A_SL_1_RBT 500m & SWIR B Stripe Nadir (bn)   | SAFESL500BNadirFile    | Swath         |
+------------------------------------------------+------------------------+---------------+
| S3A_SL_1_RBT 500m & SWIR B Stripe Oblique (bo) | SAFESL500BObliqueFile  | Swath         |
+------------------------------------------------+------------------------+---------------+
| S3A_SL_1_RBT TDI Nadir (cn)                    | SAFESL500TDINadirFile  | Swath         |
+------------------------------------------------+------------------------+---------------+
| S3A_SL_1_RBT TDI Oblique (co)                  | SAFESL500TDIObliqueFile| Swath         |
+------------------------------------------------+------------------------+---------------+
