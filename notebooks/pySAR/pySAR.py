# coding: utf-8

import datetime

__author__   = 'Alexander Myasoedov'
__email__    = 'mag@rshu.ru'
__created__  = datetime.datetime(2014, 10, 28)
__modified__ = datetime.datetime(2016, 5, 5)
__version__  = "1.0"
__status__   = "Development"



from readASAR import *

class pySAR:
    """\
    Initial class for Reading, Interpolating and Resampling data from SAR (ASAR, Sentinel-1) images
    Usage:
        s1 = readS1(iPath, fileName, pxlRes, image_size_skip, lat_crop)
            iPath    = folder
            fileName = file
            pxlRes   = image pixel resolution
            image_size_skip = max size of an image to skip if too large
            area_bbox = processing images overlaying area defined by bounding box
                        in form of [latlim, lonlim] or [latmin, latmax, lonmin, lonmax]
    """
    def __init__(self, iPath, fileName, pxlRes=800.0, image_size_skip=40*1e6, area_bbox=45, units=None, oPath=None, overwrite_nc=False, asar_conf=None, proj='EPSG:3413', numProcs=None):
        """ Reading, Interpolating, Calculating and Resampling data from SAR (ASAR, Sentinel-1) images """

        if fileName.startswith('ASA') and fileName.endswith('.N1'):

            granule_name = fileName[:-3]
            year, month, day, self.startTime = get_date_parameters(granule_name)
            _oPath = os.path.join(oPath, proj.lower().replace(':', '_'), year, month, day)
            oPath_nc = os.path.join(_oPath, granule_name+'.nc')
            mkdirs(_oPath)

            # check if nc-file exists
            if os.path.isfile(oPath_nc) and not overwrite_nc:
                logger.info("NC-file exists, skipping...")
                return
            elif os.path.isfile(oPath_nc) and overwrite_nc:
                try:
                    os.remove(oPath_nc)
                except Exception as e:
                    logger.error("Could not remove NC-file")
                    return

            logger.info("Reading ASAR file %s" %fileName)

            sigma0, lats_2, lons_2, incident_angle, polarization, lats, lons = readASAR(iPath, fileName, pxlRes, image_size_skip, area_bbox)

            self.sigma0 = sigma0
            self.lats_2 = lats_2
            self.lons_2 = lons_2
            self.lats = lats
            self.lons = lons
            self.latlim = (lats_2.min(), lats_2.max())
            self.lonlim = (lons_2.min(), lons_2.max())
            self.incident_angle = incident_angle
            self.polarization   = polarization
            self.pxlRes = pxlRes
            self.units  = units
            self.oPath_nc  = oPath_nc
            self.asar_conf = asar_conf
            self.proj = proj

            if numProcs is None:
                self.numProcs=cpu_count()-1
            else:
                self.numProcs = numProcs

    def computeSARwind(self):
        """ Computing Wind Speed from sigma0 """

        ncepGFSmodelWindSwath = addModelWind(self.startTime, self.lats_2, self.lons_2, self.numProcs)

        # calculate bearing from initial lats/lons for further wind calculation
        # Taking initial values as bearing is more accurate after
        # interpolation than vice versa
        bearing = zeros((self.lons.shape[0]-1, self.lons.shape[1]))
        for n in range(0,self.lons.shape[1]):
            col = ([self.lats[:-1,n], self.lons[:-1,n]], [self.lats[1:,n], self.lons[1:,n]])
            for m in range(0,self.lons.shape[0]-1):
                bearing[m][n] = distancelib.bearing(asarray(col[0])[:,m], asarray(col[1])[:,m])

        # interpolate to raw_counts.shape
        bearing = imresize(bearing, self.sigma0.shape)

        # NB! WINDDIR = 0 WHEN WIND BLOWS TOWARDS RADAR!
        wind_dir_model_swath_rel = 90 + bearing - ncepGFSmodelWindSwath['wind_dir']

        if self.polarization == 'H/H':
            PR = PR_Mouche(self.incident_angle, wind_dir_model_swath_rel)
            try:
                from cmod_gpu import rcs2windOpenCl
                SARwindSpeed = rcs2windOpenCl(sar=self.sigma0*PR,
                                              windir=wind_dir_model_swath_rel,
                                              theta=self.incident_angle)
            except Exception:
                from cmod_vect import rcs2windPar
                SARwindSpeed = rcs2windPar(self.sigma0*PR, cmdv=5,
                                           windir=wind_dir_model_swath_rel,
                                           theta=self.incident_angle,
                                           nprocs=self.numProcs)
        elif self.polarization == 'V/V':
            try:
                from cmod_gpu import rcs2windOpenCl
                SARwindSpeed = rcs2windOpenCl(sar=self.sigma0,
                                                 windir=wind_dir_model_swath_rel,
                                                 theta=self.incident_angle)
            except Exception:
                from cmod_vect import rcs2windPar
                SARwindSpeed = rcs2windPar(self.sigma0, cmdv=5,
                                           windir=wind_dir_model_swath_rel,
                                           theta=self.incident_angle,
                                           nprocs=self.numProcs)
        
        self.SARwindSpeed = SARwindSpeed
        self.ncepGFSmodelWindSwath = ncepGFSmodelWindSwath



    def computeRoughness(self, sigma0, incident_angle, polarization):
        """ Computing Roughness from sigma0 """
        # SARroughness = 


    def resampleSAR(self, var_name=['sigma0']):

        if 'sigma0' in var_name or 'wind_speed' in var_name or 'roughness' in var_name:
            for var in var_name:
                if 'sigma0' in var:
                    self.sigma0_res, self.swath_def, self.area_def = \
                    resample(self.sigma0, self.lats_2, self.lons_2, \
                        latlim=None, lonlim=None, pxlRes=self.pxlRes, proj=self.proj, numProcs=self.numProcs)
                elif 'wind_speed' in var:
                    self.SARwindSpeed_res, self.swath_def, self.area_def = \
                    resample(self.SARwindSpeed, self.lats_2, self.lons_2, \
                        latlim=None, lonlim=None, pxlRes=self.pxlRes, proj=self.proj, numProcs=self.numProcs)
                elif 'roughness' in var:
                    self.SARroughness_res, self.swath_def, self.area_def = \
                    resample(self.SARroughness, self.lats_2, self.lons_2, \
                        latlim=None, lonlim=None, pxlRes=self.pxlRes, proj=self.proj, numProcs=self.numProcs)      


    # Tile Generation
    def generateTiles(self, var_name=['sigma0']):

        logger.info("Creating NC tiles in %s" % self.oPath_nc)

        if 'sigma0' in var_name or 'wind_speed' in var_name or 'roughness' in var_name:
            for var in var_name:
                if 'sigma0' in var:
                    data = self.sigma0_res
                elif 'wind_speed' in var:
                    data = self.SARwindSpeed_res
                elif 'roughness' in var:
                    data = self.SARroughness_res

                data, bbox_area_def, geospatial_extent, geospatial_extent_ll = \
                    ncTilesMetainfo(data, self.lats_2, self.lons_2, \
                                    self.area_def, latlim=None, lonlim=None)
                # set Polarization formatting to comply with standards
                p = self.polarization.lower().replace('/', '')

                max_zoom_level = create_nc_tiles(data, self.oPath_nc, 'u1', var, p, configpath=self.asar_conf, logger=logger)

            write_attrib_to_nc(self.oPath_nc, bbox_area_def.area_extent, self.startTime,
                               max_zoom_level, self.pxlRes, geospatial_extent, geospatial_extent_ll, {p}, var_name)
        else:
            return

        logger.info("Wrote %sMbs" % str(round(os.path.getsize(self.oPath_nc)/1024/1024)))
        # send_redis_message(tiles_4326_r_output_dir, kml_file_4326,
        #                    nc_4326_filename, r)




        # # Reprojecting data
        # # Pixel resolution
        # # we use pxlResWind/pxlResSAR for further pyresample
        # # radius_of_influence and sigmas
        # pxlResWind = asarray(
        #     distancelib.getPixelResolution(
        #         ncepGFSmodelWind['lats_wind'],
        #         ncepGFSmodelWind['lons_wind'],
        #         ncepGFSmodelWind['lons_wind'].shape, 'km'
        #     )
        # )

        # # reproject NCEP onto S1 grid before calculations
        # # Using RectSphereBivariateSpline - Bivariate spline approximation over a rectangular mesh on a sphere
        # # as it is much more efficiant for full resolution
        # # as well as smoothes nicely the image

        # # We don't want to work with full res wind so scaling the image for about 100m resolution
        # # Adjust scale to get appropriate value
        # scale = 1

        # lts = ncepGFSmodelWind['lats_wind']
        # lns = ncepGFSmodelWind['lons_wind']

        # ncepGFSmodelWindSwath = {}
        # # check that lns increasing, if not
        # if ~all(diff(lns[0,:]) > 0):
        #     lns_ = lns[0,:]
        #     # we must start lns from -180 increasing to +180
        #     # so we reconcatenate lns array and data array, by cropping part of array and putting in front
        #     lns_ = concatenate((lns_[(lns_>=-180) & (lns_<0)], lns_[(lns_>=0)]),0)
        #     ncepGFSmodelWindSwath['wind_speed'] = (ncepGFSmodelWind['wind_speed'])
        #     ncepGFSmodelWindSwath['wind_speed'] = concatenate((ncepGFSmodelWind['wind_speed'][:,(lns_>=-180) & (lns_<0)],\
        #                                                        ncepGFSmodelWind['wind_speed'][:,(lns_>=0)]),1)
        #     ncepGFSmodelWindSwath['wind_dir'] = (ncepGFSmodelWind['wind_dir'])
        #     ncepGFSmodelWindSwath['wind_dir'] = concatenate((ncepGFSmodelWind['wind_dir'][:,(lns_>=-180) & (lns_<0)],\
        #                                                      ncepGFSmodelWind['wind_dir'][:,(lns_>=0)]),1)
        #     # only then we do interpolate

        # lts_2 = self.lats_2[::scale,::scale]
        # lns_2 = self.lons_2[::scale,::scale]

        # # RectSphereBivariateSpline uses lats and lons within the intervals (0, pi), (0, 2pi).
        # # Make sure the latitude is between 0 .. 180
        # if lts_2.min()<0:
        #     lts_2 = lts_2 + 90
        #     lts   = lts   + 90
        # if lns_2.min()<0:
        #     lns_2 = lns_2 + 180
        #     lns   = lns   + 180

        # try:
        #     ncepGFSmodelWindSwath['wind_speed'] = griddata((lts.flatten(), lns.flatten()), (ncepGFSmodelWind['wind_speed']).flatten(), (lts_2, lns_2), method='cubic')
        #     ncepGFSmodelWindSwath['wind_dir']   = griddata((lts.flatten(), lns.flatten()), (ncepGFSmodelWind['wind_dir']).flatten(), (lts_2, lns_2), method='cubic')
        #     ncepGFSmodelWindSwath['u'] = griddata((lts.flatten(), lns.flatten()), (ncepGFSmodelWind['u']).flatten(), (lts_2, lns_2), method='cubic')
        #     ncepGFSmodelWindSwath['v'] = griddata((lts.flatten(), lns.flatten()), (ncepGFSmodelWind['v']).flatten(), (lts_2, lns_2), method='cubic')
        # except:
        # #     ncepGFSmodelWindSwath['wind_speed']  = RectBivariateSpline(lts, lns, flipud(ncepGFSmodelWindSwath['wind_speed']), kx=2, ky=2)(lts_2, lns_2)
        #     pass

        # pxlResWindSwath = asarray(distancelib.getPixelResolution(lts_2, \
        #                                                     lns_2, \
        #                                                     lns_2.shape, 'km'))

        # print "Interpolated Wind cell resolution, %s km" % pxlResWindSwath

        # del lts, lns, lts_2, lns_2