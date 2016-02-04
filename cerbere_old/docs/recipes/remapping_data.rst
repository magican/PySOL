Remapping and resampling data
=============================

These examples use the :class:`Resampler` class.

Resampling gridded data on another grid
----------------------------------------



A generic resampler
-------------------

    srcmapper = srcmapperclass(url=srcproduct)
    targetmapper = targetmapperclass(url=targetproduct)
    srcfeature = srcmodelclass()
    targetfeature = targetmodelclass()
    srcfeature.load(srcmapper)
    targetfeature.load(targetmapper)
    # resampling
    if fields is None:
        fields = srcfeature.get_fieldnames()
    if radius is None:
        radius = srcfeature.get_spatial_resolution()
    print "RADIUS ", radius
    outfeature = Resampler.resample(
                        srcfeature,
                        targetfeature,
                        fields=fields,
                        radius=radius,
                        new=True,
                        add_reference=True
                        )
    # metadata
    if outfeature:
        # set global attributes
        globalattrs = {}
        source_times = outfeature.get_field('resampled_time')
        timevals = source_times.get_values()
        time0 = timevals.min()
        time1 = timevals.max()
        mintime = num2date(time0,
                                   source_times.units)
        maxtime = num2date(time1,
                                   source_times.units)
        globalattrs['time_coverage_start']\
             = mintime.strftime('%Y%m%dT%H%M%S')
        globalattrs['time_coverage_stop']\
             = maxtime.strftime('%Y%m%dT%H%M%S')
        globalattrs['spatial_colocation_radius_in_km'] = radius / 1000.
        # save file
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        fname = os.path.join(
                    outdir,
                    'resampled_' + os.path.basename(srcproduct))
        if os.path.exists(fname):
            os.remove(fname)
        outfile = ncfile.NCFile(fname, abstractmapper.WRITE_NEW, ncformat="NETCDF4")
        outfeature.save(outfile, attrs=globalattrs)
        outfile.close()
    srcmapper.close()
    targetmapper.close()