from BaseFunctions import *


# READ THE RAW_COUNTS from GRD image

def read_raw_counts(fn, fileLocation, polarization, scale, GEOgrid):


    # Note that For SLC images BitsPerSample=32 and for GRD BitsPerSample=16

    # im = Image.open(inpath + fileLocation['s1aiwgrd' + polarization])

    im = zf.read(fn[:-4] + '.SAFE' + fileLocation[fn.lower().replace("_","")[0:8] + polarization][1:])
    im = StringIO.StringIO(im) #Encode the raw data to be used by Image.open()
    im = Image.open(im)        #Open the image
    im = asarray(im)

    arrShape =  asarray([GEOgrid['numberOfLines'], GEOgrid['numberOfSamples']])
    ext, spa = _format_extent_spacing((0.,0.,arrShape[0]-1,arrShape[1]-1), scale, GEOgrid)
    if scale.max() > 1:
        im = im[ext[0]:ext[2],ext[1]:ext[3]]
        sha = (im.shape[0]/spa[0], spa[0],
               im.shape[1]/spa[1], spa[1])
        im = im.reshape(sha).mean(-1).mean(1)

    return im, ext, spa


# In[12]:

# READ the ANNOTATION

def read_anotation(fn, fileLocation, polarization):
    # open the fileLocation

    # annotation = open(inpath + fileLocation['products1aiwgrd' + polarization], "r") # Open a file in read-only mode
    # annotation = annotation.read() # read the file object

    annotation = zf.read(fn[:-4] + '.SAFE' + fileLocation['product' + fn.lower().replace("_","")[0:8] + polarization][1:])
    annotation = xmltodict.parse(annotation) # Parse the read document string

    # get geolocationGrid parameters from the Annotation Data Set Records (ADSR)
    # preallocate variables
    GEOgrid = {} # empty dict for GEOgrids

    GEOgrid['rangePixelSpacing']   = float(annotation['product']['imageAnnotation']['imageInformation']['rangePixelSpacing'])
    GEOgrid['azimuthPixelSpacing'] = float(annotation['product']['imageAnnotation']['imageInformation']['azimuthPixelSpacing'])
    GEOgrid['numberOfSamples']   = float(annotation['product']['imageAnnotation']['imageInformation']['numberOfSamples'])
    GEOgrid['numberOfLines'] = float(annotation['product']['imageAnnotation']['imageInformation']['numberOfLines'])

    geolocationGridPointList = annotation['product']['geolocationGrid']['geolocationGridPointList']
    GEOgrid['lats']  = zeros( ( int(geolocationGridPointList['@count']), 1) )
    GEOgrid['lons']  = zeros( GEOgrid['lats'].shape )
    GEOgrid['line']  = zeros( GEOgrid['lats'].shape, dtype=int )
    GEOgrid['pixel'] = zeros( GEOgrid['lats'].shape, dtype=int )
    GEOgrid['incidenceAngle'] = zeros( GEOgrid['lats'].shape )
    GEOgrid['elevationAngle'] = zeros( GEOgrid['lats'].shape )

    # read Geolocation grid points
    for n in range(int(geolocationGridPointList['@count'])):
        GEOgrid['lats'][n]  = float(geolocationGridPointList['geolocationGridPoint'][n]['latitude'])
        GEOgrid['lons'][n]  = float(geolocationGridPointList['geolocationGridPoint'][n]['longitude'])
        GEOgrid['line'][n]  = int(geolocationGridPointList['geolocationGridPoint'][n]['line'])
        GEOgrid['pixel'][n] = int(geolocationGridPointList['geolocationGridPoint'][n]['pixel'])
        GEOgrid['incidenceAngle'][n] = float(geolocationGridPointList['geolocationGridPoint'][n]['incidenceAngle'])
        GEOgrid['elevationAngle'][n] = float(geolocationGridPointList['geolocationGridPoint'][n]['elevationAngle'])


    # find zero pixel to reshape grid points to array
    ind = find(GEOgrid['pixel'] == 0)
    GEOgrid['pixel'] = reshape(GEOgrid['pixel'], (ind.size, GEOgrid['lats'].size/ind.size))
    GEOgrid['line']  = reshape(GEOgrid['line'], (ind.size, GEOgrid['lats'].size/ind.size))
    GEOgrid['lats']  = reshape(GEOgrid['lats'], (ind.size, GEOgrid['lats'].size/ind.size))
    GEOgrid['lons']  = reshape(GEOgrid['lons'], (ind.size, GEOgrid['lats'].size/ind.size))
    GEOgrid['incidenceAngle'] = reshape(GEOgrid['incidenceAngle'], (ind.size, GEOgrid['lats'].size/ind.size))
    GEOgrid['elevationAngle'] = reshape(GEOgrid['elevationAngle'], (ind.size, GEOgrid['lats'].size/ind.size))

    return GEOgrid


# In[13]:

# READ the CALIBRATION LUTs

def read_clbrtn_luts(fn, fileLocation, polarization):
    # The calibration data set contains calibration information
    # and the beta nought, sigma nought, gamma and digital
    # number (DN) Look-up Tables (LUT)s that can be used for
    # absolute product calibration.

    # We take only the calibrationVector record
    # This record holds the calibration vectors and associated fields required to
    # derive radiometrically calibrated imagery from the image MDS.

    # open the fileLocation

    # calibration = open(inpath + fileLocation['calibrations1aiwgrd' + polarization], "r") # Open a file in read-only mode
    # calibration = calibration.read() # read the file object

    calibration = zf.read(fn[:-4] + '.SAFE' + fileLocation['calibration' + fn.lower().replace("_","")[0:8] + polarization][1:])
    calibration = xmltodict.parse(calibration) # Parse the read document string

    calibrationVectorList = calibration['calibration']['calibrationVectorList']

    cLUTs = {} # empty dict for CALIBRATION LUTs

    cLUTs['line']  = zeros( ( int(calibrationVectorList['@count']), 1), dtype=int )
    cLUTs['pixel'] = zeros( (cLUTs['line'] .shape[0], int(calibrationVectorList['calibrationVector'][0]['pixel']['@count'])), dtype=int )
    cLUTs['sigmaNought']  = zeros( cLUTs['pixel'].shape, dtype=float)

    # read Calibration Vector points
    for n in range(int(calibrationVectorList['@count'])):
        cLUTs['line'][n] = int(calibrationVectorList['calibrationVector'][n]['line'])
        pixel_     = calibrationVectorList['calibrationVector'][0]['pixel']['#text']
        cLUTs['pixel'][n,:] = asarray(pixel_.split(' '), dtype=int)
        sigmaNought_     = calibrationVectorList['calibrationVector'][0]['sigmaNought']['#text']
        cLUTs['sigmaNought'][n,:] = asarray(sigmaNought_.split(' '), dtype=float)

    return cLUTs


# In[14]:

# READ the NOISE LUTs

def read_noise_luts(fn, fileLocation, polarization):
    # The L1 Noise ADS provides a LUT – with values provided in linear power –
    # that can be used to derive calibrated noise profiles which match the calibrated GRD data.

    # open the fileLocation

    noise = zf.read(fn[:-4] + '.SAFE' + fileLocation['noise' + fn.lower().replace("_","")[0:8] + polarization][1:])
    noise = xmltodict.parse(noise) # Parse the read document string

    noiseVectorList = noise['noise']['noiseVectorList']

    nLUTs = {} # empty dict for NOISE LUTs

    nLUTs['line']  = zeros( ( int(noiseVectorList['@count']), 1), dtype=int )
    nLUTs['pixel'] = zeros( (nLUTs['line'] .shape[0], int(noiseVectorList['noiseVector'][0]['pixel']['@count'])), dtype=int )
    nLUTs['noiseLut']  = zeros( nLUTs['pixel'].shape, dtype=float)

    # read Calibration Vector points
    for n in range(int(noiseVectorList['@count'])):
        nLUTs['line'][n] = int(noiseVectorList['noiseVector'][n]['line'])
        pixel_     = noiseVectorList['noiseVector'][0]['pixel']['#text']
        nLUTs['pixel'][n,:] = asarray(pixel_.split(' '), dtype=int)
        noiseLut_     = noiseVectorList['noiseVector'][0]['noiseLut']['#text']
        nLUTs['noiseLut'][n,:] = asarray(noiseLut_.split(' '), dtype=float)

    return nLUTs

def readS1_raw_counts(inpath = '/media/SOLabNFS2/tmp/sentinel-1/Svalbard-Barents/', \
fn='S1A_EW_GRDH_1SDH_20141003T133957_20141003T134057_002666_002F83_0BE7.zip', resolution=None):
    """Reading raw_counts from Sentinel-1 images"""

    global zf
    zf = zipfile.ZipFile(inpath+fn, 'r')

    manifest = zf.read(fn[:-4] + '.SAFE/manifest.safe')
    manifest = xmltodict.parse(manifest) # Parse the read document string

    # we have different XML/GeoTIFF fileLocations for vv-vh (or hh-hv) products, respectively
    # create new dict with file paths to the fileLocations
    fileLocation = {} # empty dict
    dataObject = manifest['xfdu:XFDU']['dataObjectSection']['dataObject']
    for n in range(len(dataObject)):
        if len(dataObject[n]['@ID']) > 45:
            k = dataObject[n]['@ID'][0:-45] # get the new key from @ID
        else:
            k = dataObject[n]['@ID'] # get the new key from @ID
        v = str(dataObject[n]['byteStream']['fileLocation']['@href']) # get the dict.value
        fileLocation[k] = v # assign to new dict
    #     locals()['fileLocation_'+k]=v # create local variable from 'fileLocation' dict.key


    # In[17]:

    # Get the productType/polarization

    metadataObject = manifest['xfdu:XFDU']['metadataSection']['metadataObject']
    for n in range(len(metadataObject)):
        if metadataObject[n]['@ID'] == 'generalProductInformation':
            productType = metadataObject[n]['metadataWrap']['xmlData']['s1sarl1:standAloneProductInformation']['s1sarl1:productType']
            transmitterReceiverPolarisation = metadataObject[n]['metadataWrap']['xmlData']['s1sarl1:standAloneProductInformation']['s1sarl1:transmitterReceiverPolarisation']

    polarization = []
    for p in range(0,len(transmitterReceiverPolarisation[0].lower())):
        if len(transmitterReceiverPolarisation[0].lower()) == 1:
            polarization = transmitterReceiverPolarisation.lower()
        else:
            polarization.append(transmitterReceiverPolarisation[p].lower())
    if isinstance(polarization, basestring):
        polarization = [polarization]
    print "Available polarizations: \'%s\'" %polarization

    # In[18]:

    # %%timeit -n 1 -r 1

    # set default resolution
    if resolution is None:
        resolution = 80

    raw_counts = {}

    GEOgrid = read_anotation(fn, fileLocation, polarization[0])

    # Find scale to reduce image to the specified resolution
    # arrShape =  asarray([GEOgrid['numberOfLines'], GEOgrid['numberOfSamples']])
    # scale = resolution/round(mean(asarray(distancelib.getPixelResolution(GEOgrid['lats'], \
    #                                                                      GEOgrid['lons'], \
    #                                                                      arrShape, 'km'))*1e3))
    scale = floor(resolution/asarray([GEOgrid['azimuthPixelSpacing'],GEOgrid['rangePixelSpacing']]))

    for p in polarization:
        print "Reading raw_counts: \'%s\' polarization" %p
        # READ THE RAW_COUNTS from GRD image
        try:
            raw_counts[p], _, _ = read_raw_counts(fn, fileLocation, p, scale, GEOgrid)
        except:
            print "Can't read file: ", fn[:-4] + '.SAFE' + fileLocation[fn.lower().replace("_","")[0:8] + p][1:]
            polarization.remove(p)
            pass

    # Close ZIP-file
    zf.close()

    return raw_counts, polarization, manifest, GEOgrid