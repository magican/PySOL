#!/usr/bin/env python
""" generated source for module OilSpillDetectionOp """
# 
#  * Copyright (C) 2011 by Array Systems Computing Inc. http://www.array.ca
#  *
#  * This program is free software; you can redistribute it and/or modify it
#  * under the terms of the GNU General Public License as published by the Free
#  * Software Foundation; either version 3 of the License, or (at your option)
#  * any later version.
#  * This program is distributed in the hope that it will be useful, but WITHOUT
#  * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
#  * more details.
#  *
#  * You should have received a copy of the GNU General Public License along
#  * with this program; if not, see http://www.gnu.org/licenses/
#  
# package: org.esa.nest.gpf.oceantools
import com.bc.ceres.core.ProgressMonitor

import org.esa.beam.framework.datamodel.Band

import org.esa.beam.framework.datamodel.Mask

import org.esa.beam.framework.datamodel.Product

import org.esa.beam.framework.datamodel.ProductData

import org.esa.beam.framework.gpf.Operator

import org.esa.beam.framework.gpf.OperatorException

import org.esa.beam.framework.gpf.OperatorSpi

import org.esa.beam.framework.gpf.Tile

import org.esa.beam.framework.gpf.annotations.OperatorMetadata

import org.esa.beam.framework.gpf.annotations.Parameter

import org.esa.beam.framework.gpf.annotations.SourceProduct

import org.esa.beam.framework.gpf.annotations.TargetProduct

import org.esa.beam.util.ProductUtils

import org.esa.nest.datamodel.Unit

import org.esa.nest.gpf.OperatorUtils

import org.esa.nest.gpf.TileIndex

import java.awt

import java.util.ArrayList

import java.util.HashMap

import java.util.List

# 
#  * The oil spill detection operator.
#  *
#  * The algorithm for detecting dark spots is based on adaptive thresholding. The thresholding is based on
#  * an estimate of the typical backscatter level in a large window, and the threshold is set to k decibel
#  * below the estimated local mean backscatter level. Calibrated images are used, and simple speckle filtering
#  * is applied prior to thresholding.
#  *
#  * [1] A. S. Solberg, C. Brekke and R. Solberg, "Algorithms for oil spill detection in Radarsat and ENVISAT
#  * SAR images", Geoscience and Remote Sensing Symposium, 2004. IGARSS '04. Proceedings. 2004 IEEE International,
#  * 20-24 Sept. 2004, page 4909-4912, vol.7.
#  
@OperatorMetadata(alias="Oil-Spill-Detection", category="Ocean-Tools", description="Detect oil spill.")
@SourceProduct(alias="source")
@Parameter(description="The list of source bands.", alias="sourceBands", itemAlias="band", rasterDataNodeType=Band.__class__, label="Source Bands")
@Parameter(description="Background window size", defaultValue="13", label="Background Window Size")
@Parameter(description="Threshold shift from background mean", defaultValue="2.0", label="Threshold Shift (dB)")
class OilSpillDetectionOp(Operator):
    """ generated source for class OilSpillDetectionOp """
    sourceProduct = Product()
    targetProduct = None
    sourceBandNames = None
    backgroundWindowSize = 61
    k = 2.0
    sourceImageWidth = 0
    sourceImageHeight = 0
    halfBackgroundWindowSize = 0
    kInLinearScale = 0.0
    targetBandNameToSourceBandName = HashMap()
    OILSPILLMASK_NAME = "_oil_spill_bit_msk"

    def initialize(self):
        """ generated source for method initialize """
        try:
            getMission()
            self.sourceImageWidth = sourceProduct.getSceneRasterWidth()
            self.sourceImageHeight = sourceProduct.getSceneRasterHeight()
            self.halfBackgroundWindowSize = (backgroundWindowSize - 1) / 2
            if self.k < 0:
                raise OperatorException("Threshold Shift cannot be negative")
            else:
                self.kInLinearScale = Math.pow(10.0, k / 10.0)
            self.targetProduct = Product(sourceProduct.__name__, sourceProduct.getProductType(), sourceImageWidth, sourceImageHeight)
            OperatorUtils.copyProductNodes(self.sourceProduct, self.targetProduct)
            addSelectedBands()
            addBitmasks(self.targetProduct)
        except Throwable as e:
            OperatorUtils.catchOperatorException(getId(), e)

    @classmethod
    def addBitmasks(cls, product):
        """ generated source for method addBitmasks """
        for band in product.getBands():
            if band.__name__.contains(cls.OILSPILLMASK_NAME):
                #   final BitmaskDef mask = new BitmaskDef(band.__name__+"_detection",
                #       "Oil Spill Detection", expression, Color.RED, 0.5f);
                #   product.addBitmaskDef(mask);
                mask.setDescription("Oil Spill Detection")
                mask.getImageConfig().setValue("color", Color.RED)
                mask.getImageConfig().setValue("transparency", 0.5)
                mask.getImageConfig().setValue("expression", expression)
                mask.setNoDataValue(0)
                mask.setNoDataValueUsed(True)
                product.getMaskGroup().add(mask)

    # 
    #      * Get mission from the metadata of the product.
    #      
    def getMission(self):
        """ generated source for method getMission """
    def addSelectedBands(self):
        """ generated source for method addSelectedBands """
        if self.sourceBandNames == None or self.sourceBandNames.length == 0:
            #  if user did not select any band
            for band in bands:
                if band.getUnit() != None and band.getUnit() == Unit.INTENSITY:
                    bandNameList.add(band.__name__)
            self.sourceBandNames = bandNameList.toArray([None]*len(bandNameList))
        sourceBands = [None]*sourceBandNames.length
        i = 0
        while i < sourceBandNames.length:
            if sourceBand == None:
                raise OperatorException("Source band not found: " + sourceBandName)
            sourceBands[i] = sourceBand
            i += 1
        for srcBand in sourceBands:
            if unit == None:
                raise OperatorException("band " + srcBand.__name__ + " requires a unit")
            self.targetBandNameToSourceBandName.put(targetBandName, srcBandNames)
            targetBand.setSourceImage(srcBand.getSourceImage())
            targetBandMask.setNoDataValue(0)
            targetBandMask.setNoDataValueUsed(True)
            targetBandMask.setUnit(Unit.AMPLITUDE)
            self.targetProduct.addBand(targetBandMask)

    # 
    #      * Called by the framework in order to compute a tile for the given target band.
    #      * <p>The default implementation throws a runtime exception with the message "not implemented".</p>
    #      *
    #      * @param targetBand The target band.
    #      * @param targetTile The current tile associated with the target band to be computed.
    #      * @param pm         A progress monitor which should be used to determine computation cancelation requests.
    #      * @throws org.esa.beam.framework.gpf.OperatorException
    #      *          If an error occurs during computation of the target raster.
    #      
    def computeTile(self, targetBand, targetTile, pm):
        """ generated source for method computeTile """
        try:
            while ty < maxy:
                trgIndex.calculateStride(ty)
                srcIndex.calculateStride(ty)
                while tx < maxx:
                    if v == noDataValue:
                        trgData.setElemIntAt(trgIndex.getIndex(tx), 0)
                        continue 
                    if v < threshold:
                        trgData.setElemIntAt(trgIndex.getIndex(tx), 1)
                    else:
                        trgData.setElemIntAt(trgIndex.getIndex(tx), 0)
                    tx += 1
                ty += 1
        except Throwable as e:
            OperatorUtils.catchOperatorException(getId(), e)

    # 
    #      * Compute the mean value for pixels in a given sliding window.
    #      * @param tx The x coordinate of the central point of the sliding window.
    #      * @param ty The y coordinate of the central point of the sliding window.
    #      * @param sourceTile The source image tile.
    #      * @param noDataValue the place holder for no data
    #      * @return The mena value.
    #      
    def computeBackgroundMean(self, tx, ty, sourceTile, noDataValue):
        """ generated source for method computeBackgroundMean """
        x0 = Math.max(tx - self.halfBackgroundWindowSize, 0)
        y0 = Math.max(ty - self.halfBackgroundWindowSize, 0)
        w = Math.min(tx + self.halfBackgroundWindowSize, self.sourceImageWidth - 1) - x0 + 1
        h = Math.min(ty + self.halfBackgroundWindowSize, self.sourceImageHeight - 1) - y0 + 1
        srcData = sourceTile.getDataBuffer()
        tileIndex = TileIndex(sourceTile)
        mean = 0.0
        numPixels = 0
        maxy = y0 + h
        maxx = x0 + w
        y = y0
        while y < maxy:
            tileIndex.calculateStride(y)
            while x < maxx:
                if v != noDataValue:
                    mean += v
                    numPixels += 1
                x += 1
            y += 1
        return mean / numPixels

    # 
    #      * Operator SPI.
    #      
    class Spi(OperatorSpi):
        """ generated source for class Spi """
        def __init__(self):
            """ generated source for method __init__ """
            super(Spi, self).__init__(OilSpillDetectionOp.__class__)

