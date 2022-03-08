# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 16:31:20 2021


@author: chao wang

https://code.earthengine.google.com/e78c9005738ede0559c53a5eefda6b06
In JavaScript API, it is not easy to automatic running the task one by one
this is why I convert the original javascript into python api

Updated at Nov.4, 2021
https://code.earthengine.google.com/f84d61f1456bcb4cd528b9b8bfeaa90c?noload=1

"""

import os
import sys
import ee
import math
import re

ee.Initialize()


###===========================================================
"""
Function to rename bands
"""
def getNewBandNames(prefix,bandNames_):
    #create a list with the same number of element
    seq = ee.List.sequence(1, bandNames_.length())
    
    #return the renamed band elements as list
    return seq.map(lambda bandNum: ee.String(prefix).cat(ee.Number(bandNum).int().format("%d")))
######----End Rename Bands Function--------------------------------



###===========================================================
"""
Function to calculate Principal Components
Inputs: AVIRIS-NG bands, Study Area, and Given bandNames
Output: means of AVIRIS-NG, eigenVectors, and eigenValues
"""
def CalculatePCA(InputImgCols,AOI,bandNames_,scale_define):

    #get composite image by median summary
    BRDF_ImgMedian=InputImgCols.median()

    #remove zeros
    NonZeros = BRDF_ImgMedian.gte(0).Or(BRDF_ImgMedian.lte(0))
    # Map.addLayer(NonZeros,{min:0,max:1},'mask', false)

    #update mask
    BRDF_ImgMedian = BRDF_ImgMedian.updateMask(NonZeros)

    # Mean center the data to enable a faster covariance reducer
    # and an SD stretch of the principal components.
    BRDF_meanDict = BRDF_ImgMedian.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=AOI,
        scale=scale_define,
        maxPixels=10**13,
        bestEffort=True,
        tileScale=16)

    #build the constant band image
    means = ee.Image.constant(BRDF_meanDict.values(bandNames_))

    #perform principal components with centered covariance
    BRDF_centered = BRDF_ImgMedian.subtract(means)

    # Collapse the bands of the image into a 1D array per pixel.
    BRDF_arrays = BRDF_centered.toArray()

    # Compute the covariance of the bands within the region.
    covariance = BRDF_arrays.reduceRegion(
        reducer=ee.Reducer.centeredCovariance(),
        geometry=AOI,
        scale=scale_define,
        maxPixels=10**13,
        bestEffort=True,
        tileScale=16)
    # print("covariates", covariance)

    # Get the 'array' covariance result and cast to an array.
    # This represents the band-to-band covariance within the region.
    covarArray = ee.Array(covariance.get('array'))

    # Perform an eigen analysis and slice apart the values and vectors.
    eigens = covarArray.eigen() 

    # This is a P-length vector of Eigenvalues.
    eigenValues = eigens.slice(1, 0, 1)

    # This is a PxP matrix with eigenvectors in rows.
    eigenVectors = eigens.slice(1, 1)
    
    #By performing some algebra, the proportion of variance explained (PVE) 
    #by the mth principal component is calculated using the equation:
    #It can be shown that the PVE of the mth principal component can be more 
    #simply calculated by taking the mth eigenvalue and dividing it by 
    #the number of principal components, p.
    ###PVE <- arrests.eigen$values / sum(arrests.eigen$values)
    PVE=eigenValues.divide(eigenValues.reduce(ee.Reducer.sum(), [0]).get([0,0]))

    # Convert the array image to 2D arrays for matrix computations.
    # arrayImage = arrays.toArray(1)

    # Left multiply the image array by the matrix of eigenvectors.
    # principalComponents = ee.Image(eigenVectors).matrixMultiply(arrayImage)

    # # Turn the square roots of the Eigenvalues into a P-band image.
    # sdImage = ee.Image(eigenValues.sqrt()) 
    #   .arrayProject([0]).arrayFlatten([getNewBandNames('sd')])

    # # Turn the PCs into a P-band image, normalized by SD.
    # vdpca = principalComponents
    #   # Throw out an an unneeded dimension, [[]] -> [].
    #   .arrayProject([0])
    #   # Make the one band array image a multi-band image, [] -> image.
    #   .arrayFlatten([getNewBandNames('pc')])
    #   # Normalize the PCs by their SDs.
    #   .divide(sdImage)

    return [means,ee.Image(eigenVectors),eigenValues,PVE]
######----End PCA Calculation Function--------------------------------


######=====================================================
"""
function to apply calculated Principal Components model
Inputs: each AVIRIS-NG sub-scene
Output: each sub-scene PCA
"""
def ApplyPCA2EachScene(PCA_cal_parameters):
    
    def applyPCAFun(EachSceneImg): 
        BandNameList=ee.Image(EachSceneImg).bandNames()
        #get the mean of AVIRIS-NG mosaic
        MosaicMeans=ee.Image(ee.List(PCA_cal_parameters).get(0))
    
        #Mean center the data
        EachSceneImgcentered = EachSceneImg.subtract(MosaicMeans)
    
        #Collapse the bands of the image into a 1D array per pixel
        EachScenearrays = EachSceneImgcentered.toArray()
    
        #Convert the array image to 2D arrays for matrix computations
        EachScenearrayImage = EachScenearrays.toArray(1)
    
        #apply PCA model by multply eigenVectors(matrix computations)
        EachSceneprincipalComponents = ee.Image(ee.List(PCA_cal_parameters).get(1))\
        .matrixMultiply(EachScenearrayImage)
    
        # Turn the square roots of the Eigenvalues into a P-band image.
        sdImage = ee.Image(ee.Array(ee.List(PCA_cal_parameters).get(2)).sqrt())\
        .arrayProject([0]).arrayFlatten([getNewBandNames('sd',BandNameList)])
    
        # Turn the PCs into a P-band image, normalized by SD.
        # Make the one band array image a multi-band image, [] -> image.
        EachScenePCA = EachSceneprincipalComponents.arrayProject([0])\
        .arrayFlatten([getNewBandNames('pc',BandNameList)])\
        .divide(sdImage)
    
        return EachScenePCA.copyProperties(ee.Image(EachSceneImg))
    return applyPCAFun
######----End apply PCA Function--------------------------------

#-------------------------------------------------------------------*/
##########--------End Section 1 Processing Functions----------##########
######################################################################



######=====================================================
#-------------------------------------------------------------------*/
# Section 2: Conduct PCA variable calculation
#-------------------------------------------------------------------*/

#define a Area of Interest
#statsArea = ee.Geometry.Polygon(
#        [[[-111.66100613499445, 59.04774158542581],
#          [-111.66100613499445, 58.32980878737849],
#          [-111.32592312718195, 58.32980878737849],
#          [-111.32592312718195, 59.04774158542581]]])

#statsArea = ee.Geometry.Polygon(
#        [[[-111.67081239730796, 58.995148107827866],
#          [-111.67081239730796, 58.341051440791716],
#          [-111.23135927230796, 58.341051440791716],
#          [-111.23135927230796, 58.995148107827866]]], None, False)

statsArea = ee.Geometry.Polygon(
        [[[-146.34130137143882, 66.0484878113757],
          [-146.22250757768853, 66.05715591070538],
          [-146.21528995004914, 66.04669334133754],
          [-146.15452555152902, 66.04848780792197],
          [-146.20548573465499, 66.23922433355017],
          [-146.25791441082376, 66.42997824347249],
          [-146.4515576839868, 66.4176189170355],
          [-146.39739955265156, 66.23313232814529]]], None, False)

#load BRDF corrected AVIRIS_NG Ref cols
#AVIRIS_NG_S2_Cols=ee.ImageCollection("projects/GlobalReservoirs/AVIRIS_NG_WS/AVIRIS_NG_RFL_S3_pyMR")\
#AVIRIS_NG_BRDF_Cols=ee.ImageCollection("projects/GlobalReservoirs/AVIRIS_NG_WS/AVIRIS_NG_RFL_S3_pyMR")

#AVIRIS_NG_RadNormImgCols = ee.ImageCollection("projects/GlobalReservoirs/AVIRIS_NG_WS/AVIRIS_NG_S4_RadNorm")

#For Yukon Flat, since it is not a large area, three flightlines looks seemless after BRDF
AVIRIS_NG_RadNormImgCols = ee.ImageCollection("projects/proj2019/assets/AVIRIS_NG_RFL_S3")

#######################################
#Get Unique Img Cols by distinct function
UniqueImgNameList=AVIRIS_NG_RadNormImgCols.toList(AVIRIS_NG_RadNormImgCols.size())\
.map(lambda img : ee.Image(img).get('name')).distinct()
print('UniqueImgNameList',UniqueImgNameList.getInfo())

#reduce image cols to unque image cols based on median reduce
def MedianCompositeSubFun(FileName):
    ImgWithSameNameCols=AVIRIS_NG_RadNormImgCols.filterMetadata('name','contains',FileName)
    MedianCompositeImg=ImgWithSameNameCols.median()
    return MedianCompositeImg.set({'name':FileName})

AVIRIS_NG_RadNormUniqueImgCols=ee.ImageCollection(UniqueImgNameList.map(MedianCompositeSubFun))
print('AVIRIS_NG_RadNormUniqueImgCols',AVIRIS_NG_RadNormUniqueImgCols.size().getInfo())

# Set some information about the input to be used later.

#set the output root path
#AVIRIS_NG_WS/AVIRIS_NG_RFL_S3_pyMR
#projects/GlobalReservoirs/AVIRIS_NG_WS/AVIRIS_NG_RFL_S4_MR_PCA
#exportPathRoot='projects/GlobalReservoirs/AVIRIS_NG_WS/AVIRIS_NG_S5_PCA'
exportPathRoot='projects/proj2019/assets/AVIRIS_NG_RFL_S4_PCA'
    
#get output projection
RasterProjection=ee.Image(AVIRIS_NG_RadNormImgCols.first()).projection()

#get the scale
scale=RasterProjection.nominalScale()
#    print('scale',scale.getInfo())
#scale=5

#get the crs
crs = RasterProjection.crs()

#spectral bands list
bandNames = ee.Image(AVIRIS_NG_RadNormUniqueImgCols.first()).bandNames()#.getInfo()
#bandNames=ee.List(["b16","b36","b56","b97","b241","b361"])
#bandList = bandNames.getInfo()
PCA_List=getNewBandNames('pc',bandNames)
#print('PCA_List',PCA_List.getInfo())

#define final band list to export
bandSeq = ee.List.sequence(0, 50)
   

#-------------------------------------------------------------------*/
##########--------End Parameters Preparation----------##########

#-------------------------------------------------------------------*/
# Section 2: calculate PCA parameters and Export Explained variance percentage
#-------------------------------------------------------------------*/

#calculate PCA parameters
PCA_parameters=CalculatePCA(AVIRIS_NG_RadNormUniqueImgCols,statsArea,bandNames,scale)
#print('PCA_parameters',ee.List(PCA_parameters).get(3).getInfo())

##############################
#Eigen Values
#Make a feature without geometry and set the properties to the dictionary of means.
EigenValues=ee.Array(ee.List(PCA_parameters).get(2))

#Wrap the Feature in a FeatureCollection for export.
EigenValuesFeaCols=ee.FeatureCollection(\
        ee.List.sequence(0,EigenValues.length().get([0]).subtract(1))\
        .map(lambda seqOrder : ee.Feature(None,{'EigenValues' : ee.Array(ee.List(EigenValues)).get([seqOrder,0])})))


EigenValueTableTask = ee.batch.Export.table.toDrive(
    collection=EigenValuesFeaCols,
    description='EigenValuesFeaCols',
    folder='SolarPV',
    fileNamePrefix='YukonFlatEigenValueFeaColsRadNorm')

#EigenValueTableTask.start()


##############################
#Explained Of Variance
#Make a feature without geometry and set the properties to the dictionary of means.
ExplainedOfVariance=ee.Array(ee.List(PCA_parameters).get(3))

#Wrap the Feature in a FeatureCollection for export.
ExplainedOfVarianceFeaCols=ee.FeatureCollection(\
        ee.List.sequence(0,ExplainedOfVariance.length().get([0]).subtract(1))\
        .map(lambda seqOrder : ee.Feature(None,{'PVE' : ee.Array(ee.List(ExplainedOfVariance)).get([seqOrder,0])})))
#print('ExplainedOfVarianceFeaCols',ExplainedOfVarianceFeaCols)


PVETableTask = ee.batch.Export.table.toDrive(
    collection=ExplainedOfVarianceFeaCols,
    description='ExplainedOfVarianceFeaCols',
    folder='SolarPV',
    fileNamePrefix='YukonFlatExplainedOfVarianceFeaColsRadNorm')

#PVETableTask.start()

######################################################################

#-------------------------------------------------------------------*/
# Section 3: Apply PCA calculation and Export PCA of each scene
#-------------------------------------------------------------------*/

    
#define a list to conduct calculation looply

#TargetFlightLineNameList=['ang20190715t185057'] 

#loop to export the results
for TargetFlightLineName in UniqueImgNameList.getInfo():
    print('TargetFlightLineName',TargetFlightLineName)
    
    #Apply PCA calculation to each scene
    AVIRIS_NG_PCA_Cols=AVIRIS_NG_RadNormImgCols\
        .filterMetadata('name','contains',TargetFlightLineName)\
        .map(ApplyPCA2EachScene(PCA_parameters))


    
    #filter to get the target chuck image
    ChunkRefRasterImg=ee.Image(AVIRIS_NG_PCA_Cols.\
            filterMetadata('name','contains',TargetFlightLineName).first())
    
    
    #filter to get the target chuck image
    RasterImgForPolygon=ee.Image(AVIRIS_NG_RadNormImgCols.\
            filterMetadata('name','contains',TargetFlightLineName).first()).\
            get('system:footprint')
    
    ##get the footprint of current image which is a LinearRing object
    #and extract its coordinates
    FootprintCoordinates=ee.Geometry(RasterImgForPolygon).coordinates()
#    print(FootprintCoordinates)

    #reconstruct a polygon as a study area based on the coordinates
    exportFootprintPolygon=ee.Algorithms.GeometryConstructors.Polygon(FootprintCoordinates)
    
    #set a output image filename
    ChunkFileName=TargetFlightLineName+'_radnorm_pca'
    
    #set asset id
    assetPath =exportPathRoot + '/' + ChunkFileName
    
    #select first 20 PCs
    #'EPSG:32606' Yunkon
    #'EPSG:32612' PAD
    ChunkRefRasterImgSelected=ChunkRefRasterImg\
            .reproject(ee.Projection('EPSG:32606').atScale(scale))\
            .float().select(bandSeq)
    
    taskIA2Assets_config=dict(assetId=assetPath,
        pyramidingPolicy= {'.default': 'mean'},
        region=exportFootprintPolygon,
        scale=scale,
        crs=crs,
        maxPixels=10**13)
            
    #call export to Asset API
    taskIA2assets = ee.batch.Export.image.toAsset(
        image=ChunkRefRasterImgSelected,
        description=ChunkFileName,**taskIA2Assets_config)
        
    #begin task "export image"
    taskIA2assets.start()
    


#-------------------------------------------------------------------*/
##########--------End Section 4 Export flightline----------##########
######################################################################
#After calculate PCA of each sub-scene  
#apply after get principle components of each scene
#Then apply normalized based on min-max method
#PCA_median = AVIRIS_NG_PCA_Cols.select(PCA_List).median()#.abs()
#print(PCA_median)
#
##rescale 0 to 100
#minPCA =PCA_median.reduceRegion(reducer=ee.Reducer.min(),\
#    geometry=statsArea,scale=5,bestEffort=True,tileScale=16)
#
#minPCA = ee.Dictionary(minPCA).toImage()
#
#maxPCA =PCA_median.reduceRegion(reducer=ee.Reducer.max(),\
#    geometry=statsArea,scale=5,bestEffort=True,tileScale=16)
#
#maxPCA = ee.Dictionary(maxPCA).toImage()
#
##normalized PCA by min and max method
#normPCA = PCA_median.subtract(minPCA).divide(maxPCA.subtract(minPCA))
#
##multiply 100 to scale it up
#PCA_scaled = normPCA.multiply(100)




