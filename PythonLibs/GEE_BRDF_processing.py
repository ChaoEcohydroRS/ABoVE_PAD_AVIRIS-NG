# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 16:45:07 2020
revised on Fri Feb 12 15:27:00 2021
@author: wayne
"""

import os
import sys
import ee
import math
import re

ee.Initialize()

GlobalROI = ee.Geometry.Polygon(
        [[[-111.75831480266295, 59.048459180079675],
          [-111.72576124469026, 58.67736988288882],
          [-111.6920511498481, 58.29974326145159],
          [-111.3569681420356, 58.30984464963926],
          [-111.38346813653962, 58.68530860639207],
          [-111.41181395197212, 59.058648183627376]]])

#-------------------------------------------------------------------*/
# Section 1: Loading Processing Functions
#-------------------------------------------------------------------*/
#ImgColsProcessingPkgs = require('users/Waynechao128/default:Modules/HSI/ProcessingLibs') 

#load the wavelength array
#wavelength=ImgColsProcessingPkgs.wavelength


############################################
#define function to combine footprints of each image chunk
def combChunkGeometryFun(currentImg,prePolygon):
  #because the get footprints will return a LinearRing geometry object
  #but I need a polygon as a study area before conducting reducing operation

  #get the footprint of current image which is a LinearRing object
  #and extract its coordinates
  currentCoordinates=ee.Geometry(currentImg.get('system:footprint')).coordinates()
  
  #reconstruct a polygon as a study area based on the coordinates
  currentPolygon=ee.Algorithms.GeometryConstructors.Polygon(currentCoordinates)
  
  #combine it with the previous polygon to 
  #get a footprint of one flightline since the imagery was splitted
  #into many small chunks to facility the uploading and fine-register processing
  return currentPolygon.union(prePolygon,200)
####End combFlightLineS3Fun function ################################
  

#############################################
#construct rfl image collection by combining all small image chunks 
#of each flightline into one image
def combFlightLineS3Fun(AVIRIS_NG_Cols):
    #define map function
    def CombineImgChunk(rflImg):
      
      # filter to get all chunks of each flightline
      eachFlightLineAVIRIS_NG_Cols=AVIRIS_NG_Cols.filterMetadata('name','equals',rflImg.get('name'))
      
      rflImg_median= eachFlightLineAVIRIS_NG_Cols.median()
      
      #get the geometry of the first chunk
      firstImgGeo=ee.Image(eachFlightLineAVIRIS_NG_Cols.first()).get('system:footprint')
      firstCoordinates=ee.Geometry(firstImgGeo).coordinates()
      firstPolygon=ee.Algorithms.GeometryConstructors.Polygon(firstCoordinates)
      
      #combine the footprint of all image chunk of the same flightline  
      ImgColsFootprint=eachFlightLineAVIRIS_NG_Cols.iterate(combChunkGeometryFun,firstPolygon)
            
      #add the specific rdn data into the rfl image
      return rflImg_median.set({'Fulfootprint': ImgColsFootprint,\
        'name': rflImg.get('name')})
    
    return CombineImgChunk
####End combFlightLineS3Fun function ################################
  

################################################
#Split one flightline back to original image chunks by
#by masking the one flightline with each chunk of the subimage (step 2)
def SplitBRDFCorrImgChunkFun(BRDFCorrOneFlightLineCols):
    #define map function
    def SplitImgChunk(rflImg):
      
        targetFlightLineImg=ee.Image(BRDFCorrOneFlightLineCols
        .filterMetadata('name','equals',rflImg.get('name')).first())
        
        
        #get the footprint of current image which is a LinearRing object
        #and extract its coordinates
        SubImgCoordinates=ee.Geometry(rflImg.get('system:footprint')).coordinates()
        
        #reconstruct a polygon as a study area based on the coordinates
        SubImgPolygon=ee.Algorithms.GeometryConstructors.Polygon(SubImgCoordinates)
        
        #apply clipping with the sub image extent
        SplitChunkImg=targetFlightLineImg.clip(SubImgPolygon)#.updateMask(rflImg.select([0]).mask())#

        return SplitChunkImg#.copyProperties(rflImg)
  
    return SplitImgChunk
####End SplitBRDFCorrImgChunkFun function ################################
    

################################################
##Apply Multiregression BRDF Correction to an image
def MultiRegressionBRDFAdditiveCorr(LC_MaskCols__,NumRemapClass__,imgRFL, imgRdn, scale, studyArea, bandList):
  # Parameters 
  # ----------
  # imgRFL:     reflectance images
  # imgRdn:     Orthocorrected observation parameters
  #scale:       image resoluion, here 5 meters
  #StudyArea:   ROI, here image extent 
  #bandList:    List of the target bands' name
  
  # Returns
  # -------
  # Corrected Image    

    
  #set the band list if not defined
  # if(bandList === null || bandList === undefined){
  #   bandList = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'] 
  # }
  
  #['flightHeight','sensorAzimuth','sensorZenith','sunAzimuth','sunZenith',
  # 'slope','aspect','CloudMask','ShadowMask']
  # Extract solar zenith and azimuth bands
  # to-sun-zenith (0 to 90 degrees from zenith)
  solar_zn = imgRdn.select(['sunZenith']).multiply(math.pi).divide(180)
  
  # to-sun-azimuth (0 to 360 degrees clockwise from N)
  solar_az = imgRdn.select(['sunAzimuth']).multiply(math.pi).divide(180)
  
  #to-sensor-zenith (0 to 90 degrees from zenith)
  sensor_zn = ee.Image(imgRdn.select(['sensorZenith'])).multiply(math.pi).divide(180)
  
  # #to-sensor-zenith (0 to 90 degrees from zenith)
  # #revised to -18 to 18
  # sensor_zn = ee.Image(imgRdn.select(['VZA']))
  # .multiply(math.pi).divide(180)
  
  #to-sensor-azimuth (0 to 360 degrees clockwise from N)
  #sensor_az = imgRdn.select(['sensorAzimuth']).multiply(math.pi).divide(180).cos()
  
  # ImgOrder=imgRdn.select('ImgOrder')
  

  # relative_az = sensor_az - solar_az 
  # relative_az =solar_az
  # solar_az.subtract(sensor_az)
  
  
  #Add angles to ref image as a single image
  #because the linear regression reducer (below) expects to see images
  # with all the independent values followed by all the dependent values.
  #ρUncorr=BRDF+ρAVG=f_vza*(VZA) + f_raa*(RAA) + f_sza*(SZA) + f_constant
  img_plus_angles=ee.Image([ee.Image(1).float().rename('constant'), #constant
      # ImgOrder.rename('ImgOrder'),
      sensor_zn.rename('VZA'),
      solar_az.rename('SAA'),
      # sensor_az.rename('VAA'),
      solar_zn.rename('SZA'),
      imgRFL])
  
  #cloud and shadow buffer 100 meters to get rid of the residues
  CloudShadowKernel=ee.Kernel.euclidean(100, 'meters', True)
  
  #no correction applying over the masked area  
  CloudShadowMask=imgRdn.select('ShadowMask').Or(imgRdn.select('CloudMask'))
  
  CloudShadowMaskAfterBuffer=CloudShadowMask.distance(CloudShadowKernel,False).gte(0).unmask().Not()
  
  img_plus_angles_masked = img_plus_angles.updateMask(CloudShadowMaskAfterBuffer)
  
  ################
  #Calculate BRDF coefficients for a single band
  #and than apply this coefficients for correction
  def apply_MultiRegBRDFccorr(eachBandNamePar):
    #get the water layer aside
    waterLayerMask=LC_MaskCols__.toList(LC_MaskCols__.size()).get(0)
    
    #extract water img
    WaterImg=img_plus_angles_masked.select([eachBandNamePar]).updateMask(ee.Image(waterLayerMask))
    
    
    def apply_OnEachLand(Each_LC_Mask):
        LCMaskedInputImg=img_plus_angles_masked.updateMask(Each_LC_Mask)
        
        #Calculate BRDF coefficients by linear regression 
        # Define the linear regression reducer.
        Reducer = ee.Reducer.linearRegression(
        numX=4,
        numY=1)
  
        #Fit between independent variables angles
        #and dependent variable (ref):
        #VZA: View Zenith angle
        #RAA: Relative Azimuth angle
        #SZA: Solar Zenith angle
        #BRDF = f_vza*(VZA) + f_raa*(RAA) + f_sza*(SZA) + f_constant − ρAVG.
        #ρUncorr=BRDF+ρAVG=f_vza*(VZA) + f_raa*(RAA) + f_sza*(SZA) + f_constant
        #ρCorr=ρUncorr-BRDF
        #,.add(image.select('ρAVG')
        #For example, independent = ee.Image([sin, cos, time, 1])
        brdf_coeffsByBand = LCMaskedInputImg.select(['constant',
            'VZA',
            'SAA',
            'SZA',eachBandNamePar]).reduceRegion(
                reducer=Reducer,
                geometry=studyArea,
                scale=scale,
                maxPixels=10**13,
                tileScale=16)
        
        #get the mean ref value
        ρAVG=LCMaskedInputImg.select([eachBandNamePar]).reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=studyArea,
            scale=scale,
            maxPixels=10**13,
            tileScale=4).get(eachBandNamePar)
        
        #reorignize the brdf coeffs
        #set brdf_coeffs columns' name ['SZA','RAA','VZA' and constant],
        #Corresponding to ['fvol','fgeom','fiso']  
        brdf_coeffsDictByBand=ee.Dictionary({
            'f_constant': ee.Array(brdf_coeffsByBand.get('coefficients')).get([0,0]),
            # f_order: ee.Array(brdf_coeffsByBand.get('coefficients')).get([1,0]),
            'f_vza': ee.Array(brdf_coeffsByBand.get('coefficients')).get([1,0]),
            'f_saa': ee.Array(brdf_coeffsByBand.get('coefficients')).get([2,0]),
            'f_sza': ee.Array(brdf_coeffsByBand.get('coefficients')).get([3,0])})
      
        
        # #brdf=(f_vza*(VZA) + f_raa*(RAA) + f_sza*(SZA) + f_content)− ρAVG.
        # brdf=ee.Image(ee.Number(brdf_coeffsDictByBand.get('f_vza'))).multiply(sensor_zn
        #   ).add(ee.Image(ee.Number(brdf_coeffsDictByBand.get('f_raa'))).multiply(relative_az)
        #   ).add(ee.Image(ee.Number(brdf_coeffsDictByBand.get('f_sza'))).multiply(solar_zn)
        #   ).add(ee.Image(ee.Number(brdf_coeffsDictByBand.get('f_constant')))
        #   ).subtract(ee.Image.constant(ρAVG))
          
        #brdf_nadir=(f_vza*(0) + f_raa*(RAA) + f_sza*(SZA) + f_content)− ρAVG.
        nearNadirAngle=ee.Image.constant(2).multiply(math.pi).divide(180)
        
        simulatedBRDF_nadir=ee.Image(ee.Number(brdf_coeffsDictByBand.get('f_vza'))).multiply(nearNadirAngle
            # ).add(ee.Image(ee.Number(brdf_coeffsDictByBand.get('f_order'))).multiply(LCMaskedInputImg.select('ImgOrder'))
            ).add(ee.Image(ee.Number(brdf_coeffsDictByBand.get('f_saa')))
            .multiply(LCMaskedInputImg.select('SAA'))
            ).add(ee.Image(ee.Number(brdf_coeffsDictByBand.get('f_sza')))
            .multiply(LCMaskedInputImg.select('SZA'))
            ).add(ee.Image(ee.Number(brdf_coeffsDictByBand.get('f_constant')))
            ).subtract(ee.Image.constant(ρAVG))
  
  
        #ρCorr=ρUncorr-brdf_nadir
        brdf_Corrected=LCMaskedInputImg.select([eachBandNamePar]).subtract(simulatedBRDF_nadir)
          
        #make sure the outside of the specific land cover area were masked out
        return brdf_Corrected
    
    #loop LC_MaskCols through map function
    brdf_Corrected_LC_List=LC_MaskCols__.toList(NumRemapClass__,1).map(apply_OnEachLand) #finished LC loop
    
    #add water layer into the corrected img list
    brdf_Corrected_LC_List=ee.List(brdf_Corrected_LC_List).add(WaterImg)
    
    #return the mosaic of BRDF corrected Ref of different LCs
    return ee.ImageCollection(brdf_Corrected_LC_List).mosaic()    
  
  # Calcualate brdf correction coefficients for target bands
  #and apply the correction
  img_BRDFcorr = ee.ImageCollection(bandList.map(apply_MultiRegBRDFccorr)).toBands().rename(bandList)
  
  #select the corrected bands
  # img_BRDFcorr = img_BRDFcorr.select(bandList)
  
  #added the rest bands back
  return img_BRDFcorr
#  return img_BRDFcorr.addBands(ee.Image([ee.Image.constant(1).float().rename('constant'), #constant
#      sensor_zn.rename('VZA'),
#      # ImgOrder.select('ImgOrder'),
#      # sensor_az.rename('VAA'),
#      solar_az.rename('SAA'),
#      solar_zn.rename('SZA')]))

####End brdfCorrection sub-function ################################


###################################################
# MultiRegBRDF_AdditiveCorrection_Apply
def MultiRegBRDF_AdditiveCorrection_Apply(RefTestBandsList,RdnTestBandsList,LCMaskCols_,NumRemapClass_):
    def brdfCorrFun(InputImg):
      
        imageRfl=InputImg.select(RefTestBandsList)
        imageRdn=InputImg.select(RdnTestBandsList)
        
        # VZA = imageRdn.select(['sensorZenith']).multiply(
        # imageRdn.select(['sensorAzimuth']).gt(180).subtract(0.5).divide(0.5))
        
        # FlightDir = imageRdn.select(['sensorAzimuth']).multiply(ee.Image(math.pi/180))
        # .cos().gt(0).subtract(ee.Image.constant(0.5)).divide(0.5).rename('FlightDir')
        # .reduceRegion(ee.Reducer.mode(), ee.Geometry(image.get('footprint')), 2000)
        
        # sensor_dir_value = FlightDir.get('FlightDir')
        # VZA = VZA.multiply(ee.Image.constant(ee.Number(sensor_dir_value).multiply(-1))).rename('VZA')
        
        # imageRdn=imageRdn.addBands(VZA)
      
        #call function MultiRegressionBRDFAdditiveCorr(imgRFL, imgRdn, scale, studyArea, bandList)
        brdfCorrImg = MultiRegressionBRDFAdditiveCorr(
        LCMaskCols_,
        NumRemapClass_,
        imageRfl.select(RefTestBandsList),
        imageRdn, 
        scale, 
        ee.Geometry(InputImg.get('Fulfootprint')),
        RefTestBandsList)
         
        return brdfCorrImg.copyProperties(InputImg)
    
    return brdfCorrFun

##End brdfCorrection###########################


################################################
#get the each land cover mask
def GetLandCoverMask(LCValue):
  return LandCoverRaster.eq(ee.Number(LCValue))

#-------------------------------------------------------------------*/
##########--------End Section 1 Processing Functions----------##########
######################################################################


#-------------------------------------------------------------------*/
# Section 2: Parameters Preparation
#-------------------------------------------------------------------*/
#set initial scale
scale= 5

# #projects/GlobalReservoirs/AVIRIS_NG_WS/AVIRIS_NG_RFL_S2
# #projects/GlobalReservoirs/AVIRIS_NG_WS/AVIRIS_NG_RFL_S1
# AVIRIS_NG_S2_Cols=ee.ImageCollection("projects/GlobalReservoirs/AVIRIS_NG_WS/AVIRIS_NG_RFL_S2")
# .filterMetadata('name','not_contains','ang20190715')
# # .filterMetadata('name','not_contains','ang20190716')
# # .filterMetadata('name','not_contains','ang20190715t174825')
# # .filterMetadata('name','not_contains','ang20190715t180405')
# # .filterMetadata('name','not_contains','ang20190715t181946')
# # .filterMetadata('name','not_contains','ang20190715t192202')
# # .filterMetadata('name','not_contains','ang20190715t183539')
# # .filterMetadata('name','not_contains','ang20190715t185057')


# # .filterMetadata('name','not_contains','ang20190716t182617')
# .filterMetadata('name','not_contains','ang20190716t180912')#ang20190716t182617

# # .filterMetadata('name','not_contains','ang20190716t170505')
# .filterMetadata('name','not_contains','ang20190716t160011')#ang20190716t170505
# # .filterMetadata('name','not_contains','ang20190716t172106')
# .filterMetadata('name','not_contains','ang20190716t161629')#ang20190716t172106

# # .filterMetadata('name','not_contains','ang20190716t173655')
# .filterMetadata('name','not_contains','ang20190716t163352')#ang20190716t173655#ang20190715t192202
# #bad flightline: ang20190716t163352

# # .filterMetadata('name','not_contains','ang20190716t164941')
# # .filterMetadata('name','not_contains','ang20190716t175243')#ang20190716t164941
# #ang20190716t164941/ang20190716t170505

# .filterMetadata('name','not_contains','ang20190716t164941')
# # ang20190716t172106
# .filterMetadata('name','not_contains','ang20190716t172106')
# .filterMetadata('name','not_contains','ang20190716t170505')
# #ang20190716t173655
# .filterMetadata('name','not_contains','ang20190716t173655')
# #ang20190716t182617
# .filterMetadata('name','not_contains','ang20190716t182617')
# # ang20190716t184314
# # .filterMetadata('name','contains','ang20190716t184314')
# .filterMetadata('name','not_contains','ang20190716t184314')


#   'ang20190715t174825',
#   'ang20190715t180405',
#   'ang20190715t181946',
#   'ang20190715t183539',
#   'ang20190715t185057',
#   'ang20190715t190704',
#   'ang20190715t192202'

AVIRIS_NG_S2_Cols=ee.ImageCollection("projects/GlobalReservoirs/AVIRIS_NG_WS/AVIRIS_NG_RFL_S2_Up")\
.filterMetadata('name','not_contains','ang20190716')\
.filterMetadata('name','not_contains','ang20190715t174825')\
.filterMetadata('name','not_contains','ang20190715t180405')\
.filterMetadata('name','not_contains','ang20190715t181946')\
.filterMetadata('name','not_contains','ang20190715t192202')\
.filterMetadata('name','not_contains','ang20190715t183539')\
.filterMetadata('name','not_contains','ang20190715t185057')


print('AVIRIS_NG_S2_Cols',AVIRIS_NG_S2_Cols.size().getInfo())

#Load land cover raster asset
#LandCoverRaster = #ee.Image("projects/GlobalReservoirs/AVIRIS_NG_WS/PAD_S2_LC_Jun-19_Aug-03_2019_noBlueNoDEM").updateMask(1)
LandCoverRaster = ee.Image("projects/GlobalReservoirs/AVIRIS_NG_WS/PAD_S2_LC_Jun-19_Aug-03_2019_DEM").updateMask(1)

#reclass land cover map to make it simple
# from = ee.List([0,1,2,3,4,5])
# to = ee.List([1,2,2,0,3,3])
#fromList = ee.List([0,1,2,3,4,5])
#toList = ee.List([1,2,0,3,0,4])
fromList = ee.List([1,2,3,4,5,6,7,8,9,10,11])
toList = ee.List([0,0,1,2,2,1,1,2,2,2,1])

#define remaped LC class number
NumRemapClass=3

#remap the land cover map
LandCoverRaster = LandCoverRaster.remap(fromList,toList)

#get Mask by clouds and shadows and specific land cover
LandCoverMaskCols=ee.ImageCollection(ee.List.sequence(0, NumRemapClass-1).map(GetLandCoverMask))

#-------------------------------------------------------------------*/
##########--------End Section 2 Parameters Preparation----------##########
######################################################################

#-------------------------------------------------------------------*/
# Section 3: BRDF correction
#-------------------------------------------------------------------*/

#combine the flightline
combAVIRIS_NG_Ref_Cols=AVIRIS_NG_S2_Cols.filterMetadata('system:index','contains','00000')\
.map(combFlightLineS3Fun(AVIRIS_NG_S2_Cols))
#print('combAVIRIS_NG_Ref_Cols',combAVIRIS_NG_Ref_Cols)


#Observational variables
RdnTestBandsList = ['sensorAzimuth','sensorZenith','sunAzimuth','sunZenith','CloudMask','ShadowMask']

TargetBandNames=ee.Image(AVIRIS_NG_S2_Cols.first()).bandNames()
#print('TargetBandNames',TargetBandNames)

noRefBandsList= ['flightHeight','sensorAzimuth','sensorZenith','sunAzimuth','sunZenith',
'slope','aspect','CloudMask','ShadowMask','cloudScore','viewZenith']

TargetBandNamesFiltered=TargetBandNames.removeAll(noRefBandsList)
#print('TargetBandNamesFiltered',TargetBandNamesFiltered)

#bandList = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'] 
# TestBandsList=[450,550,650,860,1580,2180]
# RefTestBandsList=["b16","b36","b56","b97","b241","b361"]

# RefTestBandsList=["b16","b36","b56"]
# print('RefTestBandsList',RefTestBandsList)

#In JavaScript API, it is not easy to automatic running the task one by one
#this is why I convert the original javascript into python api
RefTestBandsList=TargetBandNamesFiltered

# ##############################
# ####BRDF Correction
combAVIRIS_NG_RefBRDFCorrCols=combAVIRIS_NG_Ref_Cols.map(
  MultiRegBRDF_AdditiveCorrection_Apply(RefTestBandsList,RdnTestBandsList,LandCoverMaskCols,NumRemapClass))
#print('brdfCorr',combAVIRIS_NG_RefBRDFCorrCols.size().getInfo())

#-------------------------------------------------------------------*/
##########--------End Section 3 BRDF correction----------##########
######################################################################

#-------------------------------------------------------------------*/
# Section 4: Split BRDF corrected flightline and Export it
#-------------------------------------------------------------------*/

#split brdf-corrected ref of one flightline into the same extent of image chunk of Step2 
AVIRIS_NG_BRDFCorrImgChunkCols=AVIRIS_NG_S2_Cols\
.map(SplitBRDFCorrImgChunkFun(combAVIRIS_NG_RefBRDFCorrCols))
#print('split brdf-corrected ref chunk',AVIRIS_NG_BRDFCorrImgChunkCols.size().getInfo())


#set the output root path
exportPathRoot='projects/GlobalReservoirs/AVIRIS_NG_WS/AVIRIS_NG_RFL_S3_Pre2'

#get output projection
RasterProjection=ee.Image(AVIRIS_NG_S2_Cols.first()).projection()

#get the scale
scale=RasterProjection.nominalScale()
print('scale',scale.getInfo())
scale=5

#get the crs
crs = RasterProjection.crs()


#filter the original image collection to 
#limit the elements to be processed at each time
ChunkFileNameList=AVIRIS_NG_S2_Cols.toList(AVIRIS_NG_S2_Cols.size())\
.map(lambda img: ee.Image(img).get('system:index')).getInfo()
print('ChunkFileNameList',ChunkFileNameList.getInfo())
#ChunkFileNameList=['ang20190715t190704_rfl_00000']

#loop to export the results
for ChunkFileName in ChunkFileNameList:
    #filter to get the target chuck image
    ChunkRefRasterImg=ee.Image(AVIRIS_NG_BRDFCorrImgChunkCols.\
            filterMetadata('system:index','contains',ChunkFileName).first())\
            .reproject(ee.Projection('EPSG:32612').atScale(5)).float()
    
    #set a output image filename
    ChunkFileName=ChunkFileName+'_py'
    print('ChunkFileName',ChunkFileName)
    
    #get the footprint of image
    Footprint=ChunkRefRasterImg.get('footprint')
    
    #convert the footprint into polygon
    exportFootprintPolygon=ee.Algorithms.GeometryConstructors.Polygon(
    ee.Geometry(Footprint).coordinates())
    
    #set asset id
    assetPath =exportPathRoot + '/' + ChunkFileName
    
    taskIA2Assets_config=dict(assetId=assetPath,
        pyramidingPolicy= {'.default': 'mean'},
        region=exportFootprintPolygon,
        scale=scale,
        crs=crs,
        maxPixels=10**13)
            
    #call export to Asset API
    taskIA2assets = ee.batch.Export.image.toAsset(
        image=ChunkRefRasterImg,
        description=ChunkFileName,**taskIA2Assets_config)
        
    #begin task "export image"
    taskIA2assets.start()


#-------------------------------------------------------------------*/
##########--------End Section 4 Export flightline----------##########
######################################################################


