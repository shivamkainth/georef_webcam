#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#############################################################################################################
Created on Wed Oct 30 2019
Last edited on Fri Feb 11 2022

Author: Sebastian Buchelt

/******************************************************************************
 *                                                                            *
 *   This program is public software; It is distributed under the the terms   *
 *   of the Creative Commons Attribution-NonCommercial-ShareAlike 4.0         *
 *   International Public License as published by the Creative Commons        *
 *   Corporation; either version 4 of the License, or (at your option) any    *
 *   later version.                                                           *
 *                                                                            *
 ******************************************************************************/
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Name:       procedures/project2map.py
%   Purpose:    Project Images and other data into map coordinates using the 
%               output of the georef_webcam projection procedure
%   Comment:    This file contains the functions, which generate maps from  
%               images and other processed data by using the coordinate rasters 
%               and the mask generated by georef_webcam.
%
%   Overview:   project_tif: projects tif in image projection to map 
%                       coordinates
%               project_image: converts image data in other formats to tif and 
%                       then calls the project_tif function
%               project_to_image_plane: projects geotifs, which are in the same 
%                       CRS as the DEM that was used for the projection, into 
%                       the view/2-D image plane of the camera
%               project_geometry: projects shp files from image plane to map or 
%                       opposite direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#############################################################################################################
"""

###############################################################################
####### function to project data converted to tif into map coordinates ########
# input:
#       - coord_dir: directory to mask and coordinate raster layers
#       - px_size: spatial resolution of projected map
#       - pj_file: file which contains CRS information
#       - image_folder: directory, where images that should be projected are stored
#       - out_dir: output directory where projected maps are stored
#       - fill_nodata: boolean; if True, gaps in projected map are filled by interpolation
#       - file: boolean; set True, if only a single file is projected
###############################################################################
def project_tif(coord_dir, px_size, pj_file, image_folder, out_dir, fill_nodata=False, file=False):
    # import required libraries
    from osgeo.gdalconst import GA_ReadOnly
    import os, gdal, sys, math, subprocess, platform
    import numpy as np
    
    ##### get list of the images from image folder directory 
    if file:
        img_file_list = [image_folder]
    else:
        img_file_list = [os.path.join(image_folder, f) for f in os.listdir(image_folder) if (f.endswith('.tif'))]
        img_file_list.sort()
    
    ##### create directory, where output should be stored
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    
    ##### open & read coordinate raster and mask layers
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()
    input_list = ['mask.tif', 'north_raster.tif', 'east_raster.tif']
    for i in input_list:
    	datafile = os.path.join(coord_dir, i)
    	inDs = gdal.Open(datafile, GA_ReadOnly)
    	if inDs is None:
             print ('could not open ' + datafile)
             sys.exit(1)
    	canal = inDs.GetRasterBand(1)
    	if (i == 'mask.tif'):
             mask = canal.ReadAsArray().astype(np.float)
             cols = inDs.RasterXSize
             rows = inDs.RasterYSize
    	elif (i == 'north_raster.tif'):
    		north = canal.ReadAsArray().astype(np.float)
    	elif (i == 'east_raster.tif'):
    		east = canal.ReadAsArray().astype(np.float)
    
    ####################### derive extent of map array  #######################
    ##### mask georeferencing layers
    north[mask==0] = np.nan
    east[mask==0] = np.nan
    
    ##### get maximum spatial extent uf remaining image area in world coordinate system
    min_max_east = (np.nanmin(east), np.nanmax(east))
    min_max_north = (np.nanmin(north), np.nanmax(north))
    
    ##### correct extent values, so that whole extent is covered & extent is divideable by pixel_size
    half_px = px_size/2
    min_max_east = (round(min_max_east[0]/px_size)*px_size-half_px, round(min_max_east[1]/px_size)*px_size+half_px)
    min_max_north = (round(min_max_north[0]/px_size)*px_size-half_px, round(min_max_north[1]/px_size)*px_size+half_px)
    
    ##### get number of rows & cols of projected map
    row_number = int((min_max_north[1]-min_max_north[0])/px_size)
    col_number = int((min_max_east[1]-min_max_east[0])/px_size)
    
    
    ################ calculate for each image pixel the #######################
    ############## correlating position in the new array ######################
    pos_new = np.ones((rows,cols,2))*-1
    for col in range(cols):
        for row in range(rows):
            if not (math.isnan(east[row,col])) and not (math.isnan(north[row,col])):
                pos_new[row,col] = [int((east[row,col]-min_max_east[0])/px_size), int((north[row,col]-min_max_north[0])/px_size)]
    
    ################# read tif data which should be projected #################
    for img_file in img_file_list:
        inDs = gdal.Open(img_file, GA_ReadOnly)
        no_of_bands = inDs.RasterCount
        
        ##### create tif-file for projected result
        img_name = os.path.basename(img_file)               # get image name
        tif_name = img_name.split('.')[0]+'_map.tif'         # create output map name from that
        dst_ds = gdal.GetDriverByName('GTiff').Create(os.path.join(out_dir,tif_name), col_number, row_number, no_of_bands, gdal.GDT_Float64)     # create the single band raster tif-file
        
    #################### run projection for each band #########################
        for i in range(no_of_bands):
            ##### open classified image
            canal = inDs.GetRasterBand(i+1)
            noDataVal = canal.GetNoDataValue()
            image = canal.ReadAsArray()
            
            ##### create output arrays
            count_array = np.zeros((row_number, col_number))        # counts the number of image pixels assigned to each map pixel
            snow_val = np.full((row_number, col_number), -9999)     # output array for map 
            
            ##### assign image values to map pixels
            for col in range(cols):
                for row in range(rows):
                    if not (image[row, col]==noDataVal):                              # exclude no-Data pixels
                        if not -1 in pos_new[row,col]:
                            east_pos = int(pos_new[row,col,0])                  # get easting position in map array
                            north_pos = int(pos_new[row,col,1])+1               # get northing position in map array
                            # value assignment: (value*no_of_already_assigned_values + new_value_to_be_added)/new_number_of_assigned_values
                            snow_val[row_number-north_pos, east_pos] = (snow_val[row_number-north_pos, east_pos]*count_array[row_number-north_pos,east_pos]+image[row, col])/(count_array[row_number-north_pos,east_pos]+1)
                            count_array[row_number-north_pos,east_pos] +=1      # increase counter for this pixel by 1
            snow_val[count_array==0]=-9999                                      # assign noData-value, if no value was assigned before
            
            ##### write output to tif-file   
            dst_ds.GetRasterBand(i+1).WriteArray(snow_val)        # write result array into tif-file
            dst_ds.GetRasterBand(i+1).SetNoDataValue(-9999)       # define no-Data-value to -9999
        
    ########## define geotransformation & projection and save result ##########
        geotransform = (min_max_east[0], px_size, 0, min_max_north[1], 0, -px_size)     # geotransformation from spatial extent and pixel size
        dst_ds.SetGeoTransform(geotransform)                                            # set geotransfornation to tif-file
        
        # get & define projection from DEM file and save result
        dst_ds.SetProjection(pj_file)      # set projection to tif-file
        dst_ds.FlushCache()                     # write to disk
        dst_ds = None                           # save, close     
        canal = None
        inDs = None
    
    ######## optionnal: fill noData-gaps with small scale interpolation #######
    ############## using the gdal function gdal_fillnodata.py #################
    # value of fill_nodata decides to what pixel range the interpolation is applied 
        if(fill_nodata):
            for i in range(no_of_bands):
                if(platform.system()=='Windows'):
                    subprocess.Popen("gdal_fillnodata.py -md "+str(fill_nodata)+" "+os.path.join(out_dir,tif_name)+" -b " +str(i+1), shell = True,stdout=subprocess.PIPE)
                else:
                    subprocess.Popen(['gdal_fillnodata.py', '-md', str(fill_nodata), os.path.join(out_dir,tif_name), '-b', str(i+1)],stdout=subprocess.PIPE)
            
        print (img_name + ' is processed')      # print to console, that processing of file x is finished
   
################################### end #######################################
###############################################################################
        
        

###############################################################################
####### function to project data converted to tif into map coordinates ########
# input:
#       - coord_dir: directory to mask and coordinate raster layers
#       - px_size: spatial resolution of projected map
#       - pj_file: file which contains CRS information
#       - image_folder: directory, where images that should be projected are stored
#       - file_ending: file extension to look for in image folder
#       - image_file: set to filename, if only a single file is projected
#       - out_dir: output directory where projected maps are stored
#       - fill_nodata: boolean; if True, gaps in projected map are filled by interpolation
#       - data: boolean; set True, if you want to read png files via external function
###############################################################################
def project_image(coord_dir, px_size, pj_file, image_folder=None, file_ending=None, image_file = None, out_dir=None, fill_nodata=False, data=False):
    # import required libraries and submodules of georef_webcam
    import os, gdal
    from PIL import Image
    import numpy as np
    import modules.aux_results as aux_res
    
    ##### get list of the coregistered images from image folder directory, 
    # if not specified, read the master image from coord_dir
    if not (image_file == None):
        img_file_list = [image_file]
        tif_dir = os.path.join(os.path.dirname(image_file), 'image_tif')
    elif (image_folder == None):
        img_names = list(['image'+file_ending])
        img_file_list = [os.path.join(coord_dir, f)  for f in img_names]
        tif_dir = os.path.join(coord_dir, 'image_tif')
    else: 
        img_names = os.listdir(image_folder)
        img_file_list = [os.path.join(image_folder, f) for f in img_names if (f.endswith(file_ending))]
        img_file_list.sort()
        tif_dir = os.path.join(os.path.dirname(image_folder), 'image_tif')
        
    # define & create directory, where output and intermediate results should be stored
    if (out_dir == None):
        if (image_folder == None):
            out_dir = coord_dir
        else:
            out_dir = os.path.join(image_folder, 'projected')
    if not os.path.exists(tif_dir):
        os.makedirs(tif_dir)
        
    ############### convert image data to tif files ###########################
    for img_file in img_file_list:
        tif_file = os.path.join(tif_dir,os.path.splitext(os.path.basename(img_file))[0]+'.tif')
        ##### open coregistered image
        if os.path.splitext(img_file)[1]=='.png' and data:
            image = aux_res.read_png(img_file)
        else:
        	im1 =  Image.open(img_file)
        	image = np.array(im1)
        row_number, col_number, no_of_bands = image.shape
        dst_ds = gdal.GetDriverByName('GTiff').Create(tif_file, col_number, row_number, no_of_bands, gdal.GDT_Float64)   # create the 3-band raster tif-file
        for i in range(no_of_bands):
            dst_ds.GetRasterBand(i+1).WriteArray(image[:,:,i])      # write image band to the raster
            
        dst_ds.FlushCache()                     # write to disk
        dst_ds = None                           # save, close     
        
    ############### execute projection and delete intermediate tifs ###########
        project_tif(coord_dir, px_size, pj_file, tif_file, out_dir, fill_nodata=fill_nodata, file=True)
        if os.path.isfile(tif_file):
            os.remove(tif_file)
    if os.path.exists(tif_dir):
        if len(os.listdir(tif_dir))==0:
            os.rmdir(tif_dir)
################################### end #######################################
###############################################################################
            
            
            
###############################################################################
############### function to project geotifs into image plane ##################
# input:
#       - coord_dir: directory to mask and coordinate raster layers
#       - image_folder: directory, where images that should be projected are stored
#       - out_dir: output directory where projected maps are stored
#       - file: boolean; set True, if only a single file is projected
###############################################################################
def project_to_image_plane(coord_dir, image_folder, out_dir, file=False):
    # import required libraries
    from osgeo.gdalconst import GA_ReadOnly
    import os, gdal, sys, struct
    import numpy as np
    
    ##### get list of the images from image folder directory 
    if file:
        img_file_list = [image_folder]
    else:
        img_file_list = [os.path.join(image_folder, f) for f in os.listdir(image_folder) if (f.endswith('.tif'))]
        img_file_list.sort()
    
    ##### create directory, where output should be stored
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    
    ##### open & read coordinate raster and mask layers
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()
    input_list = ['mask.tif', 'north_raster.tif', 'east_raster.tif']
    for i in input_list:
    	datafile = os.path.join(coord_dir, i)
    	inDs = gdal.Open(datafile, GA_ReadOnly)
    	if inDs is None:
             print ('could not open ' + datafile)
             sys.exit(1)
    	canal = inDs.GetRasterBand(1)
    	if (i == 'mask.tif'):
             mask = canal.ReadAsArray().astype(np.float)
             cols = inDs.RasterXSize
             rows = inDs.RasterYSize
    	elif (i == 'north_raster.tif'):
    		north = canal.ReadAsArray().astype(np.float)
    	elif (i == 'east_raster.tif'):
    		east = canal.ReadAsArray().astype(np.float)
    
    
    ############### read tif data which should be projected ###################
    for img_file in img_file_list:
        inDs = gdal.Open(img_file, GA_ReadOnly)
        raster_col = inDs.RasterXSize
        raster_row = inDs.RasterYSize

        no_of_bands = inDs.RasterCount
        gt=inDs.GetGeoTransform()
        
        ##### create tif-file for projected result
        img_name = os.path.basename(img_file)               # get image name
        tif_name = img_name.split('.')[0]+'_projected.tif'         # create output map name from that
        resDs = driver.Create(os.path.join(out_dir, tif_name), cols, rows, no_of_bands, gdal.GDT_Float64)
        
    ############### run projection ############################################
        sat_proj = np.zeros((rows,cols,no_of_bands))
        for row in range(0,rows):
            for col in range(0,cols):
                if mask[row,col] and not any([np.isnan(mask[row,col]),np.isnan(east[row,col]), np.isnan(north[row,col])]):
                    ##### identify corresponding map pixel for specific image pixel
                    mx = east[row,col]
                    my = north[row,col]
                    px = int(round((mx - gt[0]) / gt[1])) #x pixel position in raster
                    py = int(round((my - gt[3]) / gt[5])) #y pixel position in raster
                    
                    ##### extract data from each band & store it in array
                    if(px>=0 and px<raster_col and py>=0 and py<raster_row):
                        for b in range(no_of_bands):                
                            ##### open classified image
                            canal = inDs.GetRasterBand(b+1)                
                            data = canal.ReadRaster(px,py,1,1,buf_type=gdal.GDT_UInt16) #Assumes 16 bit int aka 'short'
                            sat_proj[row,col,b] = struct.unpack('h' , data)[0]
                    
                    ##### set undefined pixels and pixels outside extent to noDAtaValue
                    else:
                        sat_proj[row,col] = -9999
                else:
                    sat_proj[row,col] = -9999
                    
        ##### save results as tif file
        for b in range(no_of_bands):                
            resBand = resDs.GetRasterBand(b+1)
            resBand.WriteArray(sat_proj[:,:,b])
            resBand.SetNoDataValue(-9999)       # define no-Data-value to -9999
            resBand.FlushCache()
            resBand = None
        resDs = None
################################### end #######################################
###############################################################################
        


###############################################################################
#### function to project shp files into image plane or to map coordinates #####
# input:
#       - coord_dir: directory to mask and coordinate raster layers
#       - pj_file: file which contains CRS information
#       - image_folder: directory, where images that should be projected are stored
#       - out_dir: output directory where projected maps are stored
#       - file: boolean; set True, if only a single file is projected
#       - to_map: boolean; if True, shp is projected from image view to map
#                          if False, shp is projected from map to image view
###############################################################################
def project_geometry(coord_dir, pj_file, image_folder, out_dir, file=False, to_image_view=True):
    # import required libraries
    import osgeo.ogr as ogr
    import osgeo.osr as osr
    import os, gdal, sys
    import numpy as np
    from osgeo.gdalconst import GA_ReadOnly

    ##### get list of the images from image folder directory 
    if file:
        img_file_list = [image_folder]
    else:
        img_file_list = [os.path.join(image_folder, f) for f in os.listdir(image_folder) if (f.endswith('.shp'))]
        img_file_list.sort()
    
    ##### create directory, where output should be stored
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    ##### open & read coordinate raster and mask layers
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()
    input_list = ['north_raster.tif', 'east_raster.tif']
    for i in input_list:
    	datafile = os.path.join(coord_dir, i)
    	inDs = gdal.Open(datafile, GA_ReadOnly)
    	if inDs is None:
             print ('could not open ' + datafile)
             sys.exit(1)
    	canal = inDs.GetRasterBand(1)
    	if (i == 'north_raster.tif'):
    		north = canal.ReadAsArray().astype(np.float)
    	elif (i == 'east_raster.tif'):
    		east = canal.ReadAsArray().astype(np.float)
                
    ##### get prj data for projection procedure
    if to_image_view:    
        map_pj = osr.SpatialReference()
        map_pj.ImportFromEPSG(4326)
        append = '_to_img.shp'
    else:        
        map_pj = osr.SpatialReference()
        map_pj.ImportFromWkt(pj_file)
        append = '_to_map.shp'
    
    
    ############### read shp data which should be projected ###################
    driver = ogr.GetDriverByName("ESRI Shapefile")
    for img_file in img_file_list:
        dataSource = driver.Open(img_file, 0)
        layer = dataSource.GetLayer()    
        layerDefinition = layer.GetLayerDefn()
        
        ##### create shp-file for projected result
        img_name = os.path.basename(img_file)               # get image name 
        shp_name = img_name.split('.')[0]+append         # create output map name from that
        data_source = driver.CreateDataSource(os.path.join(out_dir, shp_name))
        layer2 = data_source.CreateLayer("projected",map_pj, layer.GetGeomType())
        
        ##### Add the fields to new shp file
        for i in range(layerDefinition.GetFieldCount()):
            # get info from existing shp
            fieldName =  layerDefinition.GetFieldDefn(i).GetName()
            fieldTypeCode = layerDefinition.GetFieldDefn(i).GetType()
            fieldWidth = layerDefinition.GetFieldDefn(i).GetWidth()
            #copy info to new shp
            new_field = ogr.FieldDefn(fieldName, fieldTypeCode)
            new_field.SetWidth(fieldWidth)
            layer2.CreateField(new_field)
            
    ############### Project features ##########################################
        for feature in layer:
            ##### get single features
            feature_out = ogr.Feature(layer2.GetLayerDefn())
            geom = feature.GetGeometryRef()
            # if Polygon extract ring
            if layer.GetGeomType()==3:
                geom = geom.GetGeometryRef(0)
                ring = ogr.Geometry(ogr.wkbLinearRing)
            geom_out = ogr.Geometry(layer2.GetGeomType())
            
            ##### get position and project each point in feature
            for i in range(geom.GetPointCount()):
                data = geom.GetPoint(i)
                # or to image plane
                if to_image_view:    
                    (easting, northing, z) = data
                    if not any([np.isnan(easting),np.isnan(northing)]):
                        distance = (north-northing)**2+(east-easting)**2
                        val_1 = np.mean(np.where(distance == np.nanmin(distance))[1].tolist())
                        val_2 = -np.mean(np.where(distance == np.nanmin(distance))[0].tolist())
                    else:
                        val_1 = np.nan
                        val_2 = np.nan
                # either to map coordinates
                else:
                    (col_pos, row_pos,z) = data
                    val_1 = east[int(-(round(row_pos))),int(round(col_pos))]
                    val_2 = north[int(-(round(row_pos))),int(round(col_pos))]
                if not layer.GetGeomType()==3:
                    geom_out.AddPoint(val_1, val_2)
                # if polygon, add points to ring and add it as a Polygon to feature object
                else:
                    ring.AddPoint(val_1, val_2)
            if layer.GetGeomType()==3:
                geom_out.AddGeometry(ring)
                
            ##### store new position and copy variable values to new shp file
            feature_out.SetGeometry(geom_out)
            for i in range(layerDefinition.GetFieldCount()):
                feature_out.SetField(layerDefinition.GetFieldDefn(i).GetName(),feature.GetField(layerDefinition.GetFieldDefn(i).GetName()))
            # Create the feature in the layer (shapefile)
            layer2.CreateFeature(feature_out)
            # Dereference the feature
            feature = None
            
        # Save and close the data source
        data_source = None
################################### end #######################################
###############################################################################