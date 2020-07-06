library(lidR)
library(sf)
library(sp)
library(dplyr)
# This script converts all spatial data to CRS 3946

# This function takes the tile id (tile_id) as input and returns its extent as output
extent_from_tile_id=function(tile_id){
  xmin=as.numeric(stringr::str_sub(num,1,4))*1000
  ymin=as.numeric(stringr::str_sub(num,6,9))*1000
  res_extent=c(xmin,xmin+1000,ymin,ymin+1000) %>% extent()
  return(res_extent)
}

# Read shapefile regarding buildings
bati_og=st_read("data/cad_cadastre.cadbatiment/cad_cadastre.cadbatiment.shp") %>% 
  st_set_crs(4019) %>% 
  st_transform(crs=3946) %>% 
  mutate(bbox=purrr::map(geometry,st_bbox))%>% 
  mutate(xmin=purrr::map(bbox,"xmin"),
         xmax=purrr::map(bbox,"xmax"),
         ymin=purrr::map(bbox,"ymin"),
         ymax=purrr::map(bbox,"ymax")) %>% 
  dplyr::select(-bbox)

# Read DTM
DTM=raster("data/MNT2018_Altitude_2m/MNT2018_Altitude_2m.tif")


# This function calculates the vegetation height raster 
# It takes tile id (tile_id) as input
# Its effect corresponds to the writing of file:
# "data/vegetation_height_raster/corr_vegetation_height_{tile_id}.tif"
calc_vegetation_height=function(tile_id){
  url=paste0("https://download.data.grandlyon.com/files/grandlyon/imagerie/mnt2018/lidar/laz/",
             tile_id,".laz")
  if(!RCurl::url.exists(url)){return()}
  ##### Files' paths
  lazfile=paste0("data/lidar/",tile_id,".laz")
  rasterfile=paste0("data/vegetation_height_raster/corr_vegetation_height_",tile_id,".tif")
    if(!file.exists(lazfile)){
      download.file(url, lazfile)
    }
    las=lidR::readLAS(paste0("data/lidar/",tile_id,".laz"))
    myextent=extent(las)
    lidR::projection(las)=sp::CRS("+init=epsg:3946")
    DTM_crop=crop(DTM,myextent)
    lasveg=lidR::lasfilter(las,Classification %in% 5)
    splas=st_as_sf(lasveg@data[,c("X","Y")], coords=c("X","Y")) %>% 
      st_set_crs(3946) 
    
    splas_bbox=st_bbox(splas)
    bati=bati_og %>% 
      filter(xmin>=splas_bbox$xmin,
             xmax<=splas_bbox$xmax,
             ymin>=splas_bbox$ymin,
             ymax<=splas_bbox$ymax) %>% 
      st_buffer(dist=2)
        ####################################
        ubati=st_combine(bati)
        in_bati=st_intersects(splas,ubati,sparse=FALSE)[,1]
        lasveg@data=cbind(lasveg@data,in_bati)
        ind_in_bati=which(lasveg@data$in_bati)
        if(length(ind_in_bati)>0){
          zcorr=extract(DTM_crop,splas[ind_in_bati,])
          lasveg@data$Z[ind_in_bati]=zcorr
        }
        ############ vegetation_height raster
        gridveg <- grid_vegetation_height(lasveg, res = 0.5, p2r())
        gridveg=extend(gridveg,extent_from_tile_id(tile_id))
        gridDTM=resample(DTM_crop, gridveg)
        gridveg[is.na(gridveg[])] =gridDTM[is.na(gridveg[])] 
        gridveg=gridveg-gridDTM
        ############ substract DSM
        writeRaster(gridveg,
                    rasterfile,
                    overwrite=T)
      ############
      file.remove(lazfile)
}

# Apply function to all relevant tiles:
# Total range 1830:1868 5150:5197
library(RCurl)
listgrids=expand.grid(1831:1867,5151:5196) %>% 
  dplyr::mutate(tile_id=paste0(Var1,"_",Var2)) %>% 
  dplyr::mutate(url=paste0("https://download.data.grandlyon.com/files/grandlyon/imagerie/mnt2018/lidar/laz/",tile_id,".laz")) %>% 
  dplyr::mutate(url_exists=purrr::map(url,url.exists)) %>% 
  dplyr::filter(url_exists==TRUE)

purrr::map(listgrids$tile_id,calc_vegetation_height)

# decrease resolution of rasters by a factor of 4 (=> resolution: 2m)
all_rasters=list.files("data/vegetation_height_raster/")
agg_and_write=function(raster_name,tile_id=NA){
  if(!is.na(tile_id)){raster_name=paste0("corr_vegetation_height_",tile_id,".tif")}
  old_path=paste0("data/vegetation_height_50cm/",raster_name)
  new_path=paste0("data/vegetation_height_2m/",raster_name)
  big_raster=raster(old_path)
  small_raster=raster::aggregate(big_raster, fact=4, fun=max)
  writeRaster(small_raster,new_path, overwrite=TRUE)
}
purrr::map(all_rasters, agg_and_write)

## Merge vegetation_height rasters into larger rasters
## (aggregate based on values x, y)
listgrids =listgrids %>%
  mutate(x=floor((Var1-min(Var1))/5)+1,
         y=floor((Var2-min(Var2))/5)+1) %>% 
  mutate(mid_raster=paste0("large_raster_",x,"_",y,".tif"))
for (i in 1:8){
  for (j in 1:10){
    print(paste0("i=",i,";j=",j))
    raster_out=paste0("data/merged_rasters/large_raster_",i,"_",j,".tif")
    if(!file.exists(raster_out)){
      tile_ids=listgrids %>%
      dplyr::filter(x==i,y==j) %>% 
      select(tile_id) %>% 
      pull() 
    if(length(tile_ids>0)){
    list_rasters=paste0("data/vegetation_height_raster_2m/corr_vegetation_height_",
                        tile_ids,".tif") %>% 
      as.list() %>% 
      purrr::map(raster)
    large_raster=do.call(mosaic,c(list_rasters, fun=mean, tolerance=0.5))
    writeRaster(large_raster,
                raster_out, 
                overwrite=TRUE)
    }
    }
  }
}
## Then merge these together to create  1 global vegetation height raster (2m res.)
global_vegetation_height=do.call(mosaic,c(list.files("data/merged_rasters"),fun=mean,tolerance=0.5))
writeRaster(global_vegetation_height,
            "data/global_vegetation_height_2m.tif",
            overwrite=TRUE)
