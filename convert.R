library(av)
setwd("videos")

f <- list.files(pattern=".MP4",recursive = T)
f_avi <- list.files(pattern=".avi",recursive = T)

f <- f[!tools::file_path_sans_ext(basename(f))%in%tools::file_path_sans_ext(basename(f_avi))]

img <- "images"
dir.create(img)
for(i in f) 
{
  try_ <- try({
  if(dir.exists(img)) unlink(img,recursive = TRUE)
  av_video_images(i,destdir=img)
  imgs <- list.files(img,full.names=TRUE)
  bn <- gsub("\\.MP4|\\.mp4","",basename(i))
  cp <- dirname(i)
  av_encode_video(imgs,output = paste0(bn,".avi"),codec = "rawvideo",vfilter = "transpose=1")
  }
  ,silent = T
  )
if(inherits(try_,"try-error")) warning(paste0("video ",i," did not render"))
}


av_video_convert()
av::av_video_convert(imgs[1:20],"piv.gif")