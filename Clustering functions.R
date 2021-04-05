change_resolution <- function(image_df, x_size)
{
  ## change_resolution(image_df, x_size) subsamples an image to produce
  ## a lower resolution image. Any non-coordinate columns in the data
  ## frame are summarized with their most common value in the larger
  ## grid cell.
  ##
  ## Input:
  ## - image_df: A data frame in wide format. The x-coordinate column MUST
  ##             be named 'x' and the y-coordinate column MUST be named 'y'.
  ##             Further columns have no naming restrictions.
  ## - x_size:   The number of cells in the x-direction. The number of cells
  ##             in the vertical direction will be computed to maintain the 
  ##             perspective. There is no guarantee that the exact number
  ##             of cells in the x-direction is x_size
  ##
  ## Output:
  ## - A data frame with the same column names as image_df, but with fewer 
  ##   entries that corresponds to the reduced resolution image.
  ##
  ## Example:
  ##   library(imager)
  ##   library(dplyr)
  ##   fpath <- system.file('extdata/Leonardo_Birds.jpg',package='imager') 
  ##   im <- load.image(fpath)
  ##   im_dat<- as.data.frame(im,wide = "c") %>% rename(R = c.1, G = c.2, B = c.3) %>%
  ##            select(x,y,R,G,B)
  ##   agg_image <- change_resolution(im_dat, 50)
  
  if(!require(sp)) {
    stop("The sp packages must be installed. Run install.packages(\"sp\") and then try again.")
  }
  if(!require(dplyr)) {
    stop("The dplyr packages must be installed. Run install.packages(\"dplyr\") and then try again.")
  }
  im_dat<-image_df
  sp_dat <- im_dat 
  gridded(sp_dat) = ~x+y
  
  persp = (gridparameters(sp_dat)$cells.dim[2]/gridparameters(sp_dat)$cells.dim[1])
  y_size = floor(x_size*persp)
  orig_x_size = gridparameters(sp_dat)$cells.dim[1]
  orig_y_size = gridparameters(sp_dat)$cells.dim[2]
  
  x_res = ceiling(orig_x_size/x_size)
  y_res = ceiling(orig_y_size/y_size)
  
  gt = GridTopology(c(0.5,0.5), c(x_res, y_res),
                    c(floor(orig_x_size/x_res), floor(orig_y_size/y_res)))
  SG = SpatialGrid(gt)
  agg = aggregate(sp_dat, SG, function(x) names(which.max(table(x)))[1] )
  agg@grid@cellsize <- c(1,1)
  df <- agg %>% as.data.frame %>% rename(x = s1, y = s2)  %>% select(colnames(image_df))
  
  return(df)
  
}



process_image <- function(image_file_name, k_list){
  
  ## process_image is a function to load image and do kmeans clustering
  ## This will return a large list which will show the results of 
  ## clusterings for each of k in k_list.
  ##
  ## Input: 
  ##
  ## - image_file_name: A string ofthe name of the image file that 
  ##                    we wish to analyze by clustering.
  ## 
  ## - k_list:          A list of k values that we can use to do
  ##                    clustering. Different k-value will promote
  ##                    different clustering results as the number of
  ##                    clusters is varied. 
  ## Output:
  ##
  ## - A large list with two data.frame. One for RGB informations of
  ##   the image selected. One for clustering results. 
  ##
  ## Example: 
  ##    library(imager)
  ##    library(dplyr)
  ##    library(tidymodels)
  ##    im <- 'STA314.jpg'
  ##    k_list <- c(2:8)
  ##    process_image(im,k_list)
  
  
  image <- load.image(image_file_name)
  
  tidy_dat <- as.data.frame(image, wide='c')%>%
    rename(R=c.1, G=c.2, B=c.3)
  
  dat <- select(tidy_dat, c(-x,-y))
  kclusts <- 
    tibble(k_list)%>%
    mutate(
      kclust=map(k_list, ~kmeans(x=dat,centers= .x, nstart = 4)),
      glanced=map(kclust,glance),
      tdata=map(kclust,tidy)
    )
  
  clusterings <-
    kclusts %>%
    unnest(cols=c(glanced))
  
  return(list(tidy_dat, clusterings))
}


scree_plot<-function(cluster_info){
  
  ## scree_plot is a function to produce a scree plot showing the 
  ## relationship between k and its tot.withinss ratio, in order to
  ## help determine optimal k. 
  ##
  ## Input: 
  ##
  ## - cluster_info: A large list generated from process_image
  ##                 function that contains clustering results. 
  ##
  ## Output: 
  ## - A scree plot with two numercial variables: k values and their
  ##   tot.withinss ratio. 
  ## 
  ## Example: 
  ##
  ##    im <- 'STA314.jpg'
  ##    k_list <- c(2:8)
  ##    cluster_info <- process_image(im,k_list)
  ##    scree_plot(cluster_info)
  
  info<-cluster_info[[2]]
  
  nclust=length(info$k_list)
  
  ratio=rep(NA,nclust-1)
  
  for(kk in 2:nclust){
    
    ratio[kk-1]=info$tot.withinss[kk]/info$tot.withinss[kk-1]
  }
  plot_data<-data.frame(k=info$k_list[2:nclust],ratio)
  
  ggplot(plot_data,aes(x=k,y=ratio))+geom_line()
}


color_strips<-function(cluster_info,k){
  
  ## color_strips is a function to display k colors with their hex
  ## strings, according to k clusters. 
  ##
  ## Input: 
  ##
  ## - cluster_info: A large list generated from process_image
  ##                 function that contains clustering results. 
  ## 
  ## - k:            A number that is chosen to be the optimal number
  ##                 of clusters. 
  ## 
  ## Output: 
  ## 
  ## - A strip that displays sorted colors with their hex strings. The
  ##   number of colors equals to input k. 
  ##
  ## Example:
  ##
  ##    im <- 'STA314.jpg'
  ##    k_list <- c(2:8)
  ##    cluster_info <- process_image(im,k_list)
  ##    k = 7
  ##    color_strips(cluster_info, k)
  
  kclust<-kmeans(select(cluster_info[[1]],-x,-y),centers = k,nstart = 20 )
  centres <- tidy(kclust)
  
  centres<-centres%>%mutate(col=rgb(R,G,B))
  
  show_col(centres$col)
  
  return(centres)
}


make_pattern<-function(cluster_info,k,x_size,black_white=FALSE,background_colour=NULL){
  
  ## make_pattern is a function that produces a cross-stitch pattern
  ## of the selected image. The number of colors depends on k
  ## value. 
  ##
  ## Input: 
  ##
  ## - cluster_info: A large list generated from process_image
  ##                 function that contains clustering results. 
  ##
  ## - k:            A number that is chosen to be the optimal number
  ##                 of clusters. 
  ##
  ## - x_size :      The (approximate) total number of possible stitches 
  ##                 in the horizontal direction
  ##
  ## - black_white:  Print the pattern in 
  ##                 black and white (TRUE) or colour (FALSE,default)
  ##                             
  ## - background_colour: The colour of the background, which should
  ##                      not be stitched in the pattern.
  ##
  ## Output: 
  ## - It produces a cross_stitch pattern that can be followed, complete
  ##   with a legend that has thread colour, and a guide grid. 
  ##
  ## Example: 
  ## 
  ##    im <- 'STA314.jpg'
  ##    k_list <- c(2:8)
  ##    cluster_info <- process_image(im,k_list)
  ##    k = 7
  ##    x_size=50
  ##    make_pattern(cluster_info,k,x_size,black_white=FALSE, background_colour=NULL)
  
  kclust <- kmeans(select(cluster_info[[1]],-x,-y),centers = k,nstart = 20 )
  
  tidy_dat<-augment(kclust,cluster_info[[1]])%>%rename(cluster= .cluster)
  
  centres <- tidy(kclust)
  
  centres<-centres%>%mutate(col=rgb(R,G,B))
  
  
  centres1 <- map(centres$col,~dmc(.x,visualize = FALSE,method='euclidean'))%>%tibble%>%unnest(cols = c(.))
  
  
  im_dat <- centres1[array(tidy_dat$cluster),]%>%select(dmc,name,hex)%>%cbind(tidy_dat, .)
  
  
  reso=change_resolution(im_dat,x_size)
  
  
  ggplot(reso,aes(x,y,color=I(hex),shape=cluster))+geom_point()+scale_y_reverse()
  
}
