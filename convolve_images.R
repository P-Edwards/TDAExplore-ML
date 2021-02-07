library(optparse) 
library(utils)
library(colorspace) 
library(ggplot2)
library(gridExtra)
library(ggplotify)

MANUAL_SCALING_FOR_COLORS <- .80

option_list <- list(make_option(c("-i","--image"),type="character",default=NULL,help="Path to image file."),
                    make_option(c("-l","--list"),type="character",default=NULL,help="Path to single column .csv file listing paths to images."),
                    make_option(c("--folders"),type="character",default=NULL,help="Comma separated list of directories containing folders to convolve"),
                    make_option(c("--centersfolders"),type="character",default=NULL,help="Path to single column .csv file listing paths to images."),
                    make_option(c("--centerslist"),type="character",default=NULL,help="Path to single column .csv file listing paths to images for centers."),
                    make_option(c("-d","--directory"),type="character",default="",help="Path to directory where results should be saved."),
					          make_option(c("-t","--training"),type="character",default="",help="Path to pretrained machine learning results."),
                    make_option(c("--summaries"),type="character",default=NULL,help="Path to pretrained summaries."),
                    make_option(c("-n","--name"),type="character",default=NULL,help="Base name for saved images. Default is name stem from training RData file."),
                    make_option(c("--svm"),type="numeric",default=0,help="If non-0, convolves images using SVM classifier."),
                    make_option(c("--pca"),type="numeric",default=0,help="If non-0, convolves images given PCA dimension."),
                    make_option(c("--dot"),type="numeric",default=0,help="If non-0, convolves images using landscape norm."),
                    make_option(c("--patches"),type="numeric",default=0,help="Number of patches to use per image. Defaults to number from training RData file."),
                    make_option(c("-r","--radius"),type="numeric",default=NULL,help="Radius for patches. Defaults to radius from training file."),
                    make_option(c("-p","--proportion"),type="numeric",default=.10,help="Proportion of patches to use for samples. Defaults to proportion from training file."),
                    make_option(c("-c","--cores"),type="numeric",default=1,help="Number of cores to use for parallelized portions."), 
                    make_option(c("--threshold"),type="numeric",default=100,help="Percentage T between 0 and 100. Only the top T percent of mask values will be displayed."),
                    make_option(c("--negate"),type="numeric",default=0,help="Flip positive and negative scores."),
                    make_option(c("--separate"),type="numeric",default=0,help="If 1, output separate images for each input. If 0, masks will be in a single image file."))
opt <- parse_args(OptionParser(option_list=option_list))

real_opt <- opt

if(!is.null(opt$image)) { 
  list_of_images_to_mask <- list(opt$image)
} else if(!is.null(opt$list)) {
  list_of_images_to_mask <- as.list(read.csv(opt$list,header=FALSE)[,1])
} else if(!is.null(opt$folders)) { 
  list_of_folders <- strsplit(opt$folders,",")
  print(list_of_folders)
  file_names <- list()
  i <- 1
  for(folder in list_of_folders) {
    file_names[[i]] <- list.files(folder,full.names=TRUE)
    i <- i+1
  }
  list_of_images_to_mask <- unlist(file_names,recursive=TRUE)
}
if(!is.null(opt$centerslist)) { 
  centers_list <- as.list(read.csv(opt$centerslist,header=FALSE)[,1])
  for(i in 1:length(centers_list)) { 
    centers_list[[i]] <- as.character(centers_list[[i]])
  }
} else if(!is.null(opt$centersfolders)) { 
    list_of_folders <- strsplit(opt$centersfolders,",")
    print(list_of_folders)
    file_names <- list()
    i <- 1
    for(folder in list_of_folders) {
      file_names[[i]] <- list.files(folder,full.names=TRUE)
    }
    centers_list <- unlist(file_names,recursive=TRUE)
} else { 
    centers_list <- list()
    centers_list[1:length(list_of_images_to_mask)] <- FALSE
}

number_of_images <- length(list_of_images_to_mask)
if(opt$negate!=0) {print("NEGATIVE")}

extract_top_pca <- function(data_matrix,ml_results,types=NULL) { 
  rval <- list()
  landscapes_pca <- data_matrix%*%ml_results$landscapes_svd$v
  rval$data <- landscapes_pca[,c(ml_results$first_top_feature,ml_results$second_top_feature)]
  rval$types <- types
  rval$dim <- 2
  rval$name <- "PCA_plot_"
  rval$axis_labels <- c(paste("PCA",ml_results$first_top_feature),paste("PCA",ml_results$second_top_feature))
  return(rval)
}

top_pca_for_convolution <- function(dimension_choice=1) { 
  return(
    function(data_matrix) { 
      return(extract_top_pca(data_matrix,ml_environment$ml_results)$data[,dimension_choice])
    }
  )
}

apply_svm_rotation <- function(data_matrix,ml_results,types=NULL,classifier=1) { 
  rval <- list()
  model_matrix <- ml_results$all_svm[[classifier]]$patch_svm_model$W
  weight_vector <- model_matrix[1,(1:(dim(model_matrix)[2]-1))]
  bias_term <- model_matrix[1,dim(model_matrix)[2]]
  rval <- data_matrix%*%(weight_vector) + bias_term
  return(rval)
}

svm_rotation_for_convolution <- function(data_matrix) { return(apply_svm_rotation(data_matrix,ml_environment$ml_results,classifier=1)) }

dot_prod <- function(data_matrix) {return(diag(data_matrix%*%t(data_matrix)))}

convolution_runs <- list() 

if(opt$svm) { 
  convolution_runs$svm <- list("func"=svm_rotation_for_convolution,name="svm_")
}
if(opt$pca) { 
  convolution_runs$pca <- list("func"=top_pca_for_convolution(opt$pca),name=paste("pca_dim",opt$pca,"_",sep=""))
}
if(opt$dot) { 
  convolution_runs$dot <- list("func"=dot_prod,name="norm_")
}


ml_environment <- local({load(opt$training);environment()})

if(!is.null(opt$name)) {
  data_name_stem <- opt$name
} else {
  data_name_stem <- ml_environment$data_name_stem
}

if(!is.null(opt$directory)) { 
  image_results_directory <- opt$directory
} else {
  image_results_directory <- "results_images"
}

if(!is.null(opt$radius)) { 
  patch_radius <- opt$radius
} else { 
  patch_radius <- ml_environment$opt$radius
}

if(opt$patches) { 
  patches_per_image <- opt$patches
} else { 
  patches_per_image <- ml_environment$patches_per_image
}

if((patches_per_image %% 2) == 1) { 
  number_of_colors_for_breakpoints <- (patches_per_image %/% 2) + 1
}  else { 
  number_of_colors_for_breakpoints <- patches_per_image %/% 2
}

source(file.path(getScriptPath(),'utility_functions.R'))

number_of_points_per_patch <- max(1000,floor(.13*3*opt$radius**2))

image_summaries <- list()

if(!is.null(opt$summaries)) { 
  summaries_environment <- local({load(opt$summaries);environment()})
  if(!is.null(summaries_environment$image_summaries)) { 
    image_summaries <- summaries_environment$image_summaries
  }
} else { 
  for(i in 1:number_of_images) {
    image_summaries[[i]] <- patch_landscapes_from_image_alpha_sparse(as.character(list_of_images_to_mask[[i]]),patches_per_image,patch_center_image=centers_list[[i]],pixel_radius_for_patches = patch_radius,proportion_of_patch_sparse=.025,number_of_cores=opt$cores,max_PH_threshold=-1,return_patches=TRUE)
    if(opt$negate==1) { 
      image_summaries[[i]]$data <- -image_summaries[[i]]$data
    }
  }
}




for(run in convolution_runs) { 
  scale_max <- 1.0
  scale_min <- 0.0
  mapped_data <- NULL
  gc()

  convolution_data <- list()
  if(!is.null(opt$summaries)) { 
    convolution_data <- summaries_environment$convolution_data
  } else { 
    for(i in 1:number_of_images) { 
      convolution_data[[i]] <- weight_image_using_landscapes_and_transform(as.character(list_of_images_to_mask[[i]]),image_summaries[[i]]$data,image_summaries[[i]]$patches,patch_radius,run$func,min_weight=scale_min,max_weight=scale_max)
      convolution_data[[i]]$name <- as.character(list_of_images_to_mask[[i]])
      convolution_data[[i]]$type <- paste(run$name,unlist(strsplit(as.character(list_of_images_to_mask[[i]]),"/"))[7],sep="")
    }
  }
  z <- vector("double",number_of_images)
  a <- vector("double",number_of_images)
  w <- vector("double",number_of_images)
  for(i in 1:length(z)) { z[i] <- max(convolution_data[[i]]$mask);w[i] <- max(abs(convolution_data[[i]]$min),abs(convolution_data[[i]]$max))}
  for(i in 1:length(a)) { a[i] <- min(convolution_data[[i]]$mask)}
  max_value <- max(max(abs(z)),max(abs(a)))*MANUAL_SCALING_FOR_COLORS
  max_positive_value <- max(z) 
  min_negative_value <- min(a)

  image_grobs <- list()
  if(!opt$separate) { 
    for(i in 1:number_of_images) { 
      convolution_data[[i]]$mask[1,1] <- max_value
      convolution_data[[i]]$mask[1,2] <- -max_value
      convolution_data[[i]]$mask[which(convolution_data[[i]]$mask==0,arr.ind = TRUE)] <- -max(z)
      image_grobs[[i]] <- as.grob(function() {
        # mar=c(bottom,left,top,right)
        #par(cex=5.0,mar=c(0,0,1.0,0));
        par(mar=c(0,0,0,0));  
        image(convolution_data[[i]]$image,col=gray.colors(70000,start = 0.0),axes=FALSE);
        #image(convolution_data[[i]]$image,col=gray.colors(100000,start = 0.0),main=convolution_data[[i]]$type,axes=FALSE);
        image(convolution_data[[i]]$mask,col=divergingx_hcl(number_of_colors_for_breakpoints,palette = "RdBu",alpha=.60,rev=FALSE),add=TRUE)
      })
    }
  }
  break_points_for_histograms <- seq(-1.5,1.5,length.out=number_of_colors_for_breakpoints)
  hist_plots <- list()
  bar_colors<- list()
  for(i in 1:number_of_images) { 
    weights <- convolution_data[[i]]$weights
    # Have to set levels to full color range, otherwise they'll
    # be deleted automatically
    bar_colors <- factor(findInterval(weights,break_points_for_histograms),levels=1:number_of_colors_for_breakpoints)
    # Note drop=FALSE and nmax = number of breakpoints are required to keep histograms on same color scale
    if(i==1 || opt$separate) { 
      hist_plots[[i]] <- ggplotGrob(ggplot() + 
      geom_histogram(mapping=aes(V1,fill=bar_colors),data=as.data.frame(weights),bins=length(weights)/10,show.legend=FALSE,inherit.aes = FALSE) + 
      scale_fill_discrete_divergingx(palette="RdBu",drop=FALSE,nmax=number_of_colors_for_breakpoints) +
      geom_vline(aes(xintercept=mean(weights)),linetype="dashed",size=4.0,show.legend=FALSE) +
      ylab("Frequency") + 
      xlab("Patch value") +  
      theme_dark() +
      theme(text=element_text(size=85)) + 
      scale_x_continuous(breaks=c(-1,0,1),labels=c(paste("-1\n",unique(ml_environment$class_names)[1],sep=""),0,paste("1\n",unique(ml_environment$class_names)[2],sep="")),limits=c(min(break_points_for_histograms),max(break_points_for_histograms))))      
    } else { 
      hist_plots[[i]] <- ggplotGrob(ggplot() + 
                                    geom_histogram(mapping=aes(V1,fill=bar_colors),data=as.data.frame(weights),bins=length(weights)/10,show.legend=FALSE,inherit.aes = FALSE) + 
                                    scale_fill_discrete_divergingx(palette="RdBu",drop=FALSE,nmax=number_of_colors_for_breakpoints) +
                                    geom_vline(aes(xintercept=mean(weights)),linetype="dashed",size=4.0,show.legend=FALSE) +
                                    ylab("") + 
                                    xlab("Patch value") +
                                    theme_dark() +
                                    theme(text=element_text(size=85),plot.margin=margin(rep(0,4),"null")) +
                                    scale_x_continuous(breaks=c(-1,0,1),labels=c(paste("-1\n",unique(ml_environment$class_names)[1],sep=""),0,paste("1\n",unique(ml_environment$class_names)[2],sep="")),limits=c(min(break_points_for_histograms),max(break_points_for_histograms))))      
      }
  }

  if(!opt$separate) { 
    png(file.path(image_results_directory,paste('mask_plot_',run$type,"_",data_name_stem,'.png',sep="")),height=1.5*2048,width=(2048)*length(list_of_images_to_mask))
    layout <- matrix(0,ncol=number_of_images,nrow=3)
    for(i in 1:number_of_images) { 
        layout[,i] <- c(i,i,number_of_images+i)
    }
    grid.arrange(grobs=c(image_grobs,hist_plots),layout_matrix=layout,padding=unit(0.0,"null"))
    dev.off()
  } else { 
    for(i in 1:number_of_images) { 
      tiff(file.path(image_results_directory,paste('mask_plot_',run$type,"_",data_name_stem,i,'.tiff',sep="")),height=2048,width=2048)
      par(mar=c(0,0,0,0));  
      convolution_data[[i]]$mask[1,1] <- max_value
      convolution_data[[i]]$mask[1,2] <- -max_value
      convolution_data[[i]]$mask[which(convolution_data[[i]]$image==0,arr.ind = TRUE)] <- -max_value
      convolution_data[[i]]$image[which( (convolution_data[[i]]$mask)<((1 - opt$threshold/100)*max_positive_value)&(convolution_data[[i]]$mask >= 0),arr.ind = TRUE)] <- 0
      convolution_data[[i]]$image[which( (convolution_data[[i]]$mask)>((1 - opt$threshold/100)*min_negative_value)&(convolution_data[[i]]$mask <= 0),arr.ind = TRUE)] <- 0
      convolution_data[[i]]$mask[which( (convolution_data[[i]]$mask)<((1 - opt$threshold/100)*max_positive_value)&(convolution_data[[i]]$mask >= 0),arr.ind = TRUE)] <- -max_value
      convolution_data[[i]]$mask[which( (convolution_data[[i]]$mask)>((1 - opt$threshold/100)*min_negative_value)&(convolution_data[[i]]$mask <= 0),arr.ind = TRUE)] <- -max_value
      image(rotate(convolution_data[[i]]$image),col=gray.colors(70000,start = 0.0),axes=FALSE);
      image(rotate(convolution_data[[i]]$mask),breaks=seq(-max_value,max_value,length.out=number_of_colors_for_breakpoints+1),col=divergingx_hcl(number_of_colors_for_breakpoints,palette = "RdBu",alpha=.60,rev=FALSE),add=TRUE)
      dev.off()
      png(file.path(image_results_directory,paste('hist_plot_',run$type,data_name_stem,i,'.png',sep="")),height=1024,width=2048)
      grid.draw(hist_plots[[i]])
      dev.off()
    }
  }
}


data_file_name <- paste(data_name_stem,".RData",sep="")
save(convolution_data,opt,list_of_images_to_mask,file=file.path(opt$directory,data_file_name),version=2)