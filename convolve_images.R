library(optparse) 
library(utils)
library(colorspace) 
library(ggplot2)
library(gridExtra)
library(ggplotify)
library(grid)

MANUAL_SCALING_FOR_COLORS <- 1.0

option_list <- list(make_option(c("--training"),type="character",default="",help="Path to pretrained machine learning results."),
                    make_option(c("--image"),type="character",default=NULL,help="Path to a single image file to convolve."),
                    make_option(c("--folders"),type="character",default=NULL,help="Comma separated list of directories containing folders to convolve"),
                    make_option(c("--list"),type="character",default=NULL,help="Path to single column .csv file listing paths to images."),
                    make_option(c("--centersfolders"),type="character",default=NULL,help="Comma separated list of directories containing images from which to sample patch centers in convolution."),
                    make_option(c("--centerslist"),type="character",default=NULL,help="Path to single column .csv file listing paths to images for centers."),                    
                    make_option(c("--directory"),type="character",default=NULL,help="Path to directory where results should be saved. Default is ./tda-explore-convolutions."),
                    make_option(c("--name"),type="character",default=NULL,help="Base name for saved images. Default is name stem from training RData file."),
                    make_option(c("--svm"),type="character",default=FALSE,help="If non-0, convolves images using SVM classifier."),
                    make_option(c("--tsne"),type="character",default=FALSE,help="If non-0, convolves images using TSNE scores from input images."),
                    make_option(c("--patches"),type="numeric",default=0,help="Number of patches to use per image. Defaults to number from training RData file."),
                    make_option(c("--radius"),type="numeric",default=NULL,help="Radius for patches. Defaults to radius from training file."),
                    make_option(c("--cores"),type="numeric",default=1,help="Number of cores to use for parallelized portions."), 
                    make_option(c("--threshold"),type="numeric",default=100,help="Percentage T between 0 and 100. Only the top T percent of mask values will be displayed."),
                    make_option(c("--negate"),type="numeric",default=0,help="Flip positive and negative scores."),
                    make_option(c("--separate"),type="character",default=FALSE,help="If 1, output separate images for each input. If 0, masks will be in a single image file."))
opt <- parse_args(OptionParser(option_list=option_list))

real_opt <- opt

if(!is.null(opt$image)) { 
  list_of_images_to_mask <- list(opt$image)
} else if(!is.null(opt$list)) {
  list_of_images_to_mask <- as.list(read.csv(opt$list,header=FALSE)[,1])
} else if(!is.null(opt$folders)) { 
  list_of_folders <- strsplit(opt$folders,",")
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

if(opt$training!="") { 
  ml_environment <- local({load(opt$training);environment()})
}

apply_svm_rotation <- function(data_matrix,ml_results,types=NULL,classifier=1) { 
  rval <- list()
  model_matrix <- ml_results$svm[[classifier]]$patch_svm_model$W
  weight_vector <- model_matrix[1,(1:(dim(model_matrix)[2]-1))]
  bias_term <- model_matrix[1,dim(model_matrix)[2]]
  rval <- data_matrix%*%(weight_vector) + bias_term
  return(rval)
}

svm_rotation_for_convolution <- function(data_matrix) { return(apply_svm_rotation(data_matrix,ml_environment$training_results,classifier=1)) }

# Apply 1-dim t-sne to landscapes
# Center to 0, scale to []
apply_tsne <- function(data_matrix) { 
  output <- Rtsne::Rtsne(SparseM::as.matrix(data_matrix),dims=1,normalize=FALSE)$Y
  output <- as.vector(output,mode="double")
  # Center
  output <- output-mean(output)
  # Scale 
  output <- (output/max(abs(output)))
  output <- (output > 0)+ -(output <= 0)
  return(output)
}

convolution_runs <- list() 

if(opt$svm!=FALSE) { 
  convolution_runs$svm <- list("func"=svm_rotation_for_convolution,name="svm_","first_class"=unique(ml_environment$training_results$class_names)[2],"second_class"=unique(ml_environment$training_results$class_names)[1])
  if(opt$negate!=0) { 
    temp <- convolution_runs$svm$first_class
    convolution_runs$svm$first_class <- convolution_runs$svm$second_class
    convolution_runs$svm$second_class <- temp
  }
}

if(opt$tsne!=FALSE) { 
  convolution_runs$tsne <- list("func"=apply_tsne,name="tsne_","first_class"="Low","second_class"="High")
}


if(!is.null(opt$name)) {
  data_name_stem <- opt$name
} else {
  data_name_stem <- ml_environment$training_results$data_name_stem
}

if(!is.null(opt$directory)) { 
  image_results_directory <- opt$directory
} else {
  image_results_directory <- "tda-explore-convolutions"
}
if(!dir.exists(image_results_directory)) { 
  warning(paste("Results directory does not exist, creating ",image_results_directory,sep=""))
  dir.create(image_results_directory)
}

if(!is.null(opt$radius)) { 
  patch_radius <- opt$radius
} else { 
  patch_radius <- ml_environment$training_results$radius
}

if(opt$patches) { 
  patches_per_image <- opt$patches
} else { 
  patches_per_image <- ml_environment$training_results$patches_per_image
}

if((patches_per_image %% 2) == 1) { 
  number_of_colors_for_breakpoints <- (patches_per_image %/% 2) + 1
}  else { 
  number_of_colors_for_breakpoints <- patches_per_image %/% 2
}

image_summaries <- list()
for(i in 1:number_of_images) {
  remember <- getOption("warn")
  options(warn=-1)
  image_summaries[[i]] <- TDAExplore::patch_landscapes_from_image(as.character(list_of_images_to_mask[[i]]),patches_per_image,patch_center_image=centers_list[[i]],pixel_radius_for_patches = patch_radius,number_of_cores=opt$cores,return_patches=TRUE)
  options(warn=remember)
  if(opt$negate!=0) { 
    image_summaries[[i]]$data <- -image_summaries[[i]]$data
  }
}





for(run in convolution_runs) { 
  scale_max <- 1.0
  scale_min <- 0.0

  convolution_data <- list()
  for(i in 1:number_of_images) { 
    remember <- getOption("warn")
    options(warn=-1)
    convolution_data[[i]] <- TDAExplore::weight_image_using_landscapes_and_transform(as.character(list_of_images_to_mask[[i]]),image_summaries[[i]]$data,image_summaries[[i]]$patches,patch_radius,run$func,min_weight=scale_min,max_weight=scale_max)
    options(warn=remember)
    convolution_data[[i]]$name <- as.character(list_of_images_to_mask[[i]])
    convolution_data[[i]]$type <- paste(run$name,unlist(strsplit(as.character(list_of_images_to_mask[[i]]),"/"))[7],sep="")
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
  if(opt$separate==FALSE) { 
    for(i in 1:number_of_images) { 
      convolution_data[[i]]$mask[1,1] <- max_value
      convolution_data[[i]]$mask[1,2] <- -max_value
      convolution_data[[i]]$mask[which(convolution_data[[i]]$mask==0,arr.ind = TRUE)] <- -max_value-1
      image_grobs[[i]] <- as.grob(function() {
        # mar=c(bottom,left,top,right)
        #par(cex=5.0,mar=c(0,0,1.0,0));
        par(mar=c(0,0,0,0));  
        image(convolution_data[[i]]$image,col=c("#00000000"),zlim=c(0,max(convolution_data[[i]]$image)),axes=FALSE)
        image(convolution_data[[i]]$image,col=gray.colors(70000,start = 0.0),zlim=c(0.1,max(convolution_data[[i]]$image)),axes=FALSE,add=TRUE)    
        image(convolution_data[[i]]$mask,zlim=c(-max_value,max_value),col=divergingx_hcl(number_of_colors_for_breakpoints,palette = "RdBu",alpha=.60,rev=FALSE),add=TRUE)
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
    if(i==1 || opt$separate!=FALSE) { 
      hist_plots[[i]] <- ggplotGrob(ggplot() + 
      geom_histogram(mapping=aes(weights,fill=bar_colors),data=as.data.frame(weights),bins=length(weights)/10,show.legend=FALSE,inherit.aes = FALSE) + 
      scale_fill_discrete_divergingx(palette="RdBu",drop=FALSE,nmax=number_of_colors_for_breakpoints) +
      geom_vline(aes(xintercept=mean(weights)),linetype="dashed",size=4.0,show.legend=FALSE) +
      ylab("Frequency") + 
      xlab("Patch value") +  
      theme_dark() +
      theme(text=element_text(size=85)) + 
      scale_x_continuous(breaks=c(-1,0,1),labels=c(paste("-1\n",run$first_class,sep=""),0,paste("1\n",run$second_class,sep="")),limits=c(min(break_points_for_histograms),max(break_points_for_histograms))))      
    } else { 
      hist_plots[[i]] <- ggplotGrob(ggplot() + 
                                    geom_histogram(mapping=aes(weights,fill=bar_colors),data=as.data.frame(weights),bins=length(weights)/10,show.legend=FALSE,inherit.aes = FALSE) + 
                                    scale_fill_discrete_divergingx(palette="RdBu",drop=FALSE,nmax=number_of_colors_for_breakpoints) +
                                    geom_vline(aes(xintercept=mean(weights)),linetype="dashed",size=4.0,show.legend=FALSE) +
                                    ylab("") + 
                                    xlab("Patch value") +
                                    theme_dark() +
                                    theme(text=element_text(size=85),plot.margin=margin(rep(0,4),"null")) +
                                    scale_x_continuous(breaks=c(-1,0,1),labels=c(paste("-1\n",run$first_class,sep=""),0,paste("1\n",run$second_class,sep="")),limits=c(min(break_points_for_histograms),max(break_points_for_histograms))))      
      }
  }

  if(opt$separate==FALSE) { 
    png(file.path(image_results_directory,paste('mask_plot_',run$name,"_",data_name_stem,'.png',sep="")),height=1.5*2048,width=(2048)*length(list_of_images_to_mask),bg="transparent")
    layout <- matrix(0,ncol=number_of_images,nrow=3)
    for(i in 1:number_of_images) { 
        layout[,i] <- c(i,i,number_of_images+i)
    }
    grid.arrange(grobs=c(image_grobs,hist_plots),layout_matrix=layout,padding=unit(0.0,"null"))
    dev.off()
  } else { 
    for(i in 1:number_of_images) { 
      png(file.path(image_results_directory,paste('mask_plot_',run$name,"_",data_name_stem,i,'.png',sep="")),bg="transparent",height=2048,width=2048)
      par(mar=c(0,0,0,0));  
      convolution_data[[i]]$mask[1,1] <- max_value
      convolution_data[[i]]$mask[1,2] <- -max_value
      convolution_data[[i]]$mask[which(convolution_data[[i]]$image==0,arr.ind = TRUE)] <- -max_value-1
      convolution_data[[i]]$image[which( (convolution_data[[i]]$mask)<((1 - opt$threshold/100)*max_positive_value)&(convolution_data[[i]]$mask >= 0),arr.ind = TRUE)] <- 0
      convolution_data[[i]]$image[which( (convolution_data[[i]]$mask)>((1 - opt$threshold/100)*min_negative_value)&(convolution_data[[i]]$mask <= 0),arr.ind = TRUE)] <- 0
      convolution_data[[i]]$mask[which( (convolution_data[[i]]$mask)<((1 - opt$threshold/100)*max_positive_value)&(convolution_data[[i]]$mask >= 0),arr.ind = TRUE)] <- -max_value
      convolution_data[[i]]$mask[which( (convolution_data[[i]]$mask)>((1 - opt$threshold/100)*min_negative_value)&(convolution_data[[i]]$mask <= 0),arr.ind = TRUE)] <- -max_value
      image(convolution_data[[i]]$image,col=c("#00000000"),zlim=c(0,max(convolution_data[[i]]$image)),axes=FALSE)
      image(convolution_data[[i]]$image,col=gray.colors(70000,start = 0.0),zlim=c(0.1,max(convolution_data[[i]]$image)),axes=FALSE,add=TRUE)
      image(convolution_data[[i]]$mask,zlim=c(-max_value,max_value),col=divergingx_hcl(number_of_colors_for_breakpoints,palette = "RdBu",alpha=.60,rev=FALSE),add=TRUE)
      dev.off()
      png(file.path(image_results_directory,paste('hist_plot_',run$name,data_name_stem,i,'.png',sep="")),height=1024,width=2048)
      grid.draw(hist_plots[[i]])
      dev.off()
    }
  }
}


data_file_name <- paste(data_name_stem,".RData",sep="")
save(convolution_data,opt,list_of_images_to_mask,file=file.path(image_results_directory,data_file_name),version=2)