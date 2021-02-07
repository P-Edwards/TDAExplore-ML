option_list <- list(optparse::make_option(c("--parameters"),type="character",default=NULL,help="File name of parameters .csv file to use. Relative to current directory or absolute."),
                    optparse::make_option(c("--cores"),type="numeric",default=0,help="Number of images to process in parallel. Default is 2."),
                    optparse::make_option(c("--radius"),type="numeric",default=NULL,help="Radius for patches. Higher priority than parameters file."),
                    optparse::make_option(c("--ratio"),type="numeric",default=2,help="Patch ratio. Default is 2."),
                    optparse::make_option(c("--svm"),type="numeric",default=FALSE,help="If TRUE, 5-fold SVM cross-validaton on patch landscapes will be performed."),
                    optparse::make_option(c("--randforest"),type="numeric",default=FALSE,help="Experimental. If TRUE, PCA and random forest ranking will be computed."),
                    optparse::make_option(c("--avgsvm"),type="numeric",default=FALSE,help="Experimental. If TRUE, 5-fold SVM cross-validaton performed by first averaging patch landscapes summaries for images."),
                    optparse::make_option(c("--lower"),type="numeric",default=0,help="Experimental. Minimum percentage distance from centroid for patches."),
                    optparse::make_option(c("--upper"),type="numeric",default=1,help="Experimental. Maximum percentage distance from centroid for patches."),
                    optparse::make_option(c("--fullsave"),type="numeric",default=TRUE,help="Debugging. If flag is set to FALSE or 0, then no full save is done."),
                    optparse::make_option(c("--benchmark"),type="numeric",default=FALSE,help="Debugging. If TRUE, uses a fork based cluster type for benchmarking. Does not work on Windows."))
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
if(is.null(opt$parameters)) { 
  print("Using default parameters file, <current_working_directory>/data_parameters.csv")
  parameters_file_name <- "data_parameters.csv"
} else { 
  parameters_file_name <- opt$parameters
}

data_parameters <- utils::read.csv(parameters_file_name)


provided_parameters <- colnames(data_parameters)
# Default behavior: Use only working directory
if("image_directories" %in% provided_parameters) {
  image_type_names <- data_parameters[,"image_directories"]
} else { 
  all.files <- list.files(rec=F)
  all.files[!file.info(all.files)$isdir]
}
# Default behavior: Every directory is in different class
if("directory_classes" %in% provided_parameters) { 
  class_names <- data_parameters[,"directory_classes"]
} else { 
  print("No class names for input directories specified. Using directory names.")
  class_names <- image_type_names
}
if("patch_center_images_directories" %in% provided_parameters) { 
  patch_center_names <- data_parameters[,"patch_center_images_directories"]
  patch_center_flag <- TRUE
} else { patch_center_flag <- FALSE }


image_file_names_by_directory <- list()
patch_center_file_names_by_directory <- list()

total_number_of_files <- 0
for(i in 1:length(image_type_names)) { 
  image_file_names_by_directory[[i]] <- list.files(levels(image_type_names)[image_type_names[i]],full.names=TRUE)
  if(patch_center_flag) { 
    patch_center_file_names_by_directory[[i]] <- list.files(levels(patch_center_names)[patch_center_names[i]],full.names=TRUE)
  } else { 
    # instead sets up a vector of appropriate length that is all FALSE
    patch_center_file_names_by_directory[[i]] <- as.vector(image_file_names_by_directory[[i]],mode="logical")
    patch_center_file_names_by_directory[[i]][1:length(patch_center_file_names_by_directory[[i]])] <- FALSE
  }
  total_number_of_files <- total_number_of_files + length(image_file_names_by_directory[[i]])
}

# Default behavior: Pixel radius is 2.5% of image length.
if(!is.null(opt$radius)) { 
    radius_of_patches <- opt$radius
  } else if("radius_of_patches" %in% provided_parameters) { 
      radius_of_patches <- data_parameters[1,"radius_of_patches"]
  } else { 
      print("No pixel radius for patches specified. Using default, which is 2.5% of image length.")
      first_image <- OpenImageR::readImage(unlist(image_file_names_by_directory)[2])
      radius_of_patches <- floor(0.025*dim(first_image)[1])
    }

# Default behavior: Use 2 cores
if(opt$cores) { 
  number_of_cores <- opt$cores
} else if("number_of_cores" %in% provided_parameters) { 
  number_of_cores <- data_parameters[1,"number_of_cores"]
} else { 
  number_of_cores <- 2
}

# Default behavior: Calculation should take at most 30s
image_for_dimensions <- OpenImageR::readImage(image_file_names_by_directory[[1]][1])
pixel_area_of_images <- nrow(image_for_dimensions)*ncol(image_for_dimensions)
patch_ratio <- opt$ratio
patches_per_image <- floor(patch_ratio*pixel_area_of_images/(2*radius_of_patches)**2)
print(patches_per_image)

if("experiment_name" %in% provided_parameters) {
  experiment_name <- data_parameters[1,"experiment_name"]
} else { 
  experiment_name <- floor(ruinf(1,min=0,max=1e6))
}


if("data_results_directory" %in% provided_parameters) { 
  data_results_directory <- file.path(data_parameters[1,"data_results_directory"])
  if(!dir.exists(data_results_directory)) { 
    dir.create(data_results_directory)
  }
} else { 
  data_results_directory <- ""  
}


data_name_stem <- paste(experiment_name,patches_per_image,"patches_",radius_of_patches,"radius_",runif(1,min=1,max=1e9),sep = "")
if(opt$upper - opt$lower < 1) { 
  data_name_stem <- paste(data_name_stem,"_radialthresh")
}


i <- 1
type_vector <- vector("character",length=total_number_of_files*patches_per_image)
offset <- 0
for(image_list in image_file_names_by_directory) { 
  type_vector[(1+offset):(offset+length(image_list)*patches_per_image)] <- class_names[[i]]
  offset <- offset + length(image_list)*patches_per_image
  i <- i+1 
}


image_name_iterator <- iterators::iter(cbind(unlist(image_file_names_by_directory),unlist(patch_center_file_names_by_directory)),by="row")

if(!opt$benchmark) {
  cl <- parallel::makeCluster(number_of_cores,outfile="")
  doParallel::registerDoParallel(cl)
  } else { 
    # Forking using DoMC is picked up by SLURM
    # Only available on Linux
    doMC::registerDoMC(cores=number_of_cores)
  }
`%dopar` <- foreach::`%dopar%`
unscrambled_data <- foreach::foreach(image_name=image_name_iterator,.combine=rbind,.packages=c("TDAExplore")) %dopar% { 
	options(warn=-1)
  print(paste("Started image ", image_name))    
  if(image_name[2]!=FALSE) {
    print(paste("With image for patch centers: ",patch_center_image))
  }
	patch_summaries <- TDAExplore::patch_landscapes_from_image(image_name[1],patches_per_image,patch_center_image=image_name[2],pixel_radius_for_patches = radius_of_patches,proportion_of_patch_sparse=.025,max_PH_threshold=-1,lower_threshold=opt$lower,upper_threshold=opt$upper) 
  print(paste("Started image ", image_name))    
	patch_summaries
}
if(!opt$benchmark) {
  parallel::stopCluster(cl)
}

if(length(type_vector)>=2147483647) { 
  unscrambled_data <- as.matrix(do.call(rbind,unscrambled_data))
} else { 
  unscrambled_data <- as.matrix.csr(do.call(rbind,unscrambled_data))
}


image_file_names <- unlist(image_file_names_by_directory)

ml_results <- list()
ml_results$all_svm <- list()

if(opt$randforest) {
  print("Starting random forest and PCA")
  library(RSpectra)
  ml_results$random_forest <- list()
  # Use PCA to dimensionality reduce to top 50 principal components
  ml_results$random_forest$landscapes_svd <-svds(as.matrix(unscrambled_data),k=50,nu=0,nv=50,center=FALSE,scale=FALSE,tol=NULL)
  landscapes_pca<-unscrambled_data%*%(ml_results$landscapes_svd$v)

  # Train random forest on reduced vectors
  library(randomForest)
  ml_results$random_forest$model <- randomForest(as.matrix(landscapes_pca),factor(type_vector),ntree=500,mtry=7)

  # Use top 2 features from random forest training for summary
  ml_results$random_forest$first_top_feature <- order(-ml_results$rf_model$importance)[1]
  ml_results$random_forest$second_top_feature <- order(-ml_results$rf_model$importance)[2]
}
if(opt$svm) {
  print("Starting per-landscape SVM")
  library(LiblineaR)
  # Cross validation on PCA-rotated patch landscapes
  number_of_validation_steps <- 5
  
  #randomly shuffle the data
  shuffled_order <- sample(length(image_file_names))
  shuffled_image_names <- image_file_names[shuffled_order]

  shuffled_patch_indices <- vector(mode="double",length=length(type_vector))
  image_folds <- cut(seq(1,length(shuffled_image_names)),breaks=number_of_validation_steps,labels=FALSE)
  folds <- vector(mode="double",length=length(type_vector))
  for(i in 1:length(shuffled_order)) { 
    image_index <- shuffled_order[i]
    shuffled_patch_indices[((i-1)*patches_per_image+1):(i*patches_per_image)] <- ((image_index-1)*patches_per_image+1):(image_index*patches_per_image)
    folds[((i-1)*patches_per_image+1):(i*patches_per_image)] <- image_folds[i]
  }

  # First: Use 20% of the data to tune SVM parameters.  

  shuffled_pca <- unscrambled_data[shuffled_patch_indices,]
  shuffled_types <- type_vector[shuffled_patch_indices]
  trainIndexes <- which(folds==1,arr.ind=TRUE)

  # Quick heuristic parameter tuning
  cost <- heuristicC(as.matrix(shuffled_pca[trainIndexes,]))



  # Now do cross validation with the remaining data and selected parameters
  reduced_data <- shuffled_pca
  reduced_types <- shuffled_types

  shuffled_pca <- NULL
  shuffled_types <- NULL
  gc()

  patch_accuracies <- vector("double",number_of_validation_steps)
  image_accuracies <- vector("double",number_of_validation_steps)

  SVM_file_path <- file.path(data_results_directory,paste(data_name_stem,"_SVM_cross_validation.csv",sep=""))
  SVM_avg_per_patch_file_path <- file.path(data_results_directory,paste(data_name_stem,"_patch_feature_average_image_SVM_cross_validation.csv",sep=""))
  proportion_vector <- vector("double",length=max(reduced_types))
  for(i in 1:max(reduced_types)) { 
    proportion_vector[i] <- sum(reduced_types==i)/length(reduced_types)
  }
  write.csv(data.frame(proportion_vector,row.names=levels(class_names)),file=SVM_file_path,append=FALSE)
  write.csv(data.frame(proportion_vector,row.names=levels(class_names)),file=SVM_avg_per_patch_file_path,append=FALSE)

  reduced_types_unique <- unique(type_vector)
  for(i in 1:number_of_validation_steps) { 
    print(paste("Starting cross validation fold number ",i))
    testIndexes <- which(folds==i,arr.ind=TRUE)
    trainIndexes <- -testIndexes
    weights <- vector(mode="double",length=length(unique(reduced_types[trainIndexes])))
    names(weights) <- unique(reduced_types[trainIndexes]) 
    for(j in 1:length(weights)) { 
      weights[j] <- sum(reduced_types[trainIndexes]==names(weights)[j])
    }
    max_number <- max(weights)
    for(j in 1:length(weights)) { 
      weights[j] <- max_number/weights[j]
    }
    reduced_types_training <- reduced_types[trainIndexes]
    train_data <- reduced_data[trainIndexes,]
    if(reduced_types_training[1]!=reduced_types_unique[1]) { 
      switch_index <- min(which(reduced_types_training==reduced_types_unique[1]))
      reduced_types_training[1] <- reduced_types_unique[1]
      reduced_types_training[switch_index] <- reduced_types_unique[2]
      temporary_row <- train_data[switch_index,]
      train_data[switch_index,] <- train_data[1,]
      train_data[1,] <- temporary_row
    }
    svm_model <- LiblineaR(data=train_data,target=factor(reduced_types_training),wi=weights,cost=cost,type=2)

    ml_results$all_svm[[i]] <- list()
    ml_results$all_svm[[i]]$testing_data <- reduced_data[testIndexes,]
    ml_results$all_svm[[i]]$testing_labels <- reduced_types[testIndexes]
    
    prediction_values <- predict(svm_model,reduced_data[testIndexes,])
    actual_names <- levels(class_names)
    predicted_names <- levels(class_names)
    for(j in 1:length(actual_names)) { 
      actual_names[j] <- paste(actual_names[j],"_actual")
      predicted_names[j] <- paste(predicted_names[j],"_predicted")
    }
    #pred_table <- table(factor(prediction_values$predictions,labels=predicted_names[1:length(unique(prediction_values$predictions))],levels=unique(reduced_types)),factor(reduced_types[testIndexes],labels=actual_names,levels=unique(reduced_types)))
    pred_table <- table(factor(prediction_values$predictions,labels=predicted_names,levels=unique(reduced_types)),factor(reduced_types[testIndexes],labels=actual_names,levels=unique(reduced_types)))
    write.table(pred_table,file=SVM_file_path,sep=",",dec=".",append=TRUE)
    patch_accuracies[i] <- sum(diag(pred_table))/sum(pred_table)
    
    # Classification with average transformed image data
    # i.e. 
    # (1) Transform all patch landscapes into single numbers using model
    # (2) Average all numbers across each image
    # (3) Train and test an SVM model for classifying images
    model_matrix <- svm_model$W
    weight_vector <- model_matrix[1,(1:(dim(model_matrix)[2]-1))]
    bias_term <- model_matrix[1,dim(model_matrix)[2]]
    transformed_data <- as.matrix.csr(reduced_data%*%(weight_vector) + bias_term,nrow=nrow(reduced_data),ncol=1)
    averaged_data <- average_vectors_for_images(transformed_data,patches_per_image,shuffled_image_names,reduced_types,1,1)
    transformed_data <- averaged_data$image_weights
    transformed_types <- averaged_data$image_types
    
    costIndexes <- which(image_folds==1,arr.ind=TRUE)
    testIndexes <- which(image_folds==i,arr.ind=TRUE)
    trainIndexes <- -testIndexes

    # Quick heuristic parameter tuning
    image_cost <- heuristicC(as.matrix(transformed_data[costIndexes,]))

    train_data <- as.matrix.csr(transformed_data[trainIndexes,],ncol=1)
    test_data <- as.matrix.csr(transformed_data[testIndexes,],ncol=1)

    weights <- vector(mode="double",length=length(unique(transformed_types[trainIndexes])))
    names(weights) <- unique(transformed_types[trainIndexes]) 
    for(j in 1:length(weights)) { 
      weights[j] <- sum(transformed_types[trainIndexes]==names(weights)[j])
    }
    max_number <- max(weights)
    for(j in 1:length(weights)) { 
      weights[j] <- max_number/weights[j]
    }

    transformed_types_training <- transformed_types[trainIndexes]
    if(transformed_types_training[1]!=reduced_types_unique[1]) { 
      switch_index <- min(which(transformed_types_training==reduced_types_unique[1]))
      transformed_types_training[1] <- reduced_types_unique[1]
      transformed_types_training[switch_index] <- reduced_types_unique[2]
      temporary_row <- train_data[switch_index,]
      train_data[switch_index,] <- train_data[1,]
      train_data[1,] <- temporary_row
    }
    save.image(file="test.RData",version=2)  
    image_svm_model <- LiblineaR(data=train_data,target=factor(transformed_types_training),cost=image_cost,wi=weights,type=2)
    prediction_values <- predict(image_svm_model,test_data)
    actual_names <- levels(class_names)
    predicted_names <- levels(class_names)
    for(j in 1:length(actual_names)) { 
      actual_names[j] <- paste(actual_names[j],"_actual")
      predicted_names[j] <- paste(predicted_names[j],"_predicted")
    }
    pred_table <- table(factor(prediction_values$predictions,labels=predicted_names,levels=unique(transformed_types)),factor(transformed_types[testIndexes],labels=actual_names,levels=unique(transformed_types)))
    write.table(pred_table,file=SVM_avg_per_patch_file_path,sep=",",dec=".",append=TRUE)
    image_accuracies[i] <- sum(diag(pred_table))/sum(pred_table)
    
    ml_results$all_svm[[i]]$svm_model <- svm_model
    ml_results$all_svm[[i]]$patch_svm_model <- svm_model
    ml_results$all_svm[[i]]$image_svm_model <- image_svm_model
    ml_results$all_svm[[i]]$svm_image_training_indices <- shuffled_order[trainIndexes]
    ml_results$all_svm[[i]]$svm_patch_accuracies <- patch_accuracies
    ml_results$all_svm[[i]]$svm_image_accuracies <- image_accuracies
    ml_results$all_svm[[i]]$svm_cost <- cost 
  }

  write.table(table(patch_accuracies),file=SVM_file_path,append=TRUE,sep=",",dec=".")
  write.table(table(c(mean(patch_accuracies))),file=SVM_file_path,append=TRUE,sep=",",dec=".")
  print("Radius")
  print(radius_of_patches)
  print("Average patch accuracies")
  print(mean(patch_accuracies))
  print("Average image accuracies")
  print(mean(image_accuracies))
  write.table(table(image_accuracies),file=SVM_avg_per_patch_file_path,append=TRUE,sep=",",dec=".")
  write.table(table(c(mean(image_accuracies))),file=SVM_avg_per_patch_file_path,append=TRUE,sep=",",dec=".")

}
if(opt$avgsvm) { 
  print("Starting average per image SVM")
  library(LiblineaR)
  # Cross validation on averages of patch landscapes for each mage

  # First: Use 20% of the data to tune SVM parameters. This data will need to be discarded otherwise.
  
  # randomly shuffle the data
  averaged_data <- average_vectors_for_images(unscrambled_data,patches_per_image,image_file_names,type_vector,1,dim(unscrambled_data)[2])
  shuffled_order <- sample(nrow(averaged_data$image_weights))
  folds <- cut(seq(1,nrow(averaged_data$image_weights)),breaks=5,labels=FALSE)

  shuffled_pca <- averaged_data$image_weights[shuffled_order,]
  shuffled_types <- averaged_data$image_types[shuffled_order]
  trainIndexes <- which(folds==1,arr.ind=TRUE)

  # Quick heuristic parameter tuning
  cost <- heuristicC(shuffled_pca[trainIndexes,])



  # Now do cross validation with the remaining data and selected parameters
  number_of_validation_steps <- 5

  # Throw out tuning data
  #reduced_data <- shuffled_pca[-trainIndexes,]
  reduced_data <- shuffled_pca
  reduced_types <- shuffled_types

  shuffled_pca <- NULL
  shuffled_types <- NULL
  gc()

  folds <- cut(seq(1,nrow(reduced_data)),breaks=number_of_validation_steps,labels=FALSE)
  accuracies <- vector("double",number_of_validation_steps)

  SVM_file_path <- file.path(data_results_directory,paste(data_name_stem,"_image_average_landscape_SVM_cross_validation.csv",sep=""))
  proportion_vector <- vector("double",length=max(reduced_types))
  for(i in 1:max(reduced_types)) { 
    proportion_vector[i] <- sum(reduced_types==i)/length(reduced_types)
  }
  write.csv(data.frame(proportion_vector,row.names=levels(class_names)),file=SVM_file_path,append=FALSE)
  for(i in 1:number_of_validation_steps) { 
    print(paste("Starting run number ",i))
    testIndexes <- which(folds==i,arr.ind=TRUE)
    trainIndexes <- -testIndexes
    weights <- vector(mode="double",length=length(unique(reduced_types[trainIndexes])))
    names(weights) <- unique(reduced_types[trainIndexes]) 
    for(j in 1:length(weights)) { 
      weights[j] <- sum(reduced_types[trainIndexes]==names(weights)[j])
    }
    max_number <- max(weights)
    for(j in 1:length(weights)) { 
      weights[j] <- max_number/weights[j]
    }
    svm_model <- LiblineaR(data=reduced_data[trainIndexes,],target=factor(reduced_types[trainIndexes]),cost=cost,wi=weights,type=2)
    prediction_values <- predict(svm_model,reduced_data[testIndexes,])
    actual_names <- levels(class_names)
    predicted_names <- levels(class_names)
    for(j in 1:length(actual_names)) { 
      actual_names[j] <- paste(actual_names[j],"_actual")
      predicted_names[j] <- paste(predicted_names[j],"_predicted")
    }
    pred_table <- table(factor(prediction_values$predictions,labels=predicted_names[1:length(unique(prediction_values$predictions))]),factor(reduced_types[testIndexes],labels=actual_names))
    write.table(pred_table,file=SVM_file_path,sep=",",dec=".",append=TRUE)
    accuracies[i] <- sum(diag(pred_table))/sum(pred_table)
  }

  write.table(table(accuracies),file=SVM_file_path,append=TRUE,sep=",",dec=".")
  write.table(table(c(mean(accuracies))),file=SVM_file_path,append=TRUE,sep=",",dec=".")

  ml_results$svm_avg_model <- svm_model
  ml_results$svm_avg_accuracies <- accuracies
  ml_results$svm_avg_cost <- cost
}

  



if(opt$svm==1) { 
  save(ml_results,data_name_stem,patches_per_image,image_file_names,type_vector,class_names,opt,version=2,file=file.path(data_results_directory,paste("ML_results_and_summaries_",data_name_stem,".RData",sep="")))
} else {
  save(ml_results,unscrambled_data,data_name_stem,patches_per_image,image_file_names,type_vector,class_names,opt,version=2,file=file.path(data_results_directory,paste("ML_results_and_summaries_",data_name_stem,".RData",sep="")))
}
  
if(opt$fullsave) {
  data_file_name <- paste(data_name_stem,".RData",sep="")
  save.image(file=file.path(data_results_directory,data_file_name),version=2)
}


