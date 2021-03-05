option_list <- list(optparse::make_option(c("--parameters"),type="character",help="File name of parameters .csv file to use. Relative to current directory or absolute. Required."),
                    optparse::make_option(c("--cores"),type="numeric",default=2,help="Number of images to process in parallel. Default is 2."),
                    optparse::make_option(c("--radius"),type="numeric",default=FALSE,help="Radius for patches. Higher priority than parameters file."),
                    optparse::make_option(c("--ratio"),type="numeric",default=2,help="Patch ratio. Default is 2."),
                    optparse::make_option(c("--svm"),type="character",default=FALSE,help="If TRUE, 5-fold SVM cross-validaton on patch landscapes will be performed."),
                    optparse::make_option(c("--landscapes"),type="character",default=TRUE,help="If TRUE, constructs landscapes for images."),
                    optparse::make_option(c("--plot"),type="character",default=FALSE,help="If TRUE, save plots summarizing ML computations."),
                    optparse::make_option(c("--summaries"),type="character",default=FALSE,help="Path to RData file previously saved by this application. Images and landscapes from run in RData file will be used for plotting."),
                    optparse::make_option(c("--training"),type="character",default=FALSE,help="Path to RData file previously saved by this application. SVM model will be used for plotting."),
                    optparse::make_option(c("--imagedir"),type="character",default=FALSE,help="Path to directory in which to save plots. Default is ./tda-explore-plots."),
                    optparse::make_option(c("--lower"),type="numeric",default=0,help="Experimental. Minimum percentage distance from centroid for patches."),
                    optparse::make_option(c("--upper"),type="numeric",default=1,help="Experimental. Maximum percentage distance from centroid for patches."),
                    optparse::make_option(c("--fullsave"),type="numeric",default=FALSE,help="Debugging. If flag is set to TRUE then uses save.image at end. Default FALSE."),
                    optparse::make_option(c("--benchmark"),type="numeric",default=FALSE,help="Debugging. If TRUE, uses a fork based cluster type for benchmarking. Does not work on Windows."))
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))


replace_falsy_with_false <- function(flag_to_check) { 
  list_of_falsy_objects <- c(FALSE,"FALSE",0,"False","no","No","F","f")
  for(entry in list_of_falsy_objects) { 
    if(flag_to_check==entry) { 
      return(FALSE)
    }
  }
  return(TRUE)
}

options_to_falsify <- c("svm","landscapes","plot","fullsave","benchmark")
for(command_line_option in options_to_falsify) { 
  opt[command_line_option]  <- replace_falsy_with_false(opt[command_line_option])
}

if(opt$landscapes!=FALSE) { 
  training_results <- TDAExplore::TDAExplore(parameters=opt$parameters,number_of_cores=opt$cores,radius_of_patches=opt$radius,patch_ratio=opt$ratio,svm=opt$svm,lower=opt$lower,upper=opt$upper,verbose=TRUE,benchmark=opt$benchmark)  
  save(training_results,version=2,file=file.path(training_results$data_results_directory,paste(training_results$data_name_stem,".RData",sep="")))
}

if(opt$plot!=FALSE) { 
  if(opt$summaries!=FALSE) { 
    extract <- local({load(opt$summaries);environment()})
    summaries_results <- extract$training_results
    remove(extract)
  } else { 
    summaries_results <- training_results 
  }

  if(opt$training!=FALSE) { 
    extract <- local({load(opt$training);environment()})
    training_results <- extract$training_results
    remove(extract)
  } 
  
  if(is.null(training_results$svm)) { 
    stop("There is no SVM training from the current run or from --training tag.")
  }

  apply_patch_svm_rotation <- function(data_matrix,labels,classifier=1) { 
    rval <- list()
    model_matrix <- training_results$svm[[classifier]]$patch_svm_model$W
    weight_vector <- model_matrix[1,(1:(dim(model_matrix)[2]-1))]
    bias_term <- model_matrix[1,dim(model_matrix)[2]]
    # SparseM and transposes are to coerce this into the correct
    # data type
    rval$data <- SparseM::t(weight_vector%*%SparseM::t(data_matrix)) + bias_term
    rval$types <- labels
    rval$named_types <- factor(rval$types,levels=unique(rval$types),labels=levels(summaries_results$class_names)[unique(rval$types)])
    rval$dim <- 1
    rval$name <- paste("SVM_plot_classifier_",classifier,sep="")
    rval$axis_labels <- "SVM classifier output"
    return(rval)
  }
  apply_image_svm_rotation <- function(data_matrix,classifier=1) { 
    model_matrix <- training_results$svm[[classifier]]$image_svm_model$W
    weight_vector <- model_matrix[1,(1:(dim(model_matrix)[2]-1))]
      bias_term <- model_matrix[1,dim(model_matrix)[2]]
      return(data_matrix%*%(weight_vector)+bias_term)
  }

  identity <- function(x) {return(x)}

  average_by_image <- function(data_matrix,transformation_map,max_dim=1,post_transformation_map=identity) { 
    rval <- list()
    transformed_data <- transformation_map(data_matrix) 
    average_data <- TDAExplore:::average_vectors_for_images(transformed_data$data,summaries_results$patches_per_image,summaries_results$image_file_names,transformed_data$types,1,transformed_data$dim)
    average_data$image_weights <- post_transformation_map(average_data$image_weights)
    rval$data <- -average_data$image_weights
    rval$name <- paste(transformed_data$name,"images_",sep="")
    rval$types <- average_data$image_types
    rval$named_types <- factor(rval$types,levels=1:length(levels(summaries_results$class_names)),labels=levels(summaries_results$class_names))
    rval$dim <- transformed_data$dim
    rval$axis_labels  <- transformed_data$axis_labels
    return(rval)
  }


  plotting_data <- list() 




  number_of_classifiers <- length(training_results$svm)

  individual_results <- list()
  individual_image_results <- list()
  testing_order_image_names <- vector("character")
  if(opt$summaries==FALSE && opt$svm!=FALSE) {      
    for(i in 1:number_of_classifiers) {     
      testing_order_image_names <- c(testing_order_image_names,training_results$image_file_names[-training_results$svm[[i]]$svm_image_training_indices])
      individual_results[[i]] <- apply_patch_svm_rotation(as.matrix(training_results$svm[[i]]$testing_data),summaries_results$svm[[i]]$testing_labels,classifier=i)
      individual_results[[i]]$data <- -individual_results[[i]]$data
      individual_image_results[[i]] <- average_by_image(training_results$svm[[i]]$testing_data,function(data_matrix) {return(apply_patch_svm_rotation(data_matrix,summaries_results$svm[[i]]$testing_labels,classifier=i))},post_transformation_map=function(data_matrix) {return(apply_image_svm_rotation(data_matrix,classifier=i))})
    }
  } else { 
    individual_results[[1]] <- apply_patch_svm_rotation(summaries_results$summaries,classifier=1)
    individual_results[[1]]$data <- -individual_results[[1]]$data
    individual_image_results[[1]] <- average_by_image(summaries_results$summaries,function(data_matrix) {return(apply_patch_svm_rotation(data_matrix,summaries_results$patch_types,classifier=1))},post_transformation_map=function(data_matrix) {return(apply_image_svm_rotation(data_matrix,classifier=1))})
  }
  full_plotting_result <- individual_results[[1]]
  full_plotting_image_result <- individual_image_results[[1]]
  
  full_plotting_result$name <- "Patch_SVM_scores"
  full_plotting_image_result$name <- "Image_SVM_scores"

  full_plotting_result$data <- vector(mode="double")
  full_plotting_image_result$data <- vector(mode="double")
  full_plotting_result$types <- vector(mode="character")
  full_plotting_image_result$types <- vector(mode="character")
  full_plotting_result$which_run <- vector(mode="double")
  full_plotting_image_result$which_run <- vector(mode="double")
  for(i in 1:length(individual_image_results)) { 
    full_plotting_result$data <- c(full_plotting_result$data,individual_results[[i]]$data)
    full_plotting_result$types <- c(full_plotting_result$types,individual_results[[i]]$types)
    temp <- individual_results[[i]]$types
    temp[1:length(temp)] <- i
    full_plotting_result$which_run <- c(full_plotting_result$which_run,temp)
    full_plotting_image_result$data <- c(full_plotting_image_result$data,individual_image_results[[i]]$data)
    full_plotting_image_result$types <- c(full_plotting_image_result$types,individual_image_results[[i]]$types)
    full_plotting_image_result$which_run <- c(full_plotting_image_result$which_run,rep.int(i,length(individual_image_results[[i]]$types)))
  }
  full_plotting_result$named_types <- factor(full_plotting_result$types,1:length(levels(summaries_results$class_names)),labels=levels(summaries_results$class_names))
  full_plotting_image_result$named_types <- factor(full_plotting_image_result$types,1:length(levels(summaries_results$class_names)),labels=levels(summaries_results$class_names))
  plotting_data$patches_svm <- full_plotting_result
  plotting_data$image_svm <- full_plotting_image_result
  


  library(ggplot2)
  if(opt$imagedir!=FALSE) { 
    image_results_directory <- file.path(opt$imagedir)
  } else { 
    image_results_directory <- "tda-explore-plots"  
  }
  if(!dir.exists(image_results_directory)) { 
    dir.create(image_results_directory)
  }

  if(opt$summaries!=FALSE) { 
    plotting_name_stem <- paste(training_results$data_name_stem,"_training_applied_to_",summaries_results$data_name_stem,"_summaries",sep="")
  } else { 
    plotting_name_stem <- paste(training_results$data_name_stem,"_plots",sep="")
  }

  for(plotting_details in plotting_data) { 
    if(is.null(plotting_details)) { next }
    if(plotting_details$dim == 1) { 
      boxplot_data <- data.frame(features=plotting_details$data,types=plotting_details$named_types)
      p <- ggplot(boxplot_data,aes(x=types,y=features))+xlab("Experimental groups")+ylab("Average topological score") + geom_boxplot(fill=NA)+theme_minimal()+scale_fill_manual(values=c("blue","red","green","purple","turquoise","orange"))+scale_color_manual(values=c("blue","red","green","purple","turquoise","orange"))+theme(legend.position="none",text=element_text(size=20)) 
      if(length(plotting_details$types) < 200) {
        p <- p+geom_dotplot(binaxis='y',stackdir='center',dotsize=0.5,mapping=aes(fill=types,color=types))
      }
      p
      ggsave(file.path(image_results_directory,paste('boxplot_',plotting_details$name,plotting_name_stem,'.svg',sep="")),width=15,height=10)
    } 
  }

  if(opt$svm!=FALSE) {
    plotting_details <- plotting_data$image_svm
    classes <- levels(summaries_results$class_names)
    classification_status <- vector("character",length=length(plotting_details$data))
    classification_name <- vector("character",length=length(plotting_details$data))
    for(i in 1:length(plotting_details$data)) {     
      if(plotting_details$data[i] < 0) { 
        classification_name[i] <- classes[1]
        if(plotting_details$named_types[i] == classes[1]) {
          classification_status[i] <- "Correct"  
        } else { 
          classification_status[i] <- "Incorrect"
        }
      } else { 
        classification_name[i] <- classes[2]
        if(plotting_details$named_types[i]==classes[2]) {
         classification_status[i] <- "Correct"  
        } else {
          classification_status[i] <- "Incorrect"
        }
      }
    }
    boxplot_image_data <- data.frame(features=plotting_details$data,types=plotting_details$named_types,testtypes=factor(classification_name),labelstatus=classification_status,replicate=plotting_data$image_svm$which_run)  
    # Confusion matrices 
    boxplot_data_for_confusion <- boxplot_image_data

    boxplot_data_for_confusion[,"testtypes"] <- factor(boxplot_data_for_confusion[,"testtypes"],levels=rev(levels(boxplot_data_for_confusion[,"testtypes"])))
    actual_labels <- levels(boxplot_data_for_confusion[,"types"])
    levels(boxplot_data_for_confusion[,"types"]) <- c(paste("Actual ",actual_labels[1]),paste("Actual ",actual_labels[2]))
    predicted_labels <- levels(boxplot_data_for_confusion[,"testtypes"])
    levels(boxplot_data_for_confusion[,"testtypes"]) <- c(paste("Predicted ",predicted_labels[1]),paste("Predicted ",predicted_labels[2]))

    base_confusion_matrix <- table(boxplot_data_for_confusion[,"testtypes"],boxplot_data_for_confusion[,"types"])
    if(actual_labels[1]!=predicted_labels[1]) {
      base_confusion_matrix <- base_confusion_matrix[,c(2,1)]
      accuracy <- 100*sum(diag(base_confusion_matrix))/sum(base_confusion_matrix)
    } else { 
      accuracy <- 100*sum(diag(base_confusion_matrix))/sum(base_confusion_matrix)  
    }
    confusion_matrix <- as.data.frame(base_confusion_matrix[,c(2,1)])

    p <- ggplot(data = confusion_matrix,
           mapping = aes(x = Var1,
                         y = Var2)) +
      geom_tile(aes(fill = Freq)) +
      geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1,size=10) +
      scale_fill_gradient(low = "blue",
                          high = "red") +
      xlab(paste("Testing accuracy: ",round(accuracy,0),"%")) +
      ylab("") +
      theme(legend.position="none",panel.background = element_blank(),text=element_text(size=16,face="bold"))
    p
    ggsave(file.path(image_results_directory,paste(plotting_name_stem,"_confusion_matrix_.svg",sep="")),width=7,height=3.5)    
    image_file_names <- training_results$image_file_names
    save(plotting_data,boxplot_image_data,image_file_names,testing_order_image_names,file=file.path(image_results_directory,paste(plotting_name_stem,".RData",sep="")),version=2)
  } else { 
    save(plotting_data,file=file.path(image_results_directory,paste(plotting_name_stem,".RData",sep="")),version=2)
  }
}
  
if(opt$fullsave) {
  data_file_name <- paste(training_results$data_name_stem,".RData",sep="")
  save.image(file=file.path(data_results_directory,data_file_name),version=2)
}


