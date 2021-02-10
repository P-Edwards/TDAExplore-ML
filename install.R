
# Base packages
base_packages <- c("devtools","optparse", "LiblineaR","colorspace","ggplot2","gridExtra","ggplotify","svglite")


# Dependency packages for Windows users
windows_base_packages <- c("Rcpp","methods","foreach","doParallel","TDA","iterators","SparseM","OpenImageR","rdist")

if(.Platform$OS.type == "windows") { 
	packages_to_install <- c(base_packages,windows_base_packages)
} else { 
	cat("Warning: Several required packages rely on shared libraries you may have to install via your package manager. See the README for some details.\n\n")
	packages_to_install <- base_packages
}


new.packages <- packages_to_install[!(packages_to_install %in% installed.packages()[,"Package"])]
if(length(new.packages)) { install.packages(new.packages) }

packages_to_install <- c("TDAExplore")
new.packages <- packages_to_install[!(packages_to_install %in% installed.packages()[,"Package"])]
if(length(new.packages)) { 
	if(.Platform$OS.type== "windows") { 
		install.packages("https://dl.dropboxusercontent.com/s/4h3d8ghvhkjmthn/TDAExplore.zip",repos=NULL)
	} else { 
		devtools::install_github("P-Edwards/TDAExplore")
	}
}

# Can put copies of shell scripts somewhere in the path, if desired
copy_to_path <- utils::askYesNo("Would you like to install copies of the TDAExplore scripts into a folder accessible from the command line? Defaults require you install as admin.")
cat(" \n \n")

if(copy_to_path) { 
	if(.Platform$OS.type=="windows") { 
		current_path <- file.path(Sys.getenv("ProgramFiles"),"TDAExplore")
		ml_executable <- "ml-tda.bat"
		convolve_executable <- "convolve-tda.bat"
	} else { 
		current_path <- "/usr/local/bin"
		ml_executable <- "ml-tda"
		convolve_executable <- "convolve-tda"
	}
	new_path <- readline(prompt=paste("Install script copies to: (leave blank for default ",current_path,") "))
	if(length(new_path)>1) { 
		current_path <- new_path
	} else { 
		# Need to add TDAExplore to PATH on Windows
		if(.Platform$OS.type=="windows") { 
			shell(paste('setx PATH "',normalizePath(current_path),'"',sep=""))
		}
	}
	if(!dir.exists(normalizePath(current_path))) { dir.create(normalizePath(current_path),showWarnings=FALSE) }
	
	# Make sure the scripts are executable and copy everything
	if(!(.Platform$OS.type=="windows")) { 
		Sys.chmod(ml_executable)
		Sys.chmod(convolve_executable)
	}

	file.copy(c(ml_executable,convolve_executable,"landscapes_and_ML.R","convolve_images.R"),current_path,overwrite=TRUE)
	print("Done.")	
}