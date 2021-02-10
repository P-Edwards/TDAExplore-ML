TDAExplore-ML
-------------

Copyright (C) 2021 [Parker
Edwards](https://sites.nd.edu/parker-edwards)

Description
-----------

Command line R applications for machine learning using the R package [TDAExplore](https://github.com/P-Edwards/TDAExplore).

Version 1.0.0
-------------

External requirements
---------------------
[R 4.0](https://www.r-project.org/)


Installation
---------------------------
Install TDAExplore-ML using the included script install.R. From the command line, make sure that this package's root directory is the working directory (using `getwd()`  and `setwd()` in R). Then run:

``` sh
	R
	> source("install.R")    
```
Similarly, from RStudio, check that the working directory is this package's root directory. Then run:

```R
source("install.R")
```

On Windows you will need to add `R.exe` and `Rscript.exe` to your PATH to run the batch scripts properly. See [here](https://cran.r-project.org/bin/windows/base/rw-FAQ.html#Rcmd-is-not-found-in-my-PATH_0021) for instructions.

On Linux several of the R package dependencies of TDAExplore-ML require that your system has certain C/C++ libraries installed. The links below are to appropriate Ubuntu packages to highlight the correct library. This list is not exhaustive: recursive package requirements may requirem more libraries. If library installations throw such an error, search through the logs for the right library, install it, then rerun the TDAExplore-ML install script.

* [libtiff-dev](https://packages.ubuntu.com/search?keywords=libtiff-dev)
* [libgmp-dev](https://packages.ubuntu.com/search?keywords=libgmp-dev)

Usage
------
The package includes two command line applications: `ml-tda` and `convolve-tda`. On Windows these are batch scripts. For both, pass with the `--help` flag to see a listing of command line flags. The applications require a `.csv` parameters file with path indicated by `.csv`. A template is provided that you can copy, rename, edit, etc. 


```sh
ml-tda --parameters parameters_template.csv --svm TRUE --cores 6 
```

If you do not want to relocate the scripts or put them on your PATH, you can always call them directly:

```sh 
TDAExplore-ML/ml-convolve --parameters parameters_template.csv --svm TRUE --cores 6
```

Example data
------------
Images for the examples are distributed separately. You can download the data [here](https://drive.google.com/drive/folders/1LJSaqZTr9eVa8DcO65XQyPKtIQiy0w_D?usp=sharing) and place the folder in this project's root directory to execute the examples.


License
-------
TDAExplore is licensed under GPLv3. 
