
# Adaptable RtSNE

    # Description
      # Running bh-tSNE on data from FCS or CSV files, with additional functionality
      # Read data from FCS files, run tSNE, write back to FCS file with new tSNE values
      # Prints the 'verbose' from the tSNE analysis as well
      # Most parameters that control the tSNE algorithm are modifiable in the Rtsne function

    # Details
      # Author: Thomas Myles Ashhurst
      # Published: 2017-09-18
      # Contact: tomashhurst@gmail.com
      # Github: https://github.com/sydneycytometry/Adaptable-RtSNE
      # Website: www.sydneycytometry.org.au

    # Acknowledgements
      # tsne algorithm developed by Lauren van der Maaten https://github.com/lvdmaaten/bhtsne
      # van der Maaten, L.J.P., Hinton, G.E. (2008). Visualizing High-Dimensional Data Using t-SNE. Journal of Machine Learning Research 9(Nov):2579-2605.
      # van der Maaten, L.J.P. (2014). Accelerating t-SNE using Tree-Based Algorithms. Journal of Machine Learning Research 15(Oct):3221-3245.
      # tsne licence: https://github.com/lvdmaaten/bhtsne/blob/master/LICENSE.txt

      # Script structure adapted from https://github.com/lmweber/Rtsne-example
      # Rtsne function by Jesse Krijthe
      # Version 0.12 and 0.13 of Rtsne compiled after discussion with Jesse Krijthe, regarding adding additional input arguments into Rtsne function
      # FCS file reading and writing adapted from https://gist.github.com/yannabraham/c1f9de9b23fb94105ca5

    # Version
      # R version 3.2.3 (2015-12-10) -- "Wooden Christmas-Tree"
      # flowCore version    1.42.2
      # plyr version        1.8.4
      # Biobase version     2.36.2
      # Rtsne version       0.13        # important, v 0.12 and 0.13 use different value systems for stop_lying_iter and mom_switch_iter
      # flowViz version     1.40.0

    # Summary
      # 1. Load packages and set working directory
      # 2. Read data from FCS file or CSV file
      # 3. Select parameters to use in tSNE run
      # 4. Adjust sample scaling and size (if necessary)
      # 5. Run tSNE (also saving verbose)
      # 6. Save tSNE output as .csv
      # 7. Add tSNE parameters to data
      # 8. Write data (with tSNE parameters) to new .csv and .fcs file

      # Creating coloured plots is in a different script: github.com/sydneycytometry/tSNEplots


############################## BEGIN USER INPUT ##############################


### 1. Load packages and set working directory

    ## Install packages (if not already installed)
    if(!require('flowCore')) {install.packages('flowCore')}
    if(!require('plyr')) {install.packages('plyr')}
    if(!require("Biobase")) {install.packages("Biobase")}
    if(!require("Rtsne")) {install.packages("Rtsne")}
    if(!require("flowViz")) {
      source("https://bioconductor.org/biocLite.R") 
      biocLite("flowViz")}

    ## Load packages
    library('flowCore')
    library('plyr')
    library('Biobase')
    library('Rtsne')
    library('flowViz')

    ## Set working directory and assign as 'PrimaryDirectory'
    setwd("/Users/thomasashhurst/Desktop/SampleFolder/Adaptable-RtSNE/") # choose your working directory here
    PrimaryDirectory <- getwd()
    PrimaryDirectory

    
### 2. Chose FCS file for tSNE, read into 'data'

    ## Create name for output files
    Data_name <- "Sample_Data" # name the dataset something meaningful
    Data_name
    
    ## Show the name of .fcs or .csv files in the working directory and load data from those files
       
        ## OPTION A: FCS
        FileNames <- list.files(path=PrimaryDirectory, pattern = ".fcs")
        FileNames
        
        fcsfilename <- "sample_data.fcs" # Enter name of FCS file in between the ""

        data <- exprs(read.FCS(fcsfilename, transformation = FALSE)) 
        dim(data) # check the total number of events (first entry), and number of parameters (second entry)
        
        data_parameters <- parameters(read.FCS(fcsfilename)) 
        data_parameters
        
        
        ## OR ##

          
        ## OPTION B: CSV
        FileNames <- list.files(path=PrimaryDirectory, pattern = ".csv")
        FileNames
        
        csvfilename <- "sample_data.csv"
        
        data <- read.csv(csvfilename)
        dim(data)

    ## Review data
    head(data) # show data with headings
    unname(colnames(data))  # show reporter and marker names
    ColumnNames <- unname(colnames(data)) # assign reporter and marker names (column names) to 'ColumnNames'
    ColumnNames # Use ColumnNames (assigned above) to review the column names

    
### 3. Adjust sample scaling and size

    ## Remove duplicates (optional, may cause problems for tSNE if ignored)
      dim(data) # note the second entry, which is number of parameters
      data <- data[!duplicated(data), ]  # remove rows containing duplicate values within rounding
      dim(data) # check to see that the number of parameters has reduced
    
    
    ## IF REQUIRED, Arcsinh transformation (only if required, data may have already been transformed)
      data # note the kinds of values in the data set (i.e. are the between 1 and 10, 100 and 1000, etc)
      asinh_scale <- 5 # for CyTOF choose between 5 and 15, for flow you will need > 200.
      data <- asinh(data / asinh_scale)  # transforms all columns! including event number etc
      data # check the data has been converted 

    
    ## IF REQUIRED, subsampling (only if required)
      nsub <- 1000 # set the desired number of cells for subsampling
      set.seed(123)  # set random seed
      data <- data[sample(1:nrow(data), nsub), ]
      dim(data)

    # Note, the final tSNE parameters are added to THIS data stage (after transformation and subsampling)
    # The parameter selection in step 4 specify which parameters are used to generate the tSNE parameters
    

### 4. Select parameters for tSNE run

    ## Select the markers to use in calculation of t-SNE parameters, by including below which columns in 'data' to include
    ColumnNames
    
    ColumnNames_for_tSNE <- unname(colnames(data))[c(2,3,4,5)] # e.g. [c(11, 23, 10)] to include the markers corresponding to the column numbers 11, 23, 10
    ColumnNames_for_tSNE  # check that the column names that appear are the ones you want to analyse
    
    dim(data) # Review dimensionality of 'data' -- N events and N parameters
    data_specific <- data[, ColumnNames_for_tSNE] # Prepare data for Rtsne --  select columns to use
    dim(data_specific)     # Review dimensionality of 'data' -- N events and N parameters

    
############################## END USER INPUT (unless modifying tSNE defaults) ##############################


### 5. Run tSNE ('data' at this point is scaled, transformed) -- RUN ALL OF STEP 5
    
    ## Create 'Output_info' folder and set as working directory
    setwd(PrimaryDirectory)
    
    dir.create("Output", showWarnings = FALSE)
    setwd("Output")
    OutputDirectory <- getwd()
    
    dir.create("Output_info", showWarnings = FALSE)
    setwd("Output_info")
    getwd()

    ## Export subsampled data in CSV format (record of what was run)
    subsampled_data_name <- paste0("Output_subsampled_data_", Data_name)
    subsampled_data <- paste(subsampled_data_name, ".csv", sep = "")
    write.csv(x = data, file = subsampled_data)

    ## Run tSNE
        # Will send the verbose text to a .txt file
        verbose_name <- paste0("Output_verbose_", Data_name, ".txt")
        sink(file = verbose_name, append=TRUE, split=FALSE, type = c("output", "message"))

        # Run tSNE algorithm
        set.seed(42)                                    # default = 42 -- sets seed for reproducibility, delete or comment out for random tSNE
        tsne_out <- Rtsne(as.matrix(data_specific),     # dataset (assigned above) read as a matrix
                          dims = 2,                     # default = 2 -- only input 2 or 3 (2 = plot in 2D, 3 = plot in 3D)
                          initial_dims = 50,            # default = 50 -- number of dimensions retained in initial PCA step
                          perplexity = 30,              # default = 30
                          theta = 0.5,                  # default = 0.5 -- use 0.5 for Barnes-Hut tSNE, 0.0 for exact tSNE (takes longer)
                          check_duplicates = TRUE,      # default = TRUE, can set to FALSE if problematic
                          pca = TRUE,                   # default = TRUE -- performs a PCA prior to actual tSNE run
                          max_iter = 1000,              # default = 1000 -- default total iterations
                          verbose = TRUE,               # default = FALSE -- best to set as TRUE to receive feedback
                          is_distance = FALSE,          # default = FALSE -- experimental, using X as a distance matrix
                          Y_init = NULL,                # default = NULL -- recommend to use NULL for random initialisation
                          stop_lying_iter = 250,        # default = 250 -- number of iterations of early exaggeration
                          mom_switch_iter = 250,        # default = 250 -- number of iterations before increased momentum of spread
                          momentum = 0.5,               # default = 0.5 -- initial momentum of spread
                          final_momentum = 0.8,         # default = 0.8 -- momentum of spread at 'final_momentum'
                          eta = 200,                    # default = 200 -- learning rate
                          exaggeration_factor = 12.0    # default = 12.0 -- factor used during early exaggeration
                          )

        # turn sink off
        sink()

        # if tSNE detected duplicates, the following errors should occur: 
            # "Error in Rtsne.default(as.matrix(data_specific), dims = 2, initial_dims = 50,  : Remove duplicates before running TSNE."
            # "Error in inherits(.data, "split") : object 'tsne_out' not found"
            
            
### 6. Save tSNE output information at .csv

    ## Save tsne_out output (i.e., output settings ad raw results)
    tsne_out.df <- ldply(tsne_out, data.frame)
    output_name <- paste0("Output_parameters_", Data_name, ".csv")
    write.csv(x = tsne_out.df, file = output_name) # pretty blood good -- doesn't give row number for Y, costs, or itercosts -- but easy to figure out


### 7. Add tSNE parameters to data
    
    ## Create 'Output_info' folder and set as working directory
    setwd(OutputDirectory)
    dir.create("Output_data", showWarnings = FALSE)
    setwd("Output_data")
    getwd()
    
    ## Add Y (tSNE coordinate) values to starting data
    tsne_out_Y <- tsne_out$Y
    tsne_out_Y

    data <- cbind(data, tSNE1 = tsne_out_Y[,1])
    data <- cbind(data, tSNE2 = tsne_out_Y[,2])

  
### 8. Write data (with tSNE parameters) to .csv and .fcs
    
    ## Save data (with new tSNE parameters) as CSV
    new_file_name_csv <- paste0("Output_data_", Data_name, "_with_tSNE", ".csv")
    write.csv(x = data, file = new_file_name_csv)
    
    ## Check data and data column names
    head(data)
    dimnames(data)[[2]]
    
    ## Create FCS file metadata - column names with descriptions
    metadata <- data.frame(name=dimnames(data)[[2]],
                       desc=paste('this is column',dimnames(data)[[2]],'from your CSV')
                       )
    
    ## Double check metadata
    metadata
    
    ## Create FCS file metadata - ranges, min, and max settings
    metadata$range <- apply(apply(data,2,range),2,diff)
    metadata$minRange <- apply(data,2,min)
    metadata$maxRange <- apply(data,2,max)
    
    metadata$range
    metadata$minRange
    metadata$maxRange
    
    ## Create flowframe with tSNE data
    data.ff <- new("flowFrame",
              exprs=data,
              parameters=AnnotatedDataFrame(metadata)
              )
    
    ## Check flow frame
    data.ff
    
    ## Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
    new_file_name_fcs <- paste0("Output_data_", Data_name, "_with_tSNE", ".fcs")
    write.FCS(data.ff, new_file_name_fcs)
    
    ## Move back to PrimaryDirectory
    setwd(PrimaryDirectory)
    getwd()
    
    # note -- in the final data output, all parameters are included, but only the subsampled and/or transformed cells
    
