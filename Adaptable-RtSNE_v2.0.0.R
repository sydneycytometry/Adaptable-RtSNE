# Adaptable RtSNE v2.0.0

    # Description
      # Running bh-tSNE on data from FCS or CSV files, with additional functionality including phenograph
      # Read data from FCS files, run tSNE, write back to FCS file with new tSNE values
      # Prints the 'verbose' from the tSNE analysis as well
      # Most parameters that control the tSNE algorithm are modifiable in the Rtsne function

    # Details
      # Author: Thomas Myles Ashhurst
      # Published: 2017-10-11
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
      
      # 1 INSTALL AND LOAD PACKAGES
        # 1.1 - Install packages
        # 1.2 - Load packages
        # 1.3 - Set working directory

      # 2 USER INPUT -- DATA PREPARATION
        # 2.1 - Define options
        # 2.2 - Read files (FCS or CSV) into data
        # 2.3 - Create parameter names
        # 2.4 - Remove duplicates (optional, but recommended)
        # 2.5 - Modifications (if required, skip if not required)
        # 2.6 - Select markers for the calculation of t-SNE parameters

      # 3 END USER INPUT -- RUN ENTIRE SCRIPT
        # 3.1 - Create output directories
        # 3.2 - Save a copy of 'data' (original data) and 'data specific' (data used to calculate tSNE)
        # 3.3 - Run tSNE (if defined in 2.1)
        # 3.4 - Run Phenograph (if defined in 2.1)
        # 3.5 - Write data (with tSNE and phenograph parameters) to .csv
        # 3.6 - Write data (with tSNE and phenograph parameters) to .fcs

        # 3.7 - Create coloured tSNE plots (separate script: github.com/sydneycytometry/tSNEplots)


###################################################### 1. INSTALL AND-OR LOAD PACKAGES ###################################################### 

    ### 1.1 - Install packages (if not already installed)
        if(!require('flowCore')) {install.packages('flowCore')}
        if(!require('plyr')) {install.packages('plyr')}
        if(!require('Biobase')) {install.packages("Biobase")}
        if(!require('Rtsne')) {install.packages("Rtsne")}
        if(!require('ggplot2')) {install.packages('ggplot2')}
        if(!require('devtools')){install.packages("devtools")}
        if(!require('Rphenograph')){devtools::install_github("JinmiaoChenLab/Rphenograph")}
        if(!require('flowViz')) {source("https://bioconductor.org/biocLite.R") 
          biocLite('flowViz')}

    ### 1.2 Load packages       
        library('flowCore')
        library('plyr')
        library('Biobase')
        library('Rtsne')
        library('ggplot2')
        library('devtools')
        library("Rphenograph")
        library('flowViz')    


    ### 1.3 - Set working directory and assign as 'PrimaryDirectory'

        setwd("/Volumes/GoogleDrive/My Drive/2017 Ashhurst work/1. ASHHURST/01. R&D - Technology/[Script projects]/New/Phenograph/Sample data/")

        setwd("/Users/Tom/Desktop/tSNEplots-master 2/Input/")           # Set your working directory here (e.g. "/Users/Tom/Desktop/") -- press tab when selected after the '/' to see options
        getwd()                                                         # Check your working directory has changed correctly
        PrimaryDirectory <- getwd()                                     # Assign the working directory as 'PrimaryDirectory'
        PrimaryDirectory

    ### Alternative working directory method, with rstudio api
        
        #dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        #setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        #getwd()
        #PrimaryDirectory <- getwd()
        #PrimaryDirectory
        

###################################################### 2. BEGIN USER INPUT ###################################################### 

    ### 2.1 - Define options
        
        ## Use to list the .fcs or .csv files in the working directory
        list.files(path=PrimaryDirectory, pattern = ".fcs")     # list of FCS files
        list.files(path=PrimaryDirectory, pattern = ".csv")     # list of CSV files
        
        ## Specify options
        Data_name         <- "Sample_Data"      # Name the dataset something meaningful (e.g. "CNS_sample", "Human_PMBC")
        File_Type         <- ".csv"             # Use ".fcs" OR ".csv", in lower case
        Run_Phenograph    <- 1                  # Enter a 1 to run Phenograph, or 0 to skip
        Run_tSNE          <- 1                  # Enter a 1 to run tSNE, or 0 to skip 
        
        # Important, the only FCS/CSV file in the directory should be the one desired for analysis. If more than one are found, only the first file will be used
        
    ### 2.2 - Read files into 'data' (run all at once, no modifications required) -- if both are there, CSV will take priority
        
        if (File_Type == ".fcs"){ # if files are .fcs, will read into data
          FileNames <- list.files(path=PrimaryDirectory, pattern = ".fcs")
          FileNames
          fcsfilename <- FileNames[[1]] # Enter name of FCS file in between the ""
          data <- exprs(read.FCS(fcsfilename, transformation = FALSE)) 
          dim(data) # check the total number of events (first entry), and number of parameters (second entry)
          data_parameters <- parameters(read.FCS(fcsfilename)) 
          data_parameters}
        
        if (File_Type == ".csv"){ # if files are .csv, will read into data
          FileNames <- list.files(path=PrimaryDirectory, pattern = ".csv")
          FileNames
          FileNames[[1]]
          csvfilename <- FileNames[[1]]
          data <- read.csv(csvfilename)
          dim(data)}

    ### 2.3 - Create set of column (parameter) names
        
        head(data)                            # show data with headings
        ColumnNames <- unname(colnames(data)) # assign reporter and marker names (column names) to 'ColumnNames'
        ColumnNames                           # Use ColumnNames (assigned above) to review the column names

    ### 2.4 - Remove duplicates (optional, may cause problems for tSNE if ignored)
        
        dim(data) # note the second entry, which is number of parameters
        data <- data[!duplicated(data), ]  # remove rows containing duplicate values within rounding
        dim(data) # check to see if the number of parameters has reduced (no change means no duplicates were present)
            
    ### 2.5 - Modifications (IF REQUIRED, SKIP IF NOT REQUIRED)        
        
            ## IF REQUIRED, Arcsinh transformation (only if required, data may have already been transformed)
            data                                # note the kinds of values in the data set (i.e. are the between 1 and 10, 100 and 1000, etc)
            asinh_scale <- 5                    # for CyTOF choose between 5 and 15, for flow you will need > 200.
            data <- asinh(data / asinh_scale)   # transforms all columns! including event number etc
            data                                # check the data has been converted 
            
            
            ## IF REQUIRED, subsampling (only if required)
            dim(data)
            nsub <- 10000                               # set the desired number of cells for subsampling
            set.seed(123)                               # set random seed
            data <- data[sample(1:nrow(data), nsub), ]
            dim(data)
            
    ### 2.6 - Select the markers to use in calculation of t-SNE parameters, by including below which columns in 'data' to include
        as.matrix(ColumnNames) # view the column 'number' for each parameter
        
        ColumnNames_for_tSNE <- unname(colnames(data))[c(2,3,4,5)] # e.g. [c(11, 23, 10)] to include the markers corresponding to the column numbers 11, 23, 10
        ColumnNames_for_tSNE  # check that the column names that appear are the ones you want to analyse
        
        dim(data) # Review dimensionality of 'data' -- N events and N parameters
        data_specific <- data[, ColumnNames_for_tSNE] # Prepare data for Rtsne --  select columns to use
        dim(data_specific) # Review dimensionality of 'data' -- N events and N parameters
        
        
        # Note, the final tSNE and Phenograph parameters are added 'data' (after transformation and subsampling), not 'data_specific' stage
        # data_specific is used to specify which parametes are used to run tSNE and Phenograph
      

###################################################### 3. END USER INPUT (unless modifying tSNE defaults) ###################################################### 
        
    ### Run the rest of the script all at once

    ### 3.1 - Create output directories
        
        setwd(PrimaryDirectory)
        
        dir.create("Output", showWarnings = FALSE)
        setwd("Output")
        OutputDirectory <- getwd()
        
        setwd(OutputDirectory)
        dir.create("Output_info", showWarnings = FALSE)
        setwd(OutputDirectory)
        dir.create("Output_data", showWarnings = FALSE)

    ### 3.2 - Save a copy of 'data' and 'data specific'        
        
        setwd(OutputDirectory)
        setwd("Output_info")
        
        ## Export subsampled data in CSV format (record of what was run)
        subsampled_data_name <- paste0("Subsampled_data_", Data_name)
        subsampled_data <- paste(subsampled_data_name, ".csv", sep = "")
        write.csv(x = data, file = subsampled_data)
        
        ## Export a list of parameters used to run tSNE and Phenograph
        tSNEparametersname <- paste0("Columns_used_for_tSNE_and_phenograph_", Data_name)
        tSNEparameters <- paste(tSNEparametersname, ".csv", sep = "")
        write.csv(x = data_specific, file = tSNEparameters)      
        
    ### 3.3 - Run tSNE ('data' at this point is scaled, transformed) -- RUN ALL OF STEP 5
        
        if (Run_tSNE == 1) {
        
            ## Set working directory to "Output_info"
            setwd(OutputDirectory)
            setwd("Output_info")

            ## Turn sink on (send the verbose text from the tSNE algorithm progress updates to a .txt file)
            verbose_name <- paste0("tSNE_verbose_", Data_name, ".txt")
            sink(file = verbose_name, append=TRUE, split=FALSE, type = c("output", "message"))
    
            ## Run tSNE algorithm
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
    
            ## turn sink off
            sink()
            
            ## if tSNE detected duplicates, the following errors should occur: 
                # "Error in Rtsne.default(as.matrix(data_specific), dims = 2, initial_dims = 50,  : Remove duplicates before running TSNE."
                # "Error in inherits(.data, "split") : object 'tsne_out' not found"
            
            ## Save tsne_out output info as .csv (i.e., output settings as raw results)
            tsne_out.df <- ldply(tsne_out, data.frame)
            output_name <- paste0("tSNE_parameter_values_", Data_name, ".csv")
            write.csv(x = tsne_out.df, file = output_name) # pretty blood good -- doesn't give row number for Y, costs, or itercosts -- but easy to figure out
            
            ## Add Y (tSNE coordinate) values to starting data
            tsne_out_Y <- tsne_out$Y
            tsne_out_Y
            
            data <- cbind(data, tSNE1 = tsne_out_Y[,1])
            data <- cbind(data, tSNE2 = tsne_out_Y[,2])
            
            }
            
    ### 3.4 - Run Phenograph

        if (Run_Phenograph == 1) {
          
          setwd(OutputDirectory)
          setwd("Output_info")
          
          # set sink (send the verbose text from the tSNE algorithm progress updates to a .txt file)
          verbose_ph_name <- paste0("Phenograph_verbose_", Data_name, ".txt")
          sink(file = verbose_ph_name, append=TRUE, split=FALSE, type = c("output", "message"))
          
          # Run Phenograph
          data_ph_matrix <- as.matrix(data_specific)
          Rphenograph_out <- Rphenograph(data_ph_matrix, k = 45)
          Rphenograph_out
          
          modularity(Rphenograph_out[[2]])
          membership(Rphenograph_out[[2]])
          
          # Save phenograph cluster number results
          output_name_ph <- paste0("Phenograph_cluster_values_", Data_name, ".csv")
          write.csv(x = factor(membership(Rphenograph_out[[2]])), file = output_name_ph)
          
          # Add phenograph cluster numbers to 'data'
          as.numeric(as.character(factor(membership(Rphenograph_out[[2]])))) # converts factor to integer
          data$phenograph_cluster <- as.numeric(as.character(factor(membership(Rphenograph_out[[2]])))) # writes integer version to data
          
          # end sink
          sink()
          sink() # just in case
        }
        
        # ggplot(data, aes(x=tSNE1, y=tSNE2, col=phenograph_cluster)) + geom_point(size = 3)+theme_bw()
        # ggplot(data, aes(x=X149Sm_CD19, y=X152Gd_CD3e, col=phenograph_cluster)) + geom_point(size = 3)+theme_bw()
        # shape=phenograph_cluster

    ### 3.5 - Write data (with tSNE and phenograph parameters) to .csv
    
        ## Set working directory
        setwd(OutputDirectory)
        setwd("Output_data")
        
        ## Save data (with new tSNE parameters) as CSV
        new_file_name_csv <- paste0("Output_", Data_name, "_with_tSNE", ".csv")
        write.csv(x = data, file = new_file_name_csv)

   ### 3.6 - Write data (with tSNE and phenograph parameters) to .fcs
       
        setwd(OutputDirectory)
        setwd("Output_data")
        
        ## Check data and data column names
        head(data)
        dimnames(data)[[2]]
        
        ## Create FCS file metadata - column names with descriptions
        metadata <- data.frame(name=dimnames(data)[[2]],
                           desc=paste('column',dimnames(data)[[2]],'from dataset')
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
                  exprs=as.matrix(data), # in order to create a flow frame, data needs to be read as matrix
                  parameters=AnnotatedDataFrame(metadata)
                  )
        
        ## Check flow frame
        data.ff
        
        ## Save flowframe as .fcs file -- save data (with new tSNE parameters) as FCS
        new_file_name_fcs <- paste0("Output_ ", Data_name, "_with_tSNE", ".fcs")
        write.FCS(data.ff, new_file_name_fcs)
        
        ## Move back to PrimaryDirectory
        setwd(PrimaryDirectory)
        getwd()
        
        # note -- in the final data output, all parameters are included, but only the subsampled and/or transformed cells

    ### 3.7 - Create coloured tSNE plots
        
        ## Uses a separate script
        ## github.com/sydneycytometry/tSNEplots
