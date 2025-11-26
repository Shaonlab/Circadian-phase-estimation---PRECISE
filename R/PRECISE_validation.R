#owner:ANJOOM NIKHAT
#date:24-11-2025
#lab:Dr. Shaon Chakrabarti's lab
#program to run PRECISE to infer the phase of different samples (single or other wise)
#this code specifically validates the performance of precise so we know the actual phases of the test samples as well
#######################################################################################################################
start_time <- Sys.time() #to check the time required for implementation of PRECISE
#load the required libraries

#installing the OdeGP package latest version
#source("~/Documents/ANJOOM/time_telling_project/final_submission_documents/github_page_files/ODeGP-main/oscSetup.R")


#install the local version of OdeGP package if not already done
library(ODeGP)  #package to learn the oscillatory patterns in the training gene expression dataset
library(DEoptim)#package used of hyper parameter optimisation during ODeGP fit
library(tidyverse)
library(ggpubr) #package used for plotting
library(ggpmisc) #package used to detect the local maximas in the likelihood versus phase curve
library(circular)  #package to calculate the circular mean to get the average phase estimate mentione in the paper

gene_names <- c("bmal1","tef","nr1d1","nr1d2") #the genes that you want to train PRECISE one and also used for prediction

#TRAINING DATA INPUT
#input file that has the training phases information. The training phases are determined using any one gene with the help of wavelet transform
train_timepoint_phase_data <- read.csv("~/Documents/ANJOOM/time_telling_project/final_submission_documents/github_page_files/github_page_upload_files/training_phase_nr1d1.csv")
train_phase <- train_timepoint_phase_data$phase

#search limits for the phase outputs. We search for which training phase fits the data the best
min_train_phase <- min(train_phase)
max_train_phase <- max(train_phase)

#training data files (list of files each of which represents one time point)

#set the working directory : folder that contains the training files
setwd("~/Documents/ANJOOM/time_telling_project/final_submission_documents/github_page_files/training_data/")
train_files <- list.files()

Raw_train <- list() #we will now import each of the files corresponding to the different time points in the training set
for(i in 1:length(train_files)){
  print(train_files[[i]]) #just to check if the ordering of the time points are correct or not
  temp <- read.csv(train_files[[i]])
  temp <- temp[,gene_names] #only including the columns that have the gene mRNA counts
  
  #the name of the file when including in the list of Raw_training_files is the phase value assigned to that timepoint
  Raw_train[[as.character(train_phase)[[i]]]] <- temp 
}


#########################################################################################################
#TEST DATA INPUT 
#input file that has the test phases information. The test phases are determined using the same gene used in training data
test_timepoint_phase_data <- read.csv("~/Documents/ANJOOM/time_telling_project/final_submission_documents/github_page_files/github_page_upload_files/test_phase_nr1d1.csv")
test_phase <- test_timepoint_phase_data$phase
test_time <- as.character(test_timepoint_phase_data$time)

#test data files (list of files each of which represents one time point)

#set the working directory : folder that contains the test files
setwd("~/Documents/ANJOOM/time_telling_project/final_submission_documents/github_page_files/test_data/")
test_files <- list.files()

Raw_test <- list() #we will now import each of the files corresponding to the different time points in the testing set
for(i in 1:length(test_files)){
  print(test_files[[i]]) #just to check if the ordering of the time points are correct or not
  temp <- read.csv(test_files[[i]])
  temp <- temp[,gene_names] #only including the columns that have the gene mRNA counts
  
  #the name of the file when including in the list of Raw_training_files is the phase value assigned to that timepoint
  Raw_test[[as.character(test_phase)[[i]]]] <- temp 
}

#########################################################################################################
#VARIOUS FUNCTIONS THAT WILL REQUIRED FOR PHASE INFERENCE

#FUNCTION 1: the function to bootstrap the data from the different time points
bootstrap <- function(data,ntimes,nsamples,gene_names){
  #data : list containing files corresponding to the different time points
  #ntimes : number of final bootstapped samples required (number of iterations)
  #nsamples : number of samples to be averaged over in each iteration
  #gene_names : the genes that the user wants to utilise (oscillatory genes)
  
  Bootstrap_data <- list()
  for(num in 1:length(data)){
    df <- data[[num]]
    results_list <- vector("list", ntimes)  #list to carry the means of the sampled rows for all the genes per iteration
    
    #perform the sampling and the calculation for each iteration
    for(i in 1:ntimes){
      sampled_rows <- sample_n(df,nsamples)
      
      #when considering multiple oscillatory genes
      if(length(gene_names) > 1){mean_values <- colMeans(sampled_rows[,gene_names])}
      
      #when considering only one oscillatory gene
      if(length(gene_names) == 1){mean_values  <- mean(sampled_rows[[gene_names]])}
      
      results_list[[i]] <- mean_values
    }
    
    # Convert results to a data frame
    results_df <- as.data.frame(do.call(rbind, results_list))
    
    # Add column names to the results data frame
    colnames(results_df) <- gene_names
    
    Time <- replicate(nrow(results_df),names(data)[[num]]) 
    results_df <- cbind(Time, results_df)
    
    Bootstrap_data[[names(data)[[num]]]] <- results_df
  }
  return(Bootstrap_data)
}


#FUNCTION 2: function to perform the normalize the data across different time points
data_normalize_func <- function(data,gene_names){
  #data : bootstrapped data
  #gene_names : the genes that the user wants to utilise (oscillatory genes)
  
  #list to contain the mean and sd values of the different genes across the different time points
  #the mean and sd values will be used to normalise the data
  gene_normal_data <- list()
  for (gene in gene_names){
    gene_nfactor <- list()
    gene_val_time <- c() #list that will contain the mean expression value of gene across timepoints
    for(num in 1:length(data)){ 
      gene_val_time <- c(gene_val_time,mean(data[[num]][[gene]]))
    }
    gene_nfactor[["mean"]] <- mean(gene_val_time) #mean of a gene across time
    gene_nfactor[["sd"]] <- sd(gene_val_time) #sd of a gene across time
    
    gene_normal_data[[gene]]<-gene_nfactor
  }
  
  #normalising the spot count for the different genes for each of the samples of the different time points
  for(gene in gene_names){
    mean_time <- gene_normal_data[[gene]]$mean
    sd_time <- gene_normal_data[[gene]]$sd
    
    for(name in names(data)){
      data[[name]][gene] = (data[[name]][gene] - mean_time)/sd_time
    }
  }
  
  return(data)
}


#FUNCTION 3: posterior probability distribution for test sample of one gene
testPosteriorDist <- function (test_time,params, dataset, kernel="nonstationary",test_gene){
  #params - hyper parameters for a particular gene that was generated previously
  #dataset - data of training data for gene i after ExtractData function
  #kernel - nonstationary kernel default or else the same kernel that was used for hyper-parameter generation
  #test_gene - single test point whose time is not known, gi* is value of gene i at test point
  #test_time - time that we want to predict
  
  #calculating the mean of the posterior distribution
  mu_test = 0 #mean expression of gene i at test point according to the prior
  mu_train = 0  #mean expression of the training gene expression data under the prior
  gi = dataset$Y  #expression value of the gene i as a function in time
  matrix1 = covarianceMatrix(params , test_time, dataset$X, NULL, kernel)         # K(t*,t)
  matrix2 = covarianceMatrix(params , dataset$X, dataset$X, dataset$S, kernel)    # K(t,t)
  matrix3 = covarianceMatrix(params , test_time,test_time, NULL ,kernel)          # K(t*,t*)
  
  mean_post = as.numeric(mu_test + matrix1 %*% calcInverse(matrix2) %*% (gi - mu_train))  #mean of the posterior distribution
  
  variance_post = as.numeric(matrix3 - matrix1 %*% calcInverse(matrix2) %*% t(matrix1))   #variance of the posterior distribution 
  
  if(is.null(test_gene)){density = 0} #check this
  else{density =  dnorm(test_gene,mean = mean_post,sd = sqrt(variance_post))}
  
  res <- list(density,mean_post,variance_post)
  return(res)
}

#FUNCTION 4: total log likelihood for a test sample (considering all genes)
totalLogLikelihood <- function(params,dataset,test_time,test_data){
  #params are the list of hyper-parameters stored for the oscillatory genes
  #dataset of genes that are finally classified as oscillatory
  #test_time is the time of the sample that we want to predict
  #test_data is a vector of the gene expression of the oscillatory genes in the test data for one sample
  
  tll=0  #we will find the total log likelihood of observing a given test gene expression point
  for (i in 1:length(dataset)){
    likelihood <- testPosteriorDist(params = params[[i]], dataset = osciData_raw[[i]], kernel="nonstationary", test_time, test_data[[i]])
    ll<-log(likelihood[[1]])
    tll=tll+ll #to calculate total ll
  }
  
  return(tll)
}

#########################################################################################################
#training PRECISE
#since the random variable changes from the raw counts to normalized counts hence we accordingly update the mean and sd of the genes over different timepoints
bootstrapped_train_data <- bootstrap(data = Raw_train,ntimes = 50,nsamples = 100,gene_names = gene_names)

#normalizing the data
odeGP_train_input <- list()

for(gene in gene_names){
  gene_mean <- c()
  gene_sd <- c()
  for(num in 1:length(bootstrapped_train_data)){
    df <- bootstrapped_train_data[[num]]
    
    gene_mean <- c(gene_mean,mean(df[[gene]]))
    gene_sd <- c(gene_sd,sd(df[[gene]]))
  }
  
  X = train_phase
  
  Y = (gene_mean - mean(gene_mean))/sd(gene_mean)
  
  S = gene_sd/sd(gene_mean)
  
  df_temp <- data.frame(X = as.numeric(X),Y=Y,S=S)
  odeGP_train_input[[gene]] <- df_temp
}


#hyper parameter optimisation of the non-stationary kernel (ODeGP based fitting of the oscillatory genes as function of training phase)

threshold <- 14 #user defined
nonstat <- TRUE
kernel<-"nonstationary"  #if nonstat=TRUE then set to non-stationary or else change to kernel in use

osciData_raw<-list()      #list of normalised gene data that are actually classified as oscillatory by the GP algorithm
hyper_params_value<-list()  #list to contain the hyper parameters for non-stationary kernel for the different genes present in data

m=1
for(len in 1:length(odeGP_train_input)){
  data <- odeGP_train_input[[len]] 
  print(names(odeGP_train_input)[[len]]) #printing the gene name 
  df <- extractData(data,errorbars = TRUE) 
  results <- oscOrNot(df, nonstat, threshold=14, plotting=TRUE, Q=1)
  if(results$bayesF>threshold){ #if the gene data is oscillatory then hyper-params of the non-stationary kernel recorded
    hyper_params_value[[names(odeGP_train_input)[len]]]<-results[[2]]
    osciData_raw[[names(odeGP_train_input)[len]]]<-df #normalized data for the oscillatory genes recorded 
    m=m+1
  }  
}

#########################################################################################################
#validating PRECISE
NSAMPLES <- c(1,3,10,20,50,70,100,150)   #different sample number we want to try out (user defined)
NSAMPLES <- c(1)
test_index <- seq(2,4)  #center time samples being tested because true phase estimation suffers from edge effects

All_AVG_errors <- list()  #list to contain the errors for each time point different sample size - AVG for the local maximas of the likelihood curve
Actual_phase_AVG <-list() #list to contain the AVG predicted phase for each time point different sample sizes
Likelihood_curves <- list() #list to contain the likelihood curves for each time point different sample sizes

for (size in NSAMPLES){
  #TEST DATA GENERATION
  test_data <- bootstrap(data = Raw_test ,nsamples = size,ntimes = 50,gene_names = gene_names)
  test_data <- data_normalize_func(data = test_data, gene_names = gene_names)
  
  test_errors_AVG <- list() #list to contain the AVG error for the 50 samples of a different test time points
  predicted_phase_AVG <- list() #list to contain the actual AVG predicted phases for the 50 samples of the different test time points
  ll_curves <- list()       #list to contain the ll curves for the 50 samples of the different test time points
  
  for(index in test_index){
    tp = test_phase[[index]]
    df = test_data[[as.character(tp)]]
    
    errors_AVG <- c()  #vector to contain the AVG error for the 50 samples of current test point
    phases_test_AVG <- c() #vector to contain the AVG phase estimate for the 50 samples of current test point
    ll_phase <- list() #list to contain the ll curves for the 50 samples of the current test point
    
    for(sample in 1:nrow(df)){
      print(test_time[[index]])  
      print(size)
      print(sample)
      print(tp)
      
      test_sample <- df[sample,]
      test_sample$Time <- NULL
      test_data_sample <- as.list(test_sample)  
      
      #calculating the AVG error
      #plotting likelihood versus phase graph
      phases <- seq(min_train_phase,max_train_phase,0.05)
      sample_tll <- c()
      for(phase in phases){
        sample_tll_tp <- totalLogLikelihood(params = hyper_params_value,dataset = osciData_raw,test_time = phase,test_data = test_data_sample)
        sample_tll <- c(sample_tll,exp(sample_tll_tp))
      }
      
      df_ll <- data.frame(Phase = phases,Likelihood = sample_tll)
      
      x_ll <- df_ll$Likelihood
      
      df_ll$scaled_Likelihood <- (x_ll - min(x_ll)) / (max(x_ll) - min(x_ll))
      
      local_max_ll <- df_ll[["Phase"]][ggpmisc::find_peaks(df_ll$scaled_Likelihood, span = 5,global.threshold = 0.0,
                                                           local.threshold = 0.05)]  #function to detect the phases at local maximas in the likelihood curves
      
      max_likelihood <- max(df_ll$scaled_Likelihood) 
      max_ll_peak <- df_ll$Phase[which(df_ll$scaled_Likelihood == max(df_ll$scaled_Likelihood))] #detecting the phase corresponding to maximum likelihood 
      
      local_max_ll = union(local_max_ll,max_ll_peak) #including all the detected phase as our training data is over 48 hours here hence multiple peaks are expected
      
      phases_test_AVG <- c(phases_test_AVG,mean.circular(local_max_ll))  #actual AVG predicted phase
      
      #likelihood curve plot for the sample
      plot <- ggplot(data=df_ll,aes(x=Phase,y=Likelihood))+geom_line()+labs(x="Phase",y="Likelihood")+theme_bw()+
        theme(text = element_text(size=5))+geom_vline(xintercept = local_max_ll)
      ll_phase[[sample]] <- plot
      
      error_local_max <- c()
      for(pp in local_max_ll){
        diff_lm <- abs(pp - tp)
        if(diff_lm <= pi){
          error_phase_ll <- diff_lm*(24/(2*pi))
          error_local_max <- c(error_local_max,error_phase_ll)
        }
        if(diff_lm > pi){
          error_phase_ll <- abs(2*pi-diff_lm)*(24/(2*pi))
          error_local_max <- c(error_local_max,error_phase_ll)
        }
      }
      
      errors_AVG <- c(errors_AVG,mean(error_local_max,na.rm=TRUE))
      
    }
    test_errors_AVG[[as.character(test_time[[index]])]] <- errors_AVG
    predicted_phase_AVG[[as.character(test_time[[index]])]] <- phases_test_AVG
    ll_curves[[as.character(test_time[[index]])]] <- ll_phase
    
  }
  All_AVG_errors[[as.character(size)]] <- test_errors_AVG
  Actual_phase_AVG[[as.character(size)]] <- predicted_phase_AVG
  Likelihood_curves[[as.character(size)]] <- ll_curves
}


#plotting the AVG errors  graph
errors_test_AVG <- list()
for(tp in test_time[test_index]){
  error_timepoint <- c()
  sample_size <- c()
  for(sample in NSAMPLES){
    error_timepoint <- c(error_timepoint,All_AVG_errors[[as.character(sample)]][[as.character(tp)]])
    sample_size <- c(sample_size,replicate(50,sample))
  }
  
  df1 <- data.frame(sample_size = sample_size,Error = error_timepoint)
  errors_test_AVG[[as.character(tp)]] <- df1
}

graphs_AVG <- list()
for(i in 1:length(errors_test_AVG)){
  p <- ggplot(data=errors_test_AVG[[i]],aes(x=as.factor(sample_size),y=Error,color=as.factor(sample_size)))+
    geom_boxplot()+scale_y_continuous(breaks=seq(0,12.5,1),limits = c(0,12.5))+
    geom_jitter(shape=16, position=position_jitter(0.2),size=1)+theme_bw()+
    theme(text = element_text(size = 10),aspect.ratio = 1)+ 
    labs(title=paste("Deviations for test time ",names(errors_test_AVG)[[i]],sep=""),x="No.of cells averaged over",y="Deviation from true phase(Hrs)")+
    guides(color = guide_legend(title = "Sample Size"))
  graphs_AVG[[names(errors_test_AVG)[[i]]]] <- p
}

ggarrange (plotlist = graphs_AVG)

end_time <- Sys.time()
print(end_time-start_time)



