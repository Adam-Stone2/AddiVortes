if (!requireNamespace("randomForest", quietly = TRUE)) {
  install.packages("randomForest")
}

if (!requireNamespace("plotrix", quietly = TRUE)) {
  install.packages("plotrix")
}

if (!requireNamespace("doParallel", quietly = TRUE)) {
  install.packages("doParallel")
}

if (!requireNamespace("invgamma", quietly = TRUE)) {
  install.packages("invgamma")
}

if (!requireNamespace("FNN", quietly = TRUE)) {
  install.packages("FNN")
}

if (!requireNamespace("BayesTree", quietly = TRUE)) {
  install.packages("BayesTree")
}

if (!requireNamespace("gbm", quietly = TRUE)) {
  install.packages("gbm")
}

if (!requireNamespace("SoftBart", quietly = TRUE)) {
  install.packages("SoftBart")
}

library("plotrix")
library(parallel)
library(doParallel) 
library(parallel)
library(foreach)
library(SoftBart)
library(gbm)
library(BayesTree)
library(invgamma)
library(FNN)
library(randomForest)


figure2<-function(list_of_datasets){
  library(doParallel) 
  library(foreach)
  library(randomForest)
  
  par(mfrow =c(1,1))
  # Number of cross-validation folds
  # Number of cross-validation folds
  NumOfRep <- 20
  num_folds <- 5
  
  AddiVortes_Algorithm<-function(y,x,m = 200 ,max_iter = 1200,burn_in= 200,nu = 6,q =0.85,k = 3 ,var = 0.8 ,Omega = 3,lambda_rate = 25,YTest,XTest,IntialSigma = "Linear"){
    
    #Scaling x and y
    yScaled=(y-(max(y)+min(y))/2)/(max(y)-min(y))
    xScaled=x;
    for (i in 1:length(x[1,])){
      xScaled[,i]=(x[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
    }
    
    
    for (i in 1:length(XTest[1,])){
      XTest[,i]=(XTest[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
    }
    
    #Initialize:
    #Prediction Set (A list of vectors with the output values for each tessellation), 
    #Dimension set (A list of vectors with the covariates included in the tessellaions);
    #and Tessellation Set (A list of matrices that give the coordinates of the centers in the tessellations)
    
    Pred<-rep(list(matrix(mean(yScaled)/m)),m)
    Dim=vector(length = m)
    Tess=vector(length = m)
    for (i in 1:m){
      Dim[i]<-list(sample(1:length(x[1,]), 1))
      Tess[i]<-(list(matrix(rnorm(1,0,var))))
    }
    
    #Prepare some variables used in the backfitting algorithm
    SumOfAllTess=rep(mean(yScaled),length(yScaled))
    SigmaSquaredMu=(0.5/(k*sqrt(m)))^2
    LastTessPred=matrix
    
    #Matrices that will hold the samples from the poseterior distribution for the training samples and test samples.
    PredictionMatrix<-array(dim=c(length(y),(max_iter-burn_in)))
    TestMatrix<-array(dim=c(length(YTest),(max_iter-burn_in)))
    
    #finding lambda
    if (IntialSigma=="Naive"){ # Usually used if p is greater then n. Uses Standard deviation of y to predict sigma.
      SigmaSquaredHat=var(yScaled)
    }
    else{  # Default method using residual standard deviation from a least-squared linear regression of y, to predict sigma.
      MultiLinear<-lm(yScaled ~ xScaled)
      SigmaSquaredHat=sum(MultiLinear$residuals^2)/(length(yScaled)-length(xScaled[1,])-1)
    }
    
    #Find lambda
    lambda=1;
    lambda <- optim(par = 1,
                    fitting_function,
                    method = "Brent",
                    lower = 0.001,
                    upper = 100,
                    q=q , nu=nu, sigmaSquared_hat=SigmaSquaredHat)$par
    
    for (i in 1:max_iter){
      
      #Sample Sigma squared using all tessellations to predict the outcome variables
      SigmaSquared=SigmaSquaredCalculation(yScaled,nu,lambda,SumOfAllTess)
      
      for (j in 1:m){
        NewTessOutput<-NewTess(xScaled,j,Tess,Dim,var) #Propose new Tessellation 
        TessStar<-NewTessOutput[[1]]  
        DimStar<-NewTessOutput[[2]]
        Modification<-NewTessOutput[[3]]
        
        ResidualsOutput<-CalculateResiduals(yScaled,xScaled,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred) #Calculate the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation and the number of observations in each cell.
        R_ijOld<-ResidualsOutput[[1]]   #Old and New refer to the original and proposed tessellations
        n_ijOld<-ResidualsOutput[[2]]
        R_ijNew<-ResidualsOutput[[3]]
        n_ijNew<-ResidualsOutput[[4]]
        SumOfAllTess<-ResidualsOutput[[5]] #Keeps track of the prediction for all tessellations to help sample sigma squared.
        IndexesStar<-ResidualsOutput[[6]] #Gives the row of each observation for the cell it falls in for the proposed tessellation.
        Indexes<-ResidualsOutput[[7]]  #Gives the row of each observation for the cell it falls in for the original tessellation.
        
        if (!any(n_ijNew==0)){ #automatically rejects proposed tessellation if there exists a cell with no observations in.
          
          LOGAcceptenceProb=AlphaCalculation(xScaled,TessStar,DimStar,j,R_ijOld,n_ijOld,R_ijNew,n_ijNew,SigmaSquared,Modification,SigmaSquaredMu,Omega,lambda_rate) #Gives the log of the acceptence probability.
          
          if (log(runif(n=1, min=0, max=1))<LOGAcceptenceProb){ #Accepts the proposed tessellation is accepted then calculates the new output values for the new tessellation. 
            Tess=TessStar
            Dim=DimStar
            Pred[[j]]=NewPredSet(j,TessStar,R_ijNew,n_ijNew,SigmaSquaredMu,SigmaSquared)
            LastTessPred=Pred[[j]][IndexesStar]
          }
          else { #Rejects the proposed tesellation then calculates new output values for the original tessellation.
            Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
            LastTessPred=Pred[[j]][Indexes];
          }
        }
        else{ #Rejects the proposed tesellation then calculates new output values for the original tessellation.
          Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
          LastTessPred=Pred[[j]][Indexes];
        }
        if (j==m){ #If j equals m then adds the last tessellation output values to give a prediction.
          SumOfAllTess=SumOfAllTess+LastTessPred;
        }
      }
      
      if (i>burn_in){ #vectors that hold the predictions for each iteration after burn in.
        PredictionMatrix[,(i-burn_in)]=SumOfAllTess;
        TestMatrix[,(i-burn_in)]=TestPrediction(XTest,m,Tess,Dim,Pred);
      }
    }
    
    #finding the mean of the predition over the iterations and then unscaling the predictions.
    mean_yhat=(rowSums(PredictionMatrix)/(max_iter-burn_in))*(max(y)-min(y))+((max(y)+min(y))/2)
    mean_yhat_Test=(rowSums(TestMatrix)/(max_iter-burn_in))*(max(y)-min(y))+((max(y)+min(y))/2)
    
    return( #Returns the RMSE value for the test samples.
      data.frame(
        In_sample_RMSE = sqrt(mean((y-mean_yhat)^2)),
        Out_of_sample_RMSE = sqrt(mean((YTest-mean_yhat_Test)^2))     
      )
    )
  }
  
  
  SigmaSquaredCalculation<-function(yScaled,nu,lambda,SumOfAllTess){ #Sample sigma squared from inverse gamma distribution
    
    n=length(yScaled)
    SigmaSquared<-rinvgamma(1,shape=(nu+n)/2,rate=(nu*lambda+sum((yScaled-SumOfAllTess)^2))/2)
    
    return(SigmaSquared)
  }
  
  NewTess<-function(x,j,Tess,Dim,var){ #Propose a new tessellation
    
    p=runif(1,0,1) #Randomly sample p to decide the proposed modification to tessellation.
    
    DimStar=Dim # Let proposed dimension matrix equal original dimension matrix.
    TessStar=Tess #Similar for the tessellation matrix.
    
    if (p<0.2 & length(Dim[[j]])!=length(x[1,]) | length(Dim[[j]])==1 & p<0.4){ #Add a dimension if p is less then 0.2 or if p is less then 0.4 when there is only one dimension in the Tessellation due to adjustments (Supplementary Material).
      NumberOfCovariates=1:length(x[1,]) #Let NumberOfCovariates be a vector from 1 to the number of covariates considered.
      NumberOfCovariates=NumberOfCovariates[-Dim[[j]]] #Remove all values that for the covariates that are already in the tessellation.
      DimStar[[j]]<-c(Dim[[j]],sample(NumberOfCovariates,1)) # Uniformly sample a new covariate and add it to the dimension matrix.
      TessStar[[j]]=cbind(Tess[[j]],rnorm(length(Tess[[j]][,1]),0,var)) # Sample new coordinates from Normal distribution for the new dimension and add it to the Tessellation matrix.
      Modification="AD"}
    else if (p<0.4){ #Remove a dimension if p is less then 0.4.
      RemovedDim=sample(1:length(Dim[[j]]),1) #Uniformly sample the dimension to be removed.
      DimStar[[j]]=DimStar[[j]][-RemovedDim] #Remove the dimension from the dimesion Matrix.
      TessStar[[j]]=matrix(TessStar[[j]][,-RemovedDim],ncol=length(DimStar[[j]])) #Remove the coordinates in the Tessellation matrix corresponding to the dimension removed.
      Modification="RD"}
    else if (p<0.6 || p<0.8 & length(Tess[[j]][,1])==1){ #Add a centre if p is less then 0.6 or if p is less then 0.4 when there is only one center in the Tessellation due to adjustments (Supplementary Material).
      TessStar[[j]]=rbind(Tess[[j]],rnorm(length(Dim[[j]]),0,var)) #Add a new row of coordinates, sampled from a normal distribution, to the Tessellation matrix to add a center.
      Modification="AC"}
    else if (p<0.8){ #Add a centre if p is less then 0.8. 
      CenterRemoved=sample(1:length(TessStar[[j]][,1]),1) #Sample a row.
      TessStar[[j]]=matrix(TessStar[[j]][-CenterRemoved,],ncol=length(Dim[[j]])) #Remove row sampled.
      Modification="RC"}
    else if (p<0.9 || length(Dim[[j]])==length(x[1,])){ #Change a center if p is less then 0.9 or if the all the covariates are in the tessellation.
      TessStar[[j]][sample(1:length(TessStar[[j]][,1]),1),]=rnorm(length(Dim[[j]]),0,var) # Sample a row in the tessellaion matrix and change the coordinates of the centre by sampling from a normal distribution.
      Modification="Change"}
    else{ #Swop a dimension.
      NumberOfCovariates=1:length(x[1,])  #Let NumberOfCovariates be a vector from 1 to the number of covariates considered.
      NumberOfCovariates=NumberOfCovariates[-Dim[[j]]]  #Remove all values that for the covariates that are already in the tessellation.
      DimToChange=sample(1:length(Dim[[j]]),1) #Uniformly sample a dimension to change.
      DimStar[[j]][DimToChange]=sample(NumberOfCovariates,1) #Replace the Dimension to a new uniforly sampled covariate that is not already in the tessellaion.
      TessStar[[j]][,DimToChange]=rnorm(length(Tess[[j]][,1]),0,var) #Add new normally sampled coordinates new dimension added.
      Modification="Swop"}
    
    TessStar[[j]]<-matrix(TessStar[[j]],ncol=length(DimStar[[j]])) #Ensure the the Tessellation matrix is a "matrix" type.
    
    return(list(TessStar,DimStar,Modification)) #Return new proposed tessellation.
  }
  
  fitting_function<- function(lambda,q,nu,sigmaSquared_hat){ #function that calculates the squared difference between sigma squared hat and the inverse gamma function
    return((sigmaSquared_hat- qinvgamma(q, shape=nu/2, rate=nu*lambda/2))^2)
  }
  
  Indexes<-function(x,Tess,Dim){ #Gives the row (the center) of the tessellation that each obseravtion falls within.
    if (length(Tess[,1])==1){ #only 1 centre
      CellsForGivenTess=rep(1,length(x[,1]))
    }
    else{ #multiple
      CellsForGivenTess=knnx.index(Tess,matrix(x[,Dim],ncol = length(Dim)),1)
    }
    return(CellsForGivenTess)
  }
  
  AlphaCalculation<-function(x,Tess,Dim,j,R_ijOld,n_ijOld,R_ijNew,n_ijNew,SigmaSquared,Modification,SigmaSquaredMu,Omega,lambda_rate){ #Calculates the acceptence rate of the proposed tessellation.
    
    d=length(Dim[[j]]);
    NumCovariates=length(x[1,]);
    cStar=length(Tess[[j]][,1]);
    
    #The Log Likelihood Ratio in the acceptence ratio
    LOGlikelihoodRatio=0.5*(log(prod(n_ijOld*SigmaSquaredMu+SigmaSquared))-log(prod(n_ijNew*SigmaSquaredMu+SigmaSquared)))+((SigmaSquaredMu/(2*SigmaSquared))*(-sum((R_ijOld^2)/(n_ijOld*SigmaSquaredMu+SigmaSquared))+sum((R_ijNew^2)/(n_ijNew*SigmaSquaredMu+SigmaSquared))))
    
    #Calculating the acceptence probablity for "AD"=Adding a dimension, "RD"=Removing a dimension, "AC"=Adding a center, "RC"=Removing a center, "Change"=Changing the coordinates of a center and Swopping a dimension.
    if (Modification == "AD"){ 
      TessStructure=(dbinom(d-1,NumCovariates-1,Omega/NumCovariates))/(dbinom(d-2,NumCovariates-1,Omega/NumCovariates)*(NumCovariates-d+1))
      TransitionRatio=(NumCovariates-d+1)/d;
      AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
      
      #Adjustments.
      if (length(Dim[[j]])==1){
        AcceptenceProb=AcceptenceProb+log(1/2)
      }
      else if (length(Dim[[j]])==NumCovariates-1){
        AcceptenceProb=AcceptenceProb+log(2)}
    }
    else if (Modification == "RD"){
      TessStructure=(dbinom(d-1,NumCovariates,Omega/NumCovariates)*(NumCovariates-d))/(dbinom(d,NumCovariates,Omega/NumCovariates))
      TransitionRatio=(d+1)/(NumCovariates-d)
      AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
      
      #Adjustments.
      if (length(Dim[[j]])==NumCovariates){
        AcceptenceProb=AcceptenceProb+log(1/2)
      }
      else if (length(Dim[[j]])==2){
        AcceptenceProb=AcceptenceProb+log(2);
      }
    }
    else if(Modification == "AC"){
      TessStructure=dpois(cStar-1,lambda_rate)/dpois(cStar-2,lambda_rate)
      TransitionRatio=1/cStar;
      AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)+0.5*log(SigmaSquared)
      
      #Adjustments.
      if (cStar==1){
        AcceptenceProb=AcceptenceProb+log(1/2);
      }
    }
    else if (Modification == "RC"){
      TessStructure=dpois(cStar-1,lambda_rate)/dpois(cStar,lambda_rate);
      TransitionRatio=cStar+1;
      AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)-0.5*log(SigmaSquared)
      
      #Adjustments.
      if (cStar==2){
        AcceptenceProb=AcceptenceProb+log(2);
      }
    }
    else if (Modification == "Change"){
      TessStructure=1;
      TransitionRatio=1;
      AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio);
    }
    else {
      TessStructure=1;
      TransitionRatio=1;
      AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio);
    }
    return(AcceptenceProb)
  }
  
  CalculateResiduals<-function(y,x,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred){ #A function that calculates the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation.
    if (j==1){
      indexes=Indexes(x,Tess[[j]],Dim[[j]]);
      CurrentTessPred<-Pred[[j]][indexes]
      SumOfAllTess=SumOfAllTess-CurrentTessPred}
    else{
      indexes=Indexes(x,Tess[[j]],Dim[[j]]);
      CurrentTessPred<-Pred[[j]][indexes]
      SumOfAllTess=SumOfAllTess+LastTessPred-CurrentTessPred;
    }
    
    IndexesStar=Indexes(x,TessStar[[j]],DimStar[[j]]);
    R_j<-y-SumOfAllTess
    
    #Initializing Sizes
    
    R_ijOld=rep(0,length(Pred[[j]]))
    n_ijOld=rep(0,length(Pred[[j]]))
    
    for (i in 1:length(Pred[[j]])){
      R_ijOld[i]<-sum(R_j[indexes==i])
      n_ijOld[i]<-sum(indexes==i)
    }
    
    R_ijNew=rep(0,length(TessStar[[j]][,1]))
    n_ijNew=rep(0,length(TessStar[[j]][,1]))
    
    for (i in 1:length(TessStar[[j]][,1])){
      R_ijNew[i]<-sum(R_j[IndexesStar==i])
      n_ijNew[i]<-sum(IndexesStar==i)}
    
    return(list(R_ijOld,n_ijOld,R_ijNew,n_ijNew,SumOfAllTess,IndexesStar,indexes))
  }
  
  NewPredSet<-function(j,Tess,R_ijNew,n_ijNew,sigmaSquaredMu,SigmaSquared){ #Sampling the new output values for the new tessellation.
    PredSet=rep(0,length(Tess[[j]][,1]))
    for (i in 1:length(Tess[[j]][,1])){
      PredSet[i]=rnorm(1,sigmaSquaredMu*R_ijNew[i]/(sigmaSquaredMu*n_ijNew[i]+SigmaSquared),((SigmaSquared*sigmaSquaredMu)/(n_ijNew[i]*sigmaSquaredMu+SigmaSquared))^0.5);
    }
    return(PredSet)
  }
  
  TestPrediction<-function(x,m,Tess,Dim,Pred){ #A function that derives a prediction for the output variable give an observation, for one iteration in AddiVortes.
    Prediction=rep(0,length(x[,1]));
    for (j in 1:m){
      NewTessIndexes=Indexes(x,Tess[[j]],Dim[[j]]);
      Prediction=Prediction+Pred[[j]][NewTessIndexes]
    }
    return(Prediction)
  }
  
  createFolds <- function(data, k) {
    set.seed(123)  # For reproducibility
    fold_indices <- sample(rep(1:k, length.out = nrow(data)))
    fold_list <- split(1:nrow(data), fold_indices)
    return(fold_list)
  }
  
  # Create a data frame with all combinations of hyperparameter values
  hyperparametersBART <- expand.grid(
    m=c(50,200),
    nu_q=list(c(3,0.9),c(3,0.99),c(10,0.75)),
    k=c(1,2,3,5)
  )
  
  hyperparametersAddiVortes <- expand.grid(
    m=c(50,200),
    nu=c(6),
    q=c(0.85),
    k=c(1,3),
    sd=c(0.8,1.5),
    omega=c(3),
    lambda=c(5,25)
  )
  
  hyperparametersBoost <- expand.grid(
    depth=c(1,2,3,4),
    num_trees=c(50,100,200),
    eta=c(0.01, 0.05, 0.10, 0.25)
  )
  
  hyperparametersSoftBART <- expand.grid(
    m=c(50,200),
    k=c(1,2,3,5)
  )
  
  # Function to remove first and last columns
  remove_first_last_columns <- function(mat) {
    mat[, -c(1, ncol(mat))]
  }
  
  extract_last_column <- function(mat) {
    mat[, ncol(mat)]
  }
  
  # Apply the function to each matrix in the list
  BenchmarkX <- lapply(list_of_datasets, remove_first_last_columns)
  BenchmarkY<- lapply(list_of_datasets, extract_last_column)
  
  num_cores <- 10 # Specify the number of cores you want to use
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  for (l in 1:length(list_of_datasets)){
    set.seed(324)
    
    X<-BenchmarkX[[l]]
    Y<-BenchmarkY[[l]]
    Y<-as.numeric(as.matrix(Y))
    
    n=length(Y)
    TrainSet<-matrix(nrow=trunc(5*n/6),ncol = NumOfRep)
    TestSet<-matrix(nrow=n-trunc(5*n/6),ncol = NumOfRep)
    
    for (h in 1:NumOfRep){
      TrainSet[,h]=sort(sample.int(n,5*n/6));
      TestSetInit<-1:n
      TestSet[,h]=TestSetInit[! TestSetInit %in% TrainSet[,h]];
      
    }  
    
    #Cross validation
    folds <- createFolds(X, k = num_folds)  # 5-fold cross-validation
    
    # Initialize a matrix to store the results
    results_matrixAddiVortes <- matrix(ncol = 8, nrow = nrow(hyperparametersAddiVortes))
    rmse_values<-vector(length = num_folds)
    
    # Iterate over hyperparameter combinations
    results_matrixAddiVortes<-foreach(i = 1:nrow(hyperparametersAddiVortes), .combine = rbind) %dopar% {
      library('invgamma')
      library('FNN')
      library('expint')
      
      params <- hyperparametersAddiVortes[i, ]
      rmse_values <- numeric(num_folds)
      
      # Perform cross-validation using a loop
      for (j in 1:num_folds) {
        TrainSetCV <- unlist(folds[-j],use.names = FALSE)
        TestSetCV <- folds[[j]]
        
        m<-params$m
        nu <- params$nu
        q<- params$q
        sd<- params$sd
        omega<- params$omega
        lambda<- params$lambda
        k<- params$k
        
        # Make predictions on the test set
        rmse_values[j] <- AddiVortes_Algorithm(Y[TrainSetCV],as.matrix(X[TrainSetCV,]),m,1200,200,nu,q,k,sd,omega,lambda,Y[TestSetCV],as.matrix(X[TestSetCV,]))$Out_of_sample_RMSE
      }                   
      # Store the average RMSE from cross-validation
      return(c(m,nu,q,k,sd,omega,lambda,mean(rmse_values)))
      
    }
    
    # Find the best hyperparameter combination based on RMSE
    best_index <- which.min(results_matrixAddiVortes[, 8])
    best_hyperparametersAddiVortes <- results_matrixAddiVortes[best_index, 1:7]
    
    print("Finished AddiVortes Tuning")
    
    # Initialize a matrix to store the results
    results_matrixBART <- matrix(ncol = 5, nrow = nrow(hyperparametersBART))
    
    # Iterate over hyperparameter combinations
    results_matrixBART<-foreach (i = 1:nrow(hyperparametersBART), .combine = rbind) %dopar% {
      library('BayesTree')
      
      params <- hyperparametersBART[i, ]
      rmse_values <- numeric(num_folds)
      
      # Perform cross-validation using a loop
      for (j in 1:num_folds) {
        TrainSetCV <- unlist(folds[-j],use.names = FALSE)
        TestSetCV <- folds[[j]]
        
        m <- params$m
        nu <- params$nu_q[[1]][1]
        q <- params$nu_q[[1]][2]
        k <- params$k
        
        # Train the bart model
        BartOr<-bart(as.matrix(X[TrainSetCV,]),as.numeric(as.matrix(Y[TrainSetCV])),as.matrix(X[TestSetCV,]), ntree = m, sigdf = nu,sigquant = q,k = k)
        
        # Make predictions on the test set
        rmse_values[j] <- (mean((Y[TestSetCV]-BartOr$yhat.test.mean)**2))**0.5
      }
      
      # Store the average RMSE from cross-validation
      results_matrixBART[i,] <- cbind(params$m,params$nu_q[[1]][1],params$nu_q[[1]][2],params$k,mean(rmse_values))
      
    }
    
    print('Finished BART Tuning')
    
    # Find the best hyperparameter combination based on RMSE
    best_index <- which.min(results_matrixBART[, 5])
    best_hyperparametersBART <- results_matrixBART[best_index, 1:4]
    
    
    results_matrixSoftBART <- matrix(ncol = 3, nrow = nrow(hyperparametersSoftBART))
    
    results_matrixSoftBART<-foreach (i = 1:nrow(hyperparametersSoftBART), .combine = rbind) %dopar% {
      library(SoftBart)
      
      params <- hyperparametersSoftBART[i, ]
      rmse_values <- numeric(num_folds)
      
      # Perform cross-validation using a loop
      for (j in 1:num_folds) {
        TrainSetCV <- unlist(folds[-j],use.names = FALSE)
        TestSetCV <- folds[[j]]
        
        m <- params$m
        ##nu <- params$nu_q[[1]][1]
        ##q <- params$nu_q[[1]][2]
        k <- params$k
        
        df <- data.frame(X = X[TrainSetCV,], Y = Y[TrainSetCV])
        df_test <- data.frame(X = X[TestSetCV,], Y = Y[TestSetCV])
        
        # Train the bart model
        SoftBARTOr<-softbart_regression(Y ~ ., df, df_test, num_tree = m,k =k)
        
        # Make predictions on the test set
        rmse_values[j] <- sqrt(mean((SoftBARTOr$mu_test_mean-Y[TestSetCV])^2))
      }
      
      # Store the average RMSE from cross-validation
      results_matrixSoftBART[i,] <- cbind(params$m,params$k,mean(rmse_values))
      
    }
    
    print('Finished SoftBART Tuning')
    
    # Find the best hyperparameter combination based on RMSE
    best_index <- which.min(results_matrixSoftBART[,3])
    best_hyperparametersSoftBART <- results_matrixSoftBART[best_index, 1:2]
    
    # Initialize a matrix to store the results for each combination
    results_matrixRF <- matrix(ncol = 2, nrow = 4)
    PercentOfVariable = c(0.1,0.25,0.5,1)
    for (i in 1:4) {
      # Perform cross-validation using a loop
      fold_RMSE <- numeric(num_folds)
      for (j in 1:num_folds) {
        TrainSetCV <- unlist(folds[-j],use.names = FALSE)
        TestSetCV <- folds[[j]]
        
        mean.y.hat<-as.vector(randomForest(X[TrainSetCV,],Y[TrainSetCV],X[TestSetCV,],Y[TestSetCV],mtry=round(PercentOfVariable[i]*length(X[1,])))$test$predicted)
        fold_RMSE[j] <- sqrt(mean((Y[TestSetCV]-mean.y.hat)^2))
        
      }
      # Store the average RMSE from cross-validation
      results_matrixRF[i, 1] <- mean(fold_RMSE)
      results_matrixRF[i,2] <- PercentOfVariable[i]
    }
    ordered_indices <- order(results_matrixRF[, 1])
    ordered_results <- results_matrixRF[ordered_indices, ]
    PercentOfVariable<-ordered_results[1,2]
    
    print('Finished RF tuning')
    
    # Initialise results matrix
    results_matrixGBM <- matrix(ncol = 4, nrow = nrow(hyperparametersBoost))
    
    # Iterate over hyperparameter combinations
    results_matrixGBM <- foreach(i = 1:nrow(hyperparametersBoost), .combine = rbind) %dopar% {
      params <- hyperparametersBoost[i, ]
      rmse_values <- numeric(num_folds)
      library(gbm)
      
      # Perform cross-validation
      for (j in 1:num_folds) {
        TrainSetCV <- unlist(folds[-j], use.names = FALSE)
        TestSetCV <- folds[[j]]
        
        num_trees <- params$num_trees
        depth <- params$depth
        eta <- params$eta
        
        # Train the gbm model
        gbm_model <- gbm(
          formula = Y ~ .,
          distribution = "gaussian",
          data = data.frame(X = X[TrainSetCV, ], Y = Y[TrainSetCV]),
          n.trees = num_trees,
          interaction.depth = depth,
          shrinkage = eta,
          train.fraction = 1.0,
          verbose = FALSE
        )
        
        # Make predictions on the test set
        yhat_test <- predict(gbm_model, newdata = data.frame(X = X[TestSetCV, ]), n.trees = num_trees)
        rmse_values[j] <- sqrt(mean((Y[TestSetCV] - yhat_test) ^ 2))
      }
      
      # Store the average RMSE from cross-validation
      c(params$num_trees, params$depth, params$eta, mean(rmse_values))
    }
    
    print('Finished GBM Tuning')
    
    # Find the best hyperparameters
    best_hyperparametersBoost <- results_matrixGBM[which.min(results_matrixGBM[, 4]), 1:3]
    
    
    
    CV.AddiVortes.RMSE<-foreach(k = 1:NumOfRep, .combine = cbind) %dopar% {
      library('invgamma')
      library('FNN')
      
      AddiVortes.RMSE <- AddiVortes_Algorithm(Y[TrainSet[,k]],as.matrix(X[TrainSet[,k],]),best_hyperparametersAddiVortes[1],1200,200,best_hyperparametersAddiVortes[2],best_hyperparametersAddiVortes[3],best_hyperparametersAddiVortes[4],best_hyperparametersAddiVortes[5],best_hyperparametersAddiVortes[6],best_hyperparametersAddiVortes[7],Y[TestSet[,k]],as.matrix(X[TestSet[,k],]))$Out_of_sample_RMSE
    }
    
    Default.AddiVortes.RMSE<-foreach(k = 1:NumOfRep, .combine = cbind) %dopar% {
      library('invgamma')
      library('FNN')
      
      AddiVortes.RMSE <- AddiVortes_Algorithm(Y[TrainSet[,k]],as.matrix(X[TrainSet[,k],]),200,1200,200,6,0.85,3,0.8,3,25,Y[TestSet[,k]],as.matrix(X[TestSet[,k],]))$Out_of_sample_RMSE
    }
    
    BART.RMSE<-foreach(k = 1:NumOfRep, .combine = cbind) %dopar% {
      library('BayesTree')
      library('expint')
      
      #bartMAchine
      modelBART<-bart(as.matrix(X[TrainSet[,k],]),as.numeric(as.matrix(Y[TrainSet[,k]])),as.matrix(X[TestSet[,k],]), ntree = best_hyperparametersBART[1], sigdf = best_hyperparametersBART[2],sigquant = best_hyperparametersBART[3],k = best_hyperparametersBART[4])
      BART.RMSE<-(mean((Y[TestSet[,k]]-modelBART$yhat.test.mean)**2))**0.5  
    }
    
    
    SoftBART.RMSE<-foreach(k = 1:NumOfRep, .combine = cbind) %dopar% {
      library(SoftBart)
      
      df <- data.frame(X = X[TrainSet[,k],], Y = Y[TrainSet[,k]])
      df_test <- data.frame(X = X[TestSet[,k],], Y = Y[TestSet[,k]])
      
      # Train the bart model
      modelSoftBART<-softbart_regression(Y ~ ., df, df_test, num_tree = best_hyperparametersSoftBART[1],k =best_hyperparametersSoftBART[2])
      
      # Make predictions on the test 
      SoftBART.RMSE<-sqrt(mean((modelSoftBART$mu_test_mean-Y[TestSet[,k]])^2))
    }
    
    RForest.RMSE<-foreach(k = 1:NumOfRep, .combine = cbind) %dopar% {
      library('randomForest')
      #randomForests
      
      RandomF<-randomForest(X[TrainSet[,k],],Y[TrainSet[,k]],X[TestSet[,k],],Y[TestSet[,k]],mtry = round(PercentOfVariable*length(X[1,])))
      RForest.RMSE<-(mean((Y[TestSet[,k]]-RandomF$test$predicted)**2))**0.5
    }
    
    # Perform final evaluation with the best hyperparameters
    gbm.RMSE <- foreach(k = 1:NumOfRep, .combine = cbind) %dopar% {
      # Train the gbm model with best hyperparameters
      library(gbm)
      final_gbm_model <- gbm(
        formula = Y ~ .,
        distribution = "gaussian",
        data = data.frame(X = X[TrainSet[, k], ], Y = Y[TrainSet[, k]]),
        n.trees = best_hyperparametersBoost[1],
        interaction.depth = best_hyperparametersBoost[2],
        shrinkage = best_hyperparametersBoost[3],
        train.fraction = 1.0,
        verbose = FALSE
      )
      
      # Make predictions on the test set
      final_yhat_test <- predict(final_gbm_model, newdata = data.frame(X = X[TestSet[, k], ]), n.trees = best_hyperparametersBoost[1])
      sqrt(mean((Y[TestSet[, k]] - final_yhat_test) ^ 2))
    }
    
    
    RMSE<-rbind(CV.AddiVortes.RMSE,Default.AddiVortes.RMSE,BART.RMSE,RForest.RMSE,SoftBART.RMSE,gbm.RMSE)
    
    print(l)
    
    RRMSE<-apply(RMSE, 2, function(col) col / min(col))
    
    boxplot(RRMSE[1,],RRMSE[2,],RRMSE[3,],RRMSE[4,],RRMSE[5,],horizontal = TRUE, names = unique(c("AddiVortes CV","AddiVortes Default", "BART", "Random Forests","SoftBART")))
    
    if (l==1){
      All.RRMSE<-matrix(nrow=4,ncol=NumOfRep)
      All.RRMSE<-RRMSE
    }
    else{
      All.RRMSE<-cbind(All.RRMSE,RRMSE)
    }
  }
  
  stopCluster(cl)
  
  Smaller.RRMSE.AddiVortesCV<-All.RRMSE[1,All.RRMSE[1,]<1.5]
  Smaller.RRMSE.AddiVortesDef<-All.RRMSE[2,All.RRMSE[2,]<1.5]
  Smaller.RRMSE.BART<-All.RRMSE[3,All.RRMSE[3,]<1.5]
  Smaller.RRMSE.RF<-All.RRMSE[4,All.RRMSE[4,]<1.5]
  Smaller.RRMSE.SoftBART<-All.RRMSE[5,All.RRMSE[5,]<1.5]
  
  boxplot(as.numeric(Smaller.RRMSE.AddiVortesCV),as.numeric(Smaller.RRMSE.AddiVortesDef),as.numeric(Smaller.RRMSE.BART), as.numeric(Smaller.RRMSE.RF),horizontal = TRUE, names = unique(c("AddiVortes-CV", "AddiVortes-default","BART","Random Forests","SoftBART")))
  return(All.RRMSE)
  
}
                 
figure3<-function(){

  f = function(x){
    10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-0.5)^2+10*x[,4]+5*x[,5]
  }
  
  sigma = 1  #y = f(x) + sigma*z , z~N(0,1)
  n = 300      #number of observations
  set.seed(9)
  X=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
  Ey = f(X)
  Y=Ey+sigma*rnorm(n)
  
  n=length(Y)
  TrainSet=sort(sample.int(n,3*n/6))
  TestSet=1:n
  TestSet=TestSet[! TestSet %in% TrainSet]
  
  par(mfrow =c(1,3))
  
  AddiVortes_Algorithm_Plot<-function(y,x,m,max_iter,burn_in,nu,q,k,sd,Omega,lambda_rate,YTest,XTest,IntialSigma = "Linear"){
    
    #Scaling x and y
    yScaled=(y-(max(y)+min(y))/2)/(max(y)-min(y))
    xScaled=x;
    for (i in 1:length(x[1,])){
      xScaled[,i]=(x[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
    }
    
    
    for (i in 1:length(XTest[1,])){
      XTest[,i]=(XTest[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
    }
    
    
    #Initialize the Prediction Set, Dimension set and Tessellation Set
    
    Pred<-rep(list(matrix(mean(yScaled)/m)),m)
    Dim=vector(length = m)
    Tess=vector(length = m)
    for (i in 1:m){
      Dim[i]<-list(sample(1:length(x[1,]), 1))
      Tess[i]<-(list(matrix(rnorm(1,0,sd))))
    }
    
    #Prepare some variables used in the backfitting algorithm
    SumOfAllTess=rep(mean(yScaled),length(yScaled))
    SigmaSquaredMu=(0.5/(k*sqrt(m)))^2
    LastTessPred=matrix
    
    #Matrices that will hold the samples from the poseterior distribution for the training samples and test samples.
    PredictionMatrix<-array(dim=c(length(y),(max_iter-burn_in)))
    TestMatrix<-array(dim=c(length(YTest),(max_iter-burn_in)))
    TestMatrix<-array(dim=c(length(YTest),(max_iter-burn_in)))
    plotForSigmaSquared<-vector(length = max_iter)
    plotForRMSE<-vector(length = max_iter)
    
    #finding lambda
    if (IntialSigma=="Naive"){ # Usually used if p is greater then n. Uses Standard deviation of y to predict sigma.
      SigmaSquaredHat=var(yScaled)
    }
    else{  # Default method using residual standard deviation from a least-squared linear regression of y, to predict sigma.
      MultiLinear<-lm(yScaled ~ xScaled)
      SigmaSquaredHat=sum(MultiLinear$residuals^2)/(length(yScaled)-length(xScaled[1,])-1)
    }
    
    #Find lambda
    lambda=1;
    lambda <- optim(par = 1,
                    fitting_function,
                    method = "Brent",
                    lower = 0.001,
                    upper = 100,
                    q=q , nu=nu, sigmaSquared_hat=SigmaSquaredHat)$par
    
    for (i in 1:max_iter){
      
      #Sample Sigma squared using all tessellations to predict the outcome variables
      SigmaSquared=SigmaSquaredCalculation(yScaled,nu,lambda,SumOfAllTess)
      plotForSigmaSquared[i]=(SigmaSquared*(max(y)-min(y))^2)^0.5
      plotForRMSE[i]=(mean((yScaled-SumOfAllTess)^2))^0.5
      
      for (j in 1:m){
        NewTessOutput<-NewTess(xScaled,j,Tess,Dim,sd) #Propose new Tessellation 
        TessStar<-NewTessOutput[[1]]  
        DimStar<-NewTessOutput[[2]]
        Modification<-NewTessOutput[[3]]
        
        ResidualsOutput<-CalculateResiduals(yScaled,xScaled,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred) #Calculate the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation and the number of observations in each cell.
        R_ijOld<-ResidualsOutput[[1]]   #Old and New refer to the original and proposed tessellations
        n_ijOld<-ResidualsOutput[[2]]
        R_ijNew<-ResidualsOutput[[3]]
        n_ijNew<-ResidualsOutput[[4]]
        SumOfAllTess<-ResidualsOutput[[5]] #Keeps track of the prediction for all tessellations to help sample sigma squared.
        IndexesStar<-ResidualsOutput[[6]] #Gives the row of each observation for the cell it falls in for the proposed tessellation.
        Indexes<-ResidualsOutput[[7]]  #Gives the row of each observation for the cell it falls in for the original tessellation.
        
        if (!any(n_ijNew==0)){ #automatically rejects proposed tessellation if there exists a cell with no observations in.
          
          LOGAcceptenceProb=AlphaCalculation(xScaled,TessStar,DimStar,j,R_ijOld,n_ijOld,R_ijNew,n_ijNew,SigmaSquared,Modification,SigmaSquaredMu,Omega,lambda_rate) #Gives the log of the acceptence probability.
          
          
          if (log(runif(n=1, min=0, max=1))<LOGAcceptenceProb){ #Accepts the proposed tessellation is accepted then calculates the new output values for the new tessellation. 
            Tess=TessStar
            Dim=DimStar
            Pred[[j]]=NewPredSet(j,TessStar,R_ijNew,n_ijNew,SigmaSquaredMu,SigmaSquared)
            LastTessPred=Pred[[j]][IndexesStar]
          }
          else { #Rejects the proposed tesellation then calculates new output values for the original tessellation.
            Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
            LastTessPred=Pred[[j]][Indexes];
          }
        }
        else{ #Rejects the proposed tesellation then calculates new output values for the original tessellation.
          Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
          LastTessPred=Pred[[j]][Indexes];
        }
        if (j==m){ #If j equals m then adds the last tessellation output values to give a prediction.
          SumOfAllTess=SumOfAllTess+LastTessPred;
        }
      }
      
      if (i>burn_in){ #vectors that hold the predictions for each iteration after burn in.
        PredictionMatrix[,(i-burn_in)]=SumOfAllTess;
        TestMatrix[,(i-burn_in)]=TestPrediction(XTest,m,Tess,Dim,Pred);
      }
      if (i %% 100 == 0){
        cat(sprintf("Iteration %d out of %d", i, max_iter), "\n")
      }
    }
    
    #finding the mean of the predition over the iterations and then unscaling the predictions.
    mean_yhat=(rowSums(PredictionMatrix)/(max_iter-burn_in))*(max(y)-min(y))+((max(y)+min(y))/2)
    mean_yhat_Test=(rowSums(TestMatrix)/(max_iter-burn_in))*(max(y)-min(y))+((max(y)+min(y))/2)
    mean_yhat_Test_ColSum=(colSums(TestMatrix)/(length(TestMatrix[,1])))*(max(y)-min(y))+((max(y)+min(y))/2)
    
    LowerConfidenceTRAINValue<-vector(length=length(mean_yhat))
    UpperConfidenceTRAINValue<-vector(length=length(mean_yhat))
    
    for (i in 1:length(mean_yhat)){
      PredictionMatrix[i,]<-sort(PredictionMatrix[i,])
      
      if ((((max_iter-burn_in+1)*0.05))== round((max_iter-burn_in+1)*0.05)){
        LowerConfidenceTRAINValue[i]<-(PredictionMatrix[i,(max_iter-burn_in+1)*0.05])*(max(y)-min(y))+((max(y)+min(y))/2)
        UpperConfidenceTRAINValue[i]<-(PredictionMatrix[i,(max_iter-burn_in+1)*0.95])*(max(y)-min(y))+((max(y)+min(y))/2)
      }
      else{
        LowerConfidenceTRAINValue[i]<-((PredictionMatrix[i,trunc((max_iter-burn_in+1)*0.05)]+PredictionMatrix[i,trunc(((max_iter-burn_in+1)*0.05)+1)])/2)*(max(y)-min(y))+((max(y)+min(y))/2)
        UpperConfidenceTRAINValue[i]<-((PredictionMatrix[i,trunc((max_iter-burn_in+1)*0.95)]+PredictionMatrix[i,trunc(((max_iter-burn_in+1)*0.95)+1)])/2)*(max(y)-min(y))+((max(y)+min(y))/2)
      }
    }
    
    
    
    LowerConfidenceTESTValue<-vector(length=length(mean_yhat_Test))
    UpperConfidenceTESTValue<-vector(length=length(mean_yhat_Test))
    
    for (i in 1:length(mean_yhat_Test)){
      TestMatrix[i,]<-sort(TestMatrix[i,])
      
      if ((((max_iter-burn_in+1)*0.05))== round((max_iter-burn_in+1)*0.05)){
        LowerConfidenceTESTValue[i]<-(TestMatrix[i,(max_iter-burn_in+1)*0.05])*(max(y)-min(y))+((max(y)+min(y))/2)
        UpperConfidenceTESTValue[i]<-(TestMatrix[i,(max_iter-burn_in+1)*0.95])*(max(y)-min(y))+((max(y)+min(y))/2)
      }
      else{
        LowerConfidenceTESTValue[i]<-((TestMatrix[i,trunc((max_iter-burn_in+1)*0.05)]+TestMatrix[i,trunc(((max_iter-burn_in+1)*0.05)+1)])/2)*(max(y)-min(y))+((max(y)+min(y))/2)
        UpperConfidenceTESTValue[i]<-((TestMatrix[i,trunc((max_iter-burn_in+1)*0.95)]+TestMatrix[i,trunc(((max_iter-burn_in+1)*0.95)+1)])/2)*(max(y)-min(y))+((max(y)+min(y))/2)
      }
    }
  
    plotCI(f(x),mean_yhat, UpperConfidenceTRAINValue-mean_yhat,mean_yhat-LowerConfidenceTRAINValue,sfrac=0, scol = 'grey',ylab = "posterior intervals", xlab = "In-Sample f(x)",  cex.lab = 1.5) 
    abline(0,1)
    
    plotCI(YTest,mean_yhat_Test,UpperConfidenceTESTValue-mean_yhat_Test,mean_yhat_Test-LowerConfidenceTESTValue,sfrac=0, scol = 'grey', ylab = "posterior intervals", xlab = "Out-of-Sample f(x)", cex.lab = 1.5)
    abline(0,1)
    
    in_interval_Train <- (f(x) >= LowerConfidenceTRAINValue) & (f(x) <= UpperConfidenceTRAINValue)
    in_interval_TEST<- (YTest >= LowerConfidenceTESTValue) & (YTest <= UpperConfidenceTESTValue)
    
    #plot(plotForSigmaSquared,col=c(rep('red',burn_in),rep('black',max_iter-burn_in)),xlab="MCMC iteration",ylab="Sigma draw",type="l")
    plot(plotForSigmaSquared, type = "n", col = "black", lwd = 2, xlab="MCMC iteration",ylab="Sigma draw", cex.lab = 1.5)
    
    # Draw the first segment in red
    segments(x0 = 1:(burn_in - 1), y0 = plotForSigmaSquared[1:(burn_in - 1)],
             x1 = 2:burn_in, y1 = plotForSigmaSquared[2:burn_in], col = "red", lwd = 2)
    # Draw the second segment in black
    segments(x0 = (burn_in):(max_iter - 1), y0 = plotForSigmaSquared[(burn_in+1):(max_iter - 1)],
             x1 = (burn_in+1):max_iter, y1 = plotForSigmaSquared[(burn_in):max_iter], col = "black", lwd = 2)
    abline(1,0)
    
    return(
      data.frame(
        invertal_Train = sum(in_interval_Train)/length(in_interval_Train),
        invertal_Test = sum(in_interval_TEST)/length(in_interval_TEST),
        RMSE = sqrt(mean((YTest-mean_yhat_Test)^2))
      )
    )
  }
  
  AddiVortes_Algorithm_Plot(Y[TrainSet],X[TrainSet,],200,2000,200,6,0.95,3,0.8,3,25,f(X[TestSet,]),X[TestSet,],IntialSigma = "Linear")

}

figure4<-function(){
  f = function(x){
    10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-0.5)^2+10*x[,4]+5*x[,5]
  }
  
  sigma = 1  #y = f(x) + sigma*z , z~N(0,1)
  n = 300      #number of observations
  set.seed(9)
  X=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
  Ey = f(X)
  Y=Ey+sigma*rnorm(n)
  
  n=length(Y)
  TrainSet=sort(sample.int(n,3*n/6))
  TestSet=1:n
  TestSet=TestSet[! TestSet %in% TrainSet]
  
  par(mfrow =c(1,2))
  
  AddiVortes_Algorithm_Plot<-function(y,x,m,max_iter,nu,q,k,sd,Omega,lambda_rate,IntialSigma = "Linear"){
    
    #Scaling x and y
    yScaled=(y-(max(y)+min(y))/2)/(max(y)-min(y))
    xScaled=x;
    for (i in 1:length(x[1,])){
      xScaled[,i]=(x[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
    }
    
    #Initialize:
    #Prediction Set (A list of vectors with the output values for each tessellation), 
    #Dimension set (A list of vectors with the covariates included in the tessellaions);
    #and Tessellation Set (A list of matrices that give the coordinates of the centers in the tessellations)
    
    Pred<-rep(list(matrix(mean(yScaled)/m)),m)
    Dim=vector(length = m)
    Tess=vector(length = m)
    for (i in 1:m){
      Dim[i]<-list(sample(1:length(x[1,]), 1))
      Tess[i]<-(list(matrix(rnorm(1,0,sd))))
    }
    
    #Prepare some variables used in the backfitting algorithm
    SumOfAllTess=rep(mean(yScaled),length(yScaled))
    SigmaSquaredMu=(0.5/(k*sqrt(m)))^2
    LastTessPred=matrix
    
    AverageNumberOfCells<-vector(length = max_iter)
    AverageNumberOfDim<-vector(length = max_iter)
    
    #finding lambda
    if (IntialSigma=="Naive"){ # Usually used if p is greater then n. Uses Standard deviation of y to predict sigma.
      SigmaSquaredHat=var(yScaled)
    }
    else{  # Default method using residual standard deviation from a least-squared linear regression of y, to predict sigma.
      MultiLinear<-lm(yScaled ~ x)
      SigmaSquaredHat=sum(MultiLinear$residuals^2)/(length(yScaled)-length(xScaled[1,])-1)
    }
    
    #Find lambda
    lambda=1;
    lambda <- optim(par = 1,
                    fitting_function,
                    method = "Brent",
                    lower = 0.001,
                    upper = 100,
                    q=q , nu=nu, sigmaSquared_hat=SigmaSquaredHat)$par
    
    for (i in 1:max_iter){
      NumOfCells<-0
      NumOfDim<-0
      
      #Sample Sigma squared using all tessellations to predict the outcome variables
      SigmaSquared=SigmaSquaredCalculation(yScaled,nu,lambda,SumOfAllTess)
      
      for (j in 1:m){
        
        
        NewTessOutput<-NewTess(xScaled,j,Tess,Dim,sd) #Propose new Tessellation 
        TessStar<-NewTessOutput[[1]]  
        DimStar<-NewTessOutput[[2]]
        Modification<-NewTessOutput[[3]]
        
        ResidualsOutput<-CalculateResiduals(yScaled,xScaled,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred) #Calculate the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation and the number of observations in each cell.
        R_ijOld<-ResidualsOutput[[1]]   #Old and New refer to the original and proposed tessellations
        n_ijOld<-ResidualsOutput[[2]]
        R_ijNew<-ResidualsOutput[[3]]
        n_ijNew<-ResidualsOutput[[4]]
        SumOfAllTess<-ResidualsOutput[[5]] #Keeps track of the prediction for all tessellations to help sample sigma squared.
        IndexesStar<-ResidualsOutput[[6]] #Gives the row of each observation for the cell it falls in for the proposed tessellation.
        Indexes<-ResidualsOutput[[7]]  #Gives the row of each observation for the cell it falls in for the original tessellation.
        
        if (!any(n_ijNew==0)){ #automatically rejects proposed tessellation if there exists a cell with no observations in.
          
          LOGAcceptenceProb=AlphaCalculation(xScaled,TessStar,DimStar,j,R_ijOld,n_ijOld,R_ijNew,n_ijNew,SigmaSquared,Modification,SigmaSquaredMu,Omega,lambda_rate) #Gives the log of the acceptence probability.
          
          
          if (log(runif(n=1, min=0, max=1))<LOGAcceptenceProb){ #Accepts the proposed tessellation is accepted then calculates the new output values for the new tessellation. 
            Tess=TessStar
            Dim=DimStar
            Pred[[j]]=NewPredSet(j,TessStar,R_ijNew,n_ijNew,SigmaSquaredMu,SigmaSquared)
            LastTessPred=Pred[[j]][IndexesStar]
          }
          else { #Rejects the proposed tesellation then calculates new output values for the original tessellation.
            Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
            LastTessPred=Pred[[j]][Indexes];
          }
        }
        else{ #Rejects the proposed tesellation then calculates new output values for the original tessellation.
          Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
          LastTessPred=Pred[[j]][Indexes];
        }
        if (j==m){ #If j equals m then adds the last tessellation output values to give a prediction.
          SumOfAllTess=SumOfAllTess+LastTessPred;
        }
        NumOfCells<-NumOfCells+length(Tess[[j]][,1])
        NumOfDim<-NumOfDim+length(Tess[[j]][1,])
      }
      
      AverageNumberOfCells[i]<-NumOfCells/m
      AverageNumberOfDim[i]<-NumOfDim/m
      if (i %% 100 == 0){
        cat(sprintf("Iteration %d out of %d", i, max_iter), "\n")
      }
    }
  
    
    plot(AverageNumberOfDim,type="l",xlab="MCMC iteration", ylab="Average Number Of Dimensions per Tessellaion")
    plot(AverageNumberOfCells,type="l",xlab="MCMC iteration", ylab="Average Number Of Centers per Tessellaion")
  
    
  }
  
  AddiVortes_Algorithm_Plot(Y[TrainSet],X[TrainSet,],200,2000,6,0.85,3,0.8,3,25,IntialSigma = "Linear")
} 

figure5<- function(max_iter = 6000 , burn_in = 1000){
  SigmaSquaredCalculation<-function(yScaled,nu,lambda,SumOfAllTess){ #Sample sigma squared from inverse gamma distribution
  
  n=length(yScaled)
  SigmaSquared<-rinvgamma(1,shape=(nu+n)/2,rate=(nu*lambda+sum((yScaled-SumOfAllTess)^2))/2)
  
  return(SigmaSquared)
}

NewTess<-function(x,j,Tess,Dim,sd){ #Propose a new tessellation
  
  p=runif(1,0,1) #Randomly sample p to decide the proposed modification to tessellation.
  
  DimStar=Dim # Let proposed dimension matrix equal original dimension matrix.
  TessStar=Tess #Similar for the tessellation matrix.
  
  if (p<0.2 & length(Dim[[j]])!=length(x[1,]) | length(Dim[[j]])==1 & p<0.4){ #Add a dimension if p is less then 0.2 or if p is less then 0.4 when there is only one dimension in the Tessellation due to adjustments (Supplementary Material).
    NumberOfCovariates=1:length(x[1,]) #Let NumberOfCovariates be a vector from 1 to the number of covariates considered.
    NumberOfCovariates=NumberOfCovariates[-Dim[[j]]] #Remove all values that for the covariates that are already in the tessellation.
    DimStar[[j]]<-c(Dim[[j]],sample(NumberOfCovariates,1)) # Uniformly sample a new covariate and add it to the dimension matrix.
    TessStar[[j]]=cbind(Tess[[j]],rnorm(length(Tess[[j]][,1]),0,sd)) # Sample new coordinates from Normal distribution for the new dimension and add it to the Tessellation matrix.
    Modification="AD"}
  else if (p<0.4){ #Remove a dimension if p is less then 0.4.
    RemovedDim=sample(1:length(Dim[[j]]),1) #Uniformly sample the dimension to be removed.
    DimStar[[j]]=DimStar[[j]][-RemovedDim] #Remove the dimension from the dimesion Matrix.
    TessStar[[j]]=matrix(TessStar[[j]][,-RemovedDim],ncol=length(DimStar[[j]])) #Remove the coordinates in the Tessellation matrix corresponding to the dimension removed.
    Modification="RD"}
  else if (p<0.6 || p<0.8 & length(Tess[[j]][,1])==1){ #Add a centre if p is less then 0.6 or if p is less then 0.4 when there is only one center in the Tessellation due to adjustments (Supplementary Material).
    TessStar[[j]]=rbind(Tess[[j]],rnorm(length(Dim[[j]]),0,sd)) #Add a new row of coordinates, sampled from a normal distribution, to the Tessellation matrix to add a center.
    Modification="AC"}
  else if (p<0.8){ #Add a centre if p is less then 0.8. 
    CenterRemoved=sample(1:length(TessStar[[j]][,1]),1) #Sample a row.
    TessStar[[j]]=matrix(TessStar[[j]][-CenterRemoved,],ncol=length(Dim[[j]])) #Remove row sampled.
    Modification="RC"}
  else if (p<0.9 || length(Dim[[j]])==length(x[1,])){ #Change a center if p is less then 0.9 or if the all the covariates are in the tessellation.
    TessStar[[j]][sample(1:length(TessStar[[j]][,1]),1),]=rnorm(length(Dim[[j]]),0,sd) # Sample a row in the tessellaion matrix and change the coordinates of the centre by sampling from a normal distribution.
    Modification="Change"}
  else{ #Swop a dimension.
    NumberOfCovariates=1:length(x[1,])  #Let NumberOfCovariates be a vector from 1 to the number of covariates considered.
    NumberOfCovariates=NumberOfCovariates[-Dim[[j]]]  #Remove all values that for the covariates that are already in the tessellation.
    DimToChange=sample(1:length(Dim[[j]]),1) #Uniformly sample a dimension to change.
    DimStar[[j]][DimToChange]=sample(NumberOfCovariates,1) #Replace the Dimension to a new uniforly sampled covariate that is not already in the tessellaion.
    TessStar[[j]][,DimToChange]=rnorm(length(Tess[[j]][,1]),0,sd) #Add new normally sampled coordinates new dimension added.
    Modification="Swop"}
  
  TessStar[[j]]<-matrix(TessStar[[j]],ncol=length(DimStar[[j]])) #Ensure the the Tessellation matrix is a "matrix" type.
  
  return(list(TessStar,DimStar,Modification)) #Return new proposed tessellation.
}

fitting_function<- function(lambda,q,nu,sigmaSquared_hat){ #function that calculates the squared difference between sigma squared hat and the inverse gamma function
  return((sigmaSquared_hat- qinvgamma(q, shape=nu/2, rate=nu*lambda/2))^2)
}

Indexes<-function(x,Tess,Dim){ #Gives the row (the center) of the tessellation that each obseravtion falls within.
  if (length(Tess[,1])==1){ #only 1 centre
    CellsForGivenTess=rep(1,length(x[,1]))
  }
  else{ #multiple
    CellsForGivenTess=knnx.index(Tess,matrix(x[,Dim],ncol = length(Dim)),1)
  }
  return(CellsForGivenTess)
}

AlphaCalculation<-function(x,Tess,Dim,j,R_ijOld,n_ijOld,R_ijNew,n_ijNew,SigmaSquared,Modification,SigmaSquaredMu,Omega,lambda_rate){ #Calculates the acceptence rate of the proposed tessellation.
  
  d=length(Dim[[j]]);
  NumCovariates=length(x[1,]);
  cStar=length(Tess[[j]][,1]);

  #The Log Likelihood Ratio in the acceptence ratio
  LOGlikelihoodRatio=0.5*(log(prod(n_ijOld*SigmaSquaredMu+SigmaSquared))-log(prod(n_ijNew*SigmaSquaredMu+SigmaSquared)))+((SigmaSquaredMu/(2*SigmaSquared))*(-sum((R_ijOld^2)/(n_ijOld*SigmaSquaredMu+SigmaSquared))+sum((R_ijNew^2)/(n_ijNew*SigmaSquaredMu+SigmaSquared))))

  #Calculating the acceptence probablity for "AD"=Adding a dimension, "RD"=Removing a dimension, "AC"=Adding a center, "RC"=Removing a center, "Change"=Changing the coordinates of a center and Swopping a dimension.
  if (Modification == "AD"){ 
    TessStructure=(dbinom(d-1,NumCovariates-1,Omega/NumCovariates))/(dbinom(d-2,NumCovariates-1,Omega/NumCovariates)*(NumCovariates-d+1))
    TransitionRatio=(NumCovariates-d+1)/d;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
    
    #Adjustments.
    if (length(Dim[[j]])==1){
      AcceptenceProb=AcceptenceProb+log(1/2)
    }
    else if (length(Dim[[j]])==NumCovariates-1){
      AcceptenceProb=AcceptenceProb+log(2)}
  }
  else if (Modification == "RD"){
    TessStructure=(dbinom(d-1,NumCovariates,Omega/NumCovariates)*(NumCovariates-d))/(dbinom(d,NumCovariates,Omega/NumCovariates))
    TransitionRatio=(d+1)/(NumCovariates-d)
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
    
    #Adjustments.
    if (length(Dim[[j]])==NumCovariates){
      AcceptenceProb=AcceptenceProb+log(1/2)
    }
    else if (length(Dim[[j]])==2){
      AcceptenceProb=AcceptenceProb+log(2);
    }
  }
  else if(Modification == "AC"){
    TessStructure=dpois(cStar-1,lambda_rate)/dpois(cStar-2,lambda_rate)
    TransitionRatio=1/cStar;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)+0.5*log(SigmaSquared)

    #Adjustments.
    if (cStar==1){
      AcceptenceProb=AcceptenceProb+log(1/2);
    }
  }
  else if (Modification == "RC"){
    TessStructure=dpois(cStar-1,lambda_rate)/dpois(cStar,lambda_rate);
    TransitionRatio=cStar+1;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)-0.5*log(SigmaSquared)

    #Adjustments.
    if (cStar==2){
      AcceptenceProb=AcceptenceProb+log(2);
    }
  }
  else if (Modification == "Change"){
    TessStructure=1;
    TransitionRatio=1;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio);
  }
  else {
    TessStructure=1;
    TransitionRatio=1;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio);
  }
  return(AcceptenceProb)
}

CalculateResiduals<-function(y,x,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred){ #A function that calculates the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation.
  if (j==1){
    indexes=Indexes(x,Tess[[j]],Dim[[j]]);
    CurrentTessPred<-Pred[[j]][indexes]
    SumOfAllTess=SumOfAllTess-CurrentTessPred}
  else{
    indexes=Indexes(x,Tess[[j]],Dim[[j]]);
    CurrentTessPred<-Pred[[j]][indexes]
    SumOfAllTess=SumOfAllTess+LastTessPred-CurrentTessPred;
  }
  
  IndexesStar=Indexes(x,TessStar[[j]],DimStar[[j]]);
  R_j<-y-SumOfAllTess
  
  #Initializing Sizes
  
  R_ijOld=rep(0,length(Pred[[j]]))
  n_ijOld=rep(0,length(Pred[[j]]))
  
  for (i in 1:length(Pred[[j]])){
    R_ijOld[i]<-sum(R_j[indexes==i])
    n_ijOld[i]<-sum(indexes==i)
  }
  
  R_ijNew=rep(0,length(TessStar[[j]][,1]))
  n_ijNew=rep(0,length(TessStar[[j]][,1]))
  
  for (i in 1:length(TessStar[[j]][,1])){
    R_ijNew[i]<-sum(R_j[IndexesStar==i])
    n_ijNew[i]<-sum(IndexesStar==i)}
  
  return(list(R_ijOld,n_ijOld,R_ijNew,n_ijNew,SumOfAllTess,IndexesStar,indexes))
}

NewPredSet<-function(j,Tess,R_ijNew,n_ijNew,sigmaSquaredMu,SigmaSquared){ #Sampling the new output values for the new tessellation.
  PredSet=rep(0,length(Tess[[j]][,1]))
  for (i in 1:length(Tess[[j]][,1])){
    PredSet[i]=rnorm(1,sigmaSquaredMu*R_ijNew[i]/(sigmaSquaredMu*n_ijNew[i]+SigmaSquared),((SigmaSquared*sigmaSquaredMu)/(n_ijNew[i]*sigmaSquaredMu+SigmaSquared))^0.5);
  }
  return(PredSet)
}

TestPrediction<-function(x,m,Tess,Dim,Pred){ #A function that derives a prediction for the output variable give an observation, for one iteration in AddiVortes.
  Prediction=rep(0,length(x[,1]));
  for (j in 1:m){
    NewTessIndexes=Indexes(x,Tess[[j]],Dim[[j]]);
    Prediction=Prediction+Pred[[j]][NewTessIndexes]
  }
  return(Prediction)
}
  
  AddiVortes_Algorithm<-function(y,x,m,max_iter,burn_in,nu,q,k,sd,Omega,lambda_rate,YTest,XTest,IntialSigma = "Linear"){
  
    #Scaling x and y
    yScaled=(y-(max(y)+min(y))/2)/(max(y)-min(y))
    xScaled=x;
    for (i in 1:length(x[1,])){
      xScaled[,i]=(x[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
    }
    
    
    for (i in 1:length(XTest[1,])){
      XTest[,i]=(XTest[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
    }
    
    #Initialize:
    #Prediction Set (A list of vectors with the output values for each tessellation), 
    #Dimension set (A list of vectors with the covariates included in the tessellaions);
    #and Tessellation Set (A list of matrices that give the coordinates of the centers in the tessellations)
    
    Pred<-rep(list(matrix(mean(yScaled)/m)),m)
    Dim=vector(length = m)
    Tess=vector(length = m)
    for (i in 1:m){
      Dim[i]<-list(sample(1:length(x[1,]), 1))
      Tess[i]<-(list(matrix(rnorm(1,0,sd))))
    }
    
    #Prepare some variables used in the backfitting algorithm
    SumOfAllTess=rep(mean(yScaled),length(yScaled))
    SigmaSquaredMu=(0.5/(k*sqrt(m)))^2
    LastTessPred=matrix
    
    #Matrices that will hold the samples from the poseterior distribution for the training samples and test samples.
    PredictionMatrix<-array(dim=c(length(y),(max_iter-burn_in)))
    TestMatrix<-array(dim=c(length(YTest),(max_iter-burn_in)))
    
    #finding lambda
    if (IntialSigma=="Naive"){ # Usually used if p is greater then n. Uses Standard deviation of y to predict sigma.
      SigmaSquaredHat=var(yScaled)
    }
    else{  # Default method using residual standard deviation from a least-squared linear regression of y, to predict sigma.
      MultiLinear<-lm(yScaled ~ xScaled)
      SigmaSquaredHat=sum(MultiLinear$residuals^2)/(length(yScaled)-length(xScaled[1,])-1)
    }
    
    #Find lambda
    lambda=1;
    lambda <- optim(par = 1,
                    fitting_function,
                    method = "Brent",
                    lower = 0.001,
                    upper = 100,
                    q=q , nu=nu, sigmaSquared_hat=SigmaSquaredHat)$par
    
    for (i in 1:max_iter){
      
      #Sample Sigma squared using all tessellations to predict the outcome variables
      SigmaSquared=SigmaSquaredCalculation(yScaled,nu,lambda,SumOfAllTess)
      
      for (j in 1:m){
        NewTessOutput<-NewTess(xScaled,j,Tess,Dim,sd) #Propose new Tessellation 
        TessStar<-NewTessOutput[[1]]  
        DimStar<-NewTessOutput[[2]]
        Modification<-NewTessOutput[[3]]
        
        ResidualsOutput<-CalculateResiduals(yScaled,xScaled,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred) #Calculate the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation and the number of observations in each cell.
        R_ijOld<-ResidualsOutput[[1]]   #Old and New refer to the original and proposed tessellations
        n_ijOld<-ResidualsOutput[[2]]
        R_ijNew<-ResidualsOutput[[3]]
        n_ijNew<-ResidualsOutput[[4]]
        SumOfAllTess<-ResidualsOutput[[5]] #Keeps track of the prediction for all tessellations to help sample sigma squared.
        IndexesStar<-ResidualsOutput[[6]] #Gives the row of each observation for the cell it falls in for the proposed tessellation.
        Indexes<-ResidualsOutput[[7]]  #Gives the row of each observation for the cell it falls in for the original tessellation.
        
        if (!any(n_ijNew==0)){ #automatically rejects proposed tessellation if there exists a cell with no observations in.
          
          LOGAcceptenceProb=AlphaCalculation(xScaled,TessStar,DimStar,j,R_ijOld,n_ijOld,R_ijNew,n_ijNew,SigmaSquared,Modification,SigmaSquaredMu,Omega,lambda_rate) #Gives the log of the acceptence probability.
          
          if (log(runif(n=1, min=0, max=1))<LOGAcceptenceProb){ #Accepts the proposed tessellation is accepted then calculates the new output values for the new tessellation. 
            Tess=TessStar
            Dim=DimStar
            Pred[[j]]=NewPredSet(j,TessStar,R_ijNew,n_ijNew,SigmaSquaredMu,SigmaSquared)
            LastTessPred=Pred[[j]][IndexesStar]
          }
          else { #Rejects the proposed tesellation then calculates new output values for the original tessellation.
            Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
            LastTessPred=Pred[[j]][Indexes];
          }
        }
        else{ #Rejects the proposed tesellation then calculates new output values for the original tessellation.
          Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
          LastTessPred=Pred[[j]][Indexes];
        }
        if (j==m){ #If j equals m then adds the last tessellation output values to give a prediction.
          SumOfAllTess=SumOfAllTess+LastTessPred;
        }
      }
      
      if (i>burn_in){ #vectors that hold the predictions for each iteration after burn in.
        PredictionMatrix[,(i-burn_in)]=SumOfAllTess;
        TestMatrix[,(i-burn_in)]=TestPrediction(XTest,m,Tess,Dim,Pred);
      }
    }
    
    #finding the mean of the predition over the iterations and then unscaling the predictions.
    mean_yhat=(rowSums(PredictionMatrix)/(max_iter-burn_in))*(max(y)-min(y))+((max(y)+min(y))/2)
    mean_yhat_Test=(rowSums(TestMatrix)/(max_iter-burn_in))*(max(y)-min(y))+((max(y)+min(y))/2)
    
    return( #Returns the RMSE value for the test samples.
      data.frame(
        mean_yhat_Test
      )
    )
  }
  
  f = function(x){
    10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-0.5)^2+10*x[,4]+5*x[,5]
  }
  
  sigma = 1  #y = f(x) + sigma*z , z~N(0,1)
  n = 100      #number of observations
  set.seed(9)
  X=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
  Ey = f(X)
  Y=Ey+sigma*rnorm(n)
  
  n=length(Y)
  TrainSet=sort(sample.int(n,3*n/6))
  TestSet=1:n
  TestSet=TestSet[! TestSet %in% TrainSet]
  
  
  num_cores <- 11  # Specify the number of cores you want to use
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  xlabels<-c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10")
  
  par(mfrow = c(2, 5))
  
  # Assuming you have a dataset 'data' and a black-box regression model 'model'
  for (i in 1:10){
    library('invgamma')
    library('FNN')
    library('expint')
    
    # Create a data frame with X1 values and other variables the same as your dataset
    data_range <-rep(list(X),n)
    
    predictions<-vector(length=11)
    
    # Make predictions on data_range using the black-box regression model
    
    predictions<- foreach(j = c(1,10,20,30,40,50,60,70, 80,90,100),.combine = rbind) %dopar% {
      library('invgamma')
      library('FNN')
      library('expint')
      library("plotrix")
      
      h<-sort(X[,i])
      data_range[[j]][,i] <- rep(h[j],n)
      return(
        mean(AddiVortes_Algorithm(Y,X,200,max_iter,burn_in,6,0.85,3,0.8,3,25,Y,data_range[[j]])$mean_yhat_Test)
      )
    }
    
    h<-sort(X[,i])
    # Calculate the partial dependence for X1 by averaging predictions
    plot(h[c(1,10,20,30,40,50,60,70, 80,90,100)],predictions,type='b',ylim=c(9,18), xlab =xlabels[i],ylab = "")
    if (i == 1 | i==6) {
      mtext("Partial dependence",  side = 2, line = 2.5, )  # Adjust line parameter as needed
    }
    
  }
  stopCluster(cl)
}

figure6<-function(max_iter = 6000, burn_in = 1000){

AddiVortes_Algorithm_Plot_figure6<-function(y,x,m,max_iter,burn_in,nu,q,k,sd,Omega,lambda_rate,YTest,XTest,IntialSigma = "Linear"){

  p=length(x[1,]) #p is the number of covariates used
  
  #Scaling x and y
  yScaled=(y-(max(y)+min(y))/2)/(max(y)-min(y))
  xScaled=x;
  for (i in 1:length(x[1,])){
    xScaled[,i]=(x[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
  }

  for (i in 1:length(XTest[1,])){
    XTest[,i]=(XTest[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
  }
  
  #Initialize the Prediction Set, Dimension set and Tessellation Set
  
  Pred<-rep(list(matrix(mean(yScaled)/m)),m)
  
  Dim=vector(length = m)
  Tess=vector(length = m)
  for (i in 1:m){
    Dim[i]<-list(sample(1:length(x[1,]), 1))
    Tess[i]<-(list(matrix(rnorm(1,0,sd))))}
  
  SumOfAllTess=rep(mean(yScaled),length(yScaled));
  SigmaSquaredMu=(0.5/(k*sqrt(m)))^2;
  LastTessPred=matrix
  acceptedTess=0;
  PredictionMatrix<-array(dim=c(length(y),(max_iter-burn_in)))
  TestMatrix<-array(dim=c(length(YTest),(max_iter-burn_in)))
  CovariatesUsed<-rep(0,length(x[1,]))
  plotForSigmaSquared<-vector(length = max_iter)
  plotForRMSE<-vector(length = max_iter)
  AverageNumberOfCells<-vector(length = max_iter)
  AverageNumberOfDim<-vector(length = max_iter)
  
  #finding lambda
  if (IntialSigma=="Naive"){
    SigmaSquaredHat=var(yScaled)
  }
  else{
    MultiLinear<-lm(yScaled ~ xScaled)
    SigmaSquaredHat=sum(MultiLinear$residuals^2)/(length(yScaled)-length(xScaled[1,])-1)
  }
  lambda=1;
  lambda <- optim(par = 1,
                  fitting_function,
                  method = "Brent",
                  lower = 0.001,
                  upper = 100,
                  q=q , nu=nu, sigmaSquared_hat=SigmaSquaredHat)$par
  
  for (i in 1:max_iter){
    NumOfCells<-0
    NumOfDim<-0
    
    #Sample Whole model Sigma squared
    SigmaSquared=SigmaSquaredCalculation(yScaled,nu,lambda,SumOfAllTess)
    plotForSigmaSquared[i]=(SigmaSquared*(max(y)-min(y))^2)^0.5
    plotForRMSE[i]=(mean((yScaled-SumOfAllTess)^2))^0.5
    
    for (j in 1:m){
      NewTessOutput<-NewTess(xScaled,j,Tess,Dim,sd)
      TessStar<-NewTessOutput[[1]]
      DimStar<-NewTessOutput[[2]]
      Modification<-NewTessOutput[[3]]
      
      ResidualsOutput<-CalculateResiduals(yScaled,xScaled,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred)
      R_ijOld<-ResidualsOutput[[1]]
      n_ijOld<-ResidualsOutput[[2]]
      R_ijNew<-ResidualsOutput[[3]]
      n_ijNew<-ResidualsOutput[[4]]
      SumOfAllTess<-ResidualsOutput[[5]]
      IndexesStar<-ResidualsOutput[[6]]
      Indexes<-ResidualsOutput[[7]]
      
      if (!any(n_ijNew==0)){
        
        LOGAcceptenceProb=AlphaCalculation(xScaled,TessStar,DimStar,j,R_ijOld,n_ijOld,R_ijNew,n_ijNew,SigmaSquared,Modification,SigmaSquaredMu,Omega,lambda_rate);
        
        if (log(runif(n=1, min=0, max=1))<LOGAcceptenceProb){
          Tess=TessStar
          Dim=DimStar
          Pred[[j]]=NewPredSet(j,TessStar,R_ijNew,n_ijNew,SigmaSquaredMu,SigmaSquared)
          LastTessPred=Pred[[j]][IndexesStar]
          acceptedTess=acceptedTess+1
        }
        else {
          Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
          LastTessPred=Pred[[j]][Indexes];
        }
      }
      else{
        Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
        LastTessPred=Pred[[j]][Indexes];
      }
      if (j==m){
        SumOfAllTess=SumOfAllTess+LastTessPred;
      }
      if(i>burn_in){
        CovariatesUsed[Dim[[j]]]<-CovariatesUsed[Dim[[j]]]+1
      }
      NumOfCells<-NumOfCells+length(Tess[[j]][,1])
      NumOfDim<-NumOfDim+length(Tess[[j]][1,])
    }
    if (i %% 100 == 0){
      cat(sprintf("Iteration %d out of %d", i, max_iter), "\n")
    }
    
    AverageNumberOfCells[i]<-NumOfCells/m
    AverageNumberOfDim[i]<-NumOfDim/m
    
    if (i>burn_in){
      PredictionMatrix[,(i-burn_in)]=SumOfAllTess;
      TestMatrix[,(i-burn_in)]=TestPrediction(XTest,m,Tess,Dim,Pred);
    }
  }
  
  mean_yhat=(rowSums(PredictionMatrix)/(max_iter-burn_in))*(max(y)-min(y))+((max(y)+min(y))/2)
  mean_yhat_Test=(rowSums(TestMatrix)/(max_iter-burn_in))*(max(y)-min(y))+((max(y)+min(y))/2)

  CovariataesUsedPercentage<-CovariatesUsed/sum(CovariatesUsed)
  
  return(
    data.frame(
      CovariataesUsedPercentage
    )
  )
}


  
  f = function(x){
  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-0.5)^2+10*x[,4]+5*x[,5]
  }

  sigma = 1  #y = f(x) + sigma*z , z~N(0,1)
  n = 510      #number of observations
  set.seed(23)
  X=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
  Ey = f(X)
  Y=Ey+sigma*rnorm(n)

  n=length(Y)
  TrainSet=sort(sample.int(n,50*n/51));
  TestSet=1:n;
  TestSet=TestSet[! TestSet %in% TrainSet];
  par(mfrow=c(1,1))
  #Creating graph
  
  lineCol <- c("#e23224", "orange", "#762d92", "black", "#0928ff")
  lineType <- c(1, 2, 3, 4, 5)  # Different line types for each line
  labels <- c('m=10', 'm=20', 'm=50', 'm=100', 'm=200')  # Labels for each line
  i <- 1
    
  for (m in c(10, 20, 50, 100, 200)) {
     if (m == 10) {
      plot(AddiVortes_Algorithm_Plot_figure6(Y[TrainSet], as.matrix(X[TrainSet, ]), m, max_iter, burn_in, 6, 0.85, 3, 0.8, 3, 25, f(X[TestSet, ]), as.matrix(X[TestSet, ]))$CovariataesUsedPercentage, type = "b", col = lineCol[1], lty = lineType[1], xlim = c(1, 10),ylim=c(0,0.3), ylab = "Percentage Used",xlab="Covariate",cex.lab=1.5,lwd=2)
    } 
    else {
      lines(AddiVortes_Algorithm_Plot_figure6(Y[TrainSet], as.matrix(X[TrainSet, ]), m, max_iter, burn_in, 6, 0.85, 3, 0.8, 3, 25, f(X[TestSet, ]), as.matrix(X[TestSet, ]))$CovariataesUsedPercentage, col = lineCol[i],lty = lineType[i] , type="b",lwd=2)
    }
      
    # Add custom x-axis labels at breaks from 1 to 10
    axis(1, at = 1:10)
      i <- i + 1
    }

    
  # Add a legend to the plot
  legend("topright", legend = labels, col = lineCol, lty = lineType,lwd=2,inset = 0.1)
}

figure7<-function(max_iter = 6000, burn_in = 1000, num_of_datasets = 100){
  
  SigmaSquaredCalculation<-function(yScaled,nu,lambda,SumOfAllTess){ #Sample sigma squared from inverse gamma distribution
    
    n=length(yScaled)
    SigmaSquared<-rinvgamma(1,shape=(nu+n)/2,rate=(nu*lambda+sum((yScaled-SumOfAllTess)^2))/2)
    
    return(SigmaSquared)
  }
  
  NewTess<-function(x,j,Tess,Dim,sd){ #Propose a new tessellation
    
    p=runif(1,0,1) #Randomly sample p to decide the proposed modification to tessellation.
    
    DimStar=Dim # Let proposed dimension matrix equal original dimension matrix.
    TessStar=Tess #Similar for the tessellation matrix.
    
    if (p<0.2 & length(Dim[[j]])!=length(x[1,]) | length(Dim[[j]])==1 & p<0.4){ #Add a dimension if p is less then 0.2 or if p is less then 0.4 when there is only one dimension in the Tessellation due to adjustments (Supplementary Material).
      NumberOfCovariates=1:length(x[1,]) #Let NumberOfCovariates be a vector from 1 to the number of covariates considered.
      NumberOfCovariates=NumberOfCovariates[-Dim[[j]]] #Remove all values that for the covariates that are already in the tessellation.
      DimStar[[j]]<-c(Dim[[j]],sample(NumberOfCovariates,1)) # Uniformly sample a new covariate and add it to the dimension matrix.
      TessStar[[j]]=cbind(Tess[[j]],rnorm(length(Tess[[j]][,1]),0,sd)) # Sample new coordinates from Normal distribution for the new dimension and add it to the Tessellation matrix.
      Modification="AD"}
    else if (p<0.4){ #Remove a dimension if p is less then 0.4.
      RemovedDim=sample(1:length(Dim[[j]]),1) #Uniformly sample the dimension to be removed.
      DimStar[[j]]=DimStar[[j]][-RemovedDim] #Remove the dimension from the dimesion Matrix.
      TessStar[[j]]=matrix(TessStar[[j]][,-RemovedDim],ncol=length(DimStar[[j]])) #Remove the coordinates in the Tessellation matrix corresponding to the dimension removed.
      Modification="RD"}
    else if (p<0.6 || p<0.8 & length(Tess[[j]][,1])==1){ #Add a centre if p is less then 0.6 or if p is less then 0.4 when there is only one center in the Tessellation due to adjustments (Supplementary Material).
      TessStar[[j]]=rbind(Tess[[j]],rnorm(length(Dim[[j]]),0,sd)) #Add a new row of coordinates, sampled from a normal distribution, to the Tessellation matrix to add a center.
      Modification="AC"}
    else if (p<0.8){ #Add a centre if p is less then 0.8. 
      CenterRemoved=sample(1:length(TessStar[[j]][,1]),1) #Sample a row.
      TessStar[[j]]=matrix(TessStar[[j]][-CenterRemoved,],ncol=length(Dim[[j]])) #Remove row sampled.
      Modification="RC"}
    else if (p<0.9 || length(Dim[[j]])==length(x[1,])){ #Change a center if p is less then 0.9 or if the all the covariates are in the tessellation.
      TessStar[[j]][sample(1:length(TessStar[[j]][,1]),1),]=rnorm(length(Dim[[j]]),0,sd) # Sample a row in the tessellaion matrix and change the coordinates of the centre by sampling from a normal distribution.
      Modification="Change"}
    else{ #Swop a dimension.
      NumberOfCovariates=1:length(x[1,])  #Let NumberOfCovariates be a vector from 1 to the number of covariates considered.
      NumberOfCovariates=NumberOfCovariates[-Dim[[j]]]  #Remove all values that for the covariates that are already in the tessellation.
      DimToChange=sample(1:length(Dim[[j]]),1) #Uniformly sample a dimension to change.
      DimStar[[j]][DimToChange]=sample(NumberOfCovariates,1) #Replace the Dimension to a new uniforly sampled covariate that is not already in the tessellaion.
      TessStar[[j]][,DimToChange]=rnorm(length(Tess[[j]][,1]),0,sd) #Add new normally sampled coordinates new dimension added.
      Modification="Swop"}
    
    TessStar[[j]]<-matrix(TessStar[[j]],ncol=length(DimStar[[j]])) #Ensure the the Tessellation matrix is a "matrix" type.
    
    return(list(TessStar,DimStar,Modification)) #Return new proposed tessellation.
  }
  
  fitting_function<- function(lambda,q,nu,sigmaSquared_hat){ #function that calculates the squared difference between sigma squared hat and the inverse gamma function
    return((sigmaSquared_hat- qinvgamma(q, shape=nu/2, rate=nu*lambda/2))^2)
  }
  
  Indexes<-function(x,Tess,Dim){ #Gives the row (the center) of the tessellation that each obseravtion falls within.
    if (length(Tess[,1])==1){ #only 1 centre
      CellsForGivenTess=rep(1,length(x[,1]))
    }
    else{ #multiple
      CellsForGivenTess=knnx.index(Tess,matrix(x[,Dim],ncol = length(Dim)),k=1)
    }
    return(CellsForGivenTess)
  }
  
  AlphaCalculation<-function(x,Tess,Dim,j,R_ijOld,n_ijOld,R_ijNew,n_ijNew,SigmaSquared,Modification,SigmaSquaredMu,Omega,lambda_rate){ #Calculates the acceptence rate of the proposed tessellation.
    
    d=length(Dim[[j]]);
    NumCovariates=length(x[1,]);
    cStar=length(Tess[[j]][,1]);
    #The Log Likelihood Ratio in the acceptance ratio
    LOGlikelihoodRatio=0.5*(log(prod(n_ijOld*SigmaSquaredMu+SigmaSquared))-log(prod(n_ijNew*SigmaSquaredMu+SigmaSquared)))+((SigmaSquaredMu/(2*SigmaSquared))*(-sum((R_ijOld^2)/(n_ijOld*SigmaSquaredMu+SigmaSquared))+sum((R_ijNew^2)/(n_ijNew*SigmaSquaredMu+SigmaSquared))))
    #Calculating the acceptance probability for "AD"=Adding a dimension, "RD"=Removing a dimension, "AC"=Adding a center, "RC"=Removing a center, "Change"=Changing the coordinates of a center and Swopping a dimension.
    if (Modification == "AD"){ 
      TessStructure=(dbinom(d-1,NumCovariates-1,Omega/NumCovariates))/(dbinom(d-2,NumCovariates-1,Omega/NumCovariates)*(NumCovariates-d+1))
      TransitionRatio=(NumCovariates-d+1)/d;
      AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
      #Adjustments.
      if (length(Dim[[j]])==1){
        AcceptenceProb=AcceptenceProb+log(1/2)
      }
      else if (length(Dim[[j]])==NumCovariates-1){
        AcceptenceProb=AcceptenceProb+log(2)}
    }
    else if (Modification == "RD"){
      TessStructure=(dbinom(d-1,NumCovariates,Omega/NumCovariates)*(NumCovariates-d))/(dbinom(d,NumCovariates,Omega/NumCovariates))
      TransitionRatio=(d+1)/(NumCovariates-d)
      AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
      
      #Adjustments.
      if (length(Dim[[j]])==NumCovariates){
        AcceptenceProb=AcceptenceProb+log(1/2)
      }
      else if (length(Dim[[j]])==2){
        AcceptenceProb=AcceptenceProb+log(2);
      }
    }
    else if(Modification == "AC"){
      TessStructure=dpois(cStar-1,lambda_rate)/dpois(cStar-2,lambda_rate)
      TransitionRatio=1/cStar;
      AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)+0.5*log(SigmaSquared)
      
      #Adjustments.
      if (cStar==1){
        AcceptenceProb=AcceptenceProb+log(1/2);
      }
    }
    else if (Modification == "RC"){
      TessStructure=dpois(cStar-1,lambda_rate)/dpois(cStar,lambda_rate);
      TransitionRatio=cStar+1;
      AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)-0.5*log(SigmaSquared)
      
      #Adjustments.
      if (cStar==2){
        AcceptenceProb=AcceptenceProb+log(2);
      }
    }
    else if (Modification == "Change"){
      TessStructure=1;
      TransitionRatio=1;
      AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio);
    }
    else {
      TessStructure=1;
      TransitionRatio=1;
      AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio);
    }
    return(AcceptenceProb)
  }
  
  CalculateResiduals<-function(y,x,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred){ #A function that calculates the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation.
    if (j==1){
      indexes=Indexes(x,Tess[[j]],Dim[[j]]);
      CurrentTessPred<-Pred[[j]][indexes]
      SumOfAllTess=SumOfAllTess-CurrentTessPred}
    else{
      indexes=Indexes(x,Tess[[j]],Dim[[j]]);
      CurrentTessPred<-Pred[[j]][indexes]
      SumOfAllTess=SumOfAllTess+LastTessPred-CurrentTessPred;
    }
    
    IndexesStar=Indexes(x,TessStar[[j]],DimStar[[j]]);
    R_j<-y-SumOfAllTess
    
    #Initializing Sizes
    
    R_ijOld=rep(0,length(Pred[[j]]))
    n_ijOld=rep(0,length(Pred[[j]]))
    
    for (i in 1:length(Pred[[j]])){
      R_ijOld[i]<-sum(R_j[indexes==i])
      n_ijOld[i]<-sum(indexes==i)
    }
    
    R_ijNew=rep(0,length(TessStar[[j]][,1]))
    n_ijNew=rep(0,length(TessStar[[j]][,1]))
    
    for (i in 1:length(TessStar[[j]][,1])){
      R_ijNew[i]<-sum(R_j[IndexesStar==i])
      n_ijNew[i]<-sum(IndexesStar==i)}
    
    return(list(R_ijOld,n_ijOld,R_ijNew,n_ijNew,SumOfAllTess,IndexesStar,indexes))
  }
  
  NewPredSet<-function(j,Tess,R_ijNew,n_ijNew,sigmaSquaredMu,SigmaSquared){ #Sampling the new output values for the new tessellation.
    PredSet=rep(0,length(Tess[[j]][,1]))
    for (i in 1:length(Tess[[j]][,1])){
      PredSet[i]=rnorm(1,sigmaSquaredMu*R_ijNew[i]/(sigmaSquaredMu*n_ijNew[i]+SigmaSquared),((SigmaSquared*sigmaSquaredMu)/(n_ijNew[i]*sigmaSquaredMu+SigmaSquared))^0.5);
    }
    return(PredSet)
  }
  
  TestPrediction<-function(x,m,Tess,Dim,Pred){ #A function that derives a prediction for the output variable give an observation, for one iteration in AddiVortes.
    Prediction=rep(0,length(x[,1]));
    for (j in 1:m){
      NewTessIndexes=Indexes(x,Tess[[j]],Dim[[j]]);
      Prediction=Prediction+Pred[[j]][NewTessIndexes]
    }
    return(Prediction)
  }
  
  AddiVortes_Algorithm<-function(y,x,m,max_iter,burn_in,nu,q,k,sd,Omega,lambda_rate,YTest,XTest,IntialSigma = "Linear"){
  
  #Scaling x and y
  yScaled=(y-(max(y)+min(y))/2)/(max(y)-min(y))
  xScaled=x;
  for (i in 1:length(x[1,])){
    xScaled[,i]=(x[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
  }
  
  
  for (i in 1:length(XTest[1,])){
    XTest[,i]=(XTest[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
  }
  
  #Initialize:
  #Prediction Set (A list of vectors with the output values for each tessellation), 
  #Dimension set (A list of vectors with the covariates included in the tessellaions);
  #and Tessellation Set (A list of matrices that give the coordinates of the centers in the tessellations)
  
  Pred<-rep(list(matrix(mean(yScaled)/m)),m)
  Dim=vector(length = m)
  Tess=vector(length = m)
  for (i in 1:m){
    Dim[i]<-list(sample(1:length(x[1,]), 1))
    Tess[i]<-(list(matrix(rnorm(1,0,sd))))
  }
  
  #Prepare some variables used in the backfitting algorithm
  SumOfAllTess=rep(mean(yScaled),length(yScaled))
  SigmaSquaredMu=(0.5/(k*sqrt(m)))^2
  LastTessPred=matrix
  
  #Matrices that will hold the samples from the poseterior distribution for the training samples and test samples.
  PredictionMatrix<-array(dim=c(length(y),(max_iter-burn_in)))
  TestMatrix<-array(dim=c(length(YTest),(max_iter-burn_in)))
  
  #finding lambda
  if (IntialSigma=="Naive"){ # Usually used if p is greater then n. Uses Standard deviation of y to predict sigma.
    SigmaSquaredHat=var(yScaled)
  }
  else{  # Default method using residual standard deviation from a least-squared linear regression of y, to predict sigma.
    MultiLinear<-lm(yScaled ~ xScaled)
    SigmaSquaredHat=sum(MultiLinear$residuals^2)/(length(yScaled)-length(xScaled[1,])-1)
  }
  
  #Find lambda
  lambda=1;
  lambda <- optim(par = 1,
                  fitting_function,
                  method = "Brent",
                  lower = 0.001,
                  upper = 100,
                  q=q , nu=nu, sigmaSquared_hat=SigmaSquaredHat)$par
  
  for (i in 1:max_iter){
    
    #Sample Sigma squared using all tessellations to predict the outcome variables
    SigmaSquared=SigmaSquaredCalculation(yScaled,nu,lambda,SumOfAllTess)
    
    
    for (j in 1:m){
      NewTessOutput<-NewTess(xScaled,j,Tess,Dim,sd) #Propose new Tessellation 
      TessStar<-NewTessOutput[[1]]  
      DimStar<-NewTessOutput[[2]]
      Modification<-NewTessOutput[[3]]
      
      ResidualsOutput<-CalculateResiduals(yScaled,xScaled,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred) #Calculate the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation and the number of observations in each cell.
      R_ijOld<-ResidualsOutput[[1]]   #Old and New refer to the original and proposed tessellations
      n_ijOld<-ResidualsOutput[[2]]
      R_ijNew<-ResidualsOutput[[3]]
      n_ijNew<-ResidualsOutput[[4]]
      SumOfAllTess<-ResidualsOutput[[5]] #Keeps track of the prediction for all tessellations to help sample sigma squared.
      IndexesStar<-ResidualsOutput[[6]] #Gives the row of each observation for the cell it falls in for the proposed tessellation.
      Indexes<-ResidualsOutput[[7]]  #Gives the row of each observation for the cell it falls in for the original tessellation.
      
      if (!any(n_ijNew==0)){ #automatically rejects proposed tessellation if there exists a cell with no observations in.
        LOGAcceptenceProb=AlphaCalculation(xScaled,TessStar,DimStar,j,R_ijOld,n_ijOld,R_ijNew,n_ijNew,SigmaSquared,Modification,SigmaSquaredMu,Omega,lambda_rate) #Gives the log of the acceptence probability.
        
        if (log(runif(n=1, min=0, max=1))<LOGAcceptenceProb){ #Accepts the proposed tessellation is accepted then calculates the new output values for the new tessellation. 
          Tess=TessStar
          Dim=DimStar
          Pred[[j]]=NewPredSet(j,TessStar,R_ijNew,n_ijNew,SigmaSquaredMu,SigmaSquared)
          LastTessPred=Pred[[j]][IndexesStar]
        }
        else { #Rejects the proposed tesellation then calculates new output values for the original tessellation.
          Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
          LastTessPred=Pred[[j]][Indexes];
        }
      }
      else{ #Rejects the proposed tesellation then calculates new output values for the original tessellation.
        Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
        LastTessPred=Pred[[j]][Indexes];
      }
      if (j==m){ #If j equals m then adds the last tessellation output values to give a prediction.
        SumOfAllTess=SumOfAllTess+LastTessPred;
      }
    }
    
    if (i>burn_in){ #vectors that hold the predictions for each iteration after burn in.
      PredictionMatrix[,(i-burn_in)]=SumOfAllTess;
      TestMatrix[,(i-burn_in)]=TestPrediction(XTest,m,Tess,Dim,Pred);
    }
    if (i %% 100 == 0){
      cat(sprintf("Iteration %d out of %d", i, max_iter), "\n")
    }
  }
  
  #finding the mean of the predition over the iterations and then unscaling the predictions.
  mean_yhat=(rowSums(PredictionMatrix)/(max_iter-burn_in))*(max(y)-min(y))+((max(y)+min(y))/2)
  mean_yhat_Test=(rowSums(TestMatrix)/(max_iter-burn_in))*(max(y)-min(y))+((max(y)+min(y))/2)
  
  
  return( #Returns the RMSE value for the test samples.
    data.frame(
      RMSE = sqrt(mean((YTest-mean_yhat_Test)^2)),
      mean_yhat_Test
    )
  )
}


f = function(x){
  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-0.5)^2+10*x[,4]+5*x[,5]
}

sigma = 1  #y = f(x) + sigma*z , z~N(0,1)
n = 300      #number of observations
set.seed(24)
X=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
Ey = f(X)
Y=Ey+sigma*rnorm(n)

n=length(Y)
TrainSet=sort(sample.int(n,3*n/6))
TestSet=1:n
TestSet=TestSet[! TestSet %in% TrainSet]


hyperparametersAddiVortesRobust <- matrix(
  c(6, 6,10,
    0.95,0.85,0.75,
    1,3,3,
    1.5,0.8,0.8,
    3,3,1,
    25,25,5), ncol=6
)



rmse_values<-numeric(3)
RobustnessValues<-array(dim=c(3,num_of_datasets,8))
FullRobustnessValues<-array(dim=c(3,num_of_datasets,8))

num_cores <- 10  # Specify the number of cores you want to use
cl <- makeCluster(num_cores)
registerDoParallel(cl)

for (m in c(1,10,50,100,200,300,400,500)){
  t<-which(m==c(1,10,50,100,200,300,400,500))
  RobustnessValues[,,t]<-foreach (i = 1:num_of_datasets,.combine=cbind) %dopar% {
    library('invgamma')
    library('FNN')
    library('expint')
    
    for (i in 1:nrow(hyperparametersAddiVortesRobust)){
      params<-hyperparametersAddiVortesRobust[i,]
      nu <- params[1]
      q<- params[2]
      k<- params[3]
      sd<- params[4]
      omega<- params[5]
      lambda<- params[6]
      
      rmse_values[i] <- AddiVortes_Algorithm(Y[TrainSet],as.matrix(X[TrainSet,]),m,1200,200,nu,q,k,sd,omega,lambda,f(X[TestSet,]),as.matrix(X[TestSet,]))$RMSE
    }
    
    return(rmse_values)
  }
  print(t)
}


stopCluster(cl)

FullMeanValues<-matrix(nrow = 3,ncol = 8)
FullLowerquantileValues<-matrix(nrow = 3,ncol = 8)
FullUpperquantileValues<-matrix(nrow = 3,ncol = 8)

for (i in 1:8){
  FullMeanValues[,i]<-rowMeans(RobustnessValues[,,i])
  FullLowerquantileValues[,i]<- apply(RobustnessValues[,,i], 1, quantile, probs = c(0.05, 0.95))[1,]
  FullUpperquantileValues[,i]<- apply(RobustnessValues[,,i], 1, quantile, probs = c(0.05, 0.95))[2,]
}


plotCI(c(1,10,50,100,200,300,400,500),FullMeanValues[1,],FullUpperquantileValues[1,]-FullMeanValues[1,],FullMeanValues[1,]-FullLowerquantileValues[1,], sfrac=0, scol = 'red',pch = 16,cex = 0.75,col = 'red',xlim = c(1,500),ylim = c(1,5),lwd=2,xlab = "Number of Tessellations, m", ylab = "RMSE" )

points(c(1-5,10-5,50-5,100-5,200-5,300-5,400-5,500-5),FullMeanValues[2,], pch = 15,cex = 0.75,col = 'blue')
arrows( c(1-5,10-5,50-5,100-5,200-5,300-5,400-5,500-5),FullUpperquantileValues[2,], c(1-5,10-5,50-5,100-5,200-5,300-5,400-5,500-5),FullLowerquantileValues[2,], angle = 90, code = 3, length = 0, col = "blue",lwd= 2)

points(c(1+5,10+5,50+5,100+5,200+5,300+5,400+5,500+5),FullMeanValues[3,], pch = 17,cex =0.75, col = 'darkgreen')
arrows( c(1+5,10+5,50+5,100+5,200+5,300+5,400+5,500+5),FullUpperquantileValues[3,], c(1+5,10+5,50+5,100+5,200+5,300+5,400+5,500+5),FullLowerquantileValues[3,], angle = 90, code = 3, length = 0, col = "darkgreen",lwd = 2)

legend(
  "topright",
  legend = c("Aggressive", "Default", "Conservative"),
  col = c('red', 'blue', 'darkgreen'),
  pch = c(16, 15, 17),
  lwd = 2,
  inset = 0.1,bty = "n"
)


}

figure8<-function(){
  f = function(x){
    10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-0.5)^2+10*x[,4]+5*x[,5]
  }
  
  par(mfrow =c(3,3))
  
  AddiVortes_Algorithm_Plot<-function(y,x,m,max_iter,burn_in,nu,q,k,sd,Omega,lambda_rate,YTest,XTest,IntialSigma = "Linear",p){
    
    #Scaling x and y
    yScaled=(y-(max(y)+min(y))/2)/(max(y)-min(y))
    xScaled=x;
    for (i in 1:length(x[1,])){
      xScaled[,i]=(x[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
    }
    
    
    for (i in 1:length(XTest[1,])){
      XTest[,i]=(XTest[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
    }
    
    
    #Initialize the Prediction Set, Dimension set and Tessellation Set
    
    Pred<-rep(list(matrix(mean(yScaled)/m)),m)
    Dim=vector(length = m)
    Tess=vector(length = m)
    for (i in 1:m){
      Dim[i]<-list(sample(1:length(x[1,]), 1))
      Tess[i]<-(list(matrix(rnorm(1,0,sd))))
    }
    
    #Prepare some variables used in the backfitting algorithm
    SumOfAllTess=rep(mean(yScaled),length(yScaled))
    SigmaSquaredMu=(0.5/(k*sqrt(m)))^2
    LastTessPred=matrix
    
    #Matrices that will hold the samples from the poseterior distribution for the training samples and test samples.
    PredictionMatrix<-array(dim=c(length(y),(max_iter-burn_in)))
    TestMatrix<-array(dim=c(length(YTest),(max_iter-burn_in)))
    TestMatrix<-array(dim=c(length(YTest),(max_iter-burn_in)))
    plotForSigmaSquared<-vector(length = max_iter)
    plotForRMSE<-vector(length = max_iter)
    
    #finding lambda
    if (IntialSigma=="Naive"){ # Usually used if p is greater then n. Uses Standard deviation of y to predict sigma.
      SigmaSquaredHat=var(yScaled)
    }
    else{  # Default method using residual standard deviation from a least-squared linear regression of y, to predict sigma.
      MultiLinear<-lm(yScaled ~ xScaled)
      SigmaSquaredHat=sum(MultiLinear$residuals^2)/(length(yScaled)-length(xScaled[1,])-1)
    }
    
    #Find lambda
    lambda=1;
    lambda <- optim(par = 1,
                    fitting_function,
                    method = "Brent",
                    lower = 0.001,
                    upper = 100,
                    q=q , nu=nu, sigmaSquared_hat=SigmaSquaredHat)$par
    
    for (i in 1:max_iter){
      
      #Sample Sigma squared using all tessellations to predict the outcome variables
      SigmaSquared=SigmaSquaredCalculation(yScaled,nu,lambda,SumOfAllTess)
      plotForSigmaSquared[i]=(SigmaSquared*(max(y)-min(y))^2)^0.5
      plotForRMSE[i]=(mean((yScaled-SumOfAllTess)^2))^0.5
      
      for (j in 1:m){
        NewTessOutput<-NewTess(xScaled,j,Tess,Dim,sd) #Propose new Tessellation 
        TessStar<-NewTessOutput[[1]]  
        DimStar<-NewTessOutput[[2]]
        Modification<-NewTessOutput[[3]]
        
        ResidualsOutput<-CalculateResiduals(yScaled,xScaled,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred) #Calculate the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation and the number of observations in each cell.
        R_ijOld<-ResidualsOutput[[1]]   #Old and New refer to the original and proposed tessellations
        n_ijOld<-ResidualsOutput[[2]]
        R_ijNew<-ResidualsOutput[[3]]
        n_ijNew<-ResidualsOutput[[4]]
        SumOfAllTess<-ResidualsOutput[[5]] #Keeps track of the prediction for all tessellations to help sample sigma squared.
        IndexesStar<-ResidualsOutput[[6]] #Gives the row of each observation for the cell it falls in for the proposed tessellation.
        Indexes<-ResidualsOutput[[7]]  #Gives the row of each observation for the cell it falls in for the original tessellation.
        
        if (!any(n_ijNew==0)){ #automatically rejects proposed tessellation if there exists a cell with no observations in.
          
          LOGAcceptenceProb=AlphaCalculation(xScaled,TessStar,DimStar,j,R_ijOld,n_ijOld,R_ijNew,n_ijNew,SigmaSquared,Modification,SigmaSquaredMu,Omega,lambda_rate) #Gives the log of the acceptence probability.
          
          
          if (log(runif(n=1, min=0, max=1))<LOGAcceptenceProb){ #Accepts the proposed tessellation is accepted then calculates the new output values for the new tessellation. 
            Tess=TessStar
            Dim=DimStar
            Pred[[j]]=NewPredSet(j,TessStar,R_ijNew,n_ijNew,SigmaSquaredMu,SigmaSquared)
            LastTessPred=Pred[[j]][IndexesStar]
          }
          else { #Rejects the proposed tesellation then calculates new output values for the original tessellation.
            Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
            LastTessPred=Pred[[j]][Indexes];
          }
        }
        else{ #Rejects the proposed tesellation then calculates new output values for the original tessellation.
          Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
          LastTessPred=Pred[[j]][Indexes];
        }
        if (j==m){ #If j equals m then adds the last tessellation output values to give a prediction.
          SumOfAllTess=SumOfAllTess+LastTessPred;
        }
      }
      
      if (i>burn_in){ #vectors that hold the predictions for each iteration after burn in.
        PredictionMatrix[,(i-burn_in)]=SumOfAllTess;
        TestMatrix[,(i-burn_in)]=TestPrediction(XTest,m,Tess,Dim,Pred);
      }
      if (i %% 100 == 0){
        cat(sprintf("Iteration %d out of %d", i, max_iter), "\n")
      }
    }
    
    #finding the mean of the predition over the iterations and then unscaling the predictions.
    mean_yhat=(rowSums(PredictionMatrix)/(max_iter-burn_in))*(max(y)-min(y))+((max(y)+min(y))/2)
    mean_yhat_Test=(rowSums(TestMatrix)/(max_iter-burn_in))*(max(y)-min(y))+((max(y)+min(y))/2)
    mean_yhat_Test_ColSum=(colSums(TestMatrix)/(length(TestMatrix[,1])))*(max(y)-min(y))+((max(y)+min(y))/2)
    
    LowerConfidenceTRAINValue<-vector(length=length(mean_yhat))
    UpperConfidenceTRAINValue<-vector(length=length(mean_yhat))
    
    for (i in 1:length(mean_yhat)){
      PredictionMatrix[i,]<-sort(PredictionMatrix[i,])
      
      if ((((max_iter-burn_in+1)*0.05))== round((max_iter-burn_in+1)*0.05)){
        LowerConfidenceTRAINValue[i]<-(PredictionMatrix[i,(max_iter-burn_in+1)*0.05])*(max(y)-min(y))+((max(y)+min(y))/2)
        UpperConfidenceTRAINValue[i]<-(PredictionMatrix[i,(max_iter-burn_in+1)*0.95])*(max(y)-min(y))+((max(y)+min(y))/2)
      }
      else{
        LowerConfidenceTRAINValue[i]<-((PredictionMatrix[i,trunc((max_iter-burn_in+1)*0.05)]+PredictionMatrix[i,trunc(((max_iter-burn_in+1)*0.05)+1)])/2)*(max(y)-min(y))+((max(y)+min(y))/2)
        UpperConfidenceTRAINValue[i]<-((PredictionMatrix[i,trunc((max_iter-burn_in+1)*0.95)]+PredictionMatrix[i,trunc(((max_iter-burn_in+1)*0.95)+1)])/2)*(max(y)-min(y))+((max(y)+min(y))/2)
      }
    }
    
    LowerConfidenceTESTValue<-vector(length=length(mean_yhat_Test))
    UpperConfidenceTESTValue<-vector(length=length(mean_yhat_Test))
    
    for (i in 1:length(mean_yhat_Test)){
      TestMatrix[i,]<-sort(TestMatrix[i,])
      
      if ((((max_iter-burn_in+1)*0.05))== round((max_iter-burn_in+1)*0.05)){
        LowerConfidenceTESTValue[i]<-(TestMatrix[i,(max_iter-burn_in+1)*0.05])*(max(y)-min(y))+((max(y)+min(y))/2)
        UpperConfidenceTESTValue[i]<-(TestMatrix[i,(max_iter-burn_in+1)*0.95])*(max(y)-min(y))+((max(y)+min(y))/2)
      }
      else{
        LowerConfidenceTESTValue[i]<-((TestMatrix[i,trunc((max_iter-burn_in+1)*0.05)]+TestMatrix[i,trunc(((max_iter-burn_in+1)*0.05)+1)])/2)*(max(y)-min(y))+((max(y)+min(y))/2)
        UpperConfidenceTESTValue[i]<-((TestMatrix[i,trunc((max_iter-burn_in+1)*0.95)]+TestMatrix[i,trunc(((max_iter-burn_in+1)*0.95)+1)])/2)*(max(y)-min(y))+((max(y)+min(y))/2)
      }
    }
    
    plotCI(f(x),mean_yhat, UpperConfidenceTRAINValue-mean_yhat,mean_yhat-LowerConfidenceTRAINValue,sfrac=0, scol = 'grey',ylab = "posterior intervals", xlab = "In-Sample f(x)", pch = 16, cex.lab = 1.5, main = paste("p=", p)) 
    abline(0,1)
    
    plotCI(YTest,mean_yhat_Test,UpperConfidenceTESTValue-mean_yhat_Test,mean_yhat_Test-LowerConfidenceTESTValue,sfrac=0, scol = 'grey', ylab = "posterior intervals", xlab = "Out-of-Sample f(x)",pch = 16, cex.lab = 1.5, main = paste("p=", p))
    abline(0,1)
    
    in_interval_Train <- (f(x) >= LowerConfidenceTRAINValue) & (f(x) <= UpperConfidenceTRAINValue)
    in_interval_TEST<- (YTest >= LowerConfidenceTESTValue) & (YTest <= UpperConfidenceTESTValue)
    
    #plot(plotForSigmaSquared,col=c(rep('red',burn_in),rep('black',max_iter-burn_in)),xlab="MCMC iteration",ylab="Sigma draw",type="l")
    plot(plotForSigmaSquared, type = "n", col = "black", lwd = 2, xlab="MCMC iteration",ylab="Sigma draw", cex.lab = 1.5, main = paste("p=", p))
    
    # Draw the first segment in red
    segments(x0 = 1:(burn_in - 1), y0 = plotForSigmaSquared[1:(burn_in - 1)],
             x1 = 2:burn_in, y1 = plotForSigmaSquared[2:burn_in], col = "red", lwd = 2)
    # Draw the second segment in black
    segments(x0 = (burn_in):(max_iter - 1), y0 = plotForSigmaSquared[(burn_in+1):(max_iter - 1)],
             x1 = (burn_in+1):max_iter, y1 = plotForSigmaSquared[(burn_in):max_iter], col = "black", lwd = 2)
    abline(1,0)
    
    print( data.frame(
      invertal_Train = sum(in_interval_Train)/length(in_interval_Train),
      invertal_Test = sum(in_interval_TEST)/length(in_interval_TEST),
      RMSE = sqrt(mean((YTest-mean_yhat_Test)^2))
    ))
    
  }

  
  sigma = 1  #y = f(x) + sigma*z , z~N(0,1)
  n = 300      #number of observations
  set.seed(9)
  X=matrix(runif(n*20),n,20) #20 variables, only first 5 matter
  Ey = f(X)
  Y=Ey+sigma*rnorm(n)
  
  n=length(Y)
  TrainSet=sort(sample.int(n,3*n/6))
  TestSet=1:n
  TestSet=TestSet[! TestSet %in% TrainSet]
  
  AddiVortes_Algorithm_Plot(Y[TrainSet],X[TrainSet,],50,2000,200,6,0.85,3,0.8,3,25,f(X[TestSet,]),X[TestSet,],IntialSigma = "Linear", p=20)
  
  sigma = 1  #y = f(x) + sigma*z , z~N(0,1)
  n = 300      #number of observations
  X=matrix(runif(n*100),n,100) #100 variables, only first 5 matter
  Ey = f(X)
  Y=Ey+sigma*rnorm(n)
  
  n=length(Y)
  TrainSet=sort(sample.int(n,3*n/6))
  TestSet=1:n
  TestSet=TestSet[! TestSet %in% TrainSet]
  
  AddiVortes_Algorithm_Plot(Y[TrainSet],X[TrainSet,],50,2000,200,6,0.99,3,0.2,3,25,f(X[TestSet,]),X[TestSet,],IntialSigma = "Naive", p=100)
  
  sigma = 1  #y = f(x) + sigma*z , z~N(0,1)
  n = 300      #number of observations
  X=matrix(runif(n*1000),n,1000) #1000 variables, only first 5 matter
  Ey = f(X)
  Y=Ey+sigma*rnorm(n)
  
  n=length(Y)
  TrainSet=sort(sample.int(n,3*n/6))
  TestSet=1:n
  TestSet=TestSet[! TestSet %in% TrainSet]
  
  AddiVortes_Algorithm_Plot(Y[TrainSet],X[TrainSet,],50,2000,200,6,0.99,3,0.2,3,25,f(X[TestSet,]),X[TestSet,],IntialSigma = "Naive", p=1000)
  
}


figure9<-function(max_iter = 1200, burn_in= 200, num_of_datasets= 100){

  AddiVortes_Algorithm<-function(y,x,m = 200 ,max_iter = 1200,burn_in= 200,nu = 6,q =0.85,k = 3 ,sd = 0.8 ,Omega = 3,lambda_rate = 25,YTest,XTest,IntialSigma = "Linear"){
    
    #Scaling x and y
    yScaled=(y-(max(y)+min(y))/2)/(max(y)-min(y))
    xScaled=x;
    for (i in 1:length(x[1,])){
      xScaled[,i]=(x[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
    }
    
    
    for (i in 1:length(XTest[1,])){
      XTest[,i]=(XTest[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
    }
    
    #Initialize:
    #Prediction Set (A list of vectors with the output values for each tessellation), 
    #Dimension set (A list of vectors with the covariates included in the tessellaions);
    #and Tessellation Set (A list of matrices that give the coordinates of the centers in the tessellations)
    
    Pred<-rep(list(matrix(mean(yScaled)/m)),m)
    Dim=vector(length = m)
    Tess=vector(length = m)
    for (i in 1:m){
      Dim[i]<-list(sample(1:length(x[1,]), 1))
      Tess[i]<-(list(matrix(rnorm(1,0,sd))))
    }
    
    #Prepare some variables used in the backfitting algorithm
    SumOfAllTess=rep(mean(yScaled),length(yScaled))
    SigmaSquaredMu=(0.5/(k*sqrt(m)))^2
    LastTessPred=matrix
    
    #Matrices that will hold the samples from the poseterior distribution for the training samples and test samples.
    PredictionMatrix<-array(dim=c(length(y),(max_iter-burn_in)))
    TestMatrix<-array(dim=c(length(YTest),(max_iter-burn_in)))
    
    #finding lambda
    if (IntialSigma=="Naive"){ # Usually used if p is greater then n. Uses Standard deviation of y to predict sigma.
      SigmaSquaredHat=var(yScaled)
    }
    else{  # Default method using residual standard deviation from a least-squared linear regression of y, to predict sigma.
      MultiLinear<-lm(yScaled ~ xScaled)
      SigmaSquaredHat=sum(MultiLinear$residuals^2)/(length(yScaled)-length(xScaled[1,])-1)
    }
    
    #Find lambda
    lambda=1;
    lambda <- optim(par = 1,
                    fitting_function,
                    method = "Brent",
                    lower = 0.001,
                    upper = 100,
                    q=q , nu=nu, sigmaSquared_hat=SigmaSquaredHat)$par
    
    for (i in 1:max_iter){
      
      #Sample Sigma squared using all tessellations to predict the outcome variables
      SigmaSquared=SigmaSquaredCalculation(yScaled,nu,lambda,SumOfAllTess)
      
      for (j in 1:m){
        NewTessOutput<-NewTess(xScaled,j,Tess,Dim,sd) #Propose new Tessellation 
        TessStar<-NewTessOutput[[1]]  
        DimStar<-NewTessOutput[[2]]
        Modification<-NewTessOutput[[3]]
        
        ResidualsOutput<-CalculateResiduals(yScaled,xScaled,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred) #Calculate the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation and the number of observations in each cell.
        R_ijOld<-ResidualsOutput[[1]]   #Old and New refer to the original and proposed tessellations
        n_ijOld<-ResidualsOutput[[2]]
        R_ijNew<-ResidualsOutput[[3]]
        n_ijNew<-ResidualsOutput[[4]]
        SumOfAllTess<-ResidualsOutput[[5]] #Keeps track of the prediction for all tessellations to help sample sigma squared.
        IndexesStar<-ResidualsOutput[[6]] #Gives the row of each observation for the cell it falls in for the proposed tessellation.
        Indexes<-ResidualsOutput[[7]]  #Gives the row of each observation for the cell it falls in for the original tessellation.
        
        if (!any(n_ijNew==0)){ #automatically rejects proposed tessellation if there exists a cell with no observations in.
          
          LOGAcceptenceProb=AlphaCalculation(xScaled,TessStar,DimStar,j,R_ijOld,n_ijOld,R_ijNew,n_ijNew,SigmaSquared,Modification,SigmaSquaredMu,Omega,lambda_rate) #Gives the log of the acceptence probability.
          
          if (log(runif(n=1, min=0, max=1))<LOGAcceptenceProb){ #Accepts the proposed tessellation is accepted then calculates the new output values for the new tessellation. 
            Tess=TessStar
            Dim=DimStar
            Pred[[j]]=NewPredSet(j,TessStar,R_ijNew,n_ijNew,SigmaSquaredMu,SigmaSquared)
            LastTessPred=Pred[[j]][IndexesStar]
          }
          else { #Rejects the proposed tesellation then calculates new output values for the original tessellation.
            Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
            LastTessPred=Pred[[j]][Indexes];
          }
        }
        else{ #Rejects the proposed tesellation then calculates new output values for the original tessellation.
          Pred[[j]]=NewPredSet(j,Tess,R_ijOld,n_ijOld,SigmaSquaredMu,SigmaSquared);
          LastTessPred=Pred[[j]][Indexes];
        }
        if (j==m){ #If j equals m then adds the last tessellation output values to give a prediction.
          SumOfAllTess=SumOfAllTess+LastTessPred;
        }
      }
      
      if (i>burn_in){ #vectors that hold the predictions for each iteration after burn in.
        PredictionMatrix[,(i-burn_in)]=SumOfAllTess;
        TestMatrix[,(i-burn_in)]=TestPrediction(XTest,m,Tess,Dim,Pred);
      }
    }
    
    #finding the mean of the predition over the iterations and then unscaling the predictions.
    mean_yhat=(rowSums(PredictionMatrix)/(max_iter-burn_in))*(max(y)-min(y))+((max(y)+min(y))/2)
    mean_yhat_Test=(rowSums(TestMatrix)/(max_iter-burn_in))*(max(y)-min(y))+((max(y)+min(y))/2)
    
    return( #Returns the RMSE value for the test samples.
      data.frame(
        Out_of_sample_RMSE = sqrt(mean((YTest-mean_yhat_Test)^2))
      )
    )
  }
  
  
    SigmaSquaredCalculation<-function(yScaled,nu,lambda,SumOfAllTess){ #Sample sigma squared from inverse gamma distribution
    
      n=length(yScaled)
      SigmaSquared<-rinvgamma(1,shape=(nu+n)/2,rate=(nu*lambda+sum((yScaled-SumOfAllTess)^2))/2)
      
      return(SigmaSquared)
    }
    
    NewTess<-function(x,j,Tess,Dim,sd){ #Propose a new tessellation
      
      p=runif(1,0,1) #Randomly sample p to decide the proposed modification to tessellation.
      
      DimStar=Dim # Let proposed dimension matrix equal original dimension matrix.
      TessStar=Tess #Similar for the tessellation matrix.
      
      if (p<0.2 & length(Dim[[j]])!=length(x[1,]) | length(Dim[[j]])==1 & p<0.4){ #Add a dimension if p is less then 0.2 or if p is less then 0.4 when there is only one dimension in the Tessellation due to adjustments (Supplementary Material).
        NumberOfCovariates=1:length(x[1,]) #Let NumberOfCovariates be a vector from 1 to the number of covariates considered.
        NumberOfCovariates=NumberOfCovariates[-Dim[[j]]] #Remove all values that for the covariates that are already in the tessellation.
        DimStar[[j]]<-c(Dim[[j]],sample(NumberOfCovariates,1)) # Uniformly sample a new covariate and add it to the dimension matrix.
        TessStar[[j]]=cbind(Tess[[j]],rnorm(length(Tess[[j]][,1]),0,sd)) # Sample new coordinates from Normal distribution for the new dimension and add it to the Tessellation matrix.
        Modification="AD"}
      else if (p<0.4){ #Remove a dimension if p is less then 0.4.
        RemovedDim=sample(1:length(Dim[[j]]),1) #Uniformly sample the dimension to be removed.
        DimStar[[j]]=DimStar[[j]][-RemovedDim] #Remove the dimension from the dimesion Matrix.
        TessStar[[j]]=matrix(TessStar[[j]][,-RemovedDim],ncol=length(DimStar[[j]])) #Remove the coordinates in the Tessellation matrix corresponding to the dimension removed.
        Modification="RD"}
      else if (p<0.6 || p<0.8 & length(Tess[[j]][,1])==1){ #Add a centre if p is less then 0.6 or if p is less then 0.4 when there is only one center in the Tessellation due to adjustments (Supplementary Material).
        TessStar[[j]]=rbind(Tess[[j]],rnorm(length(Dim[[j]]),0,sd)) #Add a new row of coordinates, sampled from a normal distribution, to the Tessellation matrix to add a center.
        Modification="AC"}
      else if (p<0.8){ #Add a centre if p is less then 0.8. 
        CenterRemoved=sample(1:length(TessStar[[j]][,1]),1) #Sample a row.
        TessStar[[j]]=matrix(TessStar[[j]][-CenterRemoved,],ncol=length(Dim[[j]])) #Remove row sampled.
        Modification="RC"}
      else if (p<0.9 || length(Dim[[j]])==length(x[1,])){ #Change a center if p is less then 0.9 or if the all the covariates are in the tessellation.
        TessStar[[j]][sample(1:length(TessStar[[j]][,1]),1),]=rnorm(length(Dim[[j]]),0,sd) # Sample a row in the tessellaion matrix and change the coordinates of the centre by sampling from a normal distribution.
        Modification="Change"}
      else{ #Swop a dimension.
        NumberOfCovariates=1:length(x[1,])  #Let NumberOfCovariates be a vector from 1 to the number of covariates considered.
        NumberOfCovariates=NumberOfCovariates[-Dim[[j]]]  #Remove all values that for the covariates that are already in the tessellation.
        DimToChange=sample(1:length(Dim[[j]]),1) #Uniformly sample a dimension to change.
        DimStar[[j]][DimToChange]=sample(NumberOfCovariates,1) #Replace the Dimension to a new uniforly sampled covariate that is not already in the tessellaion.
        TessStar[[j]][,DimToChange]=rnorm(length(Tess[[j]][,1]),0,sd) #Add new normally sampled coordinates new dimension added.
        Modification="Swop"}
      
      TessStar[[j]]<-matrix(TessStar[[j]],ncol=length(DimStar[[j]])) #Ensure the the Tessellation matrix is a "matrix" type.
      
      return(list(TessStar,DimStar,Modification)) #Return new proposed tessellation.
    }
    
    fitting_function<- function(lambda,q,nu,sigmaSquared_hat){ #function that calculates the squared difference between sigma squared hat and the inverse gamma function
      return((sigmaSquared_hat- qinvgamma(q, shape=nu/2, rate=nu*lambda/2))^2)
    }
    
    Indexes<-function(x,Tess,Dim){ #Gives the row (the center) of the tessellation that each obseravtion falls within.
      if (length(Tess[,1])==1){ #only 1 centre
        CellsForGivenTess=rep(1,length(x[,1]))
      }
      else{ #multiple
        CellsForGivenTess=knnx.index(Tess,matrix(x[,Dim],ncol = length(Dim)),1)
      }
      return(CellsForGivenTess)
    }
    
    AlphaCalculation<-function(x,Tess,Dim,j,R_ijOld,n_ijOld,R_ijNew,n_ijNew,SigmaSquared,Modification,SigmaSquaredMu,Omega,lambda_rate){ #Calculates the acceptence rate of the proposed tessellation.
      
      d=length(Dim[[j]]);
      NumCovariates=length(x[1,]);
      cStar=length(Tess[[j]][,1]);
    
      #The Log Likelihood Ratio in the acceptence ratio
      LOGlikelihoodRatio=0.5*(log(prod(n_ijOld*SigmaSquaredMu+SigmaSquared))-log(prod(n_ijNew*SigmaSquaredMu+SigmaSquared)))+((SigmaSquaredMu/(2*SigmaSquared))*(-sum((R_ijOld^2)/(n_ijOld*SigmaSquaredMu+SigmaSquared))+sum((R_ijNew^2)/(n_ijNew*SigmaSquaredMu+SigmaSquared))))
    
      #Calculating the acceptence probablity for "AD"=Adding a dimension, "RD"=Removing a dimension, "AC"=Adding a center, "RC"=Removing a center, "Change"=Changing the coordinates of a center and Swopping a dimension.
      if (Modification == "AD"){ 
        TessStructure=(dbinom(d-1,NumCovariates-1,Omega/NumCovariates))/(dbinom(d-2,NumCovariates-1,Omega/NumCovariates)*(NumCovariates-d+1))
        TransitionRatio=(NumCovariates-d+1)/d;
        AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
        
        #Adjustments.
        if (length(Dim[[j]])==1){
          AcceptenceProb=AcceptenceProb+log(1/2)
        }
        else if (length(Dim[[j]])==NumCovariates-1){
          AcceptenceProb=AcceptenceProb+log(2)}
      }
      else if (Modification == "RD"){
        TessStructure=(dbinom(d-1,NumCovariates,Omega/NumCovariates)*(NumCovariates-d))/(dbinom(d,NumCovariates,Omega/NumCovariates))
        TransitionRatio=(d+1)/(NumCovariates-d)
        AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
        
        #Adjustments.
        if (length(Dim[[j]])==NumCovariates){
          AcceptenceProb=AcceptenceProb+log(1/2)
        }
        else if (length(Dim[[j]])==2){
          AcceptenceProb=AcceptenceProb+log(2);
        }
      }
      else if(Modification == "AC"){
        TessStructure=dpois(cStar-1,lambda_rate)/dpois(cStar-2,lambda_rate)
        TransitionRatio=1/cStar;
        AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)+0.5*log(SigmaSquared)
    
        #Adjustments.
        if (cStar==1){
          AcceptenceProb=AcceptenceProb+log(1/2);
        }
      }
      else if (Modification == "RC"){
        TessStructure=dpois(cStar-1,lambda_rate)/dpois(cStar,lambda_rate);
        TransitionRatio=cStar+1;
        AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)-0.5*log(SigmaSquared)
    
        #Adjustments.
        if (cStar==2){
          AcceptenceProb=AcceptenceProb+log(2);
        }
      }
      else if (Modification == "Change"){
        TessStructure=1;
        TransitionRatio=1;
        AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio);
      }
      else {
        TessStructure=1;
        TransitionRatio=1;
        AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio);
      }
      return(AcceptenceProb)
    }
    
    CalculateResiduals<-function(y,x,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred){ #A function that calculates the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation.
      if (j==1){
        indexes=Indexes(x,Tess[[j]],Dim[[j]]);
        CurrentTessPred<-Pred[[j]][indexes]
        SumOfAllTess=SumOfAllTess-CurrentTessPred}
      else{
        indexes=Indexes(x,Tess[[j]],Dim[[j]]);
        CurrentTessPred<-Pred[[j]][indexes]
        SumOfAllTess=SumOfAllTess+LastTessPred-CurrentTessPred;
      }
      
      IndexesStar=Indexes(x,TessStar[[j]],DimStar[[j]]);
      R_j<-y-SumOfAllTess
      
      #Initializing Sizes
      
      R_ijOld=rep(0,length(Pred[[j]]))
      n_ijOld=rep(0,length(Pred[[j]]))
      
      for (i in 1:length(Pred[[j]])){
        R_ijOld[i]<-sum(R_j[indexes==i])
        n_ijOld[i]<-sum(indexes==i)
      }
      
      R_ijNew=rep(0,length(TessStar[[j]][,1]))
      n_ijNew=rep(0,length(TessStar[[j]][,1]))
      
      for (i in 1:length(TessStar[[j]][,1])){
        R_ijNew[i]<-sum(R_j[IndexesStar==i])
        n_ijNew[i]<-sum(IndexesStar==i)}
      
      return(list(R_ijOld,n_ijOld,R_ijNew,n_ijNew,SumOfAllTess,IndexesStar,indexes))
    }
    
    NewPredSet<-function(j,Tess,R_ijNew,n_ijNew,sigmaSquaredMu,SigmaSquared){ #Sampling the new output values for the new tessellation.
      PredSet=rep(0,length(Tess[[j]][,1]))
      for (i in 1:length(Tess[[j]][,1])){
        PredSet[i]=rnorm(1,sigmaSquaredMu*R_ijNew[i]/(sigmaSquaredMu*n_ijNew[i]+SigmaSquared),((SigmaSquared*sigmaSquaredMu)/(n_ijNew[i]*sigmaSquaredMu+SigmaSquared))^0.5);
      }
      return(PredSet)
    }
    
    TestPrediction<-function(x,m,Tess,Dim,Pred){ #A function that derives a prediction for the output variable give an observation, for one iteration in AddiVortes.
      Prediction=rep(0,length(x[,1]));
      for (j in 1:m){
        NewTessIndexes=Indexes(x,Tess[[j]],Dim[[j]]);
        Prediction=Prediction+Pred[[j]][NewTessIndexes]
      }
      return(Prediction)
    }
    
set.seed(32)

f = function(x){
  10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-0.5)^2+10*x[,4]+5*x[,5]
}

sigma = 1  #y = f(x) + sigma*z , z~N(0,1)
n = 1000      #number of observations

#100 Friedman datasets
BenchmarkX<-list()
BenchmarkY<-list()
for (i in 1:num_of_datasets){
  BenchmarkX[[i]]=matrix(runif(n*10),n,10) #10 variables, only first 5 matter
  Ey = f(BenchmarkX[[i]])
  BenchmarkY[[i]]=Ey+sigma*rnorm(n)
}

# Number of cross-validation folds
NumOfRep <- 1
num_folds <- 5

createFolds <- function(data, k) {
  set.seed(123)  # For reproducibility
  fold_indices <- sample(rep(1:k, length.out = nrow(data)))
  fold_list <- split(1:nrow(data), fold_indices)
  return(fold_list)
}

# Create a data frame with all combinations of hyperparameter values
hyperparametersBART <- expand.grid(
  m=c(50,200),
  nu_q=list(c(3,0.9),c(3,0.99),c(10,0.75)),
  k=c(1,2,3,5)
)


hyperparametersBoost <- expand.grid(
  depth=c(1,2,3,4),
  num_trees=c(50,100,200),
  eta=c(0.01, 0.05, 0.10, 0.25)
)

hyperparametersSoftBART <- expand.grid(
  m=c(50,200),
  k=c(1,2,3,5)
)

hyperparametersAddiVortes <- expand.grid(
  m=c(50,200),
  nu=c(6),
  q=c(0.85),
  k=c(1,3),
  sd=c(0.8,1.5),
  omega=c(3),
  lambda=c(5,25)
)

X<-BenchmarkX[[1]]
Y<-BenchmarkY[[1]]

#Cross validation
folds <- createFolds(X, k = num_folds)  # 5-fold cross-validation

num_cores <- 10  # Specify the number of cores you want to use
cl3 <- makeCluster(num_cores)
registerDoParallel(cl3)


# Initialize a matrix to store the results
results_matrixAddiVortes <- matrix(ncol = 8, nrow = nrow(hyperparametersAddiVortes))
rmse_values<-vector(length = num_folds)

# Iterate over hyperparameter combinations
results_matrixAddiVortes<-foreach(i = 1:nrow(hyperparametersAddiVortes), .combine = rbind) %dopar% {
  library('invgamma')
  library('FNN')
  library('expint')
  
  params <- hyperparametersAddiVortes[i, ]
  rmse_values <- numeric(num_folds)
  
  # Perform cross-validation using a loop
  for (j in 1:num_folds) {
    TestSetCV <- unlist(folds[-j],use.names = FALSE)
    TrainSetCV <- folds[[j]]
    
    m<-params$m
    nu <- params$nu
    q<- params$q
    sd<- params$sd
    omega<- params$omega
    lambda<- params$lambda
    k<- params$k
    
    # Make predictions on the test set
    rmse_values[j] <- AddiVortes_Algorithm(Y[TrainSetCV],as.matrix(X[TrainSetCV,]),m,1200,200,nu,q,k,sd,omega,lambda,f(X[TestSetCV,]),as.matrix(X[TestSetCV,]))$Out_of_sample_RMSE
  }                   
  # Store the average RMSE from cross-validation
  return(c(m,nu,q,k,sd,omega,lambda,mean(rmse_values)))
  
}

# Find the best hyperparameter combination based on RMSE
best_index <- which.min(results_matrixAddiVortes[, 8])
best_hyperparametersAddiVortes <- results_matrixAddiVortes[best_index, 1:7]

print("Finished AddiVortes Tuning")

# Initialize a matrix to store the results
results_matrixBoost <- matrix(ncol = 4, nrow = nrow(hyperparametersBoost))

# Iterate over hyperparameter combinations
results_matrixBoost<-foreach(i = 1:nrow(hyperparametersBoost), .combine = rbind) %dopar% {
  library('gbm')
  
  params <- hyperparametersBoost[i, ]
  rmse_values <- numeric(num_folds)
  
  # Perform cross-validation using a loop
  for (j in 1:num_folds) {
    TrainSetCV <- unlist(folds[-j],use.names = FALSE)
    TestSetCV <- folds[[j]]
    
    data_frame_Train<-as.data.frame(cbind(Y[TrainSetCV],X[TrainSetCV,]))
    data_frame_Test<-as.data.frame(cbind(f(X[TestSetCV,]),X[TestSetCV,]))
    
    depth <- params$depth
    num_trees <- params$num_trees
    eta <- params$eta
    
    # Train the bart model
    gradBoost<-gbm(V1~.,data=data_frame_Train,distribution = "gaussian", interaction.depth=depth,n.trees = num_trees, shrinkage = eta)
    
    # Make predictions on the test set
    rmse_values[j] <- (mean((f(X[TestSetCV,])-predict(gradBoost,newdata=data_frame_Test))**2))**0.5
  }
  
  # Store the average RMSE from cross-validation
  results_matrixBoost[i,]<-c(depth,num_trees,eta,mean(rmse_values))
  
}

# Find the best hyperparameter combination based on RMSE
best_index <- which.min(results_matrixBoost[, 4])
best_hyperparametersBoost <- results_matrixBoost[best_index, 1:3]

print('FinishedBoostTuning')

# Initialize a matrix to store the results
results_matrixBART <- matrix(ncol = 5, nrow = nrow(hyperparametersBART))

# Iterate over hyperparameter combinations
results_matrixBART<-foreach (i = 1:nrow(hyperparametersBART), .combine = rbind) %dopar% {
  library('BayesTree')
  
  params <- hyperparametersBART[i, ]
  rmse_values <- numeric(num_folds)
  
  # Perform cross-validation using a loop
  for (j in 1:num_folds) {
    TestSetCV <- unlist(folds[-j],use.names = FALSE)
    TrainSetCV <- folds[[j]]
    
    m <- params$m
    nu <- params$nu_q[[1]][1]
    q <- params$nu_q[[1]][2]
    k <- params$k
    
    # Train the bart model
    BartOr<-bart(as.matrix(X[TrainSetCV,]),as.numeric(as.matrix(Y[TrainSetCV])),as.matrix(X[TestSetCV,]), ntree = m, sigdf = nu,sigquant = q,k = k)
    
    # Make predictions on the test set
    rmse_values[j] <- (mean((f(X[TestSetCV,])-BartOr$yhat.test.mean)**2))**0.5
  }
  
  # Store the average RMSE from cross-validation
  results_matrixBART[i,] <- cbind(params$m,params$nu_q[[1]][1],params$nu_q[[1]][2],params$k,mean(rmse_values))
  
}

print('FinishedBARTTuning')

# Find the best hyperparameter combination based on RMSE
best_index <- which.min(results_matrixBART[, 5])
best_hyperparametersBART <- results_matrixBART[best_index, 1:4]


results_matrixSoftBART <- matrix(ncol = 3, nrow = nrow(hyperparametersSoftBART))

results_matrixSoftBART<-foreach (i = 1:nrow(hyperparametersSoftBART), .combine = rbind) %dopar% {
  library(SoftBart)
  
  params <- hyperparametersSoftBART[i, ]
  rmse_values <- numeric(num_folds)
  
  # Perform cross-validation using a loop
  for (j in 1:num_folds) {
    TestSetCV <- unlist(folds[-j],use.names = FALSE)
    TrainSetCV <- folds[[j]]
    
    m <- params$m
    ##nu <- params$nu_q[[1]][1]
    ##q <- params$nu_q[[1]][2]
    k <- params$k
    
    df <- data.frame(X = X[TrainSetCV,], Y = Y[TrainSetCV])
    df_test <- data.frame(X = X[TestSetCV,], Y = f(X[TestSetCV,]))
    
    # Train the bart model
    SoftBARTOr<-softbart_regression(Y ~ ., df, df_test, num_tree = m,k =k)
    
    # Make predictions on the test set
    rmse_values[j] <- sqrt(mean((SoftBARTOr$mu_test_mean-f(X[TestSetCV,]))^2))
  }
  
  # Store the average RMSE from cross-validation
  results_matrixSoftBART[i,] <- cbind(params$m,params$k,mean(rmse_values))
  
}

print('Finished SoftBART Tuning')

# Find the best hyperparameter combination based on RMSE
best_index <- which.min(results_matrixSoftBART[,3])
best_hyperparametersSoftBART <- results_matrixSoftBART[best_index, 1:2]


# Initialize a matrix to store the results for each combination
results_matrixRF <- matrix(ncol = 2, nrow = 4)
PercentOfVariable = c(0.1,0.25,0.5,1)
for (i in 1:4) {
  # Perform cross-validation using a loop
  fold_RMSE <- numeric(num_folds)
  for (j in 1:num_folds) {
    TestSetCV <- unlist(folds[-j],use.names = FALSE)
    TrainSetCV <- folds[[j]]
    
    mean.y.hat<-as.vector(randomForest(X[TrainSetCV,],Y[TrainSetCV],X[TestSetCV,],f(X[TestSetCV,]),mtry=round(PercentOfVariable[i]*length(X[1,])))$test$predicted)
    fold_RMSE[j] <- sqrt(mean((f(X[TestSetCV,])-mean.y.hat)^2))
    
  }
  # Store the average RMSE from cross-validation
  results_matrixRF[i, 1] <- mean(fold_RMSE)
  results_matrixRF[i,2] <- PercentOfVariable[i]
}
ordered_indices <- order(results_matrixRF[, 1])
ordered_results <- results_matrixRF[ordered_indices, ]
PercentOfVariable<-ordered_results[1,2]


for (l in 1:num_of_datasets){
  
  X<-BenchmarkX[[l]]
  Y<-BenchmarkY[[l]]
  
  Y<-as.numeric(as.matrix(Y))
  
  n=length(Y)
  TrainSet<-matrix(nrow=trunc(n/6),ncol = NumOfRep)
  TestSet<-matrix(nrow=n-trunc(n/6),ncol = NumOfRep)
  
  for (h in 1:NumOfRep){
    TrainSet[,h]=sort(sample.int(n,n/6));
    TestSetInit<-1:n
    TestSet[,h]=TestSetInit[! TestSetInit %in% TrainSet[,h]];
    
  } 
  
  
  Default.AddiVortes.RMSE<-foreach(k = 1:NumOfRep, .combine = cbind) %dopar% {
    library('invgamma')
    library('FNN')
    
    AddiVortes.RMSE <- AddiVortes_Algorithm(Y[TrainSet[,k]],as.matrix(X[TrainSet[,k],]),best_hyperparametersAddiVortes[1],1200,200,best_hyperparametersAddiVortes[2],best_hyperparametersAddiVortes[3],best_hyperparametersAddiVortes[4],best_hyperparametersAddiVortes[5],best_hyperparametersAddiVortes[6],best_hyperparametersAddiVortes[7],f(X[TestSet[,k],]),as.matrix(X[TestSet[,k],]))$Out_of_sample_RMSE
  }
  
  BART.RMSE<-foreach(k = 1:NumOfRep, .combine = cbind) %dopar% {
    library('BayesTree')
    library('expint')
    
    #bartMAchine
    modelBART<-bart(as.matrix(X[TrainSet[,k],]),as.numeric(as.matrix(Y[TrainSet[,k]])),as.matrix(X[TestSet[,k],]),nskip = 1000,ndpost = 3000, ntree = best_hyperparametersBART[1], sigdf = best_hyperparametersBART[2],sigquant = best_hyperparametersBART[3],k = best_hyperparametersBART[4])
    BART.RMSE<-(mean((f(X[TestSet[,k],])-modelBART$yhat.test.mean)**2))**0.5  
  }
  
  RForest.RMSE<-foreach(k = 1:NumOfRep, .combine = cbind) %dopar% {
    library('randomForest')
    #randomForests
    
    RandomF<-randomForest(X[TrainSet[,k],],Y[TrainSet[,k]],X[TestSet[,k],],f(X[TestSet[,k],]),mtry = round(PercentOfVariable*length(X[1,])))
    RForest.RMSE<-(mean((f(X[TestSet[,k],])-RandomF$test$predicted)**2))**0.5
  }
  
  # Perform cross-validation with the bart function
  Boost.RMSE<-vector(length=NumOfRep)
  Boost.RMSE<-foreach(k = 1:NumOfRep,.combine = cbind) %dopar% {
    library('gbm')
    #Gradient Boosting
    
    data_frame_Train<-as.data.frame(cbind(Y[TrainSet],X[TrainSet,]))
    data_frame_Test<-as.data.frame(cbind(f(X[TestSet,]),X[TestSet,]))
    
    gradientBoost<-gbm(V1~.,data=data_frame_Train,distribution="gaussian",interaction.depth = best_hyperparametersBoost[1],n.trees = best_hyperparametersBoost[2],shrinkage = best_hyperparametersBoost[3])
    Boost.RMSE[k]<-(mean((f(X[TestSet,])-predict(gradientBoost,newdata=data_frame_Test))**2))**0.5
  }
  
  SoftBART.RMSE<-foreach(k = 1:NumOfRep, .combine = cbind) %dopar% {
    library(SoftBart)
    
    df <- data.frame(X = X[TrainSet[,k],], Y = Y[TrainSet[,k]])
    df_test <- data.frame(X = X[TestSet[,k],], Y=f(X[TestSet[,k],]))
    
    # Train the bart model
    modelSoftBART<-softbart_regression(Y ~ ., df, df_test, num_tree = best_hyperparametersSoftBART[1],k =best_hyperparametersSoftBART[2])
    
    # Make predictions on the test 
    SoftBART.RMSE<-sqrt(mean((modelSoftBART$mu_test_mean-f(X[TestSet[,k],]))^2))
  }
  
  RMSE<-rbind(Default.AddiVortes.RMSE,BART.RMSE,RForest.RMSE,Boost.RMSE,SoftBART.RMSE)
  
  
  print('Finished')
  print(l)
  
  boxplot(RMSE[1,],RMSE[2,],RMSE[3,],RMSE[4,],RMSE[5,],horizontal = TRUE, names = unique(c("AddiVortes Default", "BART", "Random Forests","Gradient Boosting","SoftBART")))
  
  if (l==1){
    All.RMSE<-matrix(nrow=4,ncol=NumOfRep)
    All.RMSE<-RMSE
  }
  else{
    All.RMSE<-cbind(All.RMSE,RMSE)
  }
  
}

stopCluster(cl3)

for(i in 1:5){
  All.RMSE[i,]<-as.vector(All.RMSE[i,])
}

mean(All.RMSE[2,])

boxplot(All.RMSE[1,],All.RMSE[2,],All.RMSE[3,], All.RMSE[4,], All.RMSE[5,],horizontal = TRUE, names = unique(c("AddiVortes","BART","RF", "Boosting","SoftBART")))


  return(All.RMSE)
  
  }
