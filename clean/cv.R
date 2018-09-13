########################## Cross Validation ##########################


divide<-function(K,n){
  # Divide samples into K folds randomly
  # 
  # Args:
  #   K = number of folds
  #   n = sample size
  #
  # Returns:
  #   folds = list of indices
  
  i.mix = sample(1:n)
  folds = vector(mode="list",length=K)
  temp<-split(c(1:n),1:K)
  for(i in 1:K){
    folds[[i]]<-i.mix[temp[[i]]]
  }
  return(folds)
}

cv.FASTA<-function(X,y,f, gradf, g, proxg, x0, tau1, max_iters = 100, w = 10, 
                   backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                   eps_n = 1e-15,m_X,m_W,m_G,m_I,lambda,lambda2,K=10,n=100,restart,truth){
  # Summary of MSE, prediction errors, and model size in cross validation with different regularization parameters
  #
  # Args:
  #   X = design matrix
  #   y = vector of responses
  #   f = differentialble part in loss function
  #   gradf = gradient of f
  #   g = nondifferentialble part in loss function
  #   proxg = proximal operator for g
  #   x0 = start point
  #   tau1 =  initial stepsize
  #   max_iters = maximum iterations before automatic termination
  #   w = lookback window for non-montone line search
  #   backtrack =  logical; if TRUE, perform backtracking line search
  #   recordIterates = logical; If TRUE, record iterate sequence
  #   stepsizeShrink = multplier to decrease step size
  #   eps_n = epsilon to prevent normalized residual from dividing by zero
  #   m_X = number of baseline covariates
  #   m_W = number of treatment covariates. Default is 1
  #   m_G = number of biomarkers
  #   m_I = number of interaction terms, should equals m_G
  #   K = K-folds. Default is 10
  #   n = sample size
  #   lambda_candidate = list of regularization parameters for group lasso structure
  #   lambda_candidate2 = list of regularization parameters for ridge structure
  #   restart = logical; if TRUE (default), perform adaptive restart
  #   truth = vector of true beta
  #
  # Returns:
  #   mean_pred = Mean of prediction error in K-fold 
  #   MSE_beta = Mean of MSE in K-fold
  #   SD_pred = SD of prediction error in K-fold
  #   Var_num = Variance of estimated model size in K-fold
  #   mean_num = Mean of estimated model size in K-fold
  
  folds<-divide(K,n)
  
  sol_cv<-list()
  TestErr_pred<-rep(0,K)
  TestErr_beta<-rep(0,K)
  num<-rep(0,K)
  
  for(k in 1:K){
    
    # Generating training and test datasets for cross validation
    test_X<-X[folds[[k]],]
    train_X<-X[setdiff(c(1:n),folds[[k]]),]
    test_y<-y[folds[[k]]]
    train_y<-y[setdiff(c(1:n),folds[[k]])]
    
    # Training model on training dataset
    
    sol_cv[[k]]<-FASTA(train_X,train_y,f, gradf, g, proxg, x0, tau1, max_iters = 100, w = 10, 
                        backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                        eps_n = 1e-15,m_X,m_W,m_G,m_I,lambda,lambda2,restart)
    x0<-sol_cv[[k]]$x
    
    # Test on the validation group
    TestErr_pred[k]<-norm(f0(sol_cv[[k]]$x,test_X,test_y)-f0(truth,test_X,test_y),"2")/length(test_y) # Prediction Errors in Cross Validation
    TestErr_beta[k]<-norm(sol_cv[[k]]$x-truth,"2")^2 # MSE
    num[k]<-length(which(sol_cv[[k]]$x!=0)) # estimated model size
    
  }
  return(list(Err_pred=TestErr_pred,Err_beta=TestErr_beta,start=x0,num=num))
}


opt_lambda<-function(X,y,f, gradf, g, proxg, x0, tau1, max_iters = 100, w = 10, 
                      backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                      eps_n = 1e-15,m_X,m_W,m_G,m_I,K=10,n=100,lamb_candidate,lamb_candidate2,restart,truth){
  # Summary of MSE, prediction errors, and model size in cross validation with different regularization parameters
  #
  # Args:
  #   X = design matrix
  #   y = vector of responses
  #   f = differentialble part in loss function
  #   gradf = gradient of f
  #   g = nondifferentialble part in loss function
  #   proxg = proximal operator for g
  #   x0 = start point
  #   tau1 =  initial stepsize
  #   max_iters = maximum iterations before automatic termination
  #   w = lookback window for non-montone line search
  #   backtrack =  logical; if TRUE, perform backtracking line search
  #   recordIterates = logical; If TRUE, record iterate sequence
  #   stepsizeShrink = multplier to decrease step size
  #   eps_n = epsilon to prevent normalized residual from dividing by zero
  #   m_X = number of baseline covariates
  #   m_W = number of treatment covariates. Default is 1
  #   m_G = number of biomarkers
  #   m_I = number of interaction terms, should equals m_G
  #   K = K-folds. Default is 10
  #   n = sample size
  #   lambda_candidate = list of regularization parameters for group lasso structure
  #   lambda_candidate2 = list of regularization parameters for ridge structure
  #   restart = logical; if TRUE (default), perform adaptive restart
  #   truth = vector of true beta
  #
  # Returns:
  #   mean_pred = Mean of prediction error in K-fold 
  #   MSE_beta = Mean of MSE in K-fold
  #   SD_pred = SD of prediction error in K-fold
  #   Var_num = Variance of estimated model size in K-fold
  #   mean_num = Mean of estimated model size in K-fold

  
  lamb_candidate<-lamb_candidate[order(lamb_candidate,decreasing = T)]
  lamb_candidate2<-lamb_candidate2[order(lamb_candidate2,decreasing = T)]
  
  
  TestErr_pred<-matrix(0,nrow=length(lamb_candidate),ncol=length(lamb_candidate2)) 
  MSE_beta<-matrix(0,nrow=length(lamb_candidate),ncol=length(lamb_candidate2)) 
  SDErr_pred<-matrix(0,nrow=length(lamb_candidate),ncol=length(lamb_candidate2)) 
  Mean_num<-matrix(0,nrow=length(lamb_candidate),ncol=length(lamb_candidate2)) 
  Var_num<-matrix(0,nrow=length(lamb_candidate),ncol=length(lamb_candidate2))
  
  for(i in seq_along(lamb_candidate)){
    for(j in seq_along(lamb_candidate2)){
      print(c(i,j))
    rst<-cv.FASTA(X,y,f, gradf, g, proxg, x0, tau1, max_iters = max_iters, w = 10, 
                   backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                   eps_n = 1e-15,m_X,m_W,m_G,m_I,lamb_candidate[i],lamb_candidate2[j],K,n,restart=TRUE,truth)
    cv.Err_pred<-rst$Err_pred
    cv.Err_beta<-rst$Err_beta
    cv.num<-rst$num
    x0<-rst$start
    TestErr_pred[i,j]<-mean(cv.Err_pred)
    MSE_beta[i,j]<-mean(cv.Err_beta)
    Mean_num[i,j]<-mean(cv.num)
    SDErr_pred[i,j]<-sqrt(var(cv.Err_pred)/K)
    Var_num[i,j]<-var(cv.num)
    print(c(paste("lambda 1=",lamb_candidate[i]),paste("lambda 2=",lamb_candidate2[j])))
    print(cv.Err_pred)
    print(cv.Err_beta)
    }
    TestErr_pred[i,]<-rev(TestErr_pred[i,])
    MSE_beta[i,]<-rev(MSE_beta[i,])
    SDErr_pred[i,]<-rev(SDErr_pred[i,])
    Mean_num[i,]<-rev(Mean_num[i,])
    Var_num[i,]<-rev(Var_num[i,])
  }
  
  TestErr_pred<-apply(TestErr_pred,2,rev)
  MSE_beta<-apply(MSE_beta,2,rev)
  SDErr_pred<-apply(SDErr_pred,2,rev)
  Mean_num<-apply(Mean_num,2,rev)
  Var_num<-apply(Var_num,2,rev)
  
  return(list(mean_pred=TestErr_pred,MSE_beta=MSE_beta,SD_pred=SDErr_pred,Var_num=Var_num,mean_num=Mean_num))
}


get_lambda<-function(sol_cv_glasso,k_mean,k_var=2,lamb_candidate,lamb_candidate2){
  # Get the optimal lambda and lambda2
  # 
  # Args:
  #   sol_cv_glasso = results from opt_lambda()
  #   k_mean = weight for mean_num
  #   k_var =  weight for var_num
  #   lambda_candidate = list of regularization parameters for group lasso structure
  #   lambda_candidate2 = list of regularization parameters for ridge structure
  #
  # Returns:
  #   lamb_opt1 = lambda1 optimal value
  #   lamb_opt2 = lambda2 optimal value
  
  rst<-sol_cv_glasso$mean_pred+k_mean*sol_cv_glasso$mean_num+k_var*sol_cv_glasso$Var_num # Could be customized
  lamb_loc <- which(rst==min(rst),arr.ind = T)
  lamb_opt1<-lamb_candidate[lamb_loc[1]]
  lamb_opt2<-lamb_candidate2[lamb_loc[2]]
  return(list(lamb_opt1=lamb_opt1,lamb_opt2=lamb_opt2))
}

