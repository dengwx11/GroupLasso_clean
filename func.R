############# differentiable part in loss function #############
#############           Gradient                   #############
# In this differentiable part, we combine L2 loss function for linear model 
# and ridge term into f(\theta) together
# Because they are both differentiable

## Logistic family ##
f0 <- function(beta,X,y) { -t(y)%*%(X%*%beta) + sum(log(1+exp(X%*%beta))) } # objective function
gradf0 <- function(beta,X,y) { -t(X)%*%(y-plogis(X%*%beta)) } # gradient
f <- function(beta,X,y,m_X,m_W,m_G,m_I,lambda2){
 f0(beta,X,y)+lambda2*norm(rep(0:1,c(m_X+m_W,m_G*2))*beta,'2')
}
gradf <- function(beta,X,y,m_X,m_W,m_G,m_I,lambda2) {
 gradf0(beta,X,y)+as.matrix(rep(0:1,c(m_X+m_W,m_G*2))*beta,nrow=m_X+m_W+m_G+m_I) 
}

## Gaussian family ##
f0 <- function(beta,X,y) { 0.5*norm(X%*%beta - y, "F")^2}  # objective function
gradf0 <- function(beta,X,y) { t(X)%*%(X%*%beta - y)  } # f0 gradient
f <- function(beta,X,y,m_X,m_W,m_G,m_I,lambda2) { 0.5*norm(X%*%beta - y, "F")^2+
    lambda2*norm(rep(0:1,c(m_X+m_W,m_G*2))*beta,'2')} # objective function + elastic net
gradf <- function(beta,X,y,m_X,m_W,m_G,m_I,lambda2) { t(X)%*%(X%*%beta - y) + 
    lambda2*as.matrix(rep(0:1,c(m_X+m_W,m_G*2))*beta,nrow=m_X+m_W+m_G+m_I) } # gradient f


############# Nondifferentiable part in loss function #############
#############              Subgradient                #############
# In this nondifferentiable part, we only include group lasso structure as g(\theta)

### Group Lasso ###

split_beta <- function(beta,m_X,m_W,m_G,m_I) {
  # Split coefficients into list of different parts given dimensions of X,W,G
  #
  # Args:
  #   beta =  vector of coefficients
  #   m_X = number of baseline covariates
  #   m_W = number of treatment covariates. Default is 1
  #   m_G = number of biomarkers
  #   m_I = number of biomarkers for interaction term
  #
  # Returns:
  #   A user-friendly version for beta after spliting different types of covariates
  
  if (m_X != 0 && m_W != 0) {
    ma <- rep(1:4,c(m_X,m_W,m_G,m_I))
    beta <- split(beta,ma)
    names(beta) <- c("X","W","G","I")
  } else if (m_X == 0 && m_W == 0) {
    ma <- rep(1:2,c(m_G,m_I))
    beta <- split(beta,ma)
    names(beta) <- c("G","I")
  } else if (m_X==0) {
    ma <- rep(1:3,c(m_W,m_G,m_I))
    beta <- split(beta,ma)
    names(beta) <- c("W","G","I")
  } else {
    ma <- rep(1:3,c(m_X,m_G,m_I))
    beta <- split(beta,ma)
    names(beta) <- c("X","G","I")
  }
  return(beta)
}

split_X <- function(X,m_X,m_W,m_G,m_I) {
  # Split covariates into list of different parts given dimensions of X,W,G
  #
  # Args:
  #   X =  design matrix
  #   m_X = number of baseline covariates
  #   m_W = number of treatment covariates. Default is 1
  #   m_G = number of biomarkers
  #   m_I = number of biomarkers for interaction term
  #
  # Returns:
  #   A user-friendly version for X after spliting different types of covariates
  
  ma <- rep(1:(2+m_G),c(m_X,m_W,rep(1,m_G)))
  ma <- c(ma,rep((2+1):(2+m_G),rep(1,m_I)))
  X_split <- split(X,ma)
  X_split <- lapply(X_split,function(x) x=matrix(x,nrow=n))
  return(X_split)
}

group_penalty <- function(X,m_X,m_W,m_G,m_I){
  # Penalty parameters for main and interation terms
  #
  # Args:
  #   X =  design matrix
  #   m_X = number of baseline covariates
  #   m_W = number of treatment covariates. Default is 1
  #   m_G = number of biomarkers
  #   m_I = number of biomarkers for interaction term
  #
  # Returns:
  #   Feature weights
  
  X_split <- split_X(X,m_X,m_W,m_G,m_I)
  if (m_X != 0 && m_W != 0) { ## If we have both baseline and treatment covariates
  para_I <- sapply(X_split[-1:-2], function(x) norm(as.matrix(x[,2]),'2'))
  para_mI <- sapply(X_split[-1:-2], function(x) norm(cbind(x[,1],sqrt(1.4*(1-sqrt(2/pi)))*x[,2]),'F'))
  } else if ( m_X == 0 && m_W != 0 ) { ## If we only have treatment covariate
    para_I <- sapply(X_split[-1], function(x) norm(as.matrix(x[,2]),'2'))
    para_mI <- sapply(X_split[-1], function(x) norm(cbind(x[,1],sqrt(0.64-sqrt(2/pi)*exp(-0.5))*x[,2]),'F'))
  }
  para <- list("I"=para_I,"MI"=para_mI) # I=feature weights of Interaction terms, MI=feature weights for Main and Interaction terms
  return(para)
}

g <- function(X,beta,m_X,m_W,m_G,m_I,lambda) {
  # group lasso term
  #
  # Args:
  #   X =  design matrix
  #   beta = vector of beta
  #   m_X = number of baseline covariates
  #   m_W = number of treatment covariates. Default is 1
  #   m_G = number of biomarkers
  #   m_I = number of biomarkers for interaction term
  #   lambda = regularization parameter for group lasso structure
  #
  # Returns:
  #   g() value
  
  beta <- split_beta(beta,m_X,m_W,m_G,m_I)
  para <- group_penalty(X,m_X,m_W,m_G,m_I)
  penalty <- lambda*norm(as.matrix(para$I*beta$I),'1')+lambda*sum(para$mI*sqrt(beta$G^2+beta$I^2))
  return(penalty)
}

proxg <- function(X,beta,m_X,m_W,m_G,m_I,tau,lambda) { 
  # Proximal operator for group lasso term
  #
  # Args:
  #   X =  design matrix
  #   beta = vector of beta
  #   m_X = number of baseline covariates
  #   m_W = number of treatment covariates. Default is 1
  #   m_G = number of biomarkers
  #   m_I = number of biomarkers for interaction term
  #   tau =  stepsize in proximal operator
  #   lambda = regularization parameter for group lasso structure
  #
  # Returns:
  #   minimizers in proximal operator

  beta <- split_beta(beta,m_X,m_W,m_G,m_I) ## results for baseline and treatment keep the same
  para <- group_penalty(X,m_X,m_W,m_G,m_I)
  
  beta_temp <- cbind(beta$G,beta$I)
  beta_temp <- split(beta_temp,rep(1:m_G,rep(1,m_G)))

  beta_temp <- mapply(FUN=function(x,y,z) 
  { max(0,1-lambda*tau*z/(x[1]^2+(sign(x[2])*max(abs(x[2])-lambda*tau*y,0))^2)^0.5) * c(x[1],sign(x[2])*max(abs(x[2])-lambda*tau*y,0))},beta_temp,para$I,para$MI)
  beta$G <- beta_temp[1,] ## results for prognostic effect coefficient
  beta$I <- beta_temp[2,] ## results for predictive effect coefficient
  

  
  beta <- unlist(beta,use.names = F)
  return(beta)
}

### Regular Lasso ###

glasso<- function(X,beta,m_X,m_W,m_G,m_I,lambda) {
  # Regular Lasso Term
  #
  # Args:
  #   X =  design matrix
  #   beta = vector of beta
  #   m_X = number of baseline covariates
  #   m_W = number of treatment covariates. Default is 1
  #   m_G = number of biomarkers
  #   m_I = number of biomarkers for interaction term
  #   lambda = regularization parameter for group lasso structure
  #
  # Returns:
  #   g() value
  
  beta<-split_beta(beta,m_X,m_W,m_G,m_I)
  penalty<-lambda*norm(as.matrix(beta$I),'1')+lambda*norm(as.matrix(beta$G),'1')
  return(penalty)
}
proxglasso<- function(X,beta,m_X,m_W,m_G,m_I,tau,lambda){
  # Proximal operator for Lasso Term
  #
  # Args:
  #   X =  design matrix
  #   beta = vector of beta
  #   m_X = number of baseline covariates
  #   m_W = number of treatment covariates. Default is 1
  #   m_G = number of biomarkers
  #   m_I = number of biomarkers for interaction term
  #   tau =  stepsize in proximal operator
  #   lambda = regularization parameter for group lasso structure
  #
  # Returns:
  #   minimizers in proximal operator
  
  beta<-split_beta(beta,m_X,m_W,m_G,m_I)
  beta$G<-sapply(beta$G ,function(x) sign(x)*max(0,abs(x)-lambda*tau))
  beta$I<-sapply(beta$I, function(x) sign(x)*max(0,abs(x)-lambda*tau))
  beta<-unlist(beta,use.names = F)
  return(beta)
}

