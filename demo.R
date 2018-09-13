dir<-"C:\\Users\\auz5836\\Documents\\GitHub\\GroupLasso\\clean"
source(sprintf("%s\\cv.R",dir))
source(sprintf("%s\\func.R",dir))
source(sprintf("%s\\opt.R",dir))
source(sprintf("%s\\sim.R",dir))


# Simulation Setup
n <- 100 # sample size
m_X <- 5 # dim of baseline covariates
m_W <- 1 # dim of treatment covariate
m_G <- 100 # dim of biomarkers
m_I <- m_G # dim of interaction terms
SNR <- 10 # signal to noise ratio
tau1 <- 1 # initial stepsize



# Generate X and Y
sigma <- cov_block(m_G,.3,20) # Covariance Matrix
x <- sim_X(m_X,m_W,m_G,sigma,n) # Continuous
# binprob <- runif(m_G) # SNP
# x <- sim_X_cate(m_X,m_W,m_G,sigma,n,binprob) # SNP
beta <- sim_beta_const(m_X,m_W=1,m_G,main_nonzero=0.1,inter_nonzero=0.1,both_nonzero=0.01,const=c(3,5),heir=TRUE)
y0 <- x%*%beta

# Noise with fixed SNR
noise <- rnorm(n,sd=1)
SNRmtl <- as.numeric(sqrt(var(y0)/(SNR*var(noise))))
y <- y0+SNRmtl*noise  
colnames(x) <- c(1:dim(x)[2])
truth <- which(beta!=0)


# Cross Validation
x0 <- rep(0,dim(x)[2])
lamb_candidate <- seq(1,20,2)
lamb_candidate2 <- seq(1,5,1)
sol_cv <- opt_lambda(x,y,f, gradf, g, proxg, x0, tau1, max_iters = 100, w = 10, 
                   backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5, 
                   eps_n = 1e-15,m_X,m_W,m_G,m_I,K=5,n=100, lamb_candidate, lamb_candidate2,restart=TRUE,beta)
lamb_opt <- get_lambda(sol_cv,log(n),2,lamb_candidate,lamb_candidate2)
lamb_opt1 <- lamb_opt$lamb_opt1
lamb_opt2 <- lamb_opt$lamb_opt2
sol <- FASTA(x,y,f, gradf, g, proxg, x0, tau1, max_iters = 300, w = 10,
            backtrack = TRUE, recordIterates = FALSE, stepsizeShrink = 0.5,
            eps_n = 1e-15,m_X,m_W,m_G,m_G,lamb_opt1,lamb_opt2,restart=TRUE)
glassorst <- which(sol$x!=0)
