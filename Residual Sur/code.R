
border <- function(alpha,x,y,M,bord.size,blur.sd)
{
	# This function adds a border to the scatterplot
	# to make Yhat and R0 orthogonol with "spacing" 
	# determined by alpha.

	# x=Yhat and y=R0 are original data.
	
	# M is nubmer of points to add to each border.

	# bord.size is % away from original max and min that 
	# the border is.
	
	# blur.sd  is sd of "blur" for border
	# x (Yhat) and Y (R0) are both centered and scaled. This seems to improve
	# numerical stability.
	
    range.x <- diff(range(x))
    range.y <- diff(range(y))

    min.x <- min(x)-bord.size*range.x
    max.x <- max(x)+bord.size*range.x
    min.y <- min(y)-bord.size*range.y
    max.y <- max(y)+bord.size*range.y

    range.y <- max.y-min.y
    range.x <- max.x-min.x
    
    u.a <- seq(from=0,to=1,length=M)^alpha

    left.border.x <- rnorm(M,mean=min.x,sd=blur.sd)
    left.border.y <- min.y+(range.y)*u.a
               
    right.border.x <- rnorm(M,mean=max.x,sd=blur.sd)
    right.border.y <- min.y+(range.y)*(1-u.a)
    
    bottom.border.y <- rnorm(M,mean=min.y,sd=blur.sd)
    bottom.border.x <- min.x+(range.x)*u.a
    
    top.border.y <- rnorm(M,mean=max.y,sd=blur.sd)
    top.border.x <- min.x+(range.x)*(1-u.a)
    
    border.x <- c(bottom.border.x,top.border.x,left.border.x,right.border.x)
    border.y <- c(bottom.border.y,top.border.y,left.border.y,right.border.y)

    x <- c(x,border.x)
    y <- c(y,border.y)
    
    
    list(x=as.vector(scale(x)),y=as.vector(scale(y)))
}


slope <- function(alpha,x,y,M,bord.size,blur.sd)
  {
  	# internal fuction for finding optimal alpha
    temp <- border(alpha,x,y,M,bord.size,blur.sd)
    x <- temp$x
    y <- temp$y
    as.vector(lm(y~x)$coef[2])
  }

find.bord <- function(x,y,M=100,bord.size=.05,min.a=0.0001,max.a=100, blur.sd=0)
  {
  	# function to add a border to Yhat and R0
  	
    alpha <- uniroot(slope,c(min.a,max.a), 
    tol = 10e-13,x=x,y=y,M=M,bord.size=bord.size,blur.sd=blur.sd)
    if (abs(alpha$f.root)>(.01)) print("problem: larger M? bigger range for alpha?")
    border(alpha$root,x,y,M,bord.size,blur.sd)
  }



find.regression <- function(Yhat,R0,coef.det,p,j=1,eps=10e-13)
  {
  	# start with orthogonal Yhat and R0 and find an X (matrix) and Y (vector) pair.
  	# coef.det is coefficient of determination
  	# p is number of non-intercept covariates
  	# j is the column of M to use for iterations
  	# eps is a convergence criterium
  	
    sdYhat <- sd(Yhat)
    sdR0 <- sd(R0)

    n <- length(R0)
    CD <- coef.det
    beta0 <- 1
    betas <- runif(p)
    j <- 1

    Yhat <- (sdR0/sdYhat)*sqrt(CD/(1-CD))*Yhat
    
    PR0 <- outer(R0,R0,"*")/sum(R0^2)
    Z <- rnorm(n,mean=0,sd=sdR0)
    M <- matrix(rnorm(p*n,mean=0,sd=sdYhat),n,p)
    test <- 1
    while (test>eps)
      {
        W <- cbind(rep(1,n),(diag(rep(1,n))-PR0)%*%M)
        A <- W%*%solve(t(W)%*%W)%*%t(W)
        
        term.1 <- A%*%Z
        term.2 <- PR0%*%M%*%betas
        term.3 <- apply(t(t(M)*betas)[,-j],1,sum)
        temp <- (1/betas[j])*(Yhat-beta0-term.1+term.2-term.3)
        
        new.M <- M
        new.M[,j] <- temp
        
        delta <- (diag(rep(1,n))-PR0)%*%((new.M-M))
        test <- max(abs(delta))
        
        M <- new.M
      }
    W <- cbind(rep(1,n),(diag(rep(1,n))-PR0)%*%M)
    A <- W%*%solve(t(W)%*%W)%*%t(W)
    
    epsilon <- as.vector(R0+A%*%Z)
    
    X <- (diag(rep(1,n))-PR0)%*%M
    Y <- beta0+X%*%betas+epsilon
    
    data.frame(cbind(Y=Y,X=X))
  }



temp <- scan("yufree.pbm",skip = 1)
# yufree.pbm is a bitmap file (from R's files) that was saved
# as a bitmap text file (bitsize=1) by imagemagick
# try convert 'fig.jpg'  -resize 100x100 -extent 100x100 -monochrome -compress none 'fig.pbm'

# The code below makes it into "scatterplot" format.
y <- rep(1:temp[2],each=temp[1])
x <- rep(1:temp[1],times=temp[2])
y <- -y[temp[-(1:2)]==1]
x <- x[temp[-(1:2)]==1]

# Augment the data to make x and y orthogonal.

picture <- find.bord(x,y,blur.sd=0)
Yhat <- picture$x
R0 <- picture$y

data <- find.regression(Yhat,R0,coef.det=.05,p=5,j=1,eps=10e-13)
# data$X has covariates (without intercept)
# data$Y has response.

save(data,file = 'data.RData')