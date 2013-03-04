##############################################################################################
## Programs Name:  FrSVM.R
## TODO:  	   An R algorithm, called FrSVM, which integrates protein-protein 
##		 		   interaction network information into gene selection for 
##		 		   microarry classification
##
## Autor: 		   Yupeng Cun @ B-IT, Uni-Bonn
## Contact: 	   yupeng.cun
## Create Dtate:   12/12/2011 
## Lastest Update: 14/12/2012
## Refer papers:   Yupeng Cun, Holger Fröhlich (2012) Integrating Prior Knowledge Into 
##					Prognostic Biomarker Discovery Based on Network Structure.arXiv:1212.3214 
##
## Parameters---- 
##   		x: gene expression data
##   		y: class labels
##   		d: damping factor for GeneRank, defaults value is 0.5
## 		 Gsub: Adjacency matrix of Protein-protein intersction network
##   	folds: # of -folds cross validation (CV)
## 	  repeats: $ of CV repeat times
## 	 parallel: paralle computing or not
##   	cores: cores used in parallel computing
##   	DEBUG: show more results or not
##   top.uper: the uper bound of top ranked genes
##  top.lower: the lower bound of top ranked genes
##
## Returned results---- 
##		  auc: AUC values of each test fold in CV
##     labels: original class labels
##   	 fits: SVM model for each training fold
## 		 feat: Selected features in each training folds
##############################################################################################

library(ROCR)
library(Matrix)
library(kernlab)

# source("../FrSVM.R")
# res <- frSVM.cv(x, yy, folds=10,Gsub=Gsub, repeats=10, parallel = TRUE, cores = 4, DEBUG=TRUE,d=0.85,top.uper=0.95,top.lower=0.9)
###
## fnetSVM: filter feature based on network property 
##
##

frSVM.cv <- function(x, y, folds=10,Gsub=adjacency.matrix, repeats=5, parallel = TRUE, cores = 2, DEBUG=TRUE,d=0.85, top.uper=0.95,top.lower=0.9)
{
	multicore <- ("package:multicore" %in% search())
  
	if(multicore == TRUE && parallel == TRUE)  
	{
    
		if(is.null(cores))       
			cores <- multicore:::detectCores()    
		options(cores = cores - 1)     
		
		cat("Detected ", cores," cores. Will use ", getOption("cores"), " of them.\n")    
		parallel <- TRUE  
	}  
	else  
	{
		if(parallel == TRUE)       
		cat('Package \'multicore\' not loaded. Please, load it manually prior to calling this function if you want to run classification in parallel.\n',sep='')
		cat('Will continue with sequential crossvalidation.\n', sep='')
		parallel <- FALSE
	}

	if(!is.factor(y)) stop("y must be factor!\n")
  
	if(length(levels(y)) != 2) stop('y must be factor with 2 levels.\n')
	n     <- length(y)
	folds <- trunc(folds)
  
	if(length(y) != nrow(x)) stop('y must have same length as nrow(x).\n')
	if (folds < 2) stop("folds should be greater than or equal to 2.\n")
	if (folds > n) stop("folds should be less than or equal to the number of observations.\n")

	cuts  <- cv.repeats <- list()  	  
	op =	top.uper 
	aa =    top.lower 
	
	
			
	set.seed(1234)	
	for(r in 1:repeats)
	{
		perm = sample(1:n)		
		#perm <- sample(1:n) #Sampling a random integer between 1:n
		repeat.models <- NULL 
    
		for(k in 1:folds) #randomly divide the training set in to 10 folds
		{
			tst <- perm[seq(k, n, by=folds)]  #      
			trn <- setdiff(1:n, tst)            
			cuts[[k]] <- list(trn=trn, tst=tst)    
		}    	
    
		if(DEBUG) cat('Starting classification of repeat:',r,'\n')
    
		if(parallel)	repeat.models <- mclapply(1:folds, classify.fnet, cuts=cuts, x=x, y=y,cv.repeat=r, DEBUG=DEBUG,Gsub=Gsub,d=d,op=op,aa=aa)
		else		repeat.models <-   lapply(1:folds, classify.fnet, cuts=cuts, x=x, y=y, cv.repeat=r, DEBUG=DEBUG, Gsub=Gsub,d=d,op=op,aa=aa)
	
		
    
		if(length(repeat.models) != folds)
		{      
			geterrmessage()      
			stop("One or more processes did not return. May be due to lack of memory.\n")
		}
	
		if(DEBUG)	cat('All models of repeat:',r,'have been trained.\n')
    
		cv.repeats[[r]] <- repeat.models
  
  	} 
	 
  	auc <- sapply(cv.repeats, function(cv.repeat) sapply(cv.repeat, function(model) model$auc))  
  	colnames(auc) <- paste("Repeat",1:repeats,sep="")  
  	rownames(auc) <- paste("Fold",1:folds,sep="")    	
  	
  	fits <- lapply(cv.repeats,function(cv.repeat) lapply(cv.repeat, function(model) model$model))  
  	names(fits) <- paste("Repeat",1:repeats,sep="")  
  	fits <- lapply(fits, function(x)  {names(x) = paste("Fold", 1:folds, sep = ""); x })  
  	
  	feat <- lapply(cv.repeats,function(cv.repeat) lapply(cv.repeat, function(model) model$feat))  
  	names(feat) <- paste("Repeat",1:repeats,sep="")  
  	feat <- lapply(feat, function(x)  {names(x) = paste("Fold", 1:folds, sep = ""); x })
  	
  	res <- list(feat=feat, auc=auc,fits=fits, labels=y)  
  	class(res) <- 'pathClassResult' 
  	
  	return(res)  
	#return ( cv.repeats)
}

classify.fnet <- function(fold, cuts, x, y, cv.repeat, DEBUG=DEBUG,Gsub=Gsub,d=d,op=op,aa=aa)
{
	gc() 	
	if(DEBUG) cat('starting Fold:',fold,'\n')	
  	#ad.list<- as.adjacencyList(Gsub)
	## get training and test indices
	trn <- cuts[[fold]]$trn
	tst <- cuts[[fold]]$tst		

	label 	<- sign(as.numeric(y[tst]) - 1.5) # for computing the AUC
	 	    	  	  
	fits  <- list()
  	best.boundL <- list()
  	featL  <- list()
 
	if(DEBUG) cat("Geting Gene Ranking \n")			
	ranks = getRanking(x=x[trn,], y=y[trn], Gsub=Gsub, d=d)	
		
	topRanks=ranks[which(ranks > quantile(ranks,aa))]
	topRanks2=ranks[which(ranks > quantile(ranks,op))]		
		
	ranksI=sort(topRanks,decreasing=T)
	#nodesI=sort(topNodes,decreasing=T)
	nn = length(ranksI)
	nna =length(topRanks2)-1
		
	if(DEBUG)cat('selected features(max): ', nn,'\n')
	if(DEBUG)cat('selected features(min): ', nna,'\n')
		
	i=1		
	while(nn > nna)
	{ 		  				  		
			feat.rank = names(ranksI[1:nn])
	  		#feat.node = names(topScorces)
	  		#feat = intersect(feat.rank,feat.node)
	  		feat = feat.rank
	  		featL[[i]]= feat
			xtt=x[trn,feat] 
		  	print(length(feat))	
	  		fit <- svm.fit(x=xtt, y=y[trn], Cs=10^c(-3:3), scale="scale", DEBUG=FALSE) 	  		  
	  		fits[[i]]=fit	  		  		
		    best.boundL[[i]] <-  fits[[i]]$error.bound  
			nn=nn-50
			i=i+1
	}	
		
	if(DEBUG)cat("the opitimal steps: ", i-1, "\n")
		
	best.boundLs= unlist(best.boundL)
	best.index = which(best.boundLs==min(best.boundLs))		
	n=length(best.index)						
			
	train   =	fits[[best.index[n]]]		
	feat = featL[[n]]	
	#train$w[feat]=train$w[feat]/rank(ranks[feat])
		
	xts= x[tst,feat]				
	test1    <- svm.predict(fit=train, newdata=xts, type="response")
	#print("test")
	#print(t(test1))
	## calculate the AUC					
	acutt	<- calc.auc(test1, label)
	auc     <- acutt$auc						

	
	if(DEBUG) 
	{
			cat("\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n")
			cat("=> the best AUC   is ", auc, "\t Best features length:  ", length(feat),"\n")		
			cat("===================================================\n")
	}
	if(DEBUG) {cat('Finished fold:',fold,'\n\n')}

	gc()	
	res=list(fold=fold, model=train, auc=auc,feat= feat)
		
	return(res)
}


getRanking = function(x=x, y=y,Gsub=Gsub, d=d)
{
	int= intersect(colnames(x),colnames(Gsub))
    x=x[,int]
    Gsub=Gsub[int,int]     
    x = scale(x)
    #calculate the t-scorce of each probe.
    xtt = matrix(0, ncol(x))
    yy=sign(as.numeric(y)-1.5)

    for(i in 1: ncol(x))
    	xtt[i]= abs(as.numeric(t.test(x[,i], yy,paired=T)$statistic))
    
    names(xtt)= colnames(x)	
    exprs = xtt[,1]      	
    names(exprs) = colnames(x)
    ranks   <- geneRank(W=Gsub, ex=exprs, d=d)   
    names(ranks) <- colnames(Gsub)  	
    return(ranks)
}


svm.fit = function(x, y, Cs, scale, DEBUG=FALSE){

  if(missing(x))     stop('No epxression matrix provided.')
  if(missing(y))     stop('No class-lables provided.')
  if(missing(Cs))    stop('No tuning parameter \'C\' provided.')
  if(missing(scale)) stop('Parameter \'scale\' must be in \'scale\', \'center\' or NULL.')
  if(length(levels(factor(y))) != 2) warning('y must have 2 levels.')
    
  scale.mean <- scale.std <- NULL

  if(!is.null(scale)){
    scale <- tolower(scale)

    if("center" %in% scale){
      x = scale(x,center=T)
      ## save centering coefficient
      scale.mean = attr(x,"scaled:center")
      names(scale.mean) = colnames(x)
    }
    
    if("scale" %in% scale){
      x = scale(x,center=F,scale=T)
      ## save scaling coefficient    
      scale.std = attr(x,"scaled:scale")
      names(scale.std) = colnames(x)
    }
  }

  ## this happens sometimes whenn all
  ## probe sets have exactely the same value
  ## and after centering everything is zero then.
  ## Due to this the scaling fails and produces
  ## NaNs:
  ## setting NaN columns to 0
  nan.cols     <- apply(x,2,function(y) all(is.nan(y)))
  x[,nan.cols] <- 0
  
  K <- kernelMatrix(vanilladot(), x)
  best.bound <- Inf

  for(C in Cs){
    if(DEBUG) cat('Trying C=',C,'\n')
    K2 = as.kernelMatrix(K + 1/C*diag(NROW(x)))
    fit.tmp = ksvm(K2, y, C=Inf, type="C-svc", shrinking=FALSE, tol=0.01,scaled=F)
    bound = spanbound(fit.tmp, K2, sign(as.numeric(y) - 1.5))
    if(bound < best.bound){
      model = fit.tmp
      best.bound = bound
      Cbest = C
    }
  }
  if(DEBUG) cat('Best C=',Cbest,'\n')
  
  ## alphaindex: The index of the resulting support vectors in the data matrix
  svs <- unlist(alphaindex(model))

  ## coef: The corresponding coefficients times the training labels.
  w <- abs(t(unlist(coef(model))) %*% x[svs,])
  
  
  fit <- list(fit.svm=model,  w=w, K=K, C=Cbest, xsvs=x[svs,,drop=FALSE], error.bound=best.bound, scale.mean=scale.mean, scale.std=scale.std, features=colnames(x), R=NULL)	
  return(fit)
}

svm.predict = function(fit, newdata, type="response"){

  ## do the prediction only with those genes
  ## that were use for training
  newdata <- newdata[,fit$features]

  if(!is.null(fit$scale.mean))
    newdata <- scale(newdata, center=fit$scale.mean[fit$features], scale=FALSE)
  
  if(!is.null(fit$scale.std))
    newdata <- scale(newdata, center=FALSE, scale=fit$scale.std[fit$features])
  
  Ktst        <- kernelMatrix(vanilladot(), newdata, fit$xsvs[,fit$features])
  Ktst2       <- kernelMatrix(rbfdot(sigma=0.001), newdata, fit$xsvs[,fit$features])
  ident       <- which(Ktst2 == 1)		
  Ktst[ident] <- Ktst[ident] + 1/fit$C		
  alpha       <- as.matrix(unlist(coef(fit$fit.svm)))	
  yhat        <- Ktst%*%alpha - b(fit$fit.svm)

  if(type == "class") yhat <- sign(yhat)	

  return(yhat)
}


spanbound <- function(fit, xtrn, ytrn){
  
  svindex = unlist(alphaindex(fit))	
  alpha = unlist(coef(fit))	
  pos = which(alpha > 0)
  neg = which(alpha < 0)
  alpha = abs(alpha)	
  if(class(xtrn) != "kernelMatrix"){
    yhat = predict(fit, xtrn, type="decision")
    if(param(fit)$C != Inf)
      error("span bound is only for L2-SVM!")
    K = kernelMatrix(kernelf(fit), xtrn)
  }
  else{
    K = xtrn
    yhat = K[,svindex]%*%as.matrix(unlist(coef(fit))) - b(fit)
  }
  output = ytrn*yhat	
  Cpos = Inf
  Cneg = Inf
  eps = 1e-5	
  boundpos = (alpha[pos] >= Cpos*(1-eps))
  boundneg = (alpha[neg] >= Cneg*(1-eps))
  sv1pos = svindex[pos[!boundpos]]
  sv2pos = svindex[pos[boundpos]]
  sv1neg = svindex[neg[!boundneg]]
  sv2neg = svindex[neg[boundneg]]	
  sv1 = sort(c(sv1pos, sv1neg))
  sv2 = sort(c(sv2pos, sv2neg))
  n = ncol(K)
  span = double(n)
  alpha1 = double(n)
  alpha1[svindex] = alpha 		
  if(length(sv1) > 0){ # in-bound SVs 								
    ell = length(sv1)	
    invK = chol2inv(chol(K[sv1,sv1,drop=FALSE]))
    T = -1/sum(invK)
    T2 = invK%*%as.matrix(rep(1,ell))*T
    T3 = t(as.matrix(rep(1,ell)))%*%invK
    invKSV = rbind(cbind(invK + T2%*%T3, -T2), cbind(-T*T3, T))
    tmp = diag(as.matrix(invKSV)) + 1e-10
    span[sv1] = 1./tmp[1:ell]
  }	
  else
    warning("No in-bound SVs!")			
  if(length(sv2) > 0){	# bound SVs	
    span[sv2] = diag(as.matrix(K[sv2,sv2,drop=FALSE]))
    if(length(sv1) > 0){
      V = rbind(K[sv1,sv2,drop=FALSE], rep(1,length(sv2)))			
      span[sv2] = span[sv2] - diag(t(V)%*%invKSV%*%V)			
    }				 
  }	
  loo = mean((output - alpha1*span <= 0)*1)				
  ## cat("Span bound =", loo,"\n\n")		
  loo
}



calc.auc <- function(prob,labels)
{
	
  ## this corrects a bug in ROCR:
  ## if all labels are from one group and there
  ## is no missclassification, ROCR is not able
  ## to calculate the auc
  ## patch => add a artificial prediction with prob = 0

  if(length(unique(labels)) == 1)
  {
    if(sign(labels[1]) == -1)
      labels <- c(labels,1)
    else
      labels <- c(labels,-1)
    prob <- c(prob,0)
  }

	tab.classes<- NULL
	sensitivity.classes<- NA
	specificity.classes<- NA
	labels.universe=NULL
	#pred.class = sign(prob)	
	
	pred.class = factor(2*as.numeric (prob > 0) -1)
	
	#print("original class labels:")	
	#print(labels)
	#print("predicted class labels:")	
	#print(pred.class)
	#print(t(prob))	
	
	tab.classes<-table(pred.class,labels)		
	
				
	if (nrow(tab.classes)!= ncol(tab.classes) | nrow(tab.classes)== 1 )
	{
		### 
		print("original")
		print(tab.classes)
		tabt<-.extend.to.quad.matrix (tab.classes, labels.universe=labels.universe)
		#print(tabt)

		tab2=matrix(0,ncol=2,nrow=2)		
		colnames(tab2)=c(-1,1)
		rownames(tab2)=c(-1,1)
		#print(tab2)
		tab2["-1","-1"]= tabt["-1","-1"]
		tab2["-1","1"]= tabt["-1","1"]
		tab2["1","-1"]= tabt["1","-1"]
		tab2["1","1"]= tabt["1","1"]
		tab.classes=tab2
	}
	
	print(tab.classes)
	
	
	# sensitivity = TP/ all P = TP /(TP + FN)
	sensitivity.classes<- tab.classes[2,2]/sum(tab.classes[,2])
	# specificity = TN/ all N = TN /(TN + FP)
	specificity.classes <- tab.classes[1,1]/sum(tab.classes[,1])
	# secondary diagonal
  
	sensitivity=sensitivity.classes
	specificity=specificity.classes 
  pred <- prediction(prob, labels)
  auc=unlist(performance(pred, "auc")@y.values) 
  res= list(se=sensitivity,sp=specificity,auc=auc)
  return(res)

}



geneRank <- function(W,ex,d, max.degree=Inf)
{


  W=Matrix(W)
  
  ex = abs(ex)

  ## normalize expression values
  norm_ex = ex/max(ex)

  ## try sparse Matrices in R => later
  ##w = sparse(W)
  dimW = dim(W)[1]
  if(dim(W)[2]!=dimW) stop("W must be a square matrix.")
  
  ## get the in-degree for every node
  ## from KEGG we get a directed graph
  ## thus, the column sums correspond to
  ## to the in-degree of a particular gene

  degrees = pmin(max.degree, pmax(1,colSums(W), na.rm=T))

  ## A = Identity Matrix with dimensions
  ## same as the adjacency matrix

  A=Matrix(0, nrow = dimW, ncol = dimW)
  diag(A) = 1
  
  ## produce a matrix with the degrees on
  ## the diagonal

  D1=Matrix(0, nrow = dimW, ncol = dimW)
  diag(D1) = 1.0/degrees

  ## divide the in-degrees of the gene
  ## by the overall in-degree(colSum) of the gene
  ## => kind of normalizing the in-degrees

  A = A - d*(Matrix(t(W)) %*% D1)

  ## here, we give 1-d 'for free'
  b  = (1-d) * norm_ex

  ## we want to solve:
  ## (I - d W^t D^-1)r = (1-d)ex which is the Jacobi of the PageRank
  ## where A = (I - d W^t D^-1)
  ## and   b = (1-d)ex
  ## therefore:  Ar = b
  r = as.numeric(solve(A,b))
  return(r)
  
}
