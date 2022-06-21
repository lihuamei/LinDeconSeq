solveRLMInternal<-function(S, B){
	solution<- coef(MASS::rlm(B ~ S, maxit = 100000))
	names(solution)<-colnames(S)
	return(solution)
}

solveDampenedRLM <- function(S, B, intercept){
	if (intercept) S <- cbind(1, S)
	solution<-solveOLSInternal(S, B)
	iterations <- 0; changes <- c()
	j <- findDampeningConstant(S,B,solution)
	change<-1
	while(change>.01 & iterations<1000){
		newsolution <- solveDampenedWLSj(S, B, solution, j)
		solutionAverage <- rowMeans(cbind(newsolution, matrix(solution,nrow = length(solution), ncol = 4)))
		change<-norm(as.matrix(solutionAverage - solution))
		solution<-solutionAverage
		iterations<-iterations + 1
		changes <- c(changes, change)
	}
	return(solution/sum(solution))
}

#solve WLS given a dampening constant
solveDampenedWLSj <- function(S,B,goldStandard,j){
	multiplier<-1*2^(j-1)
	sol<-goldStandard
	ws<-as.vector((1/(S%*%sol))^2)
	wsScaled<-ws/min(ws)
	wsDampened<-wsScaled
	wsDampened[which(wsScaled>multiplier)]<-multiplier
	W<-diag(wsDampened)
	D<-t(S)%*%W%*%S
	d<- t(S)%*%W%*%B
	A<-cbind(diag(dim(S)[2]))
	bzero<-c(rep(0,dim(S)[2]))
	sc <- norm(D,"2")
	BB <- B/sc; SS <- S/sc
	solution<- coef(MASS::rlm(BB ~ SS - 1, maxit = 100000))
	names(solution)<-colnames(S)
	return(solution)
}

#find a dampening constant for the weights using cross-validation
findDampeningConstant<-function(S,B,goldStandard){
	solutionsSd<-NULL
	#goldStandard is used to define the weights
	sol<-goldStandard
	ws<-as.vector((1/(S%*%sol))^2)
	wsScaled<-ws/min(ws)
	wsScaledMinusInf<-wsScaled
	#ignore infinite weights
	if(max(wsScaled)=="Inf"){
		wsScaledMinusInf<-wsScaled[-which(wsScaled=="Inf")]
	}
	#try multiple values of the dampening constant (multiplier)
	#for each, calculate the variance of the dampened weighted solution for a subset of genes
	for (j in 1:ceiling(log2(max(wsScaledMinusInf)))){
		multiplier<-1*2^(j-1)
		wsDampened<-wsScaled
		wsDampened[which(wsScaled>multiplier)]<-multiplier
		solutions<-NULL
		seeds<-c(1:100)
		for (i in 1:100){
			set.seed(seeds[i]) #make nondeterministic
			subset<-sample(length(ws),size=length(ws)*0.5) #randomly select half of gene set
			#solve dampened weighted least squares for subset
			fit = lm (B[subset] ~ -1+S[subset,],weights=wsDampened[subset])
			sol<-fit$coef*sum(goldStandard)/sum(fit$coef)
			solutions<-cbind(solutions,sol)
		}
		solutionsSd<-cbind(solutionsSd,apply(solutions,1,sd))
	}
	#choose dampening constant that results in least cross-validation variance
	j<-which.min(colMeans(solutionsSd^2))
	return(j)
}
