##
## define working directory
setwd("put directory here")
## read in data
rawdata <- read.csv('polydata.csv',header=T)

#run eigenvector spatial filtering
library(spfilteR)
library(parallel)

#define matrices for pairwise spatial similarity and phylogenetic similarity
#by definition, the diagonal elements are zero and each row is normalized to 1
#use the matrices derived in Hua et al. 2022
load("Wnb.Rdata") #if two langauges overlap
load("Wphy.Rdata") #similarity based on phylogenetic correlation
load("Wsp.Rdata") #similarity based on spatial distance-decay in correlation

#extract eigenvectors
#using SAR model
EVs.nb <- getEVs(W=Wnb,covars=NULL)
EVs.phy <- getEVs(W=Wphy,covars=NULL)
EVs.sp <- getEVs(W=Wsp,covars=NULL)
Ec.nb <- EVs.nb$moran/max(EVs.nb$moran)>=.25 #the threshold 0.25 is suggested by Griffith (2003)
Ec.phy <- EVs.phy$moran/max(EVs.phy$moran)>=.25
Ec.sp <- EVs.sp$moran/max(EVs.sp$moran)>=.25
evecs <- cbind(EVs.nb$vectors[,Ec.nb],EVs.phy$vectors[,Ec.phy],EVs.sp$vectors[,Ec.sp])

#following Koplenig 2024, using intercept-only model
#eigenvectors are selected by the objective: eigenvectors are added to the model until no more 
#	eigenvectors can further reduce the residual autocorrelation. 

########################
#polysynthesis analysis
########################
#defind variables
variables <- c("Polysynthetic", "Extended", "L1.pop", "Area", "Bordering.languages",
               "Bordering.languages.km", "Vehicularity", "Roughness", "Altitude", 
               "Small.family"
               )

data <- rawdata[,variables]

#log transform highly right skewed predictors
data$L1.pop <- log(data$L1.pop+0.5)
data$Area <- log(data$Area)
data$Bordering.languages <- sqrt(data$Bordering.languages)
data$Bordering.languages.km <- sqrt(data$Bordering.languages.km)

#scaling continuous data to help interpreting results
data$L1.pop <- scale(data$L1.pop)
data$Area <- scale(data$Area)
data$Bordering.languages <- scale(data$Bordering.languages)
data$Bordering.languages.km <- scale(data$Bordering.languages.km)
data$Roughness <- scale(data$Roughness)
data$Altitude <- scale(data$Altitude)

#define response variable
y <- data$Polysynthetic
#y <- data$Extended

#define similarity matrix to calculate Moran's I
W <- list(Wnb,Wphy,Wsp)

#define intercept
n <- length(y)
x <- as.matrix(rep(1, n))

#define objective function
objfunc <-function(y, xe, W) {
    o <- glm(y~xe-1, family=binomial(link=probit), maxit = 100)
    resid <- (y - o$fitted.values) / sqrt(o$fitted.values * (1 - o$fitted.values)) 
    test <- 0
    for (i in 1:length(W)) {
     	test <- test + spfilteR:::MI.resid(resid = resid, x = xe, W = W[[i]])$zI^2
    }
    return(test)
  }

#calculate initial I^2 values
init <- glm(y~x-1, family=binomial(link=probit), maxit = 100)
Rsquared_init <- 1-(init$deviance/init$null.deviance)
resid_init <- (y - init$fitted.values) / sqrt(init$fitted.values * (1 - init$fitted.values)) #pearson residual
zMI_init <- 0
for (i in 1:length(W)) {
     zMI_init <- zMI_init + spfilteR:::MI.resid(resid = resid_init, x = x, W = W[[i]])$zI^2
}
pMI_init <- 1-pchisq(zMI_init,df=3)

#search algorithm
oldpMI <- pMI_init
oldzMI <- zMI_init
sel <- as.logical(rep(1,dim(evecs)[2]))
sel_id <- NULL
selset <- which(sel)

for (i in which(sel)) {

	ref <- Inf
	sid <- NULL
	
	f <- function (j,y,x,W,evecs,sel_id) {
        o <- try(objfunc(y = y, xe = cbind(x, evecs[, sel_id], evecs[, j]), W = W),silent=T)
        if (class(o)=="try-error" || is.na(o)) {
        	o <- Inf
        }
        o
	}
	
	test <- unlist(mclapply(selset, f, y, x, W, evecs, sel_id, mc.cores=10))
		
	sid <- which(test==min(test))
	ref <- test[sid]
	
     # stopping rules
     pMI <- 1-pchisq(ref,df=3)
     if (pMI > oldpMI || ref < oldzMI) {
     	 oldzMI <- ref
         oldpMI <- pMI
         sel_id <- c(sel_id, selset[sid])
     } else {
     	break #adding more eigenvectors does not further reduce Moran's I
     }
      
     # remove selected eigenvectors from candidate set
     selset <- selset[!(selset %in% sel_id)]
}
sel_id
 # end selection

# filtered regression
xev <- cbind(x,evecs[, sel_id])
out <- glm(y~xev-1, family=binomial(link=probit), maxit = 100)
resid_out <- (y - out$fitted.values) / sqrt(out$fitted.values * (1 - out$fitted.values)) 
zMI_out <- 0
for (i in 1:length(W)) {
   zMI_out <- zMI_out + spfilteR:::MI.resid(resid = resid_out, x = xev, W = W[[i]])$zI^2
}
pMI_out <- 1-pchisq(zMI_out,df=3)
   
#run all possible combinations of predictors
all_vars <- variables[-c(1,2)]
nx <- length(all_vars)
combs <- lapply(1:nx, function(x) combn(all_vars, x))

model_combos <- lapply(combs, function(x) as.list(as.data.frame(x)))
model_combos <- Reduce(c, model_combos)

f <- function (j,y,xev,data) {
		x <- as.matrix(cbind(xev,data[,as.character(j)]))
        glm(y~x-1,family=binomial(link=probit), maxit = 100)
	}
res <- mclapply(model_combos, f, y, xev, data, mc.cores=10)

#calculate model probabiltiy
logp <- sapply(res,function (i) length(i$coefficients)-i$aic/2)
pp <- exp(logp-max(logp)) 
pp <- pp/sum(pp)

#extract estimated regression coefficients
beta <- lapply(1:length(res),function (j) {
	coef <- res[[j]]$coefficients[-c(1:(length(sel_id)+1))]
	names(coef) <- model_combos[[j]]
	coef
	})

#calculate bma and pip
bma <- pip <- rep(NA,nx-1)
names(bma) <- names(pip) <- all_vars
for (var in all_vars) {
	idx <- sapply(1:length(res),function (i) is.element(var,model_combos[[i]]))
	beta_var <- sapply(which(idx),function (i) beta[[i]][names(beta[[i]])==var])
	bma[var] <- sum(beta_var*pp[idx])
	pip[var] <- sum(pp[idx])
}

#write result table
result <- matrix(NA,length(res),length(all_vars)+2)
colnames(result) <- c("AIC","df",all_vars)
for (i in 1:length(res)) {
	nam <- as.character(model_combos[[i]])
	result[i,nam] <- res[[i]]$coefficients[-c(1:dim(xev)[2])]
	result[i,"AIC"] <- res[[i]]$aic
	result[i,"df"] <- length(res[[i]]$coefficients)
}

########################
#endangerment analysis
########################
library(MASS)
data <- read.csv('data.csv',header=T)
variable <- read.csv('variable.csv',header=T)

#define response variable
y <- data$EGIDS_tr
levels(y)[1:6] <- "6a"
J <- c(1:length(levels(y)))

#add polysynthetic variables
data <- cbind(data,rawdata[,c("Polysynthetic","Extended")])
#add interaction terms
Geninter <- function (var,region,nam) {
	region2 <- unique(region)
	region_nam <- gsub("\\s*\\([^\\)]+\\)","",as.character(region2))
	region_nam <- gsub(" ","_",region_nam,fixed=T)
	region_nam <- paste0(nam,"_",region_nam)
	out <- matrix(0,length(var),length(region_nam))
	colnames(out) <- region_nam
	idx <- NULL
	for (i in 1:length(region_nam)) {
		tmp <- var[region==region2[i]]
		if (sd(tmp,na.rm=T)>0) {idx <- c(idx,i)}
		out[region==region2[i],region_nam[i]] <- tmp
	}
	out[,idx,drop=F]	
}
variable2 <- variable[-((n.var-11):n.var),]
for (i in 1:dim(variable2)[1]) {
  var <- data[,which(names(data)==variable2[i,1])]
  data <- cbind(data,Geninter(var,data$lang_subregion_lang,variable[i,1]))
  add <- colnames(Geninter(var,data$lang_subregion_lang,variable2[i,1]))
  if (length(add)>1) {
  variable <- rbind(variable,data.frame(variable=add,group=variable2[i,2]))
  }
  }

#define objective function
objfunc <-function(y, xe, W) {
    o <- polr(y~xe, method="probit")
    resid <- o$fitted.values%*%J / sqrt(o$fitted.values%*%(J^2)-(o$fitted.values%*%J)^2) 
    test <- 0
    for (i in 1:length(W)) {
     	test <- test + spfilteR:::MI.resid(resid = resid, x = xe, W = W[[i]])$zI^2
    }
    return(test)
  }

#calculate initial I^2 values
init <- polr(y~1, method="probit")
resid_init <- init$fitted.values%*%J / sqrt(init$fitted.values%*%(J^2)-(init$fitted.values%*%J)^2) 
zMI_init <- 0
for (i in 1:length(W)) {
     zMI_init <- zMI_init + spfilteR:::MI.resid(resid = resid_init, x = x, W = W[[i]])$zI^2
}
pMI_init <- 1-pchisq(zMI_init,df=3)

#search algorithm
oldpMI <- pMI_init
oldzMI <- zMI_init
sel <- as.logical(rep(1,dim(evecs)[2]))
sel_id <- NULL
selset <- which(sel)

for (i in which(sel)) {

	ref <- Inf
	sid <- NULL
	
	f <- function (j,y,x,W,evecs,sel_id) {
        o <- try(objfunc(y = y, xe = cbind(x, evecs[, sel_id], evecs[, j]), W = W),silent=T)
        if (class(o)=="try-error" || is.na(o)) {
        	o <- Inf
        }
        o
	}
	
	test <- unlist(mclapply(selset, f, y, x, W, evecs, sel_id, mc.cores=10))
		
	sid <- which(test==min(test))
	ref <- test[sid]

     # stopping rules
     pMI <- 1-pchisq(ref,df=3)
     if (pMI > oldpMI || ref < oldzMI) {
     	 oldzMI <- ref
         oldpMI <- pMI
         sel_id <- c(sel_id, sid)
     } else {
     	break #adding more eigenvectors does not further reduce Moran's I
     }
      
     # remove selected eigenvectors from candidate set
     selset <- selset[!(selset %in% sel_id)]
}
sel_id
 # end selection

# filtered regression
xev <- cbind(x,evecs[, sel_id])
out <- polr(y~xev, method="probit")
resid_out <- out$fitted.values%*%J / sqrt(out$fitted.values%*%(J^2)-(out$fitted.values%*%J)^2) 
zMI_out <- 0
for (i in 1:length(W)) {
   zMI_out <- zMI_out + spfilteR:::MI.resid(resid = resid_out, x = xev, W = W[[i]])$zI^2
}
pMI_out <- 1-pchisq(zMI_out,df=3)

# model selection
#start stepwise selection procedure
X <- data[,as.character(variable[,1])]
n.var <- dim(X)[2]
loglik <- numeric(n.var)
names(loglik) <- colnames(X)[1:n.var]
tmp <- which(!is.na(rowSums(X)))
X.tmp <- X[tmp,]
y.tmp <- y[tmp]
xev.tmp <- xev[tmp,]
for (i in 1:n.var) {
	x <- cbind(xev.tmp,X.tmp[,i])
    out <- try(logLik(polr(y.tmp~x, method="probit")),silent=T)
    if (inherits(out,"try-error")) {
       			 	loglik[i] <- -Inf
    } else {
    	loglik[i] <- out
    }
}
   
library(foreach)
library(doParallel)
library(MASS)
registerDoParallel(15)
model <- vector("list",50)
var.group <- sapply(unique(variable[,2]),function (i) variable[variable[,2]==i,1])

for (aa in 1:700) {
set.seed(aa)
basemodel <- names(which(loglik==max(loglik)))
baseloglik <- max(loglik)
var.group1 <- var.group
basemodel1 <- NULL
while (!setequal(basemodel,basemodel1)) {
    basemodel1 <- basemodel
    ord <- sample(1:length(var.group),length(var.group),replace=F)
    for (j in ord) {
        vars <- as.character(var.group1[[j]])
        vars <- vars[!is.element(vars,basemodel1)]
        if (length(vars)>0) {
    		logLik <- foreach (i = 1:length(vars), .packages='MASS') %dopar% {
       			 altermodel <- c(basemodel,vars[i])
       			 x <- as.matrix(cbind(xev.tmp,X.tmp[,altermodel]))
       			 out <- try(logLik(polr(y.tmp~x, method="probit",start=c(out$coefficients,numeric(length(altermodel)),out$zeta),control=list(maxit=1000))),silent=T)
       			 if (inherits(out,"try-error")) {
       			 	out <- -Inf
       			 }
       			 out
    		}
    		logLik <- as.numeric(logLik)
    		names(logLik) <- vars
    		likratio.test <- 1-pchisq(2*(logLik-baseloglik),df=1)
    		tmp1 <- vars[which(likratio.test<=(0.05/length(likratio.test)))]
    		if (length(tmp1)>0) {
        		basemodel <- c(basemodel,tmp1)
        		x <- as.matrix(cbind(xev.tmp,X.tmp[,basemodel]))
        		baseloglik <- try(logLik(polr(y.tmp~x, method="probit",start=c(out$coefficients,numeric(length(basemodel)),out$zeta))),silent=T)
        		var <- tmp1
        		var.group1[[j]] <- setdiff(vars,var)
        	}
        }
    	vars <- basemodel
    	if (length(vars)>0) {
    		logLik <- foreach (i = 1:length(vars), .packages='MASS') %dopar% {
        		altermodel <- basemodel[-which(basemodel==vars[i])]
        		x <- as.matrix(cbind(xev.tmp,X.tmp[,altermodel]))
       			 out <- try(logLik(polr(y.tmp~x, method="probit"),start=c(out$coefficients,numeric(length(altermodel)),out$zeta),control=list(maxit=1000)),silent=T)
       			 if (inherits(out,"try-error")) {
       			 	out <- -Inf
       			 }
       			 out
    		}
    		logLik <- as.numeric(logLik)
    		names(logLik) <- vars
    		likratio.test <- 1-pchisq(2*(-logLik+baseloglik),df=1)
    		tmp1 <- vars[which(likratio.test>(0.05/length(likratio.test)))]
    		if (length(tmp1)>0) {
        		basemodel <- setdiff(basemodel,tmp1)
        		x <- as.matrix(cbind(xev.tmp,X.tmp[,basemodel]))
        		baseloglik <- try(logLik(polr(y.tmp~x, method="probit",start=c(out$coefficients,numeric(length(basemodel)),out$zeta))),silent=T)
        		var <- tmp1
        		jj <- sapply(var,function(y) which(sapply(var.group, function (x) is.element(y,x))))
        		for (ii in var) {
        		if (!is.element(ii,var.group1[[jj[ii]]])) {
        			var.group1[[jj[[ii]]]] <- c(var.group1[[jj[ii]]],ii)
        		}
        		}
    		}
    	}
    }
    print(c(length(basemodel),length(basemodel1)))
}
model[[aa]] <- basemodel
}

#remove duplicated model
for (i in 1:699) {
	for (j in (i+1):700) {
		idx <- sum(is.element(model[[i]],model[[j]]))==length(model[[i]])
		if (idx) {
			model[[i]] <- NA
		}
	}
}
model <- model[!is.na(model)]

res <- foreach (i = 1:length(model), .packages='MASS') %dopar% {
       			 j <- model[[i]]
	x <- as.matrix(cbind(xev.tmp,X.tmp[,j]))
    polr(y.tmp~x, method="probit",start=c(out$coefficients,numeric(length(j)),out$zeta),control=list(maxit=1000))
    }
    
#including surrounding models
logp <- sapply(res,function (i) -i$deviance/2)
pp <- exp(logp-max(logp)) 
pp <- pp/sum(pp) 
names(pp) <- c(1:700)
# 343          242           55           94 
#5.852797e-03 2.241671e-02 3.129696e-02 9.397084e-01 
keep <- 94

for (i in keep) {
	vars <- all_vars[!is.element(all_vars,model2[[i]])]
	#add a variable
	model <- c(model,lapply(vars,function (x) c(model2[[i]],x)))
	#remove a variable
	tmp <- lapply(model2[[i]],function (x) model2[[i]][-which(model2[[i]]==x)])
	model <- c(model,tmp)
	#add a variable to each removed model
	for (y in 1:length(tmp)) {
		model <- c(model,lapply(vars,function (x) c(tmp[[y]],x)))
	}
}

res <- foreach (i = 1:length(model), .packages='MASS') %dopar% {
       			 j <- model[[i]]
	x <- as.matrix(cbind(xev.tmp,X.tmp[,j]))
    polr(y.tmp~x, method="probit",start=c(out$coefficients,numeric(length(j)),out$zeta),control=list(maxit=1000))
    }

#calculate model probabiltiy
logp <- sapply(res,function (i) -i$deviance/2)
pp <- exp(logp-max(logp))
pp <- pp/sum(pp)

#extract coefficients
beta <- lapply(1:length(res),function (j) {
	coef <- res[[j]]$coefficients[-c(1:dim(xev)[2])]
	names(coef) <- model[[j]]
	coef
	})

#calculate bma and pip
all_vars <- as.character(variable[,1])
bma <- pip <- rep(NA,length(all_vars))
names(bma) <- names(pip) <- all_vars
for (var in all_vars) {
	idx <- sapply(1:length(res),function (i) is.element(var,model[[i]]))
	if (sum(idx)>0) {
		beta_var <- sapply(which(idx),function (i) beta[[i]][names(beta[[i]])==var])
		bma[var] <- sum(beta_var*pp[idx])
		pip[var] <- sum(pp[idx])
	}
}
bma <- bma[!is.na(bma)]
pip <- pip[!is.na(pip)]