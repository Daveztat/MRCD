###############
# Murder rate data
# Application that illustrates the MRCD subset selection and the MRCD covariance based regression analysis
# Data is taken from Khan et al. (2007), as published on http://users.ugent.be/~svaelst/software/RLARS.html 
###############

rm(list = ls())  # cleanup

#setwd!

sapply(list("MRCD.R","optimalh.R", "MRCDreg.R"),source)


# 1. Read data
#--------------

murderdata = read.table("demoUSA.txt",header=TRUE)
statenames = murderdata[,1] # first column has the names of the states
murderdata = murderdata[,-1]
rownames(murderdata) = statenames

y = matrix(murderdata[,"M"],ncol=1) # murder rate is the dependent variable
xnames = setdiff(colnames(murderdata),"M")
X = data.matrix( murderdata[,xnames])

# 2. Apply MRCD for regressions
#------------------------------

mX = cbind(y,X) # MRD regression is done on the joint data y+X, as explained in Rousseeuw et al. (2004)

# Select optimal h

start.time <- Sys.time()
seth=seq(ceiling(dim(mX)[1]/2),dim(mX)[1],1)
hresult= optimalh(mX=t(mX), seth=seth)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

plot(seth,hresult$obj, xlab="subset size h", ylab="MRCD objective function",xaxt="n")
axis(1, at = seth, las=2)

msecov=c()
for (i in 1:(length(seth)-1)){
  msecov[i]=mean((hresult$ccc[[i+1]]-hresult$ccc[[i]])^2)
}
plot(seth[-1],msecov, xlab="subset size h", ylab="Frobenius distance between K(h) and K(h-1)", xaxt="n")  
axis(1, at = seth[-1], las=2)


# MRCD regression
mrcd = MRCDreg(y=y,X=X,alpha=44/50)
estbeta = mrcd$coef
# OLS regression
OLSbeta = lm(y ~X)$coef[-1]

select = c(23)
#postscript("betas.eps",width=8,height=5,horizontal=FALSE)
par(mfrow=c(1,1))
par(mar=c(4,4,2,2))
plot(OLSbeta,estbeta,xlab="OLS coefficients",ylab="MRCD coefficients",ylim=c(-1.65,1.3))
text( OLSbeta[select],estbeta[select],colnames(X)[select],cex=0.8,pos=1)
abline(a=0,b=1,lty=3)
abline(h=0,lty=3)
abline(v=0,lty=3)
#dev.off()

data.frame( 1:25, names(OLSbeta) , OLSbeta,estbeta)[select,]

subset = mrcd$index
outliers = setdiff(1:length(y),subset)

#postscript("phones.eps",width=8,height=5,horizontal=FALSE)
par(mfrow=c(1,1))
par(mar=c(4,6,2,2))
plot(X[,"PH"],y ,col="white", xlab="Number of telephones per 100 residents" , ylab="Number of murders \n per 100,000 residents")
text( x=42, y=190 ,expression( beta[OLS]==-0.48~ vs  ~beta[MRCD]==-0.76 ))
lines(X[subset,"PH"],y[subset] ,type="p",col="gray",pch=18)
lines(X[outliers,"PH"],y[outliers] ,type="p", col="red",pch=17)
text( X[outliers,"PH"],y[outliers] ,names(X[outliers,"PH"]),cex=0.8,pos=1)
#dev.off()

  