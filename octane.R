###############
# Octane data:
# Application that illustrates the MRCD subset selection and the MRCD outlier detection analysis
###############

rm(list = ls())  # cleanup

#setwd!

sapply(list("MRCD.R","optimalh.R"),source)

library('rrcov') #for PcaHubert

# 1. Read data
#--------------

data = read.csv( paste( getwd(),"/octane.csv",sep=""),header=F)
octane = data[,2:dim(data)[2]]
dim(octane) # n x p = 39 x 226
n = dim(octane)[1]
p = dim(octane)[2]
obs = 1:n

# 2. Apply ROBPCA (robust PCA analysis)
#--------------------------------------
robpca = PcaHubert(octane, k=2, kmax=p, alpha=0.75, mcd=FALSE)
outl.robpca = obs[robpca@flag==FALSE]; outl.robpca
# Observations flagged as outliers by ROBPCA:
# 25, 26, 36, 37, 38, 39

# Plot the orthogonal distances versus the score distances
pch = rep(20,n); pch[robpca@flag==FALSE] = 17
col = rep('black',n); col[robpca@flag==FALSE] = 'red'
plot(robpca, pch=pch, col=col)


# 3. Apply MRCD for outlier detection
#------------------------------------

mX=t(octane)

# Select optimal h

start.time <- Sys.time()
#consider all possible h-values
seth=seq(ceiling(dim(mX)[2]/2),dim(mX)[2],1)
hresult=optimalh(mX=mX, seth=seth)
end.time <- Sys.time()
#computation time needed
time.taken <- end.time - start.time
time.taken


plot(seth,hresult$obj, xlab="subset size h", ylab="MRCD objective function",xaxt="n")
axis(1, at = seth, las=2)

msecov=c()
for (i in 1:(length(seth)-1)){
  msecov[i]=mean((hresult$ccc[[i+1]]-hresult$ccc[[i]])^2)
}
plot(seth[-1],msecov, xlab="subset size h", ylab="Frobenius distance between K(h) and K(20)", xaxt="n") 
axis(1, at = seth[-1], las=2)

hopt=33

result = mrcd(mX=mX,h=hopt)

dist=mahalanobis(octane, center=result$mu, cov=result$icov, inverted=TRUE)

#postscript("octanedist.eps",width=8,height=5,horizontal=FALSE)
par(mar=c(4,4,2,2))
plot(sqrt(dist), xlab='Index i', ylab='MRCD robust distance', main='', pch=pch, col=col)
#dev.off()

result$rho



 






