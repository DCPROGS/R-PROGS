# Author: David Colquhoun, University College London
#
#Two_sample-randomisation test
# exact version, via all possible samples
require(MASS)
#require(graphics)
#
#START INPUTS
#define the two samples
#test by using Cushny data
myn=10 #same for A and B
mynt=20 #total number ofobservations
mysampA=c(0.7,-1.6,-0.2,-1.2,-0.1,3.4,3.7,0.8,0.0,2)
mysampB=c(1.9,0.8,1.1,0.1,-0.1,4.4,5.5,1.6,4.6,3.4)
# Define name of output text file
outfile="rantest-exact2.txt" #name for output file
#END OF INPUTS
#
#myobsdif= -1.58
meanA=sum(mysampA/myn)
meanB=sum(mysampB/myn)
myobsdif= meanA-meanB
sdA=sd(mysampA)
sdB=sd(mysampB)
#
#initialisations
#
#Number of possible random samples
ncom=choose(mynt,myn)
mycmat=combn(mynt,myn)
# mycmat is a matrix with columns= randomssample of 10 integers,# the indices of the sample of 10 obs out of 20.  
# The ncom columns ofmycmat list all possible choices

difmean1=vector(length=ncom)  #for randist
difmean1[1:ncom]=0
#
nranlo=0 # number of diffs <= obs diff
nranhi=0 # number of diffs >= minus (obs diff)

# Do t test
myresult=t.test(mysampB,mysampA,alternative="two.sided", paired = FALSE, var.equal = TRUE,
                conf.level = 0.95)
mydiff = myresult$estimate[1]-myresult$estimate[2]
myp=myresult$p.value
mypval=myp 
myloCI=myresult$conf.int[1]
myhiCI=myresult$conf.int[2]
#
#put all obs into one array
allobs <- c(mysampA, mysampB)
totobs=sum(allobs)
#
for (myr in c(1:ncom))  #go through all possible samples
{
  rsamp1=allobs[mycmat[1:myn,myr]]
  mysumA=sum(rsamp1)
  mysumB=totobs-mysumA
  difmean1[myr]=(mysumA-mysumB)/myn
  if (difmean1[myr] <= myobsdif) nranlo=nranlo+1
  if (difmean1[myr] >= -myobsdif) nranhi=nranhi+1
}  
# mean and sd
mean1=mean(difmean1)   # for ran dist
sd1=sd(difmean1)
#
# percentiles as limits
lo95lim=quantile(difmean1, probs=0.025, names=FALSE)
hi95lim=quantile(difmean1, probs=0.975, names=FALSE)

# proportion of values below obs diff
ptoplo=nranlo/ncom
# proportion of values above minus() obs diff)
ptophi=nranhi/ncom
#two tail P
p2tail=ptoplo+ptophi
#
#plot ogram of diff between means
X11()   #next plot in new window
xlo=mean1-4*sd1
xhi=mean1+4*sd1
#ensurobs difference (red line) is in scale
if(myobsdif<=xlo) xlo=min(trunc(myobsdif-2*sd1),xlo)
if(myobsdif>=xlo) xhi=max(xhi,trunc(myobsdif+2*sd1))
xlim1=c(xlo,xhi)
#xaxt="n) supresses x axis http://www.statmethods.net/advgraphs/axes.html
mydat=hist(difmean1,breaks=40,xlim=xlim1,xaxt="n", cex.axis=1.5,main="Randomisation distribution (all possible samples)", xlab="difference between means")
axis(side=1,pos=0,at=seq(-3,3,1))  #draw x axis at y = 0
summary(mydat)
# red line for obsdif
# lines for quantiles
abline(v=myobsdif, col="red", lwd=2)
abline(v=c(lo95lim,hi95lim),col=320,lty=2, lwd=2)   #quantiles
#
#
#write results to a file (use 'cat' not 'write' to output list 
cat("Randomisation test: exact calculation all possible samples","\n","\n",file=outfile,append=FALSE)
cat("INPUTS: exact calculation: all possible samples","\n",file=outfile,append=TRUE)
cat("Total number of combinations = ",ncom,"\n",file=outfile,append=TRUE)
cat("number obs per sample = ",myn,"\n",file=outfile,append=TRUE)
cat("sample A", mysampA, "\n", file=outfile,append=TRUE)
cat("sample B", mysampB,"\n", file=outfile,append=TRUE)

cat("\n","OUTPUTS","\n",file=outfile,append=TRUE)
cat("mean for sample A= ",meanA,"\n",file=outfile,append=TRUE)
cat("mean for sample B = ",meanB,"\n",file=outfile,append=TRUE)
cat("Observed difference between means (A-B)",myobsdif, "\n",file=outfile,append=TRUE)
cat("SD for sample A) = ",sdA,"\n",file=outfile,append=TRUE)
cat("SD for sample B) = ",sdB,"\n",file=outfile,append=TRUE)
cat("mean and SD for randomisation dist = ",mean1,sd1,"\n",file=outfile,append=TRUE)
cat("quantiles for ran dist (0.025, 0.975)", lo95lim,hi95lim, "\n", file=outfile,append=TRUE)
cat("Area equal to orless than observed diff",ptoplo,"\n",file=outfile,append=TRUE)
cat("Area equal to or greater than minus observed diff",ptophi,"\n",file=outfile,append=TRUE)
cat("Two-tailed P value",p2tail,"\n","\n",file=outfile,append=TRUE)
#
cat("Result of t test","\n",file=outfile,append=TRUE)
cat("P value (2 tail)", mypval, "\n",file=outfile,append=TRUE)
cat("confidence interval", myloCI, myhiCI, "\n",file=outfile,append=TRUE)