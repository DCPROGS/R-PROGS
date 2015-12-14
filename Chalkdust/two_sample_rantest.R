# Author: David Colquhoun, University College London
#
#Two_mean_simulation
require(MASS)
#require(graphics)

#START INPUTS
mynsim=100000  #number of re-amplings to run
# Define name of output text file
outfile="rantest2.txt" #name for output file
#define the two samples
#test by using Cushny data
myn=10
mysampA=c(0.7,-1.6,-0.2,-1.2,-0.1,3.4,3.7,0.8,0.0,2)
mysampB=c(1.9,0.8,1.1,0.1,-0.1,4.4,5.5,1.6,4.6,3.4)

#END OF INPUTS
#
#myobsdif= -1.58
meanA=sum(mysampA/myn)
meanB=sum(mysampB/myn)
myobsdif= meanA-meanB
sdA=sd(mysampA)
sdB=sd(mysampB)

#initialisations
#
difmean1=vector(length=mynsim)  #for randist
difmean1[1:mynsim]=0
difmean2=vector(length=mynsim)  #for bootstrap
difmean2[1:mynsim]=0
#
nranlo=0 # number of diffs <= obs diff
nranhi=0 # number of diffs >= minus (obs diff)
nbootlo=0
nboothi=0
#
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
# set random number generator seed so sequence repeats
#set.seed(1941)
#
for (myr in c(1:mynsim))
{
  #Randomisation test
  rsamp1 <-sample(allobs,myn, replace=FALSE)
  mysumA=sum(rsamp1)
  mysumB=totobs-mysumA
  difmean1[myr]=(mysumA-mysumB)/myn
  difmean2[myr]=(mysumA-mysumB)/myn
  if (difmean1[myr] <= myobsdif) nranlo=nranlo+1
  if (difmean1[myr] >= -myobsdif) nranhi=nranhi+1
  

  #repeat for sampling with replacement (bootstrap)
  rsamp2 <-sample(allobs,myn, replace=TRUE)
  mysumA=sum(rsamp2)
  mysumB=totobs-mysumA
  difmean2[myr]=(mysumA-mysumB)/myn
  if (difmean2[myr] <= myobsdif) nbootlo=nbootlo+1
  if (difmean2[myr] >= -myobsdif) nboothi=nboothi+1
  }

# mean and sd
mean1=mean(difmean1)   # for ran dist
sd1=sd(difmean1)
mean2=mean(difmean2)   # for bootstrap
sd2=sd(difmean2)
# percentiles as limits
lo95lim=quantile(difmean1, probs=0.025, names=FALSE)
hi95lim=quantile(difmean1, probs=0.975, names=FALSE)
lo95lim2=quantile(difmean2, probs=0.025, names=FALSE)
hi95lim2=quantile(difmean2, probs=0.975, names=FALSE)

# proportion of values below obs diff
ptoplo=nranlo/mynsim
# proportion of values above minus() obs diff)
ptophi=nranhi/mynsim
#two tail P
p2tail=ptoplo+ptophi
#
#plot histogram of diff between means
X11()   #next plot in new window
xlo=mean1-4*sd1
xhi=mean1+4*sd1
#ensurobs difference (red line) is in scale
if(myobsdif<=xlo) xlo=min(trunc(myobsdif-2*sd1),xlo)
if(myobsdif>=xlo) xhi=max(xhi,trunc(myobsdif+2*sd1))
xlim1=c(xlo,xhi)
# supresses x axis http://www.statmethods.net/advgraphs/axes.html

mydat=hist(difmean1,breaks=40,xlim=xlim1,cex.axis=1.5,xaxt="n",main="Randomsation distribution (simulations)", xlab="difference between means")
axis(side=1,pos=0,at=seq(-3,3,1))  #draw x axis at y = 0
summary(mydat)
# red line for obsdif
# lines for quantiles
abline(v=myobsdif, col="red", lwd=2)
abline(v=c(lo95lim,hi95lim),col=320,lty=2, lwd=2)   #quantiles
#
#
str(mydat)
#
#write results to a file (use 'cat' not 'write' to output list 
#
cat("Randomisation distribution using random samples","\n",file=outfile,append=FALSE)
cat("INPUTS","\n",file=outfile,append=TRUE)
cat("number of resamplings = ",mynsim,"\n",file=outfile,append=TRUE)
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
cat("quantiles for ran dist (0.025, 0.975)", lo95lim,hi95lim, "\n","\n", file=outfile,append=TRUE)
#
cat("Area below observed diff",ptoplo,"\n",file=outfile,append=TRUE)
cat("Area above minus observed diff",ptophi,"\n",file=outfile,append=TRUE)
cat("Two-tailed P value",p2tail,"\n","\n",file=outfile,append=TRUE)
#
cat("Result of t test","\n",file=outfile,append=TRUE)
cat("P value (2 tail)", mypval, "\n",file=outfile,append=TRUE)
cat("confidence interval", myloCI, myhiCI, "\n",file=outfile,append=TRUE)