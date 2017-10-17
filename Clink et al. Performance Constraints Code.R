#### R code for Clink et al. 
#### Part 1 Bandwidth calculation
#### Part 2 Signal-to-noise ratio calculation
#### Part 3 Model selection code


# Set working directory to location of .wav and .csv files
setwd("/Users/denajaneclink/Desktop/Performance Constraints/Performance_constraints_wav_files")


# Load required libraries
library(tuneR)
library(seewave)
library(ggplot2)
library(zoo)
library(plotrix)
library(plyr)
library(stringr)
library(lme4)
library(usdm)
library(lme4)
library(bbmle)
library(MuMIn)
library(gridExtra)
library(lattice)


#### Part 1 Bandwidth calculation  #### 

#Set input directory to read in .wav files
input.dir <-"/Users/denajaneclink/Desktop/Performance Constraints/Performance_constraints_wav_files/TrillNotes"
L = list.files(input.dir, pattern="*.wav", full.names=FALSE)
L
filehandles <- str_split_fixed(L, pattern=".wav", n=2)[,1]

# Create an empty list that will store bandwidth values
bandwidth.list = list()


# Function to calculate bandwidth and create bandwidth figure
bandwidth.function <- function (spec, f = NULL, level = -3, plot = TRUE, colval = "black",
                                cexval = 1, fontval = 1, flab = "Frequency (Hz)", alab = "Relative amplitude (dB)",
                                type = "l", ...)
{
  if (is.null(f)) {
    if (is.vector(spec))
      stop("'f' is missing")
    else if (is.matrix(spec))
      f <- spec[nrow(spec), 1] * 2000
  }
  if (is.matrix(spec))
    spec <- spec[, 2]
  range <- c(f/2000/length(spec), f/2000)
  if (max(spec) == 1)
    stop("data must be in dB")
  if (which.max(spec) == 1)
    stop("maximal peak cannot be the first value of the spectrum")
  n0 <- length(spec)
  spec1 <- approx(spec, n = 102400)
  spec1 <- as.matrix(spec1$y)
  n1 <- nrow(spec1)
  level2 <- round(max(spec1[, 1]), 1) + level
  f0 <- which.max(spec1[, 1])
  specA <- as.matrix(spec1[1:f0, 1])
  nA <- nrow(specA)
  specB <- as.matrix(spec1[f0:nrow(spec1), 1])
  f1 <- which(round(specA, 1) == level2)
  f1khz <- ((f1[length(f1)]/n1) * (range[2] - range[1])) +
    range[1]
  f2 <- which(round(specB, 1) == level2) + (nA - 1)
  f2khz <- ((f2[1]/n1) * (range[2] - range[1])) + range[1]
  Q <- f2khz-f1khz
  results <- list(Q=Q, bdw=f2khz-f1khz)
  if (plot) {
    x <- seq(range[1], range[2], length.out = n0)
    plot(x = x, y = spec, xlab = flab, ylab = alab, type = type,
         ...)
    arrows(f1khz, level2, f2khz, level2, length = 0.1, col = colval,
           code = 3, angle = 15)
    text(paste("Bandwidth at -15 dB (Hz) =", as.character(round(Q, 2))), x = f2khz,
         y = level2, pos = 4, col = colval, cex = cexval,
         font = fontval)
    invisible(results)
  }
  return(results)
}


# Loop to calculate bandwidth
for (j in 1:length(filehandles)) { 
      filehandle <-  L[j]
      filename <- paste(input.dir, "/", filehandle, sep="")
      print(paste("processing", filehandle))
    
    # Read in wav file and downsample (which also filters frequencies of 4 kHz and above)
      wav <- readWave(filename)
    
      wav <-downsample(wav, samp.rate=8000)
      N <- length(wav@left) 
    
    # Pad with zeros to achieve length 8000.
      waveform <- c(wav@left, rep(0, 8000 - N))
    
    # Re-insert into wave object. 
      wav@left <- waveform   
      
    #Calculate the power spectral density which provdes a 1-hertz resolution
      psd <- spec(wav,f=8000, plot=T, dB="max0",wl=8000)
  
    # Calculate the moving average to smooth out the psd
      frequency <- rollapply((seq(from=1,to=nrow(psd),by=1)),width=100, by=50,median) 
      amplitude <- rollapply(psd[,2],width=100, by=50, mean)
      
    # Normalize the amplitude
      amplitude <- amplitude-max(amplitude)
      
    # Combine by column into new object that can be used to calculated bandwidth
      psd.comb<-cbind(frequency,amplitude)
      bandwidth <- bandwidth.function(psd.comb, level=-15, plot=T)$bdw
      
    # Print each value once calculated  
      print(bandwidth)
      
    # Add to the list object created earlier
      bandwidth.list[[j]]<- bandwidth
}

# Check bandwidth list to make sure values are reasonable
bandwidth.list



#### Part 2 Signal-to-noise ratio calculation #### 

# Read in data file
snr.data <- read.csv("data.for.snr.csv",header=T)

# Create an empty list to store SNR values
snr.list = list()

# Create index of calls for the loop

great.call.index <- unique(snr.data$great.call)

for (j in 1:length(great.call.index)) { 
  # Read each .wav file
    great.call.in <- great.call.index[j]
    print(paste("processing file", great.call.in)) 
    newdata <- subset(snr.data, great.call==great.call.in)
    great.call.in <- paste(great.call.in,".wav",sep="")
    wav <- readWave(great.call.in)

  # Define the trill bin start and stop times
    trill.start.time <- min(newdata$Begin.Time..s.)
    trill.last.bin <- max(newdata$End.Time..s.)
    n.bins <- as.integer(trill.last.bin-trill.start.time,length=0)
    bin.seq <- seq(from=0, to=n.bins,by=1 )
    bin.seq <- trill.start.time + bin.seq
    
  # Create a list of all the 1-s bins for SNR analysis
    bin.seq.length <- length(bin.seq)-1
    subsamps <- lapply(1:bin.seq.length, function(i) extractWave(wav, from=as.numeric(bin.seq[i]), to=as.numeric(bin.seq[i+1]), xunit = c("time"),plot=F,output="Wave"))
    
  # Calculate SNR
    for (k in 1:length(subsamps)) { 
  
  # Read in .wav file 
      wav <- subsamps[[k]]
  
  # Filter the .wav file so focus on the gibbon frequency range 
      w.dn.filt <- fir(wav, from=450, to=1800,output = "Wave") 
      w.dn.filt <- normalize(w.dn.filt, unit="16")
  
  # Now have a filtered .wav file
  
  # Create a spectrogram matrix
      w.dn.filt.spectro <- spectro(w.dn.filt,dB = NULL,plot=F,wl=441,norm=T)
      spectro.matrix <- as.matrix(w.dn.filt.spectro$amp)
      
  # Calculate mean amplitude for each column of the spectrogram 
      col.means <- colSums(spectro.matrix)
      col.means.sort <- sort(col.means)
      
  # Plot results
      plot(col.means,type="h",ylab="Sum Amplitude", xlab="Time (10 ms)")
      plot(sort(col.means),type="h",ylab="Sum Amplitude", xlab="Power Order")
      arrows(0,col.means.sort[20], 20, y1 = col.means.sort[20], length = 0.25, angle = 30,lwd=2)
      arrows(0,col.means.sort[70], 70, y1 = col.means.sort[70], length = 0.25, angle = 30,col="red",lwd=2)
  
  # Extract signal and noise values then convert to decibels
      noise <- col.means.sort[20]
      signal <- col.means.sort[70]
      snr.value <- 20*log10(signal/noise)
      snr.value
    }
    
    snr.list[[j]] <- snr.value
    
  }

### Check values in list to see if reasonable (reported in decibels)  
snr.list



#### Part 3 Model Selection code #### 

### Read in data file and check data structure
### Only 1-sec bins with >20 dB are included for analysis
data.bw.r.calc <- read.csv("/Users/denajaneclink/Desktop/performance.constraints.df.08.25.2017.csv",header=T)

#data.subset <- droplevels(subset(data.bw.r.calc, site=="SAF"))
str(data.bw.r.calc)

## Combine model variables and check for multicollinearity (or excessive correlation among explanatory variables)
## Variance inflation factor (VIF) for single explanatory variable is obtained using the r-squared value 
## of the regression of that variable against all other explanatory variables

multi.df <- cbind.data.frame(data.bw.r.calc$trill.rate,data.bw.r.calc$selection,data.bw.r.calc$playback,data.bw.r.calc$bin.new,data.bw.r.calc$trill.dur)
vif(multi.df)


## Check if bandwidth data are normally distributed
hist(data.bw.r.calc$bandwidth)

qqnorm(data.bw.r.calc$bandwidth)
qqline(data.bw.r.calc$bandwidth, col = 2)

qqnorm(log(data.bw.r.calc$bandwidth))
qqline(log(data.bw.r.calc$bandwidth), col = 2)

## Log transformation does slightly better so will log transform bandwidth data along with selection
data.bw.r.calc$bandwidth <- log(data.bw.r.calc$bandwidth)
data.bw.r.calc$selection <- log(data.bw.r.calc$selection)


### Bandwidth Model Selection ###
library(bbmle)
library(lme4)
m0 <-  lmer(bandwidth ~  (1|site/female/great.call), data=data.bw.r.calc)
m1 <- lmer(bandwidth ~ trill.rate +(1|site/female/great.call), data=data.bw.r.calc)
m2 <- lmer(bandwidth ~ bin.new + (1|site/female/great.call), data=data.bw.r.calc)
m3 <- lmer(bandwidth ~ playback +(1|site/female/great.call), data=data.bw.r.calc)
m4 <- lmer(bandwidth ~ selection +(1|site/female/great.call), data=data.bw.r.calc)
m5 <- lmer(bandwidth ~ trill.dur +(1|site/female/great.call), data=data.bw.r.calc)
m6 <- lmer(bandwidth ~ trill.rate +bin.new +(1|site/female/great.call), data=data.bw.r.calc)
m7 <- lmer(bandwidth ~ trill.rate +playback +(1|site/female/great.call), data=data.bw.r.calc)
m8 <- lmer(bandwidth ~ trill.rate +selection +(1|site/female/great.call), data=data.bw.r.calc)
m9 <- lmer(bandwidth ~ trill.rate +trill.dur +(1|site/female/great.call), data=data.bw.r.calc)
m10 <- lmer(bandwidth ~ trill.rate +bin.new +selection +(1|site/female/great.call), data=data.bw.r.calc)
m11 <- lmer(bandwidth ~ trill.rate +bin.new +playback +(1|site/female/great.call), data=data.bw.r.calc)
m12 <- lmer(bandwidth ~ trill.rate +bin.new +trill.dur +(1|site/female/great.call), data=data.bw.r.calc)



table <- print(AICctab(m0,m1,m2,m3, m4,m5,m6,m7,m8,m9,m10,m11,m12,base=T,weights=T,logLik=T, delta=T),min.weight=0.1)

## Model 6 is top model
summary(m6)


## Create coefficient plot for top model for bandwidth
summary <- summary(m6)
model.average.df <- summary$coefficients
variable.vec <- c("Intercept", "Trill Rate (Hz)","Bin")#, "SNR")
model.average.df <- as.data.frame(model.average.df)
model.average.df <- cbind(model.average.df,variable.vec)
colnames(model.average.df) <- c("Coefficient","SE", "t-value", "Variable")

##Remove intercept
model.average.df<- model.average.df[2:3,]

# Specify the width of the confidence intervals
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier

# Coefficient plot for bandwidth
#model.average.df$Variable <- factor(model.average.df$Variable, levels = model.average.df$Variable[order(model.average.df$`t-value`)])
zp1 <- ggplot(model.average.df)
zp1 <- zp1 + geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)
zp1 <- zp1 + geom_linerange(aes(x = Variable, ymin = Coefficient - SE*interval2,
                                ymax = Coefficient + SE*interval2),
                            lwd = 1, position = position_dodge(width = 1/2))
zp1 <- zp1 + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2),
                             lwd = 1/2, position = position_dodge(width = 1/2),
                             shape = 10, fill = "WHITE")
zp1 <- zp1 + coord_flip() + theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

zp1 <- zp1 + ggtitle("Bandwidth") 
zp1 <- zp1 +theme(plot.title = element_text(family = "Trebuchet MS", face="bold", size=24, hjust=0.5))
zp1 <- zp1+theme(axis.text=element_text(size=14)) +theme(axis.title.y=element_blank())+
  theme(axis.title.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=14,face="italic"))
print(zp1)

## Check to see if site-level variance is important for bandwidth
mtop.bandwidth <- lmer(bandwidth ~   trill.rate +bin.new+(1|site/female/great.call) , data=data.bw.r.calc)
mtop.without.bandwidth <- lmer(bandwidth ~ trill.rate +bin.new+(1|female/great.call) , data=data.bw.r.calc)
table <- print(AICctab(mtop.bandwidth,mtop.without.bandwidth,base=T,weights=T,logLik=T, delta=T),min.weight=0.1)




### Trill Rate Model Selection ###
m0 <-  glmer(trill.rate ~ (1|site/female/great.call), data=data.bw.r.calc, family="poisson")
m1 <- glmer(trill.rate ~ bin.new +(1|site/female/great.call), data=data.bw.r.calc, family="poisson")
m2 <- glmer(trill.rate ~ playback+ (1|site/female/great.call), data=data.bw.r.calc, family="poisson")
m3 <- glmer(trill.rate ~ selection+ (1|site/female/great.call), data=data.bw.r.calc, family="poisson")
m4 <- glmer(trill.rate ~ bin.new + selection+ playback+ (1|site/female/great.call), data=data.bw.r.calc, family="poisson")
m5 <- glmer(trill.rate ~ selection + playback +(1|site/female/great.call) , data=data.bw.r.calc, family="poisson")
m6 <- glmer(trill.rate ~ bin.new + selection + trill.dur +(1|site/female/great.call) , data=data.bw.r.calc, family="poisson")
m7 <- glmer(trill.rate ~ trill.dur + (1|site/female/great.call), data=data.bw.r.calc, family="poisson")
m8 <- glmer(trill.rate ~ bin.new  + trill.dur +(1|site/female/great.call) , data=data.bw.r.calc, family="poisson")
m9 <- glmer(trill.rate ~ bin.new + selection+(1|site/female/great.call) , data=data.bw.r.calc, family="poisson")



table.trill <- print(AICctab(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,base=T,weights=T,logLik=T, delta=T),min.weight=0.1)



## Model 8 is top model

model.average.trill <- summary(m8)
model.average.trill <-model.average.trill$coefficients
#model.average.trill <- summary(summary.trill)$coefmat.full

variable.vec <- c("Intercept","Bin","Trill Duration")#,"Playback")
model.average.trill <- as.data.frame(model.average.trill)
model.average.trill <- cbind(model.average.trill,variable.vec)
colnames(model.average.trill) <- c("Coefficient","SE", "z-value", "p-value", "Variable")
model.average.trill<- model.average.trill[2:3,]

model.average.trill$Variable <- factor(model.average.trill$Variable, levels = model.average.trill$Variable[order(model.average.trill$`p-value`)])

perf.const <-ggplot(model.average.trill)
perf.const <- perf.const + geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)
perf.const <- perf.const + geom_linerange(aes(x = Variable, ymin = Coefficient - SE*interval2,
                                              ymax = Coefficient + SE*interval2),
                                          lwd = 1, position = position_dodge(width = 1/16))
perf.const <- perf.const + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                               ymax = Coefficient + SE*interval2),
                                           lwd = 1/2, position = position_dodge(width = 1/4),
                                           shape = 10, fill = "WHITE")
perf.const <- perf.const + coord_flip() + theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
perf.const <- perf.const + ggtitle("Trill Rate")
perf.const <- perf.const + scale_x_discrete(position = "top")
#perf.const <- perf.const + scale_x_discrete(labels=c("Bin","Playback","Sequence in \n Calling Bout  "))
perf.const <- perf.const +theme(axis.text=element_text(size=12))
perf.const <- perf.const +theme(axis.title.y=element_blank())
perf.const <- perf.const +theme(plot.title = element_text(family = "Trebuchet MS", face="bold", size=24, hjust=0.5))
perf.const <- perf.const+theme(axis.text=element_text(size=14)) +theme(axis.title.y=element_blank())+
  theme(axis.title.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=14,face="italic"))

print(perf.const)  

## Is site-level variance important for trill rate?
mtop.trill <- glmer(trill.rate ~ bin.new  + trill.dur+(1|site/female/great.call), data=data.bw.r.calc, family="poisson")
mtop.trill.without <- glmer(trill.rate ~ bin.new  + trill.dur+(1|female/great.call), data=data.bw.r.calc, family="poisson")

site.trill <- print(AICctab(mtop.trill,mtop.trill.without,base=T,weights=T,logLik=T, delta=T),min.weight=0.1)



### Figure 6. Coefficient plot for bandwidth and trill rate
grid.arrange(zp1,perf.const,nrow=1)


### Figure 6. Coefficient plot for bandwidth and trill rate
grid.arrange(zp1, perf.const, nrow = 1)

model1Frame <- data.frame(Variable = rownames(summary(mtop.trill)$coef),
                          Coefficient = summary(mtop.trill)$coef[, 1],
                          SE = summary(mtop.trill)$coef[, 2],
                          modelName = "Trill Rate")
model2Frame <- data.frame(Variable = rownames(summary(mtop.bandwidth)$coef),
                          Coefficient = summary(mtop.bandwidth)$coef[, 1],
                          SE = summary(mtop.bandwidth)$coef[, 2],
                          modelName = "Bandwidth")

#Combine these data.frames
allModelFrame <- data.frame(rbind(model1Frame, model2Frame))  # etc.

allModelFrame <-allModelFrame[-c(1, 4), ] 
#variable.vec <- c("Bin","Trill Duration (s)", "Trill Rate (Hz)","Bin")
allModelFrame$Variable <- as.factor(variable.vec)


allModelFrame$Variable <- factor(allModelFrame$Variable, 
                                 levels = allModelFrame$Variable[order(allModelFrame$Coefficient)])


# Specify the width of your confidence intervals
interval1 <- -qnorm((1-0.95)/2)  # 90% multiplier
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier

# Plot
zp1 <- ggplot(allModelFrame, aes(colour = modelName))
zp1 <- zp1 + geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)+scale_colour_manual(name="Response Variable",values = c("red","blue"))
zp1 <- zp1 + geom_linerange(aes(x = Variable, ymin = Coefficient - SE*interval1,
                                ymax = Coefficient + SE*interval1),
                            lwd = 2, position = position_dodge(width = 1/2))
zp1 <- zp1 + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                 ymax = Coefficient + SE*interval2),
                             lwd = 1/2, position = position_dodge(width = 1/2),
                             shape = 21, fill = "WHITE")+
  scale_x_discrete(labels=c("Intro Duration (s)" ,"Trill Rate (Hz)", "Bin"))

zp1 <- zp1  + theme_bw() +
  coord_flip() + theme_bw()+theme(axis.text=element_text(size=18)) +theme(axis.title.y=element_blank())+theme(axis.title.x=element_text(size=18))+
  theme(axis.text.y=element_text(size=18,face="italic"))+ theme(legend.text=element_text(size=18))+ theme(legend.title =element_text(size=18, face="bold"))

print(zp1)  # The trick to these is position_dodge().


x <- data.bw.r.calc$trill.rate ## This is actually trill rate
y <- data.bw.r.calc$trill.dur


my.spar <- 1 ## this is the smoothing parameter for all spline bootstraps
sp.frame <- data.frame(x=x,y=y)
sp.resampler <- function() {
  n <- nrow(sp.frame)
  resample.rows <- sample(1:n,size=n,replace=TRUE)
  return(sp.frame[resample.rows,])
}

sp.spline.estimator <- function(data,m=300) {
  # Fit spline to data, with cross-validation to pick lambda
  fit <- smooth.spline(x=data[,1],y=data[,2],spar=my.spar)
  # Set up a grid of m evenly-spaced points on which to evaluate the spline
  eval.grid <- seq(from=min(x),to=max(x),length.out=m)
  # Slightly inefficient to re-define the same grid every time we call this,
  # but not a big overhead
  # Do the prediction and return the predicted values
  return(predict(fit,x=eval.grid)$y) # We only want the predicted values
}

sp.spline.cis <- function(B,alpha,m=300) {
  spline.main <- sp.spline.estimator(sp.frame,m=m)
  # Draw B boottrap samples, fit the spline to each
  spline.boots <- replicate(B,sp.spline.estimator(sp.resampler(),m=m))
  # Result has m rows and B columns
  cis.lower <- 2*spline.main - apply(spline.boots,1,quantile,probs=1-alpha/2)
  cis.upper <- 2*spline.main - apply(spline.boots,1,quantile,probs=alpha/2)
  return(list(main.curve=spline.main,lower.ci=cis.lower,upper.ci=cis.upper,
              x=seq(from=min(x),to=max(x),length.out=m)))
}

sp.cis <- sp.spline.cis(B=1000,alpha=0.05)

plot(x,y,xlab="Trill Duration (s)",
     ylab="Trill Rate",col= "grey80", pch =".", cex=3, ylim=c(min(y),max(y)))
#smooth.spline(y,x, cv=FALSE)
lines(x=sp.cis$x,y=sp.cis$main.curve, col="red")
lines(x=sp.cis$x,y=sp.cis$lower.ci)
lines(x=sp.cis$x,y=sp.cis$upper.ci)


str(data.bw.r.calc)








## Figure 7. Plotting site-level differences in bandwidth and trill rate
mtop.bandwidth <- lmer(bandwidth ~ trill.rate + bin.new +(1|site/female/great.call) , data=data.bw.r.calc)
mtop.bandwidth.ranef <- ranef(mtop.bandwidth,condVar=TRUE)
colnames(mtop.bandwidth.ranef[[3]]) <- c("Bandwidth")
bandwidth.re <- dotplot(mtop.bandwidth.ranef,main=F,ylab="",scales=list(x=list(cex=1),y=list(cex=1)))

bandwidth.re$site[[35]][[1]]$y <- factor(bandwidth.re$site[[35]][[1]]$y,
                                         levels = c("DK", "KB", "SAF", "CR","DV","MB","IC"))
bandwidth.re$site[[33]]$labels <- c("Dermakot", "Kinabatangan", "Kalabakan", "Crocker \n Range  ","Danum \n Valley  ","Maliau \n Basin  ","Imbak \n Canyon ")


m1.top.tr <- glmer(trill.rate ~ bin.new +(1|site/female/great.call), data=data.bw.r.calc, family="poisson")
model <- ranef(m1.top.tr,condVar=TRUE)
colnames(model[[3]]) <- c("Trill Rate") 
trill.re <- dotplot(model,main=F,ylab="",scales=list(x=list(cex=1)))
trill.re$site[[33]]$labels <- c("", "", "", "","","","")

grid.arrange(bandwidth.re$site,trill.re$site,nrow=1)



##Create figure by female
mean.bw <- aggregate(data.bw.r.calc$bandwidth, list(group = data.bw.r.calc$female),mean, na.rm=T)
mean.sd <- aggregate(data.bw.r.calc$bandwidth, list(group = data.bw.r.calc$female),sd, na.rm=T)
n.calls <- aggregate(data.bw.r.calc$bandwidth, list(group = data.bw.r.calc$female),sd, na.rm=T)
mean.std.err <-aggregate(data.bw.r.calc$bandwidth, list(group = data.bw.r.calc$female),FUN = function(x) c(mean = mean(x), sd = sd(x),
                                                                                                           n = length(x)))
trill.rate.for.ci <- aggregate(data.bw.r.calc$trill.rate, list(group = data.bw.r.calc$female), FUN = function(x) c(mean = mean(x), sd = sd(x),
                                                                                                                   n = length(x)))
myData.bw <- cbind.data.frame(mean.bw, mean.sd$x, mean.std.err$x[,3])
names(myData.bw) <- c("group","bandwidth","sd","n")
myData.bw$se <- myData.bw$sd / sqrt(myData.bw$n)

myData.trill.rate<- do.call(data.frame, trill.rate.for.ci)
myData.trill.rate$se <- myData.trill.rate$x.sd / sqrt(myData.trill.rate$x.n)

library(quantreg)
performance.constraints.df <- cbind.data.frame(myData.bw,myData.trill.rate)
performance.constraints.df <- subset(performance.constraints.df,!group=="SAFMAT_0B")
names(performance.constraints.df)<- c("group","bandwidth", "bandwidth.sd","bandwidth.n","bandwidth.se","group1","trill.rate", "trill.rate.sd","trill.rate.n","trill.rate.se")
mod <- rq(performance.constraints.df$bandwidth~ performance.constraints.df$trill.rate, data=data.bw.r.calc, tau = 0.90)
summary(rq(performance.constraints.df$bandwidth~ performance.constraints.df$trill.rate, data=data.bw.r.calc, tau = 0.90))
summ <- summary(mod,se="boot",R=10000)
summ

str(performance.constraints.df)
site.perf.df <- str_split_fixed(performance.constraints.df$group,pattern="_", n=2)[,1]



site.perf.df[site.perf.df == "MBGG"] <- "MB"
site.perf.df[site.perf.df == "DVGG"] <- "DV"
site.perf.df[site.perf.df == "SAFA"] <- "SAF"
site.perf.df[site.perf.df == "SAFB100"] <- "SAF"
site.perf.df[site.perf.df == "SAFBA"] <- "SAF"
site.perf.df[site.perf.df == "SAFBB"] <- "SAF"
site.perf.df[site.perf.df == "SAFBC"] <- "SAF"
site.perf.df[site.perf.df == "SAFBD"] <- "SAF"
site.perf.df[site.perf.df == "SAFBN"] <- "SAF"
site.perf.df[site.perf.df == "SAFC"] <- "SAF"
site.perf.df[site.perf.df == "SAFCHIL"] <- "SAF"
site.perf.df[site.perf.df == "SAFCRIP"] <- "SAF"
site.perf.df[site.perf.df == "SAFF"] <- "SAF"
site.perf.df[site.perf.df == "SAFL"] <- "SAF"
site.perf.df[site.perf.df == "SAFMAT"] <- "SAF"
site.perf.df[site.perf.df == "VJRN"] <- "SAF"
site.perf.df[site.perf.df == "VJRS"] <- "SAF"

table(site.perf.df)

### Code to create plot in ggplot
#tiff(filename = "performance.constraints.tiff",
#     width = 480, height = 480, units = "px", pointsize = 12)
trill.low <- subset(performance.constraints.df, group == "IC_07")
trill.high <- subset(performance.constraints.df, group == "DK_10")

library(ggplot2)
library(viridis)


cbbPalette <- viridis(n=7)
cbbPalette[1] <- "grey"
cbbPalette[3] <- "orange"
cbbPalette[5] <- "blue"
cbbPalette[7] <- "orange3"

ggplot(performance.constraints.df, aes(colour=site,trill.rate, bandwidth))+
#ggplot(performance.constraints.df, aes(trill.rate, bandwidth))+
  #ggplot(performance.constraints.df, aes(trill.rate, bandwidth))+
  geom_point(size=6)+
  scale_colour_manual(values=cbbPalette)+
  #scale_shape_manual(values= c(1,2,3,4,5,6,7))+
  #scale_color_viridis(discrete=T)+
  #stat_ellipse(l=0.9, linetype=2)+
  #geom_point(aes(colour=site,size=7))+
  geom_point(data=trill.low, colour=cbbPalette[4],size=12,pch=24,bg=cbbPalette[4]) +
  geom_point(data=trill.high, colour=cbbPalette[2],size=12,pch=24,bg=cbbPalette[2]) +
  #geom_text(aes(label=group))+
  geom_errorbar(width=0.001,aes(ymax = performance.constraints.df$bandwidth+performance.constraints.df$bandwidth.se, ymin=performance.constraints.df$bandwidth-performance.constraints.df$bandwidth.se))+
  geom_errorbarh(aes(xmax = performance.constraints.df$trill.rate+performance.constraints.df$trill.rate.se, xmin=performance.constraints.df$trill.rate-performance.constraints.df$trill.rate.se))+
  geom_abline(,colour="black", lwd=1.5,intercept=coef(mod)[1], slope=coef(mod)[2])+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("Bandwidth (Hz)") +xlab("Trill Rate (Hz)")+
  xlim(5,10)+   guides(color=guide_legend(title="Site"))+
  theme(text = element_text(size=20),axis.text.y= element_text(size=20),axis.text.x= element_text(size=20)) 



feature.data <- read.csv("/Users/denajaneclink/Desktop/features.all.females.04.26.2017.csv", header=T)
feature.data$call.dur <- feature.data$trill.dur + feature.data$trill.dur
plot(lda.data$unlist.trill.rate.list.~lda.data$trill.dur, xlab="Intro Duration", ylab="Trill Rate")
abline(lm(lda.data$unlist.trill.rate.list.~lda.data$trill.dur), col="red")


## Site-level differences in center frequency??

center.freq.df <- read.csv("/Users/denajaneclink/downloads/trill.data.with.center.freq.csv")
str(center.freq.df)

great.call.index <- unique(center.freq.df$great.call)

## Create standardized time column for binning
list.time.std <- sapply(1:length(great.call.index), function(j){
  grp.index<- great.call.index[j]
  data.subset <- subset(center.freq.df,(great.call==grp.index))
  start.time<- data.subset$Begin.Time..s.[1]
  data.subset$time.std.new <-  data.subset$Begin.Time..s.- start.time
  time.standard <- data.subset$time.std.new 
  time.standard <- as.data.frame(time.standard)
  do.call("rbind", time.standard)})

time.std.all <- unlist(list.time.std)
center.freq.df$time.std.new <- time.std.all

mean.cent.freq <- aggregate(center.freq.df$Center.Freq..Hz., 
                     list(great.call = center.freq.df$great.call, site=center.freq.df$site, female=center.freq.df$group,  cut(center.freq.df$time.std.new, 
                                                                      breaks=c(0, 1, 2, 3, 4, 5, 6, 7, 8,9,10,11,12,13,14), 
                                                                      right=T)),
                     mean,na.rm=F)



trill.rate<- aggregate(center.freq.df$time.std.new, 
                       list(group=center.freq.df$great.call, great.call = center.freq.df$great.call,cut(center.freq.df$time.std.new, 
                                                                                                        breaks=c(0, 1, 2, 3, 4, 5, 6, 7, 8,9,10,11,12,13,14), 
                                                                                                        right=T)),
                       length)

str(mean.cent.freq)
new.trill.data.combined <-  cbind.data.frame(mean.cent.freq$site, mean.cent.freq$female, mean.cent.freq$great.call,mean.cent.freq$x,mean.cent.freq$x,trill.rate$x)
colnames(new.trill.data.combined) <- c("site", "female", "great.call", "bin","center.freq","trill.rate")
str(new.trill.data.combined)

m0 <-  lmer(center.freq ~ trill.rate+ (1|female/great.call), data=new.trill.data.combined)
m1 <- lmer(center.freq ~ trill.rate +(1|site/female/great.call), data=new.trill.data.combined)

print(AICctab(m0,m1,base=T,weights=T,logLik=T, delta=T),min.weight=0.1)



center.freq.ranef <- ranef(m1,condVar=TRUE)
colnames(center.freq.ranef[[3]]) <- c("Center Freq")
center.freq.ranef <- dotplot(center.freq.ranef,main=F,ylab="",scales=list(x=list(cex=1),y=list(cex=1)))

center.freq.ranef$site[[35]][[1]]$y <- factor(center.freq.ranef$site[[35]][[1]]$y,
                                         levels = c("DK", "KB", "SAF", "CR","DV","MB","IC"))
center.freq.ranef$site[[33]]$labels <- c("Dermakot", "Kinabatangan", "Kalabakan", "Crocker \n Range  ","Danum \n Valley  ","Maliau \n Basin  ","Imbak \n Canyon ")
plot(center.freq.ranef$site)


cor.test(new.trill.data.combined$center.freq,new.trill.data.combined$trill.rate)

sjp.lmer(m1, type="fe.resid")




##Create figure by female
mean.cf <- aggregate(center.freq.df$Center.Freq..Hz., list(group = center.freq.df$group),mean, na.rm=T)
mean.cf.sd <- aggregate(center.freq.df$Center.Freq..Hz., list(group = center.freq.df$Center.Freq..Hz.),sd, na.rm=T)
n.calls <- aggregate(center.freq.df$Center.Freq..Hz., list(group = center.freq.df$Center.Freq..Hz.),sd, na.rm=T)
mean.cf.std.err <-aggregate(center.freq.df$Center.Freq..Hz., list(group = center.freq.df$Center.Freq..Hz.),FUN = function(x) c(mean = mean(x), sd = sd(x),
                                                                                                           n = length(x)))
trill.rate.for.ci <- aggregate(center.freq.df$trill.rate, 
                               list(group = center.freq.df$group), FUN = function(x) 
                                 c(mean = mean(x), sd = sd(x),
                                                                                                                   n = length(x)))
myData.bw <- cbind.data.frame(mean.cf, mean.cf$x, mean.cf.std.err$x[,3])
names(myData.bw) <- c("group","center.freq","sd","n", "se")
myData.bw$se <- myData.bw$sd / sqrt(myData.bw$n)

myData.trill.rate<- do.call(data.frame, trill.rate.for.ci)
myData.trill.rate$se <- myData.trill.rate$x.sd / sqrt(myData.trill.rate$x.n)

library(quantreg)
performance.constraints.df <- cbind.data.frame(myData.bw,myData.trill.rate)

names(performance.constraints.df)<- c("group","center.freq", "center.freq.sd","center.freq.n",
                                      "center.freq.se","group1","trill.rate", "trill.rate.sd",
                                      "trill.rate.n","trill.rate.se")
mod <- rq(performance.constraints.df$center.freq~ performance.constraints.df$trill.rate, data=performance.constraints.df, tau = 0.90)
summary(rq(performance.constraints.df$center.freq~ performance.constraints.df$trill.rate, data=performance.constraints.df, tau = 0.90))
summ <- summary(mod,se="boot",R=10000)
summ

str(performance.constraints.df)
site.perf.df <- str_split_fixed(performance.constraints.df$group,pattern="_", n=2)[,1]



site.perf.df[site.perf.df == "MBGG"] <- "MB"
site.perf.df[site.perf.df == "DVGG"] <- "DV"
site.perf.df[site.perf.df == "SAFA"] <- "SAF"
site.perf.df[site.perf.df == "SAFB100"] <- "SAF"
site.perf.df[site.perf.df == "SAFBA"] <- "SAF"
site.perf.df[site.perf.df == "SAFBB"] <- "SAF"
site.perf.df[site.perf.df == "SAFBC"] <- "SAF"
site.perf.df[site.perf.df == "SAFBD"] <- "SAF"
site.perf.df[site.perf.df == "SAFBN"] <- "SAF"
site.perf.df[site.perf.df == "SAFC"] <- "SAF"
site.perf.df[site.perf.df == "SAFCHIL"] <- "SAF"
site.perf.df[site.perf.df == "SAFCRIP"] <- "SAF"
site.perf.df[site.perf.df == "SAFF"] <- "SAF"
site.perf.df[site.perf.df == "SAFL"] <- "SAF"
site.perf.df[site.perf.df == "SAFMAT"] <- "SAF"
site.perf.df[site.perf.df == "VJRN"] <- "SAF"
site.perf.df[site.perf.df == "VJRS"] <- "SAF"

performance.constraints.df$site <- (site.perf.df)

### Code to create plot in ggplot
#tiff(filename = "performance.constraints.tiff",
#     width = 480, height = 480, units = "px", pointsize = 12)
trill.low <- subset(performance.constraints.df, group == "IC_07")
trill.high <- subset(performance.constraints.df, group == "DK_10")

library(ggplot2)
library(viridis)


cbbPalette <- viridis(n=7)
cbbPalette[1] <- "grey"
cbbPalette[3] <- "orange"
cbbPalette[5] <- "blue"
cbbPalette[7] <- "orange3"

str(performance.constraints.df)

ggplot(performance.constraints.df, aes(trill.rate, center.freq,colour=site))+
  #ggplot(performance.constraints.df, aes(trill.rate, bandwidth))+
  #ggplot(performance.constraints.df, aes(trill.rate, bandwidth))+
  geom_point(size=6)+
  scale_colour_manual(values=cbbPalette)+
  #scale_shape_manual(values= c(1,2,3,4,5,6,7))+
  #scale_color_viridis(discrete=T)+
  #stat_ellipse(l=0.9, linetype=2)+
  #geom_point(aes(colour=site,size=7))+
  geom_point(data=trill.low, colour=cbbPalette[4],size=12,pch=24,bg=cbbPalette[4]) +
  geom_point(data=trill.high, colour=cbbPalette[2],size=12,pch=24,bg=cbbPalette[2]) +
  #geom_text(aes(label=group))+
  geom_errorbar(width=0.001,aes(ymax = performance.constraints.df$center.freq+performance.constraints.df$center.freq.se, ymin=performance.constraints.df$center.freq-performance.constraints.df$center.freq.se))+
  geom_errorbarh(aes(xmax = performance.constraints.df$trill.rate+performance.constraints.df$trill.rate.se, xmin=performance.constraints.df$trill.rate-performance.constraints.df$trill.rate.se))+
  geom_abline(,colour="black", lwd=1.5,intercept=coef(mod)[1], slope=coef(mod)[2])+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("Bandwidth (Hz)") +xlab("Trill Rate (Hz)")+
  xlim(5,10)+   guides(color=guide_legend(title="Site"))+
  theme(text = element_text(size=20),axis.text.y= element_text(size=20),axis.text.x= element_text(size=20)) 






