#### R code for Clink et al. Evidence of vocal performance constraints in a female non-human primate 
#### Part 1 Bandwidth calculation
#### Part 2 Signal-to-noise ratio calculation
#### Part 3 Quantile Regression and Model selection code

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
library(quantreg)
library(MASS)
library(car)
library(viridis)


#### Part 1 Bandwidth calculation  #### 
#Set input directory to read in .wav files
input.dir <-"/Users/denajaneclink/Desktop/Performance Constraints/Files_for_R_Online_Supporting_Material/TrillNotes"
L = list.files(input.dir, pattern = "*.wav", full.names = FALSE)
L
filehandles <- str_split_fixed(L, pattern = ".wav", n = 2)[, 1]

# Create an empty list that will store bandwidth values
bandwidth.list = list()


# Function to calculate bandwidth and create Figure 3

bandwidth.function <- function (spec, f = NULL, level = -3, 
                                plot = TRUE, colval = "black",
                                cexval = 1, fontval = 1, 
                                flab = "Frequency (Hz)", 
                                alab = "Relative amplitude (dB)", 
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
  range <- c(f/2000/length(spec), f / 2000)
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
  f1khz <- ((f1[length(f1)] / n1) * (range[2] - range[1])) +
    range[1]
  f2 <- which(round(specB, 1) == level2) + (nA - 1)
  f2khz <- ((f2[1] / n1) * (range[2] - range[1])) + range[1]
  Q <- f2khz - f1khz
  results <- list(Q = Q, bdw = f2khz - f1khz)
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
      filename <- paste(input.dir, "/", filehandle, sep = "")
      print(paste("processing", filehandle))
    
      # Read in wav file and downsample (which also filters 
      # frequencies of 4 kHz and above)
      wav <- readWave(filename)
    
      # Bandpass filter to prevent aliasing during downsampling
      w.dn.filt <- fir(wav, from=450, to=3000,output = "Wave")
      w.dn.filt <- normalize(w.dn.filt, unit="16")
      wav@left <-w.dn.filt@left
      
      # Downsample
      wav <-downsample(wav, samp.rate=5000)
      str(wav)
      # Calculate number of samples in .wav file
      N <- length(wav@left) 
      
      # Pad with zeros to achieve length 5000.
      waveform <- c(wav@left, rep(0, 5000 - N))
      
      # Re-insert into wave object. 
      wav@left <- waveform   
      
      #Calculate the power spectral density which provdes a 1-hertz resolution
      psd <- spec(wav,f=5000, plot=T, dB="max0",wl=5000)
      
      # Calculate the moving average to smooth out the psd
      frequency <- rollapply((seq(from=1,to=nrow(psd),by=1)),width=50, by=50,median) # calculate the moving average
      amplitude <- rollapply(psd[,2],width=50, by=50, mean)
      
      # Normalize the amplitude
      amplitude <- amplitude-max(amplitude)
      
      # Combine by column into new object that can be used to calculated bandwidth
      psd.comb<-cbind(frequency,amplitude)
      bandwidth <- bandwidth.function(psd.comb, level = -15, plot = T)$bdw
      
    # Print each value once calculated  
      print(bandwidth)
      
    # Add to the list object created earlier
      bandwidth.list[[j]] <- bandwidth
}

# Check bandwidth list to make sure values are reasonable
bandwidth.list



#### Part 2 Signal-to-noise ratio calculation #### 

    # Read in data file
    snr.data <- read.csv("/Users/denajaneclink/Desktop/Performance Constraints/Files_for_R_Online_Supporting_Material/snr.data.subset.csv",header = T)
    str(snr.data)
    
    # Set input directory to the file with the .wav files
    input.dir.new <- "/Users/denajaneclink/Desktop/Performance Constraints/Performance_constraints_wav_files/long_wav_files_for_snr"
    L.snr = list.files(input.dir.new, pattern="*.wav", full.names=T)
    L.snr
    
    # Create empty list 
    noise.list.long <- list()

    # Loop to find candidate noise from individual recording of gibbon great calls
    for (j in 1:length(L.snr)) { tryCatch({
      
      # Read in two minute .wav file
      snr.file <- L.snr[[j]]
      snr.wav <- readWave(snr.file)
      print(paste("processing", snr.file))
      
      # Filter the .wav file so focus on the gibbon frequency range 
      w.dn.filt <- fir(snr.wav, from=450, to=1800,output = "Wave") 
      w.dn.filt <- normalize(w.dn.filt, unit="16")
      
      # Create 1-sec bins 
      bin.seq <- seq(from=0, to=120, by=1)
      length(bin.seq)
      
      # Create a list of all the 1-s bins for SNR analysis
      bin.seq.length <- length(bin.seq)-1
      subsamps.1sec <- lapply(1:bin.seq.length, function(i) extractWave(w.dn.filt, from=as.numeric(bin.seq[i]), to=as.numeric(bin.seq[i+1]), xunit = c("time"),plot=F,output="Wave"))
      
      # Calculate noise for each 1 sec bin for the longer recording
      noise.list <- list()
      for (k in 1:length(subsamps.1sec)) { 
        
        # Read in .wav file 
        wav <- subsamps.1sec[[k]]
        
        bin.seq.smaller <- seq(from=0, to=44100, by=441)
        bin.seq.length.small <- length(bin.seq.smaller)-1
        subsamps.short <- lapply(1:bin.seq.length.small, function(i) 
          extractWave(wav, from=as.numeric(bin.seq.smaller[i]), 
                      to=as.numeric(bin.seq.smaller[i+1]), xunit = c("samples"),
                      plot=F,output="Wave"))
        
        # Calculate sum of squares for each 10 ms long sample
        sum.square.list <- list ()
        
        for (f in 1:length(subsamps.short)){
          wav.temp <- subsamps.short[[f]]
          sum.square.list[[f]] <- sum(wav.temp@left^2)
        }
        
        #Find median value for each 10 ms clip
        noise.list[[k]] <- median(unlist(sum.square.list))
        
      }
      
      ## Now have a list with 120 values corresponding to median value for each 10 ms subsample
      ## that are candidate noise; use the 10th percentile of the distribution for noise estimate
      noise.value <- quantile(unlist(noise.list), c(.10))
      noise.value
      
      ## Combine recording id and noise estimate value 
      great.call <- str_split_fixed(snr.file, pattern="/", n=8)[[8]]
      great.call <- gsub(".wav", "", great.call)
 
      ## Add to list of noise values
      noise.list.long[[j]] <-cbind.data.frame(great.call,noise.value)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }

    ## Convert the list of noise values to a dataframe
    noise.data.frame <- do.call(rbind.data.frame, lapply(noise.list.long, data.frame, stringsAsFactors=FALSE))
    str(noise.data.frame)


    ##Set wd to location of wave files and Raven .txt files
    setwd("/Users/denajaneclink/Desktop/Performance Constraints/Files_for_R_Online_Supporting_Material/wav_and_txt_files")
    
    # Create empty list to populate with signal values
    snr.list <- list()
    
    # Create index of great calls
    great.call.index <- unique(snr.data$great.call)
    
    # Loop to calculate signal candidates
    for (j in 1:length(great.call.index)) {                                      
      great.call.in <- great.call.index[j]
      print(paste("processing file", great.call.in)) 
      newdata <- subset(snr.data, great.call==great.call.in)
      great.call.in <- paste(great.call.in,".wav",sep="")
      wav <- readWave(great.call.in)
      
      # defining the break points for 1-sec bins for the trill
      trill.start.time <- min(newdata$Begin.Time..s.)
      trill.last.bin <- max(newdata$End.Time..s.)
      n.bins <- as.integer(trill.last.bin-trill.start.time,length=0)
      bin.seq <- seq(from=0, to=n.bins,by=1 )
      bin.seq <- trill.start.time + bin.seq
      
      # create a list of all the 1-sec segments
      bin.seq.length <- length(bin.seq)-1
      subsamps <- lapply(1:bin.seq.length, function(i) extractWave(wav, from=as.numeric(bin.seq[i]), to=as.numeric(bin.seq[i+1]), xunit = c("time"),plot=F,output="Wave"))
      
      # Create an empty list for the signal candidates
      signal.list <- list()
      
      # Loop to calculate power estimate for each 10 ms time slice
      for (k in 1:length(subsamps)) { 
        # Read in .wav file 
        wav <- subsamps[[k]]
        
        # Cut 1-sec .wav files into 10 ms .wav files
        bin.seq.smaller <- seq(from=0, to=44100, by=441)
        bin.seq.length.small <- length(bin.seq.smaller)-1
        subsamps.short <- lapply(1:bin.seq.length.small, function(i) 
          extractWave(wav, from=as.numeric(bin.seq.smaller[i]), 
                      to=as.numeric(bin.seq.smaller[i+1]), xunit = c("samples"),
                      plot=F,output="Wave"))
        
        # Create empty list to store signal values
        sum.square.list <- list ()
        
        # Calculate sum of squares for each 10 ms long sample
        for (f in 1:length(subsamps.short)){
          wav.temp <- subsamps.short[[f]]
          #sum.square.list[[f]] <- rms(wav.temp@left)
          sum.square.list[[f]] <- sum(wav.temp@left^2)
        }
        
        #Find 75% quantile value for each 10 ms clip
        signal.list[[k]] <- quantile(unlist(sum.square.list), c(.75))
        
      }
      
      # Create a dataframe for the signal values
      signal.df <- do.call(rbind.data.frame, lapply(signal.list, data.frame, stringsAsFactors=FALSE))
      colnames(signal.df) <- c("Signal")
      
      # Match noise measure with signal measure
      ## Subtract estimated noise
      call.in <- great.call.index[j]
      call.in <- as.character(droplevels(call.in))
      noise.measure <- subset(noise.data.frame, great.call== call.in)
      
      # The signal value calculated above is actually the signal + noise, therefore need to subtract the noise 
      bins.minus.noise <- signal.df$Signal - noise.measure$noise.value
      
      #Finally, calculate the SNR value for each 1-sec bin
      SNR.dB <- 20 * log10(bins.minus.noise/noise.measure$noise.value)
      
      # add back to subsetted data
      # subsamp.index <- length(subsamps)+1
      data.subset.list <- lapply(1:length(SNR.dB), function(i) 
        subset(newdata,time.std.new >= i-1 & time.std.new <= i))
      
      new.data.subset.list.with.snr = list()
      for (m in 1:length(data.subset.list)) { 
        list.element <- as.data.frame(data.subset.list[m])
        snr.value <- SNR.dB[m]
        new.snr.vector <-c(rep(snr.value, nrow(list.element)))
        list.element$snr <- new.snr.vector
        new.data.subset.list.with.snr[[m]] <- list.element
      }
      
      # Add to the list created above
      newdata.snr <- do.call(rbind.data.frame, lapply(new.data.subset.list.with.snr, data.frame, stringsAsFactors=FALSE))
      #new.data.compiled <- do.call("cbind.data.frame", newdata.snr) 
      snr.list[[j]] <- newdata.snr
    }
    
    # Convert the list with SNR values to a data frame
    snr.dataframe<- do.call(rbind.data.frame, lapply(snr.list, data.frame, stringsAsFactors=FALSE))
    snr.dataframe

#### Part 3 Quantile Regression and Model Selection code #### 

### Read in data file and check data structure
### Only 1-sec bins with >20 dB are included for analysis
data.bw.r.calc <- read.csv("/Users/denajaneclink/Desktop/Performance Constraints/Files_for_R_Online_Supporting_Material/performance.constraints.df.10.09.2017.csv", header=T)
    

## Check number of sites, females and great calls

length(unique(data.bw.r.calc$site))
length(unique(data.bw.r.calc$female))
length(unique(data.bw.r.calc$great.call))


### Code to run quantile regression

# Calculate mean bandwidth and trill rate per female
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

# Combine mean bandwidth and trill rate into dataframe
performance.constraints.df <- cbind.data.frame(myData.bw,myData.trill.rate)
names(performance.constraints.df)<- c("group","bandwidth", "bandwidth.sd","bandwidth.n","bandwidth.se","group1","trill.rate", "trill.rate.sd","trill.rate.n","trill.rate.se")

# Calculate the 90% quantile regression
mod <- rq(performance.constraints.df$bandwidth~ performance.constraints.df$trill.rate, data=data.bw.r.calc, tau = 0.90)
summary(rq(performance.constraints.df$bandwidth~ performance.constraints.df$trill.rate, data=data.bw.r.calc, tau = 0.90))

# Calculate the p-value
# NOTE: as this is done using a bootstrap method the values will be slightly different each time.
summ <- summary(mod,se="boot",R=10000)
summ


### Convert file names to site for color-coding
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

performance.constraints.df$site <- as.factor(site.perf.df)

### Code to create plot in ggplot

# Identify high and low performance trills
trill.low <- subset(performance.constraints.df, group == "IC_07")
trill.high <- subset(performance.constraints.df, group == "DK_10")


# Set color palette for figure
cbbPalette <- viridis(n=7)
cbbPalette[1] <- "grey"
cbbPalette[3] <- "orange"
cbbPalette[5] <- "blue"
cbbPalette[7] <- "orange3"


#Create plot using ggplot
ggplot(performance.constraints.df, aes(trill.rate, bandwidth,colour=site))+
  geom_point(size=6)+
  scale_colour_manual(values=cbbPalette)+
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


### Model selection

## Combine model variables and check for multicollinearity 
##  (or excessive correlation among explanatory variables)
## Variance inflation factor (VIF) for single explanatory variable
##  is obtained using the r-squared value 
##  of the regression of that variable against all other explanatory variables
data.bw.r.calc <- na.omit(data.bw.r.calc)
multi.df <- cbind.data.frame(data.bw.r.calc$trill.rate,
                             data.bw.r.calc$selection,
                             data.bw.r.calc$playback,
                             data.bw.r.calc$bin.new)
vif(multi.df)


## Check if bandwidth data are normally distributed
hist(data.bw.r.calc$bandwidth)

qqnorm(data.bw.r.calc$bandwidth)
qqline(data.bw.r.calc$bandwidth, col = 2)

qqnorm(log(data.bw.r.calc$bandwidth))
qqline(log(data.bw.r.calc$bandwidth), col = 2)

## Log transformation does slightly better so will log transform bandwidth and selection data
data.bw.r.calc$bandwidth <- log(data.bw.r.calc$bandwidth)
data.bw.r.calc$selection <- log(data.bw.r.calc$selection)



### Bandwidth Model Selection ###
#If you want to compare models that differ in fixed effects terms, then you must use ordinary likelihood;
# which is why REML=F

m0.bandwidth <-  lmer(bandwidth ~  (1|site/female/great.call), 
                      data=data.bw.r.calc)
m1.bandwidth <- lmer(bandwidth ~ trill.rate +(1|site/female/great.call), 
                     data=data.bw.r.calc)
m2.bandwidth <- lmer(bandwidth ~ bin.new + (1|site/female/great.call), 
                     data=data.bw.r.calc)
m3.bandwidth <- lmer(bandwidth ~ playback +(1|site/female/great.call), 
                     data=data.bw.r.calc)
m4.bandwidth <- lmer(bandwidth ~ selection +(1|site/female/great.call), 
                     data=data.bw.r.calc)
m5.bandwidth <- lmer(bandwidth ~ trill.dur +(1|site/female/great.call), 
                     data=data.bw.r.calc)
m6.bandwidth <- lmer(bandwidth ~ trill.rate +bin.new +(1|site/female/great.call), 
                     data=data.bw.r.calc)
m7.bandwidth <- lmer(bandwidth ~ trill.rate +playback +(1|site/female/great.call), 
                     data=data.bw.r.calc)
m8.bandwidth <- lmer(bandwidth ~ trill.rate +selection +(1|site/female/great.call), 
                     data=data.bw.r.calc)
m9.bandwidth <- lmer(bandwidth ~ trill.rate +trill.dur +(1|site/female/great.call), 
                     data=data.bw.r.calc)
m10.bandwidth <- lmer(bandwidth ~ trill.rate +bin.new +selection +(1|site/female/great.call), 
                      data=data.bw.r.calc)
m11.bandwidth <- lmer(bandwidth ~ trill.rate +bin.new +playback +(1|site/female/great.call), 
                      data=data.bw.r.calc)
m12.bandwidth <- lmer(bandwidth ~ trill.rate +bin.new +trill.dur +(1|site/female/great.call), 
                      data=data.bw.r.calc)

## Compare AIC of bandwidth models
table <- print(AICctab(m0.bandwidth, m1.bandwidth, m2.bandwidth, m3.bandwidth, m4.bandwidth, m5.bandwidth, m6.bandwidth,
                       m7.bandwidth, m8.bandwidth, m9.bandwidth, m10.bandwidth, m11.bandwidth, m12.bandwidth, 
                       base = T, weights = T, logLik = T, 
                       delta = T), min.weight = 0.1)



## Model 6 is top model
summary(m6.bandwidth)


## Create coefficient plot for top model for bandwidth
summary <- summary(m6.bandwidth)
model.average.df <- summary$coefficients
variable.vec <- c("Intercept", "Trill Rate (Hz)","Bin")
model.average.df <- as.data.frame(model.average.df)
model.average.df <- cbind(model.average.df,variable.vec)
colnames(model.average.df) <- c("Coefficient","SE", "t-value", "Variable")

##Remove intercept
model.average.df<- model.average.df[2:3, ]

# Specify the width of the confidence intervals
interval2 <- -qnorm((1 - 0.95) / 2)  # 95% multiplier

# Coefficient plot for bandwidth
model.average.df$Variable <- factor(model.average.df$Variable, 
                                    levels = model.average.df$Variable[order(model.average.df$`t-value`)])
bandwidth.plot <- ggplot(model.average.df)
bandwidth.plot <- bandwidth.plot + geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)
bandwidth.plot <- bandwidth.plot + geom_linerange(aes(x = Variable,
                                ymin = Coefficient - SE * interval2,
                                ymax = Coefficient + SE * interval2),
                                lwd = 1, 
                                position = position_dodge(width = 1 / 2))
bandwidth.plot <- bandwidth.plot + geom_pointrange(aes(x = Variable, y = Coefficient, 
                                ymin = Coefficient - SE * interval2,
                                ymax = Coefficient + SE*interval2),
                                lwd = 1 / 2, 
                                position = position_dodge(width = 1 / 2),
                                shape = 10, fill = "WHITE")
bandwidth.plot <- bandwidth.plot + coord_flip() + theme_bw() +
     theme(panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())

bandwidth.plot <- bandwidth.plot + ggtitle("Bandwidth") 
bandwidth.plot <- bandwidth.plot + theme(plot.title = element_text(family = "Trebuchet MS", 
                                            face="bold", size = 24, 
                                            hjust = 0.5))
bandwidth.plot <- bandwidth.plot + theme(axis.text = element_text(size = 12)) +
   theme(axis.text=element_text(size=14)) +theme(axis.title.y=element_blank())+
  theme(axis.title.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=14,face="italic"))

print(bandwidth.plot)

## Check to see if site-level variance is important for bandwidth
mtop.bandwidth <- lmer(bandwidth ~ trill.rate +bin.new+
                         (1|site/female/great.call), 
                       data = data.bw.r.calc, REML = F)
mtop.without <- lmer(bandwidth ~ trill.rate +bin.new+
                       (1|female/great.call), 
                     data = data.bw.r.calc, REML = F)
table <- print(AICctab(mtop.bandwidth, mtop.without, 
                       base = T, weights = T, logLik = T, 
                       delta = T), min.weight = 0.1)



### Trill Rate Model Selection ###

## Check if poisson distribution is appropriate
poisson <- fitdistr(data.bw.r.calc$trill.rate, "Poisson")
qqp(data.bw.r.calc$trill.rate, "pois", poisson$estimate)


## Model selection code
m0.trill <-  glmer(trill.rate ~ (1|site/female/great.call), 
                   data=data.bw.r.calc, family="poisson")
m1.trill <- glmer(trill.rate ~ bin.new +(1|site/female/great.call), 
                  data=data.bw.r.calc, family="poisson")
m2.trill <- glmer(trill.rate ~ playback+ (1|site/female/great.call), 
                  data=data.bw.r.calc, family="poisson")
m3.trill <- glmer(trill.rate ~ selection+ (1|site/female/great.call), 
                  data=data.bw.r.calc, family="poisson")
m4.trill <- glmer(trill.rate ~ bin.new + selection+ playback+ (1|site/female/great.call), 
                  data=data.bw.r.calc, family="poisson")
m5.trill <- glmer(trill.rate ~ selection + playback +(1|site/female/great.call) , 
                  data=data.bw.r.calc, family="poisson")
m6.trill <- glmer(trill.rate ~ trill.dur + (1|site/female/great.call), 
                  data=data.bw.r.calc, family="poisson")
m7.trill <- glmer(trill.rate ~ bin.new  + trill.dur +(1|site/female/great.call), 
                  data=data.bw.r.calc, family="poisson")
m8.trill <- glmer(trill.rate ~ bin.new + selection+(1|site/female/great.call), 
                  data=data.bw.r.calc, family="poisson")


# Compare models using AIC
table.trill <- print(AICctab(m0.trill,m1.trill,m2.trill,m3.trill,m4.trill,m5.trill,m6.trill,
                             m7.trill,m8.trill,base=T,weights=T,logLik=T, delta=T),min.weight=0.1)



# Model with bin and trill duration is top model
model.average.trill <- summary(m7.trill)

## Code to create coefficient plot for trill rate
model.average.trill <-model.average.trill$coefficients


variable.vec <- c("Intercept","Bin","Trill Duration")
model.average.trill <- as.data.frame(model.average.trill)
model.average.trill <- cbind(model.average.trill,variable.vec)
colnames(model.average.trill) <- c("Coefficient","SE", "z-value", "p-value", "Variable")
model.average.trill<- model.average.trill[2:3,]

model.average.trill$Variable <- factor(model.average.trill$Variable, levels = model.average.trill$Variable[order(model.average.trill$`p-value`)])

trill.plot <-ggplot(model.average.trill)
trill.plot <- trill.plot + geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)
trill.plot <- trill.plot + geom_linerange(aes(x = Variable, ymin = Coefficient - SE*interval2,
                                              ymax = Coefficient + SE*interval2),
                                          lwd = 1, position = position_dodge(width = 1/16))
trill.plot <- trill.plot + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                               ymax = Coefficient + SE*interval2),
                                           lwd = 1/2, position = position_dodge(width = 1/4),
                                           shape = 10, fill = "WHITE")
trill.plot <- trill.plot + coord_flip() + theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
trill.plot <- trill.plot + ggtitle("Trill Rate")
trill.plot <- trill.plot + scale_x_discrete(position = "top")
#trill.plot <- trill.plot + scale_x_discrete(labels=c("Bin","Playback","Sequence in \n Calling Bout  "))
trill.plot <- trill.plot +theme(axis.text=element_text(size=12))
trill.plot <- trill.plot +theme(axis.title.y=element_blank())
trill.plot <- trill.plot +theme(plot.title = element_text(family = "Trebuchet MS", face="bold", size=24, hjust=0.5))
trill.plot <- trill.plot+theme(axis.text=element_text(size=14)) +theme(axis.title.y=element_blank())+
  theme(axis.title.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=14,face="italic"))

print(trill.plot)  

## Is site-level variance important for trill rate?
mtop.trill <- glmer(trill.rate ~ bin.new + trill.dur + (1 | site / female / great.call), 
                    data = data.bw.r.calc, family = "poisson")
mtop.trill.without <- glmer(trill.rate ~ bin.new + trill.dur + (1 | female / great.call),
                            data = data.bw.r.calc, family = "poisson")

site.trill <- print(AICctab(mtop.trill,mtop.trill.without,
                            base = T,weights = T,logLik = T,
                            delta = T), min.weight = 0.1)



### Figure 6. Combine coefficient plots for bandwidth and trill rate
grid.arrange(bandwidth.plot, trill.plot, nrow = 1)


## Figure 7. Plotting site-level differences in bandwidth and trill rate
mtop.bandwidth <- lmer(bandwidth ~ trill.rate +  bin +
                       (1 | site / female / great.call), 
                       data = data.bw.r.calc)
mtop.bandwidth.ranef <- ranef(mtop.bandwidth, condVar = TRUE)
colnames(mtop.bandwidth.ranef[[3]]) <- c("Bandwidth")
bandwidth.re <- dotplot(mtop.bandwidth.ranef,main = F,
                        ylab = "",
                        scales = list(x = list(cex = 1),
                                      y = list(cex = 1)))

bandwidth.re$site[[35]][[1]]$y <- factor(bandwidth.re$site[[35]][[1]]$y,
                                         levels = c("DK", "KB", "SAF", 
                                                    "CR", "DV", "MB", "IC"))
bandwidth.re$site[[33]]$labels <- c("Dermakot", "Kinabatangan", "Kalabakan", 
                                    "Crocker \n Range  ","Danum \n Valley  ",
                                    "Maliau \n Basin  ", "Imbak \n Canyon ")


m1.top.tr <- glmer(trill.rate ~ bin.new + trill.dur + (1 | site / female / great.call),
                   data = data.bw.r.calc, family = "poisson")
model <- ranef(m1.top.tr, condVar = TRUE)
colnames(model[[3]]) <- c("Trill Rate") 
trill.re <- dotplot(model, main = F, ylab = "",
                    scales = list(x = list(cex = 1)))
trill.re$site[[33]]$labels <- c("", "", "", "", "", "","")

# Create plot
grid.arrange(bandwidth.re$site, trill.re$site, nrow = 1)



