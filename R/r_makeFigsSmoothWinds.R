library(ggpubr)
library(Hmisc)
library(gdata)

figsPath = "/Users/leo/rice/tmp/"
paperPath = "/Users/leo/rice/res/papers/RECOMB2016/svn/WABI19-Leo/lipics-v2019-authors/figs/"



###
# Fig X mosq. 6tax timings
###
fileToRead = paste(figsPath,"final_mosq_6tax_20x/FinalTimings/r_all",sep="")
values = read.csv(fileToRead, header=TRUE)

#command line with saving?
myPlot <- ggerrorplot(values, x = "Method", y = "Time", xlab = "Method Used", ylab = "Time (s)",
                      desc_stat = "mean_sd", color = "black", 
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(14,"bold", "black")
)

myPlot <- ggpar(myPlot, x.text.angle = 90) #angles x labels

ggsave(filename=paste(paperPath,"fig_mosq_6tax_time.pdf",sep=""), plot=myPlot, width=8,height=8)



#myPlot <- ggboxplot(values, x = "Method", y = "Time", xlab = "Method Used", ylab = "Time (s)",
#                      desc_stat = "mean_sd", color = "black", 
#                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black"),angle=30
#)
#myPlot <- ggviolin(values, x = "Method", y = "Time", xlab = "Method Used", ylab = "Time (s)",
#                    desc_stat = "mean_sd", color = "black", 
#                    add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black"),angle=30
#)


###
# Fig X mosq. 6tax accuracy
###
fileToRead = paste(figsPath,"final_mosq_6tax_20x/FinalResults/r_all",sep="")
values = read.csv(fileToRead, header=TRUE)

#command line with saving?
myPlot <- ggerrorplot(values, x = "Method", y = "Accuracy", xlab = "Method Used", ylab = "Symmetric Difference",
                      desc_stat = "mean_sd", color = "black", 
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black")
)

myPlot <- ggpar(myPlot, x.text.angle = 90) #angles x labels

# THIS IS AN OPTION ALSO
# myPlot <- myPlot + coord_flip()

ggsave(filename=paste(paperPath,"fig_mosq6tax_accuracy.pdf",sep=""), plot=myPlot, width=8,height=8)


#compare vs trying with just ggplot
fileToRead = paste(figsPath,"final_mosq_6tax_20x/FinalResults/r_all",sep="")
values = read.csv(fileToRead, header=TRUE)



#command line with saving?
myPlot <- ggerrorplot(values, x = "Method", y = "Accuracy", xlab = "Method Used", ylab = "Symmetric Difference",
                      desc_stat = "mean_sd", color = "black", 
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black")
)

myPlot <- ggpar(myPlot, x.text.angle = 90) #angles x labels

ggsave(filename=paste(paperPath,"fig_mosq6tax_accuracy.pdf",sep=""), plot=myPlot, width=8,height=8)



###
# Fig - primates 11tax accuracy
###

fileToRead = paste(figsPath,"final_primates_11tax_20x/FinalResults/r_all",sep="")
values = read.csv(fileToRead, header=TRUE)

#do a lil renaming
#values[1:21,]
values [,1] = gsub("-CNTRD", "", values[,1])
values [,1] = gsub("LessMem", "", values[,1])

#sort em by mean
means = aggregate(values[, 2], list(values$Method), mean)
meansOrdered = means[order(means$x),]
meansOrderedNames = meansOrdered[,1]
list(meansOrderedNames)
meanOrderedNames = fct_inorder(rev(meansOrderedNames)) #rev does reverse order (so currently best is on top, without rev its on bottom)

#meanOrderedNames <- fct_inorder(meanOrderedNames, ordered = TRUE)

meanOrderedNames <- factor(meanOrderedNames, fct_inorder(meanOrderedNames))
#meanOrderedNames
values$Method <- factor(values$Method, meanOrderedNames)
#values$Method #WORKS?

#command line with saving?
myPlot <- ggerrorplot(values, x = "Method", y = "Accuracy", xlab = "Method Used", ylab = "Symmetric Difference",
                      desc_stat = "mean_sd", color = "black", 
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black")
)
myPlot <- ggpar(myPlot, x.text.angle = 90) #angles x labels
# THIS IS AN OPTION ALSO
myPlot <- myPlot + coord_flip()
ggsave(filename=paste(paperPath,"fig_primates_11tax_accuracy.pdf",sep=""), plot=myPlot, width=12,height=8)






###
# Fig - generate ALL ACCURACY FIGS (since the naming conventions are easy and they are all the same)
###

myFunctionMakePlot <- function(inputFile,outputFile) {

fileToRead = paste(figsPath,inputFile,sep="")
values = read.csv(fileToRead, header=TRUE)

#do a lil renaming
values [,1] = gsub("-CNTRD", "", values[,1])
values [,1] = gsub("LessMem", "", values[,1])

#sort em by mean
means = aggregate(values[, 2], list(values$Method), mean)
meansOrdered = means[order(means$x),]
meansOrderedNames = meansOrdered[,1]
meanOrderedNames = fct_inorder(rev(meansOrderedNames)) #rev does reverse order (so currently best is on top, without rev its on bottom)

#meanOrderedNames <- factor(meanOrderedNames, fct_inorder(meanOrderedNames))
values$Method <- factor(values$Method, meanOrderedNames)
#values$Method #WORKS?

#command line with saving?
myPlot <- ggerrorplot(values, x = "Method", y = "Accuracy", xlab = "Method Used", ylab = "Symmetric Difference",
                      desc_stat = "mean_sd", color = "black", 
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black")
)
myPlot <- ggpar(myPlot, x.text.angle = 90) #angles x labels
# THIS IS AN OPTION ALSO (flips it)
myPlot <- myPlot + coord_flip()
ggsave(filename=paste(paperPath,outputFile,sep=""), plot=myPlot, width=12,height=8)
}

myFunctionMakePlot("final_butterfly_8tax_20x/FinalResults/r_all","fig_butterfly_8tax_accuracy.pdf")
myFunctionMakePlot("final_butterfly_12tax_20x/FinalResults/r_all","fig_butterfly_12tax_accuracy.pdf")
myFunctionMakePlot("final_bird_8tax_20x/FinalResults/r_all","fig_bird_8tax_accuracy.pdf")
myFunctionMakePlot("final_bird_12tax_20x/FinalResults/r_all","fig_bird_12tax_accuracy.pdf")
myFunctionMakePlot("final_mosq_12tax_20x/FinalResults/r_all","fig_mosq_12tax_accuracy.pdf")
myFunctionMakePlot("final_mosq_18tax_20x/FinalResults/r_all","fig_mosq_18tax_accuracy.pdf")
myFunctionMakePlot("final_bird_4tax_20x/FinalResults/r_all","fig_bird_4tax_accuracy.pdf")
myFunctionMakePlot("final_butterfly_4tax_20x/FinalResults/r_all","fig_butterfly_4tax_accuracy.pdf")
myFunctionMakePlot("final_mosq_6tax_20x/FinalResults/r_all","fig_mosq_6tax_accuracy.pdf")
myFunctionMakePlot("final_primates_11tax_20x/FinalResults/r_all","fig_primates_11tax_accuracy.pdf")
myFunctionMakePlot("final_humans_4tax_20x/FinalResults/r_all","fig_humans_4tax_accuracy.pdf")


###
# GENERATE ALL TIMING FIGS
###
myFunctionMakeTimingPlot <- function(inputFile,outputFile) {
  
fileToRead = paste(figsPath,inputFile,sep="")
values = read.csv(fileToRead, header=TRUE)

#command line with saving?
myPlot <- ggerrorplot(values, x = "Method", y = "Time", xlab = "Method Used", ylab = "Time (s)",
                      desc_stat = "mean_sd", color = "black", 
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(14,"bold", "black")
)

myPlot <- ggpar(myPlot, x.text.angle = 90) #angles x labels

ggsave(filename=paste(paperPath,outputFile,sep=""), plot=myPlot, width=8,height=8)
}

myFunctionMakeTimingPlot("final_mosq_6tax_20x/FinalTimings/r_all","fig_mosq_6tax_time.pdf")
myFunctionMakeTimingPlot("final_mosq_12tax_20x/FinalTimings/r_all","fig_mosq_12tax_time.pdf")
myFunctionMakeTimingPlot("final_mosq_18tax_20x/FinalTimings/r_all","fig_mosq_18tax_time.pdf")
myFunctionMakeTimingPlot("final_bird_4tax_20x/FinalTimings/r_all","fig_bird_4tax_time.pdf")






###########################
# UNUSED TEST CODE (might use ggplot now instead of ggerrorplot to make it formatted just slightly better and easier)
#########################

#test code
library(tidyverse)
set.seed(0)
n <- 10
cases <- rnorm(n, log2(64), 0.25)
controls <- rnorm(n, log2(64), 0.25)
cases <- 2^(cases)
controls <- 2^(controls)
cases[1:2] <- c(110, 150) #introduce outliers
dat <- data.frame(x = factor(rep(c("Controls", "Cases"), each = n), 
                             levels = c("Controls", "Cases")),
                  Outcome = c(controls, cases))
#One option is simply to show the data points, which you can do like this:
  
dat %>% ggplot(aes(x, Outcome)) + 
geom_jitter(width = 0.05)

myPlot <- ggplot(dat,aes(x, Outcome)) + 
geom_jitter(width = 0.05)

myPlot + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#maybe i just go with this
myPlot + coord_flip()
  
  
####
# TESTING ORDER BY MEANS OF CATEGORIES
###
  
  zz <- "chrom   start
  chr2    39482
  chr1    203918
  chr1    198282
  chrX    7839028
  chr17   3874"
  Data <- read.table(text=zz, header = TRUE)
  
  #library(Hmisc)
  #library(gdata)
  
  Data$chrom  <- reorder.factor(Data$chrom , levels = c("chr1","chr2", "chr17", "chrX"))
  
  Data[order(Data$chrom), ]
  
#### MAKING A FUNCTION
  mysummary <- function(x,y) {
    cat(x)
    cat(y)
  }
  