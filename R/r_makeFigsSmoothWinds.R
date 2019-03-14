library(ggpubr)

figsPath = "/Users/leo/rice/tmp/"
paperPath = "/Users/leo/rice/res/papers/RECOMB2016/svn/ISBRA19-Leo/draft1/figs/"



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

ggsave(filename=paste(paperPath,"fig_mosq6tax_time.pdf",sep=""), plot=myPlot, width=8,height=8)



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
  