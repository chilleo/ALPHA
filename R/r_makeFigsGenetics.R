#install.packages("ggpubr")
library(ggpubr)

#demo data sets..work?

figsPath = "/Users/leo/rice/res/data/dgen/tmp/"
paperPath = "/Users/leo/rice/res/papers/RECOMB2016/svn/Genetics19-Leo/draft1/figs/"

#####
# this is what im gonna do
#####

###
# wabi fig 3 updated dgen
###
fileToRead = paste(figsPath,"fig3/fig_fig3.txt",sep="")
values = read.csv(fileToRead, header=TRUE)

#command line with saving?
#plus weird expression stuff to make migration rate not look bolded
myPlot <- ggerrorplot(values, x = "Migration.Rate", y = "p.Value", xlab = expression(paste("Migration ","Rate")), ylab = expression(paste(italic("p"),"-value")),
                      desc_stat = "mean_sd", color = "black", 
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black")
)

ggsave(filename=paste(paperPath,"fig_fig3.pdf",sep=""), plot=myPlot)

###
# wabi fig 3b(new) updated dgen
###
fileToRead = paste(figsPath,"fig3/fig_fig3.txt",sep="")
test_values = read.csv(fileToRead, header=TRUE)
test_valuesBars = aggregate(as.integer(trimws(test_values$Significant)=="TRUE"), by=list(test_values$Migration.Rate), FUN=sum)
test_vals <- test_valuesBars$x
names(test_vals) <- test_valuesBars$Group.1
pdf(paste(paperPath,"fig_fig3b.pdf",sep=""))
par(mar=c(5,6,4,1)+.1) #fix left cutoff
#par(mar=c(7,6,6,1)+.1) #added some bottom and top to fit long y axis title (not working)
barplot(test_vals,xlab = "Migration Rate", ylab="Introgression Detected (# Data Sets)",cex.axis=1.7,cex.names=1.2,cex.lab=2.2)
dev.off()
test_vals # for excel plotting


###
# wabi fig 4b updated dgen
###
fileToRead = paste(figsPath,"fig4/fig_fig4rate.txt",sep="")
values = read.csv(fileToRead, header=TRUE)

myPlot <- ggerrorplot(values, x = "Migration.Rate", y = "p.Value", xlab = "Migration Rate", ylab = "p Value",
                      desc_stat = "mean_sd", color = "black", 
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black")
)
ggsave(filename=paste(paperPath,"fig_fig4b.pdf",sep=""), plot=myPlot,width=6.5)

###
# wabi fig 4c updated dgen
###

fileToRead = paste(figsPath,"fig4/fig_fig4rate.txt",sep="")
test_values = read.csv(fileToRead, header=TRUE)
test_valuesBars = aggregate(as.integer(trimws(test_values$Significant)=="TRUE"), by=list(test_values$Migration.Rate), FUN=sum)
test_vals <- test_valuesBars$x
names(test_vals) <- test_valuesBars$Group.1
pdf(paste(paperPath,"fig_fig4c.pdf",sep=""))
par(mar=c(5,6,4,1)+.1) #fix left cutoff
barplot(test_vals,xlab = "Migration Rate", ylab="Number Of Significant p Values",cex.axis=1.7,cex.names=1.2,cex.lab=2.2)
dev.off()

###
# wabi fig 4e updated dgen
###
fileToRead = paste(figsPath,"fig4/fig_fig4time.txt",sep="")
values = read.csv(fileToRead, header=TRUE)

myPlot <- ggerrorplot(values, x = "Migration.Time", y = "p.Value", xlab = "Migration Time", ylab = "p Value",
                      desc_stat = "mean_sd", color = "black", 
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black")
)
ggsave(filename=paste(paperPath,"fig_fig4e.pdf",sep=""), plot=myPlot)

###
# wabi fig 4f updated dgen
###

fileToRead = paste(figsPath,"fig4/fig_fig4time.txt",sep="")
test_values = read.csv(fileToRead, header=TRUE)
test_valuesBars = aggregate(as.integer(trimws(test_values$Significant)=="TRUE"), by=list(test_values$Migration.Time), FUN=sum)
test_vals <- test_valuesBars$x
names(test_vals) <- test_valuesBars$Group.1
pdf(paste(paperPath,"fig_fig4f.pdf",sep=""))
par(mar=c(5,6,4,1)+.1) #fix left cutoff
barplot(test_vals,xlab = "Migration Time", ylab="Number Of Significant p Values",cex.axis=1.7,cex.names=1.2,cex.lab=2.2)
dev.off()



###
# wabi fig 7b updated dgen
###

fileToRead = paste(figsPath,"fig7/fig_fig7.txt",sep="")
values = read.csv(fileToRead, header=TRUE)

#do a lil renaming
values[,2] = gsub("Dgen", "DGEN", values[,2])
values[,2] = gsub("Dstat", "D-statistic", values[,2])

#ive decided just to show dgen p values since it might be kinda confusing with both (values[2]==" Dgen")
#values = values[1:10,] " Dgen" #only do dgen ones, should always be just 1:10 how i currently have the code working
#values = values[values[2]==" Dgen",] #more precise way to do it just in case

myPlot <- ggerrorplot(values, x = "Statistic", y = "D.or.p.Value", ylab = "Dgen p Value and D statistic D value", xlab = "Statistic Type",
                      desc_stat = "mean_sd", color = "black", 
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black")
)
ggsave(filename=paste(paperPath,"fig_fig7b.pdf",sep=""), plot=myPlot)




###
# wabi fig 7c updated dgen
###

fileToRead = paste(figsPath,"fig7/fig_fig7.txt",sep="")
test_values = read.csv(fileToRead, header=TRUE)

#do a lil renaming
test_values[,2] = gsub("Dgen", "DGEN", test_values[,2])
test_values[,2] = gsub("Dstat", "D-statistic", test_values[,2])

test_valuesBars = aggregate(as.integer(trimws(test_values$Significant)=="TRUE"), by=list(test_values$Statistic), FUN=sum)
test_vals <- test_valuesBars$x
names(test_vals) <- test_valuesBars$Group.1
pdf(paste(paperPath,"fig_fig7c.pdf",sep=""))
par(mar=c(5,6,4,1)+.1) #fix left cutoff
barplot(test_vals,xlab = "Statistic Type", ylab="Number Of Significant p Values",cex.axis=1.7,cex.names=1.2,cex.lab=2.2)
dev.off()
test_vals #to make excel



###
# wabi fig 7 NEW S3 EMBEDDING SCENARIO updated dgen
###
fileToRead = paste(figsPath,"fig7b/fig_fig7b.txt",sep="")
test_values = read.csv(fileToRead, header=TRUE)

#do a lil renaming
test_values[,2] = gsub("Dgen", "DGEN", test_values[,2])
test_values[,2] = gsub("Dstat", "D-statistic", test_values[,2])

test_valuesBars = aggregate(as.integer(trimws(test_values$Significant)=="TRUE"), by=list(test_values$Statistic), FUN=sum)
test_vals <- test_valuesBars$x
names(test_vals) <- test_valuesBars$Group.1
pdf(paste(paperPath,"fig_fig7cNewNetWithS3.pdf",sep=""))
par(mar=c(5,6,4,1)+.1) #fix left cutoff
#barplot(test_vals,xlab = "Statistic Type", ylab="Number Of Significant p Values",cex.axis=1.7,cex.names=1.2,cex.lab=2.2)
barplot(test_vals,xlab = "Statistic Type", ylab="Number Of Significant p Values",cex.axis=1.7,cex.names=1.2,cex.lab=2.2, ylim=c(0,20))
dev.off()
test_vals #to make excel




####
# plot X - 4tax violations (redone cuz of -em stop and end)
####
fileToRead = "/Users/leo/rice/res/ALPHA/wabiFinal-ALPHA-master/CommandLineFiles/fig_4taxViolations.txt"
values = read.csv(fileToRead, header=TRUE)

#added weird expression stuff to have italicized D
myPlot <- ggerrorplot(values, x = "Scenario", y = "D.Value", 
                      desc_stat = "mean_sd", color = "black", 
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black"), ylab = expression(paste(italic("D")," value")), xlab = expression(paste("Scenario",""))
)

ggsave(filename=paste(paperPath,"fig_4taxViolations.pdf",sep=""), plot=myPlot)
#ggsave(filename="fig_4taxViolations.pdf", plot=myPlot)



###################### OLD ####################



####
# plot 1 - 6taxplusone
####
#data("./fig_6taxPlusOne.txt")
#values<- read.table("/Users/leo/rice/res/ALPHA/ALPHA-master103/CommandLineFiles/test_dataPlot.txt")
fileToRead = "fig_6taxPlusOne.txt"
values = read.csv(fileToRead, header=TRUE)

# Add jittered points
#ggerrorplot(values, x = "Migration.Rate", y = "D.Value", 
#            desc_stat = "mean_sd", color = "black",
#            add = "jitter", add.params = list(color = "darkgray")
#)

#command line with saving?
myPlot <- ggerrorplot(values, x = "Migration.Rate", y = "D.Value", 
                      desc_stat = "mean_sd", color = "black", 
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black")
)

ggsave(filename="fig_6taxPlusOne.pdf", plot=myPlot)


####
# plot 2 - 5tax variable mig rate
####
#data("./fig_6taxPlusOne.txt")
#values<- read.table("/Users/leo/rice/res/ALPHA/ALPHA-master103/CommandLineFiles/test_dataPlot.txt")
fileToRead = "fig_5taxVariableM.txt"
values = read.csv(fileToRead, header=TRUE)

# Add jittered points
#ggerrorplot(values, x = "Migration.Rate", y = "D.Value", 
#            desc_stat = "mean_sd", color = "black",
#            add = "jitter", add.params = list(color = "darkgray")
#)

#command line with saving?
myPlot <- ggerrorplot(values, x = "Migration.Rate", y = "D.Value", 
                      desc_stat = "mean_sd", color = "black", 
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black")
)

ggsave(filename="fig_5taxVariableM.pdf", plot=myPlot)


####
# plot 2 - 5tax variable mig rate w inverses
####
#data("./fig_6taxPlusOne.txt")
#values<- read.table("/Users/leo/rice/res/ALPHA/ALPHA-master103/CommandLineFiles/test_dataPlot.txt")
fileToRead = "fig_inv_5taxVariableM.txt"
values = read.csv(fileToRead, header=TRUE)

# Add jittered points
#ggerrorplot(values, x = "Migration.Rate", y = "D.Value", 
#            desc_stat = "mean_sd", color = "black",
#            add = "jitter", add.params = list(color = "darkgray")
#)

#command line with saving?
myPlot <- ggerrorplot(values, x = "Migration.Rate", y = "D.Value", 
                      desc_stat = "mean_sd", color = "black",
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black")
)

ggsave(filename="fig_inv_5taxVariableM.pdf", plot=myPlot)


####
# plot 4 - 6taxplusone w Inverses
####
fileToRead = "fig_inv_6taxPlusOne.txt"
values = read.csv(fileToRead, header=TRUE)

myPlot <- ggerrorplot(values, x = "Migration.Rate", y = "D.Value", 
                      desc_stat = "mean_sd", color = "black", 
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black")
)

ggsave(filename="fig_inv_6taxPlusOne.pdf", plot=myPlot)



####
# plot 5 - 4tax violations
####
fileToRead = "fig_4taxViolations.txt"
values = read.csv(fileToRead, header=TRUE)

myPlot <- ggerrorplot(values, x = "Scenario", y = "D.Value", 
                      desc_stat = "mean_sd", color = "black", 
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black")
)

ggsave(filename="fig_4taxViolations.pdf", plot=myPlot)

####
# plot 6 - 4tax violations inv
####
fileToRead = "fig_inv_4taxViolations.txt"
values = read.csv(fileToRead, header=TRUE)

myPlot <- ggerrorplot(values, x = "Scenario", y = "D.Value", 
                      desc_stat = "mean_sd", color = "black", 
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black")
)

ggsave(filename="fig_inv_4taxViolations.pdf", plot=myPlot)


####
# plot 7 - 6tax violations
####
fileToRead = "fig_6taxViolations.txt"
values = read.csv(fileToRead, header=TRUE)

myPlot <- ggerrorplot(values, x = "Scenario", y = "D.Value", 
                      desc_stat = "mean_sd", color = "black", 
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black")
)

ggsave(filename="fig_6taxViolations.pdf", plot=myPlot)




######
# plot 8a (the new ones luay just asked for based on dfoil fig)
# (migration strength varies)
######
fileToRead = "fig_6taxStrengthTrend.txt"
values = read.csv(fileToRead, header=TRUE)

myPlot <- ggerrorplot(values, x = "Migration.Rate", y = "D.Value", 
                      desc_stat = "mean_sd", color = "black", 
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black")
)

ggsave(filename="fig_6taxStrengthTrend.pdf", plot=myPlot)


######
# plot 8a part 2 (the new ones luay just asked for based on dfoil fig)
# (migration strength varies - how many are significant bar plot)
######

#good for testing
#test_values = read.csv("/Users/leo/rice/res/ALPHA/ALPHA-master104/CommandLineFiles/fig_6taxStrengthTrend.txt", header=TRUE)
#test_valuesBars = aggregate(as.integer(trimws(test_values$Significant)=="TRUE"), by=list(test_values$Migration.Rate), FUN=sum)

#do the sum aggregations (converting the trues and falses with as.numeric)
#turns out there is a leading and trailing space that needs to be trimmed (why trimws appears), so to get to a 1 and 0 vector I have to fully do as.integer(trimws(test_values$Significant)=="TRUE")
fileToRead = "fig_6taxStrengthTrend.txt"
test_values = read.csv(fileToRead, header=TRUE)
test_valuesBars = aggregate(as.integer(trimws(test_values$Significant)=="TRUE"), by=list(test_values$Migration.Rate), FUN=sum)
test_vals <- test_valuesBars$x
names(test_vals) <- test_valuesBars$Group.1
#pdf("/Users/leo/rice/res/ALPHA/ALPHA-master104/CommandLineFiles/test_bars.pdf")

pdf("fig_bars_6taxStrengthTrend.pdf")
par(mar=c(5,6,4,1)+.1) #fix left cutoff
barplot(test_vals,xlab = "Migration Rate", ylab="Number Of Significant D Values",cex.axis=1.7,cex.names=1.2,cex.lab=2.2)
dev.off()

#trying this new way which should work
#test_vals <- 1:5
#names(test_vals) <- LETTERS[1:5]
#barplot(test_vals)
#test save for regular barplot







######
# plot 8b (the new ones luay just asked for based on dfoil fig)
# (migration TIME varies)
######

fileToRead = "fig_6taxTimeTrend.txt"
values = read.csv(fileToRead, header=TRUE)

#myPlot <- ggerrorplot(values, x = "Migration.Time", y = "D.Value", 
#                      desc_stat = "mean_sd", color = "black", title="6 Taxa D Results Varying Mig. Time",
#                      add = "jitter", add.params = list(color = "darkgray")
#)

#strip off title
myPlot <- ggerrorplot(values, x = "Migration.Time", y = "D.Value", 
                      desc_stat = "mean_sd", color = "black", 
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black")
)

ggsave(filename="fig_6taxTimeTrend.pdf", plot=myPlot)


######
# bars part of 8b
######

#do the sum aggregations (converting the trues and falses with as.numeric)
#turns out there is a leading and trailing space that needs to be trimmed (why trimws appears), so to get to a 1 and 0 vector I have to fully do as.integer(trimws(test_values$Significant)=="TRUE")
fileToRead = "fig_6taxTimeTrend.txt"
test_values = read.csv(fileToRead, header=TRUE)
test_valuesBars = aggregate(as.integer(trimws(test_values$Significant)=="TRUE"), by=list(test_values$Migration.Time), FUN=sum)
test_vals <- test_valuesBars$x
names(test_vals) <- test_valuesBars$Group.1

pdf("fig_bars_6taxTimeTrend.pdf") 
#barplot(test_vals,main="Significance With Varying Mig. Times", xlab = "Migration Time", ylab="Number Of Significant D Values") #strip off title
par(mar=c(5,6,4,1)+.1) #fix left cutoff
barplot(test_vals, xlab = "Migration Time", ylab="Number Of Significant D Values",cex.axis=1.7,cex.names=1.7,cex.lab=2.2) #strip off title
dev.off()



####
# plot 9 - 6tax violations full 6 tax dgen
####
fileToRead = "fig_6taxViolationsDgen.txt"
values = read.csv(fileToRead, header=TRUE)

myPlot <- ggerrorplot(values, x = "Scenario", y = "D.Value", 
                      desc_stat = "mean_sd", color = "black", 
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black")
)

ggsave(filename="fig_6taxViolationsDgen.pdf", plot=myPlot)


####
# plot 10 - 6tax violations full 6 tax dgen+dstat (currently just manually combining the two fig files into one)
####
fileToRead = "fig_6taxViolBoth.txt" #combine fig_6taxViolationsDgen and fig_6taxViolations
values = read.csv(fileToRead, header=TRUE)

myPlot <- ggerrorplot(values, x = "Statistic", y = "D.Value", 
                      desc_stat = "mean_sd", color = "black", 
                      add = "jitter", add.params = list(color = "darkgray"),font.x = c(24, "bold", "black"),font.y = c(24, "bold", "black"),font.tickslab = c(18,"bold", "black")
)

ggsave(filename="fig_6taxViolBoth.pdf", plot=myPlot)

