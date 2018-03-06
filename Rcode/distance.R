library(ggplot2)
library(RColorBrewer)
library(scales)
library(plyr)
library(reshape)


#Read in raw data
serCnyl50 <- read.table("SER236OtoCYS28Cnyl50", header=F)
colnames(serCnyl50) <- c("Picoseconds", "Angstroms")
provalCnyl50 <- read.table("SER236OtoCYS28Cnylproval50", heade=F)
colnames(provalCnyl50) <- c("Picoseconds", "Angstroms")

#Combine Nisin A and Nisin PV dataframe Columns
serproval <- cbind(serCnyl50, provalCnyl50)
serproval <- serproval[,-3] 
colnames(serproval) <- c("Residue", "Nisin A", "Nisin PV")

#Change dataframe structure for input to ggplot2
dfmpv <- melt(serproval,id.vars = 1)
dfmpv <- data.frame(dfmpv)
dfmpvnano <- dfmpv[,1]/1000
dfmpvnano <- data.frame(dfmpvnano)
dfmpvnano <- cbind(dfmpvnano,dfmpv)
dfmpvnano <- dfmpvnano[,-2]
colnames(dfmpvnano) <- c("Nanoseconds", "Mutation", "Angstroms")


#Create variable for ggplot2 axis text
main4=expression(Distance[ring(A)])

#Create plot of distance between cleavage point nitrogen in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 50 nanoeconds
pdf(file="serprovalDIST50nanoslineText.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmpvnano,aes(x = Nanoseconds,y = Angstroms, group=Mutation, col=Mutation)) + 
  labs(y=main4) +
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_discrete(h=c(360, 100)) +
  theme(axis.text=element_text(angle = 0,hjust = 1, size=40),
        axis.title=element_text(size=40),legend.text=element_text(size = 40), legend.title=element_text(size=40,face="bold")) +
  ggtitle("")
dev.off()

