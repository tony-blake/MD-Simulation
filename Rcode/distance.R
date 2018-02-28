library(ggplot2)
library(RColorBrewer)
library(scales)
library(plyr)



etot <- taxaplot4 <- ggplot(bigd3, aes(x = Inf14, y = Inf3, size = length, colour = Taxa)) + 
  scale_x_log10(limits=c(0.001,5000), labels=comma) +
  scale_y_log10(limits=c(0.001,2000), labels=comma) +
  xlab("Coverage (Inf14)") +
  ylab("Coverage (Inf3)") +
  geom_point(alpha=0.1, colour = 'black') +
  geom_point(shape=1) +  
  scale_colour_manual(name="Taxa",values=pcol) +
  scale_size_area(name= "Scaffold length", max_size=20) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 19)))
theme(axis.text=element_text(size=14),
      axis.title=element_text(size=18, face="bold"),legend.title=element_text(size=18,face="bold"))

pdf(file="taxaplot4.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
taxaplot4
dev.off()


summary_etot <- read.table("summary.ETOT", header=F)
colnames(summary_etot) <- c("picoseconds", "Energy")

etot <- ggplot(summary_etot, aes(x = picoseconds, y = Energy, color = gc, size = length)) + 
  scale_x_continuous(limits=c(0,503), labels = comma) +
  scale_y_continuous(limits=c(-180000,-100000), labels = comma) +
  xlab("Time (ps)") +
  ylab("Energy (kcal/mol)") +
  geom_point(alpha = 0.5) +
  scale_size_area(name= "Scaffold length", max_size=20) +
  scale_colour_gradientn(colours=c('blue','green','red')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold"))



pdf(file="etot.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(summary_etot, aes(x = picoseconds, y = Energy, color = gc, size = length)) + 
  scale_x_continuous(limits=c(0,503), labels = comma) +
  scale_y_continuous(limits=c(-180000,-100000), labels = comma) +
  xlab("Time (ps)") +
  ylab("Energy (kcal/mol)") +
  geom_point(alpha = 0.5) +
  scale_size_area(name= "Scaffold length", max_size=20) +
  scale_colour_gradientn(colours=c('blue','green','red')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold"))
dev.off()


plot4 <- ggplot(dcopy, aes(x = , y = Inf3, size = length, colour = Taxa)) + 
  scale_x_log10(limits=c(0.001,5000), labels = comma) +
  scale_y_log10(limits=c(0.001,2000), labels = comma) +
  xlab("Coverage (Inf14)") +
  ylab("Coverage (Inf3)") +
  geom_point(alpha=0.1, colour = 'black') +
  geom_point(shape=1) +  
  scale_colour_manual(name="Taxa",values=pcol) +
  scale_size_area(name= "Scaffold length", max_size=20) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 19))) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold"))



ggplot(d, aes(x = picoseconds, y = Energy, color = gc, size = length)) + 
  scale_x_log10(limits=c(0.001,5000), labels = comma) +
  scale_y_log10(limits=c(0.001,2000), labels = comma) +
  xlab("Coverage (Inf14)") +
  ylab("Coverage (Inf3)") +
  geom_point(alpha = 0.5) +
  scale_size_area(name= "Scaffold length", max_size=20) +
  scale_colour_gradientn(colours=c('blue','green','red')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold"))


library(dplyr)
library(ggplot2)
library(viridis) # for colours
library(ggseas)  # for stat_rollaplyr, rolling average on fly
library(scales)
library(ggrepel)

pdf(file="etot.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(summary_etot,aes(x = picoseconds, y = Energy, color=Energy)) +
  borders('nz', colour = NA, fill = terrain.colors(7)[7]) +
  geom_point(alpha = 0.5) +
  coord_map() +
  scale_size_area("Energy",max_size = 25000) +
  scale_color_viridis("Depth")
dev.off()

ser1 <- read.table("SER236toSER29", header=F)
colnames(ser1) <- c("Picoseconds", "Angstroms")

pdf(file="ser.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(ser1,aes(x = Picoseconds, y = Angstroms)) +
geom_line() +
expand_limits(y=0) +
ggtitle("Change in distance between centre of mass in SER29 in Nisin and centre of mass in SER236 in NSR over 860 picoseconds")
dev.off()

serO <- read.table("SER236OtoSER29O", header=F)
colnames(serO) <- c("Picoseconds", "Angstroms")

pdf(file="serO.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serO,aes(x = Picoseconds, y = Angstroms)) +
geom_line() +
expand_limits(y=0) +
ggtitle("Change in distance between Sidechain Oxygen in SER29 in Nisin and Sidechain Oxygen in SER236 in NSR over 860 picoseconds")
dev.off()


pro1 <- read.table("SER236toPRO29", heade=F)
colnames(pro1) <- c("Picoseconds", "Angstroms")

pdf(file="pro1.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(pro1,aes(x = Picoseconds, y = Angstroms)) +
geom_line() +
expand_limits(y=0) +
ggtitle("Change in distance between center of mass in PRO29 in Nisin and center of mass in SER236 in NSR over 860 picoseconds")
dev.off()



library(dplyr)

serpro <- inner_join(ser1,pro1, by="Picoseconds")

ser1$PRO.blue.SER.yellow <- 1
pro1$PRO.blue.SER.yellow <- 2

serpro1 <- rbind(ser1, pro1)

p <- ggplot(serpro1, aes(x=Picoseconds, y=Angstroms, group=group, col=group, fill=group)) +
  geom_point() +
  geom_smooth(size=1)


pdf(file="serpro6nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro1, aes(x=Picoseconds, y=Angstroms, group=group, col=group, fill=group)) +
geom_smooth(size=1) +
expand_limits(y=0) +
ggtitle("Change in distance between center of mass in PRO29 in Nisin mutant and SER29 in Nisin and center of mass in SER236 in NSR over 860 picoseconds")
dev.off()

pdf(file="serpro6nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro1, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
geom_smooth(size=1) +
expand_limits(y=0) +
scale_colour_gradientn(colours = c('yellow', 'blue')) +
theme(axis.text=element_text(size=14),
    axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
ggtitle("Change in distance between center of mass in Residue 29 in Nisin and the center of mass in SER236 in NSR over 860 picoseconds")
dev.off()



gs.pal <- colorRampPalette(c("#AF1E2D","#0147FA","#FFFF00"),bias=.1,space="rgb")

ser1$Residue <- "SER"
pro1$Residue <- "PRO"

serpro1 <- rbind(ser1, pro1)


pdf(file="serprobright.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serprobright, aes(x=Picoseconds, y=Angstroms, group=group, col=group)) +
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_discrete(h=c(240, 300)) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between center of mass in Residue 29 in Nisinand the center of mass in SER236 in NSR over 860 picoseconds")
dev.off()




pdf(file="serpro6nanos.line.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro1, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_line() +
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between center of mass in Residue 29 in Nisin and the center of mass in SER236 in NSR over 860 picoseconds")
dev.off()

pdf(file="serpro6nano.point.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro1, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_point() +
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between center of mass in Residue 29 in Nisin and the center of mass in SER236 in NSR over 860 picoseconds")
dev.off()


serCA <- read.table("SER236OtoSER29CA", header=F)
colnames(serCA) <- c("Picoseconds", "Angstroms")


proCA <- read.table("SER236OHtoPRO29CA", heade=F)
colnames(proCA) <- c("Picoseconds", "Angstroms")

serCA$PRO.blue.SER.yellow <- 1
proCA$PRO.blue.SER.yellow <- 2

serproCA <- rbind(serCA, proCA)

pdf(file="serproCA6nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproCA, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between alpha-carbon in Residue 29 in Nisin and the Sidechain O in SER236 in NSR over 6 nanoeconds")
dev.off()

pdf(file="serproCA6nanos.line.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproCA, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_line() +
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between alpha-carbon in Residue 29 in Nisin and the Sidechain O in SER236 in NSR over 6 nanoeconds")
dev.off()

pdf(file="serproCA6nano.point.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproCA, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_point() +
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between alpha-carbon in Residue 29 in Nisin and the Sidechain O in SER236 in NSR over 6 nanoeconds")
dev.off()


serC <- read.table("SER236OtoCYS28C", header=F)
colnames(serC) <- c("Picoseconds", "Angstroms")

proC <- read.table("SER236OtoCYS28Cpro", heade=F)
colnames(proC) <- c("Picoseconds", "Angstroms")


serC$PRO.blue.SER.yellow <- 1
proC$PRO.blue.SER.yellow <- 2

serproC <- rbind(serC, proC)

pdf(file="serproC19nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproC, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between cleavage point carbon in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 19 nanoeconds")
dev.off()

pdf(file="serproC19nanos.line.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproC, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_line() +
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between cleavage point carbon in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 19 nanoeconds")
dev.off()

pdf(file="serproC19nano.point.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproC, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_point() +
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between cleavage point carbon in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 19 nanoeconds")
dev.off()


protest <- read.table("SER236OHtoCYS28CpTEST2", heade=F)
colnames(protest) <- c("Picoseconds", "Angstroms")

protest$PRO.blue.SER.yellow <- 2
serprotest <- rbind(serC, protest)


pdf(file="serproCTEST219nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serprotest, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between cleavage point carbon in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 9 nanoeconds")
dev.off()

pdf(file="serproCTEST@9nanos.line.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serprotest, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_line() +
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between cleavage point carbon in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 9 nanoeconds")
dev.off()

pdf(file="serproCTEST29nano.point.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serprotest, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_point() +
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between cleavage point carbon in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 9 nanoeconds")
dev.off()



serCcorrect <- read.table("SER236OtoCYS28Ccorrect", header=F)
colnames(serCcorrect) <- c("Picoseconds", "Angstroms")

proCcorrect <- read.table("SER236OtoCYS28Cpro", heade=F)
colnames(proCcorrect) <- c("Picoseconds", "Angstroms")


serCcorrect$PRO.blue.SER.yellow <- 1
proCcorrect$PRO.blue.SER.yellow <- 2

serproCcorrect <- rbind(serCcorrect, proCcorrect)

pdf(file="serproC25nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproCcorrect, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between cleavage point carbon in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 25 nanoeconds")
dev.off()

pdf(file="serproC25nanos.line.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproCcorrect, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_line() +
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between cleavage point carbon in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 25 nanoeconds")
dev.off()

pdf(file="serproC25nano.point.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproC, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_point() +
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between cleavage point carbon in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 25 nanoeconds")
dev.off()



serCcorrect <- read.table("SER236OtoSER29Ncorrect", header=F)
colnames(serCcorrect) <- c("Picoseconds", "Angstroms")

proNcorrect <- read.table("SER236OtoPRO29Nprocorrect", heade=F)
colnames(proNcorrect) <- c("Picoseconds", "Angstroms")


serNcorrect$PRO.blue.SER.yellow <- 1
proNcorrect$PRO.blue.SER.yellow <- 2

serproNcorrect <- rbind(serNcorrect, proNcorrect)

pdf(file="serproN25nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproNcorrect, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between cleavage point nitrogen in Residue 29 in Nisin and the Sidechain O in SER236 in NSR over 25 nanoeconds")
dev.off()

pdf(file="serproC25nanos.line.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproCcorrect, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_line() +
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between cleavage point carbon in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 25 nanoeconds")
dev.off()

pdf(file="serproC25nano.point.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproC, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_point() +
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between cleavage point carbon in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 25 nanoeconds")
dev.off()



serCnyl <- read.table("SER236OtoCYS28Cnyl", header=F)
colnames(serCnyl) <- c("Picoseconds", "Angstroms")

proCnyl <- read.table("SER236OtoCYS28Cnylpro", heade=F)
colnames(proCnyl) <- c("Picoseconds", "Angstroms")


serCnyl$PRO.blue.SER.yellow <- 1
proCnyl$PRO.blue.SER.yellow <- 2

serproCnyl <- rbind(serCnyl, proCnyl)

pdf(file="serproC40nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproCnyl, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between cleavage point nitrogen in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 40 nanoeconds")
dev.off()

pdf(file="serproC40nanos.line.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproCnyl, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_line() +
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between cleavage point carbon in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 40 nanoeconds")
dev.off()

pdf(file="serproC40nano.point.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproCnyl, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_point() +
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between cleavage point carbon in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 40 nanoeconds")
dev.off()


################################################


serCnyl50 <- read.table("SER236OtoCYS28Cnyl50", header=F)
colnames(serCnyl50) <- c("Picoseconds", "Angstroms")

proCnyl50 <- read.table("SER236OtoCYS28Cnylpro50", heade=F)
colnames(proCnyl50) <- c("Picoseconds", "Angstroms")


serpro50 <- cbind(serCnyl50, proCnyl50)

serpro50 <- serpro50[,-3] 

colnames(serpro50) <- c("Residue", "SERINE", "PROLINE")


dfm50 <- melt(serprobinding50,id.vars = 1)



dfm50$value <- as.numeric(as.character(dfm50$value))

dfm50$Residue <- factor(dfm50$Residue, levels = dfm50$Residue)


serCnyl50$PRO.blue.SER.yellow <- 1
proCnyl50$PRO.blue.SER.yellow <- 2

serproCnyl50 <- rbind(serCnyl50, proCnyl50)


pdf(file="serproC50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproCnyl50, aes(x=Picoseconds, y=Angstroms, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Change in distance between cleavage point nitrogen in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 50 nanoeconds")
dev.off()


############################################ 51 ###################################

serCnyl50 <- read.table("SER236OtoCYS28Cnyl50", header=F)
colnames(serCnyl50) <- c("Picoseconds", "Angstroms")

proCnyl51 <- read.table("SER236OtoCYS28Cnylpro51", heade=F)
colnames(proCnyl50) <- c("Picoseconds", "Angstroms")


serpro51 <- cbind(serCnyl50, proCnyl51)

serpro51 <- serpro51[,-3] 

colnames(serpro51) <- c("Residue", "SERINE", "PROLINE")


dfm51 <- melt(serpro51,id.vars = 1)

colnames(dfm51) <- c("Picoseconds", "Mutation", "Angstroms")

dfm51$value <- as.numeric(as.character(dfm51$value))

dfm50$Picoseconds <- factor(dfm50$Picoseconds, levels = dfm50$Picoseconds)



pdf(file="serproDIST51nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfm51,aes(x = Picoseconds,y = Angstroms, group=Mutation, col=Mutation)) + 
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_discrete(h=c(360, 100)) +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Change in distance between cleavage point nitrogen in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 50 nanoeconds")
dev.off()



############################################ proval50 ###################################

library(reshape)

serCnyl50 <- read.table("SER236OtoCYS28Cnyl50", header=F)
colnames(serCnyl50) <- c("Picoseconds", "Angstroms")

provalCnyl50 <- read.table("SER236OtoCYS28Cnylproval50", heade=F)
colnames(provalCnyl50) <- c("Picoseconds", "Angstroms")


serproval <- cbind(serCnyl50, provalCnyl50)

serproval <- serproval[,-3] 

colnames(serproval) <- c("Residue", "SERINE", "PROLINE")


dfmpv <- melt(serproval,id.vars = 1)

colnames(dfmpv) <- c("Picoseconds", "Mutation", "Angstroms")

dfmpv$value <- as.numeric(as.character(dfmpv$value))

dfm50$Picoseconds <- factor(dfm50$Picoseconds, levels = dfm50$Picoseconds)



pdf(file="serprovalDIST50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmpv,aes(x = Picoseconds,y = Angstroms, group=Mutation, col=Mutation)) + 
  geom_smooth(size=1) +
  geom_line() +
  expand_limits(y=0) +
  scale_colour_discrete(h=c(360, 100)) +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Change in distance between cleavage point nitrogen in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 50 nanoeconds")
dev.off()



pdf(file="serprovalDIST50nanoslineText.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmpv,aes(x = Picoseconds,y = Angstroms, group=Mutation, col=Mutation)) + 
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_discrete(h=c(360, 100)) +
  theme(axis.text=element_text(angle = 90,hjust = 1, ,vjust = 0, size=40),
        axis.title=element_text(size=30,face="bold"),legend.text=element_text(size = 20), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Change in distance between cleavage point nitrogen in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 50 nanoeconds")
dev.off()


############################################ alanine 50 ###################################

library(reshape)

serCnyl50 <- read.table("SER236OtoCYS28Cnyl50", header=F)
colnames(serCnyl50) <- c("Picoseconds", "Angstroms")

alaCnyl50 <- read.table("SER236OtoCYS28Cnylala50", heade=F)
colnames(alaCnyl50) <- c("Picoseconds", "Angstroms")


ala <- cbind(serCnyl50, alaCnyl50)

ala <- ala[,-3] 

colnames(ala) <- c("Residue", "SERINE", "ALALINE")


dfmala <- melt(ala,id.vars = 1)

colnames(dfmala) <- c("Picoseconds", "Mutation", "Angstroms")

dfmala$value <- as.numeric(as.character(dfmala$value))

dfm50$Picoseconds <- factor(dfm50$Picoseconds, levels = dfm50$Picoseconds)



pdf(file="alaDIST50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmala,aes(x = Picoseconds,y = Angstroms, group=Mutation, col=Mutation)) + 
  geom_smooth(size=1) +
  geom_line() +
  expand_limits(y=0) +
  scale_colour_discrete(h=c(360, 100)) +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Change in distance between cleavage point nitrogen in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 50 nanoeconds")
dev.off()



pdf(file="alaDIST50nanosline.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmala,aes(x = Picoseconds,y = Angstroms, group=Mutation, col=Mutation)) + 
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_discrete(h=c(360, 100)) +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Change in distance between cleavage point nitrogen in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 50 nanoeconds")
dev.off()

############################################ alanine and proval and ser 50 ###################################

library(reshape)

serCnyl50 <- read.table("SER236OtoCYS28Cnyl50", header=F)
colnames(serCnyl50) <- c("Picoseconds", "Angstroms")

alaCnyl50 <- read.table("SER236OtoCYS28Cnylala50", heade=F)
colnames(alaCnyl50) <- c("Picoseconds", "Angstroms")

provalCnyl50 <- read.table("SER236OtoCYS28Cnylproval50", heade=F)
colnames(provalCnyl50) <- c("Picoseconds", "Angstroms")


alaps <- cbind(serCnyl50, provalCnyl50, alaCnyl50)

alaps <- alaps[,-c(3,5)] 

colnames(alaps) <- c("Residue", "SERINE", "PROLINE","ALALINE")


dfmalaps <- melt(alaps,id.vars = 1)

colnames(dfmalaps) <- c("Picoseconds", "Mutation", "Angstroms")

dfmalaps$value <- as.numeric(as.character(dfmalaps$value))

dfm50$Picoseconds <- factor(dfm50$Picoseconds, levels = dfm50$Picoseconds)







pdf(file="alaproserDIST50nanosline2.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmalaps,aes(x = Picoseconds,y = Angstroms, group=Mutation, col=Mutation)) + 
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_discrete(h=c(360, 100, 40)) + 
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Change in distance between cleavage point nitrogen in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 50 nanoeconds")
dev.off()






dfmproalaBE <- melt(proalabinding50,id.vars = 1)



dfmproalaBE$value <- as.numeric(as.character(dfmproalaBE$value))

dfmproalaBE$Residue <- factor(dfmproalaBE$Residue, levels = dfmproalaBE$Residue)

pdf(file="proalaBindE50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmproalaBE,aes(x = Residue,y = value)) + 
  geom_bar(aes(fill = variable),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Binding Energy of Residues")
dev.off()

########################## pico to nano #################################


dfmpvnano <- dfmpv[,1]/1000

dfmpvnano <- cbind(dfmpvnano,dfmpv)

dfmpvnano <- data.frame(dfmpvnano)

dfmpvnano <- dfmpvnano[,-2]

colnames(dfmpvnano) <- c("Nanoseconds", "Mutation", "Angstroms")


#Change in distance between cleavage point nitrogen in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 50 nanoeconds
pdf(file="serprovalDIST50nanoslineText.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmpvnano,aes(x = Nanoseconds,y = Angstroms, group=Mutation, col=Mutation)) + 
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_discrete(h=c(360, 100)) +
  theme(axis.text=element_text(angle = 0,hjust = 1, size=40),
        axis.title=element_text(size=30,face="bold"),legend.text=element_text(size = 40), legend.title=element_text(size=40,face="bold")) +
  ggtitle("a")
dev.off()


########################## pico to nano #################################


colnames(serprovalbinding50) <- c("Residue", "Nisin A", "Nisin PV")



colnames(serproval) <- c("Residue", "Nisin A", "Nisin PV")


dfmpv <- melt(serproval,id.vars = 1)

dfmpv <- data.frame(dfmpv)





dfmpvnano <- dfmpv[,1]/1000

dfmpvnano <- data.frame(dfmpvnano)

dfmpvnano <- cbind(dfmpvnano,dfmpv)

dfmpvnano <- dfmpvnano[,-2]

colnames(dfmpvnano) <- c("Nanoseconds", "Mutation", "Angstroms")

main4=expression(Distance[ring(A)])

#RMSD between cleavage point nitrogen in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 50 nanoeconds
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

