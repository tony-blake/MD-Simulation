rmsd <- read.xlsx("perresavg.BB.xlsx", sheetIndex = 1)
rmsd <- rmsd[,-c(4,5)]
rmsd.sub <- rmsd[c(1:10,30:40,204:210,288:300),]

rmsdpro <- read.xlsx("perresavg.BB.pro.xls", sheetIndex = 1, header=F)

rmsdpro <- rmsdpro[,-c(1,3)]
rmsdpro <- rmsdpro[-1,]

rmsdpro.sub <- rmsdpro[c(1:10,30:40,204:210,288:300),]

colnames(rmsd.sub) <- c("Residue", "RMSD", "std.dev")
colnames(rmsdpro.sub) <- c("Residue", "RMSD", "std.dev")


serproRMSD <- cbind(rmsd.sub, rmsdpro.sub)

serproRMSD <- serproRMSD[,-c(3,4,6)] 

colnames(serproRMSD) <- c("Residue", "SERINE.RMSD", "PROLINE.RMSD")





dfmRMSD <- melt(serproRMSD,id.vars = 1)



dfmRMSD$value <- as.numeric(as.character(dfmRMSD$value))



pdf(file="serproRMSD35nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmRMSD,aes(x = Residue,y = value)) + 
  geom_bar(aes(fill = variable),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(size=4),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("RMSD of Residues")
dev.off()


######################################################################################################

colnames(rmsd) <- c("Residue", "RMSD", "std.dev")
colnames(rmsdpro) <- c("Residue", "RMSD", "std.dev")

serpro.All.RMSD <- cbind(rmsd, rmsdpro)

serpro.All.RMSD <- serpro.All.RMSD[,-c(3,4,6)] 

colnames(serpro.All.RMSD) <- c("Residue", "SERINE.RMSD", "PROLINE.RMSD")


dfm.All.RMSD <- melt(serpro.All.RMSD,id.vars = 1)



dfm.All.RMSD$value <- as.numeric(as.character(dfm.All.RMSD$value))


pdf(file="serproAllRMSD35nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfm.All.RMSD,aes(x = Residue,y = value)) + 
  geom_bar(aes(fill = variable),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(size=4),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("RMSD of Residues")
dev.off()



######################################################################################################
rmsd.nisin <- rmsd[288:300,]
rmsdpro.nisin <- rmsdpro[288:300,]

colnames(rmsd) <- c("Residue", "RMSD", "std.dev")
colnames(rmsdpro) <- c("Residue", "RMSD", "std.dev")

serpro.nisin.RMSD <- cbind(rmsd.nisin, rmsdpro.nisin)

serpro.nisin.RMSD <- serpro.nisin.RMSD[,-c(3,4,6)] 

colnames(serpro.nisin.RMSD) <- c("Residue", "SERINE.RMSD", "PROLINE.RMSD")


dfm.nisin.RMSD <- melt(serpro.nisin.RMSD,id.vars = 1)



dfm.nisin.RMSD$value <- as.numeric(as.character(dfm.nisin.RMSD$value))

dummy1 <- serprobinding[14:26,]

dfm.nisin.RMSD$RESIDUE <- dummy1$Residue
dfm.nisin.RMSD$RESIDUE <- factor(dfm.nisin.RMSD$RESIDUE, levels = dfm.nisin.RMSD$RESIDUE)

pdf(file="serpronisinRMSD35nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfm.nisin.RMSD,aes(x = RESIDUE,y = value)) + 
  geom_bar(aes(fill = variable),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("RMSD of Residues")
dev.off()

#######################################


rmsd50ser <- read.delim("perresavg.BB.50.dat", sep="", stringsAsFactors=FALSE, header = TRUE)
rmsd50pro <- read.delim("perresavg.BB.50pro.dat", sep="", stringsAsFactors=FALSE, header = TRUE)

rmsd.nisin50ser <- rmsd50ser[288:300,]
rmsdpro.nisin50pro <- rmsd50pro[288:300,]

colnames(rmsd.nisin50ser) <- c("Residue", "RMSD", "std.dev")
colnames(rmsdpro.nisin50pro) <- c("Residue", "RMSD", "std.dev")

serpro.nisin.RMSD50 <- cbind(rmsd.nisin50ser, rmsdpro.nisin50pro)

serpro.nisin.RMSD50 <- serpro.nisin.RMSD50[,-c(3,4,6)] 

colnames(serpro.nisin.RMSD50) <- c("Residue", "SERINE.RMSD", "PROLINE.RMSD")


dfm.nisin.RMSD50 <- melt(serpro.nisin.RMSD50,id.vars = 1)



dfm.nisin.RMSD50$value <- as.numeric(as.character(dfm.nisin.RMSD$value))

dummy1 <- serprobinding[14:26,]

dfm.nisin.RMSD50$RESIDUE <- dummy1$Residue
dfm.nisin.RMSD50$RESIDUE <- factor(dfm.nisin.RMSD50$RESIDUE, levels = dfm.nisin.RMSD50$RESIDUE)

pdf(file="serpronisinRMSD50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfm.nisin.RMSD50,aes(x = RESIDUE,y = value)) + 
  geom_bar(aes(fill = variable),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("RMSD of Residues after 50 nanosecond")
dev.off()

############################################### 51 #######################################################

rmsd50ser <- read.delim("perresavg.BB.50.dat", sep="", stringsAsFactors=FALSE, header = TRUE)
rmsd51pro <- read.delim("perresavg.BB.51pro.dat", sep="", stringsAsFactors=FALSE, header = TRUE)

rmsd.nisin50ser <- rmsd50ser[288:300,]
rmsdpro.nisin51pro <- rmsd51pro[288:300,]

colnames(rmsd.nisin50ser) <- c("Residue", "RMSD", "std.dev")
colnames(rmsdpro.nisin51pro) <- c("Residue", "RMSD", "std.dev")

serpro.nisin.RMSD51 <- cbind(rmsd.nisin50ser, rmsdpro.nisin51pro)

serpro.nisin.RMSD51 <- serpro.nisin.RMSD51[,-c(3,4,6)] 

colnames(serpro.nisin.RMSD51) <- c("Residue", "SERINE.RMSD", "PROLINE.RMSD")


dfm.nisin.RMSD51 <- melt(serpro.nisin.RMSD51,id.vars = 1)



dfm.nisin.RMSD51$value <- as.numeric(as.character(dfm.nisin.RMSD51$value))

dummy1 <- serprobinding[14:26,]

dfm.nisin.RMSD51$RESIDUE <- dummy1$Residue
dfm.nisin.RMSD51$RESIDUE <- factor(dfm.nisin.RMSD51$RESIDUE, levels = dfm.nisin.RMSD51$RESIDUE)

pdf(file="serpronisinRMSD51nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfm.nisin.RMSD51,aes(x = RESIDUE,y = value)) + 
  geom_bar(aes(fill = variable),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("RMSD of Residues after 50 nanosecond")
dev.off()


################# rmsf ##################################

############################################## 51 #######################################################


RMSF51pro <- read.delim("RMSF50nisin", sep="", stringsAsFactors=FALSE, header = TRUE)


serpro.nisin.RMSF51 <- cbind(RMSF50ser, RMSF51pro)

serpro.nisin.RMSF51 <- serpro.nisin.RMSF51[,-c(3)] 

colnames(serpro.nisin.RMSF51) <- c("Residue", "SERINE.RMSF", "PROLINE.RMSF")
RMSF50ser <- read.delim("RMSF50nisinser", sep="", stringsAsFactors=FALSE, header = TRUE)

dfm.nisin.RMSF51 <- melt(serpro.nisin.RMSF51,id.vars = 1)



dfm.nisin.RMSF51$value <- as.numeric(as.character(dfm.nisin.RMSF51$value))

dummy1 <- serprobinding[14:26,]

dfm.nisin.RMSF51$RESIDUE <- dummy1$Residue
dfm.nisin.RMSF51$RESIDUE <- factor(dfm.nisin.RMSF51$RESIDUE, levels = dfm.nisin.RMSF51$RESIDUE)

pdf(file="serpronisinRMSF51nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfm.nisin.RMSF51,aes(x = RESIDUE,y = value)) + 
  geom_bar(aes(fill = variable),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("RMSF of Residues after 50 nanosecond")
dev.off()


############################################## proval #######################################################


RMSF50proval <- read.delim("RMSF50nisinproval", sep="", stringsAsFactors=FALSE, header = TRUE)


serproval.nisin.RMSF51 <- cbind(RMSF50ser, RMSF50proval)

serproval.nisin.RMSF51 <- serproval.nisin.RMSF51[,-c(3)] 

colnames(serproval.nisin.RMSF51) <- c("Residue", "SERINE.RMSF", "PROLINE.RMSF")
RMSF50ser <- read.delim("RMSF50nisinser", sep="", stringsAsFactors=FALSE, header = TRUE)

dfmpv.nisin.RMSF51 <- melt(serproval.nisin.RMSF51,id.vars = 1)



dfmpv.nisin.RMSF51$value <- as.numeric(as.character(dfmpv.nisin.RMSF51$value))

dummy1 <- serprovalbinding50[14:26,]

dfmpv.nisin.RMSF51$RESIDUE <- dummy1$Residue
dfmpv.nisin.RMSF51$RESIDUE <- factor(dfmpv.nisin.RMSF51$RESIDUE, levels = dfmpv.nisin.RMSF51$RESIDUE)

pdf(file="serprovalnisinRMSF5onanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmpv.nisin.RMSF51,aes(x = RESIDUE,y = value)) + 
  geom_bar(aes(fill = variable),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("RMSF of Residues after 50 nanosecond")
dev.off()

################################# rmsd proval #######################################################3#####

rmsd50ser <- read.delim("perresavg.BB.50.dat", sep="", stringsAsFactors=FALSE, header = TRUE)
rmsd51proval <- read.delim("perresavg.BB.50proval.dat", sep="", stringsAsFactors=FALSE, header = TRUE)

rmsd.nisin50ser <- rmsd50ser[288:300,]
rmsdproval.nisin51proval <- rmsd51proval[288:300,]

colnames(rmsd.nisin50ser) <- c("Residue", "RMSD", "std.dev")
colnames(rmsdproval.nisin51proval) <- c("Residue", "RMSD", "std.dev")

serproval.nisin.RMSD51 <- cbind(rmsd.nisin50ser, rmsdproval.nisin51proval)

serproval.nisin.RMSD51 <- serproval.nisin.RMSD51[,-c(3,4,6)] 

colnames(serproval.nisin.RMSD51) <- c("Residue", "SERINE.RMSD", "PROLINE.RMSD")

library(reshape)
dfmproval.nisin.RMSD51 <- melt(serproval.nisin.RMSD51,id.vars = 1)



dfmproval.nisin.RMSD51$value <- as.numeric(as.character(dfmproval.nisin.RMSD51$value))

dummy1 <- serprobinding[14:26,]

dfmproval.nisin.RMSD51$RESIDUE <- dummy1$Residue
dfmproval.nisin.RMSD51$RESIDUE <- factor(dfmproval.nisin.RMSD51$RESIDUE, levels = dfmproval.nisin.RMSD51$RESIDUE)

pdf(file="serprovalnisinRMSD50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmproval.nisin.RMSD51,aes(x = RESIDUE,y = value)) + 
  geom_bar(aes(fill = variable),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("RMSD of Residues after 50 nanosecond")
dev.off()

########################### RMSD proval over time  #######################################

rms_Time50ser <- read.delim("rms_vs_time.BB.50.dat", sep="", stringsAsFactors=FALSE, header = TRUE)
rms_Time50proval <- read.delim("rms_vs_time.BB.50proval.dat", sep="", stringsAsFactors=FALSE, header = TRUE)

rms_nisinser <- read.delim("rmsdnisinovertime", sep="", stringsAsFactors=FALSE, header = TRUE)
rms_nisinproval <- read.delim("rmsdpvovertime", sep="", stringsAsFactors=FALSE, header = TRUE)


library(reshape)
library(ggplot2)

colnames(rms_nisinser) <- c("Picoseconds", "Angstroms")
colnames(rms_nisinproval) <- c("Picoseconds", "Angstroms")




rms.serproval <- cbind(rms_nisinser, rms_nisinproval)

rms.serproval <- rms.serproval[,-3] 

colnames(rms.serproval) <- c("Residue", "SERINE", "PROLINE")


dfmrmspv <- melt(rms.serproval,id.vars = 1)

colnames(dfmrmspv) <- c("Picoseconds", "Mutation", "Angstroms")

dfmrmspv$value <- as.numeric(as.character(dfmrmspv$value))

dfm50$Picoseconds <- factor(dfm50$Picoseconds, levels = dfm50$Picoseconds)

library(ggplot2)

pdf(file="serprovalRMSTIME50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmrmspv,aes(x = Picoseconds,y = Angstroms, group=Mutation, col=Mutation)) + 
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_discrete(h=c(360, 100)) +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("RMS deviation of nisin from 1st picosecond over 50 nanoeconds")
dev.off()

########################### RMSD pro over time  #######################################

rms_Time50ser <- read.delim("rms_vs_time.BB.50.dat", sep="", stringsAsFactors=FALSE, header = TRUE)
rms_Time50proval <- read.delim("rms_vs_time.BB.50proval.dat", sep="", stringsAsFactors=FALSE, header = TRUE)

rms_nisinser <- read.delim("rmsdnisinovertime", sep="", stringsAsFactors=FALSE, header = TRUE)
rms_nisinpro <- read.delim("rmsdproovertime", sep="", stringsAsFactors=FALSE, header = TRUE)


library(reshape)
library(ggplot2)

colnames(rms_nisinser) <- c("Picoseconds", "Angstroms")
colnames(rms_nisinpro) <- c("Picoseconds", "Angstroms")




rms.serpro <- cbind(rms_nisinser, rms_nisinpro)

rms.serpro <- rms.serpro[,-3] 

colnames(rms.serpro) <- c("Residue", "SERINE", "PROLINE")


dfmrmspro <- melt(rms.serpro,id.vars = 1)

colnames(dfmrmspro) <- c("Picoseconds", "Mutation", "Angstroms")

dfmrmspv$value <- as.numeric(as.character(dfmrmspv$value))

dfm50$Picoseconds <- factor(dfm50$Picoseconds, levels = dfm50$Picoseconds)



pdf(file="serproRMSTIME50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmrmspro,aes(x = Picoseconds,y = Angstroms, group=Mutation, col=Mutation)) + 
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_discrete(h=c(360, 100)) +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("RMS deviation of nisin from 1st picosecond over 50 nanoeconds")
dev.off()




########################### RMSD proval over time with lines
rms_Time50ser <- read.delim("rms_vs_time.BB.50.dat", sep="", stringsAsFactors=FALSE, header = TRUE)
rms_Time50proval <- read.delim("rms_vs_time.BB.50proval.dat", sep="", stringsAsFactors=FALSE, header = TRUE)

rms_nisinser <- read.delim("rmsdnisinovertime", sep="", stringsAsFactors=FALSE, header = TRUE)
rms_nisinproval <- read.delim("rmsdpvovertime", sep="", stringsAsFactors=FALSE, header = TRUE)


library(reshape)
library(ggplot2)

colnames(rms_nisinser) <- c("Picoseconds", "Angstroms")
colnames(rms_nisinproval) <- c("Picoseconds", "Angstroms")




rms.serproval <- cbind(rms_nisinser, rms_nisinproval)

rms.serproval <- rms.serproval[,-3] 

colnames(rms.serproval) <- c("Residue", "SERINE", "PROLINE")


dfmrmspv <- melt(rms.serproval,id.vars = 1)

colnames(dfmrmspv) <- c("Picoseconds", "Mutation", "Angstroms")

dfmrmspv$value <- as.numeric(as.character(dfmrmspv$value))

dfm50$Picoseconds <- factor(dfm50$Picoseconds, levels = dfm50$Picoseconds)



pdf(file="serprovalRMSTIME50nanosline.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmrmspv,aes(x = Picoseconds,y = Angstroms, group=Mutation, col=Mutation)) + 
  geom_line() +
  expand_limits(y=0) +
  scale_colour_discrete(h=c(360, 100)) +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("RMS deviation of nisin from 1st picosecond over 50 nanoeconds")
dev.off()


########################### RMSD alanine over time  #######################################

rms_Time50ser <- read.delim("rms_vs_time.BB.50.dat", sep="", stringsAsFactors=FALSE, header = TRUE)
rms_Time50ala <- read.delim("rms_vs_time.BB.50ala.dat", sep="", stringsAsFactors=FALSE, header = TRUE)

rms_nisinser <- read.delim("rmsdnisinovertime", sep="", stringsAsFactors=FALSE, header = TRUE)
rms_nisinala <- read.delim("rmsdalaovertime", sep="", stringsAsFactors=FALSE, header = TRUE)


library(reshape)
library(ggplot2)

colnames(rms_nisinser) <- c("Picoseconds", "Angstroms")
colnames(rms_nisinala) <- c("Picoseconds", "Angstroms")




rms.ala <- cbind(rms_nisinser, rms_nisinala)

rms.ala <- rms.ala[,-3] 

colnames(rms.ala) <- c("Residue", "SERINE", "ALANINE")


dfmrmsala <- melt(rms.ala,id.vars = 1)

colnames(dfmrmsala) <- c("Picoseconds", "Mutation", "Angstroms")

dfmrmsala$value <- as.numeric(as.character(dfmrmsala$value))

dfm50$Picoseconds <- factor(dfm50$Picoseconds, levels = dfm50$Picoseconds)



pdf(file="alaRMSTIME50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmrmsala,aes(x = Picoseconds,y = Angstroms, group=Mutation, col=Mutation)) + 
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_discrete(h=c(360, 100)) +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("RMS deviation of nisin from 1st picosecond over 50 nanoeconds")
dev.off()
############################################## rmsf alanine #######################################################


RMSF50ala <- read.delim("RMSF50nisinala", sep="", stringsAsFactors=FALSE, header = TRUE)


ala.nisin.RMSF51 <- cbind(RMSF50ser, RMSF50ala)

ala.nisin.RMSF51 <- ala.nisin.RMSF51[,-c(3)] 

colnames(ala.nisin.RMSF51) <- c("Residue", "SERINE.RMSF", "ALANINE.RMSF")
RMSF50ser <- read.delim("RMSF50nisinser", sep="", stringsAsFactors=FALSE, header = TRUE)

dfmala.nisin.RMSF51 <- melt(ala.nisin.RMSF51,id.vars = 1)



dfmala.nisin.RMSF51$value <- as.numeric(as.character(dfmala.nisin.RMSF51$value))

dummy1 <- serprovalbinding50[14:26,]

dfmala.nisin.RMSF51$RESIDUE <- dummy1$Residue
dfmala.nisin.RMSF51$RESIDUE <- factor(dfmala.nisin.RMSF51$RESIDUE, levels = dfmala.nisin.RMSF51$RESIDUE)

pdf(file="alanisinRMSF5onanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmala.nisin.RMSF51,aes(x = RESIDUE,y = value)) + 
  geom_bar(aes(fill = variable),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("RMSF of Residues after 50 nanosecond")
dev.off()

############################################## rmsf alanine and proline #######################################################


RMSF50ala <- read.delim("RMSF50nisinala", sep="", stringsAsFactors=FALSE, header = TRUE)


proala.nisin.RMSF51 <- cbind(RMSF50ser, RMSF50ala, RMSF50proval)

proala.nisin.RMSF51 <- proala.nisin.RMSF51[,-c(3,5)] 

colnames(proala.nisin.RMSF51) <- c("Residue", "SERINE.RMSF", "ALANINE.RMSF", "PROLINE.RMSF")
RMSF50ser <- read.delim("RMSF50nisinser", sep="", stringsAsFactors=FALSE, header = TRUE)

dfmproala.nisin.RMSF51 <- melt(proala.nisin.RMSF51,id.vars = 1)



dfmproala.nisin.RMSF51$value <- as.numeric(as.character(dfmproala.nisin.RMSF51$value))

dummy1 <- serprovalbinding50[14:26,]

dfmproala.nisin.RMSF51$RESIDUE <- dummy1$Residue
dfmproala.nisin.RMSF51$RESIDUE <- factor(dfmproala.nisin.RMSF51$RESIDUE, levels = dfmproala.nisin.RMSF51$RESIDUE)

pdf(file="proalanisinRMSF5onanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmproala.nisin.RMSF51,aes(x = RESIDUE,y = value)) + 
  geom_bar(aes(fill = variable),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("RMSF of Residues after 50 nanosecond")
dev.off()


########################### RMSD alanine over time  #######################################

rms_Time50ser <- read.delim("rms_vs_time.BB.50.dat", sep="", stringsAsFactors=FALSE, header = TRUE)
rms_Time50ala <- read.delim("rms_vs_time.BB.50ala.dat", sep="", stringsAsFactors=FALSE, header = TRUE)

rms_nisinser <- read.delim("rmsdnisinovertime", sep="", stringsAsFactors=FALSE, header = TRUE)
rms_nisinala <- read.delim("rmsdalaovertime", sep="", stringsAsFactors=FALSE, header = TRUE)
rms_nisinproval <- read.delim("rmsdpvovertime", sep="", stringsAsFactors=FALSE, header = TRUE)

library(reshape)
library(ggplot2)

colnames(rms_nisinser) <- c("Picoseconds", "Angstroms")
colnames(rms_nisinala) <- c("Picoseconds", "Angstroms")
colnames(rms_nisinproval) <- c("Picoseconds", "Angstroms")




rms.alaps <- cbind(rms_nisinser, rms_nisinproval, rms_nisinala)

rms.alaps <- rms.alaps[,-c(3,5)] 

colnames(rms.alaps) <- c("Residue", "SERINE", "PROLINE", "ALANINE")


dfmrmsalaps <- melt(rms.alaps,id.vars = 1)

colnames(dfmrmsalaps) <- c("Picoseconds", "Mutation", "Angstroms")

dfmrmsalaps$value <- as.numeric(as.character(dfmrmsalaps$value))

dfm50$Picoseconds <- factor(dfm50$Picoseconds, levels = dfm50$Picoseconds)



pdf(file="alaproserRMSTIME50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmrmsalaps,aes(x = Picoseconds,y = Angstroms, group=Mutation, col=Mutation)) + 
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_discrete(h=c(360, 100, 40)) +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("RMS deviation of nisin from 1st picosecond over 50 nanoeconds")
dev.off()



########################## pico to nano #################################


colnames(serprovalbinding50) <- c("Residue", "Nisin A", "Nisin PV")



colnames(rms.serproval) <- c("Residue", "Nisin A", "Nisin PV")








dfmrmspvnano <- dfmrmspv[,1]/1000

dfmrmspvnano <- data.frame(dfmrmspvnano)
dfmrmspvnano <- cbind(dfmpvnano,dfmrmspv)



dfmrmspvnano <- dfmrmspvnano[,-c(2,3,4)]

colnames(dfmrmspvnano) <- c("Nanoseconds", "Mutation", "Angstroms")

 

#RMSF between cleavage point nitrogen in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 50 nanoeconds
pdf(file="serprovalRMSD50nanoslineText.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmrmspvnano,aes(x = Nanoseconds,y = Angstroms, group=Mutation, col=Mutation)) + 
  labs(y=main4) +
  geom_smooth(size=1) +
  expand_limits(y=0) +
  scale_colour_discrete(h=c(360, 100)) +
  theme(axis.text=element_text(angle = 0,hjust = 1, size=40),
        axis.title=element_text(size=40),legend.text=element_text(size = 40), legend.title=element_text(size=40,face="bold")) +
  ggtitle("")
dev.off()



main4=expression(RMSD[ring(A)])

main=expression('Lys'[22], 'Dbb'[23], 'Ala'[24], 'Dbb'[25], 'Cys'[26], 'Asn'[27], 'Cys'[28], 'Ser Pro'[29], 'ile Val'[30], 'His'[31], 'Val'[32], 'Dha'[33], 'Lys'[34])

main2=expression('Thr'[230], 'Asn'[231], 'Hie'[232], 'Lys'[233], 'Thr'[234], 'Ala'[235], 'Ser'[236], 'Ser'[237], 'Ala'[238], 'Glu'[239], 'Met'[240], 'Thr'[241], 'Phe'[242])


main5=expression('Thr'[230], 'Asn'[231], 'Hie'[232], 'Lys'[233], 'Thr'[234], 'Ala'[235], 'Ser'[236], 'Ser'[237], 'Ala'[238], 'Glu'[239], 'Met'[240], 'Thr'[241], 'Phe'[242], 'Lys'[22], 'Dbb'[23], 'Ala'[24], 'Dbb'[25], 'Cys'[26], 'Asn'[27], 'Cys'[28], 'Ser Pro'[29], 'ile Val'[30], 'His'[31], 'Val'[32], 'Dha'[33], 'Lys'[34])


################################# rmsd proval Text #######################################################3#####

rmsd50ser <- read.delim("perresavg.BB.50.dat", sep="", stringsAsFactors=FALSE, header = TRUE)
rmsd51proval <- read.delim("perresavg.BB.50proval.dat", sep="", stringsAsFactors=FALSE, header = TRUE)

rmsd.nisin50ser <- rmsd50ser[288:300,]
rmsdproval.nisin51proval <- rmsd51proval[288:300,]

colnames(rmsd.nisin50ser) <- c("Residue", "RMSD", "std.dev")
colnames(rmsdproval.nisin51proval) <- c("Residue", "RMSD", "std.dev")

serproval.nisin.RMSD51 <- cbind(rmsd.nisin50ser, rmsdproval.nisin51proval)

serproval.nisin.RMSD51 <- serproval.nisin.RMSD51[,-c(3,4,6)] 

colnames(serproval.nisin.RMSD51) <- c("Residue", "Nisin A", "Nisin PV")

library(reshape)
dfmproval.nisin.RMSD51 <- melt(serproval.nisin.RMSD51,id.vars = 1)



dfmproval.nisin.RMSD51$value <- as.numeric(as.character(dfmproval.nisin.RMSD51$value))

dummy1 <- serprobinding[14:26,]

dfmproval.nisin.RMSD51$RESIDUE <- dummy1$Residue
dfmproval.nisin.RMSD51$RESIDUE <- factor(dfmproval.nisin.RMSD51$RESIDUE, levels = dfmproval.nisin.RMSD51$RESIDUE)

colnames(dfmproval.nisin.RMSD51) <- c("Residue", "Mutation", "RMSF", "RESIDUE")

#"RMSD of Residues after 50 nanosecond"
pdf(file="serprovalnisinRMSF50nanosText.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmproval.nisin.RMSD51,aes(x = RESIDUE,y = RMSF)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  labs(y=main4) +
  scale_x_discrete(labels=main) +
  theme(axis.text.x=element_text(angle = 60, hjust = 1,size=35),axis.text.y=element_text(angle = 0, hjust = 1,size=30),
        axis.title=element_text(size=30,face="bold"),legend.text=element_text(size = 40), legend.title=element_text(size=40,face="bold")) +
  ggtitle("")
dev.off()