library(xlsx)
library(ggplot2)
library(plyr)
library(reshape)
library(dplyr)



########################### RMSD proval over time (RMSF) #######################################

#Input raw data and make dataframes
rms_nisinser <- read.delim("rmsdnisinovertime", sep="", stringsAsFactors=FALSE, header = TRUE)
rms_nisinproval <- read.delim("rmsdpvovertime", sep="", stringsAsFactors=FALSE, header = TRUE)
colnames(rms_nisinser) <- c("Picoseconds", "Angstroms")
colnames(rms_nisinproval) <- c("Picoseconds", "Angstroms")

#combine Nisin A and Nisin PV dataframes
rms.serproval <- cbind(rms_nisinser, rms_nisinproval)
rms.serproval <- rms.serproval[,-3] 
colnames(rms.serproval) <- c("Residue", "Nisin A", "Nisin PV")


#Change dataframe structure for ggplot2 input
dfmrmspv <- melt(rms.serproval,id.vars = 1)
dfmrmspvnano <- dfmrmspv[,1]/1000
dfmrmspvnano <- data.frame(dfmrmspvnano)
dfmrmspvnano <- cbind(dfmpvnano,dfmrmspv)
dfmrmspvnano <- dfmrmspvnano[,-c(2,3,4)]
colnames(dfmrmspvnano) <- c("Nanoseconds", "Mutation", "Angstroms")


#create variables to hold axix title text in ggplot2 figures
main4=expression(RMSD[ring(A)])
main=expression('Lys'[22], 'Dbb'[23], 'Ala'[24], 'Dbb'[25], 'Cys'[26], 'Asn'[27], 'Cys'[28], 'Ser Pro'[29], 'ile Val'[30], 'His'[31], 'Val'[32], 'Dha'[33], 'Lys'[34])
main2=expression('Thr'[230], 'Asn'[231], 'Hie'[232], 'Lys'[233], 'Thr'[234], 'Ala'[235], 'Ser'[236], 'Ser'[237], 'Ala'[238], 'Glu'[239], 'Met'[240], 'Thr'[241], 'Phe'[242])
main5=expression('Thr'[230], 'Asn'[231], 'Hie'[232], 'Lys'[233], 'Thr'[234], 'Ala'[235], 'Ser'[236], 'Ser'[237], 'Ala'[238], 'Glu'[239], 'Met'[240], 'Thr'[241], 'Phe'[242], 'Lys'[22], 'Dbb'[23], 'Ala'[24], 'Dbb'[25], 'Cys'[26], 'Asn'[27], 'Cys'[28], 'Ser Pro'[29], 'ile Val'[30], 'His'[31], 'Val'[32], 'Dha'[33], 'Lys'[34])



 

#Create figure that shows RMSF between cleavage point nitrogen in Residue 28 in Nisin and the Sidechain O in SER236 in NSR over 50 nanoeconds
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




################################# RMSD proval #######################################################3#####

#Input raw data and make dataframes
rmsd50ser <- read.delim("perresavg.BB.50.dat", sep="", stringsAsFactors=FALSE, header = TRUE)
rmsd51proval <- read.delim("perresavg.BB.50proval.dat", sep="", stringsAsFactors=FALSE, header = TRUE)
rmsd.nisin50ser <- rmsd50ser[288:300,]
rmsdproval.nisin51proval <- rmsd51proval[288:300,]
colnames(rmsd.nisin50ser) <- c("Residue", "RMSD", "std.dev")
colnames(rmsdproval.nisin51proval) <- c("Residue", "RMSD", "std.dev")


#Combine datframes for Nisin A and Nisin PV
serproval.nisin.RMSD51 <- cbind(rmsd.nisin50ser, rmsdproval.nisin51proval)
serproval.nisin.RMSD51 <- serproval.nisin.RMSD51[,-c(3,4,6)] 
colnames(serproval.nisin.RMSD51) <- c("Residue", "Nisin A", "Nisin PV")

#modify datframe for input
dfmproval.nisin.RMSD51 <- melt(serproval.nisin.RMSD51,id.vars = 1)
dfmproval.nisin.RMSD51$value <- as.numeric(as.character(dfmproval.nisin.RMSD51$value))

#Use column from dataframe in Bindingenergy.R script "serprobinding" to list residue names
dummy1 <- serprobinding[14:26,]
dfmproval.nisin.RMSD51$RESIDUE <- dummy1$Residue
dfmproval.nisin.RMSD51$RESIDUE <- factor(dfmproval.nisin.RMSD51$RESIDUE, levels = dfmproval.nisin.RMSD51$RESIDUE)
colnames(dfmproval.nisin.RMSD51) <- c("Residue", "Mutation", "RMSF", "RESIDUE")

#Create plot of RMSD of Residues after 50 nanosecond
pdf(file="serprovalnisinRMSF50nanosText.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmproval.nisin.RMSD51,aes(x = RESIDUE,y = RMSF)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  labs(y=main4) +
  scale_x_discrete(labels=main) +
  theme(axis.text.x=element_text(angle = 60, hjust = 1,size=35),axis.text.y=element_text(angle = 0, hjust = 1,size=30),
        axis.title=element_text(size=30,face="bold"),legend.text=element_text(size = 40), legend.title=element_text(size=40,face="bold")) +
  ggtitle("")
dev.off()
