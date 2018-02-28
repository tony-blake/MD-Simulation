library(xlsx)
library(ggplot2)
library(plyr)
library(reshape)
library(dplyr)


all.bridge.avg <- read.table("all.bridge.avg.dat", header =T, fill=T)
all.UU.avg <- read.table("all.UU.avg.dat", header =F, fill=T)

bbhbond <- read.table("bbhbond", header=T)



FINAL_DECOMP_MMPBSA <- read.xlsx("FINAL_DECOMP_MMPBSA.xlsx", sheetIndex = 1, header=F)
FINAL_DECOMP_MMPBSA.PRO <- read.xlsx("FINAL_DECOMP_MMPBSA.dat.pro.xlsx", sheetIndex = 1, header=F)

bindE <- FINAL_DECOMP_MMPBSA[1:26,c("X2","X18")]

bindEPRO <- FINAL_DECOMP_MMPBSA.PRO[9:34,c("X2","X18")]

bindE$PRO.blue.SER.yellow <- 1
bindEPRO$PRO.blue.SER.yellow <- 2


bindEserpro <- rbind(bindE, bindEPRO)
colnames(bindEserpro) <- c("Residue", "Binding.Energy", "PRO.blue.SER.yellow")

pdf(file="serproBindE25nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(bindEserpro, aes(x=Residue, y=Binding.Energy, group=PRO.blue.SER.yellow, col=PRO.blue.SER.yellow)) +
  geom_bar(colour="black", stat="identity",
           position=position_dodge(),
           size=.3) +                
  expand_limits(y=0) +
  scale_colour_gradientn(colours = c('yellow', 'blue')) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Binding Energies of Residues")
dev.off()



 serprobinding <- cbind(bindE, bindEPRO)
 
 serprobinding <- serprobinding[,-3] 
 
 colnames(serprobinding) <- c("Residue", "SERINE", "PROLINE")
 
 dfm <- melt(serprobinding,id.vars = 1)
 dfm <- data.frame(lapply(dfm, function(x) if(class(x)=="character") trimws(x) else(x)), stringsAsFactors=F)
 
 dfm <- data.frame(lapply(dfm, trimws))
 
pdf(file="serproBindE25nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfm,aes(x = Residue,y = value)) + 
   geom_bar(aes(fill = variable),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(size=4),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Binding Energies of Residues")
dev.off()

dfm$value <- as.numeric(as.character(dfm$value))
dat1 <- subset(dfm, value >= 0)
dat2 <- subset(dfm, value < 0)


pdf(file="serproBindE25nanos2.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot() + 
  geom_bar(data = dat1, aes(x=Residue, y=value, fill=variable),stat = "identity") +
  geom_bar(data = dat2, aes(x=Residue, y=value, fill=variable),stat = "identity") +
  theme(axis.text=element_text(size=4),
        axis.title=element_text(size=18,face="bold"),legend.text=element_text(size = 14), legend.title=element_text(size=18,face="bold")) +
  ggtitle("Binding Energies of Residues")
dev.off()




FINAL_DECOMP_MMPBSAser2 <- read.xlsx("FINAL_DECOMP_MMPBSAser2.xlsx", sheetIndex = 1, header=F)
FINAL_DECOMP_MMPBSA35.PRO <- read.xlsx("FINAL_DECOMP_MMPBSA35.xlsx", sheetIndex = 1, header=F)

bindEser2 <- FINAL_DECOMP_MMPBSAser2[9:34,c("X2","X18")]

bindE35PRO <- FINAL_DECOMP_MMPBSA35.PRO[9:34,c("X2","X18")]

bindEser2$PRO.blue.SER.yellow <- 1
bindE35PRO$PRO.blue.SER.yellow <- 2


bindE35serpro <- rbind(bindEser2, bindE35PRO)
colnames(bindE35serpro) <- c("Residue", "Binding.Energy", "PRO.blue.SER.yellow")




serprobindingtest <- cbind(bindEser2, bindE35PRO)

serprobindingtest <- serprobindingtest[,-3] 

colnames(serprobindingtest) <- c("Residue", "SERINE", "PROLINE")





dfmtest <- melt(serprobindingtest,id.vars = 1)



dfmtest$value <- as.numeric(as.character(dfmtest$value))

dfmtest$Residue <- factor(dfmtest$Residue, levels = dfmtest$Residue)

pdf(file="serproBindE35nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmtest,aes(x = Residue,y = value)) + 
  geom_bar(aes(fill = variable),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Binding Energy of Residues")
dev.off()
################################################################################################

FINAL_DECOMP_MMPBSA50.SER <- read.delim("FINAL_DECOMP_MMPBSAser50.dat", sep="", header=F)
FINAL_DECOMP_MMPBSA35.PRO <- read.delim("FINAL_DECOMP_MMPBSA50.dat", sep="", header=F)


FINAL_DECOMP_MMPBSA50.SER <- read.xlsx("FINAL_DECOMP_MMPBSAser50.xlsx", sheetIndex = 1, header=F)
FINAL_DECOMP_MMPBSA35.PRO <- read.xlsx("FINAL_DECOMP_MMPBSA50.xlsx", sheetIndex = 1, header=F)


bindEser50 <- FINAL_DECOMP_MMPBSA50.SER[9:34,c("X2","X18")]

bindE50PRO <- FINAL_DECOMP_MMPBSA35.PRO[9:34,c("X2","X18")]

bindEser2$PRO.blue.SER.yellow <- 1
bindE35PRO$PRO.blue.SER.yellow <- 2


bindE35serpro <- rbind(bindEser2, bindE35PRO)
colnames(bindE35serpro) <- c("Residue", "Binding.Energy", "PRO.blue.SER.yellow")




serprobinding50 <- cbind(bindEser50, bindE50PRO)

serprobinding50 <- serprobinding50[,-3] 

colnames(serprobinding50) <- c("Residue", "SERINE", "PROLINE")






dfm50 <- melt(serprobinding50,id.vars = 1)



dfm50$value <- as.numeric(as.character(dfm50$value))

dfm50$Residue <- factor(dfm50$Residue, levels = dfm50$Residue)

pdf(file="serproBindE50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfm50,aes(x = Residue,y = value)) + 
  geom_bar(aes(fill = variable),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Binding Energy of Residues")
dev.off()


######################################### 51 #########################################

FINAL_DECOMP_MMPBSA50.SER <- read.xlsx("FINAL_DECOMP_MMPBSAser50.xlsx", sheetIndex = 1, header=F)
FINAL_DECOMP_MMPBSA51.PRO <- read.xlsx("FINAL_DECOMP_MMPBSA51.xlsx", sheetIndex = 1, header=F)


bindEser50 <- FINAL_DECOMP_MMPBSA50.SER[9:34,c("X2","X18")]

bindE51PRO <- FINAL_DECOMP_MMPBSA51.PRO[9:34,c("X2","X18")]

bindEser2$PRO.blue.SER.yellow <- 1
bindE35PRO$PRO.blue.SER.yellow <- 2


bindE35serpro <- rbind(bindEser2, bindE35PRO)
colnames(bindE35serpro) <- c("Residue", "Binding.Energy", "PRO.blue.SER.yellow")




serprobinding51 <- cbind(bindEser50, bindE51PRO)

serprobinding51 <- serprobinding51[,-3] 

colnames(serprobinding51) <- c("Residue", "SERINE", "PROLINE")






dfm51 <- melt(serprobinding51,id.vars = 1)



dfm51$value <- as.numeric(as.character(dfm51$value))

dfm51$Residue <- factor(dfm51$Residue, levels = dfm51$Residue)

pdf(file="serproBindE51nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfm51,aes(x = Residue,y = value)) + 
  geom_bar(aes(fill = variable),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Binding Energy of Residues")
dev.off()

############################################ gfpro2 #############################################


FINAL_DECOMP_MMPBSA50.SER <- read.xlsx("FINAL_DECOMP_MMPBSAser50.xlsx", sheetIndex = 1, header=F)
FINAL_DECOMP_MMPBSA51TEST.PRO <- read.xlsx("FINAL_DECOMP_MMPBSA51TEST.xlsx", sheetIndex = 1, header=F)


bindEser50 <- FINAL_DECOMP_MMPBSA50.SER[9:34,c("X2","X18")]

bindE51PROTEST <- FINAL_DECOMP_MMPBSA51TEST.PRO[9:34,c("X2","X18")]

bindEser2$PRO.blue.SER.yellow <- 1
bindE35PRO$PRO.blue.SER.yellow <- 2


bindE35serpro <- rbind(bindEser2, bindE51PROTEST)
colnames(bindE35serpro) <- c("Residue", "Binding.Energy", "PRO.blue.SER.yellow")




serprobinding51TEST <- cbind(bindEser50, bindE51PROTEST)

serprobinding51TEST <- serprobinding51TEST[,-3] 

colnames(serprobinding51TEST) <- c("Residue", "SERINE", "PROLINE")






dfm51TEST <- melt(serprobinding51TEST,id.vars = 1)



dfm51TEST$value <- as.numeric(as.character(dfm51TEST$value))

dfm51TEST$Residue <- factor(dfm51TEST$Residue, levels = dfm51TEST$Residue)

pdf(file="serproBindETEST51nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfm51TEST,aes(x = Residue,y = value)) + 
  geom_bar(aes(fill = variable),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Binding Energy of Residues")
dev.off()

###################################### proval #############################################

FINAL_DECOMP_MMPBSA50PROVAL <- read.xlsx("FINAL_DECOMP_MMPBSA50PROVAL.xlsx", sheetIndex = 1, header=F)

bindE51PROVAL <- FINAL_DECOMP_MMPBSA50PROVAL[9:34,c("X2","X18")]

serprovalbinding50 <- cbind(bindEser50, bindE51PROVAL)

serprovalbinding50 <- serprovalbinding50[,-3] 

colnames(serprovalbinding50) <- c("Residue", "SERINE", "PROLINE")






dfmprovalBE <- melt(serprovalbinding50,id.vars = 1)



dfmprovalBE$value <- as.numeric(as.character(dfmprovalBE$value))

dfmprovalBE$Residue <- factor(dfmprovalBE$Residue, levels = dfmprovalBE$Residue)

pdf(file="serprovalBindE50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmprovalBE,aes(x = Residue,y = value)) + 
  geom_bar(aes(fill = variable),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Binding Energy of Residues")
dev.off()

###################################### alanine #############################################

FINAL_DECOMP_MMPBSA50ALA <- read.xlsx("FINAL_DECOMP_ALA50_MMPBSA.xlsx", sheetIndex = 1, header=F)

bindE50ALA <- FINAL_DECOMP_MMPBSA50ALA[9:34,c("X2","X18")]

alabinding50 <- cbind(bindEser50, bindE50ALA)

alabinding50 <- alabinding50[,-3] 

colnames(alabinding50) <- c("Residue", "SERINE", "ALANINE")






dfmalaBE <- melt(alabinding50,id.vars = 1)



dfmalaBE$value <- as.numeric(as.character(dfmalaBE$value))

dfmalaBE$Residue <- factor(dfmalaBE$Residue, levels = dfmalaBE$Residue)

pdf(file="alaBindE50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmalaBE,aes(x = Residue,y = value)) + 
  geom_bar(aes(fill = variable),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Binding Energy of Residues")
dev.off()


###################################### proline alanine #############################################

FINAL_DECOMP_MMPBSA50ALA <- read.xlsx("FINAL_DECOMP_ALA50_MMPBSA.xlsx", sheetIndex = 1, header=F)

bindE50ALA <- FINAL_DECOMP_MMPBSA50ALA[9:34,c("X2","X18")]

proalabinding50 <- cbind(bindEser50, bindE50ALA, bindE51PROVAL)

proalabinding50 <- proalabinding50[,-3] 
proalabinding50 <- proalabinding50[,-4] 

colnames(proalabinding50) <- c("Residue", "SERINE", "ALANINE", "PROLINE")






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


###################################### proval #############################################

FINAL_DECOMP_MMPBSA50PROVAL <- read.xlsx("FINAL_DECOMP_MMPBSA50PROVAL.xlsx", sheetIndex = 1, header=F)

bindE51PROVAL <- FINAL_DECOMP_MMPBSA50PROVAL[9:34,c("X2","X18")]

serprovalbinding50 <- cbind(bindEser50, bindE51PROVAL)

serprovalbinding50 <- serprovalbinding50[,-3] 

colnames(serprovalbinding50) <- c("Residue", "Nisin A", "Nisin PV")






dfmprovalBE <- melt(serprovalbinding50,id.vars = 1)



dfmprovalBE$value <- as.numeric(as.character(dfmprovalBE$value))

dfmprovalBE$Residue <- factor(dfmprovalBE$Residue, levels = dfmprovalBE$Residue)

colnames(dfmprovalBE) <- c("Residue", "Mutation", "value")

pdf(file="serprovalBindE50nanosTEXT.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmprovalBE,aes(x = Residue,y = value)) +
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
 labs(y=main4) +
  scale_x_discrete(labels=main5) +
  theme(axis.text.x=element_text(angle = 90, hjust = 1,size=10),axis.text.y=element_text(angle = 0, hjust = 1,size=30),
        axis.title=element_text(size=20,face="bold"),legend.text=element_text(size = 40), legend.title=element_text(size=40,face="bold")) +
  ggtitle("Binding Energy of Residues")
dev.off()


main4=expression('Net Energy/Res'[kcal/mol])




main=expression('Lys'[22], 'Dbb'[23], 'Ala'[24], 'Dbb'[25], 'Cys'[26], 'Asn'[27], 'Cys'[28], 'Ser Pro'[29], 'ile Val'[30], 'His'[31], 'Val'[32], 'Dha'[33], 'Lys'[34])

main2=expression('Thr'[230], 'Asn'[231], 'Hie'[232], 'Lys'[233], 'Thr'[234], 'Ala'[235], 'Ser'[236], 'Ser'[237], 'Ala'[238], 'Glu'[239], 'Met'[240], 'Thr'[241], 'Phe'[242])


main5=expression('Thr'[230], 'Asn'[231], 'Hie'[232], 'Lys'[233], 'Thr'[234], 'Ala'[235], 'Ser'[236], 'Ser'[237], 'Ala'[238], 'Glu'[239], 'Met'[240], 'Thr'[241], 'Phe'[242], 'Lys'[22], 'Dbb'[23], 'Ala'[24], 'Dbb'[25], 'Cys'[26], 'Asn'[27], 'Cys'[28], 'Ser Pro'[29], 'ile Val'[30], 'His'[31], 'Val'[32], 'Dha'[33], 'Lys'[34])

###################################### proval split #############################################

FINAL_DECOMP_MMPBSA50PROVAL <- read.xlsx("FINAL_DECOMP_MMPBSA50PROVAL.xlsx", sheetIndex = 1, header=F)

bindE51PROVAL <- FINAL_DECOMP_MMPBSA50PROVAL[9:34,c("X2","X18")]

serprovalbinding50 <- cbind(bindEser50, bindE51PROVAL)

serprovalbinding50 <- serprovalbinding50[,-3] 

colnames(serprovalbinding50) <- c("Residue", "SERINE", "PROLINE")






dfmprovalBE <- melt(serprovalbinding50,id.vars = 1)



dfmprovalBE$value <- as.numeric(as.character(dfmprovalBE$value))

dfmprovalBE$Residue <- factor(dfmprovalBE$Residue, levels = dfmprovalBE$Residue)

dfmprovalsplit1 <- dfmprovalBE[c(14:26,40:52),]
dfmprovalsplit2 <- dfmprovalBE[c(1:13,27:39),]

pdf(file="serprovalBindE50nanosTEXT1.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmprovalsplit1,aes(x = Residue,y = value)) +
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  labs(y=main4) +
  scale_x_discrete(labels=main) +
  theme(axis.text.x=element_text(angle = 60, hjust = 1,size=35),axis.text.y=element_text(angle = 0, hjust = 1,size=30),
        axis.title=element_text(size=40,face="bold"),legend.text=element_text(size = 40), legend.title=element_text(size=40,face="bold")) +
  ggtitle("")
dev.off()



pdf(file="serprovalBindE50nanosTEXT2.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmprovalsplit2,aes(x = Residue,y = value)) +
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  labs(y=main4) +
  scale_x_discrete(labels=main2) +
  theme(axis.text.x=element_text(angle = 60, hjust = 1,size=35),axis.text.y=element_text(angle = 0, hjust = 1,size=30),
        axis.title=element_text(size=25,face="bold"),legend.text=element_text(size = 40), legend.title=element_text(size=40,face="bold")) +
  ggtitle("")
dev.off()
