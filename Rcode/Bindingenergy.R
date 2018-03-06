library(xlsx)
library(ggplot2)
library(plyr)
library(reshape)
library(dplyr)


#Read in raw data
FINAL_DECOMP_MMPBSA50.SER <- read.xlsx("FINAL_DECOMP_MMPBSAser50.xlsx", sheetIndex = 1, header=F)
FINAL_DECOMP_MMPBSA50PROVAL <- read.xlsx("FINAL_DECOMP_MMPBSA50PROVAL.xlsx", sheetIndex = 1, header=F)

#Extract relevant columns
bindEser50 <- FINAL_DECOMP_MMPBSA50.SER[9:34,c("X2","X18")]
bindE51PROVAL <- FINAL_DECOMP_MMPBSA50PROVAL[9:34,c("X2","X18")]

#Combine extracted column data forNisin-A and Nisin-PV
serprovalbinding50 <- cbind(bindEser50, bindE51PROVAL)
serprovalbinding50 <- serprovalbinding50[,-3] 

#Change dataframe structure for ggplot2 input format
dfmprovalBE <- melt(serprovalbinding50,id.vars = 1)
dfmprovalBE$value <- as.numeric(as.character(dfmprovalBE$value))
dfmprovalBE$Residue <- factor(dfmprovalBE$Residue, levels = dfmprovalBE$Residue)


#Split datframe into values for Nisin residues and NSR residues repsectively
dfmprovalsplit1 <- dfmprovalBE[c(14:26,40:52),]
dfmprovalsplit2 <- dfmprovalBE[c(1:13,27:39),]

#Encode figure axis text in variable for ggplot2 functions
main4=expression('Net Energy/Res'[kcal/mol])
main=expression('Lys'[22], 'Dbb'[23], 'Ala'[24], 'Dbb'[25], 'Cys'[26], 'Asn'[27], 'Cys'[28], 'Ser Pro'[29], 'ile Val'[30], 'His'[31], 'Val'[32], 'Dha'[33], 'Lys'[34])
main2=expression('Thr'[230], 'Asn'[231], 'Hie'[232], 'Lys'[233], 'Thr'[234], 'Ala'[235], 'Ser'[236], 'Ser'[237], 'Ala'[238], 'Glu'[239], 'Met'[240], 'Thr'[241], 'Phe'[242])
main5=expression('Thr'[230], 'Asn'[231], 'Hie'[232], 'Lys'[233], 'Thr'[234], 'Ala'[235], 'Ser'[236], 'Ser'[237], 'Ala'[238], 'Glu'[239], 'Met'[240], 'Thr'[241], 'Phe'[242], 'Lys'[22], 'Dbb'[23], 'Ala'[24], 'Dbb'[25], 'Cys'[26], 'Asn'[27], 'Cys'[28], 'Ser Pro'[29], 'ile Val'[30], 'His'[31], 'Val'[32], 'Dha'[33], 'Lys'[34])


#Create bar plot figure for binding energies of nisin residues
pdf(file="serprovalBindE50nanosTEXT1.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmprovalsplit1,aes(x = Residue,y = value)) +
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  labs(y=main4) +
  scale_x_discrete(labels=main) +
  theme(axis.text.x=element_text(angle = 60, hjust = 1,size=35),axis.text.y=element_text(angle = 0, hjust = 1,size=30),
        axis.title=element_text(size=40,face="bold"),legend.text=element_text(size = 40), legend.title=element_text(size=40,face="bold")) +
  ggtitle("")
dev.off()


#Create bar plot figure for binding energies of NSR residues
pdf(file="serprovalBindE50nanosTEXT2.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(dfmprovalsplit2,aes(x = Residue,y = value)) +
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  labs(y=main4) +
  scale_x_discrete(labels=main2) +
  theme(axis.text.x=element_text(angle = 60, hjust = 1,size=35),axis.text.y=element_text(angle = 0, hjust = 1,size=30),
        axis.title=element_text(size=25,face="bold"),legend.text=element_text(size = 40), legend.title=element_text(size=40,face="bold")) +
  ggtitle("")
dev.off()
