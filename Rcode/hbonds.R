





#Input raw data and create dataframes
hbondtable51proval <- read.delim("All.UU.avg.all10proval.dat", sep="", stringsAsFactors=FALSE)
hbondtable50sernisin <- read.delim("hbondsnisinSER2", sep="", stringsAsFactors=FALSE, header = FALSE)
hbondtable51provalnisin <- read.delim("hbondprovalsnisinPROVAL2", sep="", stringsAsFactors=FALSE, header = FALSE)
colnames(hbondtable51provalnisin) <- colnames(hbondtable50pro)
colnames(hbondtable50sernisin) <- colnames(hbondtable50pro)
hbondtable50sernisin$Mutation <- "Nisin A"
hbondtable51provalnisin$Mutation <- "Nisin PV"

#Combine dataframes for Nisin A and Nisin PV and format for input to ggplot2
serproval.nisin.hbonds51 <- rbind(hbondtable50sernisin, hbondtable51provalnisin)
serproval.nisin.hbonds.temp51 <- within(serproval.nisin.hbonds51, bond <- paste(X.Acceptor,DonorH,sep='-'))
serproval.nisin.hbonds.neo51 <- serproval.nisin.hbonds.temp51[,c("bond", "Frac", "Mutation")]
serproval.nisin.hbonds.neo51$bond <- factor(serproval.nisin.hbonds.neo51$bond, levels = serproval.nisin.hbonds.neo51$bond)

#Create variable to store y-axis title text
main4=expression('Hbond Occupancy')

#Create plot to show hbond occupancies of Residues for every residue 
pdf(file="serprovalnisinHBONDS51nanosText.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproval.nisin.hbonds.neo51,aes(x = bond,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  labs(y=main4) +
  theme(axis.text.x=element_text(angle = 60, hjust = 1,size=5),axis.text.y=element_text(angle = 0, hjust = 1,size=30),
        axis.title=element_text(size=30,face="bold"),legend.text=element_text(size = 40), legend.title=element_text(size=40,face="bold")) +
  ggtitle("")
dev.off()

#Create new dataframe to only show residues whose hbond occupancy is greater than 20% 
serproval.nisin.hbonds.neo51$Frac > 0.2
serproval.nisin.hbonds.neo51.greatertwenty <- serproval.nisin.hbonds.neo51[serproval.nisin.hbonds.neo51$Frac > 0.2,]
serproval.nisin.hbonds.neo51.greatertwenty$bond <- factor(serproval.nisin.hbonds.neo51.greatertwenty$bond, levels = serproval.nisin.hbonds.neo51.greatertwenty$bond)

#Create variable to store y-axis title text
main4=expression('Occupancy')

#Create plot to show hbond occupancies of Residues for residues with occupancies > 20% 
pdf(file="serprovalnisinHBONDS51nanosText20.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproval.nisin.hbonds.neo51.greatertwenty,aes(x = bond,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  labs(y=main4) +
  theme(axis.text.x=element_text(angle = 60, hjust = 1,size=20),axis.text.y=element_text(angle = 0, hjust = 1,size=40),
        axis.title=element_text(size=40,face="bold"),legend.text=element_text(size = 40), legend.title=element_text(size=40,face="bold")) +
  ggtitle("")
dev.off()
