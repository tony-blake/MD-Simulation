hbondtable50pro <- read.delim("All.UU.avg.all112pro.dat", sep="", stringsAsFactors=FALSE)

hbondtable50sernisin <- read.delim("hbondsnisinSER", sep="", stringsAsFactors=FALSE, header = FALSE)

hbondtable50pronisin <- read.delim("hbondsnisinPRO", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(hbondtable50pronisin) <- colnames(hbondtable50pro)

hbondtable50sernisin$Mutation <- "SERINE"
hbondtable50pronisin$Mutation <- "PROLINE"

serpro.nisin.hbonds <- rbind(hbondtable50sernisin, hbondtable50pronisin)

serpro.nisin.hbonds.temp <- within(serpro.nisin.hbonds, x <- paste(X.Acceptor,DonorH,sep='-'))

serpro.nisin.hbonds.neo <- serpro.nisin.hbonds.temp[,c("x", "Frac", "Mutation")]


serpro.nisin.hbonds.neo$x <- factor(serpro.nisin.hbonds.neo$x, levels = serpro.nisin.hbonds.neo$x)

pdf(file="serpronisinHBONDS50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro.nisin.hbonds.neo,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=2),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues")
dev.off()


############################ nisin core ###############################

corehbondtable50sernisin <- read.delim("hbondcoreser", sep="", stringsAsFactors=FALSE, header = FALSE)

corehbondtable50pronisin <- read.delim("hbondcorepro", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(corehbondtable50pronisin) <- colnames(hbondtable50pro)
colnames(corehbondtable50sernisin) <- colnames(hbondtable50pro)

corehbondtable50sernisin$Mutation <- "SERINE"
corehbondtable50pronisin$Mutation <- "PROLINE"

serpro.nisincore.hbonds <- rbind(corehbondtable50sernisin, corehbondtable50pronisin)

serpro.nisincore.hbonds.temp <- within(serpro.nisincore.hbonds, x <- paste(X.Acceptor,DonorH,sep='-'))

serpro.nisincore.hbonds.neo <- serpro.nisincore.hbonds.temp[,c("x", "Frac", "Mutation")]


serpro.nisincore.hbonds.neo$x <- factor(serpro.nisincore.hbonds.neo$x, levels = serpro.nisincore.hbonds.neo$x)

pdf(file="serpronisincoreHBONDS50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro.nisincore.hbonds.neo,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues CYS28-SER29(PRO29)-ILE30")
dev.off()


############################ nisin C-Terminus ###############################

ctermhbondtable50sernisin <- read.delim("hbondCTERMINUSser", sep="", stringsAsFactors=FALSE, header = FALSE)

ctermhbondtable50pronisin <- read.delim("hbondCTERMINUSpro", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(ctermhbondtable50sernisin) <- colnames(hbondtable50pro)
colnames(ctermhbondtable50pronisin) <- colnames(hbondtable50pro)

ctermhbondtable50sernisin$Mutation <- "SERINE"
ctermhbondtable50pronisin$Mutation <- "PROLINE"

serpro.nisincterm.hbonds <- rbind(ctermhbondtable50sernisin, ctermhbondtable50pronisin)

serpro.nisincterm.hbonds.temp <- within(serpro.nisincterm.hbonds, x <- paste(X.Acceptor,DonorH,sep='-'))

serpro.nisincterm.hbonds.neo <- serpro.nisincterm.hbonds.temp[,c("x", "Frac", "Mutation")]


serpro.nisincterm.hbonds.neo$x <- factor(serpro.nisincterm.hbonds.neo$x, levels = serpro.nisincterm.hbonds.neo$x)

pdf(file="serpronisinCTERMINUSHBONDS50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro.nisincterm.hbonds.neo,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues HIE31-VAL32-DHA33")
dev.off()

################################ nisin Lyseine 300 ################################################

cLYS300hbondtable50sernisin <- read.delim("hbondLYS300ser", sep="", stringsAsFactors=FALSE, header = FALSE)

cLYS300hbondtable50pronisin <- read.delim("hbondLYS300pro", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(cLYS300hbondtable50sernisin) <- colnames(hbondtable50pro)
colnames(cLYS300hbondtable50pronisin ) <- colnames(hbondtable50pro)

cLYS300hbondtable50sernisin$Mutation <- "SERINE"
cLYS300hbondtable50pronisin$Mutation <- "PROLINE"

serpro.nisincLYS300.hbonds <- rbind(cLYS300hbondtable50sernisin, cLYS300hbondtable50pronisin) 
serpro.nisincLYS300.hbonds.temp <- within(serpro.nisincLYS300.hbonds, x <- paste(X.Acceptor,DonorH,sep='-'))

serpro.nisincLYS300.hbonds.neo <- serpro.nisincLYS300.hbonds.temp[,c("x", "Frac", "Mutation")]


serpro.nisincLYS300.hbonds.neo$x <- factor(serpro.nisincLYS300.hbonds.neo$x, levels = serpro.nisincLYS300.hbonds.neo$x)

pdf(file="serpronisinCLYS300HBONDS50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro.nisincLYS300.hbonds.neo,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues LYS34 ")
dev.off()

################################ nisin ring E ################################################

ringEhbondtable50sernisin <- read.delim("hbondsringE", sep="", stringsAsFactors=FALSE, header = FALSE)

ringEhbondtable50pronisin <- read.delim("hbondsringEpro", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(ringEhbondtable50sernisin) <- colnames(hbondtable50pro)
colnames(ringEhbondtable50pronisin) <- colnames(hbondtable50pro)

ringEhbondtable50sernisin$Mutation <- "SERINE"
ringEhbondtable50pronisin$Mutation <- "PROLINE"

serpro.nisinRingE.hbonds <- rbind(ringEhbondtable50sernisin, ringEhbondtable50pronisin) 
serpro.nisinRingE.hbonds.temp <- within(serpro.nisinRingE.hbonds, x <- paste(X.Acceptor,DonorH,sep='-'))

serpro.nisinRingE.hbonds.neo <- serpro.nisinRingE.hbonds.temp[,c("x", "Frac", "Mutation")]


serpro.nisinRingE.hbonds.neo$x <- factor(serpro.nisinRingE.hbonds.neo$x, levels = serpro.nisinRingE.hbonds.neo$x)

pdf(file="serpronisinRingEHBONDS50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro.nisinRingE.hbonds.neo,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues DBB25-CYS26-ASN27")
dev.off()


################################ nisin ring D ################################################

ringDhbondtable50sernisin <- read.delim("hbondsringDser", sep="", stringsAsFactors=FALSE, header = FALSE)

ringDhbondtable50pronisin <- read.delim("hbondsringDpro", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(ringDhbondtable50sernisin) <- colnames(hbondtable50pro)
colnames(ringDhbondtable50pronisin) <- colnames(hbondtable50pro)

ringDhbondtable50sernisin$Mutation <- "SERINE"
ringDhbondtable50pronisin$Mutation <- "PROLINE"

serpro.nisinRingD.hbonds <- rbind(ringDhbondtable50sernisin, ringDhbondtable50pronisin) 
serpro.nisinRingD.hbonds.temp <- within(serpro.nisinRingD.hbonds, x <- paste(X.Acceptor,DonorH,sep='-'))

serpro.nisinRingD.hbonds.neo <- serpro.nisinRingD.hbonds.temp[,c("x", "Frac", "Mutation")]


serpro.nisinRingD.hbonds.neo$x <- factor(serpro.nisinRingD.hbonds.neo$x, levels = serpro.nisinRingD.hbonds.neo$x)

pdf(file="serpronisinRingDHBONDS50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro.nisinRingD.hbonds.neo,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues LYS22-DBB23-ALA24")
dev.off()



########################################## pro 51 ##############################################

hbondtable51pro <- read.delim("All.UU.avg.all113pro.dat", sep="", stringsAsFactors=FALSE)

hbondtable50sernisin <- read.delim("hbondsnisinSER", sep="", stringsAsFactors=FALSE, header = FALSE)

hbondtable51pronisin <- read.delim("hbondsnisinPRO51", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(hbondtable51pronisin) <- colnames(hbondtable50pro)

hbondtable50sernisin$Mutation <- "SERINE"
hbondtable51pronisin$Mutation <- "PROLINE"

serpro.nisin.hbonds51 <- rbind(hbondtable50sernisin, hbondtable51pronisin)

serpro.nisin.hbonds.temp51 <- within(serpro.nisin.hbonds51, x <- paste(X.Acceptor,DonorH,sep='-'))

serpro.nisin.hbonds.neo51 <- serpro.nisin.hbonds.temp51[,c("x", "Frac", "Mutation")]


serpro.nisin.hbonds.neo51$x <- factor(serpro.nisin.hbonds.neo51$x, levels = serpro.nisin.hbonds.neo51$x)

pdf(file="serpronisinHBONDS51nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro.nisin.hbonds.neo51,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=2),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues")
dev.off()


############################ nisin core ###############################

corehbondtable50sernisin <- read.delim("hbondcoreser", sep="", stringsAsFactors=FALSE, header = FALSE)

corehbondtable51pronisin <- read.delim("hbondcorepro", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(corehbondtable51pronisin) <- colnames(hbondtable50pro)
colnames(corehbondtable50sernisin) <- colnames(hbondtable50pro)

corehbondtable50sernisin$Mutation <- "SERINE"
corehbondtable51pronisin$Mutation <- "PROLINE"

serpro51.nisincore.hbonds <- rbind(corehbondtable50sernisin, corehbondtable51pronisin)

serpro51.nisincore.hbonds.temp <- within(serpro51.nisincore.hbonds, x <- paste(X.Acceptor,DonorH,sep='-'))

serpro51.nisincore.hbonds.neo <- serpro51.nisincore.hbonds.temp[,c("x", "Frac", "Mutation")]


serpro51.nisincore.hbonds.neo$x <- factor(serpro51.nisincore.hbonds.neo$x, levels = serpro51.nisincore.hbonds.neo$x)

pdf(file="serpronisincoreHBONDS51nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro51.nisincore.hbonds.neo,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues CYS28-SER29(PRO29)-ILE30")
dev.off()


############################ nisin C-Terminus ###############################

ctermhbondtable50sernisin <- read.delim("hbondCTERMINUSser", sep="", stringsAsFactors=FALSE, header = FALSE)

ctermhbondtable51pronisin <- read.delim("hbondCTERMINUSpro", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(ctermhbondtable50sernisin) <- colnames(hbondtable50pro)
colnames(ctermhbondtable51pronisin) <- colnames(hbondtable51pro)

ctermhbondtable50sernisin$Mutation <- "SERINE"
ctermhbondtable51pronisin$Mutation <- "PROLINE"

serpro51.nisincterm.hbonds <- rbind(ctermhbondtable50sernisin, ctermhbondtable51pronisin)

serpro51.nisincterm.hbonds.temp <- within(serpro51.nisincterm.hbonds, x <- paste(X.Acceptor,DonorH,sep='-'))

serpro51.nisincterm.hbonds.neo <- serpro51.nisincterm.hbonds.temp[,c("x", "Frac", "Mutation")]


serpro51.nisincterm.hbonds.neo$x <- factor(serpro.nisincterm.hbonds.neo$x, levels = serpro.nisincterm.hbonds.neo$x)

pdf(file="serpronisinCTERMINUSHBONDS51nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro51.nisincterm.hbonds.neo,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues HIE31-VAL32-DHA33")
dev.off()

################################ nisin Lyseine 300 ################################################

cLYS300hbondtable50sernisin <- read.delim("hbondLYS300ser", sep="", stringsAsFactors=FALSE, header = FALSE)

cLYS300hbondtable50pronisin <- read.delim("hbondLYS300pro", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(cLYS300hbondtable50sernisin) <- colnames(hbondtable50pro)
colnames(cLYS300hbondtable50pronisin ) <- colnames(hbondtable50pro)

cLYS300hbondtable50sernisin$Mutation <- "SERINE"
cLYS300hbondtable50pronisin$Mutation <- "PROLINE"

serpro.nisincLYS300.hbonds <- rbind(cLYS300hbondtable50sernisin, cLYS300hbondtable50pronisin) 
serpro.nisincLYS300.hbonds.temp <- within(serpro.nisincLYS300.hbonds, x <- paste(X.Acceptor,DonorH,sep='-'))

serpro.nisincLYS300.hbonds.neo <- serpro.nisincLYS300.hbonds.temp[,c("x", "Frac", "Mutation")]


serpro.nisincLYS300.hbonds.neo$x <- factor(serpro.nisincLYS300.hbonds.neo$x, levels = serpro.nisincLYS300.hbonds.neo$x)

pdf(file="serpronisinCLYS300HBONDS50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro.nisincLYS300.hbonds.neo,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues LYS34 ")
dev.off()

################################ nisin ring E ################################################

ringEhbondtable50sernisin <- read.delim("hbondsringE", sep="", stringsAsFactors=FALSE, header = FALSE)

ringEhbondtable50pronisin <- read.delim("hbondsringEpro", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(ringEhbondtable50sernisin) <- colnames(hbondtable50pro)
colnames(ringEhbondtable50pronisin) <- colnames(hbondtable50pro)

ringEhbondtable50sernisin$Mutation <- "SERINE"
ringEhbondtable50pronisin$Mutation <- "PROLINE"

serpro.nisinRingE.hbonds <- rbind(ringEhbondtable50sernisin, ringEhbondtable50pronisin) 
serpro.nisinRingE.hbonds.temp <- within(serpro.nisinRingE.hbonds, x <- paste(X.Acceptor,DonorH,sep='-'))

serpro.nisinRingE.hbonds.neo <- serpro.nisinRingE.hbonds.temp[,c("x", "Frac", "Mutation")]


serpro.nisinRingE.hbonds.neo$x <- factor(serpro.nisinRingE.hbonds.neo$x, levels = serpro.nisinRingE.hbonds.neo$x)

pdf(file="serpronisinRingEHBONDS50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro.nisinRingE.hbonds.neo,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues DBB25-CYS26-ASN27")
dev.off()


################################ nisin ring D ################################################

ringDhbondtable50sernisin <- read.delim("hbondsringDser", sep="", stringsAsFactors=FALSE, header = FALSE)

ringDhbondtable50pronisin <- read.delim("hbondsringDpro", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(ringDhbondtable50sernisin) <- colnames(hbondtable50pro)
colnames(ringDhbondtable50pronisin) <- colnames(hbondtable50pro)

ringDhbondtable50sernisin$Mutation <- "SERINE"
ringDhbondtable50pronisin$Mutation <- "PROLINE"

serpro.nisinRingD.hbonds <- rbind(ringDhbondtable50sernisin, ringDhbondtable50pronisin) 
serpro.nisinRingD.hbonds.temp <- within(serpro.nisinRingD.hbonds, x <- paste(X.Acceptor,DonorH,sep='-'))

serpro.nisinRingD.hbonds.neo <- serpro.nisinRingD.hbonds.temp[,c("x", "Frac", "Mutation")]


serpro.nisinRingD.hbonds.neo$x <- factor(serpro.nisinRingD.hbonds.neo$x, levels = serpro.nisinRingD.hbonds.neo$x)

pdf(file="serpronisinRingDHBONDS50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro.nisinRingD.hbonds.neo,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues LYS22-DBB23-ALA24")
dev.off()

######################################## 51 ##########################################################

############################ nisin core ###############################

corehbondtable50sernisin <- read.delim("hbondcoreser", sep="", stringsAsFactors=FALSE, header = FALSE)

corehbondtable50pronisin <- read.delim("hbondcorepro", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(corehbondtable50pronisin) <- colnames(hbondtable50pro)
colnames(corehbondtable50sernisin) <- colnames(hbondtable50pro)

corehbondtable50sernisin$Mutation <- "SERINE"
corehbondtable50pronisin$Mutation <- "PROLINE"

serpro.nisincore.hbonds <- rbind(corehbondtable50sernisin, corehbondtable50pronisin)

serpro.nisincore.hbonds.temp <- within(serpro.nisincore.hbonds, x <- paste(X.Acceptor,DonorH,sep='-'))

serpro.nisincore.hbonds.neo <- serpro.nisincore.hbonds.temp[,c("x", "Frac", "Mutation")]


serpro.nisincore.hbonds.neo$x <- factor(serpro.nisincore.hbonds.neo$x, levels = serpro.nisincore.hbonds.neo$x)

pdf(file="serpronisincoreHBONDS50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro.nisincore.hbonds.neo,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues CYS28-SER29(PRO29)-ILE30")
dev.off()


############################ nisin C-Terminus ###############################

ctermhbondtable50sernisin <- read.delim("hbondCTERMINUSser", sep="", stringsAsFactors=FALSE, header = FALSE)

ctermhbondtable50pronisin <- read.delim("hbondCTERMINUSpro", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(ctermhbondtable50sernisin) <- colnames(hbondtable50pro)
colnames(ctermhbondtable50pronisin) <- colnames(hbondtable50pro)

ctermhbondtable50sernisin$Mutation <- "SERINE"
ctermhbondtable50pronisin$Mutation <- "PROLINE"

serpro.nisincterm.hbonds <- rbind(ctermhbondtable50sernisin, ctermhbondtable50pronisin)

serpro.nisincterm.hbonds.temp <- within(serpro.nisincterm.hbonds, x <- paste(X.Acceptor,DonorH,sep='-'))

serpro.nisincterm.hbonds.neo <- serpro.nisincterm.hbonds.temp[,c("x", "Frac", "Mutation")]


serpro.nisincterm.hbonds.neo$x <- factor(serpro.nisincterm.hbonds.neo$x, levels = serpro.nisincterm.hbonds.neo$x)

pdf(file="serpronisinCTERMINUSHBONDS50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro.nisincterm.hbonds.neo,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues HIE31-VAL32-DHA33")
dev.off()

################################ nisin Lyseine 300 ################################################

cLYS300hbondtable50sernisin <- read.delim("hbondLYS300ser", sep="", stringsAsFactors=FALSE, header = FALSE)

cLYS300hbondtable51pronisin <- read.delim("hbondLYS300pro", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(cLYS300hbondtable50sernisin) <- colnames(hbondtable50pro)
colnames(cLYS300hbondtable51pronisin ) <- colnames(hbondtable51pro)

cLYS300hbondtable50sernisin$Mutation <- "SERINE"
cLYS300hbondtable51pronisin$Mutation <- "PROLINE"

serpro51.nisincLYS300.hbonds <- rbind(cLYS300hbondtable50sernisin, cLYS300hbondtable51pronisin) 
serpro51.nisincLYS300.hbonds.temp <- within(serpro51.nisincLYS300.hbonds, x <- paste(X.Acceptor,DonorH,sep='-'))

serpro51.nisincLYS300.hbonds.neo <- serpro51.nisincLYS300.hbonds.temp[,c("x", "Frac", "Mutation")]


serpro51.nisincLYS300.hbonds.neo$x <- factor(serpro51.nisincLYS300.hbonds.neo$x, levels = serpro51.nisincLYS300.hbonds.neo$x)

pdf(file="serpronisinCLYS300HBONDS51nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro51.nisincLYS300.hbonds.neo,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues LYS34 ")
dev.off()

################################ nisin ring E ################################################

ringEhbondtable50sernisin <- read.delim("hbondsringE", sep="", stringsAsFactors=FALSE, header = FALSE)

ringEhbondtable51pronisin <- read.delim("hbondsringEpro", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(ringEhbondtable50sernisin) <- colnames(hbondtable50pro)
colnames(ringEhbondtable51pronisin) <- colnames(hbondtable50pro)

ringEhbondtable50sernisin$Mutation <- "SERINE"
ringEhbondtable51pronisin$Mutation <- "PROLINE"

serpro51.nisinRingE.hbonds <- rbind(ringEhbondtable50sernisin, ringEhbondtable51pronisin) 
serpro51.nisinRingE.hbonds.temp <- within(serpro51.nisinRingE.hbonds, x <- paste(X.Acceptor,DonorH,sep='-'))

serpro51.nisinRingE.hbonds.neo <- serpro51.nisinRingE.hbonds.temp[,c("x", "Frac", "Mutation")]


serpro51.nisinRingE.hbonds.neo$x <- factor(serpro51.nisinRingE.hbonds.neo$x, levels = serpro51.nisinRingE.hbonds.neo$x)

pdf(file="serpronisinRingEHBONDS51nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro51.nisinRingE.hbonds.neo,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues DBB25-CYS26-ASN27")
dev.off()


################################ nisin ring D ################################################

ringDhbondtable50sernisin <- read.delim("hbondsringDser", sep="", stringsAsFactors=FALSE, header = FALSE)

ringDhbondtable51pronisin <- read.delim("hbondsringDpro", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(ringDhbondtable50sernisin) <- colnames(hbondtable50pro)
colnames(ringDhbondtable51pronisin) <- colnames(hbondtable51pro)

ringDhbondtable50sernisin$Mutation <- "SERINE"
ringDhbondtable51pronisin$Mutation <- "PROLINE"

serpro51.nisinRingD.hbonds <- rbind(ringDhbondtable50sernisin, ringDhbondtable51pronisin) 
serpro51.nisinRingD.hbonds.temp <- within(serpro51.nisinRingD.hbonds, x <- paste(X.Acceptor,DonorH,sep='-'))

serpro51.nisinRingD.hbonds.neo <- serpro51.nisinRingD.hbonds.temp[,c("x", "Frac", "Mutation")]


serpro51.nisinRingD.hbonds.neo$x <- factor(serpro51.nisinRingD.hbonds.neo$x, levels = serpro51.nisinRingD.hbonds.neo$x)

pdf(file="serpronisinRingDHBONDS51nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro51.nisinRingD.hbonds.neo,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues LYS22-DBB23-ALA24")
dev.off()


########################################## proval 51 ##############################################

hbondtable51proval <- read.delim("All.UU.avg.all10proval.dat", sep="", stringsAsFactors=FALSE)

hbondtable50sernisin <- read.delim("hbondsnisinSER", sep="", stringsAsFactors=FALSE, header = FALSE)

hbondtable51provalnisin <- read.delim("hbondprovalsnisinPROVAL", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(hbondtable51provalnisin) <- colnames(hbondtable50pro)

colnames(hbondtable50sernisin) <- colnames(hbondtable50pro)

hbondtable50sernisin$Mutation <- "SERINE"
hbondtable51provalnisin$Mutation <- "PROLINE"

serproval.nisin.hbonds51 <- rbind(hbondtable50sernisin, hbondtable51provalnisin)

serproval.nisin.hbonds.temp51 <- within(serproval.nisin.hbonds51, x <- paste(X.Acceptor,DonorH,sep='-'))

serproval.nisin.hbonds.neo51 <- serproval.nisin.hbonds.temp51[,c("x", "Frac", "Mutation")]


serproval.nisin.hbonds.neo51$x <- factor(serproval.nisin.hbonds.neo51$x, levels = serproval.nisin.hbonds.neo51$x)

pdf(file="serprovalnisinHBONDS51nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproval.nisin.hbonds.neo51,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=2),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues")
dev.off()

library(ggplot2)


############################ nisin core ###############################

corehbondtable50sernisin <- read.delim("hbondcoreser", sep="", stringsAsFactors=FALSE, header = FALSE)

corehbondtable51provalnisin <- read.delim("hbondcoreproval", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(corehbondtable51provalnisin) <- colnames(hbondtable50pro)
colnames(corehbondtable50sernisin) <- colnames(hbondtable50pro)

corehbondtable50sernisin$Mutation <- "SERINE"
corehbondtable51provalnisin$Mutation <- "PROLINE"

serproval51.nisincore.hbonds <- rbind(corehbondtable50sernisin, corehbondtable51provalnisin)

serproval51.nisincore.hbonds.temp <- within(serproval51.nisincore.hbonds, x <- paste(X.Acceptor,DonorH,sep='-'))

serproval51.nisincore.hbonds.neo <- serproval51.nisincore.hbonds.temp[,c("x", "Frac", "Mutation")]


serproval51.nisincore.hbonds.neo$x <- factor(serproval51.nisincore.hbonds.neo$x, levels = serproval51.nisincore.hbonds.neo$x)

pdf(file="serprovalnisincoreHBONDS51nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproval51.nisincore.hbonds.neo,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues CYS28-SER29(PRO29)-ILE30(VAL30")
dev.off()

########################################## alanine 50 ##############################################

hbondtable51proval <- read.delim("All.UU.avg.all10proval.dat", sep="", stringsAsFactors=FALSE)

hbondtable50sernisin <- read.delim("hbondsnisinSER", sep="", stringsAsFactors=FALSE, header = FALSE)

hbondtable51provalnisin <- read.delim("hbondprovalsnisinPROVAL", sep="", stringsAsFactors=FALSE, header = FALSE)


hbondtable50ala <- read.delim("hbondalasnisinALA", sep="", stringsAsFactors=FALSE)

colnames(hbondtable50ala) <- colnames(hbondtable50pro)

colnames(hbondtable50sernisin) <- colnames(hbondtable50pro)

colnames(hbondtable51provalnisin) <- colnames(hbondtable50pro)



hbondtable50sernisin$Mutation <- "SERINE"
hbondtable51provalnisin$Mutation <- "PROLINE"
hbondtable50ala$Mutation <- "ALANINE"


ala.nisin.hbonds51 <- rbind(hbondtable50sernisin, hbondtable51provalnisin, hbondtable50ala)

ala.nisin.hbonds.temp51 <- within(ala.nisin.hbonds51, x <- paste(X.Acceptor,DonorH,sep='-'))

ala.nisin.hbonds.neo51 <- ala.nisin.hbonds.temp51[,c("x", "Frac", "Mutation")]


ala.nisin.hbonds.neo51$x <- factor(ala.nisin.hbonds.neo51$x, levels = ala.nisin.hbonds.neo51$x)

pdf(file="alanisinHBONDS51nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(ala.nisin.hbonds.neo51,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=2),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues")
dev.off()

library(ggplot2)


############################ nisin core ###############################

corehbondtable50sernisin <- read.delim("hbondcoreser", sep="", stringsAsFactors=FALSE, header = FALSE)

corehbondtable51provalnisin <- read.delim("hbondcoreproval", sep="", stringsAsFactors=FALSE, header = FALSE)

corehbondtable51ala <- read.delim("hbondalacore", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(corehbondtable51provalnisin) <- colnames(hbondtable50pro)
colnames(corehbondtable50sernisin) <- colnames(hbondtable50pro)
colnames(corehbondtable51ala) <- colnames(hbondtable50pro)


corehbondtable50sernisin$Mutation <- "SERINE"
corehbondtable51provalnisin$Mutation <- "PROLINE"
corehbondtable51ala$Mutation <- "ALANINE"

ala.nisincore.hbonds <- rbind(corehbondtable50sernisin, corehbondtable51provalnisin, corehbondtable51ala)

ala.nisincore.hbonds.temp <- within(ala.nisincore.hbonds, x <- paste(X.Acceptor,DonorH,sep='-'))

ala.nisincore.hbonds.neo <- ala.nisincore.hbonds.temp[,c("x", "Frac", "Mutation")]


ala.nisincore.hbonds.neo$x <- factor(ala.nisincore.hbonds.neo$x, levels = ala.nisincore.hbonds.neo$x)

pdf(file="alanisincoreHBONDS51nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(ala.nisincore.hbonds.neo,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues CYS28-SER29(PRO29)-ILE30(VAL30")
dev.off()

############################ nisin C-Terminus ###############################

ctermhbondtable50sernisin <- read.delim("hbondCTERMINUSser", sep="", stringsAsFactors=FALSE, header = FALSE)

ctermhbondtable51proval <- read.delim("hbondprovalCTERMINUS", sep="", stringsAsFactors=FALSE, header = FALSE)

ctermhbondtable50ala <- read.delim("hbondalaCTERMINUS", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(ctermhbondtable50sernisin) <- colnames(hbondtable50pro)
colnames(ctermhbondtable51proval) <- colnames(hbondtable51pro)
colnames(ctermhbondtable50ala) <- colnames(hbondtable51pro)

ctermhbondtable50sernisin$Mutation <- "SERINE"
ctermhbondtable51proval$Mutation <- "PROLINE"
ctermhbondtable50ala$Mutation <- "ALANINE"

ala.nisincterm.hbonds <- rbind(ctermhbondtable50sernisin, ctermhbondtable51proval, ctermhbondtable50ala)

ala.nisincterm.hbonds.temp <- within(ala.nisincterm.hbonds, x <- paste(X.Acceptor,DonorH,sep='-'))

ala.nisincterm.hbonds.neo <- ala.nisincterm.hbonds.temp[,c("x", "Frac", "Mutation")]


ala.nisincterm.hbonds.neo$x <- factor(ala.nisincterm.hbonds.neo$x, levels = ala.nisincterm.hbonds.neo$x)

pdf(file="alaCTERMINUSHBONDS51nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(ala.nisincterm.hbonds.neo,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues HIE31-VAL32-DHA33")
dev.off()

################################ nisin Lyseine 300 ################################################

cLYS300hbondtable50sernisin <- read.delim("hbondLYS300ser", sep="", stringsAsFactors=FALSE, header = FALSE)

cLYS300hbondtable50provalnisin <- read.delim("hbondprovalLYS300", sep="", stringsAsFactors=FALSE, header = FALSE)

cLYS300hbondtable50ala <- read.delim("hbondalaLYS300", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(cLYS300hbondtable50sernisin) <- colnames(hbondtable50pro)
colnames(cLYS300hbondtable50provalnisin) <- colnames(hbondtable50pro)
colnames(cLYS300hbondtable50ala) <- colnames(hbondtable50pro)

cLYS300hbondtable50sernisin$Mutation <- "SERINE"
cLYS300hbondtable50provalnisin$Mutation <- "PROLINE"
cLYS300hbondtable50ala$Mutation <- "ALANINE"

ala.nisincLYS300.hbonds <- rbind(cLYS300hbondtable50sernisin, cLYS300hbondtable50provalnisin, cLYS300hbondtable50ala) 
ala.nisincLYS300.hbonds.temp <- within(ala.nisincLYS300.hbonds, x <- paste(X.Acceptor,DonorH,sep='-'))

ala.nisincLYS300.hbonds.neo <- ala.nisincLYS300.hbonds.temp[,c("x", "Frac", "Mutation")]


ala.nisincLYS300.hbonds.neo$x <- factor(ala.nisincLYS300.hbonds.neo$x, levels = ala.nisincLYS300.hbonds.neo$x)

pdf(file="alaCLYS300HBONDS50nanos.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serpro.nisincLYS300.hbonds.neo,aes(x = x,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  theme(axis.text=element_text(angle = 90, hjust = 1,size=8),
        axis.title=element_text(size=10,face="bold"),legend.text=element_text(size = 10), legend.title=element_text(size=10,face="bold")) +
  ggtitle("Hbonds of Residues LYS34 ")
dev.off()






############################ nisin core Text ###############################

corehbondtable50sernisin <- read.delim("hbondcoreser", sep="", stringsAsFactors=FALSE, header = FALSE)

corehbondtable51provalnisin <- read.delim("hbondprovalcore", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(corehbondtable51provalnisin) <- colnames(hbondtable50pro)
colnames(corehbondtable50sernisin) <- colnames(hbondtable50pro)

corehbondtable50sernisin$Mutation <- "Nisin A"
corehbondtable51provalnisin$Mutation <- "Nisin PV"

serproval51.nisincore.hbonds <- rbind(corehbondtable50sernisin, corehbondtable51provalnisin)

serproval51.nisincore.hbonds.temp <- within(serproval51.nisincore.hbonds, bond <- paste(X.Acceptor,DonorH,sep='-'))

serproval51.nisincore.hbonds.neo <- serproval51.nisincore.hbonds.temp[,c("bond", "Frac", "Mutation")]


serproval51.nisincore.hbonds.neo$bond <- factor(serproval51.nisincore.hbonds.neo$bond, levels = serproval51.nisincore.hbonds.neo$bond)

main4=expression("Hbond Occupancy (out of 1)")

#Hbonds of Residues CYS28-SER29(PRO29)-ILE30(VAL30
pdf(file="serprovalnisincoreHBONDS51nanosText.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproval51.nisincore.hbonds.neo,aes(x = bond,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  labs(y=main4) +
  theme(axis.text.x=element_text(angle = 60, hjust = 1,size=8),axis.text.y=element_text(angle = 0, hjust = 1,size=30),
        axis.title=element_text(size=30,face="bold"),legend.text=element_text(size = 40), legend.title=element_text(size=40,face="bold")) +
  ggtitle("")
dev.off()



########################################## proval 51 Text ##############################################

hbondtable51proval <- read.delim("All.UU.avg.all10proval.dat", sep="", stringsAsFactors=FALSE)

hbondtable50sernisin <- read.delim("hbondsnisinSER2", sep="", stringsAsFactors=FALSE, header = FALSE)

hbondtable51provalnisin <- read.delim("hbondprovalsnisinPROVAL2", sep="", stringsAsFactors=FALSE, header = FALSE)

colnames(hbondtable51provalnisin) <- colnames(hbondtable50pro)

colnames(hbondtable50sernisin) <- colnames(hbondtable50pro)

hbondtable50sernisin$Mutation <- "Nisin A"
hbondtable51provalnisin$Mutation <- "Nisin PV"

serproval.nisin.hbonds51 <- rbind(hbondtable50sernisin, hbondtable51provalnisin)

serproval.nisin.hbonds.temp51 <- within(serproval.nisin.hbonds51, bond <- paste(X.Acceptor,DonorH,sep='-'))

serproval.nisin.hbonds.neo51 <- serproval.nisin.hbonds.temp51[,c("bond", "Frac", "Mutation")]


serproval.nisin.hbonds.neo51$bond <- factor(serproval.nisin.hbonds.neo51$bond, levels = serproval.nisin.hbonds.neo51$bond)

#Hbonds of Residues

pdf(file="serprovalnisinHBONDS51nanosText.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproval.nisin.hbonds.neo51,aes(x = bond,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  labs(y=main4) +
  theme(axis.text.x=element_text(angle = 60, hjust = 1,size=5),axis.text.y=element_text(angle = 0, hjust = 1,size=30),
        axis.title=element_text(size=30,face="bold"),legend.text=element_text(size = 40), legend.title=element_text(size=40,face="bold")) +
  ggtitle("")
dev.off()

serproval.nisin.hbonds.neo51$Frac > 0.1
serproval.nisin.hbonds.neo51.greaterten <- serproval.nisin.hbonds.neo51[serproval.nisin.hbonds.neo51$Frac > 0.1,]
serproval.nisin.hbonds.neo51.greaterten$bond <- factor(serproval.nisin.hbonds.neo51.greaterten$bond, levels = serproval.nisin.hbonds.neo51.greaterten$bond)



main4=expression('Hbond Occupancy')



pdf(file="serprovalnisinHBONDS51nanosText10.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproval.nisin.hbonds.neo51.greaterten,aes(x = bond,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  labs(y=main4) +
  theme(axis.text.x=element_text(angle = 60, hjust = 1,size=10),axis.text.y=element_text(angle = 0, hjust = 1,size=40),
        axis.title=element_text(size=40,face="bold"),legend.text=element_text(size = 40), legend.title=element_text(size=40,face="bold")) +
  ggtitle("")
dev.off()


serproval.nisin.hbonds.neo51$Frac > 0.2
serproval.nisin.hbonds.neo51.greatertwenty <- serproval.nisin.hbonds.neo51[serproval.nisin.hbonds.neo51$Frac > 0.2,]
serproval.nisin.hbonds.neo51.greatertwenty$bond <- factor(serproval.nisin.hbonds.neo51.greatertwenty$bond, levels = serproval.nisin.hbonds.neo51.greatertwenty$bond)



main4=expression('Occupancy')



pdf(file="serprovalnisinHBONDS51nanosText20.pdf",paper="A4r",width=16,height=8.5,onefile = FALSE,useDingbats = FALSE)
ggplot(serproval.nisin.hbonds.neo51.greatertwenty,aes(x = bond,y = Frac)) + 
  geom_bar(aes(fill = Mutation),width=0.4, position = position_dodge(width=0.5),stat="identity") +
  labs(y=main4) +
  theme(axis.text.x=element_text(angle = 60, hjust = 1,size=20),axis.text.y=element_text(angle = 0, hjust = 1,size=40),
        axis.title=element_text(size=40,face="bold"),legend.text=element_text(size = 40), legend.title=element_text(size=40,face="bold")) +
  ggtitle("")
dev.off()