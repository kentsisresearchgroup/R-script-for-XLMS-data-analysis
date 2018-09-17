rm(list=ls())
ls()
getwd()
setwd("Z:\\04_projects\\XLMS_proteomics\\BSA_XLMS\\20170713_xlBSA_replicates\\20170830_Rcleanscript_xlBSAreplicates_fullanalysis")
library(dplyr)
library(plyr)
library(VennDiagram)
library(gridExtra)

#Score cutoffs
s1=0  #ETDMS2 humandb
s2=0  #HCDMS2 humandb
s3=0  #EThcDMS2 humandb
s4=0  #HCDMS3 humandb
s5=0  #ETDMS2_hek humandb
s6=0  #HCDMS2_hek humandb
s7=0  #EThcDMS2_hek humandb
s8=0  #HCDMS3_hek humandb
s9=0 #ETDMS2 bsadb
s10=0 #HCDMS2 bsadb
s11=0 #EThcDMS2 bsadb
s12=0 #HCDMS3 bsadb
s13=0 #ETDMS2_hek bsadb
s14=0 #HCDMS2_hek bsadb
s15=0 #EThcDMS2_hek bsadb
s16=0 #HCDMS3_hek bsadb

##ETD cross-linked BSA human database
ETD1_xlBSA_hdb=read.table("170629_xlBSA1_hdb_CIDMS2_ETDMS2_Crosslinks.txt", fill=TRUE, header=TRUE)
ETD2_xlBSA_hdb=read.table("170629_xlBSA2_hdb_CIDMS2_ETDMS2_Crosslinks.txt", fill=TRUE, header=TRUE)
ETD3_xlBSA_hdb=read.table("170629_xlBSA3_hdb_CIDMS2_ETDMS2_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
ETD1_xlBSA_hdb$Replicate=c(1)
ETD2_xlBSA_hdb$Replicate=c(2)
ETD3_xlBSA_hdb$Replicate=c(3)
##Compile all replicates
ETD_xlBSA_hdb=rbind(ETD1_xlBSA_hdb,ETD2_xlBSA_hdb,ETD3_xlBSA_hdb)
write.table(ETD_xlBSA_hdb,"20170830_ETDMS2_xlBSA_hdb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
ETD_xlBSA_hdb$PositionA=pmax(ETD_xlBSA_hdb$Position.A,ETD_xlBSA_hdb$Position.B)
ETD_xlBSA_hdb$PositionB=pmin(ETD_xlBSA_hdb$Position.A,ETD_xlBSA_hdb$Position.B)
##Keep BSA cross-links, duplicates kept
ETD_xlBSA_hdb_bsaonly<- ETD_xlBSA_hdb[grepl("P02769", ETD_xlBSA_hdb[["Accession.A"]]) & grepl("P02769", ETD_xlBSA_hdb[["Accession.B"]]), ]
write.table(ETD_xlBSA_hdb_bsaonly,"20170830_ETDMS2_xlBSA_hdb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
ETD_xlBSA_hdb_nonbsa<- ETD_xlBSA_hdb[!grepl("P02769", ETD_xlBSA_hdb[["Accession.A"]]) | !grepl("P02769", ETD_xlBSA_hdb[["Accession.B"]]), ]
write.table(ETD_xlBSA_hdb_nonbsa,"20170830_ETDMS2_xlBSA_hdb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
ETD_xlBSA_hdb_dup=ETD_xlBSA_hdb_bsaonly[!duplicated(ETD_xlBSA_hdb_bsaonly[c("PositionA","PositionB")]),]
write.table(ETD_xlBSA_hdb_dup,"20170830_ETDMS2_xlBSA_hdb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
ETD_xlBSA_hdb_dup$Accession.A=as.character(ETD_xlBSA_hdb_dup$Accession.A)
ETD_xlBSA_hdb_dup$Accession.A<-replace(ETD_xlBSA_hdb_dup$Accession.A,ETD_xlBSA_hdb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
ETD_xlBSA_hdb_dup$Accession.B=as.character(ETD_xlBSA_hdb_dup$Accession.B)
ETD_xlBSA_hdb_dup$Accession.B<-replace(ETD_xlBSA_hdb_dup$Accession.B,ETD_xlBSA_hdb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
ETD_xlBSA_hdb_xinet<-data.frame(ETD_xlBSA_hdb_dup$Max.XlinkX.Score,as.factor(ETD_xlBSA_hdb_dup$Accession.A),as.factor(ETD_xlBSA_hdb_dup$Accession.B),
ETD_xlBSA_hdb_dup$PositionA,ETD_xlBSA_hdb_dup$PositionB)
colnames(ETD_xlBSA_hdb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(ETD_xlBSA_hdb_xinet,"20170830_ETDMS2_xlBSA_hdb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
ETD1_xlBSA_hdb$PositionA=pmax(ETD1_xlBSA_hdb$Position.A,ETD1_xlBSA_hdb$Position.B)
ETD1_xlBSA_hdb$PositionB=pmin(ETD1_xlBSA_hdb$Position.A,ETD1_xlBSA_hdb$Position.B)
ETD2_xlBSA_hdb$PositionA=pmax(ETD2_xlBSA_hdb$Position.A,ETD2_xlBSA_hdb$Position.B)
ETD2_xlBSA_hdb$PositionB=pmin(ETD2_xlBSA_hdb$Position.A,ETD2_xlBSA_hdb$Position.B)
ETD3_xlBSA_hdb$PositionA=pmax(ETD3_xlBSA_hdb$Position.A,ETD3_xlBSA_hdb$Position.B)
ETD3_xlBSA_hdb$PositionB=pmin(ETD3_xlBSA_hdb$Position.A,ETD3_xlBSA_hdb$Position.B)
#Filter XlinkX Score >80
ETD1_xlBSA_hdb=ETD1_xlBSA_hdb[ETD1_xlBSA_hdb[,"Max.XlinkX.Score"]>s1,]
ETD2_xlBSA_hdb=ETD2_xlBSA_hdb[ETD2_xlBSA_hdb[,"Max.XlinkX.Score"]>s1,]
ETD3_xlBSA_hdb=ETD3_xlBSA_hdb[ETD3_xlBSA_hdb[,"Max.XlinkX.Score"]>s1,]
#Remove Duplicates
ETD1_xlBSA_hdb_dup=ETD1_xlBSA_hdb[!duplicated(ETD1_xlBSA_hdb[c("PositionA","PositionB")]),]
ETD2_xlBSA_hdb_dup=ETD2_xlBSA_hdb[!duplicated(ETD2_xlBSA_hdb[c("PositionA","PositionB")]),]
ETD3_xlBSA_hdb_dup=ETD3_xlBSA_hdb[!duplicated(ETD3_xlBSA_hdb[c("PositionA","PositionB")]),]
#ETD Number of BSA crosslinks
ETD1_xlBSA_hdb_bsaonly<- ETD1_xlBSA_hdb_dup[grepl("P02769", ETD1_xlBSA_hdb_dup[["Accession.A"]]) & grepl("P02769", ETD1_xlBSA_hdb_dup[["Accession.B"]]), ]
ETD2_xlBSA_hdb_bsaonly<- ETD2_xlBSA_hdb_dup[grepl("P02769", ETD2_xlBSA_hdb_dup[["Accession.A"]]) & grepl("P02769", ETD2_xlBSA_hdb_dup[["Accession.B"]]), ]
ETD3_xlBSA_hdb_bsaonly<- ETD3_xlBSA_hdb_dup[grepl("P02769", ETD3_xlBSA_hdb_dup[["Accession.A"]]) & grepl("P02769", ETD3_xlBSA_hdb_dup[["Accession.B"]]), ]
#ETD Number of all crosslinks, no duplicates
ETD_xlBSA_hdb_crosslinks=c(nrow(ETD1_xlBSA_hdb_dup),nrow(ETD2_xlBSA_hdb_dup),nrow(ETD3_xlBSA_hdb_dup))
#ETD Number of BSA crosslinks, no duplicates
ETD_xlBSA_hdb_bsaonly=c(nrow(ETD1_xlBSA_hdb_bsaonly),nrow(ETD2_xlBSA_hdb_bsaonly),nrow(ETD3_xlBSA_hdb_bsaonly))
#ETD Number of non-BSA crosslinks
ETD_xlBSA_hdb_nonbsa=(ETD_xlBSA_hdb_crosslinks-ETD_xlBSA_hdb_bsaonly)
#Output for Number of crosslinks
ETD_xlBSA_hdb=data.frame(ETD_xlBSA_hdb_crosslinks,ETD_xlBSA_hdb_bsaonly,ETD_xlBSA_hdb_nonbsa)
colnames(ETD_xlBSA_hdb)=c("ETD_xlBSA_hdb_Total","ETD_xlBSA_hdb_bsacrosslinks","ETD_xlBSA_hdb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(ETD1_xlBSA_hdb_bsaonly,ETD2_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(ETD1_xlBSA_hdb_bsaonly,ETD3_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(ETD2_xlBSA_hdb_bsaonly,ETD3_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB")))
tiff('ETD_xlBSA_hdb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(ETD1_xlBSA_hdb_bsaonly),area2=nrow(ETD2_xlBSA_hdb_bsaonly),area3=nrow(ETD3_xlBSA_hdb_bsaonly),
n12=nrow(merge(ETD1_xlBSA_hdb_bsaonly,ETD2_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(ETD1_xlBSA_hdb_bsaonly,ETD3_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(ETD2_xlBSA_hdb_bsaonly,ETD3_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(ETD_xlBSA_hdb_dup)-(nrow(ETD1_xlBSA_hdb_bsaonly)+nrow(ETD2_xlBSA_hdb_bsaonly)+nrow(ETD3_xlBSA_hdb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("ETD_xlBSA_hdb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()

##HCD cross-linked BSA human database
HCD1_xlBSA_hdb=read.table("170629_xlBSA1_hdb_CIDMS2_HCDMS2_Crosslinks.txt", fill=TRUE, header=TRUE)
HCD2_xlBSA_hdb=read.table("170629_xlBSA2_hdb_CIDMS2_HCDMS2_Crosslinks.txt", fill=TRUE, header=TRUE)
HCD3_xlBSA_hdb=read.table("170629_xlBSA3_hdb_CIDMS2_HCDMS2_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
HCD1_xlBSA_hdb$Replicate=c(1)
HCD2_xlBSA_hdb$Replicate=c(2)
HCD3_xlBSA_hdb$Replicate=c(3)
##Compile all replicates
HCD_xlBSA_hdb=rbind(HCD1_xlBSA_hdb,HCD2_xlBSA_hdb,HCD3_xlBSA_hdb)
write.table(HCD_xlBSA_hdb,"20170830_HCDMS2_xlBSA_hdb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
HCD_xlBSA_hdb$PositionA=pmax(HCD_xlBSA_hdb$Position.A,HCD_xlBSA_hdb$Position.B)
HCD_xlBSA_hdb$PositionB=pmin(HCD_xlBSA_hdb$Position.A,HCD_xlBSA_hdb$Position.B)
##Keep BSA cross-links, duplicates kept
HCD_xlBSA_hdb_bsaonly<- HCD_xlBSA_hdb[grepl("P02769", HCD_xlBSA_hdb[["Accession.A"]]) & grepl("P02769", HCD_xlBSA_hdb[["Accession.B"]]), ]
write.table(HCD_xlBSA_hdb_bsaonly,"20170830_HCDMS2_xlBSA_hdb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
HCD_xlBSA_hdb_nonbsa<- HCD_xlBSA_hdb[!grepl("P02769", HCD_xlBSA_hdb[["Accession.A"]]) | !grepl("P02769", HCD_xlBSA_hdb[["Accession.B"]]), ]
write.table(HCD_xlBSA_hdb_nonbsa,"20170830_HCDMS2_xlBSA_hdb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
HCD_xlBSA_hdb_dup=HCD_xlBSA_hdb_bsaonly[!duplicated(HCD_xlBSA_hdb_bsaonly[c("PositionA","PositionB")]),]
write.table(HCD_xlBSA_hdb_dup,"20170830_HCDMS2_xlBSA_hdb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
HCD_xlBSA_hdb_dup$Accession.A=as.character(HCD_xlBSA_hdb_dup$Accession.A)
HCD_xlBSA_hdb_dup$Accession.A<-replace(HCD_xlBSA_hdb_dup$Accession.A,HCD_xlBSA_hdb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
HCD_xlBSA_hdb_dup$Accession.B=as.character(HCD_xlBSA_hdb_dup$Accession.B)
HCD_xlBSA_hdb_dup$Accession.B<-replace(HCD_xlBSA_hdb_dup$Accession.B,HCD_xlBSA_hdb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
HCD_xlBSA_hdb_xinet<-data.frame(HCD_xlBSA_hdb_dup$Max.XlinkX.Score,as.factor(HCD_xlBSA_hdb_dup$Accession.A),as.factor(HCD_xlBSA_hdb_dup$Accession.B),
HCD_xlBSA_hdb_dup$PositionA,HCD_xlBSA_hdb_dup$PositionB)
colnames(HCD_xlBSA_hdb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(HCD_xlBSA_hdb_xinet,"20170830_HCDMS2_xlBSA_hdb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
HCD1_xlBSA_hdb$PositionA=pmax(HCD1_xlBSA_hdb$Position.A,HCD1_xlBSA_hdb$Position.B)
HCD1_xlBSA_hdb$PositionB=pmin(HCD1_xlBSA_hdb$Position.A,HCD1_xlBSA_hdb$Position.B)
HCD2_xlBSA_hdb$PositionA=pmax(HCD2_xlBSA_hdb$Position.A,HCD2_xlBSA_hdb$Position.B)
HCD2_xlBSA_hdb$PositionB=pmin(HCD2_xlBSA_hdb$Position.A,HCD2_xlBSA_hdb$Position.B)
HCD3_xlBSA_hdb$PositionA=pmax(HCD3_xlBSA_hdb$Position.A,HCD3_xlBSA_hdb$Position.B)
HCD3_xlBSA_hdb$PositionB=pmin(HCD3_xlBSA_hdb$Position.A,HCD3_xlBSA_hdb$Position.B)
#Filter XlinkX Score >80
HCD1_xlBSA_hdb=HCD1_xlBSA_hdb[HCD1_xlBSA_hdb[,"Max.XlinkX.Score"]>s2,]
HCD2_xlBSA_hdb=HCD2_xlBSA_hdb[HCD2_xlBSA_hdb[,"Max.XlinkX.Score"]>s2,]
HCD3_xlBSA_hdb=HCD3_xlBSA_hdb[HCD3_xlBSA_hdb[,"Max.XlinkX.Score"]>s2,]
#Remove Duplicates
HCD1_xlBSA_hdb_dup=HCD1_xlBSA_hdb[!duplicated(HCD1_xlBSA_hdb[c("PositionA","PositionB")]),]
HCD2_xlBSA_hdb_dup=HCD2_xlBSA_hdb[!duplicated(HCD2_xlBSA_hdb[c("PositionA","PositionB")]),]
HCD3_xlBSA_hdb_dup=HCD3_xlBSA_hdb[!duplicated(HCD3_xlBSA_hdb[c("PositionA","PositionB")]),]
#HCD Number of BSA crosslinks
HCD1_xlBSA_hdb_bsaonly<- HCD1_xlBSA_hdb_dup[grepl("P02769", HCD1_xlBSA_hdb_dup[["Accession.A"]]) & grepl("P02769", HCD1_xlBSA_hdb_dup[["Accession.B"]]), ]
HCD2_xlBSA_hdb_bsaonly<- HCD2_xlBSA_hdb_dup[grepl("P02769", HCD2_xlBSA_hdb_dup[["Accession.A"]]) & grepl("P02769", HCD2_xlBSA_hdb_dup[["Accession.B"]]), ]
HCD3_xlBSA_hdb_bsaonly<- HCD3_xlBSA_hdb_dup[grepl("P02769", HCD3_xlBSA_hdb_dup[["Accession.A"]]) & grepl("P02769", HCD3_xlBSA_hdb_dup[["Accession.B"]]), ]
#HCD Number of all crosslinks, no duplicates
HCD_xlBSA_hdb_crosslinks=c(nrow(HCD1_xlBSA_hdb_dup),nrow(HCD2_xlBSA_hdb_dup),nrow(HCD3_xlBSA_hdb_dup))
#HCD Number of BSA crosslinks, no duplicates
HCD_xlBSA_hdb_bsaonly=c(nrow(HCD1_xlBSA_hdb_bsaonly),nrow(HCD2_xlBSA_hdb_bsaonly),nrow(HCD3_xlBSA_hdb_bsaonly))
#HCD Number of non-BSA crosslinks
HCD_xlBSA_hdb_nonbsa=(HCD_xlBSA_hdb_crosslinks-HCD_xlBSA_hdb_bsaonly)
#Output for Number of crosslinks
HCD_xlBSA_hdb=data.frame(HCD_xlBSA_hdb_crosslinks,HCD_xlBSA_hdb_bsaonly,HCD_xlBSA_hdb_nonbsa)
colnames(HCD_xlBSA_hdb)=c("HCD_xlBSA_hdb_Total","HCD_xlBSA_hdb_bsacrosslinks","HCD_xlBSA_hdb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(HCD1_xlBSA_hdb_bsaonly,HCD2_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(HCD1_xlBSA_hdb_bsaonly,HCD3_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(HCD2_xlBSA_hdb_bsaonly,HCD3_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB")))
tiff('HCD_xlBSA_hdb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(HCD1_xlBSA_hdb_bsaonly),area2=nrow(HCD2_xlBSA_hdb_bsaonly),area3=nrow(HCD3_xlBSA_hdb_bsaonly),
n12=nrow(merge(HCD1_xlBSA_hdb_bsaonly,HCD2_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(HCD1_xlBSA_hdb_bsaonly,HCD3_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(HCD2_xlBSA_hdb_bsaonly,HCD3_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(HCD_xlBSA_hdb_dup)-(nrow(HCD1_xlBSA_hdb_bsaonly)+nrow(HCD2_xlBSA_hdb_bsaonly)+nrow(HCD3_xlBSA_hdb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("HCD_xlBSA_hdb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()

##ETHCD cross-linked BSA+hek293 human database
ETHCD1_xlBSA_hdb=read.table("170629_xlBSA1_hdb_CIDMS2_EThcDMS2_Crosslinks.txt", fill=TRUE, header=TRUE)
ETHCD2_xlBSA_hdb=read.table("170629_xlBSA2_hdb_CIDMS2_EThcDMS2_Crosslinks.txt", fill=TRUE, header=TRUE)
ETHCD3_xlBSA_hdb=read.table("170629_xlBSA3_hdb_CIDMS2_EThcDMS2_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
ETHCD1_xlBSA_hdb$Replicate=c(1)
ETHCD2_xlBSA_hdb$Replicate=c(2)
ETHCD3_xlBSA_hdb$Replicate=c(3)
##Compile all replicates
ETHCD_xlBSA_hdb=rbind(ETHCD1_xlBSA_hdb,ETHCD2_xlBSA_hdb,ETHCD3_xlBSA_hdb)
write.table(ETHCD_xlBSA_hdb,"20170830_ETHCDMS2_xlBSA_hdb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
ETHCD_xlBSA_hdb$PositionA=pmax(ETHCD_xlBSA_hdb$Position.A,ETHCD_xlBSA_hdb$Position.B)
ETHCD_xlBSA_hdb$PositionB=pmin(ETHCD_xlBSA_hdb$Position.A,ETHCD_xlBSA_hdb$Position.B)
##Keep BSA cross-links, duplicates kept
ETHCD_xlBSA_hdb_bsaonly<- ETHCD_xlBSA_hdb[grepl("P02769", ETHCD_xlBSA_hdb[["Accession.A"]]) & grepl("P02769", ETHCD_xlBSA_hdb[["Accession.B"]]), ]
write.table(ETHCD_xlBSA_hdb_bsaonly,"20170830_ETHCDMS2_xlBSA_hdb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
ETHCD_xlBSA_hdb_nonbsa<- ETHCD_xlBSA_hdb[!grepl("P02769", ETHCD_xlBSA_hdb[["Accession.A"]]) | !grepl("P02769", ETHCD_xlBSA_hdb[["Accession.B"]]), ]
write.table(ETHCD_xlBSA_hdb_nonbsa,"20170830_ETHCDMS2_xlBSA_hdb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
ETHCD_xlBSA_hdb_dup=ETHCD_xlBSA_hdb_bsaonly[!duplicated(ETHCD_xlBSA_hdb_bsaonly[c("PositionA","PositionB")]),]
write.table(ETHCD_xlBSA_hdb_dup,"20170830_ETHCDMS2_xlBSA_hdb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
ETHCD_xlBSA_hdb_dup$Accession.A=as.character(ETHCD_xlBSA_hdb_dup$Accession.A)
ETHCD_xlBSA_hdb_dup$Accession.A<-replace(ETHCD_xlBSA_hdb_dup$Accession.A,ETHCD_xlBSA_hdb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
ETHCD_xlBSA_hdb_dup$Accession.B=as.character(ETHCD_xlBSA_hdb_dup$Accession.B)
ETHCD_xlBSA_hdb_dup$Accession.B<-replace(ETHCD_xlBSA_hdb_dup$Accession.B,ETHCD_xlBSA_hdb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
ETHCD_xlBSA_hdb_xinet<-data.frame(ETHCD_xlBSA_hdb_dup$Max.XlinkX.Score,as.factor(ETHCD_xlBSA_hdb_dup$Accession.A),as.factor(ETHCD_xlBSA_hdb_dup$Accession.B),
ETHCD_xlBSA_hdb_dup$PositionA,ETHCD_xlBSA_hdb_dup$PositionB)
colnames(ETHCD_xlBSA_hdb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(ETHCD_xlBSA_hdb_xinet,"20170830_ETHCDMS2_xlBSA_hdb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
ETHCD1_xlBSA_hdb$PositionA=pmax(ETHCD1_xlBSA_hdb$Position.A,ETHCD1_xlBSA_hdb$Position.B)
ETHCD1_xlBSA_hdb$PositionB=pmin(ETHCD1_xlBSA_hdb$Position.A,ETHCD1_xlBSA_hdb$Position.B)
ETHCD2_xlBSA_hdb$PositionA=pmax(ETHCD2_xlBSA_hdb$Position.A,ETHCD2_xlBSA_hdb$Position.B)
ETHCD2_xlBSA_hdb$PositionB=pmin(ETHCD2_xlBSA_hdb$Position.A,ETHCD2_xlBSA_hdb$Position.B)
ETHCD3_xlBSA_hdb$PositionA=pmax(ETHCD3_xlBSA_hdb$Position.A,ETHCD3_xlBSA_hdb$Position.B)
ETHCD3_xlBSA_hdb$PositionB=pmin(ETHCD3_xlBSA_hdb$Position.A,ETHCD3_xlBSA_hdb$Position.B)
#Filter XlinkX Score >80
ETHCD1_xlBSA_hdb=ETHCD1_xlBSA_hdb[ETHCD1_xlBSA_hdb[,"Max.XlinkX.Score"]>s3,]
ETHCD2_xlBSA_hdb=ETHCD2_xlBSA_hdb[ETHCD2_xlBSA_hdb[,"Max.XlinkX.Score"]>s3,]
ETHCD3_xlBSA_hdb=ETHCD3_xlBSA_hdb[ETHCD3_xlBSA_hdb[,"Max.XlinkX.Score"]>s3,]
#Remove Duplicates
ETHCD1_xlBSA_hdb_dup=ETHCD1_xlBSA_hdb[!duplicated(ETHCD1_xlBSA_hdb[c("PositionA","PositionB")]),]
ETHCD2_xlBSA_hdb_dup=ETHCD2_xlBSA_hdb[!duplicated(ETHCD2_xlBSA_hdb[c("PositionA","PositionB")]),]
ETHCD3_xlBSA_hdb_dup=ETHCD3_xlBSA_hdb[!duplicated(ETHCD3_xlBSA_hdb[c("PositionA","PositionB")]),]
#ETHCD Number of BSA crosslinks
ETHCD1_xlBSA_hdb_bsaonly<- ETHCD1_xlBSA_hdb_dup[grepl("P02769", ETHCD1_xlBSA_hdb_dup[["Accession.A"]]) & grepl("P02769", ETHCD1_xlBSA_hdb_dup[["Accession.B"]]), ]
ETHCD2_xlBSA_hdb_bsaonly<- ETHCD2_xlBSA_hdb_dup[grepl("P02769", ETHCD2_xlBSA_hdb_dup[["Accession.A"]]) & grepl("P02769", ETHCD2_xlBSA_hdb_dup[["Accession.B"]]), ]
ETHCD3_xlBSA_hdb_bsaonly<- ETHCD3_xlBSA_hdb_dup[grepl("P02769", ETHCD3_xlBSA_hdb_dup[["Accession.A"]]) & grepl("P02769", ETHCD3_xlBSA_hdb_dup[["Accession.B"]]), ]
#ETHCD Number of all crosslinks, no duplicates
ETHCD_xlBSA_hdb_crosslinks=c(nrow(ETHCD1_xlBSA_hdb_dup),nrow(ETHCD2_xlBSA_hdb_dup),nrow(ETHCD3_xlBSA_hdb_dup))
#ETHCD Number of BSA crosslinks, no duplicates
ETHCD_xlBSA_hdb_bsaonly=c(nrow(ETHCD1_xlBSA_hdb_bsaonly),nrow(ETHCD2_xlBSA_hdb_bsaonly),nrow(ETHCD3_xlBSA_hdb_bsaonly))
#ETHCD Number of non-BSA crosslinks
ETHCD_xlBSA_hdb_nonbsa=(ETHCD_xlBSA_hdb_crosslinks-ETHCD_xlBSA_hdb_bsaonly)
#Output for Number of crosslinks
ETHCD_xlBSA_hdb=data.frame(ETHCD_xlBSA_hdb_crosslinks,ETHCD_xlBSA_hdb_bsaonly,ETHCD_xlBSA_hdb_nonbsa)
colnames(ETHCD_xlBSA_hdb)=c("ETHCD_xlBSA_hdb_Total","ETHCD_xlBSA_hdb_bsacrosslinks","ETHCD_xlBSA_hdb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(ETHCD1_xlBSA_hdb_bsaonly,ETHCD2_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(ETHCD1_xlBSA_hdb_bsaonly,ETHCD3_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(ETHCD2_xlBSA_hdb_bsaonly,ETHCD3_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB")))
tiff('ETHCD_xlBSA_hdb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(ETHCD1_xlBSA_hdb_bsaonly),area2=nrow(ETHCD2_xlBSA_hdb_bsaonly),area3=nrow(ETHCD3_xlBSA_hdb_bsaonly),
n12=nrow(merge(ETHCD1_xlBSA_hdb_bsaonly,ETHCD2_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(ETHCD1_xlBSA_hdb_bsaonly,ETHCD3_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(ETHCD2_xlBSA_hdb_bsaonly,ETHCD3_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(ETHCD_xlBSA_hdb_dup)-(nrow(ETHCD1_xlBSA_hdb_bsaonly)+nrow(ETHCD2_xlBSA_hdb_bsaonly)+nrow(ETHCD3_xlBSA_hdb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("ETHCD_xlBSA_hdb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()

##HCDMS3 cross-linked BSA human database
HCDMS31_xlBSA_hdb=read.table("170629_xlBSA1_hdb_CIDMS2_HCDMS3_Crosslinks.txt", fill=TRUE, header=TRUE)
HCDMS32_xlBSA_hdb=read.table("170629_xlBSA2_hdb_CIDMS2_HCDMS3_Crosslinks.txt", fill=TRUE, header=TRUE)
HCDMS33_xlBSA_hdb=read.table("170629_xlBSA3_hdb_CIDMS2_HCDMS3_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
HCDMS31_xlBSA_hdb$Replicate=c(1)
HCDMS32_xlBSA_hdb$Replicate=c(2)
HCDMS33_xlBSA_hdb$Replicate=c(3)
##Compile all replicates
HCDMS3_xlBSA_hdb=rbind(HCDMS31_xlBSA_hdb,HCDMS32_xlBSA_hdb,HCDMS33_xlBSA_hdb)
write.table(HCDMS3_xlBSA_hdb,"20170830_HCDMS3_xlBSA_hdb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
HCDMS3_xlBSA_hdb$PositionA=pmax(HCDMS3_xlBSA_hdb$Position.A,HCDMS3_xlBSA_hdb$Position.B)
HCDMS3_xlBSA_hdb$PositionB=pmin(HCDMS3_xlBSA_hdb$Position.A,HCDMS3_xlBSA_hdb$Position.B)
##Keep BSA cross-links, duplicates kept
HCDMS3_xlBSA_hdb_bsaonly<- HCDMS3_xlBSA_hdb[grepl("P02769", HCDMS3_xlBSA_hdb[["Accession.A"]]) & grepl("P02769", HCDMS3_xlBSA_hdb[["Accession.B"]]), ]
write.table(HCDMS3_xlBSA_hdb_bsaonly,"20170830_HCDMS3_xlBSA_hdb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
HCDMS3_xlBSA_hdb_nonbsa<- HCDMS3_xlBSA_hdb[!grepl("P02769", HCDMS3_xlBSA_hdb[["Accession.A"]]) | !grepl("P02769", HCDMS3_xlBSA_hdb[["Accession.B"]]), ]
write.table(HCDMS3_xlBSA_hdb_nonbsa,"20170830_HCDMS3_xlBSA_hdb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
HCDMS3_xlBSA_hdb_dup=HCDMS3_xlBSA_hdb_bsaonly[!duplicated(HCDMS3_xlBSA_hdb_bsaonly[c("PositionA","PositionB")]),]
write.table(HCDMS3_xlBSA_hdb_dup,"20170830_HCDMS3_xlBSA_hdb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
HCDMS3_xlBSA_hdb_dup$Accession.A=as.character(HCDMS3_xlBSA_hdb_dup$Accession.A)
HCDMS3_xlBSA_hdb_dup$Accession.A<-replace(HCDMS3_xlBSA_hdb_dup$Accession.A,HCDMS3_xlBSA_hdb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
HCDMS3_xlBSA_hdb_dup$Accession.B=as.character(HCDMS3_xlBSA_hdb_dup$Accession.B)
HCDMS3_xlBSA_hdb_dup$Accession.B<-replace(HCDMS3_xlBSA_hdb_dup$Accession.B,HCDMS3_xlBSA_hdb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
HCDMS3_xlBSA_hdb_xinet<-data.frame(HCDMS3_xlBSA_hdb_dup$Max.XlinkX.Score,as.factor(HCDMS3_xlBSA_hdb_dup$Accession.A),as.factor(HCDMS3_xlBSA_hdb_dup$Accession.B),
HCDMS3_xlBSA_hdb_dup$PositionA,HCDMS3_xlBSA_hdb_dup$PositionB)
colnames(HCDMS3_xlBSA_hdb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(HCDMS3_xlBSA_hdb_xinet,"20170830_HCDMS3_xlBSA_hdb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
HCDMS31_xlBSA_hdb$PositionA=pmax(HCDMS31_xlBSA_hdb$Position.A,HCDMS31_xlBSA_hdb$Position.B)
HCDMS31_xlBSA_hdb$PositionB=pmin(HCDMS31_xlBSA_hdb$Position.A,HCDMS31_xlBSA_hdb$Position.B)
HCDMS32_xlBSA_hdb$PositionA=pmax(HCDMS32_xlBSA_hdb$Position.A,HCDMS32_xlBSA_hdb$Position.B)
HCDMS32_xlBSA_hdb$PositionB=pmin(HCDMS32_xlBSA_hdb$Position.A,HCDMS32_xlBSA_hdb$Position.B)
HCDMS33_xlBSA_hdb$PositionA=pmax(HCDMS33_xlBSA_hdb$Position.A,HCDMS33_xlBSA_hdb$Position.B)
HCDMS33_xlBSA_hdb$PositionB=pmin(HCDMS33_xlBSA_hdb$Position.A,HCDMS33_xlBSA_hdb$Position.B)
#Filter XlinkX Score >80
HCDMS31_xlBSA_hdb=HCDMS31_xlBSA_hdb[HCDMS31_xlBSA_hdb[,"Max.XlinkX.Score"]>s4,]
HCDMS32_xlBSA_hdb=HCDMS32_xlBSA_hdb[HCDMS32_xlBSA_hdb[,"Max.XlinkX.Score"]>s4,]
HCDMS33_xlBSA_hdb=HCDMS33_xlBSA_hdb[HCDMS33_xlBSA_hdb[,"Max.XlinkX.Score"]>s4,]
#Remove Duplicates
HCDMS31_xlBSA_hdb_dup=HCDMS31_xlBSA_hdb[!duplicated(HCDMS31_xlBSA_hdb[c("PositionA","PositionB")]),]
HCDMS32_xlBSA_hdb_dup=HCDMS32_xlBSA_hdb[!duplicated(HCDMS32_xlBSA_hdb[c("PositionA","PositionB")]),]
HCDMS33_xlBSA_hdb_dup=HCDMS33_xlBSA_hdb[!duplicated(HCDMS33_xlBSA_hdb[c("PositionA","PositionB")]),]
#HCDMS3 Number of BSA crosslinks
HCDMS31_xlBSA_hdb_bsaonly<- HCDMS31_xlBSA_hdb_dup[grepl("P02769", HCDMS31_xlBSA_hdb_dup[["Accession.A"]]) & grepl("P02769", HCDMS31_xlBSA_hdb_dup[["Accession.B"]]), ]
HCDMS32_xlBSA_hdb_bsaonly<- HCDMS32_xlBSA_hdb_dup[grepl("P02769", HCDMS32_xlBSA_hdb_dup[["Accession.A"]]) & grepl("P02769", HCDMS32_xlBSA_hdb_dup[["Accession.B"]]), ]
HCDMS33_xlBSA_hdb_bsaonly<- HCDMS33_xlBSA_hdb_dup[grepl("P02769", HCDMS33_xlBSA_hdb_dup[["Accession.A"]]) & grepl("P02769", HCDMS33_xlBSA_hdb_dup[["Accession.B"]]), ]
#HCDMS3 Number of all crosslinks, no duplicates
HCDMS3_xlBSA_hdb_crosslinks=c(nrow(HCDMS31_xlBSA_hdb_dup),nrow(HCDMS32_xlBSA_hdb_dup),nrow(HCDMS33_xlBSA_hdb_dup))
#HCDMS3 Number of BSA crosslinks, no duplicates
HCDMS3_xlBSA_hdb_bsaonly=c(nrow(HCDMS31_xlBSA_hdb_bsaonly),nrow(HCDMS32_xlBSA_hdb_bsaonly),nrow(HCDMS33_xlBSA_hdb_bsaonly))
#HCDMS3 Number of non-BSA crosslinks
HCDMS3_xlBSA_hdb_nonbsa=(HCDMS3_xlBSA_hdb_crosslinks-HCDMS3_xlBSA_hdb_bsaonly)
#Output for Number of crosslinks
HCDMS3_xlBSA_hdb=data.frame(HCDMS3_xlBSA_hdb_crosslinks,HCDMS3_xlBSA_hdb_bsaonly,HCDMS3_xlBSA_hdb_nonbsa)
colnames(HCDMS3_xlBSA_hdb)=c("HCDMS3_xlBSA_hdb_Total","HCDMS3_xlBSA_hdb_bsacrosslinks","HCDMS3_xlBSA_hdb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(HCDMS31_xlBSA_hdb_bsaonly,HCDMS32_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(HCDMS31_xlBSA_hdb_bsaonly,HCDMS33_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(HCDMS32_xlBSA_hdb_bsaonly,HCDMS33_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB")))
tiff('HCDMS3_xlBSA_hdb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(HCDMS31_xlBSA_hdb_bsaonly),area2=nrow(HCDMS32_xlBSA_hdb_bsaonly),area3=nrow(HCDMS33_xlBSA_hdb_bsaonly),
n12=nrow(merge(HCDMS31_xlBSA_hdb_bsaonly,HCDMS32_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(HCDMS31_xlBSA_hdb_bsaonly,HCDMS33_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(HCDMS32_xlBSA_hdb_bsaonly,HCDMS33_xlBSA_hdb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(HCDMS3_xlBSA_hdb_dup)-(nrow(HCDMS31_xlBSA_hdb_bsaonly)+nrow(HCDMS32_xlBSA_hdb_bsaonly)+nrow(HCDMS33_xlBSA_hdb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("HCDMS3_xlBSA_hdb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()


##ETD cross-linked BSA+hek293t human database
ETD1_xlBSA_293_hdb=read.table("170629_xlBSA1_hdb_CIDMS2_ETDMS2_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
ETD2_xlBSA_293_hdb=read.table("170629_xlBSA2_hdb_CIDMS2_ETDMS2_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
ETD3_xlBSA_293_hdb=read.table("170629_xlBSA3_hdb_CIDMS2_ETDMS2_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
ETD1_xlBSA_293_hdb$Replicate=c(1)
ETD2_xlBSA_293_hdb$Replicate=c(2)
ETD3_xlBSA_293_hdb$Replicate=c(3)
##Compile all replicates
ETD_xlBSA_293_hdb=rbind(ETD1_xlBSA_293_hdb,ETD2_xlBSA_293_hdb,ETD3_xlBSA_293_hdb)
write.table(ETD_xlBSA_293_hdb,"20170830_ETDMS2_xlBSA_293_hdb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
ETD_xlBSA_293_hdb$PositionA=pmax(ETD_xlBSA_293_hdb$Position.A,ETD_xlBSA_293_hdb$Position.B)
ETD_xlBSA_293_hdb$PositionB=pmin(ETD_xlBSA_293_hdb$Position.A,ETD_xlBSA_293_hdb$Position.B)
##Keep BSA cross-links, duplicates kept
ETD_xlBSA_293_hdb_bsaonly<- ETD_xlBSA_293_hdb[grepl("P02769", ETD_xlBSA_293_hdb[["Accession.A"]]) & grepl("P02769", ETD_xlBSA_293_hdb[["Accession.B"]]), ]
write.table(ETD_xlBSA_293_hdb_bsaonly,"20170830_ETDMS2_xlBSA_293_hdb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
ETD_xlBSA_293_hdb_nonbsa<- ETD_xlBSA_293_hdb[!grepl("P02769", ETD_xlBSA_293_hdb[["Accession.A"]]) | !grepl("P02769", ETD_xlBSA_293_hdb[["Accession.B"]]), ]
write.table(ETD_xlBSA_293_hdb_nonbsa,"20170830_ETDMS2_xlBSA_293_hdb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
ETD_xlBSA_293_hdb_dup=ETD_xlBSA_293_hdb_bsaonly[!duplicated(ETD_xlBSA_293_hdb_bsaonly[c("PositionA","PositionB")]),]
write.table(ETD_xlBSA_293_hdb_dup,"20170830_ETDMS2_xlBSA_293_hdb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
ETD_xlBSA_293_hdb_dup$Accession.A=as.character(ETD_xlBSA_293_hdb_dup$Accession.A)
ETD_xlBSA_293_hdb_dup$Accession.A<-replace(ETD_xlBSA_293_hdb_dup$Accession.A,ETD_xlBSA_293_hdb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
ETD_xlBSA_293_hdb_dup$Accession.B=as.character(ETD_xlBSA_293_hdb_dup$Accession.B)
ETD_xlBSA_293_hdb_dup$Accession.B<-replace(ETD_xlBSA_293_hdb_dup$Accession.B,ETD_xlBSA_293_hdb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
ETD_xlBSA_293_hdb_xinet<-data.frame(ETD_xlBSA_293_hdb_dup$Max.XlinkX.Score,as.factor(ETD_xlBSA_293_hdb_dup$Accession.A),as.factor(ETD_xlBSA_293_hdb_dup$Accession.B),
ETD_xlBSA_293_hdb_dup$PositionA,ETD_xlBSA_293_hdb_dup$PositionB)
colnames(ETD_xlBSA_293_hdb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(ETD_xlBSA_293_hdb_xinet,"20170830_ETDMS2_xlBSA_293_hdb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
ETD1_xlBSA_293_hdb$PositionA=pmax(ETD1_xlBSA_293_hdb$Position.A,ETD1_xlBSA_293_hdb$Position.B)
ETD1_xlBSA_293_hdb$PositionB=pmin(ETD1_xlBSA_293_hdb$Position.A,ETD1_xlBSA_293_hdb$Position.B)
ETD2_xlBSA_293_hdb$PositionA=pmax(ETD2_xlBSA_293_hdb$Position.A,ETD2_xlBSA_293_hdb$Position.B)
ETD2_xlBSA_293_hdb$PositionB=pmin(ETD2_xlBSA_293_hdb$Position.A,ETD2_xlBSA_293_hdb$Position.B)
ETD3_xlBSA_293_hdb$PositionA=pmax(ETD3_xlBSA_293_hdb$Position.A,ETD3_xlBSA_293_hdb$Position.B)
ETD3_xlBSA_293_hdb$PositionB=pmin(ETD3_xlBSA_293_hdb$Position.A,ETD3_xlBSA_293_hdb$Position.B)
#Filter XlinkX Score >80
ETD1_xlBSA_293_hdb=ETD1_xlBSA_293_hdb[ETD1_xlBSA_293_hdb[,"Max.XlinkX.Score"]>s5,]
ETD2_xlBSA_293_hdb=ETD2_xlBSA_293_hdb[ETD2_xlBSA_293_hdb[,"Max.XlinkX.Score"]>s5,]
ETD3_xlBSA_293_hdb=ETD3_xlBSA_293_hdb[ETD3_xlBSA_293_hdb[,"Max.XlinkX.Score"]>s5,]
#Remove Duplicates
ETD1_xlBSA_293_hdb_dup=ETD1_xlBSA_293_hdb[!duplicated(ETD1_xlBSA_293_hdb[c("PositionA","PositionB")]),]
ETD2_xlBSA_293_hdb_dup=ETD2_xlBSA_293_hdb[!duplicated(ETD2_xlBSA_293_hdb[c("PositionA","PositionB")]),]
ETD3_xlBSA_293_hdb_dup=ETD3_xlBSA_293_hdb[!duplicated(ETD3_xlBSA_293_hdb[c("PositionA","PositionB")]),]
#ETD Number of BSA crosslinks
ETD1_xlBSA_293_hdb_bsaonly<- ETD1_xlBSA_293_hdb_dup[grepl("P02769", ETD1_xlBSA_293_hdb_dup[["Accession.A"]]) & grepl("P02769", ETD1_xlBSA_293_hdb_dup[["Accession.B"]]), ]
ETD2_xlBSA_293_hdb_bsaonly<- ETD2_xlBSA_293_hdb_dup[grepl("P02769", ETD2_xlBSA_293_hdb_dup[["Accession.A"]]) & grepl("P02769", ETD2_xlBSA_293_hdb_dup[["Accession.B"]]), ]
ETD3_xlBSA_293_hdb_bsaonly<- ETD3_xlBSA_293_hdb_dup[grepl("P02769", ETD3_xlBSA_293_hdb_dup[["Accession.A"]]) & grepl("P02769", ETD3_xlBSA_293_hdb_dup[["Accession.B"]]), ]
#ETD Number of all crosslinks, no duplicates
ETD_xlBSA_293_hdb_crosslinks=c(nrow(ETD1_xlBSA_293_hdb_dup),nrow(ETD2_xlBSA_293_hdb_dup),nrow(ETD3_xlBSA_293_hdb_dup))
#ETD Number of BSA crosslinks, no duplicates
ETD_xlBSA_293_hdb_bsaonly=c(nrow(ETD1_xlBSA_293_hdb_bsaonly),nrow(ETD2_xlBSA_293_hdb_bsaonly),nrow(ETD3_xlBSA_293_hdb_bsaonly))
#ETD Number of non-BSA crosslinks
ETD_xlBSA_293_hdb_nonbsa=(ETD_xlBSA_293_hdb_crosslinks-ETD_xlBSA_293_hdb_bsaonly)
#Output for Number of crosslinks
ETD_xlBSA_293_hdb=data.frame(ETD_xlBSA_293_hdb_crosslinks,ETD_xlBSA_293_hdb_bsaonly,ETD_xlBSA_293_hdb_nonbsa)
colnames(ETD_xlBSA_293_hdb)=c("ETD_xlBSA_293_hdb_Total","ETD_xlBSA_293_hdb_bsacrosslinks","ETD_xlBSA_293_hdb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(ETD1_xlBSA_293_hdb_bsaonly,ETD2_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(ETD1_xlBSA_293_hdb_bsaonly,ETD3_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(ETD2_xlBSA_293_hdb_bsaonly,ETD3_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB")))
tiff('ETD_xlBSA_293_hdb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(ETD1_xlBSA_293_hdb_bsaonly),area2=nrow(ETD2_xlBSA_293_hdb_bsaonly),area3=nrow(ETD3_xlBSA_293_hdb_bsaonly),
n12=nrow(merge(ETD1_xlBSA_293_hdb_bsaonly,ETD2_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(ETD1_xlBSA_293_hdb_bsaonly,ETD3_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(ETD2_xlBSA_293_hdb_bsaonly,ETD3_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(ETD_xlBSA_293_hdb_dup)-(nrow(ETD1_xlBSA_293_hdb_bsaonly)+nrow(ETD2_xlBSA_293_hdb_bsaonly)+nrow(ETD3_xlBSA_293_hdb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("ETD_xlBSA_293_hdb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()

##HCD cross-linked BSA+hek293t human database
HCD1_xlBSA_293_hdb=read.table("170629_xlBSA1_hdb_CIDMS2_HCDMS2_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
HCD2_xlBSA_293_hdb=read.table("170629_xlBSA1_hdb_CIDMS2_HCDMS2_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
HCD3_xlBSA_293_hdb=read.table("170629_xlBSA3_hdb_CIDMS2_HCDMS2_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
HCD1_xlBSA_293_hdb$Replicate=c(1)
HCD2_xlBSA_293_hdb$Replicate=c(2)
HCD3_xlBSA_293_hdb$Replicate=c(3)
##Compile all replicates
HCD_xlBSA_293_hdb=rbind(HCD1_xlBSA_293_hdb,HCD2_xlBSA_293_hdb,HCD3_xlBSA_293_hdb)
write.table(HCD_xlBSA_293_hdb,"20170830_HCDMS2_xlBSA_293_hdb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
HCD_xlBSA_293_hdb$PositionA=pmax(HCD_xlBSA_293_hdb$Position.A,HCD_xlBSA_293_hdb$Position.B)
HCD_xlBSA_293_hdb$PositionB=pmin(HCD_xlBSA_293_hdb$Position.A,HCD_xlBSA_293_hdb$Position.B)
##Keep BSA cross-links, duplicates kept
HCD_xlBSA_293_hdb_bsaonly<- HCD_xlBSA_293_hdb[grepl("P02769", HCD_xlBSA_293_hdb[["Accession.A"]]) & grepl("P02769", HCD_xlBSA_293_hdb[["Accession.B"]]), ]
write.table(HCD_xlBSA_293_hdb_bsaonly,"20170830_HCDMS2_xlBSA_293_hdb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
HCD_xlBSA_293_hdb_nonbsa<- HCD_xlBSA_293_hdb[!grepl("P02769", HCD_xlBSA_293_hdb[["Accession.A"]]) | !grepl("P02769", HCD_xlBSA_293_hdb[["Accession.B"]]), ]
write.table(HCD_xlBSA_293_hdb_nonbsa,"20170830_HCDMS2_xlBSA_293_hdb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
HCD_xlBSA_293_hdb_dup=HCD_xlBSA_293_hdb_bsaonly[!duplicated(HCD_xlBSA_293_hdb_bsaonly[c("PositionA","PositionB")]),]
write.table(HCD_xlBSA_293_hdb_dup,"20170830_HCDMS2_xlBSA_293_hdb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
HCD_xlBSA_293_hdb_dup$Accession.A=as.character(HCD_xlBSA_293_hdb_dup$Accession.A)
HCD_xlBSA_293_hdb_dup$Accession.A<-replace(HCD_xlBSA_293_hdb_dup$Accession.A,HCD_xlBSA_293_hdb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
HCD_xlBSA_293_hdb_dup$Accession.B=as.character(HCD_xlBSA_293_hdb_dup$Accession.B)
HCD_xlBSA_293_hdb_dup$Accession.B<-replace(HCD_xlBSA_293_hdb_dup$Accession.B,HCD_xlBSA_293_hdb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
HCD_xlBSA_293_hdb_xinet<-data.frame(HCD_xlBSA_293_hdb_dup$Max.XlinkX.Score,as.factor(HCD_xlBSA_293_hdb_dup$Accession.A),as.factor(HCD_xlBSA_293_hdb_dup$Accession.B),
HCD_xlBSA_293_hdb_dup$PositionA,HCD_xlBSA_293_hdb_dup$PositionB)
colnames(HCD_xlBSA_293_hdb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(HCD_xlBSA_293_hdb_xinet,"20170830_HCDMS2_xlBSA_293_hdb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
HCD1_xlBSA_293_hdb$PositionA=pmax(HCD1_xlBSA_293_hdb$Position.A,HCD1_xlBSA_293_hdb$Position.B)
HCD1_xlBSA_293_hdb$PositionB=pmin(HCD1_xlBSA_293_hdb$Position.A,HCD1_xlBSA_293_hdb$Position.B)
HCD2_xlBSA_293_hdb$PositionA=pmax(HCD2_xlBSA_293_hdb$Position.A,HCD2_xlBSA_293_hdb$Position.B)
HCD2_xlBSA_293_hdb$PositionB=pmin(HCD2_xlBSA_293_hdb$Position.A,HCD2_xlBSA_293_hdb$Position.B)
HCD3_xlBSA_293_hdb$PositionA=pmax(HCD3_xlBSA_293_hdb$Position.A,HCD3_xlBSA_293_hdb$Position.B)
HCD3_xlBSA_293_hdb$PositionB=pmin(HCD3_xlBSA_293_hdb$Position.A,HCD3_xlBSA_293_hdb$Position.B)
#Filter XlinkX Score >80
HCD1_xlBSA_293_hdb=HCD1_xlBSA_293_hdb[HCD1_xlBSA_293_hdb[,"Max.XlinkX.Score"]>s6,]
HCD2_xlBSA_293_hdb=HCD2_xlBSA_293_hdb[HCD2_xlBSA_293_hdb[,"Max.XlinkX.Score"]>s6,]
HCD3_xlBSA_293_hdb=HCD3_xlBSA_293_hdb[HCD3_xlBSA_293_hdb[,"Max.XlinkX.Score"]>s6,]
#Remove Duplicates
HCD1_xlBSA_293_hdb_dup=HCD1_xlBSA_293_hdb[!duplicated(HCD1_xlBSA_293_hdb[c("PositionA","PositionB")]),]
HCD2_xlBSA_293_hdb_dup=HCD2_xlBSA_293_hdb[!duplicated(HCD2_xlBSA_293_hdb[c("PositionA","PositionB")]),]
HCD3_xlBSA_293_hdb_dup=HCD3_xlBSA_293_hdb[!duplicated(HCD3_xlBSA_293_hdb[c("PositionA","PositionB")]),]
#HCD Number of BSA crosslinks
HCD1_xlBSA_293_hdb_bsaonly<- HCD1_xlBSA_293_hdb_dup[grepl("P02769", HCD1_xlBSA_293_hdb_dup[["Accession.A"]]) & grepl("P02769", HCD1_xlBSA_293_hdb_dup[["Accession.B"]]), ]
HCD2_xlBSA_293_hdb_bsaonly<- HCD2_xlBSA_293_hdb_dup[grepl("P02769", HCD2_xlBSA_293_hdb_dup[["Accession.A"]]) & grepl("P02769", HCD2_xlBSA_293_hdb_dup[["Accession.B"]]), ]
HCD3_xlBSA_293_hdb_bsaonly<- HCD3_xlBSA_293_hdb_dup[grepl("P02769", HCD3_xlBSA_293_hdb_dup[["Accession.A"]]) & grepl("P02769", HCD3_xlBSA_293_hdb_dup[["Accession.B"]]), ]
#HCD Number of all crosslinks, no duplicates
HCD_xlBSA_293_hdb_crosslinks=c(nrow(HCD1_xlBSA_293_hdb_dup),nrow(HCD2_xlBSA_293_hdb_dup),nrow(HCD3_xlBSA_293_hdb_dup))
#HCD Number of BSA crosslinks, no duplicates
HCD_xlBSA_293_hdb_bsaonly=c(nrow(HCD1_xlBSA_293_hdb_bsaonly),nrow(HCD2_xlBSA_293_hdb_bsaonly),nrow(HCD3_xlBSA_293_hdb_bsaonly))
#HCD Number of non-BSA crosslinks
HCD_xlBSA_293_hdb_nonbsa=(HCD_xlBSA_293_hdb_crosslinks-HCD_xlBSA_293_hdb_bsaonly)
#Output for Number of crosslinks
HCD_xlBSA_293_hdb=data.frame(HCD_xlBSA_293_hdb_crosslinks,HCD_xlBSA_293_hdb_bsaonly,HCD_xlBSA_293_hdb_nonbsa)
colnames(HCD_xlBSA_293_hdb)=c("HCD_xlBSA_293_hdb_Total","HCD_xlBSA_293_hdb_bsacrosslinks","HCD_xlBSA_293_hdb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(HCD1_xlBSA_293_hdb_bsaonly,HCD2_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(HCD1_xlBSA_293_hdb_bsaonly,HCD3_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(HCD2_xlBSA_293_hdb_bsaonly,HCD3_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB")))
tiff('HCD_xlBSA_293_hdb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(HCD1_xlBSA_293_hdb_bsaonly),area2=nrow(HCD2_xlBSA_293_hdb_bsaonly),area3=nrow(HCD3_xlBSA_293_hdb_bsaonly),
n12=nrow(merge(HCD1_xlBSA_293_hdb_bsaonly,HCD2_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(HCD1_xlBSA_293_hdb_bsaonly,HCD3_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(HCD2_xlBSA_293_hdb_bsaonly,HCD3_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(HCD_xlBSA_293_hdb_dup)-(nrow(HCD1_xlBSA_293_hdb_bsaonly)+nrow(HCD2_xlBSA_293_hdb_bsaonly)+nrow(HCD3_xlBSA_293_hdb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("HCD_xlBSA_293_hdb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()

##ETHCD cross-linked BSA+hek293t human database
ETHCD1_xlBSA_293_hdb=read.table("170629_xlBSA1_hdb_CIDMS2_EThcDMS2_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
ETHCD2_xlBSA_293_hdb=read.table("170629_xlBSA2_hdb_CIDMS2_EThcDMS2_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
ETHCD3_xlBSA_293_hdb=read.table("170629_xlBSA3_hdb_CIDMS2_EThcDMS2_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
ETHCD1_xlBSA_293_hdb$Replicate=c(1)
ETHCD2_xlBSA_293_hdb$Replicate=c(2)
ETHCD3_xlBSA_293_hdb$Replicate=c(3)
##Compile all replicates
ETHCD_xlBSA_293_hdb=rbind(ETHCD1_xlBSA_293_hdb,ETHCD2_xlBSA_293_hdb,ETHCD3_xlBSA_293_hdb)
write.table(ETHCD_xlBSA_293_hdb,"20170830_ETHCDMS2_xlBSA_293_hdb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
ETHCD_xlBSA_293_hdb$PositionA=pmax(ETHCD_xlBSA_293_hdb$Position.A,ETHCD_xlBSA_293_hdb$Position.B)
ETHCD_xlBSA_293_hdb$PositionB=pmin(ETHCD_xlBSA_293_hdb$Position.A,ETHCD_xlBSA_293_hdb$Position.B)
##Keep BSA cross-links, duplicates kept
ETHCD_xlBSA_293_hdb_bsaonly<- ETHCD_xlBSA_293_hdb[grepl("P02769", ETHCD_xlBSA_293_hdb[["Accession.A"]]) & grepl("P02769", ETHCD_xlBSA_293_hdb[["Accession.B"]]), ]
write.table(ETHCD_xlBSA_293_hdb_bsaonly,"20170830_ETHCDMS2_xlBSA_293_hdb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
ETHCD_xlBSA_293_hdb_nonbsa<- ETHCD_xlBSA_293_hdb[!grepl("P02769", ETHCD_xlBSA_293_hdb[["Accession.A"]]) | !grepl("P02769", ETHCD_xlBSA_293_hdb[["Accession.B"]]), ]
write.table(ETHCD_xlBSA_293_hdb_nonbsa,"20170830_ETHCDMS2_xlBSA_293_hdb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
ETHCD_xlBSA_293_hdb_dup=ETHCD_xlBSA_293_hdb_bsaonly[!duplicated(ETHCD_xlBSA_293_hdb_bsaonly[c("PositionA","PositionB")]),]
write.table(ETHCD_xlBSA_293_hdb_dup,"20170830_ETHCDMS2_xlBSA_293_hdb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
ETHCD_xlBSA_293_hdb_dup$Accession.A=as.character(ETHCD_xlBSA_293_hdb_dup$Accession.A)
ETHCD_xlBSA_293_hdb_dup$Accession.A<-replace(ETHCD_xlBSA_293_hdb_dup$Accession.A,ETHCD_xlBSA_293_hdb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
ETHCD_xlBSA_293_hdb_dup$Accession.B=as.character(ETHCD_xlBSA_293_hdb_dup$Accession.B)
ETHCD_xlBSA_293_hdb_dup$Accession.B<-replace(ETHCD_xlBSA_293_hdb_dup$Accession.B,ETHCD_xlBSA_293_hdb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
ETHCD_xlBSA_293_hdb_xinet<-data.frame(ETHCD_xlBSA_293_hdb_dup$Max.XlinkX.Score,as.factor(ETHCD_xlBSA_293_hdb_dup$Accession.A),as.factor(ETHCD_xlBSA_293_hdb_dup$Accession.B),
ETHCD_xlBSA_293_hdb_dup$PositionA,ETHCD_xlBSA_293_hdb_dup$PositionB)
colnames(ETHCD_xlBSA_293_hdb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(ETHCD_xlBSA_293_hdb_xinet,"20170830_ETHCDMS2_xlBSA_293_hdb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
ETHCD1_xlBSA_293_hdb$PositionA=pmax(ETHCD1_xlBSA_293_hdb$Position.A,ETHCD1_xlBSA_293_hdb$Position.B)
ETHCD1_xlBSA_293_hdb$PositionB=pmin(ETHCD1_xlBSA_293_hdb$Position.A,ETHCD1_xlBSA_293_hdb$Position.B)
ETHCD2_xlBSA_293_hdb$PositionA=pmax(ETHCD2_xlBSA_293_hdb$Position.A,ETHCD2_xlBSA_293_hdb$Position.B)
ETHCD2_xlBSA_293_hdb$PositionB=pmin(ETHCD2_xlBSA_293_hdb$Position.A,ETHCD2_xlBSA_293_hdb$Position.B)
ETHCD3_xlBSA_293_hdb$PositionA=pmax(ETHCD3_xlBSA_293_hdb$Position.A,ETHCD3_xlBSA_293_hdb$Position.B)
ETHCD3_xlBSA_293_hdb$PositionB=pmin(ETHCD3_xlBSA_293_hdb$Position.A,ETHCD3_xlBSA_293_hdb$Position.B)
#Filter XlinkX Score >80
ETHCD1_xlBSA_293_hdb=ETHCD1_xlBSA_293_hdb[ETHCD1_xlBSA_293_hdb[,"Max.XlinkX.Score"]>s7,]
ETHCD2_xlBSA_293_hdb=ETHCD2_xlBSA_293_hdb[ETHCD2_xlBSA_293_hdb[,"Max.XlinkX.Score"]>s7,]
ETHCD3_xlBSA_293_hdb=ETHCD3_xlBSA_293_hdb[ETHCD3_xlBSA_293_hdb[,"Max.XlinkX.Score"]>s7,]
#Remove Duplicates
ETHCD1_xlBSA_293_hdb_dup=ETHCD1_xlBSA_293_hdb[!duplicated(ETHCD1_xlBSA_293_hdb[c("PositionA","PositionB")]),]
ETHCD2_xlBSA_293_hdb_dup=ETHCD2_xlBSA_293_hdb[!duplicated(ETHCD2_xlBSA_293_hdb[c("PositionA","PositionB")]),]
ETHCD3_xlBSA_293_hdb_dup=ETHCD3_xlBSA_293_hdb[!duplicated(ETHCD3_xlBSA_293_hdb[c("PositionA","PositionB")]),]
#ETHCD Number of BSA crosslinks
ETHCD1_xlBSA_293_hdb_bsaonly<- ETHCD1_xlBSA_293_hdb_dup[grepl("P02769", ETHCD1_xlBSA_293_hdb_dup[["Accession.A"]]) & grepl("P02769", ETHCD1_xlBSA_293_hdb_dup[["Accession.B"]]), ]
ETHCD2_xlBSA_293_hdb_bsaonly<- ETHCD2_xlBSA_293_hdb_dup[grepl("P02769", ETHCD2_xlBSA_293_hdb_dup[["Accession.A"]]) & grepl("P02769", ETHCD2_xlBSA_293_hdb_dup[["Accession.B"]]), ]
ETHCD3_xlBSA_293_hdb_bsaonly<- ETHCD3_xlBSA_293_hdb_dup[grepl("P02769", ETHCD3_xlBSA_293_hdb_dup[["Accession.A"]]) & grepl("P02769", ETHCD3_xlBSA_293_hdb_dup[["Accession.B"]]), ]
#ETHCD Number of all crosslinks, no duplicates
ETHCD_xlBSA_293_hdb_crosslinks=c(nrow(ETHCD1_xlBSA_293_hdb_dup),nrow(ETHCD2_xlBSA_293_hdb_dup),nrow(ETHCD3_xlBSA_293_hdb_dup))
#ETHCD Number of BSA crosslinks, no duplicates
ETHCD_xlBSA_293_hdb_bsaonly=c(nrow(ETHCD1_xlBSA_293_hdb_bsaonly),nrow(ETHCD2_xlBSA_293_hdb_bsaonly),nrow(ETHCD3_xlBSA_293_hdb_bsaonly))
#ETHCD Number of non-BSA crosslinks
ETHCD_xlBSA_293_hdb_nonbsa=(ETHCD_xlBSA_293_hdb_crosslinks-ETHCD_xlBSA_293_hdb_bsaonly)
#Output for Number of crosslinks
ETHCD_xlBSA_293_hdb=data.frame(ETHCD_xlBSA_293_hdb_crosslinks,ETHCD_xlBSA_293_hdb_bsaonly,ETHCD_xlBSA_293_hdb_nonbsa)
colnames(ETHCD_xlBSA_293_hdb)=c("ETHCD_xlBSA_293_hdb_Total","ETHCD_xlBSA_293_hdb_bsacrosslinks","ETHCD_xlBSA_293_hdb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(ETHCD1_xlBSA_293_hdb_bsaonly,ETHCD2_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(ETHCD1_xlBSA_293_hdb_bsaonly,ETHCD3_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(ETHCD2_xlBSA_293_hdb_bsaonly,ETHCD3_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB")))
tiff('ETHCD_xlBSA_293_hdb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(ETHCD1_xlBSA_293_hdb_bsaonly),area2=nrow(ETHCD2_xlBSA_293_hdb_bsaonly),area3=nrow(ETHCD3_xlBSA_293_hdb_bsaonly),
n12=nrow(merge(ETHCD1_xlBSA_293_hdb_bsaonly,ETHCD2_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(ETHCD1_xlBSA_293_hdb_bsaonly,ETHCD3_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(ETHCD2_xlBSA_293_hdb_bsaonly,ETHCD3_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(ETHCD_xlBSA_293_hdb_dup)-(nrow(ETHCD1_xlBSA_293_hdb_bsaonly)+nrow(ETHCD2_xlBSA_293_hdb_bsaonly)+nrow(ETHCD3_xlBSA_293_hdb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("ETHCD_xlBSA_293_hdb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()

##HCDMS3 cross-linked BSA+hek293t human database
HCDMS31_xlBSA_293_hdb=read.table("170629_xlBSA1_hdb_CIDMS2_HCDMS3_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
HCDMS32_xlBSA_293_hdb=read.table("170629_xlBSA2_hdb_CIDMS2_HCDMS3_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
HCDMS33_xlBSA_293_hdb=read.table("170629_xlBSA3_hdb_CIDMS2_HCDMS3_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
HCDMS31_xlBSA_293_hdb$Replicate=c(1)
HCDMS32_xlBSA_293_hdb$Replicate=c(2)
HCDMS33_xlBSA_293_hdb$Replicate=c(3)
##Compile all replicates
HCDMS3_xlBSA_293_hdb=rbind(HCDMS31_xlBSA_293_hdb,HCDMS32_xlBSA_293_hdb,HCDMS33_xlBSA_293_hdb)
write.table(HCDMS3_xlBSA_293_hdb,"20170830_HCDMS3_xlBSA_293_hdb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
HCDMS3_xlBSA_293_hdb$PositionA=pmax(HCDMS3_xlBSA_293_hdb$Position.A,HCDMS3_xlBSA_293_hdb$Position.B)
HCDMS3_xlBSA_293_hdb$PositionB=pmin(HCDMS3_xlBSA_293_hdb$Position.A,HCDMS3_xlBSA_293_hdb$Position.B)
##Keep BSA cross-links, duplicates kept
HCDMS3_xlBSA_293_hdb_bsaonly<- HCDMS3_xlBSA_293_hdb[grepl("P02769", HCDMS3_xlBSA_293_hdb[["Accession.A"]]) & grepl("P02769", HCDMS3_xlBSA_293_hdb[["Accession.B"]]), ]
write.table(HCDMS3_xlBSA_293_hdb_bsaonly,"20170830_HCDMS3_xlBSA_293_hdb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
HCDMS3_xlBSA_293_hdb_nonbsa<- HCDMS3_xlBSA_293_hdb[!grepl("P02769", HCDMS3_xlBSA_293_hdb[["Accession.A"]]) | !grepl("P02769", HCDMS3_xlBSA_293_hdb[["Accession.B"]]), ]
write.table(HCDMS3_xlBSA_293_hdb_nonbsa,"20170830_HCDMS3_xlBSA_293_hdb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
HCDMS3_xlBSA_293_hdb_dup=HCDMS3_xlBSA_293_hdb_bsaonly[!duplicated(HCDMS3_xlBSA_293_hdb_bsaonly[c("PositionA","PositionB")]),]
write.table(HCDMS3_xlBSA_293_hdb_dup,"20170830_HCDMS3_xlBSA_293_hdb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
HCDMS3_xlBSA_293_hdb_dup$Accession.A=as.character(HCDMS3_xlBSA_293_hdb_dup$Accession.A)
HCDMS3_xlBSA_293_hdb_dup$Accession.A<-replace(HCDMS3_xlBSA_293_hdb_dup$Accession.A,HCDMS3_xlBSA_293_hdb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
HCDMS3_xlBSA_293_hdb_dup$Accession.B=as.character(HCDMS3_xlBSA_293_hdb_dup$Accession.B)
HCDMS3_xlBSA_293_hdb_dup$Accession.B<-replace(HCDMS3_xlBSA_293_hdb_dup$Accession.B,HCDMS3_xlBSA_293_hdb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
HCDMS3_xlBSA_293_hdb_xinet<-data.frame(HCDMS3_xlBSA_293_hdb_dup$Max.XlinkX.Score,as.factor(HCDMS3_xlBSA_293_hdb_dup$Accession.A),as.factor(HCDMS3_xlBSA_293_hdb_dup$Accession.B),
HCDMS3_xlBSA_293_hdb_dup$PositionA,HCDMS3_xlBSA_293_hdb_dup$PositionB)
colnames(HCDMS3_xlBSA_293_hdb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(HCDMS3_xlBSA_293_hdb_xinet,"20170830_HCDMS3_xlBSA_293_hdb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
HCDMS31_xlBSA_293_hdb$PositionA=pmax(HCDMS31_xlBSA_293_hdb$Position.A,HCDMS31_xlBSA_293_hdb$Position.B)
HCDMS31_xlBSA_293_hdb$PositionB=pmin(HCDMS31_xlBSA_293_hdb$Position.A,HCDMS31_xlBSA_293_hdb$Position.B)
HCDMS32_xlBSA_293_hdb$PositionA=pmax(HCDMS32_xlBSA_293_hdb$Position.A,HCDMS32_xlBSA_293_hdb$Position.B)
HCDMS32_xlBSA_293_hdb$PositionB=pmin(HCDMS32_xlBSA_293_hdb$Position.A,HCDMS32_xlBSA_293_hdb$Position.B)
HCDMS33_xlBSA_293_hdb$PositionA=pmax(HCDMS33_xlBSA_293_hdb$Position.A,HCDMS33_xlBSA_293_hdb$Position.B)
HCDMS33_xlBSA_293_hdb$PositionB=pmin(HCDMS33_xlBSA_293_hdb$Position.A,HCDMS33_xlBSA_293_hdb$Position.B)
#Filter XlinkX Score >80
HCDMS31_xlBSA_293_hdb=HCDMS31_xlBSA_293_hdb[HCDMS31_xlBSA_293_hdb[,"Max.XlinkX.Score"]>s8,]
HCDMS32_xlBSA_293_hdb=HCDMS32_xlBSA_293_hdb[HCDMS32_xlBSA_293_hdb[,"Max.XlinkX.Score"]>s8,]
HCDMS33_xlBSA_293_hdb=HCDMS33_xlBSA_293_hdb[HCDMS33_xlBSA_293_hdb[,"Max.XlinkX.Score"]>s8,]
#Remove Duplicates
HCDMS31_xlBSA_293_hdb_dup=HCDMS31_xlBSA_293_hdb[!duplicated(HCDMS31_xlBSA_293_hdb[c("PositionA","PositionB")]),]
HCDMS32_xlBSA_293_hdb_dup=HCDMS32_xlBSA_293_hdb[!duplicated(HCDMS32_xlBSA_293_hdb[c("PositionA","PositionB")]),]
HCDMS33_xlBSA_293_hdb_dup=HCDMS33_xlBSA_293_hdb[!duplicated(HCDMS33_xlBSA_293_hdb[c("PositionA","PositionB")]),]
#HCDMS3 Number of BSA crosslinks
HCDMS31_xlBSA_293_hdb_bsaonly<- HCDMS31_xlBSA_293_hdb_dup[grepl("P02769", HCDMS31_xlBSA_293_hdb_dup[["Accession.A"]]) & grepl("P02769", HCDMS31_xlBSA_293_hdb_dup[["Accession.B"]]), ]
HCDMS32_xlBSA_293_hdb_bsaonly<- HCDMS32_xlBSA_293_hdb_dup[grepl("P02769", HCDMS32_xlBSA_293_hdb_dup[["Accession.A"]]) & grepl("P02769", HCDMS32_xlBSA_293_hdb_dup[["Accession.B"]]), ]
HCDMS33_xlBSA_293_hdb_bsaonly<- HCDMS33_xlBSA_293_hdb_dup[grepl("P02769", HCDMS33_xlBSA_293_hdb_dup[["Accession.A"]]) & grepl("P02769", HCDMS33_xlBSA_293_hdb_dup[["Accession.B"]]), ]
#HCDMS3 Number of all crosslinks, no duplicates
HCDMS3_xlBSA_293_hdb_crosslinks=c(nrow(HCDMS31_xlBSA_293_hdb_dup),nrow(HCDMS32_xlBSA_293_hdb_dup),nrow(HCDMS33_xlBSA_293_hdb_dup))
#HCDMS3 Number of BSA crosslinks, no duplicates
HCDMS3_xlBSA_293_hdb_bsaonly=c(nrow(HCDMS31_xlBSA_293_hdb_bsaonly),nrow(HCDMS32_xlBSA_293_hdb_bsaonly),nrow(HCDMS33_xlBSA_293_hdb_bsaonly))
#HCDMS3 Number of non-BSA crosslinks
HCDMS3_xlBSA_293_hdb_nonbsa=(HCDMS3_xlBSA_293_hdb_crosslinks-HCDMS3_xlBSA_293_hdb_bsaonly)
#Output for Number of crosslinks
HCDMS3_xlBSA_293_hdb=data.frame(HCDMS3_xlBSA_293_hdb_crosslinks,HCDMS3_xlBSA_293_hdb_bsaonly,HCDMS3_xlBSA_293_hdb_nonbsa)
colnames(HCDMS3_xlBSA_293_hdb)=c("HCDMS3_xlBSA_293_hdb_Total","HCDMS3_xlBSA_293_hdb_bsacrosslinks","HCDMS3_xlBSA_293_hdb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(HCDMS31_xlBSA_293_hdb_bsaonly,HCDMS32_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(HCDMS31_xlBSA_293_hdb_bsaonly,HCDMS33_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(HCDMS32_xlBSA_293_hdb_bsaonly,HCDMS33_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB")))
tiff('HCDMS3_xlBSA_293_hdb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(HCDMS31_xlBSA_293_hdb_bsaonly),area2=nrow(HCDMS32_xlBSA_293_hdb_bsaonly),area3=nrow(HCDMS33_xlBSA_293_hdb_bsaonly),
n12=nrow(merge(HCDMS31_xlBSA_293_hdb_bsaonly,HCDMS32_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(HCDMS31_xlBSA_293_hdb_bsaonly,HCDMS33_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(HCDMS32_xlBSA_293_hdb_bsaonly,HCDMS33_xlBSA_293_hdb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(HCDMS3_xlBSA_293_hdb_dup)-(nrow(HCDMS31_xlBSA_293_hdb_bsaonly)+nrow(HCDMS32_xlBSA_293_hdb_bsaonly)+nrow(HCDMS33_xlBSA_293_hdb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("HCDMS3_xlBSA_293_hdb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()

##ETD cross-linked BSA BSA database
ETD1_xlBSA_bsadb=read.table("170629_xlBSA1_bsadb_CIDMS2_ETDMS2_Crosslinks.txt", fill=TRUE, header=TRUE)
ETD2_xlBSA_bsadb=read.table("170629_xlBSA2_bsadb_CIDMS2_ETDMS2_Crosslinks.txt", fill=TRUE, header=TRUE)
ETD3_xlBSA_bsadb=read.table("170629_xlBSA3_bsadb_CIDMS2_ETDMS2_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
ETD1_xlBSA_bsadb$Replicate=c(1)
ETD2_xlBSA_bsadb$Replicate=c(2)
ETD3_xlBSA_bsadb$Replicate=c(3)
##Compile all replicates
ETD_xlBSA_bsadb=rbind(ETD1_xlBSA_bsadb,ETD2_xlBSA_bsadb,ETD3_xlBSA_bsadb)
write.table(ETD_xlBSA_bsadb,"20170830_ETDMS2_xlBSA_bsadb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
ETD_xlBSA_bsadb$PositionA=pmax(ETD_xlBSA_bsadb$Position.A,ETD_xlBSA_bsadb$Position.B)
ETD_xlBSA_bsadb$PositionB=pmin(ETD_xlBSA_bsadb$Position.A,ETD_xlBSA_bsadb$Position.B)
##Keep BSA cross-links, duplicates kept
ETD_xlBSA_bsadb_bsaonly<- ETD_xlBSA_bsadb[grepl("P02769", ETD_xlBSA_bsadb[["Accession.A"]]) & grepl("P02769", ETD_xlBSA_bsadb[["Accession.B"]]), ]
write.table(ETD_xlBSA_bsadb_bsaonly,"20170830_ETDMS2_xlBSA_bsadb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
ETD_xlBSA_bsadb_nonbsa<- ETD_xlBSA_bsadb[!grepl("P02769", ETD_xlBSA_bsadb[["Accession.A"]]) | !grepl("P02769", ETD_xlBSA_bsadb[["Accession.B"]]), ]
write.table(ETD_xlBSA_bsadb_nonbsa,"20170830_ETDMS2_xlBSA_bsadb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
ETD_xlBSA_bsadb_dup=ETD_xlBSA_bsadb_bsaonly[!duplicated(ETD_xlBSA_bsadb_bsaonly[c("PositionA","PositionB")]),]
write.table(ETD_xlBSA_bsadb_dup,"20170830_ETDMS2_xlBSA_bsadb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
ETD_xlBSA_bsadb_dup$Accession.A=as.character(ETD_xlBSA_bsadb_dup$Accession.A)
ETD_xlBSA_bsadb_dup$Accession.A<-replace(ETD_xlBSA_bsadb_dup$Accession.A,ETD_xlBSA_bsadb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
ETD_xlBSA_bsadb_dup$Accession.B=as.character(ETD_xlBSA_bsadb_dup$Accession.B)
ETD_xlBSA_bsadb_dup$Accession.B<-replace(ETD_xlBSA_bsadb_dup$Accession.B,ETD_xlBSA_bsadb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
ETD_xlBSA_bsadb_xinet<-data.frame(ETD_xlBSA_bsadb_dup$Max.XlinkX.Score,as.factor(ETD_xlBSA_bsadb_dup$Accession.A),as.factor(ETD_xlBSA_bsadb_dup$Accession.B),
ETD_xlBSA_bsadb_dup$PositionA,ETD_xlBSA_bsadb_dup$PositionB)
colnames(ETD_xlBSA_bsadb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(ETD_xlBSA_bsadb_xinet,"20170830_ETDMS2_xlBSA_bsadb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
ETD1_xlBSA_bsadb$PositionA=pmax(ETD1_xlBSA_bsadb$Position.A,ETD1_xlBSA_bsadb$Position.B)
ETD1_xlBSA_bsadb$PositionB=pmin(ETD1_xlBSA_bsadb$Position.A,ETD1_xlBSA_bsadb$Position.B)
ETD2_xlBSA_bsadb$PositionA=pmax(ETD2_xlBSA_bsadb$Position.A,ETD2_xlBSA_bsadb$Position.B)
ETD2_xlBSA_bsadb$PositionB=pmin(ETD2_xlBSA_bsadb$Position.A,ETD2_xlBSA_bsadb$Position.B)
ETD3_xlBSA_bsadb$PositionA=pmax(ETD3_xlBSA_bsadb$Position.A,ETD3_xlBSA_bsadb$Position.B)
ETD3_xlBSA_bsadb$PositionB=pmin(ETD3_xlBSA_bsadb$Position.A,ETD3_xlBSA_bsadb$Position.B)
#Filter XlinkX Score >80
ETD1_xlBSA_bsadb=ETD1_xlBSA_bsadb[ETD1_xlBSA_bsadb[,"Max.XlinkX.Score"]>s9,]
ETD2_xlBSA_bsadb=ETD2_xlBSA_bsadb[ETD2_xlBSA_bsadb[,"Max.XlinkX.Score"]>s9,]
ETD3_xlBSA_bsadb=ETD3_xlBSA_bsadb[ETD3_xlBSA_bsadb[,"Max.XlinkX.Score"]>s9,]
#Remove Duplicates
ETD1_xlBSA_bsadb_dup=ETD1_xlBSA_bsadb[!duplicated(ETD1_xlBSA_bsadb[c("PositionA","PositionB")]),]
ETD2_xlBSA_bsadb_dup=ETD2_xlBSA_bsadb[!duplicated(ETD2_xlBSA_bsadb[c("PositionA","PositionB")]),]
ETD3_xlBSA_bsadb_dup=ETD3_xlBSA_bsadb[!duplicated(ETD3_xlBSA_bsadb[c("PositionA","PositionB")]),]
#ETD Number of BSA crosslinks
ETD1_xlBSA_bsadb_bsaonly<- ETD1_xlBSA_bsadb_dup[grepl("P02769", ETD1_xlBSA_bsadb_dup[["Accession.A"]]) & grepl("P02769", ETD1_xlBSA_bsadb_dup[["Accession.B"]]), ]
ETD2_xlBSA_bsadb_bsaonly<- ETD2_xlBSA_bsadb_dup[grepl("P02769", ETD2_xlBSA_bsadb_dup[["Accession.A"]]) & grepl("P02769", ETD2_xlBSA_bsadb_dup[["Accession.B"]]), ]
ETD3_xlBSA_bsadb_bsaonly<- ETD3_xlBSA_bsadb_dup[grepl("P02769", ETD3_xlBSA_bsadb_dup[["Accession.A"]]) & grepl("P02769", ETD3_xlBSA_bsadb_dup[["Accession.B"]]), ]
#ETD Number of all crosslinks, no duplicates
ETD_xlBSA_bsadb_crosslinks=c(nrow(ETD1_xlBSA_bsadb_dup),nrow(ETD2_xlBSA_bsadb_dup),nrow(ETD3_xlBSA_bsadb_dup))
#ETD Number of BSA crosslinks, no duplicates
ETD_xlBSA_bsadb_bsaonly=c(nrow(ETD1_xlBSA_bsadb_bsaonly),nrow(ETD2_xlBSA_bsadb_bsaonly),nrow(ETD3_xlBSA_bsadb_bsaonly))
#ETD Number of non-BSA crosslinks
ETD_xlBSA_bsadb_nonbsa=(ETD_xlBSA_bsadb_crosslinks-ETD_xlBSA_bsadb_bsaonly)
#Output for Number of crosslinks
ETD_xlBSA_bsadb=data.frame(ETD_xlBSA_bsadb_crosslinks,ETD_xlBSA_bsadb_bsaonly,ETD_xlBSA_bsadb_nonbsa)
colnames(ETD_xlBSA_bsadb)=c("ETD_xlBSA_bsadb_Total","ETD_xlBSA_bsadb_bsacrosslinks","ETD_xlBSA_bsadb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(ETD1_xlBSA_bsadb_bsaonly,ETD2_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(ETD1_xlBSA_bsadb_bsaonly,ETD3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(ETD2_xlBSA_bsadb_bsaonly,ETD3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB")))
tiff('ETD_xlBSA_bsadb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(ETD1_xlBSA_bsadb_bsaonly),area2=nrow(ETD2_xlBSA_bsadb_bsaonly),area3=nrow(ETD3_xlBSA_bsadb_bsaonly),
n12=nrow(merge(ETD1_xlBSA_bsadb_bsaonly,ETD2_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(ETD1_xlBSA_bsadb_bsaonly,ETD3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(ETD2_xlBSA_bsadb_bsaonly,ETD3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(ETD_xlBSA_bsadb_dup)-(nrow(ETD1_xlBSA_bsadb_bsaonly)+nrow(ETD2_xlBSA_bsadb_bsaonly)+nrow(ETD3_xlBSA_bsadb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("ETD_xlBSA_bsadb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()

##HCD cross-linked BSA BSA database
HCD1_xlBSA_bsadb=read.table("170629_xlBSA1_bsadb_CIDMS2_HCDMS2_Crosslinks.txt", fill=TRUE, header=TRUE)
HCD2_xlBSA_bsadb=read.table("170629_xlBSA2_bsadb_CIDMS2_HCDMS2_Crosslinks.txt", fill=TRUE, header=TRUE)
HCD3_xlBSA_bsadb=read.table("170629_xlBSA3_bsadb_CIDMS2_HCDMS2_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
HCD1_xlBSA_bsadb$Replicate=c(1)
HCD2_xlBSA_bsadb$Replicate=c(2)
HCD3_xlBSA_bsadb$Replicate=c(3)
##Compile all replicates
HCD_xlBSA_bsadb=rbind(HCD1_xlBSA_bsadb,HCD2_xlBSA_bsadb,HCD3_xlBSA_bsadb)
write.table(HCD_xlBSA_bsadb,"20170830_HCDMS2_xlBSA_bsadb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
HCD_xlBSA_bsadb$PositionA=pmax(HCD_xlBSA_bsadb$Position.A,HCD_xlBSA_bsadb$Position.B)
HCD_xlBSA_bsadb$PositionB=pmin(HCD_xlBSA_bsadb$Position.A,HCD_xlBSA_bsadb$Position.B)
##Keep BSA cross-links, duplicates kept
HCD_xlBSA_bsadb_bsaonly<- HCD_xlBSA_bsadb[grepl("P02769", HCD_xlBSA_bsadb[["Accession.A"]]) & grepl("P02769", HCD_xlBSA_bsadb[["Accession.B"]]), ]
write.table(HCD_xlBSA_bsadb_bsaonly,"20170830_HCDMS2_xlBSA_bsadb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
HCD_xlBSA_bsadb_nonbsa<- HCD_xlBSA_bsadb[!grepl("P02769", HCD_xlBSA_bsadb[["Accession.A"]]) | !grepl("P02769", HCD_xlBSA_bsadb[["Accession.B"]]), ]
write.table(HCD_xlBSA_bsadb_nonbsa,"20170830_HCDMS2_xlBSA_bsadb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
HCD_xlBSA_bsadb_dup=HCD_xlBSA_bsadb_bsaonly[!duplicated(HCD_xlBSA_bsadb_bsaonly[c("PositionA","PositionB")]),]
write.table(HCD_xlBSA_bsadb_dup,"20170830_HCDMS2_xlBSA_bsadb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
HCD_xlBSA_bsadb_dup$Accession.A=as.character(HCD_xlBSA_bsadb_dup$Accession.A)
HCD_xlBSA_bsadb_dup$Accession.A<-replace(HCD_xlBSA_bsadb_dup$Accession.A,HCD_xlBSA_bsadb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
HCD_xlBSA_bsadb_dup$Accession.B=as.character(HCD_xlBSA_bsadb_dup$Accession.B)
HCD_xlBSA_bsadb_dup$Accession.B<-replace(HCD_xlBSA_bsadb_dup$Accession.B,HCD_xlBSA_bsadb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
HCD_xlBSA_bsadb_xinet<-data.frame(HCD_xlBSA_bsadb_dup$Max.XlinkX.Score,as.factor(HCD_xlBSA_bsadb_dup$Accession.A),as.factor(HCD_xlBSA_bsadb_dup$Accession.B),
HCD_xlBSA_bsadb_dup$PositionA,HCD_xlBSA_bsadb_dup$PositionB)
colnames(HCD_xlBSA_bsadb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(HCD_xlBSA_bsadb_xinet,"20170830_HCDMS2_xlBSA_bsadb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
HCD1_xlBSA_bsadb$PositionA=pmax(HCD1_xlBSA_bsadb$Position.A,HCD1_xlBSA_bsadb$Position.B)
HCD1_xlBSA_bsadb$PositionB=pmin(HCD1_xlBSA_bsadb$Position.A,HCD1_xlBSA_bsadb$Position.B)
HCD2_xlBSA_bsadb$PositionA=pmax(HCD2_xlBSA_bsadb$Position.A,HCD2_xlBSA_bsadb$Position.B)
HCD2_xlBSA_bsadb$PositionB=pmin(HCD2_xlBSA_bsadb$Position.A,HCD2_xlBSA_bsadb$Position.B)
HCD3_xlBSA_bsadb$PositionA=pmax(HCD3_xlBSA_bsadb$Position.A,HCD3_xlBSA_bsadb$Position.B)
HCD3_xlBSA_bsadb$PositionB=pmin(HCD3_xlBSA_bsadb$Position.A,HCD3_xlBSA_bsadb$Position.B)
#Filter XlinkX Score >80
HCD1_xlBSA_bsadb=HCD1_xlBSA_bsadb[HCD1_xlBSA_bsadb[,"Max.XlinkX.Score"]>s10,]
HCD2_xlBSA_bsadb=HCD2_xlBSA_bsadb[HCD2_xlBSA_bsadb[,"Max.XlinkX.Score"]>s10,]
HCD3_xlBSA_bsadb=HCD3_xlBSA_bsadb[HCD3_xlBSA_bsadb[,"Max.XlinkX.Score"]>s10,]
#Remove Duplicates
HCD1_xlBSA_bsadb_dup=HCD1_xlBSA_bsadb[!duplicated(HCD1_xlBSA_bsadb[c("PositionA","PositionB")]),]
HCD2_xlBSA_bsadb_dup=HCD2_xlBSA_bsadb[!duplicated(HCD2_xlBSA_bsadb[c("PositionA","PositionB")]),]
HCD3_xlBSA_bsadb_dup=HCD3_xlBSA_bsadb[!duplicated(HCD3_xlBSA_bsadb[c("PositionA","PositionB")]),]
#HCD Number of BSA crosslinks
HCD1_xlBSA_bsadb_bsaonly<- HCD1_xlBSA_bsadb_dup[grepl("P02769", HCD1_xlBSA_bsadb_dup[["Accession.A"]]) & grepl("P02769", HCD1_xlBSA_bsadb_dup[["Accession.B"]]), ]
HCD2_xlBSA_bsadb_bsaonly<- HCD2_xlBSA_bsadb_dup[grepl("P02769", HCD2_xlBSA_bsadb_dup[["Accession.A"]]) & grepl("P02769", HCD2_xlBSA_bsadb_dup[["Accession.B"]]), ]
HCD3_xlBSA_bsadb_bsaonly<- HCD3_xlBSA_bsadb_dup[grepl("P02769", HCD3_xlBSA_bsadb_dup[["Accession.A"]]) & grepl("P02769", HCD3_xlBSA_bsadb_dup[["Accession.B"]]), ]
#HCD Number of all crosslinks, no duplicates
HCD_xlBSA_bsadb_crosslinks=c(nrow(HCD1_xlBSA_bsadb_dup),nrow(HCD2_xlBSA_bsadb_dup),nrow(HCD3_xlBSA_bsadb_dup))
#HCD Number of BSA crosslinks, no duplicates
HCD_xlBSA_bsadb_bsaonly=c(nrow(HCD1_xlBSA_bsadb_bsaonly),nrow(HCD2_xlBSA_bsadb_bsaonly),nrow(HCD3_xlBSA_bsadb_bsaonly))
#HCD Number of non-BSA crosslinks
HCD_xlBSA_bsadb_nonbsa=(HCD_xlBSA_bsadb_crosslinks-HCD_xlBSA_bsadb_bsaonly)
#Output for Number of crosslinks
HCD_xlBSA_bsadb=data.frame(HCD_xlBSA_bsadb_crosslinks,HCD_xlBSA_bsadb_bsaonly,HCD_xlBSA_bsadb_nonbsa)
colnames(HCD_xlBSA_bsadb)=c("HCD_xlBSA_bsadb_Total","HCD_xlBSA_bsadb_bsacrosslinks","HCD_xlBSA_bsadb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(HCD1_xlBSA_bsadb_bsaonly,HCD2_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(HCD1_xlBSA_bsadb_bsaonly,HCD3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(HCD2_xlBSA_bsadb_bsaonly,HCD3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB")))
tiff('HCD_xlBSA_bsadb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(HCD1_xlBSA_bsadb_bsaonly),area2=nrow(HCD2_xlBSA_bsadb_bsaonly),area3=nrow(HCD3_xlBSA_bsadb_bsaonly),
n12=nrow(merge(HCD1_xlBSA_bsadb_bsaonly,HCD2_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(HCD1_xlBSA_bsadb_bsaonly,HCD3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(HCD2_xlBSA_bsadb_bsaonly,HCD3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(HCD_xlBSA_bsadb_dup)-(nrow(HCD1_xlBSA_bsadb_bsaonly)+nrow(HCD2_xlBSA_bsadb_bsaonly)+nrow(HCD3_xlBSA_bsadb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("HCD_xlBSA_bsadb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()

##ETHCD cross-linked BSA+hek293 BSA database
ETHCD1_xlBSA_bsadb=read.table("170629_xlBSA1_bsadb_CIDMS2_EThcDMS2_Crosslinks.txt", fill=TRUE, header=TRUE)
ETHCD2_xlBSA_bsadb=read.table("170629_xlBSA2_bsadb_CIDMS2_EThcDMS2_Crosslinks.txt", fill=TRUE, header=TRUE)
ETHCD3_xlBSA_bsadb=read.table("170629_xlBSA3_bsadb_CIDMS2_EThcDMS2_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
ETHCD1_xlBSA_bsadb$Replicate=c(1)
ETHCD2_xlBSA_bsadb$Replicate=c(2)
ETHCD3_xlBSA_bsadb$Replicate=c(3)
##Compile all replicates
ETHCD_xlBSA_bsadb=rbind(ETHCD1_xlBSA_bsadb,ETHCD2_xlBSA_bsadb,ETHCD3_xlBSA_bsadb)
write.table(ETHCD_xlBSA_bsadb,"20170830_ETHCDMS2_xlBSA_bsadb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
ETHCD_xlBSA_bsadb$PositionA=pmax(ETHCD_xlBSA_bsadb$Position.A,ETHCD_xlBSA_bsadb$Position.B)
ETHCD_xlBSA_bsadb$PositionB=pmin(ETHCD_xlBSA_bsadb$Position.A,ETHCD_xlBSA_bsadb$Position.B)
##Keep BSA cross-links, duplicates kept
ETHCD_xlBSA_bsadb_bsaonly<- ETHCD_xlBSA_bsadb[grepl("P02769", ETHCD_xlBSA_bsadb[["Accession.A"]]) & grepl("P02769", ETHCD_xlBSA_bsadb[["Accession.B"]]), ]
write.table(ETHCD_xlBSA_bsadb_bsaonly,"20170830_ETHCDMS2_xlBSA_bsadb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
ETHCD_xlBSA_bsadb_nonbsa<- ETHCD_xlBSA_bsadb[!grepl("P02769", ETHCD_xlBSA_bsadb[["Accession.A"]]) | !grepl("P02769", ETHCD_xlBSA_bsadb[["Accession.B"]]), ]
write.table(ETHCD_xlBSA_bsadb_nonbsa,"20170830_ETHCDMS2_xlBSA_bsadb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
ETHCD_xlBSA_bsadb_dup=ETHCD_xlBSA_bsadb_bsaonly[!duplicated(ETHCD_xlBSA_bsadb_bsaonly[c("PositionA","PositionB")]),]
write.table(ETHCD_xlBSA_bsadb_dup,"20170830_ETHCDMS2_xlBSA_bsadb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
ETHCD_xlBSA_bsadb_dup$Accession.A=as.character(ETHCD_xlBSA_bsadb_dup$Accession.A)
ETHCD_xlBSA_bsadb_dup$Accession.A<-replace(ETHCD_xlBSA_bsadb_dup$Accession.A,ETHCD_xlBSA_bsadb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
ETHCD_xlBSA_bsadb_dup$Accession.B=as.character(ETHCD_xlBSA_bsadb_dup$Accession.B)
ETHCD_xlBSA_bsadb_dup$Accession.B<-replace(ETHCD_xlBSA_bsadb_dup$Accession.B,ETHCD_xlBSA_bsadb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
ETHCD_xlBSA_bsadb_xinet<-data.frame(ETHCD_xlBSA_bsadb_dup$Max.XlinkX.Score,as.factor(ETHCD_xlBSA_bsadb_dup$Accession.A),as.factor(ETHCD_xlBSA_bsadb_dup$Accession.B),
ETHCD_xlBSA_bsadb_dup$PositionA,ETHCD_xlBSA_bsadb_dup$PositionB)
colnames(ETHCD_xlBSA_bsadb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(ETHCD_xlBSA_bsadb_xinet,"20170830_ETHCDMS2_xlBSA_bsadb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
ETHCD1_xlBSA_bsadb$PositionA=pmax(ETHCD1_xlBSA_bsadb$Position.A,ETHCD1_xlBSA_bsadb$Position.B)
ETHCD1_xlBSA_bsadb$PositionB=pmin(ETHCD1_xlBSA_bsadb$Position.A,ETHCD1_xlBSA_bsadb$Position.B)
ETHCD2_xlBSA_bsadb$PositionA=pmax(ETHCD2_xlBSA_bsadb$Position.A,ETHCD2_xlBSA_bsadb$Position.B)
ETHCD2_xlBSA_bsadb$PositionB=pmin(ETHCD2_xlBSA_bsadb$Position.A,ETHCD2_xlBSA_bsadb$Position.B)
ETHCD3_xlBSA_bsadb$PositionA=pmax(ETHCD3_xlBSA_bsadb$Position.A,ETHCD3_xlBSA_bsadb$Position.B)
ETHCD3_xlBSA_bsadb$PositionB=pmin(ETHCD3_xlBSA_bsadb$Position.A,ETHCD3_xlBSA_bsadb$Position.B)
#Filter XlinkX Score >80
ETHCD1_xlBSA_bsadb=ETHCD1_xlBSA_bsadb[ETHCD1_xlBSA_bsadb[,"Max.XlinkX.Score"]>s11,]
ETHCD2_xlBSA_bsadb=ETHCD2_xlBSA_bsadb[ETHCD2_xlBSA_bsadb[,"Max.XlinkX.Score"]>s11,]
ETHCD3_xlBSA_bsadb=ETHCD3_xlBSA_bsadb[ETHCD3_xlBSA_bsadb[,"Max.XlinkX.Score"]>s11,]
#Remove Duplicates
ETHCD1_xlBSA_bsadb_dup=ETHCD1_xlBSA_bsadb[!duplicated(ETHCD1_xlBSA_bsadb[c("PositionA","PositionB")]),]
ETHCD2_xlBSA_bsadb_dup=ETHCD2_xlBSA_bsadb[!duplicated(ETHCD2_xlBSA_bsadb[c("PositionA","PositionB")]),]
ETHCD3_xlBSA_bsadb_dup=ETHCD3_xlBSA_bsadb[!duplicated(ETHCD3_xlBSA_bsadb[c("PositionA","PositionB")]),]
#ETHCD Number of BSA crosslinks
ETHCD1_xlBSA_bsadb_bsaonly<- ETHCD1_xlBSA_bsadb_dup[grepl("P02769", ETHCD1_xlBSA_bsadb_dup[["Accession.A"]]) & grepl("P02769", ETHCD1_xlBSA_bsadb_dup[["Accession.B"]]), ]
ETHCD2_xlBSA_bsadb_bsaonly<- ETHCD2_xlBSA_bsadb_dup[grepl("P02769", ETHCD2_xlBSA_bsadb_dup[["Accession.A"]]) & grepl("P02769", ETHCD2_xlBSA_bsadb_dup[["Accession.B"]]), ]
ETHCD3_xlBSA_bsadb_bsaonly<- ETHCD3_xlBSA_bsadb_dup[grepl("P02769", ETHCD3_xlBSA_bsadb_dup[["Accession.A"]]) & grepl("P02769", ETHCD3_xlBSA_bsadb_dup[["Accession.B"]]), ]
#ETHCD Number of all crosslinks, no duplicates
ETHCD_xlBSA_bsadb_crosslinks=c(nrow(ETHCD1_xlBSA_bsadb_dup),nrow(ETHCD2_xlBSA_bsadb_dup),nrow(ETHCD3_xlBSA_bsadb_dup))
#ETHCD Number of BSA crosslinks, no duplicates
ETHCD_xlBSA_bsadb_bsaonly=c(nrow(ETHCD1_xlBSA_bsadb_bsaonly),nrow(ETHCD2_xlBSA_bsadb_bsaonly),nrow(ETHCD3_xlBSA_bsadb_bsaonly))
#ETHCD Number of non-BSA crosslinks
ETHCD_xlBSA_bsadb_nonbsa=(ETHCD_xlBSA_bsadb_crosslinks-ETHCD_xlBSA_bsadb_bsaonly)
#Output for Number of crosslinks
ETHCD_xlBSA_bsadb=data.frame(ETHCD_xlBSA_bsadb_crosslinks,ETHCD_xlBSA_bsadb_bsaonly,ETHCD_xlBSA_bsadb_nonbsa)
colnames(ETHCD_xlBSA_bsadb)=c("ETHCD_xlBSA_bsadb_Total","ETHCD_xlBSA_bsadb_bsacrosslinks","ETHCD_xlBSA_bsadb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(ETHCD1_xlBSA_bsadb_bsaonly,ETHCD2_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(ETHCD1_xlBSA_bsadb_bsaonly,ETHCD3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(ETHCD2_xlBSA_bsadb_bsaonly,ETHCD3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB")))
tiff('ETHCD_xlBSA_bsadb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(ETHCD1_xlBSA_bsadb_bsaonly),area2=nrow(ETHCD2_xlBSA_bsadb_bsaonly),area3=nrow(ETHCD3_xlBSA_bsadb_bsaonly),
n12=nrow(merge(ETHCD1_xlBSA_bsadb_bsaonly,ETHCD2_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(ETHCD1_xlBSA_bsadb_bsaonly,ETHCD3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(ETHCD2_xlBSA_bsadb_bsaonly,ETHCD3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(ETHCD_xlBSA_bsadb_dup)-(nrow(ETHCD1_xlBSA_bsadb_bsaonly)+nrow(ETHCD2_xlBSA_bsadb_bsaonly)+nrow(ETHCD3_xlBSA_bsadb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("ETHCD_xlBSA_bsadb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()

##HCDMS3 cross-linked BSA BSA database
HCDMS31_xlBSA_bsadb=read.table("170629_xlBSA1_bsadb_CIDMS2_HCDMS3_Crosslinks.txt", fill=TRUE, header=TRUE)
HCDMS32_xlBSA_bsadb=read.table("170629_xlBSA2_bsadb_CIDMS2_HCDMS3_Crosslinks.txt", fill=TRUE, header=TRUE)
HCDMS33_xlBSA_bsadb=read.table("170629_xlBSA3_bsadb_CIDMS2_HCDMS3_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
HCDMS31_xlBSA_bsadb$Replicate=c(1)
HCDMS32_xlBSA_bsadb$Replicate=c(2)
HCDMS33_xlBSA_bsadb$Replicate=c(3)
##Compile all replicates
HCDMS3_xlBSA_bsadb=rbind(HCDMS31_xlBSA_bsadb,HCDMS32_xlBSA_bsadb,HCDMS33_xlBSA_bsadb)
write.table(HCDMS3_xlBSA_bsadb,"20170830_HCDMS3_xlBSA_bsadb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
HCDMS3_xlBSA_bsadb$PositionA=pmax(HCDMS3_xlBSA_bsadb$Position.A,HCDMS3_xlBSA_bsadb$Position.B)
HCDMS3_xlBSA_bsadb$PositionB=pmin(HCDMS3_xlBSA_bsadb$Position.A,HCDMS3_xlBSA_bsadb$Position.B)
##Keep BSA cross-links, duplicates kept
HCDMS3_xlBSA_bsadb_bsaonly<- HCDMS3_xlBSA_bsadb[grepl("P02769", HCDMS3_xlBSA_bsadb[["Accession.A"]]) & grepl("P02769", HCDMS3_xlBSA_bsadb[["Accession.B"]]), ]
write.table(HCDMS3_xlBSA_bsadb_bsaonly,"20170830_HCDMS3_xlBSA_bsadb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
HCDMS3_xlBSA_bsadb_nonbsa<- HCDMS3_xlBSA_bsadb[!grepl("P02769", HCDMS3_xlBSA_bsadb[["Accession.A"]]) | !grepl("P02769", HCDMS3_xlBSA_bsadb[["Accession.B"]]), ]
write.table(HCDMS3_xlBSA_bsadb_nonbsa,"20170830_HCDMS3_xlBSA_bsadb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
HCDMS3_xlBSA_bsadb_dup=HCDMS3_xlBSA_bsadb_bsaonly[!duplicated(HCDMS3_xlBSA_bsadb_bsaonly[c("PositionA","PositionB")]),]
write.table(HCDMS3_xlBSA_bsadb_dup,"20170830_HCDMS3_xlBSA_bsadb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
HCDMS3_xlBSA_bsadb_dup$Accession.A=as.character(HCDMS3_xlBSA_bsadb_dup$Accession.A)
HCDMS3_xlBSA_bsadb_dup$Accession.A<-replace(HCDMS3_xlBSA_bsadb_dup$Accession.A,HCDMS3_xlBSA_bsadb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
HCDMS3_xlBSA_bsadb_dup$Accession.B=as.character(HCDMS3_xlBSA_bsadb_dup$Accession.B)
HCDMS3_xlBSA_bsadb_dup$Accession.B<-replace(HCDMS3_xlBSA_bsadb_dup$Accession.B,HCDMS3_xlBSA_bsadb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
HCDMS3_xlBSA_bsadb_xinet<-data.frame(HCDMS3_xlBSA_bsadb_dup$Max.XlinkX.Score,as.factor(HCDMS3_xlBSA_bsadb_dup$Accession.A),as.factor(HCDMS3_xlBSA_bsadb_dup$Accession.B),
HCDMS3_xlBSA_bsadb_dup$PositionA,HCDMS3_xlBSA_bsadb_dup$PositionB)
colnames(HCDMS3_xlBSA_bsadb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(HCDMS3_xlBSA_bsadb_xinet,"20170830_HCDMS3_xlBSA_bsadb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
HCDMS31_xlBSA_bsadb$PositionA=pmax(HCDMS31_xlBSA_bsadb$Position.A,HCDMS31_xlBSA_bsadb$Position.B)
HCDMS31_xlBSA_bsadb$PositionB=pmin(HCDMS31_xlBSA_bsadb$Position.A,HCDMS31_xlBSA_bsadb$Position.B)
HCDMS32_xlBSA_bsadb$PositionA=pmax(HCDMS32_xlBSA_bsadb$Position.A,HCDMS32_xlBSA_bsadb$Position.B)
HCDMS32_xlBSA_bsadb$PositionB=pmin(HCDMS32_xlBSA_bsadb$Position.A,HCDMS32_xlBSA_bsadb$Position.B)
HCDMS33_xlBSA_bsadb$PositionA=pmax(HCDMS33_xlBSA_bsadb$Position.A,HCDMS33_xlBSA_bsadb$Position.B)
HCDMS33_xlBSA_bsadb$PositionB=pmin(HCDMS33_xlBSA_bsadb$Position.A,HCDMS33_xlBSA_bsadb$Position.B)
#Filter XlinkX Score >80
HCDMS31_xlBSA_bsadb=HCDMS31_xlBSA_bsadb[HCDMS31_xlBSA_bsadb[,"Max.XlinkX.Score"]>s12,]
HCDMS32_xlBSA_bsadb=HCDMS32_xlBSA_bsadb[HCDMS32_xlBSA_bsadb[,"Max.XlinkX.Score"]>s12,]
HCDMS33_xlBSA_bsadb=HCDMS33_xlBSA_bsadb[HCDMS33_xlBSA_bsadb[,"Max.XlinkX.Score"]>s12,]
#Remove Duplicates
HCDMS31_xlBSA_bsadb_dup=HCDMS31_xlBSA_bsadb[!duplicated(HCDMS31_xlBSA_bsadb[c("PositionA","PositionB")]),]
HCDMS32_xlBSA_bsadb_dup=HCDMS32_xlBSA_bsadb[!duplicated(HCDMS32_xlBSA_bsadb[c("PositionA","PositionB")]),]
HCDMS33_xlBSA_bsadb_dup=HCDMS33_xlBSA_bsadb[!duplicated(HCDMS33_xlBSA_bsadb[c("PositionA","PositionB")]),]
#HCDMS3 Number of BSA crosslinks
HCDMS31_xlBSA_bsadb_bsaonly<- HCDMS31_xlBSA_bsadb_dup[grepl("P02769", HCDMS31_xlBSA_bsadb_dup[["Accession.A"]]) & grepl("P02769", HCDMS31_xlBSA_bsadb_dup[["Accession.B"]]), ]
HCDMS32_xlBSA_bsadb_bsaonly<- HCDMS32_xlBSA_bsadb_dup[grepl("P02769", HCDMS32_xlBSA_bsadb_dup[["Accession.A"]]) & grepl("P02769", HCDMS32_xlBSA_bsadb_dup[["Accession.B"]]), ]
HCDMS33_xlBSA_bsadb_bsaonly<- HCDMS33_xlBSA_bsadb_dup[grepl("P02769", HCDMS33_xlBSA_bsadb_dup[["Accession.A"]]) & grepl("P02769", HCDMS33_xlBSA_bsadb_dup[["Accession.B"]]), ]
#HCDMS3 Number of all crosslinks, no duplicates
HCDMS3_xlBSA_bsadb_crosslinks=c(nrow(HCDMS31_xlBSA_bsadb_dup),nrow(HCDMS32_xlBSA_bsadb_dup),nrow(HCDMS33_xlBSA_bsadb_dup))
#HCDMS3 Number of BSA crosslinks, no duplicates
HCDMS3_xlBSA_bsadb_bsaonly=c(nrow(HCDMS31_xlBSA_bsadb_bsaonly),nrow(HCDMS32_xlBSA_bsadb_bsaonly),nrow(HCDMS33_xlBSA_bsadb_bsaonly))
#HCDMS3 Number of non-BSA crosslinks
HCDMS3_xlBSA_bsadb_nonbsa=(HCDMS3_xlBSA_bsadb_crosslinks-HCDMS3_xlBSA_bsadb_bsaonly)
#Output for Number of crosslinks
HCDMS3_xlBSA_bsadb=data.frame(HCDMS3_xlBSA_bsadb_crosslinks,HCDMS3_xlBSA_bsadb_bsaonly,HCDMS3_xlBSA_bsadb_nonbsa)
colnames(HCDMS3_xlBSA_bsadb)=c("HCDMS3_xlBSA_bsadb_Total","HCDMS3_xlBSA_bsadb_bsacrosslinks","HCDMS3_xlBSA_bsadb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(HCDMS31_xlBSA_bsadb_bsaonly,HCDMS32_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(HCDMS31_xlBSA_bsadb_bsaonly,HCDMS33_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(HCDMS32_xlBSA_bsadb_bsaonly,HCDMS33_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB")))
tiff('HCDMS3_xlBSA_bsadb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(HCDMS31_xlBSA_bsadb_bsaonly),area2=nrow(HCDMS32_xlBSA_bsadb_bsaonly),area3=nrow(HCDMS33_xlBSA_bsadb_bsaonly),
n12=nrow(merge(HCDMS31_xlBSA_bsadb_bsaonly,HCDMS32_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(HCDMS31_xlBSA_bsadb_bsaonly,HCDMS33_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(HCDMS32_xlBSA_bsadb_bsaonly,HCDMS33_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(HCDMS3_xlBSA_bsadb_dup)-(nrow(HCDMS31_xlBSA_bsadb_bsaonly)+nrow(HCDMS32_xlBSA_bsadb_bsaonly)+nrow(HCDMS33_xlBSA_bsadb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("HCDMS3_xlBSA_bsadb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()


##ETD cross-linked BSA+hek293t BSA database
ETD1_xlBSA_293_bsadb=read.table("170629_xlBSA1_bsadb_CIDMS2_ETDMS2_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
ETD2_xlBSA_293_bsadb=read.table("170629_xlBSA2_bsadb_CIDMS2_ETDMS2_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
ETD3_xlBSA_293_bsadb=read.table("170629_xlBSA3_bsadb_CIDMS2_ETDMS2_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
ETD1_xlBSA_293_bsadb$Replicate=c(1)
ETD2_xlBSA_293_bsadb$Replicate=c(2)
ETD3_xlBSA_293_bsadb$Replicate=c(3)
##Compile all replicates
ETD_xlBSA_293_bsadb=rbind(ETD1_xlBSA_293_bsadb,ETD2_xlBSA_293_bsadb,ETD3_xlBSA_293_bsadb)
write.table(ETD_xlBSA_293_bsadb,"20170830_ETDMS2_xlBSA_293_bsadb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
ETD_xlBSA_293_bsadb$PositionA=pmax(ETD_xlBSA_293_bsadb$Position.A,ETD_xlBSA_293_bsadb$Position.B)
ETD_xlBSA_293_bsadb$PositionB=pmin(ETD_xlBSA_293_bsadb$Position.A,ETD_xlBSA_293_bsadb$Position.B)
##Keep BSA cross-links, duplicates kept
ETD_xlBSA_293_bsadb_bsaonly<- ETD_xlBSA_293_bsadb[grepl("P02769", ETD_xlBSA_293_bsadb[["Accession.A"]]) & grepl("P02769", ETD_xlBSA_293_bsadb[["Accession.B"]]), ]
write.table(ETD_xlBSA_293_bsadb_bsaonly,"20170830_ETDMS2_xlBSA_293_bsadb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
ETD_xlBSA_293_bsadb_nonbsa<- ETD_xlBSA_293_bsadb[!grepl("P02769", ETD_xlBSA_293_bsadb[["Accession.A"]]) | !grepl("P02769", ETD_xlBSA_293_bsadb[["Accession.B"]]), ]
write.table(ETD_xlBSA_293_bsadb_nonbsa,"20170830_ETDMS2_xlBSA_293_bsadb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
ETD_xlBSA_293_bsadb_dup=ETD_xlBSA_293_bsadb_bsaonly[!duplicated(ETD_xlBSA_293_bsadb_bsaonly[c("PositionA","PositionB")]),]
write.table(ETD_xlBSA_293_bsadb_dup,"20170830_ETDMS2_xlBSA_293_bsadb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
ETD_xlBSA_293_bsadb_dup$Accession.A=as.character(ETD_xlBSA_293_bsadb_dup$Accession.A)
ETD_xlBSA_293_bsadb_dup$Accession.A<-replace(ETD_xlBSA_293_bsadb_dup$Accession.A,ETD_xlBSA_293_bsadb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
ETD_xlBSA_293_bsadb_dup$Accession.B=as.character(ETD_xlBSA_293_bsadb_dup$Accession.B)
ETD_xlBSA_293_bsadb_dup$Accession.B<-replace(ETD_xlBSA_293_bsadb_dup$Accession.B,ETD_xlBSA_293_bsadb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
ETD_xlBSA_293_bsadb_xinet<-data.frame(ETD_xlBSA_293_bsadb_dup$Max.XlinkX.Score,as.factor(ETD_xlBSA_293_bsadb_dup$Accession.A),as.factor(ETD_xlBSA_293_bsadb_dup$Accession.B),
ETD_xlBSA_293_bsadb_dup$PositionA,ETD_xlBSA_293_bsadb_dup$PositionB)
colnames(ETD_xlBSA_293_bsadb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(ETD_xlBSA_293_bsadb_xinet,"20170830_ETDMS2_xlBSA_293_bsadb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
ETD1_xlBSA_293_bsadb$PositionA=pmax(ETD1_xlBSA_293_bsadb$Position.A,ETD1_xlBSA_293_bsadb$Position.B)
ETD1_xlBSA_293_bsadb$PositionB=pmin(ETD1_xlBSA_293_bsadb$Position.A,ETD1_xlBSA_293_bsadb$Position.B)
ETD2_xlBSA_293_bsadb$PositionA=pmax(ETD2_xlBSA_293_bsadb$Position.A,ETD2_xlBSA_293_bsadb$Position.B)
ETD2_xlBSA_293_bsadb$PositionB=pmin(ETD2_xlBSA_293_bsadb$Position.A,ETD2_xlBSA_293_bsadb$Position.B)
ETD3_xlBSA_293_bsadb$PositionA=pmax(ETD3_xlBSA_293_bsadb$Position.A,ETD3_xlBSA_293_bsadb$Position.B)
ETD3_xlBSA_293_bsadb$PositionB=pmin(ETD3_xlBSA_293_bsadb$Position.A,ETD3_xlBSA_293_bsadb$Position.B)
#Filter XlinkX Score >80
ETD1_xlBSA_293_bsadb=ETD1_xlBSA_293_bsadb[ETD1_xlBSA_293_bsadb[,"Max.XlinkX.Score"]>s13,]
ETD2_xlBSA_293_bsadb=ETD2_xlBSA_293_bsadb[ETD2_xlBSA_293_bsadb[,"Max.XlinkX.Score"]>s13,]
ETD3_xlBSA_293_bsadb=ETD3_xlBSA_293_bsadb[ETD3_xlBSA_293_bsadb[,"Max.XlinkX.Score"]>s13,]
#Remove Duplicates
ETD1_xlBSA_293_bsadb_dup=ETD1_xlBSA_293_bsadb[!duplicated(ETD1_xlBSA_293_bsadb[c("PositionA","PositionB")]),]
ETD2_xlBSA_293_bsadb_dup=ETD2_xlBSA_293_bsadb[!duplicated(ETD2_xlBSA_293_bsadb[c("PositionA","PositionB")]),]
ETD3_xlBSA_293_bsadb_dup=ETD3_xlBSA_293_bsadb[!duplicated(ETD3_xlBSA_293_bsadb[c("PositionA","PositionB")]),]
#ETD Number of BSA crosslinks
ETD1_xlBSA_293_bsadb_bsaonly<- ETD1_xlBSA_293_bsadb_dup[grepl("P02769", ETD1_xlBSA_293_bsadb_dup[["Accession.A"]]) & grepl("P02769", ETD1_xlBSA_293_bsadb_dup[["Accession.B"]]), ]
ETD2_xlBSA_293_bsadb_bsaonly<- ETD2_xlBSA_293_bsadb_dup[grepl("P02769", ETD2_xlBSA_293_bsadb_dup[["Accession.A"]]) & grepl("P02769", ETD2_xlBSA_293_bsadb_dup[["Accession.B"]]), ]
ETD3_xlBSA_293_bsadb_bsaonly<- ETD3_xlBSA_293_bsadb_dup[grepl("P02769", ETD3_xlBSA_293_bsadb_dup[["Accession.A"]]) & grepl("P02769", ETD3_xlBSA_293_bsadb_dup[["Accession.B"]]), ]
#ETD Number of all crosslinks, no duplicates
ETD_xlBSA_293_bsadb_crosslinks=c(nrow(ETD1_xlBSA_293_bsadb_dup),nrow(ETD2_xlBSA_293_bsadb_dup),nrow(ETD3_xlBSA_293_bsadb_dup))
#ETD Number of BSA crosslinks, no duplicates
ETD_xlBSA_293_bsadb_bsaonly=c(nrow(ETD1_xlBSA_293_bsadb_bsaonly),nrow(ETD2_xlBSA_293_bsadb_bsaonly),nrow(ETD3_xlBSA_293_bsadb_bsaonly))
#ETD Number of non-BSA crosslinks
ETD_xlBSA_293_bsadb_nonbsa=(ETD_xlBSA_293_bsadb_crosslinks-ETD_xlBSA_293_bsadb_bsaonly)
#Output for Number of crosslinks
ETD_xlBSA_293_bsadb=data.frame(ETD_xlBSA_293_bsadb_crosslinks,ETD_xlBSA_293_bsadb_bsaonly,ETD_xlBSA_293_bsadb_nonbsa)
colnames(ETD_xlBSA_293_bsadb)=c("ETD_xlBSA_293_bsadb_Total","ETD_xlBSA_293_bsadb_bsacrosslinks","ETD_xlBSA_293_bsadb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(ETD1_xlBSA_293_bsadb_bsaonly,ETD2_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(ETD1_xlBSA_293_bsadb_bsaonly,ETD3_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(ETD2_xlBSA_293_bsadb_bsaonly,ETD3_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB")))
tiff('ETD_xlBSA_293_bsadb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(ETD1_xlBSA_293_bsadb_bsaonly),area2=nrow(ETD2_xlBSA_293_bsadb_bsaonly),area3=nrow(ETD3_xlBSA_293_bsadb_bsaonly),
n12=nrow(merge(ETD1_xlBSA_293_bsadb_bsaonly,ETD2_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(ETD1_xlBSA_293_bsadb_bsaonly,ETD3_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(ETD2_xlBSA_293_bsadb_bsaonly,ETD3_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(ETD_xlBSA_293_bsadb_dup)-(nrow(ETD1_xlBSA_293_bsadb_bsaonly)+nrow(ETD2_xlBSA_293_bsadb_bsaonly)+nrow(ETD3_xlBSA_293_bsadb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("ETD_xlBSA_293_bsadb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()

##HCD cross-linked BSA+hek293t BSA database
HCD1_xlBSA_293_bsadb=read.table("170629_xlBSA1_bsadb_CIDMS2_HCDMS2_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
HCD2_xlBSA_293_bsadb=read.table("170629_xlBSA1_bsadb_CIDMS2_HCDMS2_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
HCD3_xlBSA_293_bsadb=read.table("170629_xlBSA3_bsadb_CIDMS2_HCDMS2_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
HCD1_xlBSA_293_bsadb$Replicate=c(1)
HCD2_xlBSA_293_bsadb$Replicate=c(2)
HCD3_xlBSA_293_bsadb$Replicate=c(3)
##Compile all replicates
HCD_xlBSA_293_bsadb=rbind(HCD1_xlBSA_293_bsadb,HCD2_xlBSA_293_bsadb,HCD3_xlBSA_293_bsadb)
write.table(HCD_xlBSA_293_bsadb,"20170830_HCDMS2_xlBSA_293_bsadb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
HCD_xlBSA_293_bsadb$PositionA=pmax(HCD_xlBSA_293_bsadb$Position.A,HCD_xlBSA_293_bsadb$Position.B)
HCD_xlBSA_293_bsadb$PositionB=pmin(HCD_xlBSA_293_bsadb$Position.A,HCD_xlBSA_293_bsadb$Position.B)
##Keep BSA cross-links, duplicates kept
HCD_xlBSA_293_bsadb_bsaonly<- HCD_xlBSA_293_bsadb[grepl("P02769", HCD_xlBSA_293_bsadb[["Accession.A"]]) & grepl("P02769", HCD_xlBSA_293_bsadb[["Accession.B"]]), ]
write.table(HCD_xlBSA_293_bsadb_bsaonly,"20170830_HCDMS2_xlBSA_293_bsadb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
HCD_xlBSA_293_bsadb_nonbsa<- HCD_xlBSA_293_bsadb[!grepl("P02769", HCD_xlBSA_293_bsadb[["Accession.A"]]) | !grepl("P02769", HCD_xlBSA_293_bsadb[["Accession.B"]]), ]
write.table(HCD_xlBSA_293_bsadb_nonbsa,"20170830_HCDMS2_xlBSA_293_bsadb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
HCD_xlBSA_293_bsadb_dup=HCD_xlBSA_293_bsadb_bsaonly[!duplicated(HCD_xlBSA_293_bsadb_bsaonly[c("PositionA","PositionB")]),]
write.table(HCD_xlBSA_293_bsadb_dup,"20170830_HCDMS2_xlBSA_293_bsadb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
HCD_xlBSA_293_bsadb_dup$Accession.A=as.character(HCD_xlBSA_293_bsadb_dup$Accession.A)
HCD_xlBSA_293_bsadb_dup$Accession.A<-replace(HCD_xlBSA_293_bsadb_dup$Accession.A,HCD_xlBSA_293_bsadb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
HCD_xlBSA_293_bsadb_dup$Accession.B=as.character(HCD_xlBSA_293_bsadb_dup$Accession.B)
HCD_xlBSA_293_bsadb_dup$Accession.B<-replace(HCD_xlBSA_293_bsadb_dup$Accession.B,HCD_xlBSA_293_bsadb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
HCD_xlBSA_293_bsadb_xinet<-data.frame(HCD_xlBSA_293_bsadb_dup$Max.XlinkX.Score,as.factor(HCD_xlBSA_293_bsadb_dup$Accession.A),as.factor(HCD_xlBSA_293_bsadb_dup$Accession.B),
HCD_xlBSA_293_bsadb_dup$PositionA,HCD_xlBSA_293_bsadb_dup$PositionB)
colnames(HCD_xlBSA_293_bsadb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(HCD_xlBSA_293_bsadb_xinet,"20170830_HCDMS2_xlBSA_293_bsadb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
HCD1_xlBSA_293_bsadb$PositionA=pmax(HCD1_xlBSA_293_bsadb$Position.A,HCD1_xlBSA_293_bsadb$Position.B)
HCD1_xlBSA_293_bsadb$PositionB=pmin(HCD1_xlBSA_293_bsadb$Position.A,HCD1_xlBSA_293_bsadb$Position.B)
HCD2_xlBSA_293_bsadb$PositionA=pmax(HCD2_xlBSA_293_bsadb$Position.A,HCD2_xlBSA_293_bsadb$Position.B)
HCD2_xlBSA_293_bsadb$PositionB=pmin(HCD2_xlBSA_293_bsadb$Position.A,HCD2_xlBSA_293_bsadb$Position.B)
HCD3_xlBSA_293_bsadb$PositionA=pmax(HCD3_xlBSA_293_bsadb$Position.A,HCD3_xlBSA_293_bsadb$Position.B)
HCD3_xlBSA_293_bsadb$PositionB=pmin(HCD3_xlBSA_293_bsadb$Position.A,HCD3_xlBSA_293_bsadb$Position.B)
#Filter XlinkX Score >80
HCD1_xlBSA_293_bsadb=HCD1_xlBSA_293_bsadb[HCD1_xlBSA_293_bsadb[,"Max.XlinkX.Score"]>s14,]
HCD2_xlBSA_293_bsadb=HCD2_xlBSA_293_bsadb[HCD2_xlBSA_293_bsadb[,"Max.XlinkX.Score"]>s14,]
HCD3_xlBSA_293_bsadb=HCD3_xlBSA_293_bsadb[HCD3_xlBSA_293_bsadb[,"Max.XlinkX.Score"]>s14,]
#Remove Duplicates
HCD1_xlBSA_293_bsadb_dup=HCD1_xlBSA_293_bsadb[!duplicated(HCD1_xlBSA_293_bsadb[c("PositionA","PositionB")]),]
HCD2_xlBSA_293_bsadb_dup=HCD2_xlBSA_293_bsadb[!duplicated(HCD2_xlBSA_293_bsadb[c("PositionA","PositionB")]),]
HCD3_xlBSA_293_bsadb_dup=HCD3_xlBSA_293_bsadb[!duplicated(HCD3_xlBSA_293_bsadb[c("PositionA","PositionB")]),]
#HCD Number of BSA crosslinks
HCD1_xlBSA_293_bsadb_bsaonly<- HCD1_xlBSA_293_bsadb_dup[grepl("P02769", HCD1_xlBSA_293_bsadb_dup[["Accession.A"]]) & grepl("P02769", HCD1_xlBSA_293_bsadb_dup[["Accession.B"]]), ]
HCD2_xlBSA_293_bsadb_bsaonly<- HCD2_xlBSA_293_bsadb_dup[grepl("P02769", HCD2_xlBSA_293_bsadb_dup[["Accession.A"]]) & grepl("P02769", HCD2_xlBSA_293_bsadb_dup[["Accession.B"]]), ]
HCD3_xlBSA_293_bsadb_bsaonly<- HCD3_xlBSA_293_bsadb_dup[grepl("P02769", HCD3_xlBSA_293_bsadb_dup[["Accession.A"]]) & grepl("P02769", HCD3_xlBSA_293_bsadb_dup[["Accession.B"]]), ]
#HCD Number of all crosslinks, no duplicates
HCD_xlBSA_293_bsadb_crosslinks=c(nrow(HCD1_xlBSA_293_bsadb_dup),nrow(HCD2_xlBSA_293_bsadb_dup),nrow(HCD3_xlBSA_293_bsadb_dup))
#HCD Number of BSA crosslinks, no duplicates
HCD_xlBSA_293_bsadb_bsaonly=c(nrow(HCD1_xlBSA_293_bsadb_bsaonly),nrow(HCD2_xlBSA_293_bsadb_bsaonly),nrow(HCD3_xlBSA_293_bsadb_bsaonly))
#HCD Number of non-BSA crosslinks
HCD_xlBSA_293_bsadb_nonbsa=(HCD_xlBSA_293_bsadb_crosslinks-HCD_xlBSA_293_bsadb_bsaonly)
#Output for Number of crosslinks
HCD_xlBSA_293_bsadb=data.frame(HCD_xlBSA_293_bsadb_crosslinks,HCD_xlBSA_293_bsadb_bsaonly,HCD_xlBSA_293_bsadb_nonbsa)
colnames(HCD_xlBSA_293_bsadb)=c("HCD_xlBSA_293_bsadb_Total","HCD_xlBSA_293_bsadb_bsacrosslinks","HCD_xlBSA_293_bsadb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(HCD1_xlBSA_293_bsadb_bsaonly,HCD2_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(HCD1_xlBSA_293_bsadb_bsaonly,HCD3_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(HCD2_xlBSA_293_bsadb_bsaonly,HCD3_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB")))
tiff('HCD_xlBSA_293_bsadb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(HCD1_xlBSA_293_bsadb_bsaonly),area2=nrow(HCD2_xlBSA_293_bsadb_bsaonly),area3=nrow(HCD3_xlBSA_293_bsadb_bsaonly),
n12=nrow(merge(HCD1_xlBSA_293_bsadb_bsaonly,HCD2_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(HCD1_xlBSA_293_bsadb_bsaonly,HCD3_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(HCD2_xlBSA_293_bsadb_bsaonly,HCD3_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(HCD_xlBSA_293_bsadb_dup)-(nrow(HCD1_xlBSA_293_bsadb_bsaonly)+nrow(HCD2_xlBSA_293_bsadb_bsaonly)+nrow(HCD3_xlBSA_293_bsadb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("HCD_xlBSA_293_bsadb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()

##ETHCD cross-linked BSA+hek293t BSA database
ETHCD1_xlBSA_293_bsadb=read.table("170629_xlBSA1_bsadb_CIDMS2_EThcDMS2_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
ETHCD2_xlBSA_293_bsadb=read.table("170629_xlBSA2_bsadb_CIDMS2_EThcDMS2_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
ETHCD3_xlBSA_293_bsadb=read.table("170629_xlBSA3_bsadb_CIDMS2_EThcDMS2_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
ETHCD1_xlBSA_293_bsadb$Replicate=c(1)
ETHCD2_xlBSA_293_bsadb$Replicate=c(2)
ETHCD3_xlBSA_293_bsadb$Replicate=c(3)
##Compile all replicates
ETHCD_xlBSA_293_bsadb=rbind(ETHCD1_xlBSA_293_bsadb,ETHCD2_xlBSA_293_bsadb,ETHCD3_xlBSA_293_bsadb)
write.table(ETHCD_xlBSA_293_bsadb,"20170830_ETHCDMS2_xlBSA_293_bsadb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
ETHCD_xlBSA_293_bsadb$PositionA=pmax(ETHCD_xlBSA_293_bsadb$Position.A,ETHCD_xlBSA_293_bsadb$Position.B)
ETHCD_xlBSA_293_bsadb$PositionB=pmin(ETHCD_xlBSA_293_bsadb$Position.A,ETHCD_xlBSA_293_bsadb$Position.B)
##Keep BSA cross-links, duplicates kept
ETHCD_xlBSA_293_bsadb_bsaonly<- ETHCD_xlBSA_293_bsadb[grepl("P02769", ETHCD_xlBSA_293_bsadb[["Accession.A"]]) & grepl("P02769", ETHCD_xlBSA_293_bsadb[["Accession.B"]]), ]
write.table(ETHCD_xlBSA_293_bsadb_bsaonly,"20170830_ETHCDMS2_xlBSA_293_bsadb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
ETHCD_xlBSA_293_bsadb_nonbsa<- ETHCD_xlBSA_293_bsadb[!grepl("P02769", ETHCD_xlBSA_293_bsadb[["Accession.A"]]) | !grepl("P02769", ETHCD_xlBSA_293_bsadb[["Accession.B"]]), ]
write.table(ETHCD_xlBSA_293_bsadb_nonbsa,"20170830_ETHCDMS2_xlBSA_293_bsadb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
ETHCD_xlBSA_293_bsadb_dup=ETHCD_xlBSA_293_bsadb_bsaonly[!duplicated(ETHCD_xlBSA_293_bsadb_bsaonly[c("PositionA","PositionB")]),]
write.table(ETHCD_xlBSA_293_bsadb_dup,"20170830_ETHCDMS2_xlBSA_293_bsadb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
ETHCD_xlBSA_293_bsadb_dup$Accession.A=as.character(ETHCD_xlBSA_293_bsadb_dup$Accession.A)
ETHCD_xlBSA_293_bsadb_dup$Accession.A<-replace(ETHCD_xlBSA_293_bsadb_dup$Accession.A,ETHCD_xlBSA_293_bsadb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
ETHCD_xlBSA_293_bsadb_dup$Accession.B=as.character(ETHCD_xlBSA_293_bsadb_dup$Accession.B)
ETHCD_xlBSA_293_bsadb_dup$Accession.B<-replace(ETHCD_xlBSA_293_bsadb_dup$Accession.B,ETHCD_xlBSA_293_bsadb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
ETHCD_xlBSA_293_bsadb_xinet<-data.frame(ETHCD_xlBSA_293_bsadb_dup$Max.XlinkX.Score,as.factor(ETHCD_xlBSA_293_bsadb_dup$Accession.A),as.factor(ETHCD_xlBSA_293_bsadb_dup$Accession.B),
ETHCD_xlBSA_293_bsadb_dup$PositionA,ETHCD_xlBSA_293_bsadb_dup$PositionB)
colnames(ETHCD_xlBSA_293_bsadb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(ETHCD_xlBSA_293_bsadb_xinet,"20170830_ETHCDMS2_xlBSA_293_bsadb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
ETHCD1_xlBSA_293_bsadb$PositionA=pmax(ETHCD1_xlBSA_293_bsadb$Position.A,ETHCD1_xlBSA_293_bsadb$Position.B)
ETHCD1_xlBSA_293_bsadb$PositionB=pmin(ETHCD1_xlBSA_293_bsadb$Position.A,ETHCD1_xlBSA_293_bsadb$Position.B)
ETHCD2_xlBSA_293_bsadb$PositionA=pmax(ETHCD2_xlBSA_293_bsadb$Position.A,ETHCD2_xlBSA_293_bsadb$Position.B)
ETHCD2_xlBSA_293_bsadb$PositionB=pmin(ETHCD2_xlBSA_293_bsadb$Position.A,ETHCD2_xlBSA_293_bsadb$Position.B)
ETHCD3_xlBSA_293_bsadb$PositionA=pmax(ETHCD3_xlBSA_293_bsadb$Position.A,ETHCD3_xlBSA_293_bsadb$Position.B)
ETHCD3_xlBSA_293_bsadb$PositionB=pmin(ETHCD3_xlBSA_293_bsadb$Position.A,ETHCD3_xlBSA_293_bsadb$Position.B)
#Filter XlinkX Score >80
ETHCD1_xlBSA_293_bsadb=ETHCD1_xlBSA_293_bsadb[ETHCD1_xlBSA_293_bsadb[,"Max.XlinkX.Score"]>s15,]
ETHCD2_xlBSA_293_bsadb=ETHCD2_xlBSA_293_bsadb[ETHCD2_xlBSA_293_bsadb[,"Max.XlinkX.Score"]>s15,]
ETHCD3_xlBSA_293_bsadb=ETHCD3_xlBSA_293_bsadb[ETHCD3_xlBSA_293_bsadb[,"Max.XlinkX.Score"]>s15,]
#Remove Duplicates
ETHCD1_xlBSA_293_bsadb_dup=ETHCD1_xlBSA_293_bsadb[!duplicated(ETHCD1_xlBSA_293_bsadb[c("PositionA","PositionB")]),]
ETHCD2_xlBSA_293_bsadb_dup=ETHCD2_xlBSA_293_bsadb[!duplicated(ETHCD2_xlBSA_293_bsadb[c("PositionA","PositionB")]),]
ETHCD3_xlBSA_293_bsadb_dup=ETHCD3_xlBSA_293_bsadb[!duplicated(ETHCD3_xlBSA_293_bsadb[c("PositionA","PositionB")]),]
#ETHCD Number of BSA crosslinks
ETHCD1_xlBSA_293_bsadb_bsaonly<- ETHCD1_xlBSA_293_bsadb_dup[grepl("P02769", ETHCD1_xlBSA_293_bsadb_dup[["Accession.A"]]) & grepl("P02769", ETHCD1_xlBSA_293_bsadb_dup[["Accession.B"]]), ]
ETHCD2_xlBSA_293_bsadb_bsaonly<- ETHCD2_xlBSA_293_bsadb_dup[grepl("P02769", ETHCD2_xlBSA_293_bsadb_dup[["Accession.A"]]) & grepl("P02769", ETHCD2_xlBSA_293_bsadb_dup[["Accession.B"]]), ]
ETHCD3_xlBSA_293_bsadb_bsaonly<- ETHCD3_xlBSA_293_bsadb_dup[grepl("P02769", ETHCD3_xlBSA_293_bsadb_dup[["Accession.A"]]) & grepl("P02769", ETHCD3_xlBSA_293_bsadb_dup[["Accession.B"]]), ]
#ETHCD Number of all crosslinks, no duplicates
ETHCD_xlBSA_293_bsadb_crosslinks=c(nrow(ETHCD1_xlBSA_293_bsadb_dup),nrow(ETHCD2_xlBSA_293_bsadb_dup),nrow(ETHCD3_xlBSA_293_bsadb_dup))
#ETHCD Number of BSA crosslinks, no duplicates
ETHCD_xlBSA_293_bsadb_bsaonly=c(nrow(ETHCD1_xlBSA_293_bsadb_bsaonly),nrow(ETHCD2_xlBSA_293_bsadb_bsaonly),nrow(ETHCD3_xlBSA_293_bsadb_bsaonly))
#ETHCD Number of non-BSA crosslinks
ETHCD_xlBSA_293_bsadb_nonbsa=(ETHCD_xlBSA_293_bsadb_crosslinks-ETHCD_xlBSA_293_bsadb_bsaonly)
#Output for Number of crosslinks
ETHCD_xlBSA_293_bsadb=data.frame(ETHCD_xlBSA_293_bsadb_crosslinks,ETHCD_xlBSA_293_bsadb_bsaonly,ETHCD_xlBSA_293_bsadb_nonbsa)
colnames(ETHCD_xlBSA_293_bsadb)=c("ETHCD_xlBSA_293_bsadb_Total","ETHCD_xlBSA_293_bsadb_bsacrosslinks","ETHCD_xlBSA_293_bsadb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(ETHCD1_xlBSA_293_bsadb_bsaonly,ETHCD2_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(ETHCD1_xlBSA_293_bsadb_bsaonly,ETHCD3_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(ETHCD2_xlBSA_293_bsadb_bsaonly,ETHCD3_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB")))
tiff('ETHCD_xlBSA_293_bsadb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(ETHCD1_xlBSA_293_bsadb_bsaonly),area2=nrow(ETHCD2_xlBSA_293_bsadb_bsaonly),area3=nrow(ETHCD3_xlBSA_293_bsadb_bsaonly),
n12=nrow(merge(ETHCD1_xlBSA_293_bsadb_bsaonly,ETHCD2_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(ETHCD1_xlBSA_293_bsadb_bsaonly,ETHCD3_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(ETHCD2_xlBSA_293_bsadb_bsaonly,ETHCD3_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(ETHCD_xlBSA_293_bsadb_dup)-(nrow(ETHCD1_xlBSA_293_bsadb_bsaonly)+nrow(ETHCD2_xlBSA_293_bsadb_bsaonly)+nrow(ETHCD3_xlBSA_293_bsadb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("ETHCD_xlBSA_293_bsadb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()

##HCDMS3 cross-linked BSA+hek293t BSA database
HCDMS31_xlBSA_293_bsadb=read.table("170629_xlBSA1_bsadb_CIDMS2_HCDMS3_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
HCDMS32_xlBSA_293_bsadb=read.table("170629_xlBSA2_bsadb_CIDMS2_HCDMS3_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
HCDMS33_xlBSA_293_bsadb=read.table("170629_xlBSA3_bsadb_CIDMS2_HCDMS3_hek293_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
HCDMS31_xlBSA_293_bsadb$Replicate=c(1)
HCDMS32_xlBSA_293_bsadb$Replicate=c(2)
HCDMS33_xlBSA_293_bsadb$Replicate=c(3)
##Compile all replicates
HCDMS3_xlBSA_293_bsadb=rbind(HCDMS31_xlBSA_293_bsadb,HCDMS32_xlBSA_293_bsadb,HCDMS33_xlBSA_293_bsadb)
write.table(HCDMS3_xlBSA_293_bsadb,"20170830_HCDMS3_xlBSA_293_bsadb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
HCDMS3_xlBSA_293_bsadb$PositionA=pmax(HCDMS3_xlBSA_293_bsadb$Position.A,HCDMS3_xlBSA_293_bsadb$Position.B)
HCDMS3_xlBSA_293_bsadb$PositionB=pmin(HCDMS3_xlBSA_293_bsadb$Position.A,HCDMS3_xlBSA_293_bsadb$Position.B)
##Keep BSA cross-links, duplicates kept
HCDMS3_xlBSA_293_bsadb_bsaonly<- HCDMS3_xlBSA_293_bsadb[grepl("P02769", HCDMS3_xlBSA_293_bsadb[["Accession.A"]]) & grepl("P02769", HCDMS3_xlBSA_293_bsadb[["Accession.B"]]), ]
write.table(HCDMS3_xlBSA_293_bsadb_bsaonly,"20170830_HCDMS3_xlBSA_293_bsadb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
HCDMS3_xlBSA_293_bsadb_nonbsa<- HCDMS3_xlBSA_293_bsadb[!grepl("P02769", HCDMS3_xlBSA_293_bsadb[["Accession.A"]]) | !grepl("P02769", HCDMS3_xlBSA_293_bsadb[["Accession.B"]]), ]
write.table(HCDMS3_xlBSA_293_bsadb_nonbsa,"20170830_HCDMS3_xlBSA_293_bsadb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
HCDMS3_xlBSA_293_bsadb_dup=HCDMS3_xlBSA_293_bsadb_bsaonly[!duplicated(HCDMS3_xlBSA_293_bsadb_bsaonly[c("PositionA","PositionB")]),]
write.table(HCDMS3_xlBSA_293_bsadb_dup,"20170830_HCDMS3_xlBSA_293_bsadb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
HCDMS3_xlBSA_293_bsadb_dup$Accession.A=as.character(HCDMS3_xlBSA_293_bsadb_dup$Accession.A)
HCDMS3_xlBSA_293_bsadb_dup$Accession.A<-replace(HCDMS3_xlBSA_293_bsadb_dup$Accession.A,HCDMS3_xlBSA_293_bsadb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
HCDMS3_xlBSA_293_bsadb_dup$Accession.B=as.character(HCDMS3_xlBSA_293_bsadb_dup$Accession.B)
HCDMS3_xlBSA_293_bsadb_dup$Accession.B<-replace(HCDMS3_xlBSA_293_bsadb_dup$Accession.B,HCDMS3_xlBSA_293_bsadb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
HCDMS3_xlBSA_293_bsadb_xinet<-data.frame(HCDMS3_xlBSA_293_bsadb_dup$Max.XlinkX.Score,as.factor(HCDMS3_xlBSA_293_bsadb_dup$Accession.A),as.factor(HCDMS3_xlBSA_293_bsadb_dup$Accession.B),
HCDMS3_xlBSA_293_bsadb_dup$PositionA,HCDMS3_xlBSA_293_bsadb_dup$PositionB)
colnames(HCDMS3_xlBSA_293_bsadb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(HCDMS3_xlBSA_293_bsadb_xinet,"20170830_HCDMS3_xlBSA_293_bsadb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
HCDMS31_xlBSA_293_bsadb$PositionA=pmax(HCDMS31_xlBSA_293_bsadb$Position.A,HCDMS31_xlBSA_293_bsadb$Position.B)
HCDMS31_xlBSA_293_bsadb$PositionB=pmin(HCDMS31_xlBSA_293_bsadb$Position.A,HCDMS31_xlBSA_293_bsadb$Position.B)
HCDMS32_xlBSA_293_bsadb$PositionA=pmax(HCDMS32_xlBSA_293_bsadb$Position.A,HCDMS32_xlBSA_293_bsadb$Position.B)
HCDMS32_xlBSA_293_bsadb$PositionB=pmin(HCDMS32_xlBSA_293_bsadb$Position.A,HCDMS32_xlBSA_293_bsadb$Position.B)
HCDMS33_xlBSA_293_bsadb$PositionA=pmax(HCDMS33_xlBSA_293_bsadb$Position.A,HCDMS33_xlBSA_293_bsadb$Position.B)
HCDMS33_xlBSA_293_bsadb$PositionB=pmin(HCDMS33_xlBSA_293_bsadb$Position.A,HCDMS33_xlBSA_293_bsadb$Position.B)
#Filter XlinkX Score >80
HCDMS31_xlBSA_293_bsadb=HCDMS31_xlBSA_293_bsadb[HCDMS31_xlBSA_293_bsadb[,"Max.XlinkX.Score"]>s16,]
HCDMS32_xlBSA_293_bsadb=HCDMS32_xlBSA_293_bsadb[HCDMS32_xlBSA_293_bsadb[,"Max.XlinkX.Score"]>s16,]
HCDMS33_xlBSA_293_bsadb=HCDMS33_xlBSA_293_bsadb[HCDMS33_xlBSA_293_bsadb[,"Max.XlinkX.Score"]>s16,]
#Remove Duplicates
HCDMS31_xlBSA_293_bsadb_dup=HCDMS31_xlBSA_293_bsadb[!duplicated(HCDMS31_xlBSA_293_bsadb[c("PositionA","PositionB")]),]
HCDMS32_xlBSA_293_bsadb_dup=HCDMS32_xlBSA_293_bsadb[!duplicated(HCDMS32_xlBSA_293_bsadb[c("PositionA","PositionB")]),]
HCDMS33_xlBSA_293_bsadb_dup=HCDMS33_xlBSA_293_bsadb[!duplicated(HCDMS33_xlBSA_293_bsadb[c("PositionA","PositionB")]),]
#HCDMS3 Number of BSA crosslinks
HCDMS31_xlBSA_293_bsadb_bsaonly<-HCDMS31_xlBSA_293_bsadb_dup[grepl("P02769", HCDMS31_xlBSA_293_bsadb_dup[["Accession.A"]]) & grepl("P02769", HCDMS31_xlBSA_293_bsadb_dup[["Accession.B"]]), ]
HCDMS32_xlBSA_293_bsadb_bsaonly<-HCDMS32_xlBSA_293_bsadb_dup[grepl("P02769", HCDMS32_xlBSA_293_bsadb_dup[["Accession.A"]]) & grepl("P02769", HCDMS32_xlBSA_293_bsadb_dup[["Accession.B"]]), ]
HCDMS33_xlBSA_293_bsadb_bsaonly<-HCDMS33_xlBSA_293_bsadb_dup[grepl("P02769", HCDMS33_xlBSA_293_bsadb_dup[["Accession.A"]]) & grepl("P02769", HCDMS33_xlBSA_293_bsadb_dup[["Accession.B"]]), ]
#HCDMS3 Number of all crosslinks, no duplicates
HCDMS3_xlBSA_293_bsadb_crosslinks=c(nrow(HCDMS31_xlBSA_293_bsadb_dup),nrow(HCDMS32_xlBSA_293_bsadb_dup),nrow(HCDMS33_xlBSA_293_bsadb_dup))
#HCDMS3 Number of BSA crosslinks, no duplicates
HCDMS3_xlBSA_293_bsadb_bsaonly=c(nrow(HCDMS31_xlBSA_293_bsadb_bsaonly),nrow(HCDMS32_xlBSA_293_bsadb_bsaonly),nrow(HCDMS33_xlBSA_293_bsadb_bsaonly))
#HCDMS3 Number of non-BSA crosslinks
HCDMS3_xlBSA_293_bsadb_nonbsa=(HCDMS3_xlBSA_293_bsadb_crosslinks-HCDMS3_xlBSA_293_bsadb_bsaonly)
#Output for Number of crosslinks
HCDMS3_xlBSA_293_bsadb=data.frame(HCDMS3_xlBSA_293_bsadb_crosslinks,HCDMS3_xlBSA_293_bsadb_bsaonly,HCDMS3_xlBSA_293_bsadb_nonbsa)
colnames(HCDMS3_xlBSA_293_bsadb)=c("HCDMS3_xlBSA_293_bsadb_Total","HCDMS3_xlBSA_293_bsadb_bsacrosslinks","HCDMS3_xlBSA_293_bsadb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(HCDMS31_xlBSA_293_bsadb_bsaonly,HCDMS32_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(HCDMS31_xlBSA_293_bsadb_bsaonly,HCDMS33_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(HCDMS32_xlBSA_293_bsadb_bsaonly,HCDMS33_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB")))
tiff('HCDMS3_xlBSA_293_bsadb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(HCDMS31_xlBSA_293_bsadb_bsaonly),area2=nrow(HCDMS32_xlBSA_293_bsadb_bsaonly),area3=nrow(HCDMS33_xlBSA_293_bsadb_bsaonly),
n12=nrow(merge(HCDMS31_xlBSA_293_bsadb_bsaonly,HCDMS32_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(HCDMS31_xlBSA_293_bsadb_bsaonly,HCDMS33_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(HCDMS32_xlBSA_293_bsadb_bsaonly,HCDMS33_xlBSA_293_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(HCDMS3_xlBSA_293_bsadb_dup)-(nrow(HCDMS31_xlBSA_293_bsadb_bsaonly)+nrow(HCDMS32_xlBSA_293_bsadb_bsaonly)+nrow(HCDMS33_xlBSA_293_bsadb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("HCDMS3_xlBSA_293_bsadb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()

Crosslinks=data.frame(ETD_xlBSA_hdb,HCD_xlBSA_hdb,ETHCD_xlBSA_hdb,HCDMS3_xlBSA_hdb,ETD_xlBSA_293_hdb,HCD_xlBSA_293_hdb,ETHCD_xlBSA_293_hdb,HCDMS3_xlBSA_293_hdb,
ETD_xlBSA_bsadb,HCD_xlBSA_bsadb,ETHCD_xlBSA_bsadb,HCDMS3_xlBSA_bsadb,ETD_xlBSA_293_bsadb,HCD_xlBSA_293_bsadb,ETHCD_xlBSA_293_bsadb,HCDMS3_xlBSA_293_bsadb)
write.table(Crosslinks,"20170830_Crosslinkscompiled.txt",sep="\t",row.names=FALSE)