rm(list=ls())
ls()
getwd()
setwd("Z:\\04_projects\\XLMS_proteomics\\BSA_XLMS\\20171009_xlBSA_proteases_replicates\\xlBSA_proteases_replicates_BsaDB")
library(dplyr)
library(plyr)
library(VennDiagram)
library(gridExtra)

#Score cutoffs
s1=0  #Chymotrypsin BSA db, BSA only
s2=0  #Trypsin BSA db, BSA only
s3=0  #GluC BSA db, BSA only
s4=0  #Chymotrypsin BSA db, BSA+hek293t
s5=0  #Trypsin BSA db, BSA+hek293t
s6=0  #GluC BSA db, BSA+hek293t

##Chymotrypsin cross-linked BSA, BSA database
CHYMO1_xlBSA_bsadb=read.table("170920_xlBSA_chymotrypsin_1_Crosslinks.txt", fill=TRUE, header=TRUE)
CHYMO2_xlBSA_bsadb=read.table("170920_xlBSA_chymotrypsin_2_Crosslinks.txt", fill=TRUE, header=TRUE)
CHYMO3_xlBSA_bsadb=read.table("170920_xlBSA_chymotrypsin_3_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
CHYMO1_xlBSA_bsadb$Replicate=c(1)
CHYMO2_xlBSA_bsadb$Replicate=c(2)
CHYMO3_xlBSA_bsadb$Replicate=c(3)
##Compile all replicates
CHYMO_xlBSA_bsadb=rbind(CHYMO1_xlBSA_bsadb,CHYMO2_xlBSA_bsadb,CHYMO3_xlBSA_bsadb)
write.table(CHYMO_xlBSA_bsadb,"20171010_CHYMO_xlBSA_bsadb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
CHYMO_xlBSA_bsadb$PositionA=pmax(CHYMO_xlBSA_bsadb$Position.A,CHYMO_xlBSA_bsadb$Position.B)
CHYMO_xlBSA_bsadb$PositionB=pmin(CHYMO_xlBSA_bsadb$Position.A,CHYMO_xlBSA_bsadb$Position.B)
##Keep BSA cross-links, duplicates kept
CHYMO_xlBSA_bsadb_bsaonly<- CHYMO_xlBSA_bsadb[grepl("P02769", CHYMO_xlBSA_bsadb[["Accession.A"]]) & grepl("P02769", CHYMO_xlBSA_bsadb[["Accession.B"]]), ]
write.table(CHYMO_xlBSA_bsadb_bsaonly,"20171010_CHYMO_xlBSA_bsadb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
CHYMO_xlBSA_bsadb_nonbsa<- CHYMO_xlBSA_bsadb[!grepl("P02769", CHYMO_xlBSA_bsadb[["Accession.A"]]) | !grepl("P02769", CHYMO_xlBSA_bsadb[["Accession.B"]]), ]
write.table(CHYMO_xlBSA_bsadb_nonbsa,"20171010_CHYMO_xlBSA_bsadb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
CHYMO_xlBSA_bsadb_dup=CHYMO_xlBSA_bsadb_bsaonly[!duplicated(CHYMO_xlBSA_bsadb_bsaonly[c("PositionA","PositionB")]),]
write.table(CHYMO_xlBSA_bsadb_dup,"20171010_CHYMO_xlBSA_bsadb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
CHYMO_xlBSA_bsadb_dup$Accession.A=as.character(CHYMO_xlBSA_bsadb_dup$Accession.A)
CHYMO_xlBSA_bsadb_dup$Accession.A<-replace(CHYMO_xlBSA_bsadb_dup$Accession.A,CHYMO_xlBSA_bsadb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
CHYMO_xlBSA_bsadb_dup$Accession.B=as.character(CHYMO_xlBSA_bsadb_dup$Accession.B)
CHYMO_xlBSA_bsadb_dup$Accession.B<-replace(CHYMO_xlBSA_bsadb_dup$Accession.B,CHYMO_xlBSA_bsadb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
CHYMO_xlBSA_bsadb_xinet<-data.frame(CHYMO_xlBSA_bsadb_dup$Max.XlinkX.Score,as.factor(CHYMO_xlBSA_bsadb_dup$Accession.A),as.factor(CHYMO_xlBSA_bsadb_dup$Accession.B),
CHYMO_xlBSA_bsadb_dup$PositionA,CHYMO_xlBSA_bsadb_dup$PositionB)
colnames(CHYMO_xlBSA_bsadb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(CHYMO_xlBSA_bsadb_xinet,"20171010_CHYMO_xlBSA_bsadb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
CHYMO1_xlBSA_bsadb$PositionA=pmax(CHYMO1_xlBSA_bsadb$Position.A,CHYMO1_xlBSA_bsadb$Position.B)
CHYMO1_xlBSA_bsadb$PositionB=pmin(CHYMO1_xlBSA_bsadb$Position.A,CHYMO1_xlBSA_bsadb$Position.B)
CHYMO2_xlBSA_bsadb$PositionA=pmax(CHYMO2_xlBSA_bsadb$Position.A,CHYMO2_xlBSA_bsadb$Position.B)
CHYMO2_xlBSA_bsadb$PositionB=pmin(CHYMO2_xlBSA_bsadb$Position.A,CHYMO2_xlBSA_bsadb$Position.B)
CHYMO3_xlBSA_bsadb$PositionA=pmax(CHYMO3_xlBSA_bsadb$Position.A,CHYMO3_xlBSA_bsadb$Position.B)
CHYMO3_xlBSA_bsadb$PositionB=pmin(CHYMO3_xlBSA_bsadb$Position.A,CHYMO3_xlBSA_bsadb$Position.B)
#Filter XlinkX Score >80
CHYMO1_xlBSA_bsadb=CHYMO1_xlBSA_bsadb[CHYMO1_xlBSA_bsadb[,"Max.XlinkX.Score"]>s1,]
CHYMO2_xlBSA_bsadb=CHYMO2_xlBSA_bsadb[CHYMO2_xlBSA_bsadb[,"Max.XlinkX.Score"]>s1,]
CHYMO3_xlBSA_bsadb=CHYMO3_xlBSA_bsadb[CHYMO3_xlBSA_bsadb[,"Max.XlinkX.Score"]>s1,]
#Remove Duplicates
CHYMO1_xlBSA_bsadb_dup=CHYMO1_xlBSA_bsadb[!duplicated(CHYMO1_xlBSA_bsadb[c("PositionA","PositionB")]),]
CHYMO2_xlBSA_bsadb_dup=CHYMO2_xlBSA_bsadb[!duplicated(CHYMO2_xlBSA_bsadb[c("PositionA","PositionB")]),]
CHYMO3_xlBSA_bsadb_dup=CHYMO3_xlBSA_bsadb[!duplicated(CHYMO3_xlBSA_bsadb[c("PositionA","PositionB")]),]
#CHYMO Number of BSA crosslinks
CHYMO1_xlBSA_bsadb_bsaonly<- CHYMO1_xlBSA_bsadb_dup[grepl("P02769", CHYMO1_xlBSA_bsadb_dup[["Accession.A"]]) & grepl("P02769", CHYMO1_xlBSA_bsadb_dup[["Accession.B"]]), ]
CHYMO2_xlBSA_bsadb_bsaonly<- CHYMO2_xlBSA_bsadb_dup[grepl("P02769", CHYMO2_xlBSA_bsadb_dup[["Accession.A"]]) & grepl("P02769", CHYMO2_xlBSA_bsadb_dup[["Accession.B"]]), ]
CHYMO3_xlBSA_bsadb_bsaonly<- CHYMO3_xlBSA_bsadb_dup[grepl("P02769", CHYMO3_xlBSA_bsadb_dup[["Accession.A"]]) & grepl("P02769", CHYMO3_xlBSA_bsadb_dup[["Accession.B"]]), ]
#CHYMO Number of all crosslinks, no duplicates
CHYMO_xlBSA_bsadb_crosslinks=c(nrow(CHYMO1_xlBSA_bsadb_dup),nrow(CHYMO2_xlBSA_bsadb_dup),nrow(CHYMO3_xlBSA_bsadb_dup))
#CHYMO Number of BSA crosslinks, no duplicates
CHYMO_xlBSA_bsadb_bsaonly=c(nrow(CHYMO1_xlBSA_bsadb_bsaonly),nrow(CHYMO2_xlBSA_bsadb_bsaonly),nrow(CHYMO3_xlBSA_bsadb_bsaonly))
#CHYMO Number of non-BSA crosslinks
CHYMO_xlBSA_bsadb_nonbsa=(CHYMO_xlBSA_bsadb_crosslinks-CHYMO_xlBSA_bsadb_bsaonly)
#Output for Number of crosslinks
CHYMO_xlBSA_bsadb=data.frame(CHYMO_xlBSA_bsadb_crosslinks,CHYMO_xlBSA_bsadb_bsaonly,CHYMO_xlBSA_bsadb_nonbsa)
colnames(CHYMO_xlBSA_bsadb)=c("CHYMO_xlBSA_bsadb_Total","CHYMO_xlBSA_bsadb_bsacrosslinks","CHYMO_xlBSA_bsadb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(CHYMO1_xlBSA_bsadb_bsaonly,CHYMO2_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(CHYMO1_xlBSA_bsadb_bsaonly,CHYMO3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(CHYMO2_xlBSA_bsadb_bsaonly,CHYMO3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB")))
tiff('CHYMO_xlBSA_bsadb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(CHYMO1_xlBSA_bsadb_bsaonly),area2=nrow(CHYMO2_xlBSA_bsadb_bsaonly),area3=nrow(CHYMO3_xlBSA_bsadb_bsaonly),
n12=nrow(merge(CHYMO1_xlBSA_bsadb_bsaonly,CHYMO2_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(CHYMO1_xlBSA_bsadb_bsaonly,CHYMO3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(CHYMO2_xlBSA_bsadb_bsaonly,CHYMO3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(CHYMO_xlBSA_bsadb_dup)-(nrow(CHYMO1_xlBSA_bsadb_bsaonly)+nrow(CHYMO2_xlBSA_bsadb_bsaonly)+nrow(CHYMO3_xlBSA_bsadb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("CHYMO_xlBSA_bsadb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()

##Trypsin cross-linked BSA, BSA database
TRYP1_xlBSA_bsadb=read.table("170920_xlBSA_trypsin_1_Crosslinks.txt", fill=TRUE, header=TRUE)
TRYP2_xlBSA_bsadb=read.table("170920_xlBSA_trypsin_2_Crosslinks.txt", fill=TRUE, header=TRUE)
TRYP3_xlBSA_bsadb=read.table("170920_xlBSA_trypsin_3_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
TRYP1_xlBSA_bsadb$Replicate=c(1)
TRYP2_xlBSA_bsadb$Replicate=c(2)
TRYP3_xlBSA_bsadb$Replicate=c(3)
##Compile all replicates
TRYP_xlBSA_bsadb=rbind(TRYP1_xlBSA_bsadb,TRYP2_xlBSA_bsadb,TRYP3_xlBSA_bsadb)
write.table(TRYP_xlBSA_bsadb,"20171010_TRYP_xlBSA_bsadb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
TRYP_xlBSA_bsadb$PositionA=pmax(TRYP_xlBSA_bsadb$Position.A,TRYP_xlBSA_bsadb$Position.B)
TRYP_xlBSA_bsadb$PositionB=pmin(TRYP_xlBSA_bsadb$Position.A,TRYP_xlBSA_bsadb$Position.B)
##Keep BSA cross-links, duplicates kept
TRYP_xlBSA_bsadb_bsaonly<- TRYP_xlBSA_bsadb[grepl("P02769", TRYP_xlBSA_bsadb[["Accession.A"]]) & grepl("P02769", TRYP_xlBSA_bsadb[["Accession.B"]]), ]
write.table(TRYP_xlBSA_bsadb_bsaonly,"20171010_TRYP_xlBSA_bsadb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
TRYP_xlBSA_bsadb_nonbsa<- TRYP_xlBSA_bsadb[!grepl("P02769", TRYP_xlBSA_bsadb[["Accession.A"]]) | !grepl("P02769", TRYP_xlBSA_bsadb[["Accession.B"]]), ]
write.table(TRYP_xlBSA_bsadb_nonbsa,"20171010_TRYP_xlBSA_bsadb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
TRYP_xlBSA_bsadb_dup=TRYP_xlBSA_bsadb_bsaonly[!duplicated(TRYP_xlBSA_bsadb_bsaonly[c("PositionA","PositionB")]),]
write.table(TRYP_xlBSA_bsadb_dup,"20171010_TRYP_xlBSA_bsadb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
TRYP_xlBSA_bsadb_dup$Accession.A=as.character(TRYP_xlBSA_bsadb_dup$Accession.A)
TRYP_xlBSA_bsadb_dup$Accession.A<-replace(TRYP_xlBSA_bsadb_dup$Accession.A,TRYP_xlBSA_bsadb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
TRYP_xlBSA_bsadb_dup$Accession.B=as.character(TRYP_xlBSA_bsadb_dup$Accession.B)
TRYP_xlBSA_bsadb_dup$Accession.B<-replace(TRYP_xlBSA_bsadb_dup$Accession.B,TRYP_xlBSA_bsadb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
TRYP_xlBSA_bsadb_xinet<-data.frame(TRYP_xlBSA_bsadb_dup$Max.XlinkX.Score,as.factor(TRYP_xlBSA_bsadb_dup$Accession.A),as.factor(TRYP_xlBSA_bsadb_dup$Accession.B),
TRYP_xlBSA_bsadb_dup$PositionA,TRYP_xlBSA_bsadb_dup$PositionB)
colnames(TRYP_xlBSA_bsadb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(TRYP_xlBSA_bsadb_xinet,"20171010_TRYP_xlBSA_bsadb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
TRYP1_xlBSA_bsadb$PositionA=pmax(TRYP1_xlBSA_bsadb$Position.A,TRYP1_xlBSA_bsadb$Position.B)
TRYP1_xlBSA_bsadb$PositionB=pmin(TRYP1_xlBSA_bsadb$Position.A,TRYP1_xlBSA_bsadb$Position.B)
TRYP2_xlBSA_bsadb$PositionA=pmax(TRYP2_xlBSA_bsadb$Position.A,TRYP2_xlBSA_bsadb$Position.B)
TRYP2_xlBSA_bsadb$PositionB=pmin(TRYP2_xlBSA_bsadb$Position.A,TRYP2_xlBSA_bsadb$Position.B)
TRYP3_xlBSA_bsadb$PositionA=pmax(TRYP3_xlBSA_bsadb$Position.A,TRYP3_xlBSA_bsadb$Position.B)
TRYP3_xlBSA_bsadb$PositionB=pmin(TRYP3_xlBSA_bsadb$Position.A,TRYP3_xlBSA_bsadb$Position.B)
#Filter XlinkX Score >80
TRYP1_xlBSA_bsadb=TRYP1_xlBSA_bsadb[TRYP1_xlBSA_bsadb[,"Max.XlinkX.Score"]>s2,]
TRYP2_xlBSA_bsadb=TRYP2_xlBSA_bsadb[TRYP2_xlBSA_bsadb[,"Max.XlinkX.Score"]>s2,]
TRYP3_xlBSA_bsadb=TRYP3_xlBSA_bsadb[TRYP3_xlBSA_bsadb[,"Max.XlinkX.Score"]>s2,]
#Remove Duplicates
TRYP1_xlBSA_bsadb_dup=TRYP1_xlBSA_bsadb[!duplicated(TRYP1_xlBSA_bsadb[c("PositionA","PositionB")]),]
TRYP2_xlBSA_bsadb_dup=TRYP2_xlBSA_bsadb[!duplicated(TRYP2_xlBSA_bsadb[c("PositionA","PositionB")]),]
TRYP3_xlBSA_bsadb_dup=TRYP3_xlBSA_bsadb[!duplicated(TRYP3_xlBSA_bsadb[c("PositionA","PositionB")]),]
#TRYP Number of BSA crosslinks
TRYP1_xlBSA_bsadb_bsaonly<- TRYP1_xlBSA_bsadb_dup[grepl("P02769", TRYP1_xlBSA_bsadb_dup[["Accession.A"]]) & grepl("P02769", TRYP1_xlBSA_bsadb_dup[["Accession.B"]]), ]
TRYP2_xlBSA_bsadb_bsaonly<- TRYP2_xlBSA_bsadb_dup[grepl("P02769", TRYP2_xlBSA_bsadb_dup[["Accession.A"]]) & grepl("P02769", TRYP2_xlBSA_bsadb_dup[["Accession.B"]]), ]
TRYP3_xlBSA_bsadb_bsaonly<- TRYP3_xlBSA_bsadb_dup[grepl("P02769", TRYP3_xlBSA_bsadb_dup[["Accession.A"]]) & grepl("P02769", TRYP3_xlBSA_bsadb_dup[["Accession.B"]]), ]
#TRYP Number of all crosslinks, no duplicates
TRYP_xlBSA_bsadb_crosslinks=c(nrow(TRYP1_xlBSA_bsadb_dup),nrow(TRYP2_xlBSA_bsadb_dup),nrow(TRYP3_xlBSA_bsadb_dup))
#TRYP Number of BSA crosslinks, no duplicates
TRYP_xlBSA_bsadb_bsaonly=c(nrow(TRYP1_xlBSA_bsadb_bsaonly),nrow(TRYP2_xlBSA_bsadb_bsaonly),nrow(TRYP3_xlBSA_bsadb_bsaonly))
#TRYP Number of non-BSA crosslinks
TRYP_xlBSA_bsadb_nonbsa=(TRYP_xlBSA_bsadb_crosslinks-TRYP_xlBSA_bsadb_bsaonly)
#Output for Number of crosslinks
TRYP_xlBSA_bsadb=data.frame(TRYP_xlBSA_bsadb_crosslinks,TRYP_xlBSA_bsadb_bsaonly,TRYP_xlBSA_bsadb_nonbsa)
colnames(TRYP_xlBSA_bsadb)=c("TRYP_xlBSA_bsadb_Total","TRYP_xlBSA_bsadb_bsacrosslinks","TRYP_xlBSA_bsadb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(TRYP1_xlBSA_bsadb_bsaonly,TRYP2_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(TRYP1_xlBSA_bsadb_bsaonly,TRYP3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(TRYP2_xlBSA_bsadb_bsaonly,TRYP3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB")))
tiff('TRYP_xlBSA_bsadb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(TRYP1_xlBSA_bsadb_bsaonly),area2=nrow(TRYP2_xlBSA_bsadb_bsaonly),area3=nrow(TRYP3_xlBSA_bsadb_bsaonly),
n12=nrow(merge(TRYP1_xlBSA_bsadb_bsaonly,TRYP2_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(TRYP1_xlBSA_bsadb_bsaonly,TRYP3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(TRYP2_xlBSA_bsadb_bsaonly,TRYP3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(TRYP_xlBSA_bsadb_dup)-(nrow(TRYP1_xlBSA_bsadb_bsaonly)+nrow(TRYP2_xlBSA_bsadb_bsaonly)+nrow(TRYP3_xlBSA_bsadb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("TRYP_xlBSA_bsadb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()

##GluC cross-linked BSA, BSA database
GLUC1_xlBSA_bsadb=read.table("170920_xlBSA_gluC_1_Crosslinks.txt", fill=TRUE, header=TRUE)
GLUC2_xlBSA_bsadb=read.table("170920_xlBSA_gluC_2_Crosslinks.txt", fill=TRUE, header=TRUE)
GLUC3_xlBSA_bsadb=read.table("170920_xlBSA_gluC_3_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
GLUC1_xlBSA_bsadb$Replicate=c(1)
GLUC2_xlBSA_bsadb$Replicate=c(2)
GLUC3_xlBSA_bsadb$Replicate=c(3)
##Compile all replicates
GLUC_xlBSA_bsadb=rbind(GLUC1_xlBSA_bsadb,GLUC2_xlBSA_bsadb,GLUC3_xlBSA_bsadb)
write.table(GLUC_xlBSA_bsadb,"20171010_GLUC_xlBSA_bsadb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
GLUC_xlBSA_bsadb$PositionA=pmax(GLUC_xlBSA_bsadb$Position.A,GLUC_xlBSA_bsadb$Position.B)
GLUC_xlBSA_bsadb$PositionB=pmin(GLUC_xlBSA_bsadb$Position.A,GLUC_xlBSA_bsadb$Position.B)
##Keep BSA cross-links, duplicates kept
GLUC_xlBSA_bsadb_bsaonly<- GLUC_xlBSA_bsadb[grepl("P02769", GLUC_xlBSA_bsadb[["Accession.A"]]) & grepl("P02769", GLUC_xlBSA_bsadb[["Accession.B"]]), ]
write.table(GLUC_xlBSA_bsadb_bsaonly,"20171010_GLUC_xlBSA_bsadb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
GLUC_xlBSA_bsadb_nonbsa<- GLUC_xlBSA_bsadb[!grepl("P02769", GLUC_xlBSA_bsadb[["Accession.A"]]) | !grepl("P02769", GLUC_xlBSA_bsadb[["Accession.B"]]), ]
write.table(GLUC_xlBSA_bsadb_nonbsa,"20171010_GLUC_xlBSA_bsadb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
GLUC_xlBSA_bsadb_dup=GLUC_xlBSA_bsadb_bsaonly[!duplicated(GLUC_xlBSA_bsadb_bsaonly[c("PositionA","PositionB")]),]
write.table(GLUC_xlBSA_bsadb_dup,"20171010_GLUC_xlBSA_bsadb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
GLUC_xlBSA_bsadb_dup$Accession.A=as.character(GLUC_xlBSA_bsadb_dup$Accession.A)
GLUC_xlBSA_bsadb_dup$Accession.A<-replace(GLUC_xlBSA_bsadb_dup$Accession.A,GLUC_xlBSA_bsadb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
GLUC_xlBSA_bsadb_dup$Accession.B=as.character(GLUC_xlBSA_bsadb_dup$Accession.B)
GLUC_xlBSA_bsadb_dup$Accession.B<-replace(GLUC_xlBSA_bsadb_dup$Accession.B,GLUC_xlBSA_bsadb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
GLUC_xlBSA_bsadb_xinet<-data.frame(GLUC_xlBSA_bsadb_dup$Max.XlinkX.Score,as.factor(GLUC_xlBSA_bsadb_dup$Accession.A),as.factor(GLUC_xlBSA_bsadb_dup$Accession.B),
GLUC_xlBSA_bsadb_dup$PositionA,GLUC_xlBSA_bsadb_dup$PositionB)
colnames(GLUC_xlBSA_bsadb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(GLUC_xlBSA_bsadb_xinet,"20171010_GLUC_xlBSA_bsadb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
GLUC1_xlBSA_bsadb$PositionA=pmax(GLUC1_xlBSA_bsadb$Position.A,GLUC1_xlBSA_bsadb$Position.B)
GLUC1_xlBSA_bsadb$PositionB=pmin(GLUC1_xlBSA_bsadb$Position.A,GLUC1_xlBSA_bsadb$Position.B)
GLUC2_xlBSA_bsadb$PositionA=pmax(GLUC2_xlBSA_bsadb$Position.A,GLUC2_xlBSA_bsadb$Position.B)
GLUC2_xlBSA_bsadb$PositionB=pmin(GLUC2_xlBSA_bsadb$Position.A,GLUC2_xlBSA_bsadb$Position.B)
GLUC3_xlBSA_bsadb$PositionA=pmax(GLUC3_xlBSA_bsadb$Position.A,GLUC3_xlBSA_bsadb$Position.B)
GLUC3_xlBSA_bsadb$PositionB=pmin(GLUC3_xlBSA_bsadb$Position.A,GLUC3_xlBSA_bsadb$Position.B)
#Filter XlinkX Score >80
GLUC1_xlBSA_bsadb=GLUC1_xlBSA_bsadb[GLUC1_xlBSA_bsadb[,"Max.XlinkX.Score"]>s3,]
GLUC2_xlBSA_bsadb=GLUC2_xlBSA_bsadb[GLUC2_xlBSA_bsadb[,"Max.XlinkX.Score"]>s3,]
GLUC3_xlBSA_bsadb=GLUC3_xlBSA_bsadb[GLUC3_xlBSA_bsadb[,"Max.XlinkX.Score"]>s3,]
#Remove Duplicates
GLUC1_xlBSA_bsadb_dup=GLUC1_xlBSA_bsadb[!duplicated(GLUC1_xlBSA_bsadb[c("PositionA","PositionB")]),]
GLUC2_xlBSA_bsadb_dup=GLUC2_xlBSA_bsadb[!duplicated(GLUC2_xlBSA_bsadb[c("PositionA","PositionB")]),]
GLUC3_xlBSA_bsadb_dup=GLUC3_xlBSA_bsadb[!duplicated(GLUC3_xlBSA_bsadb[c("PositionA","PositionB")]),]
#GLUC Number of BSA crosslinks
GLUC1_xlBSA_bsadb_bsaonly<- GLUC1_xlBSA_bsadb_dup[grepl("P02769", GLUC1_xlBSA_bsadb_dup[["Accession.A"]]) & grepl("P02769", GLUC1_xlBSA_bsadb_dup[["Accession.B"]]), ]
GLUC2_xlBSA_bsadb_bsaonly<- GLUC2_xlBSA_bsadb_dup[grepl("P02769", GLUC2_xlBSA_bsadb_dup[["Accession.A"]]) & grepl("P02769", GLUC2_xlBSA_bsadb_dup[["Accession.B"]]), ]
GLUC3_xlBSA_bsadb_bsaonly<- GLUC3_xlBSA_bsadb_dup[grepl("P02769", GLUC3_xlBSA_bsadb_dup[["Accession.A"]]) & grepl("P02769", GLUC3_xlBSA_bsadb_dup[["Accession.B"]]), ]
#GLUC Number of all crosslinks, no duplicates
GLUC_xlBSA_bsadb_crosslinks=c(nrow(GLUC1_xlBSA_bsadb_dup),nrow(GLUC2_xlBSA_bsadb_dup),nrow(GLUC3_xlBSA_bsadb_dup))
#GLUC Number of BSA crosslinks, no duplicates
GLUC_xlBSA_bsadb_bsaonly=c(nrow(GLUC1_xlBSA_bsadb_bsaonly),nrow(GLUC2_xlBSA_bsadb_bsaonly),nrow(GLUC3_xlBSA_bsadb_bsaonly))
#GLUC Number of non-BSA crosslinks
GLUC_xlBSA_bsadb_nonbsa=(GLUC_xlBSA_bsadb_crosslinks-GLUC_xlBSA_bsadb_bsaonly)
#Output for Number of crosslinks
GLUC_xlBSA_bsadb=data.frame(GLUC_xlBSA_bsadb_crosslinks,GLUC_xlBSA_bsadb_bsaonly,GLUC_xlBSA_bsadb_nonbsa)
colnames(GLUC_xlBSA_bsadb)=c("GLUC_xlBSA_bsadb_Total","GLUC_xlBSA_bsadb_bsacrosslinks","GLUC_xlBSA_bsadb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(GLUC1_xlBSA_bsadb_bsaonly,GLUC2_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(GLUC1_xlBSA_bsadb_bsaonly,GLUC3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(GLUC2_xlBSA_bsadb_bsaonly,GLUC3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB")))
tiff('GLUC_xlBSA_bsadb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(GLUC1_xlBSA_bsadb_bsaonly),area2=nrow(GLUC2_xlBSA_bsadb_bsaonly),area3=nrow(GLUC3_xlBSA_bsadb_bsaonly),
n12=nrow(merge(GLUC1_xlBSA_bsadb_bsaonly,GLUC2_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(GLUC1_xlBSA_bsadb_bsaonly,GLUC3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(GLUC2_xlBSA_bsadb_bsaonly,GLUC3_xlBSA_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(GLUC_xlBSA_bsadb_dup)-(nrow(GLUC1_xlBSA_bsadb_bsaonly)+nrow(GLUC2_xlBSA_bsadb_bsaonly)+nrow(GLUC3_xlBSA_bsadb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("GLUC_xlBSA_bsadb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()

##Chymotrypsin cross-linked BSA+hek293T, BSA database
CHYMO1_xlBSA_hek293t_bsadb=read.table("170920_xlBSA_hek293t_chymotrypsin_1_Crosslinks.txt", fill=TRUE, header=TRUE)
CHYMO2_xlBSA_hek293t_bsadb=read.table("170920_xlBSA_hek293t_chymotrypsin_2_Crosslinks.txt", fill=TRUE, header=TRUE)
CHYMO3_xlBSA_hek293t_bsadb=read.table("170920_xlBSA_hek293t_chymotrypsin_3_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
CHYMO1_xlBSA_hek293t_bsadb$Replicate=c(1)
CHYMO2_xlBSA_hek293t_bsadb$Replicate=c(2)
CHYMO3_xlBSA_hek293t_bsadb$Replicate=c(3)
##Compile all replicates
CHYMO_xlBSA_hek293t_bsadb=rbind(CHYMO1_xlBSA_hek293t_bsadb,CHYMO2_xlBSA_hek293t_bsadb,CHYMO3_xlBSA_hek293t_bsadb)
write.table(CHYMO_xlBSA_hek293t_bsadb,"20171010_CHYMO_xlBSA_hek293t_bsadb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
CHYMO_xlBSA_hek293t_bsadb$PositionA=pmax(CHYMO_xlBSA_hek293t_bsadb$Position.A,CHYMO_xlBSA_hek293t_bsadb$Position.B)
CHYMO_xlBSA_hek293t_bsadb$PositionB=pmin(CHYMO_xlBSA_hek293t_bsadb$Position.A,CHYMO_xlBSA_hek293t_bsadb$Position.B)
##Keep BSA cross-links, duplicates kept
CHYMO_xlBSA_hek293t_bsadb_bsaonly<- CHYMO_xlBSA_hek293t_bsadb[grepl("P02769", CHYMO_xlBSA_hek293t_bsadb[["Accession.A"]]) & grepl("P02769", CHYMO_xlBSA_hek293t_bsadb[["Accession.B"]]), ]
write.table(CHYMO_xlBSA_hek293t_bsadb_bsaonly,"20171010_CHYMO_xlBSA_hek293t_bsadb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
CHYMO_xlBSA_hek293t_bsadb_nonbsa<- CHYMO_xlBSA_hek293t_bsadb[!grepl("P02769", CHYMO_xlBSA_hek293t_bsadb[["Accession.A"]]) | !grepl("P02769", CHYMO_xlBSA_hek293t_bsadb[["Accession.B"]]), ]
write.table(CHYMO_xlBSA_hek293t_bsadb_nonbsa,"20171010_CHYMO_xlBSA_hek293t_bsadb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
CHYMO_xlBSA_hek293t_bsadb_dup=CHYMO_xlBSA_hek293t_bsadb_bsaonly[!duplicated(CHYMO_xlBSA_hek293t_bsadb_bsaonly[c("PositionA","PositionB")]),]
write.table(CHYMO_xlBSA_hek293t_bsadb_dup,"20171010_CHYMO_xlBSA_hek293t_bsadb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
CHYMO_xlBSA_hek293t_bsadb_dup$Accession.A=as.character(CHYMO_xlBSA_hek293t_bsadb_dup$Accession.A)
CHYMO_xlBSA_hek293t_bsadb_dup$Accession.A<-replace(CHYMO_xlBSA_hek293t_bsadb_dup$Accession.A,CHYMO_xlBSA_hek293t_bsadb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
CHYMO_xlBSA_hek293t_bsadb_dup$Accession.B=as.character(CHYMO_xlBSA_hek293t_bsadb_dup$Accession.B)
CHYMO_xlBSA_hek293t_bsadb_dup$Accession.B<-replace(CHYMO_xlBSA_hek293t_bsadb_dup$Accession.B,CHYMO_xlBSA_hek293t_bsadb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
CHYMO_xlBSA_hek293t_bsadb_xinet<-data.frame(CHYMO_xlBSA_hek293t_bsadb_dup$Max.XlinkX.Score,as.factor(CHYMO_xlBSA_hek293t_bsadb_dup$Accession.A),as.factor(CHYMO_xlBSA_hek293t_bsadb_dup$Accession.B),
CHYMO_xlBSA_hek293t_bsadb_dup$PositionA,CHYMO_xlBSA_hek293t_bsadb_dup$PositionB)
colnames(CHYMO_xlBSA_hek293t_bsadb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(CHYMO_xlBSA_hek293t_bsadb_xinet,"20171010_CHYMO_xlBSA_hek293t_bsadb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
CHYMO1_xlBSA_hek293t_bsadb$PositionA=pmax(CHYMO1_xlBSA_hek293t_bsadb$Position.A,CHYMO1_xlBSA_hek293t_bsadb$Position.B)
CHYMO1_xlBSA_hek293t_bsadb$PositionB=pmin(CHYMO1_xlBSA_hek293t_bsadb$Position.A,CHYMO1_xlBSA_hek293t_bsadb$Position.B)
CHYMO2_xlBSA_hek293t_bsadb$PositionA=pmax(CHYMO2_xlBSA_hek293t_bsadb$Position.A,CHYMO2_xlBSA_hek293t_bsadb$Position.B)
CHYMO2_xlBSA_hek293t_bsadb$PositionB=pmin(CHYMO2_xlBSA_hek293t_bsadb$Position.A,CHYMO2_xlBSA_hek293t_bsadb$Position.B)
CHYMO3_xlBSA_hek293t_bsadb$PositionA=pmax(CHYMO3_xlBSA_hek293t_bsadb$Position.A,CHYMO3_xlBSA_hek293t_bsadb$Position.B)
CHYMO3_xlBSA_hek293t_bsadb$PositionB=pmin(CHYMO3_xlBSA_hek293t_bsadb$Position.A,CHYMO3_xlBSA_hek293t_bsadb$Position.B)
#Filter XlinkX Score >80
CHYMO1_xlBSA_hek293t_bsadb=CHYMO1_xlBSA_hek293t_bsadb[CHYMO1_xlBSA_hek293t_bsadb[,"Max.XlinkX.Score"]>s4,]
CHYMO2_xlBSA_hek293t_bsadb=CHYMO2_xlBSA_hek293t_bsadb[CHYMO2_xlBSA_hek293t_bsadb[,"Max.XlinkX.Score"]>s4,]
CHYMO3_xlBSA_hek293t_bsadb=CHYMO3_xlBSA_hek293t_bsadb[CHYMO3_xlBSA_hek293t_bsadb[,"Max.XlinkX.Score"]>s4,]
#Remove Duplicates
CHYMO1_xlBSA_hek293t_bsadb_dup=CHYMO1_xlBSA_hek293t_bsadb[!duplicated(CHYMO1_xlBSA_hek293t_bsadb[c("PositionA","PositionB")]),]
CHYMO2_xlBSA_hek293t_bsadb_dup=CHYMO2_xlBSA_hek293t_bsadb[!duplicated(CHYMO2_xlBSA_hek293t_bsadb[c("PositionA","PositionB")]),]
CHYMO3_xlBSA_hek293t_bsadb_dup=CHYMO3_xlBSA_hek293t_bsadb[!duplicated(CHYMO3_xlBSA_hek293t_bsadb[c("PositionA","PositionB")]),]
#CHYMO Number of BSA crosslinks
CHYMO1_xlBSA_hek293t_bsadb_bsaonly<- CHYMO1_xlBSA_hek293t_bsadb_dup[grepl("P02769", CHYMO1_xlBSA_hek293t_bsadb_dup[["Accession.A"]]) & grepl("P02769", CHYMO1_xlBSA_hek293t_bsadb_dup[["Accession.B"]]), ]
CHYMO2_xlBSA_hek293t_bsadb_bsaonly<- CHYMO2_xlBSA_hek293t_bsadb_dup[grepl("P02769", CHYMO2_xlBSA_hek293t_bsadb_dup[["Accession.A"]]) & grepl("P02769", CHYMO2_xlBSA_hek293t_bsadb_dup[["Accession.B"]]), ]
CHYMO3_xlBSA_hek293t_bsadb_bsaonly<- CHYMO3_xlBSA_hek293t_bsadb_dup[grepl("P02769", CHYMO3_xlBSA_hek293t_bsadb_dup[["Accession.A"]]) & grepl("P02769", CHYMO3_xlBSA_hek293t_bsadb_dup[["Accession.B"]]), ]
#CHYMO Number of all crosslinks, no duplicates
CHYMO_xlBSA_hek293t_bsadb_crosslinks=c(nrow(CHYMO1_xlBSA_hek293t_bsadb_dup),nrow(CHYMO2_xlBSA_hek293t_bsadb_dup),nrow(CHYMO3_xlBSA_hek293t_bsadb_dup))
#CHYMO Number of BSA crosslinks, no duplicates
CHYMO_xlBSA_hek293t_bsadb_bsaonly=c(nrow(CHYMO1_xlBSA_hek293t_bsadb_bsaonly),nrow(CHYMO2_xlBSA_hek293t_bsadb_bsaonly),nrow(CHYMO3_xlBSA_hek293t_bsadb_bsaonly))
#CHYMO Number of non-BSA crosslinks
CHYMO_xlBSA_hek293t_bsadb_nonbsa=(CHYMO_xlBSA_hek293t_bsadb_crosslinks-CHYMO_xlBSA_hek293t_bsadb_bsaonly)
#Output for Number of crosslinks
CHYMO_xlBSA_hek293t_bsadb=data.frame(CHYMO_xlBSA_hek293t_bsadb_crosslinks,CHYMO_xlBSA_hek293t_bsadb_bsaonly,CHYMO_xlBSA_hek293t_bsadb_nonbsa)
colnames(CHYMO_xlBSA_hek293t_bsadb)=c("CHYMO_xlBSA_hek293t_bsadb_Total","CHYMO_xlBSA_hek293t_bsadb_bsacrosslinks","CHYMO_xlBSA_hek293t_bsadb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(CHYMO1_xlBSA_hek293t_bsadb_bsaonly,CHYMO2_xlBSA_hek293t_bsadb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(CHYMO1_xlBSA_hek293t_bsadb_bsaonly,CHYMO3_xlBSA_hek293t_bsadb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(CHYMO2_xlBSA_hek293t_bsadb_bsaonly,CHYMO3_xlBSA_hek293t_bsadb_bsaonly,by=c("PositionA","PositionB")))
tiff('CHYMO_xlBSA_hek293t_bsadb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(CHYMO1_xlBSA_hek293t_bsadb_bsaonly),area2=nrow(CHYMO2_xlBSA_hek293t_bsadb_bsaonly),area3=nrow(CHYMO3_xlBSA_hek293t_bsadb_bsaonly),
n12=nrow(merge(CHYMO1_xlBSA_hek293t_bsadb_bsaonly,CHYMO2_xlBSA_hek293t_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(CHYMO1_xlBSA_hek293t_bsadb_bsaonly,CHYMO3_xlBSA_hek293t_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(CHYMO2_xlBSA_hek293t_bsadb_bsaonly,CHYMO3_xlBSA_hek293t_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(CHYMO_xlBSA_hek293t_bsadb_dup)-(nrow(CHYMO1_xlBSA_hek293t_bsadb_bsaonly)+nrow(CHYMO2_xlBSA_hek293t_bsadb_bsaonly)+nrow(CHYMO3_xlBSA_hek293t_bsadb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("CHYMO_xlBSA_hek293t_bsadb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()

##Trypsin cross-linked BSA+hek293T, BSA database
TRYP1_xlBSA_hek293t_bsadb=read.table("170920_xlBSA_hek293t_trypsin_1_Crosslinks.txt", fill=TRUE, header=TRUE)
TRYP2_xlBSA_hek293t_bsadb=read.table("170920_xlBSA_hek293t_trypsin_2_Crosslinks.txt", fill=TRUE, header=TRUE)
TRYP3_xlBSA_hek293t_bsadb=read.table("170920_xlBSA_hek293t_trypsin_3_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
TRYP1_xlBSA_hek293t_bsadb$Replicate=c(1)
TRYP2_xlBSA_hek293t_bsadb$Replicate=c(2)
TRYP3_xlBSA_hek293t_bsadb$Replicate=c(3)
##Compile all replicates
TRYP_xlBSA_hek293t_bsadb=rbind(TRYP1_xlBSA_hek293t_bsadb,TRYP2_xlBSA_hek293t_bsadb,TRYP3_xlBSA_hek293t_bsadb)
write.table(TRYP_xlBSA_hek293t_bsadb,"20171010_TRYP_xlBSA_hek293t_bsadb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
TRYP_xlBSA_hek293t_bsadb$PositionA=pmax(TRYP_xlBSA_hek293t_bsadb$Position.A,TRYP_xlBSA_hek293t_bsadb$Position.B)
TRYP_xlBSA_hek293t_bsadb$PositionB=pmin(TRYP_xlBSA_hek293t_bsadb$Position.A,TRYP_xlBSA_hek293t_bsadb$Position.B)
##Keep BSA cross-links, duplicates kept
TRYP_xlBSA_hek293t_bsadb_bsaonly<- TRYP_xlBSA_hek293t_bsadb[grepl("P02769", TRYP_xlBSA_hek293t_bsadb[["Accession.A"]]) & grepl("P02769", TRYP_xlBSA_hek293t_bsadb[["Accession.B"]]), ]
write.table(TRYP_xlBSA_hek293t_bsadb_bsaonly,"20171010_TRYP_xlBSA_hek293t_bsadb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
TRYP_xlBSA_hek293t_bsadb_nonbsa<- TRYP_xlBSA_hek293t_bsadb[!grepl("P02769", TRYP_xlBSA_hek293t_bsadb[["Accession.A"]]) | !grepl("P02769", TRYP_xlBSA_hek293t_bsadb[["Accession.B"]]), ]
write.table(TRYP_xlBSA_hek293t_bsadb_nonbsa,"20171010_TRYP_xlBSA_hek293t_bsadb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
TRYP_xlBSA_hek293t_bsadb_dup=TRYP_xlBSA_hek293t_bsadb_bsaonly[!duplicated(TRYP_xlBSA_hek293t_bsadb_bsaonly[c("PositionA","PositionB")]),]
write.table(TRYP_xlBSA_hek293t_bsadb_dup,"20171010_TRYP_xlBSA_hek293t_bsadb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
TRYP_xlBSA_hek293t_bsadb_dup$Accession.A=as.character(TRYP_xlBSA_hek293t_bsadb_dup$Accession.A)
TRYP_xlBSA_hek293t_bsadb_dup$Accession.A<-replace(TRYP_xlBSA_hek293t_bsadb_dup$Accession.A,TRYP_xlBSA_hek293t_bsadb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
TRYP_xlBSA_hek293t_bsadb_dup$Accession.B=as.character(TRYP_xlBSA_hek293t_bsadb_dup$Accession.B)
TRYP_xlBSA_hek293t_bsadb_dup$Accession.B<-replace(TRYP_xlBSA_hek293t_bsadb_dup$Accession.B,TRYP_xlBSA_hek293t_bsadb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
TRYP_xlBSA_hek293t_bsadb_xinet<-data.frame(TRYP_xlBSA_hek293t_bsadb_dup$Max.XlinkX.Score,as.factor(TRYP_xlBSA_hek293t_bsadb_dup$Accession.A),as.factor(TRYP_xlBSA_hek293t_bsadb_dup$Accession.B),
TRYP_xlBSA_hek293t_bsadb_dup$PositionA,TRYP_xlBSA_hek293t_bsadb_dup$PositionB)
colnames(TRYP_xlBSA_hek293t_bsadb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(TRYP_xlBSA_hek293t_bsadb_xinet,"20171010_TRYP_xlBSA_hek293t_bsadb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
TRYP1_xlBSA_hek293t_bsadb$PositionA=pmax(TRYP1_xlBSA_hek293t_bsadb$Position.A,TRYP1_xlBSA_hek293t_bsadb$Position.B)
TRYP1_xlBSA_hek293t_bsadb$PositionB=pmin(TRYP1_xlBSA_hek293t_bsadb$Position.A,TRYP1_xlBSA_hek293t_bsadb$Position.B)
TRYP2_xlBSA_hek293t_bsadb$PositionA=pmax(TRYP2_xlBSA_hek293t_bsadb$Position.A,TRYP2_xlBSA_hek293t_bsadb$Position.B)
TRYP2_xlBSA_hek293t_bsadb$PositionB=pmin(TRYP2_xlBSA_hek293t_bsadb$Position.A,TRYP2_xlBSA_hek293t_bsadb$Position.B)
TRYP3_xlBSA_hek293t_bsadb$PositionA=pmax(TRYP3_xlBSA_hek293t_bsadb$Position.A,TRYP3_xlBSA_hek293t_bsadb$Position.B)
TRYP3_xlBSA_hek293t_bsadb$PositionB=pmin(TRYP3_xlBSA_hek293t_bsadb$Position.A,TRYP3_xlBSA_hek293t_bsadb$Position.B)
#Filter XlinkX Score >80
TRYP1_xlBSA_hek293t_bsadb=TRYP1_xlBSA_hek293t_bsadb[TRYP1_xlBSA_hek293t_bsadb[,"Max.XlinkX.Score"]>s5,]
TRYP2_xlBSA_hek293t_bsadb=TRYP2_xlBSA_hek293t_bsadb[TRYP2_xlBSA_hek293t_bsadb[,"Max.XlinkX.Score"]>s5,]
TRYP3_xlBSA_hek293t_bsadb=TRYP3_xlBSA_hek293t_bsadb[TRYP3_xlBSA_hek293t_bsadb[,"Max.XlinkX.Score"]>s5,]
#Remove Duplicates
TRYP1_xlBSA_hek293t_bsadb_dup=TRYP1_xlBSA_hek293t_bsadb[!duplicated(TRYP1_xlBSA_hek293t_bsadb[c("PositionA","PositionB")]),]
TRYP2_xlBSA_hek293t_bsadb_dup=TRYP2_xlBSA_hek293t_bsadb[!duplicated(TRYP2_xlBSA_hek293t_bsadb[c("PositionA","PositionB")]),]
TRYP3_xlBSA_hek293t_bsadb_dup=TRYP3_xlBSA_hek293t_bsadb[!duplicated(TRYP3_xlBSA_hek293t_bsadb[c("PositionA","PositionB")]),]
#TRYP Number of BSA crosslinks
TRYP1_xlBSA_hek293t_bsadb_bsaonly<- TRYP1_xlBSA_hek293t_bsadb_dup[grepl("P02769", TRYP1_xlBSA_hek293t_bsadb_dup[["Accession.A"]]) & grepl("P02769", TRYP1_xlBSA_hek293t_bsadb_dup[["Accession.B"]]), ]
TRYP2_xlBSA_hek293t_bsadb_bsaonly<- TRYP2_xlBSA_hek293t_bsadb_dup[grepl("P02769", TRYP2_xlBSA_hek293t_bsadb_dup[["Accession.A"]]) & grepl("P02769", TRYP2_xlBSA_hek293t_bsadb_dup[["Accession.B"]]), ]
TRYP3_xlBSA_hek293t_bsadb_bsaonly<- TRYP3_xlBSA_hek293t_bsadb_dup[grepl("P02769", TRYP3_xlBSA_hek293t_bsadb_dup[["Accession.A"]]) & grepl("P02769", TRYP3_xlBSA_hek293t_bsadb_dup[["Accession.B"]]), ]
#TRYP Number of all crosslinks, no duplicates
TRYP_xlBSA_hek293t_bsadb_crosslinks=c(nrow(TRYP1_xlBSA_hek293t_bsadb_dup),nrow(TRYP2_xlBSA_hek293t_bsadb_dup),nrow(TRYP3_xlBSA_hek293t_bsadb_dup))
#TRYP Number of BSA crosslinks, no duplicates
TRYP_xlBSA_hek293t_bsadb_bsaonly=c(nrow(TRYP1_xlBSA_hek293t_bsadb_bsaonly),nrow(TRYP2_xlBSA_hek293t_bsadb_bsaonly),nrow(TRYP3_xlBSA_hek293t_bsadb_bsaonly))
#TRYP Number of non-BSA crosslinks
TRYP_xlBSA_hek293t_bsadb_nonbsa=(TRYP_xlBSA_hek293t_bsadb_crosslinks-TRYP_xlBSA_hek293t_bsadb_bsaonly)
#Output for Number of crosslinks
TRYP_xlBSA_hek293t_bsadb=data.frame(TRYP_xlBSA_hek293t_bsadb_crosslinks,TRYP_xlBSA_hek293t_bsadb_bsaonly,TRYP_xlBSA_hek293t_bsadb_nonbsa)
colnames(TRYP_xlBSA_hek293t_bsadb)=c("TRYP_xlBSA_hek293t_bsadb_Total","TRYP_xlBSA_hek293t_bsadb_bsacrosslinks","TRYP_xlBSA_hek293t_bsadb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(TRYP1_xlBSA_hek293t_bsadb_bsaonly,TRYP2_xlBSA_hek293t_bsadb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(TRYP1_xlBSA_hek293t_bsadb_bsaonly,TRYP3_xlBSA_hek293t_bsadb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(TRYP2_xlBSA_hek293t_bsadb_bsaonly,TRYP3_xlBSA_hek293t_bsadb_bsaonly,by=c("PositionA","PositionB")))
tiff('TRYP_xlBSA_hek293t_bsadb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(TRYP1_xlBSA_hek293t_bsadb_bsaonly),area2=nrow(TRYP2_xlBSA_hek293t_bsadb_bsaonly),area3=nrow(TRYP3_xlBSA_hek293t_bsadb_bsaonly),
n12=nrow(merge(TRYP1_xlBSA_hek293t_bsadb_bsaonly,TRYP2_xlBSA_hek293t_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(TRYP1_xlBSA_hek293t_bsadb_bsaonly,TRYP3_xlBSA_hek293t_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(TRYP2_xlBSA_hek293t_bsadb_bsaonly,TRYP3_xlBSA_hek293t_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(TRYP_xlBSA_hek293t_bsadb_dup)-(nrow(TRYP1_xlBSA_hek293t_bsadb_bsaonly)+nrow(TRYP2_xlBSA_hek293t_bsadb_bsaonly)+nrow(TRYP3_xlBSA_hek293t_bsadb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("TRYP_xlBSA_hek293t_bsadb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()

##GluC cross-linked BSA+hek293T, BSA database
GLUC1_xlBSA_hek293t_bsadb=read.table("170920_xlBSA_hek293t_gluC_1_Crosslinks.txt", fill=TRUE, header=TRUE)
GLUC2_xlBSA_hek293t_bsadb=read.table("170920_xlBSA_hek293t_gluC_2_Crosslinks.txt", fill=TRUE, header=TRUE)
GLUC3_xlBSA_hek293t_bsadb=read.table("170920_xlBSA_hek293t_gluC_3_Crosslinks.txt", fill=TRUE, header=TRUE)
#Identify replicate number
GLUC1_xlBSA_hek293t_bsadb$Replicate=c(1)
GLUC2_xlBSA_hek293t_bsadb$Replicate=c(2)
GLUC3_xlBSA_hek293t_bsadb$Replicate=c(3)
##Compile all replicates
GLUC_xlBSA_hek293t_bsadb=rbind(GLUC1_xlBSA_hek293t_bsadb,GLUC2_xlBSA_hek293t_bsadb,GLUC3_xlBSA_hek293t_bsadb)
write.table(GLUC_xlBSA_hek293t_bsadb,"20171010_GLUC_xlBSA_hek293t_bsadb.txt",sep="\t",row.names=FALSE)
##Reposition, PositionA>PositionB
GLUC_xlBSA_hek293t_bsadb$PositionA=pmax(GLUC_xlBSA_hek293t_bsadb$Position.A,GLUC_xlBSA_hek293t_bsadb$Position.B)
GLUC_xlBSA_hek293t_bsadb$PositionB=pmin(GLUC_xlBSA_hek293t_bsadb$Position.A,GLUC_xlBSA_hek293t_bsadb$Position.B)
##Keep BSA cross-links, duplicates kept
GLUC_xlBSA_hek293t_bsadb_bsaonly<- GLUC_xlBSA_hek293t_bsadb[grepl("P02769", GLUC_xlBSA_hek293t_bsadb[["Accession.A"]]) & grepl("P02769", GLUC_xlBSA_hek293t_bsadb[["Accession.B"]]), ]
write.table(GLUC_xlBSA_hek293t_bsadb_bsaonly,"20171010_GLUC_xlBSA_hek293t_bsadb_BSAonly.txt",sep="\t",row.names=FALSE)
##Non-BSA cross-links, duplicates kept
GLUC_xlBSA_hek293t_bsadb_nonbsa<- GLUC_xlBSA_hek293t_bsadb[!grepl("P02769", GLUC_xlBSA_hek293t_bsadb[["Accession.A"]]) | !grepl("P02769", GLUC_xlBSA_hek293t_bsadb[["Accession.B"]]), ]
write.table(GLUC_xlBSA_hek293t_bsadb_nonbsa,"20171010_GLUC_xlBSA_hek293t_bsadb_nonBSAonly.txt",sep="\t",row.names=FALSE)
##Remove duplicates after keeping BSA cross-links
GLUC_xlBSA_hek293t_bsadb_dup=GLUC_xlBSA_hek293t_bsadb_bsaonly[!duplicated(GLUC_xlBSA_hek293t_bsadb_bsaonly[c("PositionA","PositionB")]),]
write.table(GLUC_xlBSA_hek293t_bsadb_dup,"20171010_GLUC_xlBSA_hek293t_bsadb_BSAonly_nodup.txt",sep="\t",row.names=FALSE)
#For XiNET export
GLUC_xlBSA_hek293t_bsadb_dup$Accession.A=as.character(GLUC_xlBSA_hek293t_bsadb_dup$Accession.A)
GLUC_xlBSA_hek293t_bsadb_dup$Accession.A<-replace(GLUC_xlBSA_hek293t_bsadb_dup$Accession.A,GLUC_xlBSA_hek293t_bsadb_dup$Accession.A=="P02769","sp|P02769|ALBU_BOVIN ")
GLUC_xlBSA_hek293t_bsadb_dup$Accession.B=as.character(GLUC_xlBSA_hek293t_bsadb_dup$Accession.B)
GLUC_xlBSA_hek293t_bsadb_dup$Accession.B<-replace(GLUC_xlBSA_hek293t_bsadb_dup$Accession.B,GLUC_xlBSA_hek293t_bsadb_dup$Accession.B=="P02769","sp|P02769|ALBU_BOVIN ")
GLUC_xlBSA_hek293t_bsadb_xinet<-data.frame(GLUC_xlBSA_hek293t_bsadb_dup$Max.XlinkX.Score,as.factor(GLUC_xlBSA_hek293t_bsadb_dup$Accession.A),as.factor(GLUC_xlBSA_hek293t_bsadb_dup$Accession.B),
GLUC_xlBSA_hek293t_bsadb_dup$PositionA,GLUC_xlBSA_hek293t_bsadb_dup$PositionB)
colnames(GLUC_xlBSA_hek293t_bsadb_xinet)=c("Score","Protein1","Protein2","LinkPos1","LinkPos2")
write.csv(GLUC_xlBSA_hek293t_bsadb_xinet,"20171010_GLUC_xlBSA_hek293t_bsadb_xinet.csv",row.names=FALSE)
###For individual replicates
#Reposition, PositionA>PositionB
GLUC1_xlBSA_hek293t_bsadb$PositionA=pmax(GLUC1_xlBSA_hek293t_bsadb$Position.A,GLUC1_xlBSA_hek293t_bsadb$Position.B)
GLUC1_xlBSA_hek293t_bsadb$PositionB=pmin(GLUC1_xlBSA_hek293t_bsadb$Position.A,GLUC1_xlBSA_hek293t_bsadb$Position.B)
GLUC2_xlBSA_hek293t_bsadb$PositionA=pmax(GLUC2_xlBSA_hek293t_bsadb$Position.A,GLUC2_xlBSA_hek293t_bsadb$Position.B)
GLUC2_xlBSA_hek293t_bsadb$PositionB=pmin(GLUC2_xlBSA_hek293t_bsadb$Position.A,GLUC2_xlBSA_hek293t_bsadb$Position.B)
GLUC3_xlBSA_hek293t_bsadb$PositionA=pmax(GLUC3_xlBSA_hek293t_bsadb$Position.A,GLUC3_xlBSA_hek293t_bsadb$Position.B)
GLUC3_xlBSA_hek293t_bsadb$PositionB=pmin(GLUC3_xlBSA_hek293t_bsadb$Position.A,GLUC3_xlBSA_hek293t_bsadb$Position.B)
#Filter XlinkX Score >80
GLUC1_xlBSA_hek293t_bsadb=GLUC1_xlBSA_hek293t_bsadb[GLUC1_xlBSA_hek293t_bsadb[,"Max.XlinkX.Score"]>s6,]
GLUC2_xlBSA_hek293t_bsadb=GLUC2_xlBSA_hek293t_bsadb[GLUC2_xlBSA_hek293t_bsadb[,"Max.XlinkX.Score"]>s6,]
GLUC3_xlBSA_hek293t_bsadb=GLUC3_xlBSA_hek293t_bsadb[GLUC3_xlBSA_hek293t_bsadb[,"Max.XlinkX.Score"]>s6,]
#Remove Duplicates
GLUC1_xlBSA_hek293t_bsadb_dup=GLUC1_xlBSA_hek293t_bsadb[!duplicated(GLUC1_xlBSA_hek293t_bsadb[c("PositionA","PositionB")]),]
GLUC2_xlBSA_hek293t_bsadb_dup=GLUC2_xlBSA_hek293t_bsadb[!duplicated(GLUC2_xlBSA_hek293t_bsadb[c("PositionA","PositionB")]),]
GLUC3_xlBSA_hek293t_bsadb_dup=GLUC3_xlBSA_hek293t_bsadb[!duplicated(GLUC3_xlBSA_hek293t_bsadb[c("PositionA","PositionB")]),]
#GLUC Number of BSA crosslinks
GLUC1_xlBSA_hek293t_bsadb_bsaonly<- GLUC1_xlBSA_hek293t_bsadb_dup[grepl("P02769", GLUC1_xlBSA_hek293t_bsadb_dup[["Accession.A"]]) & grepl("P02769", GLUC1_xlBSA_hek293t_bsadb_dup[["Accession.B"]]), ]
GLUC2_xlBSA_hek293t_bsadb_bsaonly<- GLUC2_xlBSA_hek293t_bsadb_dup[grepl("P02769", GLUC2_xlBSA_hek293t_bsadb_dup[["Accession.A"]]) & grepl("P02769", GLUC2_xlBSA_hek293t_bsadb_dup[["Accession.B"]]), ]
GLUC3_xlBSA_hek293t_bsadb_bsaonly<- GLUC3_xlBSA_hek293t_bsadb_dup[grepl("P02769", GLUC3_xlBSA_hek293t_bsadb_dup[["Accession.A"]]) & grepl("P02769", GLUC3_xlBSA_hek293t_bsadb_dup[["Accession.B"]]), ]
#GLUC Number of all crosslinks, no duplicates
GLUC_xlBSA_hek293t_bsadb_crosslinks=c(nrow(GLUC1_xlBSA_hek293t_bsadb_dup),nrow(GLUC2_xlBSA_hek293t_bsadb_dup),nrow(GLUC3_xlBSA_hek293t_bsadb_dup))
#GLUC Number of BSA crosslinks, no duplicates
GLUC_xlBSA_hek293t_bsadb_bsaonly=c(nrow(GLUC1_xlBSA_hek293t_bsadb_bsaonly),nrow(GLUC2_xlBSA_hek293t_bsadb_bsaonly),nrow(GLUC3_xlBSA_hek293t_bsadb_bsaonly))
#GLUC Number of non-BSA crosslinks
GLUC_xlBSA_hek293t_bsadb_nonbsa=(GLUC_xlBSA_hek293t_bsadb_crosslinks-GLUC_xlBSA_hek293t_bsadb_bsaonly)
#Output for Number of crosslinks
GLUC_xlBSA_hek293t_bsadb=data.frame(GLUC_xlBSA_hek293t_bsadb_crosslinks,GLUC_xlBSA_hek293t_bsadb_bsaonly,GLUC_xlBSA_hek293t_bsadb_nonbsa)
colnames(GLUC_xlBSA_hek293t_bsadb)=c("GLUC_xlBSA_hek293t_bsadb_Total","GLUC_xlBSA_hek293t_bsadb_bsacrosslinks","GLUC_xlBSA_hek293t_bsadb_nonbsacrosslinks")
#Plotting Venn Diagram
n12=nrow(merge(GLUC1_xlBSA_hek293t_bsadb_bsaonly,GLUC2_xlBSA_hek293t_bsadb_bsaonly,by=c("PositionA","PositionB")))
n13=nrow(merge(GLUC1_xlBSA_hek293t_bsadb_bsaonly,GLUC3_xlBSA_hek293t_bsadb_bsaonly,by=c("PositionA","PositionB")))
n23=nrow(merge(GLUC2_xlBSA_hek293t_bsadb_bsaonly,GLUC3_xlBSA_hek293t_bsadb_bsaonly,by=c("PositionA","PositionB")))
tiff('GLUC_xlBSA_hek293t_bsadb_replicates.tiff',width=10, height=10, units="in", res=300)
g=draw.triple.venn(area1=nrow(GLUC1_xlBSA_hek293t_bsadb_bsaonly),area2=nrow(GLUC2_xlBSA_hek293t_bsadb_bsaonly),area3=nrow(GLUC3_xlBSA_hek293t_bsadb_bsaonly),
n12=nrow(merge(GLUC1_xlBSA_hek293t_bsadb_bsaonly,GLUC2_xlBSA_hek293t_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n13=nrow(merge(GLUC1_xlBSA_hek293t_bsadb_bsaonly,GLUC3_xlBSA_hek293t_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n23=nrow(merge(GLUC2_xlBSA_hek293t_bsadb_bsaonly,GLUC3_xlBSA_hek293t_bsadb_bsaonly,by=c("PositionA","PositionB"))),
n123=nrow(GLUC_xlBSA_hek293t_bsadb_dup)-(nrow(GLUC1_xlBSA_hek293t_bsadb_bsaonly)+nrow(GLUC2_xlBSA_hek293t_bsadb_bsaonly)+nrow(GLUC3_xlBSA_hek293t_bsadb_bsaonly)-n12-n13-n23),
category=c("Replicate1","Replicate2","Replicate3"),lwd=2,
fill=c("blue","red","green"),cex=3,cat.cex=2,cat.fontface=2,euler.d=FALSE,scaled=FALSE)
grid.arrange(gTree(children=g),top=textGrob("GLUC_xlBSA_hek293t_bsadb_Replicates",gp=gpar(fontsize=32,fontface=2)))
dev.off()