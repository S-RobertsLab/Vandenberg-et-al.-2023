library(plyr)
OriginalIntersect <- read.csv(file = choose.files())

WT_Intersect <- OriginalIntersect[OriginalIntersect$Genotype=="WT",]
Rad30_Intersect <- OriginalIntersect[OriginalIntersect$Genotype=="rad30D",]


WT_<- WT_Intersect[!duplicated(WT_Intersect$Start_Position),]
rad30_<- Rad30_Intersect[!duplicated(Rad30_Intersect$Start_Position),]

WT_0 <- WT_[WT_$Photoreactivation_Time_Min_=="0",]
Rad_0 <- rad30_[rad30_$Photoreactivation_Time_Min_=="0",]


WT_Rad30_Intersect <- OriginalIntersect[!duplicated(OriginalIntersect$Start_Position),]
WT_Rad30_Intersect <- WT_Rad30_Intersect[!WT_Rad30_Intersect$Genotype=="rad16D",]
WT_Rad30_Intersect <- WT_Rad30_Intersect[!WT_Rad30_Intersect$Photoreactivation_Time_Min_=="40",]

####################################
##Separate WT and Rad30 context files and then do intersect
## Merge by chr position
## WT 
WT_Intersect <- read.csv(file = choose.files())
WT_Context <- read.csv(file = choose.files())
colnames(WT_Context)[1]<- "Chr_Position"

Merged_WT <- merge(WT_Intersect,WT_Context, by="Chr_Position")

##RAD30
Rad30_Intersect <- read.csv(file = choose.files())
Rad30_Context <- read.csv(file = choose.files())
colnames(Rad30_Context)[1]<- "Chr_Position"

Merged_Rad30 <- merge(Rad30_Intersect,Rad30_Context, by="Chr_Position")

## Get sub type

Merged_Rad30$MutationTri <- Merged_Rad30$Trinucleotide
stri_sub(Merged_Rad30$MutationTri,2,2) <- substr(Merged_Rad30$Mutation,1,1)

Merged_Rad30$Ref_Sub <- stri_sub(Merged_Rad30$Trinucleotide,2,2)
Merged_Rad30$Mut_Sub <- stri_sub(Merged_Rad30$MutationTri,2,2)

Merged_Rad30$Substitution <- Merged_Rad30$Ref_Sub
stri_sub(Merged_Rad30$Substitution,2,2) <- substr(Merged_Rad30$Mut_Sub,1,1)

Merged_Rad30$Ref_Sub=NULL
Merged_Rad30$Mut_Sub=NULL
##Counts
Rad30_Counts <- ddply(Merged_Rad30,.(Merged_Rad30$Strand,Merged_Rad30$Trinucleotide,Merged_Rad30$MutationTri,Merged_Rad30$Substitution),nrow)
names(Rad30_Counts) <- c("Strand","Reference","Mutation","Substitution","Frequency")

#Separate by sub type and combine by rev comp

Rad30_TA <- Rad30_Counts[Rad30_Counts$Substitution=="TA",]
Rad30_AT <- Rad30_Counts[Rad30_Counts$Substitution=="AT",]

sum(Rad30_AT$Frequency)+sum(Rad30_TA$Frequency)
sum(WT_AT$Frequency)+ sum(WT_TA$Frequency)
WT_TA<- na.omit(WT_TA)

write.csv(Rad30_Counts,file = "Rad30_TranAsymmetryCounts.csv")

#### Fixed
#### Counts then combine

Rad30 <- read.csv(file = choose.files())
Rad26 <- read.csv(file = choose.files())
WT <- read.csv(file = choose.files())

WT_Counts <- ddply(WT,.(WT$Strand,WT$Reference,WT$Mutation,WT$Trinucleotide),nrow)
names(WT_Counts) <- c("Strand","Ref_Sub","Mut_Sub","Reference","Frequency")

WT_Counts$Ref_Sub <- stri_sub(WT_Counts$Reference,2,2)
WT_Counts$Mut_Sub <- stri_sub(WT_Counts$Mutation,2,2)

WT_Counts$Substitution <- WT_Counts$Ref_Sub
stri_sub(WT_Counts$Substitution,2,2) <- substr(WT_Counts$Mut_Sub,1,1)

WT_Counts$Ref_Sub=NULL
WT_Counts$Mut_Sub=NULL

TG_ <- WT_Counts[WT_Counts$Substitution == "TG",]
AC_ <- WT_Counts[WT_Counts$Substitution == "AC",]
AT_ <- WT_Counts[WT_Counts$Substitution == "AT",]
TA_ <- WT_Counts[WT_Counts$Substitution == "TA",]
CG_ <- WT_Counts[WT_Counts$Substitution == "CG",]
GC_ <- WT_Counts[WT_Counts$Substitution == "GC",]
AG_ <- WT_Counts[WT_Counts$Substitution == "AG",]
TC_ <- WT_Counts[WT_Counts$Substitution == "TC",]
CT_ <- WT_Counts[WT_Counts$Substitution == "CT",]
GA_ <- WT_Counts[WT_Counts$Substitution == "GA",]
CA_ <- WT_Counts[WT_Counts$Substitution == "CA",]
GT_ <- WT_Counts[WT_Counts$Substitution == "GT",]


AC_$oldStrand <- AC_$Strand
AC_$NewStrand <- AC_$Strand
AC_[AC_$NewStrand == "+", "NewStrand"] <- ".-"
AC_[AC_$NewStrand == "-", "NewStrand"] <- "+"
AC_[AC_$NewStrand == ".-", "NewStrand"] <- "-"
AC_$Strand <- AC_$NewStrand
AC_$NewStrand=NULL
AC_$oldStrand=NULL

AC_ <- na.omit(AC_)

TG_$RevComp <- sapply(TG_$Reference, function(x) as.character(reverseComplement(DNAString(x))))

AC_$RevComp <- AC_$Reference

TG_ACstrand <- merge(TG_,AC_, by=c("RevComp","Strand"), all=TRUE)
View(TG_ACstrand)
write.csv(TG_ACstrand,file="WT_TG_AC.csv")

##AT/TA##
AT_ <- SingleStrandBias[SingleStrandBias$Substitution == "AT",]
TA_ <- SingleStrandBias[SingleStrandBias$Substitution == "TA",]

AT_$oldStrand <- AT_$Strand
AT_$NewStrand <- AT_$Strand
AT_[AT_$NewStrand == "+", "NewStrand"] <- ".-"
AT_[AT_$NewStrand == "-", "NewStrand"] <- "+"
AT_[AT_$NewStrand == ".-", "NewStrand"] <- "-"
AT_$Strand <- AT_$NewStrand
AT_$NewStrand=NULL
AT_$oldStrand=NULL

AT_ <- na.omit(AT_)
TA_ <- na.omit(TA_)

TA_$RevComp <- sapply(TA_$Reference, function(x) as.character(reverseComplement(DNAString(x))))

AT_$RevComp <- AT_$Reference

TA_ATstrand <- merge(TA_,AT_, by=c("RevComp","Strand"), all=TRUE)

write.csv(TA_ATstrand,file="WT_TA_AT.csv")

##CG/GC##
CG_ <- SingleStrandBias[SingleStrandBias$Substitution == "CG",]
GC_ <- SingleStrandBias[SingleStrandBias$Substitution == "GC",]

GC_$oldStrand <- GC_$Strand
GC_$NewStrand <- GC_$Strand
GC_[GC_$NewStrand == "+", "NewStrand"] <- ".-"
GC_[GC_$NewStrand == "-", "NewStrand"] <- "+"
GC_[GC_$NewStrand == ".-", "NewStrand"] <- "-"
GC_$Strand <- GC_$NewStrand
GC_$NewStrand=NULL
GC_$oldStrand=NULL

CG_ <- na.omit(CG_)
GC_ <- na.omit(CG_)

CG_$RevComp <- sapply(CG_$Reference, function(x) as.character(reverseComplement(DNAString(x))))

GC_$RevComp <- GC_$Reference

GC_CGstrand <- merge(GC_,CG_, by=c("RevComp","Strand"), all=TRUE)

write.csv(GC_CGstrand,file="WT_GC_CG.csv")

##CT/GA##
CT_ <- SingleStrandBias[SingleStrandBias$Substitution == "CT",]
GA_ <- SingleStrandBias[SingleStrandBias$Substitution == "GA",]

GA_$oldStrand <- GA_$Strand
GA_$NewStrand <- GA_$Strand
GA_[GA_$NewStrand == "+", "NewStrand"] <- ".-"
GA_[GA_$NewStrand == "-", "NewStrand"] <- "+"
GA_[GA_$NewStrand == ".-", "NewStrand"] <- "-"
GA_$Strand <- GA_$NewStrand
AC_$NewStrand=NULL
AC_$oldStrand=NULL

CT_ <- na.omit(CT_)
GA_ <- na.omit(GA_)

GA_$RevComp <- sapply(GA_$Reference, function(x) as.character(reverseComplement(DNAString(x))))

CT_$RevComp <- CT_$Reference

CT_GAstrand <- merge(CT_,GA_, by=c("RevComp","Strand"), all=TRUE)

write.csv(CT_GAstrand,file="WT_GA_CT.csv")

##AG/TC##
AG_ <- SingleStrandBias[SingleStrandBias$Substitution == "AG",]
TC_ <- SingleStrandBias[SingleStrandBias$Substitution == "TC",]

AG_$oldStrand <- AG_$Strand
AG_$NewStrand <- AG_$Strand
AG_[AG_$NewStrand == "+", "NewStrand"] <- ".-"
AG_[AG_$NewStrand == "-", "NewStrand"] <- "+"
AG_[AG_$NewStrand == ".-", "NewStrand"] <- "-"
AG_$Strand <- AG_$NewStrand
AG_$NewStrand=NULL
AG_$oldStrand=NULL

AG_ <- na.omit(AG_)
TC_ <- na.omit(TC_)

TC_$RevComp <- sapply(TC_$Reference, function(x) as.character(reverseComplement(DNAString(x))))

AG_$RevComp <- AG_$Reference

TC_AGstrand <- merge(TC_,AG_, by=c("RevComp","Strand"), all=TRUE)
View(AG_TCstrand)

write.csv(TC_AGstrand,file="WT_AG_TC_.csv")

##CA/GT##
CA_ <- SingleStrandBias[SingleStrandBias$Substitution == "CA",]
GT_ <- SingleStrandBias[SingleStrandBias$Substitution == "GT",]

GT_$oldStrand <- GT_$Strand
GT_$NewStrand <- GT_$Strand
GT_[GT_$NewStrand == "+", "NewStrand"] <- ".-"
GT_[GT_$NewStrand == "-", "NewStrand"] <- "+"
GT_[GT_$NewStrand == ".-", "NewStrand"] <- "-"
GT_$Strand <- GT_$NewStrand

GT_$NewStrand=NULL
GT_$oldStrand=NULL
CA_ <- na.omit(CA_)
GT_ <- na.omit(GT_)

CA_$RevComp <- sapply(CA_$Reference, function(x) as.character(reverseComplement(DNAString(x))))

GT_$RevComp <- GT_$Reference

GT_CAstrand <- merge(GT_,CA_, by=c("RevComp","Strand"), all=TRUE)
View(TG_ACstrand)
write.csv(GT_CAstrand,file="WT_GT_CA_.csv")

#########DOUBLE Strand Bias#####################

WT_Counts <- ddply(WT,.(WT$Strand,WT$Substitution,WT$Context),nrow)
names(WT_Counts) <- c("Strand","Substitution","Reference","Frequency")


WT <- WTStrandBias
rad <- RadStrandBias

WTStrandBias <- read.delim(file = choose.files())
RadStrandBias <- read.delim(file = choose.files())


CCTT <- WT_Counts[WT_Counts$Substitution == "CCTT",]
GGAA <- WT_Counts[WT_Counts$Substitution == "AATG",]

GGAA$oldStrand <- GGAA$Strand
GGAA$NewStrand <- GGAA$Strand
GGAA[GGAA$NewStrand == "+", "NewStrand"] <- ".-"
GGAA[GGAA$NewStrand == "-", "NewStrand"] <- "+"
GGAA[GGAA$NewStrand == ".-", "NewStrand"] <- "-"
GGAA$Strand <- GGAA$NewStrand

GGAA$NewStrand=NULL
GGAA$oldStrand=NULL
CCTT <- na.omit(CCTT)
GGAA <- na.omit(GGAA)

CCTT$RevComp <- sapply(CCTT$Reference, function(x) as.character(reverseComplement(DNAString(x))))

GGAA$RevComp <- GGAA$Reference

GGAAstrand <- merge(GGAA,CCTT, by=c("RevComp","Strand"), all=TRUE)
View(TG_ACstrand)
write.csv(GGAAstrand,file="WT_CCTT.csv")

#RAD30
rad_Counts <- ddply(rad,.(rad$Strand,rad$Substitution,rad$Context),nrow)
names(rad_Counts) <- c("Strand","Substitution","Reference","Frequency")



CCTT <- rad_Counts[rad_Counts$Substitution == "CCTT",]
GGAA <- rad_Counts[rad_Counts$Substitution == "GGAA",]

GGAA$oldStrand <- GGAA$Strand
GGAA$NewStrand <- GGAA$Strand
GGAA[GGAA$NewStrand == "+", "NewStrand"] <- ".-"
GGAA[GGAA$NewStrand == "-", "NewStrand"] <- "+"
GGAA[GGAA$NewStrand == ".-", "NewStrand"] <- "-"
GGAA$Strand <- GGAA$NewStrand

GGAA$NewStrand=NULL
GGAA$oldStrand=NULL
CCTT <- na.omit(CCTT)
GGAA <- na.omit(GGAA)

CCTT$RevComp <- sapply(CCTT$Reference, function(x) as.character(reverseComplement(DNAString(x))))

GGAA$RevComp <- GGAA$Reference

GGAAstrand <- merge(GGAA,CCTT, by=c("RevComp","Strand"), all=TRUE)
View(TG_ACstrand)
write.csv(GGAAstrand,file="rad_CCTT.csv")
