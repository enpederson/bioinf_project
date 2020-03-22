#!/usr/bin/R

#notes:
#dark matter definition: E=-100
# teasing out lowest common ancestor between reads by traversing phylogeny tables
# kraken? unmatched reads
# correlate by quality scores and lengths (filter)
# QueryID correspond to fastq IDs
# repeats? (HMM) repeat database? unique k-mers (low complexity DUST MASKER)
# tblastx (lower error rate?)
# pick all the same genes at first from SRA
# 
# SRA codes:
# Swamp: SRX7820012, SRX7820011, SRX7820010 
# OHara: SRX7820009, SRX7820008, SRX7820007

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# pkgs <- rownames(installed.packages())
# BiocManager::install(pkgs, type = "source", checkBuilt = TRUE)
# BiocManager::install("ShortRead")
library(rBLAST)
library(ShortRead)
library(ggplot2)
library(taxonomizr)
library(dplyr)
library(forcats)
library(tibble)
library(tidyr)
library(stringr)
library(ape)
library(Biostrings)
library(DECIPHER)

options(scipen=999, digits=1) 


Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/share/pkg.7/blast+/2.7.1/install/bin", sep=":"))
Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/share/pkg.7/sratoolkit/2.9.2/install/bin/", sep=":"))

#access local BLAST database
bl <- blast(db="/projectnb/ct-shbioinf/blast/nt.fa") 

#read Swamp SRRs:
srrS1=c('SRX7820012')
srrS2=c('SRX7820011')
srrS3=c('SRX7820010')
#load fastqs
system(paste('fastq-dump', srrS1, sep=' ')) 
system(paste('fastq-dump', srrS2, sep=' '))
system(paste('fastq-dump', srrS3, sep=' '))
#read fastqs
dnaS1 = readFastq('.', pattern=srrS1)
dnaS2 = readFastq('.', pattern=srrS2)
dnaS3 = readFastq('.', pattern=srrS3)

#combine Swamp reads
S123 = c(dnaS1,dnaS2,dnaS3)
dnaMergeS <- do.call('append', S123)

#readS1 = sread(dnaS1, id=id(dnaS1))
#readS2 = sread(dnaS2, id=id(dnaS2))
#readS3 = sread(dnaS3, id=id(dnaS3))
readMergeS = sread(dnaMergeS, id=id(dnaMergeS))


#read O'Hara SRRs:
srrOH1=c('SRX7820009')
srrOH2=c('SRX7820008')
srrOH3=c('SRX7820007')
#load fastqs
system(paste('fastq-dump', srrOH1, sep=' ')) 
system(paste('fastq-dump', srrOH2, sep=' '))
system(paste('fastq-dump', srrOH3, sep=' '))
#read fastqs
dnaOH1 = readFastq('.', pattern=srrOH1)
dnaOH2 = readFastq('.', pattern=srrOH2)
dnaOH3 = readFastq('.', pattern=srrOH3)

OH123 = c(dnaOH1,dnaOH2,dnaOH3)
dnaMergeOH <- do.call('append', OH123)

#readOH1 = sread(dnaOH1, id=id(dnaOH1))
#readOH2 = sread(dnaOH2, id=id(dnaOH2))
#readOH3 = sread(dnaOH3, id=id(dnaOH3))
readMergeOH = sread(dnaMergeOH, id=id(dnaMergeOH))



#This will take a long time. Run once and save environment!
#clS <- predict(bl, readMergeS, BLAST_args = '-num_threads 12 -evalue 1e-100') #swamp
#clOH <- predict(bl, readMergeOH, BLAST_args = '-num_threads 12 -evalue 1e-100') #O'Hara

#load the taxonomy database files 'nodes' and 'names'
taxaNodes<-read.nodes.sql("/projectnb/ct-shbioinf/taxonomy/data/nodes.dmp")
taxaNames<-read.names.sql("/projectnb/ct-shbioinf/taxonomy/data/names.dmp")
# Search the taxonomy by accession ID #

accidS = as.character(clS$SubjectID) # accession IDs of BLAST hits
accidOH = as.character(clOH$SubjectID)

idS<-accessionToTaxa(accidS,'/projectnb/ct-shbioinf/taxonomy/data/accessionTaxa.sql')
idOH<-accessionToTaxa(accidOH,'/projectnb/ct-shbioinf/taxonomy/data/accessionTaxa.sql')

# displays the taxonomic names from each ID #
taxlistS=getTaxonomy(idS, taxaNodes, taxaNames)
taxlistOH=getTaxonomy(idOH, taxaNodes, taxaNames)

cltaxS=cbind(clS,taxlistS)
cltaxOH=cbind(clOH,taxlistOH)

#memory reduction
# remove(taxlistS)
# remove(taxlistOH)

### Saved data as /data/Env1.RData

topgenera = cltaxS %>%
  group_by(genus) %>% 
  #top_n(n=1, wt=(Bits)) %>%
  #group_by(species) %>% 
  #group_by(genus) %>%
  summarize(count=n()) #%>% 
#arrange(desc(count))
head(arrange(topgenera, desc(count)),n=10)

cltaxS_name = rownames_to_column(cltaxS,var="readID")
HitsQIDS = cltaxS_name %>%
  group_by(QueryID)
  
HitsQIDNAS = cltaxS_name %>%
  filter(is.na(genus))  %>%
  group_by(QueryID)
  
# filter if QueryID group has only genus=NA
HitsQID_group_S = cltaxS %>%
  filter(!is.na(genus)) %>% # get rid of NA genera
  group_by(QueryID) %>% #group
  summarize(count=n())

#Find missing QueryIDs
lenS = lengths(HitsQID_group_S) 
seq1S = seq(1:lenS[1])

#make artificial list with length of HitsQID_group_S
regexp <- "[[:digit:]]+"
seq2S = str_extract(HitsQID_group_S[[1,]], regexp)

missingintS = setdiff(seq1S,seq2S) #QueryIDs that have no matches
missingS = paste("Query_", sep="", missingintS)
NAreadsS = readMergeS[missingintS] #READS OBJECT WITH NO BLAST HITS


HitsNAOH = cltaxOH %>%
  filter(is.na(genus)) %>%
  arrange(desc(Bits))

# filter if QueryID group has only genus=NA
HitsQID_group_OH = cltaxOH %>%
  filter(!is.na(genus)) %>% # get rid of NA genera
  group_by(QueryID) %>% #group
  summarize(count=n())

#Find missing QueryIDs
lenOH = lengths(HitsQID_group_OH) 
seq1OH = seq(1:lenOH[1])

#make artificial list with length of HitsQID_group_S
regexp <- "[[:digit:]]+"
seq2OH = str_extract(HitsQID_group_OH[[1,]], regexp)

missingintOH = setdiff(seq1OH,seq2OH) #QueryIDs that have no matches
missingOH = paste("Query_", sep="", missingintOH)
NAreadsOH = readMergeS[missingintOH]

#Complexity Analysis
complexS = dustyScore(NAreadsS)
complexdfS = as.data.frame(complexS)
complexdfS_full = cbind(NAreadsS,complexdfS)

ggplot(data = complexdfS, aes(x = complexS)) +
  ggtitle("Complexity Scores of Dark Matter from Swamp samples") +
  xlab("DUST complexity score") +
  geom_histogram(bins = 250) +
  scale_x_log10()

#sequester out low and high complexity reads
complexdfS_high = complexdfS_full %>%
  filter(complexS >= 5000)

complexdfS_low = complexdfS_full %>%
  filter(complexS < 5000)


complexOH = dustyScore(NAreadsOH)
complexdfOH = as.data.frame(complexOH)

ggplot(data = complexdfOH, aes(x = complexOH)) +
  ggtitle("Complexity Scores of Dark Matter from O'Hara samples") +
  xlab("DUST complexity score") +
  geom_histogram(bins = 250) +
  scale_x_log10()

#Multiple alignments ############

pause()
# DNA_alignment_S <- AlignSeqs(NAreadsS, processors=12)
# writeXStringSet(DNA_alignment_S, file="./outputs/Alignments_S.fasta")

DNA_alignment_high_S <- AlignSeqs(complexdfS_high, processors=12)
DNA_alignment_high_staggered_S = StaggerAlignment(DNA_alignment_high_S, processors=12)
writeXStringSet(DNA_alignment_high_staggered_S, file="./outputs/Alignments_high_staggered_S.fasta")

###########
pause()
DNA_alignment_low_S <- AlignSeqs(complexdfS_low, processors=12)
DNA_alignment_low_staggered_S = StaggerAlignment(DNA_alignment_low_S, processors=12)
writeXStringSet(DNA_alignment_low_staggered_S, file="./outputs/Alignments_low_staggered_S.fasta")


##########
pause()
DNA_alignment_OH <- AlignSeqs(NAreadsOH, processors=12)
writeXStringSet(DNA_alignment_OH, file="./outputs/Alignments_OH.fasta")

IUPAC_CODE_MAP # list of ambiguity codes for reference
threshold = 0.05
DNA_consensus_S <- ConsensusSequence(DNA_alignment_S,
                  threshold = 0.05,
                  ambiguity = TRUE,
                  noConsensusChar = "+",
                  minInformation = 1 - threshold,
                  ignoreNonBases = FALSE,
                  includeTerminalGaps = FALSE)
DNA_consensus_S

