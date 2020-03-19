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
clS <- predict(bl, readMergeS, BLAST_args = '-num_threads 12 -evalue 1e-100') #swamp
clOH <- predict(bl, readMergeOH, BLAST_args = '-num_threads 12 -evalue 1e-100') #O'Hara


