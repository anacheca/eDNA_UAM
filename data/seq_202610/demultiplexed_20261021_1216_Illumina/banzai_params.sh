#!/usr/bin/env bash


################################################################################
# LIBRARY PREP
################################################################################
# CS_TAGS  : samples carry CS1/CS2 tags inside the reads (demultiplexed here)
# ILLUMINA : samples only have i5/i7 indices, so the sequencer already split
#            the reads into one R1/R2 pair per sample, Fwd primer on R1
LIBRARY_PREP="ILLUMINA"

################################################################################
# INPUT
################################################################################
# What is the file path to the directory containing all of the libraries/reads?

PARENT_DIR="/home/moncho/Projects/Novogene_X208SC26088993/raw"

# Where is the sequencing metadata file? (SEE FORMATTING GUIDELINES IN README!)
SEQUENCING_METADATA="${PARENT_DIR}"/../metadata_caban.csv

################################################################################
# OUTPUT
################################################################################
# This script will generate a directory (folder) containing the output of the script.
# Where do you want this new folder to go?/home/moncho/Projects/Novogene_X208SC26088993
OUTPUT_DIRECTORY="/home/moncho/Projects/Novogene_X208SC26088993/pipeline_output"


################################################################################
# METADATA DETAILS
################################################################################
# Specify columns for raw sequencing files:
COLNAME_FILE1="file1"
COLNAME_FILE2="file2"

# One row per sample and locus. Sample + locus must be unique.
# Sample names: letters, digits, . _ - only. Locus names: no underscores.
COLNAME_SAMPLE_ID="Sample_name"

################################################################################
# PRIMER REMOVAL
################################################################################
# Specify the primers used to generate these amplicons.
# They are read from SEQUENCING_METADATA: give the column names of the
# forward and reverse primers and of the locus.
COLNAME_PRIMER1="PrimerF"
COLNAME_PRIMER2="PrimerR"
COLNAME_LOCUS="Locus"

# Read pairs where either read is shorter than this after primer removal are
# discarded (removes primer dimers). Reads shorter than LENR1 / LENR2 are lost
# later anyway when DADA2 truncates them.
MIN_LENGTH=50

################################################################################
# USE HASH
################################################################################
# Should the sequence ID after dereplication be the output of a hash algorithm?

USE_HASH="YES"

################################################################################
# CLUSTER OTUs: USING DADA2, Unoise and secondary swarm
################################################################################

SEARCH_ASVs="YES"

SEARCH_Unoise="NO"

SECONDARY_SWARM="YES"
### Truncation lengths after primer removal (DADA2 truncLen, UNOISE cutadapt -l / -L)
# For MiFish with 2x150 reads: R1 keeps ~129 bp after the primer, R2 ~123 bp,
# and the amplicon is ~165-185 bp, so 120 + 110 leaves ~50 bp of overlap.
LENR1=100
LENR2=95

################################################################################
# REANALYSIS
################################################################################
# Would you like to pick up where a previous analysis left off?

# Reuse the primer-trimmed reads of a previous run of this pipeline?
# Point DEMULT_OUTPUT to that run's output folder (it must contain
# sample_trans.tmp, summary.csv and the noprimers folder; run it with HOARD="YES").

ALREADY_DEMULTIPLEXED="NO"
DEMULT_OUTPUT=""


################################################################################
# GENERAL SETTINGS
################################################################################
# Would you like to save every single intermediate file as we go? YES | NO
# recommendation: NO, unless testing or troubleshooting
HOARD="NO"

