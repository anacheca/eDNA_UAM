#!/usr/bin/env bash


################################################################################
# INPUT
################################################################################
# What is the file path to the directory containing all of the libraries/reads?

PARENT_DIR="/mnt/c/Users/AM.5074073/OneDrive - UAM/Documentos/Medits_test_git/data/raw_data_tfg"

# Where is the sequencing metadata file? (SEE FORMATTING GUIDELINES IN README!)
SEQUENCING_METADATA="${PARENT_DIR}"/metadata_modificado.csv

################################################################################
# OUTPUT
################################################################################
# This script will generate a directory (folder) containing the output of the script.
# Where do you want this new folder to go?
OUTPUT_DIRECTORY="/mnt/c/Users/AM.5074073/OneDrive - UAM/Documentos/Medits_test_git/pipeline_output_tfg"


################################################################################
# METADATA DETAILS
################################################################################
# Specify columns for raw sequencing files:
COLNAME_FILE1="file1"
COLNAME_FILE2="file2"

# MUST be unique for each row!
COLNAME_SAMPLE_ID="Sample"


# Your metadata must have a column corresponding to the subfolders containing the raw reads.
# In order to make this flexible across both multiple and single library preps, you must include this even if you only sequenced one library (sorry!).
COLNAME_ID1_NAME="pri_index_name"

################################################################################
# DEMULTIPLEXING
################################################################################


# Specify the nucleotide sequences that differentiate multiplexed samples
# (sometimes, confusingly referred to as "tags" or "barcodes")
# these are the secondary index -- the primary index added with the sequencing adapters should not be in the sequence data
# You can grab these from the file specified above (SEQUENCING_METADATA) by specifying the column name of index sequences.
COLNAME_ID2_SEQ="seq_index_seq_1"
COLNAME_ID2_SEQ_R="seq_index_seq_2"
COLNAME_ID2_NAME="seq_index_name_1"
COLNAME_ID2_NAME_R="seq_index_name_2"


################################################################################
# PRIMER REMOVAL
################################################################################
# Specify the primers used to generate these amplicons.
# As with the multiplex indexes, Banzai will grab these from the file SEQUENCING_METADATA.
# You must indicate the column names of the forward and reverse primers
COLNAME_PRIMER1="primerF_seq"
COLNAME_PRIMER2="primerR_seq"
COLNAME_LOCUS="locus"

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
### Parameters for DADA2
LENR1=190
LENR2=160

################################################################################
# REANALYSIS
################################################################################
# Would you like to pick up where a previous analysis left off?

# Have you already demultiplexed your reads into .1 and .2 pairs per sample.
# Point towards the output folder (must include files: sample_trans.tmp,
# barcodes.fasta, summary.csv and pcr_primers.fasta; and the folder /demultiplexed
# so the pipeline can cp all necessary files

ALREADY_DEMULTIPLEXED="NO"
DEMULT_OUTPUT=""


################################################################################
# GENERAL SETTINGS
################################################################################
# Would you like to save every single intermediate file as we go? YES | NO
# recommendation: NO, unless testing or troubleshooting
HOARD="YES"

