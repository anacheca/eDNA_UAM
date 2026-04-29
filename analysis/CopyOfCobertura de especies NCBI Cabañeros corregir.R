############################################################
# 0. LIBRERÍAS
############################################################

library(tidyverse)
library(readxl)
library(rentrez)
library(stringr)
library(Biostrings)

############################################################
# 1. CARGA DE ESPECIES
############################################################

species_df <- read_excel("~/eDNA_UAM/data/Cabañeros lista de especies.xlsx") %>%
  select(Genero_especie) %>%
  filter(!is.na(Genero_especie)) %>%
  mutate(Genero_especie = str_trim(Genero_especie))

############################################################
# 2. OBTENER taxID (NCBI)
############################################################

get_taxid <- function(sp){
  res <- entrez_search(db="taxonomy", term=sp, retmax=1)
  if(length(res$ids)==0) return(NA)
  res$ids[[1]]
}

species_df <- species_df %>%
  mutate(taxID = map_chr(Genero_especie, get_taxid))

#Revisar cuántas especies perdemos al hacer el filter
# Hay 6 sin taxid, añadir manualment (5 subspecies)

missing_taxid <-  species_df %>%
  filter(is.na(taxID))

species_df <- species_df %>%
  filter(!is.na(taxID))

############################################################
# 3. BUSCAR ACCESIONES 12S (ROBUSTO Y SEGURO)
############################################################

search_12S_accessions <- function(taxid){
  
  if(is.na(taxid)) return(character(0))
  
  query <- paste0(
    "txid", taxid,
    "[ORGN] AND (12S OR 12S ribosomal OR rRNA) OR mitochondrion"
  )
  
  res <- tryCatch(
    entrez_search(
      db = "nucleotide",
      term = query,
      retmax = 2000
    ),
    error = function(e) return(NULL)
  )
  
  if(is.null(res) || length(res$ids) == 0) return(character(0))
  
  as.character(res$ids)
}

species_df <- species_df %>%
  mutate(
    accessions_12S = map(taxID, search_12S_accessions),
    n_12S = map_int(accessions_12S, length)
  )

############################################################
# 4. LIMPIEZA + DESCARGA DE SECUENCIAS (VERSION SEGURA + DEBUG)
############################################################

# 🔴 LIMPIEZA DE ACCESIONES
species_df$accessions_12S <- lapply(species_df$accessions_12S, function(x){
  if(is.null(x)) return(character(0))
  x <- x[!is.na(x)]
  as.character(x)
})

# 🔴 EVITAR BLOQUEOS DE NCBI
options(timeout = 60)

############################################################
# 🔥 FUNCIÓN FINAL (FILTRA IDS PROBLEMÁTICOS + SOLO 12S ÚTIL)
############################################################

fetch_sequences_12S_primers <- function(ids){
  
  if(is.null(ids)) return(NULL)
  
  ids <- as.character(ids)
  ids <- ids[!is.na(ids)]
  
  if(length(ids) == 0) return(NULL)
  
  # 🔴 LIMITE REAL (evita colapsos)
  ids <- head(ids, 30)
  
  seq_list <- list()
  
  for(id in ids){
    
    # 🔴 FILTRO CLAVE: elimina IDs problemáticos tipo contigs raros
    if(grepl("^[0-9]{10,}$", id)) next
    
    cat("probando ID:", id, "\n")
    
    seqs <- tryCatch(
      entrez_fetch(
        db = "nucleotide",
        id = id,
        rettype = "fasta",
        retmode = "text"
      ),
      error = function(e) return(NULL)
    )
    
    # 🔴 control de errores NCBI
    if(is.null(seqs)) next
    if(!is.character(seqs)) next
    if(nchar(seqs) == 0) next
    
    seqs <- gsub("\r", "", seqs)
    
    # 🔴 asegurar FASTA válido
    if(!grepl(">", seqs)) next
    
    # 🔴 conversión segura Biostrings
    dna <- tryCatch(
      Biostrings::readDNAStringSet(textConnection(seqs)),
      error = function(e) return(NULL)
    )
    
    if(!is.null(dna)){
      seq_list[[length(seq_list) + 1]] <- dna
    }
  }
  
  if(length(seq_list) == 0) return(NULL)
  
  do.call(c, seq_list)
}

############################################################
# 🔴 DESCARGA FINAL (SIN map, EVITA BLOQUEOS)
############################################################

seqs_12S <- vector("list", nrow(species_df))

for(i in seq_len(nrow(species_df))){
  
  cat("\nESPECIE", i, "de", nrow(species_df), "\n")
  
  res <- fetch_sequences_12S_primers(
    species_df$accessions_12S[[i]]
  )
  
  # 🔴 CLAVE: convertir NULL en lista vacía explícita
  if(is.null(res)){
    seqs_12S[[i]] <- Biostrings::DNAStringSet()
  } else {
    seqs_12S[[i]] <- res
  }
}

# 🔴 SEGURIDAD EXTRA: forzar estructura constante
stopifnot(length(seqs_12S) == nrow(species_df))

species_df$seqs_12S <- seqs_12S

############################################################
# 5. VARIABILIDAD GENÉTICA (CORREGIDA)
############################################################

calc_variability <- function(dna){
  
  if(is.null(dna)) return(NA)
  if(!inherits(dna,"DNAStringSet")) return(NA)
  if(length(dna) < 2) return(NA)
  
  lengths <- width(dna)
  
  sd(lengths, na.rm=TRUE)
}

species_df <- species_df %>%
  mutate(var_length_12S = map_dbl(seqs_12S, calc_variability))

############################################################
# 6. PRIMERS (SOLO 12S)
############################################################

FWD <- "GGGTTGGTAAATTTCGTGCCAGC"
REV <- "CATAGTGGGGTATCTAATCCCAGTTTG"

check_primers <- function(dna){
  
  if(is.null(dna)) return(FALSE)
  
  fwd <- vcountPattern(FWD, dna, max.mismatch=3)
  rev <- vcountPattern(REV, dna, max.mismatch=3)
  
  any(fwd > 0 & rev > 0)
}

species_df <- species_df %>%
  mutate(amplicon_12S = map_lgl(seqs_12S, check_primers))

############################################################
# 7. COBERTURA
############################################################

species_df <- species_df %>%
  mutate(
    coverage_12S = case_when(
      n_12S == 0 ~ "0",
      n_12S <= 5 ~ "baja",
      n_12S <= 20 ~ "media",
      TRUE ~ "alta"
    )
  )

############################################################
# 8. TABLA FINAL
############################################################

final_table <- species_df %>%
  select(
    species = Genero_especie,
    taxID,
    n_12S,
    coverage_12S,
    var_length_12S,
    amplicon_12S
  )

write_csv(final_table, "Report_12S.csv")

############################################################
# 9. RESUMEN GLOBAL
############################################################

summary_table <- final_table %>%
  group_by(coverage_12S, amplicon_12S) %>%
  summarise(n_species = n(), .groups="drop")

write_csv(summary_table, "Summary_report.csv")
