############################################################
# 🔹 0. LIBRERÍAS
############################################################

library(tidyverse)
library(readxl)
library(rentrez)
library(stringr)
library(ape)

species <- read_excel("~/eDNA_UAM/data/Cabañeros lista de especies.xlsx") %>%
  select(Genero_especie) %>%                 # columna correcta
  filter(!is.na(Genero_especie)) %>%         # quitar NA
  mutate(Genero_especie = str_trim(Genero_especie))


############################################################
# 🔹 3. CREAR CARPETAS
############################################################

dir.create("data/COI", recursive = TRUE, showWarnings = FALSE)
dir.create("data/12S", recursive = TRUE, showWarnings = FALSE)


############################################################
# 🔹 4. FUNCIÓN DESCARGA NCBI
############################################################

download_mt <- function(sp, gene, folder){
  
  cat("Descargando:", sp, "|", gene, "\n")
  
  query <- paste0('"', sp, '"[ORGN] AND ', gene)
  
  res <- entrez_search(db = "nucleotide", term = query, retmax = 200)
  
  if(length(res$ids) == 0){
    cat("  ⚠️ sin resultados\n")
    return(NULL)
  }
  
  seqs <- entrez_fetch(
    db = "nucleotide",
    id = res$ids,
    rettype = "fasta",
    retmode = "text"
  )
  
  safe <- str_replace_all(sp, " ", "_")
  
  writeLines(seqs, file.path(folder, paste0(safe, ".fasta")))
  
  Sys.sleep(0.3)   # evita bloqueo NCBI
}


############################################################
# 🔹 5. LISTA DE GENES
############################################################

COI_query  <- "COI[Gene]"
T12S_query <- "12S[TITL]"


############################################################
# 🔹 6. DESCARGA (SEPARADO Y CORRECTO)
############################################################

walk(species$Genero_especie,
     ~download_mt(.x, COI_query, "data/COI"))

walk(species$Genero_especie,
     ~download_mt(.x, T12S_query, "data/12S"))

############################################################
# 🔹 6. PRIMERS
############################################################

COI_FWD <- "GGGTTGGTAAATTTCGTGCCAGC"
COI_REV <- "CATAGTGGGGTATCTAATCCCAGTTTG"

S12_FWD <- COI_FWD
S12_REV <- COI_REV

############################################################
# 🔹 7. FUNCIÓN ROBUSTA PRIMERS
############################################################

check_primers <- function(file, FWD, REV){
  
  seqs <- tryCatch(
    ape::read.dna(file, format="fasta"),
    error = function(e) return(NULL)
  )
  
  # ❌ archivo vacío o roto
  if(is.null(seqs) || NROW(seqs) == 0){
    return(tibble(
      species = str_remove(basename(file), "\\.fasta"),
      n_seq = 0,
      amplifica = FALSE
    ))
  }
  
  # 🔥 CONVERSIÓN SEGURA DNAbin → strings
  seqs_char <- as.character(seqs)
  
  tibble(
    species = str_remove(basename(file), "\\.fasta"),
    n_seq = length(seqs_char),
    amplifica = any(
      str_detect(seqs_char, FWD) &
        str_detect(seqs_char, REV)
    )
  )
}

############################################################
# 🔹 8. ANALISIS COI
############################################################

files_coi <- list.files("data/COI", full.names=TRUE)

results_coi <- map_dfr(files_coi,
                       ~check_primers(.x, COI_FWD, COI_REV)) %>%
  mutate(marker = "COI")

############################################################
# 🔹 9. ANALISIS 12S
############################################################

files_12s <- list.files("data/12S", full.names=TRUE)

results_12s <- map_dfr(files_12s,
                       ~check_primers(.x, S12_FWD, S12_REV)) %>%
  mutate(marker = "12S")

############################################################
# 🔹 10. RESULTADO FINAL
############################################################

results <- bind_rows(results_coi, results_12s)

############################################################
# 🔹 11. RESUMEN
############################################################

results %>%
  count(marker, amplifica)

############################################################
# 🔹 12. PLOT FINAL
############################################################

ggplot(results, aes(x=marker, fill=amplifica)) +
  geom_bar(position="fill") +
  facet_wrap(~marker)
