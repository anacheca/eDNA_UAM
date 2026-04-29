############################################################
# 0. LIBRERÍAS
############################################################

library(tidyverse)
library(readxl)
library(rentrez)
library(stringr)

############################################################
# 1. CARGA DE ESPECIES (tu lista del parque nacional)
############################################################

# Leemos la lista de especies del Excel
species <- read_excel("~/eDNA_UAM/data/Cabañeros lista de especies.xlsx") %>%
  
  # Nos quedamos solo con el nombre especie (género + especie)
  select(Genero_especie) %>%
  
  # Eliminamos valores vacíos
  filter(!is.na(Genero_especie)) %>%
  
  # Limpiamos espacios extra
  mutate(Genero_especie = str_trim(Genero_especie))

############################################################
# 2. FUNCIÓN DE BÚSQUEDA EN NCBI
############################################################

# Esta función busca cuántas secuencias hay por especie y marcador
search_ncbi <- function(sp, query){
  
  # Construimos la query tipo GenBank
  term <- paste0('"', sp, '"[ORGN] AND ', query)
  
  # Buscamos en base de datos nucleotide
  res <- entrez_search(
    db = "nucleotide",
    term = term,
    retmax = 0   # SOLO queremos el número de resultados (más rápido)
  )
  
  # Devolvemos tabla tipo "coverage"
  tibble(
    species = sp,
    counts = res$count
  )
}

############################################################
# 3. DEFINIMOS LOS MARCADORES
############################################################

# COI (barcode animal)
COI_query <- "COI[Gene]"

# 12S (muy usado en eDNA vertebrados)
T12S_query <- "12S[TITL]"

############################################################
# 4. BÚSQUEDA EN NCBI
############################################################

# COI: para cada especie
COI_raw <- map_dfr(
  species$Genero_especie,
  ~search_ncbi(.x, COI_query)
)

# 12S: para cada especie
T12S_raw <- map_dfr(
  species$Genero_especie,
  ~search_ncbi(.x, T12S_query)
)

############################################################
# 5. FORMATEO DE DATOS
############################################################

# Añadimos etiqueta de locus
COI_raw <- COI_raw %>%
  mutate(locus = "COI")

T12S_raw <- T12S_raw %>%
  mutate(locus = "12S")

############################################################
# 6. UNIMOS TODO
############################################################

to.plot <- bind_rows(COI_raw, T12S_raw)

############################################################
# 7. CLASIFICAMOS COBERTURA
############################################################

to.plot <- to.plot %>%
  
  # Creamos categorías tipo "+2", "1", "0"
  mutate(
    coverage = case_when(
      counts == 0 ~ "0",
      counts == 1 ~ "1",
      counts > 1  ~ "+2"
    )
  ) %>%
  
  # Ordenamos niveles (como fct_relevel del original)
  mutate(coverage = factor(coverage, levels = c("0", "1", "+2")))

############################################################
# 8. VISUALIZACIÓN
############################################################

ggplot(to.plot, aes(x = locus)) +
  
  # Barras por nivel de cobertura
  geom_bar(aes(fill = coverage)) +
  
  # Cada especie en un panel (como faceting en Myctophidae)
  facet_wrap(~species, scales = "free") +
  
  # Estilo tipo script original
  coord_flip() +
  
  theme(
    panel.spacing = unit(0, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(size = 6),
    axis.title = element_blank(),
    axis.text.x = element_blank()
  ) +
  
  labs(
    title = "COI vs 12S coverage en especies de Cabañeros",
    x = "",
    y = ""
  )

