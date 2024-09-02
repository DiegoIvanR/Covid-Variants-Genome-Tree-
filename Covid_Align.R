# Respuestas escritas al final del script
# Diego Iván Rodríguez Núñez A01644772
# Luis Ignacio 

library(Biostrings)  
library(seqinr)  
library(adegenet)  
library(ape)  
library(DECIPHER)  
library(viridis)  
library(ggplot2)  
library(phangorn)
library(ggtree)
# PP754897 - USA 
# PP756396 - USA
# PP298099 - MX 
# OR856797 - Brazil 
# PP754055 - Russia
# OQ059020 - Russia
# MT637143 - Russia
# MT320538 - France
# MT318827 - Germany
# LC529905 - Japan
# MT039890 - South Korea
# MW672391 - Austria
# ON506959 - Italy
# MT873892 - UK
# MT007544 - Australia
# MT670023 - Chile
# MT499210 - Poland
# MW309438 - Canada
# OL672836 - Belgium
# ON671727 - Zealand
# LR883967 - Spain
# OK180974 - Czech Republic
# MZ145384 - Serbia
# OM643294 - Romania

genbank_ids <- c("PP754897", "PP756396", "PP298099", "OR856797", "PP754055", "OQ059020",
                 "MT637143","MT320538", "MT318827", "LC529905", "MT039890", "MW672391",
                 "ON506959", "MT873892", "MT007544", "MT670023", "MT499210", "MW309438",
                 "OL672836", "ON671727", "LR883967", "OK180974", "MZ145384", "OM643294")

sequences <- read.GenBank(genbank_ids)
output_file <- "output_file.fasta"
write.dna(sequences, file = output_file, format = "fasta")
sequence_lengths <- sapply(sequences, length)


not_align <- readDNAStringSet("output_file.fasta", format = "fasta")
not_align <- OrientNucleotides(not_align)
align <- AlignSeqs(not_align)

BrowseSeqs(align, highlight = 0)
writeXStringSet(align, file = "covid_align.fasta")
fasta_file <- "covid_align.fasta"
dna_alignment <- read.dna(fasta_file, format = "fasta", as.matrix = TRUE)

numSecuencias <- 24
numBases <- 5
count_bases <- function(seq) {
  alphabetFrequency(seq, baseOnly = TRUE, as.prob = FALSE)
}
base_counts <- count_bses(not_align[1])
datos <- data.frame(bases=rep(c('A', 'C', 'G', 'T', 'other'), times=numSecuencias),
                    virus=c(rep(genbank_ids[1], times=numBases),
                            rep(genbank_ids[2], times=numBases),
                            rep(genbank_ids[3], times=numBases),
                            rep(genbank_ids[4], times=numBases),
                            rep(genbank_ids[5], times=numBases),
                            rep(genbank_ids[6], times=numBases),
                            rep(genbank_ids[7], times=numBases),
                            rep(genbank_ids[8], times=numBases),
                            rep(genbank_ids[9], times=numBases),
                            rep(genbank_ids[10], times=numBases),
                            rep(genbank_ids[11], times=numBases),
                            rep(genbank_ids[12], times=numBases),
                            rep(genbank_ids[13], times=numBases),
                            rep(genbank_ids[14], times=numBases),
                            rep(genbank_ids[15], times=numBases),
                            rep(genbank_ids[16], times=numBases),
                            rep(genbank_ids[17], times=numBases),
                            rep(genbank_ids[18], times=numBases),
                            rep(genbank_ids[19], times=numBases),
                            rep(genbank_ids[20], times=numBases),
                            rep(genbank_ids[21], times=numBases),
                            rep(genbank_ids[22], times=numBases),
                            rep(genbank_ids[23], times=numBases),
                            rep(genbank_ids[24], times=numBases)
                    ),
                    tots=c(count_bases(not_align[1]),
                           count_bases(not_align[2]),
                           count_bases(not_align[3]),
                           count_bases(not_align[4]),
                           count_bases(not_align[5]),
                           count_bases(not_align[6]),
                           count_bases(not_align[7]),
                           count_bases(not_align[8]),
                           count_bases(not_align[9]),
                           count_bases(not_align[10]),
                           count_bases(not_align[11]),
                           count_bases(not_align[12]),
                           count_bases(not_align[13]),
                           count_bases(not_align[14]),
                           count_bases(not_align[15]),
                           count_bases(not_align[16]),
                           count_bases(not_align[17]),
                           count_bases(not_align[18]),
                           count_bases(not_align[19]),
                           count_bases(not_align[20]),
                           count_bases(not_align[21]),
                           count_bases(not_align[22]),
                           count_bases(not_align[23]),
                           count_bases(not_align[24])
                    ))

bases <- ggplot(datos, aes(x = bases, y = tots, fill = bases))+
  geom_col() + #geom_bar(stat = "identity")
  facet_wrap(~ virus, ncol = 4) + 
  scale_fill_brewer(palette="Set1") +
  ggtitle("Composicion de nucleótidos de cada virus") +
  labs(x="Bases", y="Total de bases", fill="virus") +
  theme_dark()

hamming_dist <- dist.dna(dna_alignment, model = "raw")
hamming_matrix <- as.matrix(hamming_dist)
n <- nrow(hamming_matrix)
sequence_names <- rownames(hamming_matrix)
distance_data <- expand.grid(Secuencia1 = sequence_names, Secuencia2 = sequence_names)
distance_data$Distancia <- as.vector(hamming_matrix)

distancia_plot <- ggplot(distance_data, aes(x = Secuencia1, y = Secuencia2, fill = Distancia)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +  # Gradiente de gris
  labs(x = "Secuencia", y = "Secuencia", fill = "Distancia") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotar etiquetas

#ggsave("distancia.png", distancia_plot)

# 12: Crear el árbol filogenético con ape
tree <- nj(hamming_dist)
ggtree_plot <- ggtree(tree) + 
  geom_tiplab() +
  aes()

# Imprimir las longitudes de las secuencias
cat("Longitudes de las secuencias:\n")
for (i in seq_along(sequence_lengths)) {
  cat(sprintf("Secuencia %d: %d bases\n", i, sequence_lengths[i]))
}
distancia_plot

bases # Coincide en las variantes de los virus que la base timina
# es la que más abunda, seguida de adenina, guanina, y al final,
# citosina.
ggtree_plot #  árbol filogenético que muestra las relaciones 
# evolutivas inferidas entre las diferentes variantes del
# SARS-CoV-2 basadas en sus características físicas o genéticas.
# Dentro de la gráfica se puede apreciar tres grupos comunes,
# que se identifican como las regiones del planeta, principalmente
# Asia, Europa y Américas, por lo que si se puede intuir que hay
# una relación entre los grupos de variantes en sí, pero que a su 
# vez difieren de las otras regiones.
