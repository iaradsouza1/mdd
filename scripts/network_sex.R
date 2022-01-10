
# PPI networks for each group ---------------------------------------------

library(biomaRt)
library(tidyverse)
library(igraph)
library(tidygraph)
library(scatterpie)
library(ggraph)
library(magrittr)

# ------------------ GET GENES AND INTERACTIONS -------------------------
# Load diff genes table
load("results/diff_exp/diff_df.rda")
gwas_intersections <- read_csv("results/tables/gwas_intersection.csv")

diff_df %<>% 
  separate(group, into = c("region", "sex"), sep = "_")

gwas_genes <- unique(gwas_intersections$gene_name)

# FUNCTIONS ---------------------------------------------------------------

# Get interaction network from string database
get_map  <- function(ids) {
  s_map <- read.table(
    text = RCurl::postForm(
      "https://string-db.org/api/tsv/network"
      ,identifiers = ids
      ,echo_query  = "1"
      ,required_score = "0"
      ,show_query_node_labels = "1"
      ,species     = "9606"
    )[[1]], sep = "\t", header = T)
}

# Function to combine scores from different channels
combinescores <- function(dat, evidences = "all", confLevel = 0.4){
  
  if(evidences[1] == "all"){
    edat <- dat[, -c(1, 2, ncol(dat))]
  } else {
    if(!all(evidences %in% colnames(dat))) {
      stop("NOTE: one or more 'evidences' not listed in 'dat' colnames!")
    }
    edat <- dat[, evidences]
  }
  edat <- edat/1000
  
  edat <- 1 - edat
  sc <- apply(X = edat, MARGIN = 1, FUN = function(x) 1 - prod(x))
  dat <- cbind(dat[, c(1, 2)], combined_score = sc)
  idx <- dat$combined_score >= confLevel
  dat <- dat[idx,]
  return(dat)
  
}

# FEMALE ------------------------------------------------------------------

# Select the female altered genes
diff_df_female <- diff_df %>% filter(sex == "female")

# Get interaction network
ids_female <- paste0(unique(diff_df_female$hgnc_symbol), collapse = "%0d")
s_map_female <- get_map(ids_female)

# Select channels (from experimental data, coexpression, and database) and combine their scores
s_map_female %>% 
  mutate(across(contains("score"), function(i) round(i * 1000, digits = 0))) %>% 
  dplyr::select(contains("preferredName"), contains("score")) %>% 
  combinescores(evidences = c("escore", "ascore", "dscore"), confLevel = 0.4) %>% 
  unique() -> int_female

# RedeR
# Here we select proper coordinates to our network
# This a manual task, we choose the best layout to the network
g <- graph_from_edgelist(as.matrix(int_female[,1:2]), directed = F)
rdp <- RedeR::RedPort()
calld(rdp)
addGraph(rdp, g)

# Using RedeR interface, we save the network layout and used it in the next steps

# Plot network layoout ----------------------------------------------------

# nodes <- read_tsv("data/nodes1")
# edges <- read_delim("data/edges1", delim = "\t")

nodes <- read_tsv("~/Área de trabalho/redes/v3/node_female1")
edges <- read_delim("~/Área de trabalho/redes/v3/edge_female1", delim = "\t")

nodes %<>% 
  mutate(gwas = ifelse(alias %in% gwas_genes, "gwas", "not_gwas"))

g <- graph_from_data_frame(edges, nodes, directed = F)

V(g)$x <- nodes$x
V(g)$y <- nodes$y
V(g)$size <- 1
V(g)$a <- ifelse(V(g)$alias %in% unique(diff_df$hgnc_symbol[diff_df$type == "DGE"]), 1, 0)
V(g)$b <- ifelse(V(g)$alias %in% unique(diff_df$hgnc_symbol[diff_df$type == "DTE"]), 1, 0)
V(g)$c <- ifelse(V(g)$alias %in% unique(diff_df$hgnc_symbol[diff_df$type == "DTU"]), 1, 0)

ggraph(g, x = x, y = y) + 
  geom_edge_link0(edge_color = "gray", alpha = 0.7, width = 0.4) +
  geom_scatterpie(
    cols = c("a", "b", "c"),
    data = as_data_frame(g, "vertices") %>% filter(gwas == "gwas"),
    colour = NA,
    n = 5,
    pie_scale = 1.5,
    show.legend = F) +
  geom_scatterpie(
    cols = c("a", "b", "c"),
    data = as_data_frame(g, "vertices") %>% filter(gwas == "not_gwas"),
    colour = NA,
    pie_scale = 0.22,
    show.legend = F
  ) +
  geom_node_text(aes(label = alias), size = 1, nudge_x = 20, nudge_y = 20) + 
  scale_fill_manual(values = c("#0ac80aff", "#4f4affff", "#ff822fff")) +
  coord_fixed() +
  theme_graph() -> p

# ggsave(p, filename = "~/Área de trabalho/redes/v3/rede2.pdf",
#        height = 20, width = 20, device = cairo_pdf, dpi = 300, units = "cm")

svg(filename = "~/Área de trabalho/redes/v3/rede_female.svg", height = 10, width = 10)
print(p)
dev.off()

# Percentage of total genes in the network: 41,17%%
n_distinct(nodes$alias) / n_distinct(diff_df_female$hgnc_symbol)

# MALE --------------------------------------------------------------------

diff_df_male <- diff_df %>% filter(sex == "male")

# Get interaction network
ids_male <- paste0(unique(diff_df_male$hgnc_symbol), collapse = "%0d")
s_map_male <- get_map(ids_male)

# Select channels (from experimental data, coexpression, and database) and combine their scores
s_map_male %>% 
  mutate(across(contains("score"), function(i) round(i * 1000, digits = 0))) %>% 
  dplyr::select(contains("preferredName"), contains("score")) %>% 
  combinescores(evidences = c("escore", "ascore", "dscore"), confLevel = 0.4) %>% 
  unique() -> int_male

# RedeR
g <- graph_from_edgelist(as.matrix(int_male[,1:2]), directed = F)
rdp <- RedeR::RedPort()
calld(rdp)
addGraph(rdp, g)

# Using REDER, we selected the best visualization for our network.
# Nodes and edges coordinates are saved as nodes18 and edges18.

# Plot network layout ----------------------------------------------------

# nodes <- read_tsv("data/nodes1")
# edges <- read_delim("data/edges1", delim = "\t")

nodes <- read_tsv("~/Área de trabalho/redes/v3/nodes_male1")
edges <- read_delim("~/Área de trabalho/redes/v3/edges_male1", delim = "\t")

nodes %<>% 
  mutate(gwas = ifelse(alias %in% gwas_genes, "gwas", "not_gwas"))

g <- graph_from_data_frame(edges, nodes, directed = F)

V(g)$x <- nodes$x
V(g)$y <- nodes$y
V(g)$size <- 1
V(g)$a <- ifelse(V(g)$alias %in% unique(diff_df$hgnc_symbol[diff_df$type == "DGE"]), 1, 0)
V(g)$b <- ifelse(V(g)$alias %in% unique(diff_df$hgnc_symbol[diff_df$type == "DTE"]), 1, 0)
V(g)$c <- ifelse(V(g)$alias %in% unique(diff_df$hgnc_symbol[diff_df$type == "DTU"]), 1, 0)

ggraph(g, x = x, y = y) + 
  geom_edge_link0(edge_color = "gray", alpha = 0.7, width = 0.4) +
  geom_scatterpie(
    cols = c("a", "b", "c"),
    data = as_data_frame(g, "vertices") %>% filter(gwas == "gwas"),
    colour = NA,
    n = 5,
    pie_scale = 1.5,
    show.legend = F) +
    geom_scatterpie(
      cols = c("a", "b", "c"),
      data = as_data_frame(g, "vertices") %>% filter(gwas == "not_gwas"),
      colour = NA,
      pie_scale = 0.22,
      show.legend = F
  ) +
  geom_node_text(aes(label = alias), size = 1, nudge_x = 20, nudge_y = 20) + 
  #geom_node_label(aes(label = alias)) + 
  scale_fill_manual(values = c("#0ac80aff", "#4f4affff", "#ff822fff")) +
  coord_fixed() +
  theme_graph() -> p

# ggsave(p, filename = "~/Área de trabalho/redes/v3/rede2.pdf",
#        height = 20, width = 20, device = cairo_pdf, dpi = 300, units = "cm")

svg(filename = "~/Área de trabalho/redes/v3/rede_male.svg", height = 10, width = 10)
print(p)
dev.off()

# Percentage of total genes in the network: 40,18%%
n_distinct(nodes$alias) / n_distinct(diff_df_male$hgnc_symbol)
