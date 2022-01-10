
# PPI networks for each group ---------------------------------------------

library(biomaRt)
library(tidyverse)
library(igraph)
library(tidygraph)
library(scatterpie)
library(ggraph)
library(magrittr)
library(RedeR)

# ------------------ GET GENES AND INTERACTIONS -------------------------
# Load diff genes table
load("results/diff_exp/diff_df.rda")
gwas_intersections <- read_csv("results/tables/gwas_intersection.csv")

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

# Get interaction network
identifiers <- paste0(unique(diff_df$hgnc_symbol), collapse = "%0d")
s_map <- get_map(identifiers)

# Select channels (from experimental data, coexpression, and database) and combine their scores
s_map %>% 
  mutate(across(contains("score"), function(i) round(i * 1000, digits = 0))) %>% 
  dplyr::select(contains("preferredName"), contains("score")) %>% 
  combinescores(evidences = c("escore", "ascore", "dscore"), confLevel = 0.4) %>% 
  unique() -> int

# RedeR

# Here we select proper coordinates to our network
g <- graph_from_edgelist(as.matrix(int[,1:2]), directed = F)
rdp <- RedeR::RedPort()
calld(rdp)
addGraph(rdp, g)

# Using REDER, we selected the best visualization for our network.
# Nodes and edges coordinates are saved as nodes2 and edges2 in 'results/networks' directory.

# Plot network layoout ----------------------------------------------------

nodes <- read_tsv("results/networks/nodes2")
edges <- read_delim("results/networks/edges2", delim = "\t")


# Mark nodes as GWAS
nodes %<>% 
  mutate(gwas = ifelse(alias %in% gwas_genes, "gwas", "not_gwas"))

# Plot
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
    pie_scale = 1,
  show.legend = F) +
  geom_scatterpie(
    cols = c("a", "b", "c"),
    data = as_data_frame(g, "vertices") %>% filter(gwas == "not_gwas"),
    colour = NA,
    pie_scale = 0.18,
    show.legend = F
  ) +
  geom_node_text(aes(label = alias), size = 1.5, nudge_x = 10, nudge_y = 10) + 
  #geom_node_label(aes(label = alias)) + 
  scale_fill_manual(values = c("#0ac80aff", "#4f4affff", "#ff822fff")) +
  coord_fixed() +
  theme_graph() -> p

# Save 
svg(filename = "results/plots_paper/rede2.svg", height = 20, width = 20)
print(p)
dev.off()

# Percentage of total genes in the network: 49,43%%
n_distinct(nodes$alias) / n_distinct(diff_df$hgnc_symbol)

