
# Create individual networks ----------------------------------------------

genes <- c("HNRNPA3", "DDX39B")
i <- "aINS_female"

diff_df_temp <- diff_df %>% 
  select(hgnc_symbol, type) %>% 
  filter(hgnc_symbol %in% genes) %>% 
  distinct()

int_temp <- int %>% 
  filter(preferredName_A %in% diff_df_temp$hgnc_symbol, preferredName_B %in% diff_df_temp$hgnc_symbol)

g <- graph_from_edgelist(as.matrix(int_temp[,c(1,2)]))

lay <- create_layout(g, layout = "nicely")
V(g)$x <- lay[,"x"]
V(g)$y <- lay[,"y"]
V(g)$a <- ifelse(V(g)$name %in% unique(diff_df$hgnc_symbol[diff_df$type == "DGE"]), 1, 0)
V(g)$b <- ifelse(V(g)$name %in% unique(diff_df$hgnc_symbol[diff_df$type == "DTE"]), 1, 0)
V(g)$c <- ifelse(V(g)$name %in% unique(diff_df$hgnc_symbol[diff_df$type == "DTU"]), 1, 0)

ggraph(g, x = x, y = y) +
  geom_edge_link0(edge_color = "gray", alpha = 0.7, width = 0.2) +
  geom_scatterpie(
    data = as_data_frame(g, "vertices"),
    cols = c("a", "b", "c"),
    color = NA,
    pie_scale = 4,
    show.legend = F) +
  geom_node_text(aes(label = name), size = 3,  nudge_x = 0.03, nudge_y = 0.03) + 
  scale_fill_manual(values = c("#0ac80aff", "#4f4affff", "#ff822fff")) +
  coord_fixed() + 
  theme_graph() 



# -------------------------------------------------------------------------

load("results/networks/graphs_by_group_wo_nbor.rda")

make_graph_by_group <- function(genes, group_name, diff_df, int) {
  
  if(!is.null(genes)) {
    
    diff_df_temp <- diff_df %>% 
      select(hgnc_symbol, type, group) %>% 
      filter(hgnc_symbol %in% genes, group %in% group_name) %>% 
      distinct()
    
    int_temp <- int %>% 
      filter(preferredName_A %in% diff_df_temp$hgnc_symbol, preferredName_B %in% diff_df_temp$hgnc_symbol)
    
    g <- graph_from_edgelist(as.matrix(int_temp[,c(1,2)]))
    
    lay <- create_layout(g, layout = "nicely")
    V(g)$x <- lay[,"x"]
    V(g)$y <- lay[,"y"]
    V(g)$a <- ifelse(V(g)$name %in% unique(diff_df$hgnc_symbol[diff_df$type == "DGE"]), 1, 0)
    V(g)$b <- ifelse(V(g)$name %in% unique(diff_df$hgnc_symbol[diff_df$type == "DTE"]), 1, 0)
    V(g)$c <- ifelse(V(g)$name %in% unique(diff_df$hgnc_symbol[diff_df$type == "DTU"]), 1, 0)
    
    # buffer <- c(-2, 2)
    # xlims <- ceiling(range(V(g)$x)) + buffer
    # ylims <- ceiling(range(V(g)$y)) + buffer
    
    xlims <- ylims <- c(-50, 50)
    
    ggraph(g, x = x, y = y) +
      geom_edge_link0(edge_color = "gray", alpha = 0.7, width = 0.2) +
      geom_scatterpie(
        data = as_data_frame(g, "vertices"),
        cols = c("a", "b", "c"),
        color = NA,
        pie_scale = 0.5,
        show.legend = F) +
      geom_node_text(aes(label = name), size = 2,  nudge_x = 0.03, nudge_y = 0.03) + 
      scale_fill_manual(values = c("#0ac80aff", "#4f4affff", "#ff822fff")) +
      coord_fixed() +
      scale_x_continuous(limits = xlims) +
      scale_y_continuous(limits = ylims) +
      theme_graph() +
      theme(
        panel.border = element_rect(
          colour = "#161616",
          fill = NA,
          linewidth = 1
        )
      ) -> plot_g
    
    return(plot_g)
    
  } else {
    
    return(NULL)
  
  }
  
}

imap(graphs_by_group, ~ {
  
  make_graph_by_group(.x$genes_from_graph, .y, diff_df, int)

}) %>% 
  discard(is.null) -> ls_graphs

design <- "ABCDEFGH"

wrap_plots(
  ls_graphs,
  design = design,
  nrow = 2,
  ncol = 4
)
