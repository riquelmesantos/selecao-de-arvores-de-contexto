library(ggraph)
library(tidygraph)
library(igraph)
library(stringi)

plot_arvore <- function(arvore, label = T){
  arvore <- stri_reverse(arvore)
  arvore <- arvore[order(unlist(arvore), decreasing = FALSE)]
  init_g = graph.empty(n = 1, directed = TRUE)
  V(init_g)$name <- ""
  g = init_g
  for (word in arvore) {
    subwords <- stri_reverse(stri_sub(word, 1, 1:nchar(word)))
    subg <- graph.lattice(length(subwords) + 1, directed = TRUE)
    V(subg)$name  <-  c("", subwords)
    g = igraph::union(g, subg)
  }
  
  g_tbl <- as_tbl_graph(g)
  # Sem labels
  if(!label) return(ggraph(g_tbl, layout = 'tree', circular = FALSE) + 
    geom_edge_link(color = "black", width = 0.5) + 
    geom_node_point(color = "black", size = 1.5) + 
    theme_void())
  
  # Com labels
  folhas <- degree(g, mode = "out") == 0
  g_tbl <- g_tbl %>% mutate(label = ifelse(folhas, name, NA))
  return(
    ggraph(g_tbl, layout = 'tree', circular = FALSE) + 
    geom_edge_link(color = "black", width = 1.5) + 
    geom_node_point(color = "black", size = 3) + 
    
    geom_node_text(aes(label = label),
                   size = 6,
                   vjust = 1.8,
                   hjust = 0.5,
                   family = "serif",
                   na.rm = TRUE) +
    theme_void()+
    theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
    )
}

# plot_arvore(
#   estimacoes_bic_limpa[[1]]$hat_tau,
#   label = F
# )
plot(plot_arvore(estimacoes_contexto_limpa[[i]]$hat_tau, T))
# for(i in 17:19){
#   pdf(file = sprintf("desenhos/com_label/arvore_%d.pdf", i), width = 7, height = 7)
#   plot(plot_arvore(estimacoes_contexto_limpa[[i]]$hat_tau, T))
#   dev.off()
# }



# pdf("arvore.pdf", width = 10, height = 10)
# plot_arvore(
#   estimacoes_contexto_limpa[[1]]$hat_tau,
#   label = T
# )
# dev.off()















