library(ggraph)
library(tidygraph)
library(igraph)
library(stringi)
setwd('C:/Users/rique/Desktop/PIPGES/Dissertacao/Códigos/Comparacao_Estimadores')

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
    geom_node_point(color = "black", size = 0.4) + 
    theme_void())
  
  # Com labels
  folhas <- degree(g, mode = "out") == 0
  g_tbl <- g_tbl %>% mutate(label = ifelse(folhas, name, NA))
  return(
    ggraph(g_tbl, layout = 'tree', circular = FALSE) + 
    geom_edge_link(color = "black", width = 1) + 
    geom_node_point(color = "black", size = 3) + 
    
    geom_node_text(aes(label = label),
                   size = 3,
                   vjust = 1.5,
                   hjust = 0.5,
                   family = "serif",
                   na.rm = TRUE) +
    theme_void()+
    theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))+
    theme_void()
    )
}

# Carrega as arvores selecionadas
selecoes_contexto <- readRDS("analise_dados_reais/estimacoes_contexto_limpa.rds")
selecoes_bic <- readRDS("estimacoes_bic_limpa.rds")
selecoes_galves <- readRDS("analise_dados_reais/estimacoes_galves_limpa.rds")

length(selecoes_galves)
for(i in 1:11){
  pdf(file = sprintf("analise_dados_reais/desenhos/arvore_%d_galves.pdf", i), width = 7, height = 7)
  plot(plot_arvore(selecoes_bic[[4]]$hat_tau, T))
  dev.off()
}

plot(plot_arvore(estimacoes_contexto_limpa[[6]]$hat_tau, T))








