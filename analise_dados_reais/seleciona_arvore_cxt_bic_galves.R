setwd('C:/Users/rique/Desktop/PIPGES/Dissertacao/Códigos/Comparacao_Estimadores')
source("algoritmos/algoritmo_contexto.R")
source("algoritmos/algoritmo_bic.R")
source("algoritmos/algoritmo_galves.R")

sequencia_binaria <- function(cadeia, bin) {
  breaks <- seq(0, floor(max(cadeia)) + 1, by = bin)
  categorias <- cut(cadeia[,1], breaks = breaks, right = TRUE, labels = FALSE, include.lowest = TRUE)
  binario <- rep(0, length(breaks) - 1)
  binario[unique(na.omit(categorias))] <- 1
  return(binario)
}

dados_neuronio <- read.table("Conjunto de dados/cell_30.dat")
sequencia_vet <- sequencia_binaria(dados_neuronio, 0.029)
sequencia <- paste0(sequencia_vet, collapse = '')
log(nchar(sequencia)) # tamanho da árcore maximal

### Contexto ----
alfabeto <- c('0','1')
limiar <- log(nchar(sequencia))
tam_max <- log(nchar(sequencia))

estimacoes_contexto  <- list()
c <- 0.1
estimacoes_contexto[[1]] <- tau_c(sequencia, alfabeto, c*limiar, tam_max)
estimacoes_contexto[[1]]$constante <- c
for(i in 2:300){
  c <- c(1:300/10)[i]
  estimacoes_contexto[[i]] <- tau_c_calculado(estimacoes_contexto[[1]]$probabilidades, 
                                              estimacoes_contexto[[1]]$deltas, 
                                              estimacoes_contexto[[1]]$internos, 
                                              c*limiar)
  estimacoes_contexto[[i]]$constante <- c
  print(c)
  print(estimacoes_contexto[[i]]$hat_tau)
  if(length(estimacoes_contexto[[i]]$hat_tau) ==  2) break
}

estimacoes_contexto_limpa <- list()
arvores <- c()
cont <- 0
for(i in 1:length(estimacoes_contexto)){
  arvore_selecionada <- paste(estimacoes_contexto[[i]]$hat_tau, collapse = '-')
  if (arvore_selecionada %in% arvores) next
  arvores <- c(arvores, arvore_selecionada)
  cont <- cont + 1
  estimacoes_contexto_limpa[[cont]] <- estimacoes_contexto[[i]]
}

saveRDS(estimacoes_contexto_limpa, file = "analise_dados_reais/dados/selecoes_contexto_neuro_30.rds")


### BIC ----

estimacoes_bic <- list()
c <- 0.1
estimacoes_bic[[1]] <- tau_bic(sequencia, alfabeto, tam_max, c)
estimacoes_bic[[1]]$constante <- c

for(i in 2:300){
  c <- c(1:300/10)[i]
  print(c)
  estimacoes_bic[[i]] <- tau_bic_calculado(sequencia, alfabeto,
                                           estimacoes_bic[[1]]$log_P_ws, 
                                           estimacoes_bic[[1]]$probabilidades, 
                                           estimacoes_bic[[1]]$folhas, 
                                           c)
  estimacoes_bic[[i]]$constante <- c
  print(estimacoes_bic[[i]]$hat_tau)
  if(length(estimacoes_bic[[i]]$hat_tau) == 2) break
}

estimacoes_bic_limpa <- list()
arvores <- c()
cont <- 0
for(i in 1:length(estimacoes_bic)){
  arvore_selecionada <- paste(estimacoes_bic[[i]]$hat_tau, collapse = '-')
  if (arvore_selecionada %in% arvores) next
  arvores <- c(arvores, arvore_selecionada)
  cont <- cont + 1
  estimacoes_bic_limpa[[cont]] <- estimacoes_bic[[i]]
}
saveRDS(estimacoes_bic_limpa, file = "analise_dados_reais/dados/selecoes_bic_neuro_30.rds")

### Galves ----
tam_max <- 8

estimacoes_galves <- list()
c <- 0.05
estimacoes_galves[[1]] <- tau_g(sequencia, alfabeto, c, tam_max)
estimacoes_galves[[1]]$constante <- c

for(i in 2:20){
  c <- c(1:20/20)[i]
  print(c)
  estimacoes_galves[[i]] <- tau_g_calculado(estimacoes_galves[[1]]$deltas,
                                            estimacoes_galves[[1]]$nos,
                                            estimacoes_galves[[1]]$probabilidades, alfabeto, c)
  estimacoes_galves[[i]]$constante <- c
  print(estimacoes_galves[[i]]$hat_tau)
  if(length(estimacoes_galves[[i]]$hat_tau) == 0) break
}

estimacoes_galves_limpa <- list()
arvores <- c()
cont <- 0
for(i in 1:length(estimacoes_galves)){
  arvore_selecionada <- paste(estimacoes_galves[[i]]$hat_tau, collapse = '-')
  if (arvore_selecionada %in% arvores) next
  arvores <- c(arvores, arvore_selecionada)
  cont <- cont + 1
  estimacoes_galves_limpa[[cont]] <- estimacoes_galves[[i]]
}

saveRDS(estimacoes_galves_limpa, file = "analise_dados_reais/dados/selecoes_galves_neuro_30.rds")

