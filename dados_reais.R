library(igraph)
library(stringi)
library(BCT)
setwd('C:/Users/rique/Desktop/PIPGES/Dissertacao/Códigos/Comparacao_Estimadores')
source("algoritmos/algoritmo_contexto.R")
source("algoritmos/algoritmo_bic.R")
source("algoritmos/algoritmo_galves.R")

log_perda <- function(seq, probs) {
  ll <- 0
  tam_max <- max(nchar(rownames(probs)))
  for (i in (tam_max+1):nchar(seq)) {
    passado <- stri_sub(seq, from = i-tam_max, to = i-1)
    contexto <- intersect(sufixos(passado), rownames(probs))
    if (length(contexto) == 0) return('Não foi encontrado contexto para a sequência', passado)
    contexto <- contexto[nchar(contexto) == max(nchar(contexto))]
    
    prox_simbolo <- stri_sub(seq, from = i, to = i)
    prob <- probs[contexto, prox_simbolo]
    if(prob == 0) prob <- 1e-10 
    
    ll <- ll + log(prob)
  }
  return(-ll / (nchar(seq) - tam_max))
}

zero_um_perda <- function(seq, probs) {
  erros <- 0
  tam_max <- max(nchar(rownames(probs)))
  for (i in (tam_max+1):nchar(seq)) {
    passado <- stri_sub(seq, from = i - tam_max, to = i - 1)
    contexto <- intersect(sufixos(passado), rownames(probs))
    if (length(contexto) == 0) return('Não foi encontrado contexto para a sequência', passado)
    contexto <- contexto[nchar(contexto) == max(nchar(contexto))]
    
    verdadeiro <- stri_sub(seq, from = i, to = i)
    predito <- colnames(probs)[which.max(probs[contexto, ])]
    if (predito != verdadeiro) erros <- erros + 1
  }
  return(erros/ (nchar(seq) - tam_max))  
}

sequencia_binaria <- function(cadeia, bin) {
  breaks <- seq(0, floor(max(cadeia)) + 1, by = bin)
  categorias <- cut(cadeia[,1], breaks = breaks, right = TRUE, labels = FALSE, include.lowest = TRUE)
  binario <- rep(0, length(breaks) - 1)
  binario[unique(na.omit(categorias))] <- 1
  return(binario)
}

dados_neuronio <- read.table("Conjunto de dados/cell_2.dat")
diferenca <- diff(dados_neuronio[,1])
plot(diferenca, type = 'l')


# plota diferenças ordenadas
diferenca_ordenada <- sort(diferenca, decreasing = F)
plot(diferenca_ordenada, type = 'l')
# plota linha vertical indicando a média das diferenças
abline(h = mean(diferenca_ordenada), col = 'red', lty = 2)
nrow(dados_neuronio) - nrow(unique(dados_neuronio))

sequencia_vet <- sequencia_binaria(dados_neuronio, 0.05)
sequencia <- paste0(sequencia_vet, collapse = '')
prop_treino <- 0.7
estimacoes_contexto <- stri_sub(sequencia, to = nchar(sequencia) * prop_treino)
sequencia_teste <- stri_sub(sequencia, from = nchar(sequencia) * prop_treino + 1)
log(nchar(sequencia_treino))
plot(as.numeric(strsplit(sequencia, split = '')[[1]]))

### Contexto ----
alfabeto <- c('0','1')
limiar <- log(nchar(sequencia_treino))
tam_max <- log(nchar(sequencia_treino))

estimacoes_contexto  <- list()
c <- 0.1
estimacoes_contexto[[1]] <- tau_c(sequencia_treino, alfabeto, c*limiar, tam_max)
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

lls_contexto <- matrix(0, nrow = length(estimacoes_contexto_limpa), ncol = 3)
for(i in 1:length(estimacoes_contexto_limpa)){
  lls_contexto[i,1] <- paste(estimacoes_contexto_limpa[[i]]$hat_tau, collapse = '-')
  probs <- estimacoes_contexto_limpa[[i]]$probabilidades_arvore
  lls_contexto[i,2] <- log_perda(sequencia_teste, probs)
  lls_contexto[i,3] <- estimacoes_contexto_limpa[[i]]$constante
}
View(lls_contexto)

### BIC ----
alfabeto <- c('0','1')
tam_max <- log(nchar(sequencia_treino))
estimacoes_bic <- list()
c <- 0.1
estimacoes_bic[[1]] <- tau_bic(sequencia_treino, alfabeto, tam_max, c)
estimacoes_bic[[1]]$constante <- c
for(i in 2:300){
  c <- c(1:300/10)[i]
  print(c)
  estimacoes_bic[[i]] <- tau_bic_calculado(sequencia_treino, alfabeto,
                                           estimacoes_bic[[1]]$log_P_ws, 
                                           estimacoes_bic[[1]]$probabilidades, 
                                           estimacoes_bic[[1]]$folhas, 
                                           c)
  estimacoes_bic[[i]]$constante <- c
  print(estimacoes_bic[[i]]$hat_tau)
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

lls_bic <- matrix(0, nrow = length(estimacoes_bic_limpa), ncol = 3)

for(i in 1:length(estimacoes_bic_limpa)){
  lls_bic[i,1] <- paste(estimacoes_bic_limpa[[i]]$hat_tau, collapse = '-')
  probs <- estimacoes_bic_limpa[[i]]$probabilidades_arvore
  lls_bic[i,2] <- log_perda(sequencia_teste, probs)
  lls_bic[i,3] <- estimacoes_bic_limpa[[i]]$constante
}
View(lls_bic)

### Galves ----
alfabeto <- c('0','1')
tam_max <- log(nchar(sequencia_treino))
estimacoes_galves <- list()

for(i in 1:20){
  c <- c(1:20/20)[i]
  print(c)
  estimacoes_galves[[i]] <- tau_g(sequencia_treino, alfabeto, c, tam_max)
  estimacoes_galves[[i]]$constante <- c
  print(estimacoes_galves[[i]]$hat_tau)
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

lls_galves <- matrix(0, nrow = length(estimacoes_contexto_limpa), ncol = 2)

for(i in 1:length(estimacoes_contexto_limpa)){
  lls_galves[i,1] <- paste(estimacoes_galves_limpa[[i]]$hat_tau, collapse = '-')
  lls_galves[i,2] <- log_perda(sequencia_teste, estimacoes_galves_limpa[[i]]$probabilidades)
}

### BCT ----
alfabeto <- c('0','1')
tam_max <- log(nchar(sequencia_treino))
estimacoes_bct <- list()

for(i in 1:100){
  c <- c(1:100/100)[i]
  print(c)
  estimacoes_bct[[i]] <- BCT(sequencia_treino, tam_max, c)
  estimacoes_bct[[i]]$constante <- c
  print(estimacoes_bct[[i]]$Contexts)
}

estimacoes_bct_limpa <- list()
arvores <- c()
cont <- 0
for(i in 1:length(estimacoes_bct)){
  arvore_selecionada <- paste(stri_reverse(estimacoes_bct[[i]]$Contexts), collapse = '-')
  if (arvore_selecionada %in% arvores) next
  arvores <- c(arvores, arvore_selecionada)
  cont <- cont + 1
  estimacoes_bct_limpa[[cont]] <- estimacoes_bct[[i]]
  estimacoes_bct_limpa[[cont]]$constante <- estimacoes_bct[[i]]$constante
}

lls_bct <- matrix(0, nrow = length(estimacoes_bct_limpa), ncol = 3)

for(i in 1:length(estimacoes_bct_limpa)){
  lls_bct[i,1] <- paste(estimacoes_bct_limpa[[i]]$Contexts, collapse = '-')
  lls_bct[i,2] <- log_loss(sequencia, tam_max, nchar(sequencia) * prop_treino, estimacoes_bct_limpa[[i]]$constante)[ceiling(nchar(sequencia) * (1-prop_treino))]
  lls_bct[i,3] <- zero_one_loss(sequencia, tam_max, nchar(sequencia) * prop_treino, estimacoes_bct_limpa[[i]]$constante)[ceiling(nchar(sequencia) (1-prop_treino))]
}










