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

log_like <- function(cont, probs, alfabeto){
  contextos <- rownames(probs)
  ll = 0
  for(w in contextos){
    for(a in alfabeto){
      if (is.na(probs[w, a]) || probs[w, a] == 0) next
      ll <- ll + (cont[w, a] * log(probs[w, a]))
    }
  }
  return(ll)
}

residuos_deviance <- function(seq, probs) {
  residuos <- numeric(length(seq))
  tam_max <- max(nchar(rownames(probs)))
  
  for (i in (tam_max+1):length(seq)) {
    passado <- seq[(i-tam_max):(i-1)] 
    contexto <- intersect(sufixos(paste(passado, collapse = '')), rownames(probs))
    
    if (length(contexto) == 0) return('Não foi encontrado contexto para a sequência', passado)
    
    contexto <- contexto[nchar(contexto) == max(nchar(contexto))]
    prox_simbolo <- seq[i]
    p <-  probs[contexto, prox_simbolo]
    p <- min(max(p, 1e-10), 1 - 1e-10)
    
    if(seq[i] == 1)  residuos[i] <- sign(1 - p) * sqrt(-2 * log(p))
    else if(seq[i] == 0) residuos[i] <- sign(0 - p) * sqrt(-2 * log(1 - p))
    else stop("Símbolo desconhecido: ", seq[i])
  }
  return(residuos)
}

residuo_pearson <- function(seq, probs) {
  residuos <- numeric(length(seq))
  tam_max <- max(nchar(rownames(probs)))
  
  for (i in (tam_max+1):length(seq)) {
    passado <- seq[(i-tam_max):(i-1)] 
    contexto <- intersect(sufixos(paste(passado, collapse = '')), rownames(probs))
    
    if (length(contexto) == 0) return('Não foi encontrado contexto para a sequência', passado)
    
    contexto <- contexto[nchar(contexto) == max(nchar(contexto))]
    prox_simbolo <- seq[i]
    p <-  probs[contexto, prox_simbolo]
    p <- min(max(p, 1e-10), 1 - 1e-10)
    
    if(seq[i] == '1')  residuos[i] <- (1 - p) / sqrt(p * (1 - p))
    else if(seq[i] == '0') residuos[i] <- (0 - p) / sqrt(p * (1 - p))
    else stop("Símbolo desconhecido: ", seq[i])
  }
  return(residuos)
}

sequencia_binaria <- function(cadeia, bin) {
  breaks <- seq(0, floor(max(cadeia)) + 1, by = bin)
  categorias <- cut(cadeia[,1], breaks = breaks, right = TRUE, labels = FALSE, include.lowest = TRUE)
  binario <- rep(0, length(breaks) - 1)
  binario[unique(na.omit(categorias))] <- 1
  return(binario)
}

block_bootstrap <- function(sequencia, block_length, n, B) {
  sequencia <- strsplit(sequencia, '')[[1]]
  n_blocks <- ceiling(n / block_length)
  blocks <- matrix(0, nrow = n_blocks, ncol = block_length)
  samples <- list()
  
  for(i in 1:B){
    for (j in 1:n_blocks) {
      start <- sample(1:(length(sequencia) - block_length + 1), 1)
      blocks[j, ] <- sequencia[start:(start + block_length - 1)]
    }
   samples[[i]] <- paste(as.vector(t(blocks)), collapse = '')
  }
  return(samples)
}


dados_neuronio <- read.table("Conjunto de dados/cell_30.dat")
sequencia_vet <- sequencia_binaria(dados_neuronio, 0.029)
sequencia <- paste0(sequencia_vet, collapse = '')
log(nchar(sequencia))

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
# saveRDS(estimacoes_contexto_limpa, file = "estimacoes_contexto_limpa.rds")

lls_contexto <- matrix(0, nrow = length(estimacoes_contexto_limpa), ncol = 4)


residuos_deviance_contexto <- list()
for(i in 1:length(estimacoes_contexto_limpa)){
  print(i)
  lls_contexto[,1][i] <- paste(estimacoes_contexto_limpa[[i]]$hat_tau, collapse = '-')
  probs <- estimacoes_contexto_limpa[[i]]$probabilidades_arvore
  lls_contexto[i,2] <- log_perda(sequencia, probs)
  lls_contexto[i,3] <- estimacoes_contexto_limpa[[i]]$constante
}


View(lls_contexto)

### BIC ----
alfabeto <- c('0','1')
tam_max <- log(nchar(sequencia))
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
saveRDS(estimacoes_bic_limpa, file = "estimacoes_bic_limpa.rds")

lls_bic <- matrix(0, nrow = length(estimacoes_bic_limpa), ncol = 3)

for(i in 1:length(estimacoes_bic_limpa)){
  lls_bic[i,1] <- paste(estimacoes_bic_limpa[[i]]$hat_tau, collapse = '-')
  probs <- estimacoes_bic_limpa[[i]]$probabilidades_arvore
  lls_bic[i,2] <- log_perda(sequencia, probs)
  lls_bic[i,3] <- estimacoes_bic_limpa[[i]]$constante
}
View(lls_bic)

residuos_deviance_bic <- list()
residuos_pearson_bic <- list()
for(i in 1:length(estimacasoes_bic_limpa)){
  print(i)
  residuos_deviance_bic[[i]] <- residuos_deviance(sequencia_vet, estimacoes_bic_limpa[[i]]$probabilidades_arvore)
  hist(abs(residuos_deviance_bic[[i]]), main = paste("Resíduos de deviance - BIC -", paste(estimacoes_bic_limpa[[i]]$hat_tau, collapse = '-')), xlab = "Resíduos", ylab = "Frequência")
  residuos_pearson_bic[[i]] <- residuo_pearson(sequencia_vet, estimacoes_bic_limpa[[i]]$probabilidades_arvore)
  hist(abs(residuos_pearson_bic[[i]]), main = paste("Resíduos de Pearson - BIC -", paste(estimacoes_bic_limpa[[i]]$hat_tau, collapse = '-')), xlab = "Resíduos", ylab = "Frequência")
}

sbdr <- sapply(1:12, function(i) mean(residuos_deviance_bic[[i]][residuos_deviance_bic[[i]] < 0]))
plot(residuos_deviance_bic[[1]][1:1000])
plot(residuos_pearson_bic[[1]][1:1000])


### Galves ----
alfabeto <- c('0','1')
tam_max <- 8
estimacoes_galves <- list()
c <- 0.05
estimacoes_galves[[1]] <- tau_g(sequencia_treino, alfabeto, c, tam_max)
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
# saveRDS(estimacoes_galves_limpa, file = "estimacoes_galves_limpa.rds")
lls_galves <- matrix(0, nrow = length(estimacoes_galves_limpa), ncol = 2)

for(i in 1:length(estimacoes_galves_limpa)){
  lls_galves[i,1] <- paste(estimacoes_galves_limpa[[i]]$hat_tau, collapse = '-')
  lls_galves[i,2] <- log_perda(sequencia_teste, estimacoes_galves_limpa[[i]]$probabilidades)
}
View(lls_galves)




