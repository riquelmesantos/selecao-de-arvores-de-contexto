setwd('C:/Users/rique/Desktop/PIPGES/Dissertacao/Códigos/Comparacao_Estimadores')
source('funcoes_auxiliares.R')
selecoes_contexto <- readRDS("analise_dados_reais/dados/selecoes_contexto_neuro_30.rds")
selecoes_bic <- readRDS("analise_dados_reais/dados/selecoes_bic_neuro_30.rds")
selecoes_galves <- readRDS("analise_dados_reais/dados/selecoes_galves_neuro_30.rds")

log_like <- function(cont, probs, alfabeto){
  contextos <- rownames(probs)
  ll = 0
  for(w in contextos){
    for(a in alfabeto){
      if (is.na(probs[w, a]) && !is.na(probs['lambda', a])) ll <- ll + (cont[w, a] * probs['lambda', a])
      if (is.na(probs[w, a]) || probs[w, a] == 0) next
      ll <- ll + (cont[w, a] * log(probs[w, a]))
    }
  }
  return(ll)
}

aic <- function(ll, k) -2*ll + 2*k

bic <- function(ll, k, seq) -2*ll + log(length(seq))*k

residuo_pearson <- function(seq, probs) {
  residuos <- numeric(length(seq))
  tam_max <- max(nchar(rownames(probs)))
  
  for (i in (tam_max+1):length(seq)) {
    passado <- seq[(i-tam_max):(i-1)] 
    contexto <- intersect(sufixos(paste(passado, collapse = '')), rownames(probs))
    
    if (length(contexto) == 0) stop('Não foi encontrado contexto para a sequência: ', paste(passado, collapse=''))
    
    contexto <- contexto[nchar(contexto) == max(nchar(contexto))]
    prox_simbolo <- as.character(seq[i])
    p <-  probs[contexto, prox_simbolo]
    p <- min(max(p, 1e-10), 1 - 1e-10)
    
    if(seq[i] == '1')  residuos[i] <- (1 - p) / sqrt(p * (1 - p))
    else if(seq[i] == '0') residuos[i] <- (0 - p) / sqrt(p * (1 - p))
    else stop("Símbolo desconhecido: ", seq[i])
  }
  return(residuos)
}

residuo_person_deviance <- function(seq, probs) {
  residuo_pearson_ <- numeric(length(seq))
  residuo_deviance_ <- numeric(length(seq))
  tam_max <- max(nchar(rownames(probs)))
  
  for (i in (tam_max+1):length(seq)) {
    passado <- paste(seq[(i-tam_max):(i-1)] , collapse = '')
    contexto <- intersect(sufixos(passado), rownames(probs))
    
    if (length(contexto) == 0) stop('Não foi encontrado contexto para a sequência: ', paste(passado, collapse=''))
    
    contexto <- contexto[nchar(contexto) == max(nchar(contexto))]
    prox_simbolo <- as.character(seq[i])
    p <-  probs[contexto, prox_simbolo]
    p <- min(max(p, 1e-10), 1 - 1e-10)
    
    if(seq[i] == 1){
      residuo_pearson_[i] <- (1 - p) / sqrt(p * (1 - p))
      residuo_deviance_[i] <- sign(1 - p) * sqrt(-2 * log(p))
    }  
    else if(seq[i] == 0){ 
      residuo_pearson_[i] <- (0 - p) / sqrt(p * (1 - p))
      residuo_deviance_[i] <- sign(0 - p) * sqrt(-2 * log(1 - p))
      }
    else stop("Símbolo desconhecido: ", seq[i])
  }
  return(list(residuos_pearson = residuo_pearson_, 
              residuos_deviance = residuo_deviance_))
}


nomes <- rep(c('contexto', 'bic', 'galves'), 
               c(length(selecoes_contexto), length(selecoes_bic), length(selecoes_galves)))
modelos_selecionados <- c(selecoes_contexto, selecoes_bic, selecoes_galves)
for(i in 1:length(modelos_selecionados)) modelos_selecionados[[i]]$metodo <- nomes[i]
arvores_contexto <- sapply(1:length(selecoes_contexto), function(i) paste(sort(selecoes_contexto[[i]]$hat_tau), collapse = '-'))
arvores_bic <- sapply(1:length(selecoes_bic), function(i) paste(sort(selecoes_bic[[i]]$hat_tau), collapse = '-'))
arvores_galves <- sapply(1:length(selecoes_galves), function(i) paste(sort(selecoes_galves[[i]]$hat_tau), collapse = '-'))

modelos <- list()
arvores <- c()
cont <- 0

for(i in 1:length(modelos_selecionados)){
  arvore_selecionada <- paste(sort(modelos_selecionados[[i]]$hat_tau), collapse = ' ')
  if (arvore_selecionada %in% arvores) next
  arvores <- c(arvores, arvore_selecionada)
  cont <- cont + 1
  modelos[[cont]] <- modelos_selecionados[[i]]
}

conts_geral <- matrix(0, ncol = 2, dimnames = list(NULL, c(0,1)))
for(i in 1:length(modelos)){
  print(i)
  probs <- modelos[[i]]$probabilidades_arvore
  arvore <- modelos[[i]]$hat_tau
  conts <- contagens(setdiff(arvore, rownames(conts_geral)), alfabeto, sequencia)
  conts_geral <- rbind(conts_geral, conts)
  k <- length(arvore)
  modelos[[i]]$logLik <- log_like(conts_geral, probs, alfabeto)
  modelos[[i]]$aic <- aic(modelos[[i]]$logLik, k)
  modelos[[i]]$bic <- bic(modelos[[i]]$logLik, k, sequencia_vet)
  modelos[[i]]$residuo <- residuo_person_deviance(sequencia_vet, probs)
}

aics <- sapply(1:17, function(i) modelos[[i]]$aic)
bics <- sapply(1:17, function(i) modelos[[i]]$bic)
lls <- sapply(1:17, function(i) modelos[[i]]$logLik)

plot(aics, type = 'b', col = 'blue', pch = 19, 
     xlab = 'Modelos', ylab = 'AIC', main = 'Comparação de Modelos (AIC)')
plot(bics, type = 'b', col = 'red', pch = 19, 
     xlab = 'Modelos', ylab = 'BIC', main = 'Comparação de Modelos (BIC)')
plot(lls, type = 'b', col = 'green', pch = 19, 
     xlab = 'Modelos', ylab = 'Log-Likelihood', main = 'Comparação de Modelos (Log-Likelihood)')
order(aics)
order(bics)
order(lls, decreasing = T)

residuos_pearson_ <- lapply(1:17, function(i) modelos[[i]]$residuo$residuos_pearson)
residuos_pearson_0 <- lapply(1:17, function(i) modelos[[i]]$residuo$residuos_pearson[sequencia_vet == 0])
residuos_pearson_1 <- lapply(1:17, function(i) modelos[[i]]$residuo$residuos_pearson[sequencia_vet == 1])
residuos_deviance_ <- lapply(1:17, function(i) modelos[[i]]$residuo$residuos_deviance)
residuos_deviance_0 <- lapply(1:17, function(i) modelos[[i]]$residuo$residuos_deviance[sequencia_vet == 0])
residuos_deviance_1 <- lapply(1:17, function(i) modelos[[i]]$residuo$residuos_deviance[sequencia_vet == 1])
order(sapply(1:17, function(i) mean(abs(residuos_deviance_0[[i]]))) +  0.6*sapply(1:17, function(i) mean(abs(residuos_deviance_1[[i]]))))
order(sapply(1:17, function(i) mean(residuos_pearson_0[[i]]^2)) + sapply(1:17, function(i) mean(residuos_pearson_1[[i]]^2)))
order(sapply(1:17, function(i) mean(residuos_pearson_0[[i]]^2)))

modelos[[1]]$hat_tau


hist(residuos_pearson_[[9]])
hist(residuos_pearson_[[5]])
plot(residuos_pearson_[[5]][1:10000])

plot_arvore(modelos[[15]]$hat_tau, T)
