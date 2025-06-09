# Importa funções auxiliares ----
setwd("C:/Users/rique/Desktop/PIPGES/Dissertacao/Códigos/Comparacao_Estimadores")
source("funcoes_auxiliares.R")

## Teste qui-quadrado para homogeneidade
#\chi^2 = \sum_{a=0}^{m}\sum_{b=0}^{m} N_{bw}(a) \frac{(\hat{p}(a|bw) - \hat{p}(a|w))^2}{\hat{p}(a|w)} 
deltaH <- function(w, alfabeto, probs, counts){
  soma <- 0
  for(a in alfabeto){
    paw <- probs[w, a]
    if (is.nan(paw) || paw == 0) next
    for(b in alfabeto){
      pabw <- probs[b %s+% w, a]
      if (is.nan(pabw) || pabw == 0) next
      soma <- soma + counts[b %s+% w, a] * log(pabw / paw)
    }
  }
  return(soma)
}

delta_boot <- function(w, alfabeto, amostras, probs, cont){
  B <- length(amostras$boot)
  deltas <- rep(0, B)
  
  for(i in 1:B){
    nos_contar <- c(w, alfabeto %s+% w)
    cont_boot <- contagens(nos_contar, alfabeto, amostras$boot[[i]])
    probs_boot <- (cont_boot ) / (rowSums(cont_boot))
    deltas[i] <- deltaH(w, alfabeto, probs_boot, cont_boot)
  }
  
  obs <- deltaH(w, alfabeto, probs, cont)
  pvalue <- sum(deltas >= obs) / B 
  
  output <- list(
    mean = mean(deltas),
    deltas = deltas,
    obs = obs,
    pvalue = pvalue)
  return(output)
}

par_boot <- function(amostra, alfabeto, arvore, n, B, probs = NULL){
  amostras <- list('obs' = amostra, 'boot' = NULL, 'transicoes' = NULL)
  if(is.null(probs)){
    conts <- contagens(arvore, alfabeto, amostra)
    probs <- conts / rowSums(conts)
  }
  else {
    amostras$transicoes <- probs
    amostras$boot <- lapply(1:B, function(x) simula(probs, n))
  }
  return(amostras)
}

tau_boot <- function(amostra, limiar, tam_max, B, mensagens = F){
  alfabeto <- as.character(alfabeto)
  n <- nchar(amostra)
  
  # Gera o conjunto nos de todas as sequencias observadas até tam_max na amostra
  Vn <- calcula_Vn(amostra, tam_max)
  
  # Filtra os nós que não tem todos os irmãos, para garantir que a árvore maxial é completa
  sem_irmaos <- Filter(\(w) any(!irmaos(w, alfabeto) %in% Vn), Vn)
  nos_remover <- unlist(lapply(sem_irmaos, function(w) filhos(w, alfabeto)))
  
  nos <- setdiff(Vn, c(sem_irmaos, nos_remover))
  
  # Identifica as folhas e nós internos da árvore máximal
  folhas <- arvore_teste <- Filter(\(w) any(!filhos(w, alfabeto) %in% nos), nos)
  internos <- setdiff(nos, folhas)
  
  # Calcula as contagens dos nós
  cont <- contagens(nos, alfabeto, amostra, folhas)
  
  # Estima as probabilidades de transição
  
  probs <- cont / rowSums(cont)
  
  # Inicializa a indicadora C e o vetor de restantes
  restantes <- rev(internos)
  
  deltas <- setNames(vector("list", length(internos)), internos)
  deltas_medio <- setNames(rep(NA, length(internos)), internos)
  pvalues <- setNames(rep(NA, length(internos)), internos)
  
  
  while(length(restantes) > 0){
    w <- restantes[1]
    if(mensagens) cat("Testando nó: ", w, "\n")
    filhos_w <- filhos(w, alfabeto)
    arvore_teste <- setdiff(arvore_teste, filhos_w)
    arvore_teste <- c(arvore_teste, w)
    
    if(mensagens) cat('arvore teste:', arvore_teste)
    amostras_boot <- par_boot(amostra, alfabeto, arvore_teste, n, B, probs[arvore_teste,])
    delta <- delta_boot(w, alfabeto, amostras_boot, probs, cont)
    deltas[[w]] <- delta$deltas
    deltas_medio[w] <- delta$mean
    pvalues[w] <- delta$pvalue
  
    if(delta$pvalue > (limiar / length(arvore_teste))){
      if(mensagens) cat("Rejeita H0 para w = ", w,'p-value:', delta$pvalue, ". Realiza-se a poda, removemos seus filhos\n")
      restantes <- setdiff(restantes, w)
    }
    else{
      if(mensagens) cat("Não rejeita H0 para w = ",w,'p-value:', delta$pvalue, w, ". Mantemos o nó e seus filhos \n")
      arvore_teste <- setdiff(arvore_teste, w)
      arvore_teste <- c(arvore_teste, filhos_w)
      restantes <- setdiff(restantes, sufixos(w))
    }
    if(mensagens) cat('restantes: ', restantes, "\n")
  }
  
  return(list(
    'hat_tau' = arvore_teste,
    'probabilidades' = probs[arvore_teste,],
    'deltas' = deltas,
    'deltas_medio' = deltas_medio,
    'pvalues' = pvalues
  ))
}

transicoes <- matrix(c(0.1, 0.9,
                       0.7, 0.3,
                       0.4, 0.6,
                       0.9, 0.1), 
                     nrow = 4, byrow = TRUE, 
                     dimnames = list( c('00', '10', '01', '11'), c('0', '1')))
alfabeto <- c('0', '1')
N <- 5000
tam_max <- 5
limiar <- 0.05
B <- 1000

set.seed(999)
amostra <- simula(transicoes, N)
tau <- tau_boot(amostra, limiar, tam_max, B)
print(tau$probabilidades)
print(tau$pvalues)

set.seed(998)
amostra <- simula(transicoes, N)
tau <- tau_boot(amostra, limiar, tam_max, B)
print(tau$probabilidades)
print(tau$pvalues)

set.seed(997)
amostra <- simula(transicoes, N)
tau <- tau_boot(amostra, limiar, tam_max, B)
print(tau$probabilidades)
print(tau$pvalues)

set.seed(996)
amostra <- simula(transicoes, N)
tau <- tau_boot(amostra, limiar, tam_max, B)
print(tau$probabilidades)
print(tau$pvalues)

set.seed(995)
amostra <- simula(transicoes, N)
tau <- tau_boot(amostra, limiar, tam_max, B)
print(tau$probabilidades)
print(tau$pvalues)

set.seed(994)
amostra <- simula(transicoes, N)
tau <- tau_boot(amostra, limiar, tam_max, B)
print(tau$probabilidades)
print(tau$pvalues)


          