# Importa funções auxiliares ----
setwd("C:/Users/rique/Desktop/PIPGES/Dissertacao/Códigos/Comparacao_Estimadores")
source("funcoes_auxiliares.R")

# Probabilidades ponderadas
log_P_w <- function(w, alfabeto, probs, cont) {
  probs <- probs 
  fatores <- sapply(alfabeto, function(a)
    cont[w, a] * log(probs[w, a]))
  soma <- sum(fatores)
  return(soma)
}

# Função de penalização
pen <- function(amostra, alfabeto, c = 1) {
  a <- length(alfabeto)
  n <- nchar(amostra)
  return(c * (a - 1) * log(n))
}

# Seleção via penalização da verossimilhança
tau_bic <- function(amostra, alfabeto, tam_max, c = 1){
  arvore <- c()
  nos <- calcula_Vn(amostra, tam_max)
  cont <- contagens(nos, alfabeto, amostra)
  contagem_zero <- rownames(cont[rowSums(cont == 0) > 0, ])
  nos <- c(setdiff(nos, contagem_zero))
  folhas <- Filter(\(w) any(!filhos(w, alfabeto) %in% nos), nos)
  internos <- setdiff(nos, folhas)
  
  probs <- cont / rowSums(cont)
  penalizacao <- pen(amostra, alfabeto, c)
  
  log_Vw <- setNames(rep(-Inf, length(nos)), nos)
  Xw <- setNames(rep(0, length(nos)), nos)
  log_P_ws <- sapply(nos, function(w) log_P_w(w, alfabeto, probs, cont))
  log_Vw[folhas] <- sapply(folhas, function(w) log_P_ws[w] - penalizacao)
  
  for(w in rev(internos)){
    filhos_w <- intersect(filhos(w, alfabeto), nos)
    log_Vw_w <-  log_P_ws[w]- penalizacao
    log_Vw[w] <- log_Vw_w
    soma_logs_filhos <- sum(log_Vw[filhos_w])
    
    if(soma_logs_filhos > log_Vw_w){
      log_Vw[w] <- soma_logs_filhos
      Xw[w] <- 1
    }
  }
  
  for(w in rev(nos)){
    nw <- nchar(w)
    pais <- setdiff(sufixos(w), w)
    if((nw == 1 && Xw[w] == 0) || (Xw[w] == 0 && all(Xw[pais] == 1)))
      arvore <- c(arvore, w)
  }
  return(list('hat_tau' = arvore, 
              'probabilidades_arvore' = probs[arvore,],
              'log_P_ws' = log_P_ws,
              'probabilidades' = probs,
              'folhas' = folhas))
}

tau_bic_calculado <- function(amostra, alfabeto, log_P_ws, probs, folhas, c = 1){
  nos <- names(log_P_ws)
  internos <- setdiff(nos, folhas)
  penalizacao <- pen(amostra, alfabeto, c)
  Xw <- setNames(rep(0, length(nos)), nos)
  log_Vw <- setNames(rep(-Inf, length(nos)), nos)
  log_Vw[folhas] <- sapply(folhas, function(w) log_P_ws[w] - penalizacao)
  
  for(w in rev(internos)){
    filhos_w <- intersect(filhos(w, alfabeto), nos)
    log_Vw_w <-  log_P_ws[w]- penalizacao
    log_Vw[w] <- log_Vw_w
    soma_logs_filhos <- sum(log_Vw[filhos_w])
    if(soma_logs_filhos > log_Vw_w){
      log_Vw[w] <- soma_logs_filhos
      Xw[w] <- 1
    }
  }
  arvore <- c()
  for(w in rev(nos)){
    nw <- nchar(w)
    pais <- setdiff(sufixos(w), w)
    if((nw == 1 && Xw[w] == 0) || (Xw[w] == 0 && all(Xw[pais] == 1)))
      arvore <- c(arvore, w)
  }
  return(list('hat_tau' = arvore, 
              'probabilidades_arvore' = probs[arvore,],
              'log_P_ws' = log_P_ws,
              'probabilidades' = probs,
              'internos' = internos))
}



# ## Testes ----
# ### Teste 1: Cadeia de Markov com ordem 2 em alfabeto binário ----
# transicoes <- matrix(c(0.1, 0.9,
#                        0.7, 0.3,
#                        0.4, 0.6,
#                        0.9, 0.1),
#                      nrow = 4, byrow = TRUE,
#                      dimnames = list( c('00', '10', '01', '11'), c('0', '1')))
# alfabeto <- c('0', '1')
# tam_max <- 10
# N <- 500
# 
# set.seed(999)
# amostra <- simula(transicoes, N)
# arvore_selecionada <- tau_bic(amostra, alfabeto, tam_max)   
# print(arvore_selecionada$hat_tau)   # Estimou corretamente!
# print(arvore_selecionada$probabilidades)
# 
# set.seed(998)
# amostra <- simula(transicoes, N)
# arvore_selecionada <- tau_bic(amostra, alfabeto, tam_max)   
# print(arvore_selecionada$hat_tau)   # Estimou corretamente!
# print(arvore_selecionada$probabilidades)
# 
# set.seed(997)
# amostra <- simula(transicoes, N)
# arvore_selecionada <- tau_bic(amostra, alfabeto, tam_max)   
# print(arvore_selecionada$hat_tau)   # Errou pois podou um nó, talvez a penalização tenha sido excessiva
# print(arvore_selecionada$probabilidades)








