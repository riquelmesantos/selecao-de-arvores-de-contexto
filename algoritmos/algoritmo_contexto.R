# Importa funções auxiliares ----
setwd("C:/Users/rique/Desktop/PIPGES/Dissertacao/Códigos/Comparacao_Estimadores")
source("funcoes_auxiliares.R")

# Medida de discrepância ----
# \Delta_R(w) = \sum_{b \in \mathcal{A}} \sum_{a \in \mathcal{A}} N(bwa) \log\left[\frac{\hat{p}(a \mid bw)}{\hat p(a \mid w)}\right]. 

deltaR <- function(w, alfabeto, probs, counts){
  soma <- 0
  for (a in alfabeto) {
    paw <- probs[w,a]
    if (is.nan(paw) || paw == 0) next
    
    for (b in alfabeto) {
      N_bwa = counts[b %s+% w, a]
      pabw <- probs[b %s+%  w, a]
      
      if (is.nan(pabw) || pabw == 0) next
      soma <- soma + N_bwa * log(pabw/paw)
    }
  }
  return(soma)
}

# Algoritmo Contexto ---- 
tau_c <- function(amostra, alfabeto, limiar, tam_max) {
  alfabeto <- as.character(alfabeto)
  n <- nchar(amostra)
  
  # Gera o conjunto nos de todas as palavras observadas até tam_max na amostra
  nos <- calcula_Vn(amostra, tam_max)
  
  # Identifica as folhas e nós internos da árvore máxima
  folhas <- Filter(\(w) any(!filhos(w, alfabeto) %in% nos), nos)
  internos <- setdiff(nos, folhas)
  
  # Calcula as contagens dos nós
  cont <- contagens(nos, alfabeto, amostra, folhas)
  
  # Estima as probabilidades de transição
  probs <- cont / rowSums(cont)
  
  # Inicializa a matriz indicadora C
  C <- setNames(rep(0, length(nos)), nos)
  deltas <- setNames(rep(NA, length(nos)), nos)
  
  # Avalia, de forma bottom-up, quais contextos devem ser mantidos com base na medida deltaR
  for (w in rev(internos)) {
    delta <- deltaR(w, alfabeto, probs, cont)
    deltas[w] <- delta
    filhos <- paste0(alfabeto, w)
    if (any(C[filhos] == 1)) { # se algum filho for relevante, o nó pai também é
      C[w] <- 1  
      next
    }
    if (delta > limiar) C[w] <- 1  
  }
  
  # Constrói a árvore final a partir da função indicadora C
  arvore <- c()
  for (w in nos) {
    pai_w <- pai(w)
    if ((pai_w == "lambda" && C[w] == 0) || (C[w] == 0 && C[pai_w] == 1)) {
      arvore <- c(arvore, w)
    }
  }
  
  # Retorna as folhas da árvore e as probabilidades associadas
  return(list(
    'hat_tau' = arvore,
    'probabilidades_arvore' = probs[arvore,],
    'deltas' = deltas,
    'probabilidades' = probs,
    'internos' = internos
  ))
}

tau_c_calculado <- function(probs, deltas, internos, limiar){
  nos <- names(deltas)
  C <- setNames(rep(0, length(nos)), nos)
  
  for (w in rev(internos)) {
    delta <-  deltas[w] 
    filhos <- paste0(alfabeto, w)
    if (any(C[filhos] == 1)) { 
      C[w] <- 1  
      next
    }
    if (delta > limiar) C[w] <- 1  
  }
  
  arvore <- c()
  for (w in nos) {
    pai_w <- pai(w)
    if ((pai_w == "lambda" && C[w] == 0) || (C[w] == 0 && C[pai_w] == 1)) {
      arvore <- c(arvore, w)
    }
  }
  return(list(
    'hat_tau' = arvore,
    'probabilidades_arvore' = probs[arvore,],
    'deltas' = deltas,
    'probabilidades' = probs,
    'internos' = internos
  ))
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
# N <- 500
# limiar <- log(N)
# 
# set.seed(999)
# amostra <- simula(transicoes, N)
# arvore_selecionada <- tau_c(amostra, alfabeto, limiar, 10)   
# print(arvore_selecionada)   # Estimou corretamente!
# 
# set.seed(998)
# amostra <- simula(transicoes, N)
# arvore_selecionada <- tau_c(amostra, alfabeto, limiar, 10)   
# print(arvore_selecionada)   # Estimou corretamente!
# 
# set.seed(997)
# amostra <- simula(transicoes, N)
# arvore_selecionada <- tau_c(amostra, alfabeto, limiar, 10)   
# print(arvore_selecionada)   # Estimou corretamente!
# 
# ### Teste 2: Cadeia com memória de comprimento variável \folhas = {1, 10, 100, 1000, 0000} ----
# transicoes <- matrix(c(0.6, 0.4,
#                        0.4, 0.6,
#                        0.3, 0.7,
#                        0.2, 0.8,
#                        0.1, 0.9),
#                      nrow = 5, byrow = TRUE, 
#                      dimnames = list( c('1', '10', '100', '1000', '0000'), c('0', '1')))
# alfabeto <- c('0', '1')
# N <- 30000
# limiar <- log(N)
# 
# set.seed(999)
# amostra <- simula(transicoes, N)
# arvore_selecionada <- tau_c(amostra, alfabeto, limiar, 10)   
# print(arvore_selecionada)   # Estimou corretamente!
# 
# set.seed(998)
# amostra <- simula(transicoes, N)
# arvore_selecionada <- tau_c(amostra, alfabeto, limiar, 10)   
# print(arvore_selecionada)   # Estimou corretamente!
# 
# set.seed(997)
# amostra <- simula(transicoes, N)
# arvore_selecionada <- tau_c(amostra, alfabeto, limiar, 10)   
# print(arvore_selecionada)   # Errou por um ramo a menos, faltou crescer o 000...
# 
# ### Teste 3: Cadeia com memória de comprimento variável \folhas = {00, 10, 20, 11, 21, 02, 12, 22, 001, 101, 201} ----
# transicoes <- matrix(c(0.1, 0.3, 0.6,
#                        0.2, 0.5, 0.3,
#                        0.7, 0.2, 0.1,
#                        0.4, 0.4, 0.2,
#                        0.8, 0.1, 0.1,
#                        0.3, 0.4, 0.3,
#                        0.1, 0.2, 0.7,
#                        0.3, 0.6, 0.1,
#                        0.1, 0.2, 0.7,
#                        0.4, 0.5, 0.1,
#                        0.5, 0.3, 0.2),
#                      nrow = 11, byrow = TRUE, 
#                      dimnames = list(
#                        c('00', '10', '20', '11', '21', '02', '12', '22', '001', '101', '201'), 
#                        c('0', '1', '2')))
# 
# alfabeto <- c('0', '1', '2')
# N <- 9000
# limiar <- log(N)
# 
# set.seed(999)
# amostra <- simula(transicoes, N)
# arvore_selecionada <- tau_c(amostra, alfabeto, limiar, 10)   
# print(arvore_selecionada)   # Estimou corretamente!
# 
# set.seed(998)
# amostra <- simula(transicoes, N)
# arvore_selecionada <- tau_c(amostra, alfabeto, limiar, 10)   
# print(arvore_selecionada)   # Estimou corretamente!
# 
# set.seed(997)
# amostra <- simula(transicoes, N)
# arvore_selecionada <- tau_c(amostra, alfabeto, limiar, 10)   
# print(arvore_selecionada)  # Errou por alguns ramos a mais, cresceu o 02, 102, 2102, 02102, 102102, 1102102 provavelmente porque não apareceram essas sequencias o suficiente na amostra
# contagens(arvore_selecionada$hat_tau, alfabeto, amostra)
# 
# 


