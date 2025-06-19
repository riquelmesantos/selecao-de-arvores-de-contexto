# Importa funções auxiliares ----
setwd("C:/Users/rique/Desktop/PIPGES/Dissertacao/Códigos/Comparacao_Estimadores")
source("funcoes_auxiliares.R")

# Medida de discrepância ----
# \Delta_G(by) = \max_{a\in \ca A}\{|\hat p' (a\mid by) - \hat p'(a\mid y)|\}

deltaG <- function(bw, probs) {
  w <- pai(bw)
  if (w == '') paw = probs['lambda', ]
  else paw = probs[w, ]
  pabw = probs[bw, ]
  delta <- max(abs(paw - pabw))
  return(delta)
}

tau_g <- function(amostra, alfabeto, limiar, tam_max) {
  alfabeto <- as.character(alfabeto)
  
  # Gera o conjunto nos de todas as palavras observadas até tam_max na amostra
  nos <- calcula_Vn(amostra, tam_max+1)
  
  # Identifica as folhas e nós internos da árvore máxima
  folhas <- Filter(\(w) any(!filhos(w, alfabeto) %in% nos), nos)
  internos <- unique(unlist(sapply(folhas, function(w) sufixos(w, l_min = 2))))
  internos <- internos[order(nchar(internos))]
  nos_admissiveis <-  c(folhas, internos)
  nos_admissiveis <- nos_admissiveis[order(nchar(nos_admissiveis))]
  
  cont <- contagens(nos_admissiveis, alfabeto, amostra, folhas, lambda = T)
  somas <- rowSums(cont)
  # Estima as probabilidades de transição suavizadas (Laplace)
  probs <- (cont + 1) / (somas + length(alfabeto))
  
  # Inicializa a função indicadora para cada nó
  C <- setNames(rep(0, length(nos)), nos)
  # Aponta discrepâncias nas folhas
  deltas <- sapply(nos_admissiveis, function(w) deltaG(w, probs))
  
  for(w in rev(nos_admissiveis)){
    if(C[w] == 1) next
    filhos_w <- filhos(w, alfabeto)
    irmaos_w <- intersect(irmaos(w, alfabeto), nos_admissiveis)
    if(all(filhos_w %in% nos_admissiveis)){
      if(all(C[filhos_w] == 1)) next
      else if(any(deltas[irmaos_w] > limiar)){
        C[w] <- 1
        next
      }
    }
    if(any(deltas[irmaos_w] > limiar)) C[w] <- 1
  }
  arvore <- names(C[C==1])
  
  # Retorna a árvore selecionada e as probabilidades suavizadas associadas
  return(list('hat_tau' = arvore, 
              'probabilidades_arvore' = probs[c('lambda',arvore),],
              'probabilidades' = probs,
              'deltas' = deltas,
              'nos' = nos))
}

tau_g_calculado <- function(deltas, nos, probs, alfabeto, limiar) {
  nos_admissiveis <- names(deltas)
  C <- setNames(rep(0, length(nos)), nos)
  
  for(w in rev(nos_admissiveis)){
    if(C[w] == 1) next
    filhos_w <- filhos(w, alfabeto)
    irmaos_w <- intersect(irmaos(w, alfabeto), nos_admissiveis)
    if(all(filhos_w %in% nos_admissiveis)){
      if(all(C[filhos_w] == 1)) next
      else if(any(deltas[irmaos_w] > limiar)){
        C[w] <- 1
        next
      }
    }
    if(any(deltas[irmaos_w] > limiar)) C[w] <- 1
  }
  arvore <- names(C[C==1])
  
  # Retorna a árvore selecionada e as probabilidades suavizadas associadas
  return(list('hat_tau' = arvore, 
              'probabilidades_arvore' = probs[c('lambda',arvore),],
              'probabilidades' = probs,
              'deltas' = deltas,
              'nos' = nos))
}

# # 
# # # ## Testes ----
# # # ### Teste 1: Cadeia de Markov com ordem 2 em alfabeto binário ----
# transicoes <- matrix(c(0.1, 0.9,
#                        0.7, 0.3,
#                        0.4, 0.6,
#                        0.9, 0.1),
#                      nrow = 4, byrow = TRUE,
#                      dimnames = list( c('00', '10', '01', '11'), c('0', '1')))
# 
# alfabeto <- c('0', '1')
# limiar <- 0.2
# tam_max <- 5
# N <- 8000
# 
# set.seed(999)
# amostra <- simula(transicoes, N)
# arvore_selecionada <- tau_g(amostra, alfabeto, limiar, tam_max)
# print(arvore_selecionada$hat_tau)   # Estimou corretamente!
# print(arvore_selecionada$probabilidades)
# 
# set.seed(998)
# amostra <- simula(transicoes, N)
# arvore_selecionada <- tau_g(amostra, alfabeto, limiar, tam_max)
# print(arvore_selecionada$hat_tau)   # Estimou corretamente!
# print(arvore_selecionada$probabilidades)
# 
# set.seed(997)
# amostra <- simula(transicoes, N)
# arvore_selecionada <- tau_g(amostra, alfabeto, limiar, tam_max)
# print(arvore_selecionada$hat_tau)   ## Estimou corretamente!
# print(arvore_selecionada$probabilidades)
# 
# ### Teste 2: Cadeia com memória de comprimento variável \tau = {1, 10, 100, 1000, 0000} ----
# transicoes <- matrix(c(0.6, 0.4,
#                        0.4, 0.6,
#                        0.3, 0.7,
#                        0.2, 0.8,
#                        0.1, 0.9),
#                      nrow = 5, byrow = TRUE, 
#                      dimnames = list( c('1', '10', '100', '1000', '0000'), c('0', '1')))
# alfabeto <- c('0', '1')
# limiar <- 0.3
# N <- 100000
# tam_max <- 8
# 
# set.seed(999)
# amostra <- simula(transicoes, N)
# arvore_selecionada <- tau_g(amostra, alfabeto, limiar, 8)   
# print(arvore_selecionada$hat_tau)   # Estimou errado, mas manteve umma estrutura semelhante
# print(arvore_selecionada$probabilidades)
# 
# set.seed(998)
# amostra <- simula(transicoes, N)
# arvore_selecionada <- tau_g(amostra, alfabeto, limiar, 8)   
# print(arvore_selecionada$hat_tau)   # Não conseguiu estimar
# print(arvore_selecionada$probabilidades)
# 
# set.seed(997)
# amostra <- simula(transicoes, N)
# arvore_selecionada <- tau_g(amostra, alfabeto, limiar, 8)   
# print(arvore_selecionada$hat_tau)   # Estimou errado, mas manteve umma estrutura semelhante
# print(arvore_selecionada$probabilidades)
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
# limiar <- 0.5
# N <- 100000
# 
# amostra <- simula(transicoes, N)
# arvore_selecionada <- tau_g(amostra, alfabeto, limiar, 8)   
# print(arvore_selecionada$hat_tau)   # Estimou errado...
# print(arvore_selecionada$probabilidades)
# 
# set.seed(998)
# amostra <- simula(transicoes, N)
# arvore_selecionada <- tau_g(amostra, alfabeto, limiar, 8)   
# print(arvore_selecionada$hat_folhas)   # Estimou errado...
# print(arvore_selecionada$probabilidades)
# 
# set.seed(997)
# amostra <- simula(transicoes, N)
# arvore_selecionada <- tau_g(amostra, alfabeto, limiar, 8)   
# print(arvore_selecionada$hat_folhas)   # Estimou errado
# print(arvore_selecionada$probabilidades)


