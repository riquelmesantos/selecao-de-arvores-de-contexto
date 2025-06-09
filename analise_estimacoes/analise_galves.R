#  Importa funções auxiliares e pacotes necessários ----
wd <- "C:/Users/rique/Desktop/PIPGES/Dissertacao/Códigos/Comparacao_Estimadores"
setwd(wd)
source("funcoes_auxiliares.R")
library(dplyr, quietly = T)
library(tidyverse, quietly = T)

# Importa os dados da seleções ----
# tau_1_galves <- readRDS("estimacoes/result_galves_tau_1.RDS")
# tau_2_galves <- readRDS("estimacoes/result_galves_tau_2.RDS")
 tau_3_galves <- readRDS("estimacoes/result_galves_tau_3_maiores.RDS")
# tau_4_galves <- readRDS("estimacoes/result_galves_tau_4.RDS")

# Função para selecionar a árvore utilizando outros limiares ----
tau_g_deltas <- function(deltas, limiar, alfabeto, folhas_internos, tau, tam_max = NULL) {
  if(is.null(tam_max)) tam_max <- Inf
  deltas <- deltas[nchar(names(deltas)) <= tam_max]
  deltas <-  deltas[!is.na(deltas)]
  nos <- names(deltas)
  internos <- c(unique(sapply(tau, function(w) pai(w))), tau)
  externos <- setdiff(nos, c(internos, folhas_internos, tau))
  
  seleciona <- sapply(folhas_internos, function(w) any(deltas[filhos(w, alfabeto)] > limiar))
  nao_seleciona <- sapply(externos, function(w) all(deltas[w] <= limiar))
  if(all(seleciona) && all(nao_seleciona[!is.na(nao_seleciona)])) return(1)
  return(0)
}

tau_g_deltas_ <- function(deltas, alfabeto, limiar) {
  alfabeto <- as.character(alfabeto)
  nos <- names(deltas)
  
  C <- setNames(rep(0, length(nos)), nos)
  for(w in rev(nos)){
    if(C[w] == 1) next
    filhos_w <- filhos(w, alfabeto)
    irmaos_w <- intersect(irmaos(w, alfabeto), nos)
    if(all(filhos_w %in% nos)){
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
  return(arvore)
}


# Árvores modelo---
tau_1 <- c('00', '10', '01', '11')
tau_2 <- c('1', '10', '100', '1000', '0000')
tau_3 <- c('2', '00', '10', '20', '11', '21', '001', '101', '201')
tau_4 <- c('0', '1', '02', '12', '22', '03', '13', '33', '032', '232', '332', '023', '123', '223', '323', '0132', '1132', '2132', '3132')

# Análise das árvores estimadas para diferentes constantes ----
limiares <- 1:40/40

modelos <- list(
  #list(nome = "tau_1", result_selecoes = tau_1_galves, folhas = c('0', '1'), alfabeto = c('0', '1'), tau = tau_1)#,
  # list(nome = "tau_2", result_selecoes = tau_2_galves, folhas = c('0', '00', '000'), alfabeto = c('0', '1'), tau = tau_2),
  list(nome = "tau_3", result_selecoes = tau_3_galves, folhas = c('lambda','0', '1', '01'), alfabeto = c('0', '1', '2'), tau = tau_3)#,
  # list(nome = "tau_4", result_selecoes = tau_4_galves, folhas = c('2', '3', '32', '23', '132'), alfabeto = c('0', '1', '2', '3'), tau = tau_4)
)

for(modelo in modelos){
  nome_modelo <- modelo$nome
  result_selecoes <- modelo$result_selecoes
  folhas_internos <- modelo$folhas
  alfabeto <- modelo$alfabeto
  tau <- modelo$tau
  
  acertos_contantes <- matrix(0, nrow = length(limiares),
                              ncol = length(result_selecoes),
                              dimnames = list('c=' %s+% limiares, names(result_selecoes)))
  
  for(limiar in limiares){
    acertos_reteste <- matrix(0, nrow = 100, ncol = length(result_selecoes), dimnames = list(1:100, names(result_selecoes)))
    for(tamanho in names(result_selecoes)){
      for(i in 1:100){
        tam_max <- 6
        deltas <- result_selecoes[[tamanho]][[i]]$deltas
        acertos_reteste[i, tamanho] <- tau_g_deltas(deltas, limiar, alfabeto, folhas_internos, tau, tam_max)
      }
    }
    acertos_contantes['c=' %s+% limiar, ] <- colSums(acertos_reteste)
    print(sprintf("Modelo: %s, Constante: %.2f, Acertos: %s", nome_modelo, limiar, paste(acertos_contantes['c=' %s+% limiar, ], collapse = ", ")))
  }
  
  saveRDS(acertos_contantes, file = sprintf('analise_estimacoes/dados/acertos_%s_cs_galves_fixo_6.RDS', nome_modelo))
}




