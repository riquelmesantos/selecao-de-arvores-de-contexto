#  Importa funções auxiliares e pacotes necessários ----
wd <- "~/Riquelme/Comparacao_algoritmos"
source("funcoes_auxiliares.R")
library(dplyr, quietly = T)
library(tidyverse, quietly = T)

# Importa os dados da seleções ----
tau_1_bic <- readRDS("estimacoes/result_bic_tau_1.RDS")
tau_2_bic <- readRDS("estimacoes/result_bic_tau_2.RDS")
tau_3_bic <- readRDS("estimacoes/result_bic_tau_3.RDS")
tau_4_bic <- readRDS("estimacoes/result_bic_tau_4.RDS")

# Função para selecionar a árvore utilizando outros limiares ----
pen <- function(amostra, alfabeto, c = 1){
  a <- length(alfabeto)
  n <- nchar(amostra)
  return(c*(a-1)*log(n))}

calculate_Vw <- function(w, alfabeto, penalizacao, P_ws, tau){
  exp(-penalizacao) * P_ws[w]
}

tau_bic_Pws <- function(log_P_ws, alfabeto, c = 1){
  arvore <- c()
  nos <- names(log_P_ws)
  folhas <- Filter(\(w) any(!filhos(w, alfabeto) %in% nos), nos)
  internos <- setdiff(nos, folhas)
  
  log_Vw <- setNames(rep(-Inf, length(nos)), nos)
  Xw <- setNames(rep(0, length(nos)), nos)
  
  penalizacao <- pen(amostra, alfabeto, c)
  
  log_Vw <- setNames(rep(-Inf, length(nos)), internos) 
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
  return(setequal(tau, arvore))
}

# Árvores modelo---
tau_1 <- c('00', '10', '01', '11')
tau_2 <- c('1', '10', '100', '1000', '0000')
tau_3 <- c('2', '00', '10', '20', '11', '21', '001', '101', '201')
tau_4 <- c('0', '1', '02', '12', '22', '03', '13', '33', '032', '232', '332', '023', '123', '223', '323', '0132', '1132', '2132', '3132')

# Análise das árvores estimadas para diferentes constantes ----
constantes <- 0:30/10

modelos <- list(
  list(nome = "tau_1", result_selecoes = tau_1_bic, folhas = c('0', '1'), alfabeto = c('0', '1'), tau = tau_1),
  list(nome = "tau_2", result_selecoes = tau_2_bic, folhas = c('0', '00', '000'), alfabeto = c('0', '1'), tau = tau_2),
  list(nome = "tau_3", result_selecoes = tau_3_bic, folhas = c('0', '01', '1'), alfabeto = c('0', '1', '2'), tau = tau_3),
  list(nome = "tau_4", result_selecoes = tau_4_bic, folhas = c('2', '3', '32', '23', '132'), alfabeto = c('0', '1', '2', '3'), tau = tau_4)
)

for(modelo in modelos){
  nome_modelo <- modelo$nome
  result_selecoes <- modelo$result_selecoes
  folhas_internos <- modelo$folhas
  alfabeto <- modelo$alfabeto
  tau <- modelo$tau
  
  acertos_contantes <- matrix(0, nrow = length(constantes),
                              ncol = length(result_selecoes),
                              dimnames = list('c=' %s+% constantes, names(result_selecoes)))
  
  for(constante in constantes){
    acertos_reteste <- matrix(0, nrow = 100, ncol = length(result_selecoes), dimnames = list(1:100, names(result_selecoes)))
    for(tamanho in names(result_selecoes)){
      for(i in 1:100){
        log_P_ws <- result_selecoes[[tamanho]][[i]]$log_P_ws
        acertos_reteste[i, tamanho] <- tau_bic_Pws(log_P_ws,  alfabeto, constante)
      }
    }
    acertos_contantes['c=' %s+% constante, ] <- colSums(acertos_reteste)
    print(sprintf("Modelo: %s, Constante: %.2f, Acertos: %s", nome_modelo, constante, paste(acertos_contantes['c=' %s+% constante, ], collapse = ", ")))
  }
  
  saveRDS(acertos_contantes, file = sprintf('analise_estimacoes/dados/acertos_%s_cs_bic.RDS', nome_modelo))
}
