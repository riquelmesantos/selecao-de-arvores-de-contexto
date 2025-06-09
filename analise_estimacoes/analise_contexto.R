#  Importa funções auxiliares e pacotes necessários ----
wd <- "~/Riquelme/Comparacao_algoritmos"
source("funcoes_auxiliares.R")
library(dplyr, quietly = T)
library(tidyverse, quietly = T)

# Importa os dados da seleções ----
tau_1_contexto <- readRDS("estimacoes/result_contexto_tau_1.RDS")
tau_2_contexto <- readRDS("estimacoes/result_contexto_tau_2.RDS")
tau_3_contexto <- readRDS("estimacoes/result_contexto_tau_3.RDS")
tau_4_contexto <- readRDS("estimacoes/result_contexto_tau_4.RDS")

# Função para selecionar a árvore utilizando outros limiares ----
tau_c_deltas <- function(deltas, limiar, folhas_internos) {
  deltas <-  deltas[!is.na(deltas)]
  nos <- names(deltas)
  externos <- setdiff(nos, folhas_internos)
  if(all(deltas[externos] < limiar) && all(deltas[folhas_internos] > limiar)) {
    return(1)
  }
  else return(0)
}

# Árvores modelo---
tau_1 <- c('00', '10', '01', '11')
tau_2 <- c('1', '10', '100', '1000', '0000')
tau_3 <- c('2', '00', '10', '20', '11', '21', '001', '101', '201')
tau_4 <- c('0', '1', '02', '12', '22', '03', '13', '33', '032', '232', '332', '023', '123', '223', '323', '0132', '1132', '2132', '3132')

# Análise das árvores estimadas para diferentes constantes ----
constantes <- 1:30/10

modelos <- list(
  list(nome = "tau_1", result_selecoes = tau_1_contexto, folhas = c('0', '1')),
  list(nome = "tau_2", result_selecoes = tau_2_contexto, folhas = c('0', '00', '000')),
  list(nome = "tau_3", result_selecoes = tau_3_contexto, folhas = c('0', '01', '1')),
  list(nome = "tau_4", result_selecoes = tau_4_contexto, folhas = c('2', '3', '32', '23', '132'))
)

for(modelo in modelos){
  nome_modelo <- modelo$nome
  result_selecoes <- modelo$result_selecoes
  folhas_internos <- modelo$folhas
  
  acertos_contantes <- matrix(0, nrow = length(constantes),
                              ncol = length(result_selecoes),
                              dimnames = list('c=' %s+% constantes, names(result_selecoes)))
  
  for(constante in constantes){
    acertos_reteste <- matrix(0, nrow = 100, ncol = length(result_selecoes), dimnames = list(1:100, names(result_selecoes)))
    for(tamanho in names(result_selecoes)){
      for(i in 1:100){
        n <- as.numeric(stri_sub(tamanho, from = 3))
        limiar <- constante * log(n)
        deltas <- result_selecoes[[tamanho]][[i]]$deltas
        acertos_reteste[i, tamanho] <- tau_c_deltas(deltas, limiar, folhas_internos)
      }
    }
    acertos_contantes['c=' %s+% constante, ] <- colSums(acertos_reteste)
    print(sprintf("Modelo: %s, Constante: %.2f, Acertos: %s", nome_modelo, constante, paste(acertos_contantes['c=' %s+% constante, ], collapse = ", ")))
  }
  
  saveRDS(acertos_contantes, file = sprintf('analise_estimacoes/dados/acertos_%s_cs_contexto.RDS', nome_modelo))
}




