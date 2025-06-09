setwd("~/Riquelme/Comparacao_algoritmos")
source("algoritmos/algoritmo_bic.R")

# ÁRVORE 1 =====================================================================
amostras_tau_1 <- readRDS('amostras_geradas/amostras_tau_1.RDS')

result_bic_tau_1 <- setNames(lapply(colnames(amostras_tau_1), function(x) list()), colnames(amostras_tau_1))
alfabeto <- c('0', '1')

for(coluna in colnames(amostras_tau_1)){
  for(linha in rownames(amostras_tau_1)){
    amostra <- amostras_tau_1[linha, coluna]
    c <- 1
    tam_max <- floor(log(nchar(amostra)))
    result_bic_tau_1[[coluna]][[linha]] <- tau_bic(amostra, alfabeto, tam_max, c)
    cat('Árvore 1', coluna, linha, '\n')
  }
}

saveRDS(result_bic_tau_1, file = 'estimacoes/result_bic_tau_1.RDS')

# ÁRVORE 2 =====================================================================
amostras_tau_2 <- readRDS('amostras_geradas/amostras_tau_2.RDS')

result_bic_tau_2 <- setNames(lapply(colnames(amostras_tau_2), function(x) list()), colnames(amostras_tau_2))
alfabeto <- c('0', '1')

for(coluna in colnames(amostras_tau_2)){
  for(linha in rownames(amostras_tau_2)){
    amostra <- amostras_tau_2[linha, coluna]
    c <- 1
    tam_max <- floor(log(nchar(amostra)))
    result_bic_tau_2[[coluna]][[linha]] <- tau_bic(amostra, alfabeto, tam_max, c)
    cat('Árvore 2', coluna, linha, '\n')
  }
}

saveRDS(result_bic_tau_2, file = 'estimacoes/result_bic_tau_2.RDS')

# ÁRVORE 3 =====================================================================
amostras_tau_3 <- readRDS('amostras_geradas/amostras_tau_3.RDS')

result_bic_tau_3 <- setNames(lapply(colnames(amostras_tau_3), function(x) list()), colnames(amostras_tau_3))
alfabeto <- c('0', '1', '2')
for(coluna in colnames(amostras_tau_3)){
  for(linha in rownames(amostras_tau_3)){
    amostra <- amostras_tau_3[linha, coluna]
    c <- 1
    tam_max <- floor(log(nchar(amostra)))
    result_bic_tau_3[[coluna]][[linha]] <- tau_bic(amostra, alfabeto, tam_max, c)
    cat('Árvore 3', coluna, linha, '\n')
  }
}

saveRDS(result_bic_tau_3, file = 'estimacoes/result_bic_tau_3.RDS')

# ÁRVORE 4 =====================================================================
amostras_tau_4 <- readRDS('amostras_geradas/amostras_tau_4.RDS')

result_bic_tau_4 <- setNames(lapply(colnames(amostras_tau_4), function(x) list()), colnames(amostras_tau_4))
alfabeto <- c('0', '1', '2', '3')

for(coluna in colnames(amostras_tau_4)){
  for(linha in rownames(amostras_tau_4)){
    amostra <- amostras_tau_4[linha, coluna]
    c <- 1
    tam_max <- floor(log(nchar(amostra)))
    result_bic_tau_4[[coluna]][[linha]] <- tau_bic(amostra, alfabeto, tam_max, c)
    cat('Árvore 4', coluna, linha, '\n')
  }
}

saveRDS(result_bic_tau_4, file = 'estimacoes/result_bic_tau_4.RDS')
