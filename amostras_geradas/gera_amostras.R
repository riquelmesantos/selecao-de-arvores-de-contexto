setwd("C:/Users/rique/Desktop/PIPGES/Dissertacao/Códigos/Comparacao_Estimadores")
source("funcoes_auxiliares.R")

gera_amostras <- function(n_amostras = 100, tam_amostras, min_ocorr, probs, alfabeto, seeds = 1:10e10){
  # Define os contextos (tau) a partir do nome das linhas da matriz de probabilidades
  tau <- rownames(probs)
  filhos_taus <- setdiff(unique(as.vector(sapply(tau, function(w) filhos(w, alfabeto)))),tau)
  
  # Cria uma matriz para armazenar as amostras geradas.
  amostras <- matrix(0, ncol = length(tam_amostras), 
                     nrow = n_amostras, 
                     dimnames = list(1:n_amostras, paste0('n=', tam_amostras)))
  amostra_gerada <- 1
  n_seed <- 1
  # Loop até gerar todas as amostras desejadas
  while(amostra_gerada < n_amostras + 1){
    # Simula uma sequência com comprimento igual ao maior tamanho de amostra desejado
    set.seed(seeds[n_seed])
    amostra <- stri_sub(simula(probs, max(tam_amostras) + 1000*length(alfabeto)), from = 1000*length(alfabeto)+1)
    n_seed <- n_seed + 1
    # Começa com o maior tamanho de amostra e vai descendo
    i <- 1
    while(i <= length(tam_amostras)){
      # Fatia a amostra simulada até o tamanho desejado
      fatiado <- stri_sub(amostra, to = tam_amostras[i])
      cat(n_seed, amostra_gerada, 'fatiou uma amostra de tamanho', nchar(fatiado))
      # Conta quantas vezes cada contexto (tau) aparece na fatia
      cont <- rowSums(contagens(tau, alfabeto, fatiado))
      cont_filhos <- rowSums(contagens(filhos_taus, alfabeto, fatiado))
      # Verifica se todos os contextos têm pelo menos a ocorrência mínima requerida
      if(all(min_ocorr[i] <= cont) && all(1 <= cont_filhos)){
        amostras[amostra_gerada, i] <- fatiado   # Se sim, salva a fatia correspondente na matriz
        cat(' e aceitou. \n')
        i <- i + 1
      } 
      else {
        cat(' e descartou. \n')
        break
      }
    }
    print(amostra_gerada)
    # Se a amostra foi aceita em todos os níveis de tamanho, passa pra próxima
    if(all(amostras[amostra_gerada, ] != 0)) 
      amostra_gerada <- amostra_gerada + 1
  }
  
  # Retorna a matriz de amostras geradas
  return(amostras)
}

probs_1 <- matrix(
  c(0.1, 0.9, 
    0.7, 0.3, 
    0.4, 0.6, 
    0.9, 0.1),
  nrow = 4,
  byrow = TRUE,
  dimnames = list(c('00', '10', '01', '11'), c('0', '1'))
)

probs_1
n_amostras <- 100
tam_amostras_1 <- c(200, 400, 600, 800, 1000, 1200, 1400, 1600)
min_ocorr_1 <- floor(tam_amostras_1 / log(tam_amostras_1) / 4)
seeds <- 1:10e10
alfabeto_1 <- c('0', '1')

amostras_1 <- gera_amostras(n_amostras, tam_amostras_1, min_ocorr_1, probs_1, alfabeto_1)
saveRDS(amostras_1, file = 'amostras_geradas/amostras_tau_1.RDS')


probs_2 <- matrix(
  c(0.6, 0.4, 
    0.4, 0.6, 
    0.3, 0.7, 
    0.2, 0.8, 
    0.1, 0.9),
  nrow = 5,
  byrow = TRUE,
  dimnames = list(c('1', '10', '100', '1000', '0000'), c('0', '1'))
)

probs_2
n_amostras <- 100
tam_amostras_2 <- c(10000, 14000, 20000, 24000, 28000, 32000, 36000, 40000)
min_ocorr_2 <- floor(tam_amostras_2 / log(tam_amostras_2) / 16)
seeds <- 1:10e10
alfabeto_2 <- c('0', '1')

amostras_2 <- gera_amostras(n_amostras, tam_amostras_2, min_ocorr_2, probs_2, alfabeto_2)
saveRDS(amostras_2, file = 'amostras_geradas/amostras_tau_2.RDS')

probs_3 <- matrix(
  c(0.3,0.6,0.1,
    0.1,0.3,0.6,
    0.2,0.5,0.3,
    0.7,0.2,0.1,
    0.4,0.4,0.2,
    0.8,0.1,0.1,
    0.1,0.2,0.7,
    0.4,0.5,0.1,
    0.5,0.3,0.2),
  nrow = 9,
  byrow = TRUE,
  dimnames = list(
    c('2', '00', '10', '20', '11', '21', '001', '101', '201'),
    c('0', '1', '2')
  )
)

probs_3
n_amostras <- 100
tam_amostras_3 <- c(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000)
min_ocorr_3 <- floor(tam_amostras_3 / log(tam_amostras_3)/ 9)
seeds <- 1:10e10
alfabeto_3 <- c('0', '1', '2')

amostras_3 <- gera_amostras(n_amostras, tam_amostras_3, min_ocorr_3, probs_3, alfabeto_3)
saveRDS(amostras_3, file = 'amostras_geradas/amostras_tau_3.RDS')

probs_4 <- matrix(
  c(0.2, 0.3, 0.1, 0.4,
    0.1, 0.1, 0.6, 0.2,
    0.4, 0.3, 0.2, 0.1,
    0.1, 0.2, 0.4, 0.3,
    0.3, 0.1, 0.2, 0.4,
    0.1, 0.2, 0.1, 0.6,
    0.2, 0.4, 0.3, 0.1,
    0.5, 0.1, 0.2, 0.2,
    0.2, 0.2, 0.5, 0.1,
    0.1, 0.3, 0.4, 0.2,
    0.6, 0.1, 0.1, 0.2,
    0.7, 0.1, 0.1, 0.1,
    0.1, 0.5, 0.2, 0.2,
    0.3, 0.3, 0.3, 0.1,
    0.2, 0.2, 0.1, 0.5,
    0.1, 0.7, 0.1, 0.1,
    0.1, 0.1, 0.2, 0.6,
    0.3, 0.3, 0.1, 0.3,
    0.2, 0.2, 0.5, 0.1),
  nrow = 19,
  byrow = T,
  dimnames = list(
    c('0', '1', '02', '12', '22', '03', '13', '33', '032', '232', '332', '023', '123', '223', '323', '0132', '1132', '2132', '3132'),
    c('0', '1', '2', '3')
  ))
probs_4
rowSums(probs_4)
n_amostras <- 100
tam_amostras_4 <- c(5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000)
min_ocorr_4 <- floor(tam_amostras_4 / log(tam_amostras_4) / 128)
seeds <- 1:10e10
alfabeto_4 <- c('0', '1', '2','3')

amostras_4 <- gera_amostras(n_amostras, tam_amostras_4, min_ocorr_4, probs_4, alfabeto_4)
saveRDS(amostras_4, file = 'amostras_geradas/amostras_tau_4.RDS')
