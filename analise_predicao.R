# Sequência de teste
amostra_teste <- "001110110101100100010011001001001001001100111011001100110010010100101100010011001001100110011001001001100111001100100101010010110101100010011001001001100011001001100101101000110011101001001111011001101"
seq_teste <- strsplit(amostra_teste, "")[[1]]
# Transições estimadas
transicao_estimada <- matrix(c(
  0.9, 0.1,
  0.6, 0.4,
  0.3, 0.7,
  0.1, 0.9
), nrow=4, byrow=TRUE)
rownames(transicao_estimada) <- c("00", "01", "10", "11")
colnames(transicao_estimada) <- c("0", "1")

# Função pra calcular log-verossimilhança
log_verossim <- function(seq, trans_est) {
  ll <- 0
  for (i in 3:length(seq)) {
    contexto <- paste0(seq[i-2], seq[i-1])
    simbolo <- seq[i]
    prob <- trans_est[contexto, simbolo]
    if (prob > 0) {
      ll <- ll + log(prob)
    } else {
      # Se a transição nunca foi observada no treino, penaliza hard
      ll <- ll + log(1e-10)  # pseudo-count mínimo
    }
  }
  return(ll)
}

log_likelihood <- log_verossim(seq_teste, transicao_estimada)
cat("Log-verossimilhança na sequência de teste:", log_likelihood, "\n")

