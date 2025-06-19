library(ggplot2)
library(dplyr)
library(ggtext)
library(showtext)
showtext_auto()

sequencia_binaria <- function(cadeia, bin) {
  breaks <- seq(0, floor(max(cadeia)) + 1, by = bin)
  categorias <- cut(cadeia[,1], breaks = breaks, right = TRUE, labels = FALSE, include.lowest = TRUE)
  binario <- rep(0, length(breaks) - 1)
  binario[unique(na.omit(categorias))] <- 1
  return(binario)
}

# 
# vars <- c()
# 
# for(i in 1:34){
#   dados_neuronio <- read.table("Conjunto de dados/cell_"%s+%i%s+%".dat")
#   diferenca <- diff(dados_neuronio[,1])
#   nrow(dados_neuronio) - nrow(unique(dados_neuronio))
#   sequencia_vet <- sequencia_binaria(dados_neuronio, 0.01)
#   sequencia_fatiada <- split(sequencia_vet, ceiling(seq_along(sequencia_vet)/600))
#   somas <- sapply(sequencia_fatiada, sum)
#   vars[i] <- var(somas)
#   barplot(somas,  main = as.character(i))
# }

# Registrar a fonte serf (se instalada no sistema)
font_add(family = "serf", regular = "C:/Windows/Fonts/times.ttf")
setwd('C:/Users/rique/Desktop/PIPGES/Dissertacao/Códigos/Comparacao_Estimadores')
dados_neuronio <- read.table("Conjunto de dados/cell_32.dat")
diferenca <- diff(dados_neuronio$V1)
df_diferenca <- data.frame(index = 1:length(diferenca), valor = diferenca)

sequencia_vet <- sequencia_binaria(dados_neuronio, 0.01)
sequencia_fatiada <- split(sequencia_vet, ceiling(seq_along(sequencia_vet)/600))
somas <- sapply(sequencia_fatiada, sum)

df_somas <- data.frame(index = 1:length(somas), valor = somas)

g1 <- ggplot(df_diferenca, aes(x = index, y = valor)) +
  geom_segment(aes(xend = index, yend = 0), color = "#4F71CC") +
  annotate("segment",
           x = 0, xend = max(df_diferenca$index),
           y = mean(diferenca), yend = mean(diferenca),
           linetype = "dashed", color = "#D62828", size = 1) +
  annotate("text",
           x = 2000,  # posição horizontal ajustável
           y = mean(diferenca) + 0.4,           # um pouco acima da linha
           label = "média", color = "#D62828",
           size = 5, family = "serf") +
  labs(
    x = "Ordem do disparo", y = "Intervalo (ms)",
    family = 'serif'
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16, family = "serf"),
    aspect.ratio = 1
  )

plot(g1)

g2 <- ggplot(df_diferenca, aes(x = valor)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#4F71C9", color = "black") +
  geom_density(color = "#D62828", size = 1) +
  scale_y_continuous(limits = c(0, 1.2)) +
  labs(#title = "Histograma e densidade dos intervalos entre disparos",
       x = "Tempo entre disparos", y = "Densidade") +
  theme_minimal()+
  theme(text = element_text(size = 16,  family="serf"),
        aspect.ratio = 1)
plot(g2)

bins <- 1:2000 / 1000 
disparos <- sapply(bins, function(x) sum(sequencia_binaria(dados_neuronio, x)))
tamanhos <- sapply(bins, function(x) length(sequencia_binaria(dados_neuronio, x)))

df_bin <- data.frame(bin = bins, 
                     proporcao = disparos / length(dados_neuronio$V1),
                     tamanho = tamanhos)


g3 <- ggplot(df_bin, aes(x = bin, y = proporcao)) +
  geom_line(color = "#4F71CC", size = 1.5) +
  annotate("segment", x = 0, xend = 0.029, y = 0.85, yend = 0.85,
           linetype = "dashed", color = "#D62828", size = 1) +
  annotate("segment", x = 0.029, xend = 0.029, y = 0, yend = 0.85,
           linetype = "dashed", color = "#D62828", size = 1) +
  geom_point(data = data.frame(x = 0.029, y = 0.85),
             mapping = aes(x = x, y = y),
             color = "#D62828", size = 3)+
  labs(
   # title = expression("Proporção de disparos capturados por tamanho de "*italic("bin")),
    x = expression("Tamanho do "*italic("bin")*" (s)"),
    y = "Proporção de disparos capturados"
  ) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(
    breaks = c(0.029, seq(0.4, 2, by = 0.2)),
    labels = c("<span style='color:#D62828;'>0.029</span>", paste0(seq(0.4, 2, by = 0.2)))
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16, family = "serf"),
    axis.text.x = element_markdown(),
    aspect.ratio = 1
  )
plot(g3)

g4 <- ggplot(df_bin[c(0:100),], aes(x = bin, y = tamanho)) +
  geom_line(color = "#4F71CC", size = 1.5) +
  annotate("segment", x = 0, xend = 0.029, y = 372724, yend = 372724,
           linetype = "dashed", color = "#D62828", size = 1) +
  annotate("segment", x = 0.029, xend = 0.029, y = 0, yend = 372724,
           linetype = "dashed", color = "#D62828", size = 1) +
  geom_point(data = data.frame(x = 0.029, y = 372724),
             mapping = aes(x = x, y = y),
             color = "#D62828", size = 3)+
  labs(
    # title = expression("Proporção de disparos capturados por tamanho de "*italic("bin")),
    x = expression("Tamanho do "*italic("bin")*" (s)"),
    y = "Tamanho da cadeia binária"
  ) +
  scale_x_continuous(
    breaks = c(0, 0.029, 0.05, 0.075,0.1),
    labels = c(0, "<span style='color:#D62828;'>0.029</span>", 0.5, 0.75,1)
    ) +
  scale_y_continuous(
    breaks = c(3e+05, seq(2e+6, 1e+7, by=2e+6))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16, family = "serf"),
    axis.text.x = element_markdown(),
    aspect.ratio = 1
  )
plot(g4)



ggsave("analise_dados_reais/intervalos_disparos.pdf", plot = g1, width = 5, height = 5, dpi = 300)
ggsave("analise_dados_reais/histograma_disparos.pdf", plot = g2, width = 5, height = 5, dpi = 300)
ggsave("analise_dados_reais/proporcao_disparos.pdf", plot = g3, width = 5, height = 5, dpi = 300)
ggsave("analise_dados_reais/tamanho_bin.pdf", plot = g4, width = 5, height = 5, dpi = 300)












