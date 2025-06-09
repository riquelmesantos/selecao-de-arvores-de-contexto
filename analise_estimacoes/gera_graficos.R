# Importa pacotes e arquivos
wd <- "C:/Users/rique/Desktop/PIPGES/Dissertacao/Códigos/Comparacao_Estimadores"
setwd(wd)
library(dplyr, quietly = T)
library(tidyverse, quietly = T)
library(ggplot2)
library(latex2exp)

acertos_contantes <- readRDS('analise_estimacoes/dados/acertos_tau_3_cs_galves_fixo_4.RDS')
colnames(acertos_contantes)[c(1,6)] <- c('n=100000', 'n=200000')
df_const <- rownames_to_column(as.data.frame(acertos_contantes), var = "c")
df_long <- pivot_longer(df_const, cols = -c, names_to = "n", values_to = 'n_acertos')
df_long$n <- factor(df_long$n, levels = unique(df_long$n))
df_long$c <- as.numeric(gsub("c=", "", df_long$c))

cores <- colorRampPalette(c("#265F01", "#509F20", "#47BF89", "#00B84D", 
                            "#00B2B5", "#0E7ACC", "#4F71CC", "#8B64CC"))(8)

grafico <- ggplot(df_long, aes(x = c, y = n_acertos, color = n, group = n)) +
  geom_line(linewidth = 1.2, alpha = 0.6) +
  ggtitle(TeX(r"(Altura máxima fixa em 4)"))+
  scale_x_continuous(
    name = TeX(r"(Valor de $\delta$)"),
    breaks = seq(0, 1, by = 0.2)) +
  scale_color_manual(
    name = element_blank(),
    values = cores) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(y = "Número de acertos") +
  theme_minimal() +
  theme(legend.position = "bottom",
        text = element_text(size = 24,  family="serif"),
        aspect.ratio = 1)

print(grafico)
ggsave("analise_estimacoes/graficos/Proporcao_acertos_galves_tau_3_fixo_4.pdf",
       plot = grafico, width = 7, height = 8, dpi = 300)


