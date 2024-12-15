library("ggplot2")
library("ggrepel")
library("pheatmap")

exprs <- read.csv("D:/McGill/phd course/single-cell/final/GSE115978/PEX3_exhaustion.csv")

cor = cor.test(exprs$Exhaustion_score, exprs$PEX3.percentage)$estimate
p_value = cor.test(exprs$Exhaustion_score, exprs$PEX3.percentage)$p.value

labels = data.frame(
  x = rep(0, 2),
  y = c(1.2, 1.1),
  labels = rep( c( paste0('r = ', round(cor, 2)),
                   paste0('p.value = ', formatC(p_value))), 4)
)

add_label = geom_text(data = labels, aes(x = x, y = y, label = labels), size = 4, hjust = 0)

ggplot(data = exprs, aes(x = PEX3.percentage, y = Exhaustion_score)) +
  geom_point(size = 3.5, color = '#9B70F9', shape = 1, stroke = 1) +
  theme_bw() +
  xlab('PEX3.percentage') +
  ylab('Exhaustion_score') +
  geom_smooth(method = 'lm',color = '#FD6F71',fill = '#FD6F71',alpha = 0.25)+
  add_label +
  ggtitle("Melanoma Patients")
