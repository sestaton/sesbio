library(ggplot2)

data <- read.table("test_11_Lines_cov.tsv",header=T,sep="\t")
#head(data)
#  Line          Ref   Pos  Cov
#1 Hopi scaffold_500 29101 1604
#2 Hopi scaffold_500 29102 1482
#3 Hopi scaffold_500 29103 1262
#4 Hopi scaffold_500 29104 1147
#5 Hopi scaffold_500 29105   66
#6 Hopi scaffold_500 29106   45
p <- ggplot(data, aes(pos, cov)) + geom_point()
p + facet_grid(. ~ Ref)
