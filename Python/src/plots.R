# Title     : Genetic Algorithm Plot
# Objective : Plot the generated data of the Genetic Algorithm
# Created by: filipefalcao
# Created on: 2019-04-02

# Libs
library(ggplot2)
library(ggpubr)

# Read input data
input_data <- read.csv(file = "msa-ga/data/output_data.csv", header = TRUE, sep = ",")
mutation_data <- read.csv(file = "msa-ga/data/output_data2.csv", header = TRUE, sep = ",")

# Plots 1 e 2
plot1 <- ggplot(input_data, aes(x = gens, y = result)) + 
  geom_line() +
  geom_point() +
  ggtitle("Gerações vs. Função de Avaliação") +
  xlab("Número de Gerações") + ylab("Função de Avaliação")

plot2 <- ggplot(input_data, aes(x = chromosomes, y = result)) + 
  geom_line(linetype = "dashed") +
  geom_point() +
  ggtitle("Cromosomos vs. Função de Avaliação") +
  xlab("Número de Cromosomos") + 
  theme(axis.title.y = element_blank())

# Arrange
ggarrange(plot1, plot2)

# Plot 3
plot3 <- ggplot(mutation_data, aes(x = mutation, y = result)) + 
  geom_line() +
  geom_point() +
  ggtitle("Taxa de Mutação vs. Função de Avaliação") +
  xlab("Taxa de Mutação") + ylab("Função de Avaliação") +
  scale_y_continuous(limits = c(450, NA))
