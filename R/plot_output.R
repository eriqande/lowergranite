rm(list=ls())      # Removes current objects that may be in workspace

#Example data
irrigation.df = data.frame(Region = rep(c("Africa","Latin America","North America","Europe"),4), Year = factor(c(rep(1980,4),
                                                                                                                 rep(1990,4), rep(2000,4), rep(2007,4))), Area = c(9.3,12.7,21.2,18.8,11.0,15.5,21.6,25.3,13.2,17.3,
                                                                                                                                                                   23.3,26.7,13.6,17.3,23.8,26.3))

#Example plot
qplot(data = irrigation.df, x = Area, y = Region, colour = Year, main = "Irrigation")

setwd("C:/Users/mackerman/My Documents/GitHub/lowergranite/")
library(ggplot2)
library(gridExtra)

knownSum  <- read.csv("known_summary.csv", header = TRUE)
infSum    <- read.csv("inf_summary.csv",   header = TRUE)

Stock       <- colnames(knownSum)[2:length(colnames(knownSum))]
Simulation  <- c(rep("Known",length(Stock)),rep("Inferred",length(Stock)))
PercentBias <- as.vector(as.numeric(c(knownSum[3,2:(length(Stock)+1)],infSum[3,2:(length(Stock)+1)])))
PctHalfCI   <- as.vector(as.numeric(c(knownSum[10,2:(length(Stock)+1)],infSum[10,2:(length(Stock)+1)])))
CICover     <- as.vector(as.numeric(c(knownSum[8,2:(length(Stock)+1)],infSum[8,2:(length(Stock)+1)])))

data = data.frame(Stock = rep(Stock,2), Simulation = Simulation, PercentBias = PercentBias, 
                  PctHalfCI = PctHalfCI, CICover = CICover)

p1 <- qplot(data = data, x = PercentBias, y = Stock, shape = Simulation, size = I(5)) + 
  scale_x_continuous(breaks = seq(-18,12,4)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  geom_vline(xintercept = 0, linetype = "dotdash") +
  scale_shape_manual(values=c(1,17)) +
  theme(axis.text.x = element_text(size = 12, color = 'black')) +
  theme(axis.text.y = element_text(size = 12, color = 'black')) +
  theme(panel.grid.minor.x = element_line(NA)) +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(1,0,1,1),"lines"))

p2 <- qplot(data = data, x = CICover, y = Stock, shape = Simulation, size = I(5)) + 
  scale_x_continuous(breaks = seq(0,1,0.1)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  geom_vline(xintercept = 0.9, linetype = "dotdash") +
  scale_shape_manual(values=c(1,17)) +
  theme(axis.text.x = element_text(size = 12, color = 'black')) +
  theme(axis.text.y = element_blank()) +
  theme(panel.grid.minor.x = element_line(NA)) +
  theme(axis.title.y = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,-0.80),"lines"))

grid.arrange(p1,p2,ncol=2)
