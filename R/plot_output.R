# Function for plotting results contained in sumry table from 2 simulations
#
#'

rm(list=ls())

library(ggplot2)
library(gridExtra)

SUMRY.DIR <- "S:/Eagle Fish Genetics Lab/Ackerman/LGR/Kirk/ms/finalSimulations/lowergranite_big_run_CincoDeMayo_2014/"
setwd(SUMRY.DIR)
big_run <- slurp_all_results(dir(pattern = "^[CS].*"))

k <- big_run$SH11SIMPOP_StockAge.xlsx_9_500_500_FALSE$sumrys
i <- big_run$SH11SIMPOP_StockAge.xlsx_9_500_500_TRUE$sumrys

#plot_output <- function(
#  SUMRY.DIR = getwd(),
#  x = NULL,
#  y = NULL)

#k <- read.csv("SH11SIMPOP_StockAge.xlsx_9_500_500_FALSE/summary.csv", header = TRUE, row.names = 1)
#i <- read.csv("SH11SIMPOP_StockAge.xlsx_9_500_500_TRUE/summary.csv", header = TRUE, row.names = 1)

# Create data frame of results to plot
Groups       <- colnames(k)
Simulation   <- c(rep("Known",length(Groups)),rep("Inferred",length(Groups)))
Truth        <- as.vector(as.numeric(c(k[1,1:length(Groups)],i[1,1:length(Groups)])))
PercentBias  <- as.vector(as.numeric(c(k[3,1:length(Groups)],i[3,1:length(Groups)])))
PctHalfCI    <- as.vector(as.numeric(c(k[10,1:length(Groups)],i[10,1:length(Groups)])))
CICover      <- as.vector(as.numeric(c(k[8,1:length(Groups)],i[8,1:length(Groups)])))

data = data.frame(Group = rep(Groups,2), Simulation = Simulation, TruthNum = Truth,
                  PercentBias = PercentBias, PctHalfCI = PctHalfCI, CICover = CICover)

data = data[data$TruthNum != 0,]
data$Truth = format(data$TruthNum, digits = 1)
data <- droplevels(data)
data$Zero = 0

grpNames   <- as.character(data$Group[data$Simulation=="Inferred"])
ord        <- order(data$CICover[data$Simulation=="Inferred"], decreasing = FALSE)
data$Group <- factor(data$Group, levels = grpNames[ord])

# Truth Text Plot
p1 <- ggplot(data = data, aes(x = Zero, y = Group)) +
  xlab("Truth") +
  #theme(panel.border = element_rect(colour = "black")) +
  geom_text(data = data, aes(x = Zero, label = Truth), size=4.3) +
  geom_hline(yintercept = 45.6, linetype = "solid") +
  #geom_hline(yintercept = length(Groups) + 0.6, linetype = "solid") +
  theme(axis.text.x = element_blank()) + 
  theme(axis.text.y = element_text(size = 11, color = 'black')) +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(1,5.5,1.5,1),"lines")) +
  theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(panel.background = element_blank())
  
# Pct Bias
p2 <- qplot(data = data, x = PercentBias, y = Group, shape = Simulation, size = I(4)) +
  scale_x_continuous(limits = c(-25,25), breaks = seq(-25,25,5)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  geom_vline(xintercept = 0, linetype = "solid") +
  geom_vline(xintercept = 10, linetype = "dotdash") +
  geom_vline(xintercept = -10, linetype = "dotdash") +
  scale_shape_manual(values=c(1,16)) +
  theme(axis.text.x = element_text(size = 12, color = 'black')) +
  theme(axis.text.y = element_text(size = 12, color = 'black')) +
  theme(panel.grid.minor.x = element_line(NA)) +
  theme(legend.position = "none") +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(plot.margin = unit(c(1,6,1,-6.5),"lines")) + 
  #theme(axis.line = element_line(colour = "black")) 
  theme(panel.grid.major = element_blank())
  #theme(panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(colour = "black")) +
  #theme(panel.background = element_blank())

# CI Coverage Plot
p3 <- qplot(data = data, x = CICover, y = Group, shape = Simulation, size = I(4)) +
  #xlim(0,1) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  geom_vline(xintercept = 0.9, linetype = "dotdash") +
  scale_shape_manual(values=c(1,16)) +
  theme(axis.text.x = element_text(size = 12, color = 'black')) +
  theme(axis.text.y = element_blank()) +
  theme(panel.grid.minor.x = element_line(NA)) +
  theme(axis.title.y = element_blank()) +
  theme(plot.margin = unit(c(1,1,1,-7),"lines")) + 
  #theme(axis.line = element_line(colour = "black")) +
  theme(panel.grid.major = element_blank())
  #theme(panel.grid.minor = element_blank()) +
  #theme(panel.border = element_blank()) +
  #theme(panel.background = element_blank())

grid.arrange(p1,p2,p3,ncol=3)
