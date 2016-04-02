#_________________________________________________________________+
#                                                                 |
# prob_dens_plot.r - create a probability density plot            |
#_________________________________________________________________+
#                                                                 |
# Description: Take the percent divergence in LTR sequences and   |
# create a density plot of divergence values for each             |
# superfamily                                                     |
#                                                                 |
# Author: S. Evan Staton                                          |
# Contact: statonse<at>uga.edu                                    |
# Started: 4.29.10                                               |
# Updated:                                                        |
#                                                                 |
# Suggested                                                  | 
#                                                                 |                                                                                                                                                         |
#                                                              |  
#                                                                 |      
# Requires: 
#_________________________________________________________________+
#TODO: 
#
rm(list=ls())    #this clears the workspace to make sure no leftover variables are floating around
graphics.off()   #close all graphics windows

#library(lattice)

# just using ggplot2 now
library(ggplot2)    #creating the plot layers
library(gridExtra)  #specify the layout of the plots
#library(ggExtra)    #align the plots, as they are on the same scale
library(reshape)

#-------------------------------+
#DIVERGENCE VALUES FOR ALL LTRs |
#-------------------------------+
allLTRdivIN <- read.csv("/Users/spencerstaton/Desktop/Sunflower_Repeat_Project_901/HA_BACs_11-08/RESULTS_11-15/LTRdigest_clustal_alignments/RL_divergence_values.csv",header=TRUE,sep=",")
allLTRdiv <- sort(allLTRdivIN[,1])

#---------------------------------+
#DIVERGENCE VALUES FOR GYPSY LTRs |
#---------------------------------+
gypsyLTRdiv <- sort(allLTRdivIN[,2])

#---------------------------------+
#DIVERGENCE VALUES FOR COPIA LTRs |
#---------------------------------+
copiaLTRdiv <- sort(allLTRdivIN[,16])

#-----------------------------------+
#GET THE MEAN AND SD OF DIVER. VALS |
#-----------------------------------+
#get some basic statistics
#allmn  <- mean(allLTRdiv)
#alldev <- sd(allLTRdiv)

#gypm   <- mean(gypsyLTRdiv)
#gypdev <- sd(gypsyLTRdiv)

#copm   <- mean(copiaLTRdiv)
#copdev <- sd(copiaLTRdiv)

#------------------------+
# CREATE THE DOTCHART    |
#------------------------+
LTR <- read.csv("/Users/spencerstaton/Desktop/Sunflower_Repeat_Project_901/HA_BACs_11-08/RESULTS_11-15/LTRdigest_clustal_alignments/RL_divergence_values3.csv",sep=",",header=TRUE)

#preserve the order of elements
LTR$family <- factor(LTR$family, levels=unique(as.character(LTR$family)))

gpl <- ggplot(LTR, aes(x = divergence, y = LTR$family, colour = superfamily)) + geom_point() 
      
#gpl + scale_colour_manual(values = c("red","blue")) + facet_grid(superfamily ~ ., scales = "free") +
#  scale_x_continuous(breaks=NA) + theme_set(theme_bw()) + ylab(NULL) + xlab(NULL)

# use lattice
#dotplot(x = family ~ divergence, data = LTR, group = superfamily, auto.key = TRUE)

#-------------------------+
# CREATE THE DENSITY PLOT |
#-------------------------+
SUPERFAMS <- read.csv("/Users/spencerstaton/Desktop/Sunflower_Repeat_Project_901/HA_BACs_11-08/RESULTS_11-15/LTRdigest_clustal_alignments/RL_divergence_values4.csv",sep=",",header=TRUE)

#qplot(superfamily, data = SUPERFAMS, geom = "density")

gpldens <- ggplot(SUPERFAMS, aes(divergence,colour = superfamily)) + geom_density() #+ geom_density(linetype = 0, alpha = 0.6) 
           
#gpldens + scale_colour_manual(values = c("black","red","blue")) + theme_set(theme_bw()) + ylab(NULL) + xlab(NULL)

#-------------------------+
# CREATE THE GRID         | 
#-------------------------+
#leg <- ggplotGrob(gpl + opts(keep = "legend_box"))
## one needs to provide the legend with a well-defined width
#legend=gTree(children=gList(leg), cl="legendGrob")
#widthDetails.legendGrob <- function(x) unit(2, "cm")

#postscript("labeled_dotchart_dens_plot.ps")

grid.arrange(# create the dotchart
             gpl
             + facet_grid(superfamily ~ ., scales = "free")
             + scale_colour_manual(values = c("red","blue"))
             #+ scale_x_continuous(breaks=NA) # does indeed remove the breaks
             + xlab(NULL) + ylab(NULL)
             + theme_set(theme_grey())
             + opts(legend.position="none"), #,axis.text.x = theme_blank(), axis.text.y = theme_blank()), 
                    #axis.ticks = theme_blank()),
             # create the density plot
             gpldens
             + scale_colour_manual(values = c("black","red","blue"))
             + xlab(NULL) + ylab(NULL) 
             + theme_set(theme_grey())
             + opts(legend.position="none")) #,axis.text.y = theme_blank())) #,
             #legend=legend,`
             #main ="LTR retrotransposon insertion ages")#,
             #left = "LTR retrotransposon families")

#postscript("aligned_plots.ps")
#pdf("aligned_plots.pdf")

# align.plots just plots over the current layer so clear the graphics device ...
#graphics.off()

# and render the aligned plots
#align.plots(gpl,gpldens)

