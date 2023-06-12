## load libraries
library("dplyr")
library("ggplot2")

## Read args from command line
args = commandArgs(trailingOnly=TRUE)

## Uncomment For debugging only
## Comment for production mode only
## args[1] <- "./parallel_plot.PCA_df.tsv" # formated file
## args[2] <- "./snps_noks.eval" ## eval file


## Passing args to named objects
pcadf_file <- args[1]
signifpc_file <- args[2]
out_file <- "./pca_1_2.png"

eval <- read.csv(file = signifpc_file, header = F, sep = "\t") %>%
  rename(eigenvalue=V1) %>%
  mutate(PC=paste0("PC",row_number()), Variance=eigenvalue/sum(eigenvalue)) %>%
  mutate(Percentage_variance=round(100*(Variance), digits = 2))

####
## Automate PCA generation ----
## read data
pca_data.df <- read.table(file = pcadf_file,
                          header = T,
                          sep = "\t", stringsAsFactors = T) %>%
  mutate(PC1=PC1*(-1))

my_pal <- c("Chane" = "#FF6F00B2",
            "Colombian" = "#C71000B2",
            "Karitiana" = "#008EA0B2",
            "Maya" = "#8A4198B2", 
            "Mixe" = "#5A9599B2",
            "Mixtec" = "#FF6348B2",
            "Pima" = "#84D7E1B2",
            "Quechua" = "#FF95A8B2",
            "Surui" = "#3D3B25B2",
            "Zapotec" = "#ADE2D0B2")
## Uncomment if you want to ommit populations
my_pal <- my_pal[! names(my_pal) %in% c('Karitiana','Surui')]
fig_pal <- c("Anzick-1" = 1,
             "Ayayema" = 2,
             "Lovelock2" = 3,
             "Lovelock3" = 4,
             "PC537" = 5,
             "SpiritCave" = 6,
             "Sumidouro5" = 7,
             "USR1" = 8,
             "Aconcagua" = 9,
             "BigBar" = 10,
             "CAA001" = 11,
             "CDE001" = 12,
             "CK-13" = 13,
             "EFC019" = 14,
             "Kennewick" = 15,
             "Lovelock1" = 16,
             "Lovelock4" = 17,
             "PuntaSantaAna" = 18,
             "Sumidouro4" = 19,
             "Sumidouro6" = 20,
             "Sumidouro7" = 21,
             "Sumidouro8" = 22,
             "TrailCreek" = 23)
## Uncomment if you want to omit low cov
fig_pal <- fig_pal[! names(fig_pal) %in% c("Aconcagua",
                                           "BigBar",
                                           "CAA001",
                                           "CDE001",
                                           "CK-13",
                                           "EFC019",
                                           "Kennewick",
                                           "Lovelock1",
                                           "Lovelock4",
                                           "PuntaSantaAna",
                                           "Sumidouro4",
                                           "Sumidouro6",
                                           "Sumidouro7",
                                           "Sumidouro8",
                                           "TrailCreek")]
                                    
my_plot.p <- ggplot(data = pca_data.df, mapping = aes(x = PC1, y = PC2)) +
  geom_point(data = ~subset(., tag != "Ancient"), mapping = aes(color = tag), shape = 20, size = 3) +
  geom_point(data = ~subset(., tag == "Ancient"), mapping = aes(shape = sample), fill = "#93d9a0", size = 3.5) +
  scale_color_manual(name = "Group", values = my_pal) +
  scale_shape_manual(name = "Ancient Sample", values = fig_pal) +
  xlab(paste(eval$PC[1]," (",eval$Percentage_variance[1],"%)")) + 
  ylab(paste(eval$PC[2]," (",eval$Percentage_variance[2],"%)")) +
  theme_linedraw() +
  theme(
        legend.title = element_blank(),
        legend.text = element_text(family = "Helvetica Neue Light", size = 10),
        axis.text = element_text(family = "Helvetica Neue Light", size = 9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black"),
        title = element_text(family = "Helvetica Neue Light", size = 10),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent', colour = NA),
        legend.spacing.y = unit(0, 'cm'),
        legend.key = element_rect(fill = "transparent")) 


ggsave(filename = out_file, plot = my_plot.p, device = "png", dpi = 300, width = 8, height = 8, units = "in", bg = "transparent")
