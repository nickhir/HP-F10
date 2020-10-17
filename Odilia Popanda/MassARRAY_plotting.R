library(dplyr)
library(ggpubr)
library(ggplot2)
library(pheatmap)
library(viridis)

# lig4 aci ----
lig4_aci <- data.frame(Genotype=c("WT","DKO"), 
                       methylation=c(78,11),
                       sd=c(9.43,5.12)
                       )


# svg("lig4 aci.svg")
ggplot(lig4_aci, aes(x = Genotype, y = methylation)) +
    geom_col(width=0.5, fill="black") +
    geom_errorbar(data=lig4_aci, aes(x=Genotype, ymin=methylation-sd, ymax=methylation+sd),
                  width=0.35, size=1.2)+
    ylab("Methylation [%]") +
    scale_fill_manual(values = c("black", "#818181")) +
    theme_classic(base_size = 17)+
    grids(linetype="dashed")+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100), breaks = c(0,20,40,60,80,100))
# dev.off()
    



# lig4 hyp ----
lig4_hyp <- data.frame(Genotype=c("WT","DKO"), 
                       methylation=c(91,5),
                       sd=c(2,8)
)

svg("lig4 hyp.svg")
ggplot(lig4_hyp, aes(x = Genotype, y = methylation)) +
    geom_col(width=0.5, fill="black") +
    geom_errorbar(data=lig4_hyp, aes(x=Genotype, ymin=methylation-1, ymax=methylation+sd),
                  width=0.35, size=1.2)+
    ylab("Methylation [%]") +
    scale_fill_manual(values = c("black", "#818181")) +
    theme_classic(base_size = 17)+
    grids(linetype="dashed")+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100), breaks = c(0,20,40,60,80,100))
dev.off()
















# line1 aci ----
line1 <- data.frame(Genotype=c("WT","DKO"), 
                       methylation=c(65,44),
                       sd=c(3,4)
)

svg("line1.svg")
ggplot(line1, aes(x = Genotype, y = methylation)) +
    geom_col(width=0.5, fill="black") +
    geom_errorbar(data=line1, aes(x=Genotype, ymin=methylation-sd, ymax=methylation+sd),
                  width=0.35, size=1.2)+
    ylab("Methylation [%]") +
    scale_fill_manual(values = c("black", "#818181")) +
    theme_classic(base_size = 17)+
    grids(linetype="dashed")+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100), breaks = c(0,20,40,60,80,100))
dev.off()



# line1 aci ----
sat <- data.frame(Genotype=c("WT","DKO"), 
                    methylation=c(92,7),
                    sd=c(1,3)
)

# svg("sat.svg")
ggplot(sat, aes(x = Genotype, y = methylation)) +
    geom_col(width=0.5, fill="black") +
    geom_errorbar(data=sat, aes(x=Genotype, ymin=methylation-sd, ymax=methylation+sd),
                  width=0.35, size=1.2)+
    ylab("Methylation [%]") +
    theme_classic(base_size = 17)+
    grids(linetype="dashed")+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100), breaks = c(0,20,40,60,80,100))
dev.off()


# individual_lig4aci ----


individual_lig4aci <- data.frame(position = c("CpG_2","CpG_3", "CpG_4.5"),
                         WT = c(71,71,91),
                         DKO = c(7,7,18))

long_df <- melt(individual_lig4aci, id.vars = "position")

colnames(long_df) <- c("position", "Genotype", "Methylation")


# svg("individual lig4 aci.svg")
ggplot(data=long_df, aes(x=position, y=Methylation, group=Genotype, color=Genotype, shape=Genotype))+
    geom_point(size=3)+
    geom_line(size = 1)+
    theme_classic(base_size = 17)+
    grids(linetype="dashed")+
    xlab("")+
    ylab("Methylation [%]")+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100), breaks = c(0,20,40,60,80,100))

# dev.off()


# individual_lig4hyp ----
individual_lig4hyp <- data.frame(position = c("CpG_2","CpG_3", "CpG_4.5"),
                                 WT = c(92,92,88),
                                 DKO = c(0,0,16))

long_df <- melt(individual_lig4hyp, id.vars = "position")

colnames(long_df) <- c("position", "Genotype", "Methylation")


# svg("individual_lig4hyp.svg")
ggplot(data=long_df, aes(x=position, y=Methylation, group=Genotype, color=Genotype, shape=Genotype))+
    geom_point(size=3)+
    geom_line(size = 1)+
    theme_classic(base_size = 17)+
    grids(linetype="dashed")+
    xlab("")+
    ylab("Methylation [%]")+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100), breaks = c(0,20,40,60,80,100))
# dev.off()


# individual Line1 ----
individual_line1 <- data.frame(position = c("CpG_1","CpG_2", "CpG_3"),
                               WT = c(63,69,63),
                               DKO = c(38,45,48))


long_df <- melt(individual_line1, id.vars = "position")

colnames(long_df) <- c("position", "Genotype", "Methylation")


# svg("individual_line1.svg")
ggplot(data=long_df, aes(x=position, y=Methylation, group=Genotype, color=Genotype, shape=Genotype))+
    geom_point(size=3)+
    geom_line(size = 1)+
    theme_classic(base_size = 17)+
    grids(linetype="dashed")+
    xlab("")+
    ylab("Methylation [%]")+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100), breaks = c(0,20,40,60,80,100))
# dev.off()




# individual sat ----
individual_sat <- data.frame(position = c("CpG_3","CpG_4"),
                               WT = c(91,92),
                               DKO = c(4,9))


long_df <- melt(individual_sat, id.vars = "position")

colnames(long_df) <- c("position", "Genotype", "Methylation")


# svg("individual_sat.svg")
ggplot(data=long_df, aes(x=position, y=Methylation, group=Genotype, color=Genotype, shape=Genotype))+
    geom_point(size=3)+
    geom_line(size = 1)+
    theme_classic(base_size = 17)+
    grids(linetype="dashed")+
    xlab("")+
    ylab("Methylation [%]")+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100), breaks = c(0,20,40,60,80,100))
# dev.off()




# plot heatmap ----

matrix <- data.frame(row.names = c("0%", "20%", "40%", "60%", "80%", "100%", "WT", "DKO"),
                     LIG4_CpG2 = c(18,20,35,46,70,84,71,7),
                     LIG4_CpG3 = c(18,20,35,46,70,84,71,7),
                     LIG4_CpG4.5 = c(4,30,29,55,71,93,91,18),
                     Line1_CpG1 = c(2,20,31,49,60,69,63,38),
                     Line1_CpG2= c(8,23,34,50,59,74,69,45),
                     Line1_CpG3 = c(11,23,34,48,51,59,63,48),
                     SAT_CpG3 = c(0,31,57,71,77,100,91,4),
                     SAT_CpG4 = c(10, 30,47,63,78,86,92,9)
                     )

# svg("heatmap.svg", height=13, width = 13)
pheatmap(mat=matrix, 
         color = plasma(10),
         border_color =  "NA",
         fontsize=15,
         cluster_rows = F, 
         cluster_cols = F,
         angle=45,
         na_col = "black",
         cellheight = 30,
         cellwidth = 60
         )

# dev.off()














