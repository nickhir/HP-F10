library(ggplot2)
library(readxl)
library(dplyr)
library(ggpubr)
# read in excel file which contains qPCR results (provided by Supervisor)
setwd("data/")
excel <- read_excel(path = "Course_MCiPqPCR9_9_20.xls")

# select only columns that are relevant for us (Pos, Sample, Ct)
# also select rows that are relevant for us. In this case all rows containing A (standard), B(SNRPN)
# and H (spike in)

standard_gene <- excel %>%
    select(c("Pos", "Sample", "Ct")) %>%
    filter(grepl("A\\d+", Pos))

standard_spike <- excel %>%
    select(c("Pos", "Sample", "Ct")) %>%
    filter(grepl("F\\d+", Pos))

gene <- excel %>%
    select(c("Pos", "Sample", "Ct")) %>%
    filter(grepl("B\\d+", Pos))

spike <- excel %>%
    select(c("Pos", "Sample", "Ct")) %>%
    filter(grepl("J\\d+", Pos))

# for standard we have to convert column Sample to double instead of chr
standard_gene$Sample <- as.double(standard_gene$Sample)


# calculate regression equation
lm_eqn <- lm(standard_gene$Ct ~ standard_gene$Sample)

metrics <- list(
    a = format(unname(coef(lm_eqn)[1]), digits = 2),
    b = format(unname(coef(lm_eqn)[2]), digits = 2),
    r2 = format(summary(lm_eqn)$r.squared, digits = 3)
)

eq <- substitute(italic(C[T]) == a * " " * b  * italic(x) * ", " ~  ~ italic(R) ^ 2 ~ "=" ~ r2,
               metrics)

equation <- as.character(as.expression(eq))


# plot the standard (row A)

# svg("images/StandardCurve.svg")
ggplot(standard_gene, aes(x = Sample, y = Ct)) +
    geom_smooth(method = "lm", se = T) +
    geom_point() +
    xlab("Log Concentration") +
    ylab(bquote(C[T])) +
    geom_text(
        x = -0.6,
        y = 28.21,
        label = equation,
        parse = T,
        size = 6.7
    ) +
    theme_classic(base_size = 14.5)+
    scale_y_continuous(limits = c(20, 30))+
    grids(linetype="dashed")+
    ggtitle("Standard SNRPN")+
    theme(plot.title = element_text(hjust = 0.5))

# dev.off()

# use the regression equation (CT-20.967/-2.945 = Conc) to calculate the log concentration of our gene and spike in
gene["LogConcentration"] <- (gene$Ct - 20.967) / -2.945

spike["LogConcentration"] <- (spike$Ct - 20.967) / -2.945

# calculate relative concentrattion

gene["RelativeConcentration"] <- 10 ^ gene$LogConcentration

spike["RelativeConcentration"] <- 10 ^ spike$LogConcentration

# calculate average for technical duplicates

gene <- gene %>%
    group_by(Sample) %>%
    mutate(AverageRelConc = mean(RelativeConcentration, na.rm = T)) %>%
    mutate(sd = sd(RelativeConcentration, na.rm = T))

spike <- spike %>%
    group_by(Sample) %>%
    mutate(AverageRelConc = mean(RelativeConcentration, na.rm = T)) %>%
    mutate(sd = sd(RelativeConcentration, na.rm = T))






# remove genotype from sample name and move it to new column
gene$Sample <- sub("wt\\d+|dko\\d+", "", gene$Sample)
gene["Genotype"] <- c(rep("WT", 10), rep("DKO", 10))

spike$Sample <- sub("wt\\d+|dko\\d+", "", spike$Sample)
spike["Genotype"] <- c(rep("WT", 10), rep("DKO", 10))

# rearrange order of the plot
gene$Genotype <- factor(gene$Genotype, levels = c("WT", "DKO"))

spike$Genotype <- factor(spike$Genotype, levels = c("WT", "DKO"))




# plot the average Rel concentration for the gene

# svg("images/SNRP.svg")
ggplot(gene, aes(x = Sample, y = AverageRelConc, fill = Genotype)) +
    geom_errorbar(
        aes(ymin = AverageRelConc - sd, ymax = AverageRelConc + sd),
        position = position_dodge(0.45),
        width = 0.3,
        size = 1.23
    ) +
    geom_bar(stat = "identity",
             position = position_dodge(),
             width = 0.45) +
    ylab("Relative Concentration") +
    scale_fill_manual(values = c("black", "#818181")) +
    scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, 0.5),
        breaks = seq(0, 0.5, by = 0.1)
    )+
    theme_classic(base_size = 14.5)+
    grids(linetype="dashed")+
    ggtitle("Enrichment SNRPN")+
    theme(plot.title = element_text(hjust = 0.5))

# dev.off()

# plot the average Rel concentration for the spike

# svg("images/spike.svg")
ggplot(spike, aes(x = Sample, y = AverageRelConc, fill = Genotype)) +
    geom_errorbar(
        aes(ymin = AverageRelConc - sd, ymax = AverageRelConc + sd),
        position = position_dodge(0.45),
        width = 0.3,
        size = 1.23
    ) +
    geom_bar(stat = "identity",
             position = position_dodge(),
             width = 0.45) +
    ylab("Relative Concentration") +
    scale_fill_manual(values = c("black", "#818181")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5))+
    theme_classic(base_size = 14.5)+
    grids(linetype="dashed")+
    ggtitle("Enrichment SNRPN")+
    theme(plot.title = element_text(hjust = 0.5))
# dev.off()

