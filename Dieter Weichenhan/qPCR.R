library(ggplot2)
library(readxl)
library(dplyr)
library(ggpubr)

# read in excel file which contains qPCR results (provided by Supervisor)
setwd("data/")
excel <- read_excel(path = "Course_MCiPqPCR9_9_20.xls")

# select only columns that are relevant for us (Pos, Sample, Ct)
# also select rows that are relevant for us. In this case all rows containing A (standard), B(SNRPN)
# and F (spike in standard), J (spike in)

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

standard_spike$Sample <- as.double(standard_spike$Sample)


# calculate regression equation
lm_eqn_gene <- lm(standard_gene$Ct ~ standard_gene$Sample)

lm_eqn_spike <- lm(standard_spike$Ct ~ standard_spike$Sample)


# isolate metrics for the regression of the gene
metrics <- list(
    a = format(unname(coef(lm_eqn_gene)[1]), digits = 2),
    b = format(unname(coef(lm_eqn_gene)[2]), digits = 2),
    r2 = format(summary(lm_eqn_gene)$r.squared, digits = 3)
)

eq <- substitute(italic(C[T]) == a * " " * b  * italic(x) * ", " ~  ~ italic(R) ^ 2 ~ "=" ~ r2,
               metrics)

equation <- as.character(as.expression(eq))


# plot the standard (row A)

# svg("StandardCurve_SNRPN.# svg")
ggplot(standard_gene, aes(x = Sample, y = Ct)) +
    geom_smooth(method = "lm", se = T) +
    geom_point() +
    xlab("Log Concentration") +
    ylab(bquote(C[T])) +
    geom_text(
        x = -0.6,
        y = 27.21,
        label = equation,
        parse = T,
        size = 6.7
    ) +
    theme_classic(base_size = 14.5)+
    scale_y_continuous(limits = c(20, 28))+
    grids(linetype="dashed")+
    ggtitle("Standard curve SNRPN")+
    theme(plot.title = element_text(hjust = 0.5))
# dev.off()


metrics_spike <- list(
    a = format(unname(coef(lm_eqn_spike)[1]), digits = 2),
    b = format(unname(coef(lm_eqn_spike)[2]), digits = 2),
    r2 = format(summary(lm_eqn_spike)$r.squared, digits = 3)
)

eq_spike <- substitute(italic(C[T]) == a * " " * b  * italic(x) * ", " ~  ~ italic(R) ^ 2 ~ "=" ~ r2,
                 metrics_spike)

equation_spike <- as.character(as.expression(eq_spike))


# plot the standard for the spike

# svg("StandardCurve_spike.# svg")
ggplot(standard_spike, aes(x = Sample, y = Ct)) +
    geom_smooth(method = "lm", se = T) +
    geom_point() +
    xlab("Log Concentration") +
    ylab(bquote(C[T])) +
    geom_text(
        x = -0.6,
        y = 27.21,
        label = equation_spike,
        parse = T,
        size = 6.7
    ) +
    theme_classic(base_size = 14.5)+
    scale_y_continuous(limits = c(17, 29), breaks = seq(17,29, by=3))+
    grids(linetype="dashed")+
    ggtitle("Standard curve Spike DNA")+
    theme(plot.title = element_text(hjust = 0.5))

 # dev.off()






# use the regression equation to calculate the log concentration of our gene and spike in
gene["LogConcentration"] <- (gene$Ct - 20.967) / -2.945

spike["LogConcentration"] <- (spike$Ct - 18.095) / -3.237

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

# svg("SNRP.# svg")
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
    ggtitle("SNRPN")+
    theme(plot.title = element_text(hjust = 0.5))

# dev.off()

# plot the average Rel concentration for the spike

# svg("spike.# svg")
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
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.5))+
    theme_classic(base_size = 14.5)+
    grids(linetype="dashed")+
    ggtitle("Spike DNA")+
    theme(plot.title = element_text(hjust = 0.5))
# dev.off()

