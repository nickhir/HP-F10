library(ggplot2)
library(readxl)
library(dplyr)


# read in excel file which contains qPCR results (provided by Supervisor)
setwd("C:/Users/Nick/Desktop/Studium/Master/2. Semester/HP-F10/5_Dieter Weichenhan/")
excel <- read_excel(path = "Course_MCiPqPCR9_9_20.xls")

# select only columns that are relevant for us (Pos, Sample, Ct)
# also select rows that are relevant for us. In this case all rows containing A (standard), B(SNRPN) 
# and H (spike in)

standard <- excel %>%
    select(c("Pos", "Sample", "Ct")) %>%
    filter(grepl("A\\d+", Pos))

gene <- excel %>%
    select(c("Pos", "Sample", "Ct")) %>%
    filter(grepl("B\\d+", Pos))

spike <- excel %>%
    select(c("Pos", "Sample", "Ct")) %>%
    filter(grepl("H\\d+", Pos))

# for standard we have to convert column Sample to double instead of chr
standard$Sample <- as.double(standard$Sample)


# calculate regression equation
lm_eqn <- lm(standard$Ct ~ standard$Sample)

metrics <- list(a = format(unname(coef(lm_eqn)[1]), digits = 2),
                b = format(unname(coef(lm_eqn)[2]), digits = 2),
                r2 = format(summary(lm_eqn)$r.squared, digits = 3))

eq <- substitute(italic(C[T]) == a * " "* b  *italic(x)*", "~~italic(R)^2~"="~r2, metrics)
equation <- as.character(as.expression(eq))


# plot the standard (row A)

#svg("StandardCurve.svg")
ggplot(standard, aes(x=Sample, y=Ct))+
    geom_smooth(method = "lm", se = T)+
    geom_point()+
    xlab("Log Concentration")+
    ylab(bquote(C[T]))+
    geom_text(x = -0.6, y=26.21, label=equation, parse = T, size=6.7)
#dev.off()

# use the regression equation (CT-20.967/-2.945 = Conc) to calculate the log concentration of our gene and spike in
gene["LogConcentration"] <- (gene$Ct-20.967)/-2.945

spike["LogConcentration"] <- (spike$Ct-20.967)/-2.945

# calculate relative concentrattion

gene["RelativeConcentration"] <- 10^gene$LogConcentration

spike["RelativeConcentration"] <- 10^spike$LogConcentration

# calculate average for technical duplicates

gene <- gene %>% 
            group_by(Sample) %>%
            mutate(AverageRelConc = mean(RelativeConcentration, na.rm = T))

spike <- spike %>%
    group_by(Sample) %>%
    mutate(AverageRelConc = mean(RelativeConcentration, na.rm = T))



# remove genotype from sample name and move it to new column
gene$Sample <- sub("wt\\d+|dko\\d+", "", gene$Sample)
gene["Genotype"] <- c(rep("WT",10), rep("DKO",10))

spike$Sample <- sub("wt\\d+|dko\\d+", "", spike$Sample)
spike["Genotype"] <- c(rep("WT",10), rep("DKO",10))

# plot the average Rel concentration for the gene
ggplot(gene, aes(x=Sample, y=AverageRelConc, fill=Genotype))+
    geom_bar(stat="identity", position = "dodge")+
    ylab("Relative Concentration")+
    scale_fill_grey(start=0.23, end=0.67)

# plot the average Rel concentration for the spike
ggplot(spike, aes(x=Sample, y=AverageRelConc, fill=Genotype))+
    geom_bar(stat="identity", position = "dodge")+
    ylab("Relative Concentration")+
    scale_fill_grey(start=0.23, end=0.67)

