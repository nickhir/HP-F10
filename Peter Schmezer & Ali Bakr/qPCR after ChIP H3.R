library(readxl)
library(dplyr)
library(ggplot2)


# read in data
data <- read_excel("./data/qPCR results_ct values (Practical_2020).xlsx", skip=1)

# we are only interested in the position and the Ct value
data <- data %>%
    select(c("Pos", "Cp"))

# further select only the positions that we are interested in.

# select the desired values for the 10% input ----
WT.input.DSB3 <- data %>%
    filter(Pos %in% c("J1", "J2", "J3")) %>%
    select("Cp") 

KO.input.DSB3 <- data %>%
    filter(Pos %in% c("J5", "J6", "J7")) %>%
    select("Cp")

WT.input.DSB4 <- data %>%
    filter(Pos %in% c("K1", "K2", "K3")) %>%
    select("Cp")

KO.input.DSB4 <- data %>%
    filter(Pos %in% c("K5", "K6", "K7")) %>%
    select("Cp")

WT.input.DSBIII <- data %>%
    filter(Pos %in% c("L1", "L2", "L3")) %>%
    select("Cp")

KO.input.DSBIII <- data %>%
    filter(Pos %in% c("L5", "L6", "L7")) %>%
    select("Cp")

WT.input.DSBIV <- data %>%
    filter(Pos %in% c("M1", "M2", "M3")) %>%
    select("Cp")

KO.input.DSBIV <- data %>%
    filter(Pos %in% c("M5", "M6", "M7")) %>%
    select("Cp")

WT.input.DSBV <- data %>%
    filter(Pos %in% c("N1", "N2", "N3")) %>%
    select("Cp")

KO.input.DSBV <- data %>%
    filter(Pos %in% c("N5", "N6", "N7")) %>%
    select("Cp")


# select the desired values for the h3me36me3 ----
WT.h3me36me3.DSB3 <- data %>%
    filter(Pos %in% c("J9", "J10", "J11")) %>%
    select("Cp")

KO.h3me36me3.DSB3 <- data %>%
    filter(Pos %in% c("J13", "J14", "J15")) %>%
    select("Cp")

WT.h3me36me3.DSB4 <- data %>%
    filter(Pos %in% c("K9", "K10", "K11")) %>%
    select("Cp")

KO.h3me36me3.DSB4 <- data %>%
    filter(Pos %in% c("K13", "K14", "K15")) %>%
    select("Cp")

WT.h3me36me3.DSBIII <- data %>%
    filter(Pos %in% c("L9", "L10", "L11")) %>%
    select("Cp")

KO.h3me36me3.DSBIII <- data %>%
    filter(Pos %in% c("L13", "L14", "L15")) %>%
    select("Cp")

WT.h3me36me3.DSBIV <- data %>%
    filter(Pos %in% c("M9", "M10", "M11")) %>%
    select("Cp")

KO.h3me36me3.DSBIV <- data %>%
    filter(Pos %in% c("M13", "M14", "M15")) %>%
    select("Cp")

WT.h3me36me3.DSBV <- data %>%
    filter(Pos %in% c("N9", "N10", "N11")) %>%
    select("Cp")

KO.h3me36me3.DSBV <- data %>%
    filter(Pos %in% c("N13", "N14", "N15")) %>%
    select("Cp")





# for the input sample, calculate the average and substract 3.322 (number of cycles which would be needed to reach 100%) ----

WT.input.DSB3.corrected <- mean(WT.input.DSB3$Cp)-3.322
WT.input.DSB4.corrected <- mean(WT.input.DSB4$Cp)-3.322
WT.input.DSBIII.corrected <- mean(WT.input.DSBIII$Cp)-3.322
WT.input.DSBIV.corrected <- mean(WT.input.DSBIV$Cp)-3.322
WT.input.DSBV.corrected <- mean(WT.input.DSBV$Cp)-3.322

KO.input.DSB3.corrected <- mean(KO.input.DSB3$Cp)-3.322
KO.input.DSB4.corrected <- mean(KO.input.DSB4$Cp)-3.322
KO.input.DSBIII.corrected <- mean(KO.input.DSBIII$Cp)-3.322
KO.input.DSBIV.corrected <- mean(KO.input.DSBIV$Cp)-3.322
KO.input.DSBV.corrected <- mean(KO.input.DSBV$Cp)-3.322

# caculate the enrichment for h3me36me3 for each individual well with the following equation 100*2^-Ct sample /2^-Ct input. 
# For Ct input use the above calculated average ----

WT.h3me36me3.DSB3 <- WT.h3me36me3.DSB3 %>%
    mutate(enrichment = 100*(2^-WT.h3me36me3.DSB3$Cp/2^-WT.input.DSB3.corrected)) 

WT.h3me36me3.DSB4 <- WT.h3me36me3.DSB4 %>%
    mutate(enrichment = 100*(2^-WT.h3me36me3.DSB4$Cp/2^-WT.input.DSB4.corrected))

WT.h3me36me3.DSBIII <- WT.h3me36me3.DSBIII %>%
    mutate(enrichment = 100*(2^-WT.h3me36me3.DSBIII$Cp/2^-WT.input.DSBIII.corrected))

WT.h3me36me3.DSBIV <- WT.h3me36me3.DSBIV %>%
    mutate(enrichment = 100*(2^-WT.h3me36me3.DSBIV$Cp/2^-WT.input.DSBIV.corrected))

WT.h3me36me3.DSBV <- WT.h3me36me3.DSBV %>%
    mutate(enrichment = 100*(2^-WT.h3me36me3.DSBV$Cp/2^-WT.input.DSBV.corrected))


KO.h3me36me3.DSB3 <- KO.h3me36me3.DSB3 %>%
    mutate(enrichment = 100*(2^-KO.h3me36me3.DSB3$Cp/2^-KO.input.DSB3.corrected))

KO.h3me36me3.DSB4 <- KO.h3me36me3.DSB4 %>%
    mutate(enrichment = 100*(2^-KO.h3me36me3.DSB4$Cp/2^-KO.input.DSB4.corrected))

KO.h3me36me3.DSBIII <- KO.h3me36me3.DSBIII %>%
    mutate(enrichment = 100*(2^-KO.h3me36me3.DSBIII$Cp/2^-KO.input.DSBIII.corrected))

KO.h3me36me3.DSBIV <- KO.h3me36me3.DSBIV %>%
    mutate(enrichment = 100*(2^-KO.h3me36me3.DSBIV$Cp/2^-KO.input.DSBIV.corrected))

KO.h3me36me3.DSBV <- KO.h3me36me3.DSBV %>%
    mutate(enrichment = 100*(2^-KO.h3me36me3.DSBV$Cp/2^-KO.input.DSBV.corrected))






# average the calculated enrichments for each primer ----

WT.h3me36me3.DSB3.enrichment <- mean(WT.h3me36me3.DSB3$enrichment)
WT.h3me36me3.DSB4.enrichment <- mean(WT.h3me36me3.DSB4$enrichment)
WT.h3me36me3.DSBIII.enrichment <- mean(WT.h3me36me3.DSBIII$enrichment)
WT.h3me36me3.DSBIV.enrichment <- mean(WT.h3me36me3.DSBIV$enrichment)
WT.h3me36me3.DSBV.enrichment <- mean(WT.h3me36me3.DSBV$enrichment)

KO.h3me36me3.DSB3.enrichment <- mean(KO.h3me36me3.DSB3$enrichment)
KO.h3me36me3.DSB4.enrichment <- mean(KO.h3me36me3.DSB4$enrichment)
KO.h3me36me3.DSBIII.enrichment <- mean(KO.h3me36me3.DSBIII$enrichment)
KO.h3me36me3.DSBIV.enrichment <- mean(KO.h3me36me3.DSBIV$enrichment)
KO.h3me36me3.DSBV.enrichment <- mean(KO.h3me36me3.DSBV$enrichment)

# calculated enrichments standard derivation for each primer ----
WT.h3me36me3.DSB3.sd <- sd(WT.h3me36me3.DSB3$enrichment)
WT.h3me36me3.DSB4.sd <- sd(WT.h3me36me3.DSB4$enrichment)
WT.h3me36me3.DSBIII.sd <- sd(WT.h3me36me3.DSBIII$enrichment)
WT.h3me36me3.DSBIV.sd <- sd(WT.h3me36me3.DSBIV$enrichment)
WT.h3me36me3.DSBV.sd <- sd(WT.h3me36me3.DSBV$enrichment)

KO.h3me36me3.DSB3.sd <- sd(KO.h3me36me3.DSB3$enrichment)
KO.h3me36me3.DSB4.sd <- sd(KO.h3me36me3.DSB4$enrichment)
KO.h3me36me3.DSBIII.sd <- sd(KO.h3me36me3.DSBIII$enrichment)
KO.h3me36me3.DSBIV.sd <- sd(KO.h3me36me3.DSBIV$enrichment)
KO.h3me36me3.DSBV.sd <- sd(KO.h3me36me3.DSBV$enrichment)
 

# generate the data frame which will be used for plotting ---- 

overview <- data.frame(
    Genotype = c(rep("WT",5), rep("KO",5)),
    primer = rep(c("DSB-3", "DSB-4", "DSB-III", "DSB-IV", "DSB-V"),2),
    # now the actual enrichments. First WT 3,4,III,IV,V then  KO 3,4,III,IV,V
    enrichment = c(WT.h3me36me3.DSB3.enrichment, WT.h3me36me3.DSB4.enrichment, WT.h3me36me3.DSBIII.enrichment, WT.h3me36me3.DSBIV.enrichment, WT.h3me36me3.DSBV.enrichment,
                   KO.h3me36me3.DSB3.enrichment, KO.h3me36me3.DSB4.enrichment, KO.h3me36me3.DSBIII.enrichment, KO.h3me36me3.DSBIV.enrichment, KO.h3me36me3.DSBV.enrichment),
    sd = c(WT.h3me36me3.DSB3.sd, WT.h3me36me3.DSB4.sd, WT.h3me36me3.DSBIII.sd, WT.h3me36me3.DSBIV.sd, WT.h3me36me3.DSBV.sd,
           KO.h3me36me3.DSB3.sd, KO.h3me36me3.DSB4.sd, KO.h3me36me3.DSBIII.sd, KO.h3me36me3.DSBIV.sd, KO.h3me36me3.DSBV.sd)
)

# specify order of the plot
overview$Genotype <- factor(overview$Genotype, levels = c("WT", "KO"))


# plot the dataframe ---- 

# svg("qRT PCR H3.svg")
ggplot(data = overview, aes(x=primer, y=enrichment, fill=Genotype))+
    geom_errorbar(aes(ymin = enrichment-sd, ymax=enrichment+sd),
                  position=position_dodge(0.45), width=0.3, size=1.23)+
    geom_bar(stat="identity", position=position_dodge(), width=0.45)+
    ylab("Enrichment of input [%]")+
    xlab("")+
    scale_fill_manual(values = c("black", "#919191"))+
    theme_classic(base_size = 14.5)+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 37), breaks = seq(0,35,by=5))
# dev.off()
