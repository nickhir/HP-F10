library(readxl)
library(dplyr)
library(ggplot2)


# read in data
data <- read_excel("./data/qPCR results_ct values (Practical_2020).xlsx", skip=1)

# we are only interested in the position and the Ct value
data <- data %>%
    select(c("Pos", "Cp"))

# further select only the positions that we are interested in. For me: WT+4OHT vs ID3-KO+4OHT for both input and H4K16ac

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


# select the desired values for the H4K16ac ----
WT.H4K16ac.DSB3 <- data %>%
    filter(Pos %in% c("J17", "J18", "J19")) %>%
    select("Cp")

KO.H4K16ac.DSB3 <- data %>%
    filter(Pos %in% c("J21", "J22", "J23")) %>%
    select("Cp")

WT.H4K16ac.DSB4 <- data %>%
    filter(Pos %in% c("K17", "K18", "K19")) %>%
    select("Cp")

KO.H4K16ac.DSB4 <- data %>%
    filter(Pos %in% c("K21", "K22", "K23")) %>%
    select("Cp")

WT.H4K16ac.DSBIII <- data %>%
    filter(Pos %in% c("L17", "L18", "L19")) %>%
    select("Cp")

KO.H4K16ac.DSBIII <- data %>%
    filter(Pos %in% c("L21", "L22", "L23")) %>%
    select("Cp")

WT.H4K16ac.DSBIV <- data %>%
    filter(Pos %in% c("M17", "M18", "M19")) %>%
    select("Cp")

KO.H4K16ac.DSBIV <- data %>%
    filter(Pos %in% c("M21", "M22", "M23")) %>%
    select("Cp")

WT.H4K16ac.DSBV <- data %>%
    filter(Pos %in% c("N17", "N18", "N19")) %>%
    select("Cp")

KO.H4K16ac.DSBV <- data %>%
    filter(Pos %in% c("N21", "N22", "N23")) %>%
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

# caculate the enrichment for H4K16ac for each individual well with the following equation 100*2^-Ct sample /2^-Ct input. 
# For Ct input use the above calculated average ----

WT.H4K16ac.DSB3 <- WT.H4K16ac.DSB3 %>%
    mutate(enrichment = 100*(2^-WT.H4K16ac.DSB3$Cp/2^-WT.input.DSB3.corrected)) 

WT.H4K16ac.DSB4 <- WT.H4K16ac.DSB4 %>%
    mutate(enrichment = 100*(2^-WT.H4K16ac.DSB4$Cp/2^-WT.input.DSB4.corrected))

WT.H4K16ac.DSBIII <- WT.H4K16ac.DSBIII %>%
    mutate(enrichment = 100*(2^-WT.H4K16ac.DSBIII$Cp/2^-WT.input.DSBIII.corrected))

WT.H4K16ac.DSBIV <- WT.H4K16ac.DSBIV %>%
    mutate(enrichment = 100*(2^-WT.H4K16ac.DSBIV$Cp/2^-WT.input.DSBIV.corrected))

WT.H4K16ac.DSBV <- WT.H4K16ac.DSBV %>%
    mutate(enrichment = 100*(2^-WT.H4K16ac.DSBV$Cp/2^-WT.input.DSBV.corrected))


KO.H4K16ac.DSB3 <- KO.H4K16ac.DSB3 %>%
    mutate(enrichment = 100*(2^-KO.H4K16ac.DSB3$Cp/2^-KO.input.DSB3.corrected))

KO.H4K16ac.DSB4 <- KO.H4K16ac.DSB4 %>%
    mutate(enrichment = 100*(2^-KO.H4K16ac.DSB4$Cp/2^-KO.input.DSB4.corrected))

KO.H4K16ac.DSBIII <- KO.H4K16ac.DSBIII %>%
    mutate(enrichment = 100*(2^-KO.H4K16ac.DSBIII$Cp/2^-KO.input.DSBIII.corrected))

KO.H4K16ac.DSBIV <- KO.H4K16ac.DSBIV %>%
    mutate(enrichment = 100*(2^-KO.H4K16ac.DSBIV$Cp/2^-KO.input.DSBIV.corrected))

KO.H4K16ac.DSBV <- KO.H4K16ac.DSBV %>%
    mutate(enrichment = 100*(2^-KO.H4K16ac.DSBV$Cp/2^-KO.input.DSBV.corrected))






# average the calculated enrichments for each primer ----

WT.H4K16ac.DSB3.enrichment <- mean(WT.H4K16ac.DSB3$enrichment)
WT.H4K16ac.DSB4.enrichment <- mean(WT.H4K16ac.DSB4$enrichment)
WT.H4K16ac.DSBIII.enrichment <- mean(WT.H4K16ac.DSBIII$enrichment)
WT.H4K16ac.DSBIV.enrichment <- mean(WT.H4K16ac.DSBIV$enrichment)
WT.H4K16ac.DSBV.enrichment <- mean(WT.H4K16ac.DSBV$enrichment)

KO.H4K16ac.DSB3.enrichment <- mean(KO.H4K16ac.DSB3$enrichment)
KO.H4K16ac.DSB4.enrichment <- mean(KO.H4K16ac.DSB4$enrichment)
KO.H4K16ac.DSBIII.enrichment <- mean(KO.H4K16ac.DSBIII$enrichment)
KO.H4K16ac.DSBIV.enrichment <- mean(KO.H4K16ac.DSBIV$enrichment)
KO.H4K16ac.DSBV.enrichment <- mean(KO.H4K16ac.DSBV$enrichment)

# calculated enrichments standard derivation for each primer ----
WT.H4K16ac.DSB3.sd <- sd(WT.H4K16ac.DSB3$enrichment)
WT.H4K16ac.DSB4.sd <- sd(WT.H4K16ac.DSB4$enrichment)
WT.H4K16ac.DSBIII.sd <- sd(WT.H4K16ac.DSBIII$enrichment)
WT.H4K16ac.DSBIV.sd <- sd(WT.H4K16ac.DSBIV$enrichment)
WT.H4K16ac.DSBV.sd <- sd(WT.H4K16ac.DSBV$enrichment)

KO.H4K16ac.DSB3.sd <- sd(KO.H4K16ac.DSB3$enrichment)
KO.H4K16ac.DSB4.sd <- sd(KO.H4K16ac.DSB4$enrichment)
KO.H4K16ac.DSBIII.sd <- sd(KO.H4K16ac.DSBIII$enrichment)
KO.H4K16ac.DSBIV.sd <- sd(KO.H4K16ac.DSBIV$enrichment)
KO.H4K16ac.DSBV.sd <- sd(KO.H4K16ac.DSBV$enrichment)
 

# generate the data frame which will be used for plotting ---- 

overview <- data.frame(
    Genotype = c(rep("WT",5), rep("KO",5)),
    primer = rep(c("DSB-3", "DSB-4", "DSB-III", "DSB-IV", "DSB-V"),2),
    # now the actual enrichments. First WT 3,4,III,IV,V then  KO 3,4,III,IV,V
    enrichment = c(WT.H4K16ac.DSB3.enrichment, WT.H4K16ac.DSB4.enrichment, WT.H4K16ac.DSBIII.enrichment, WT.H4K16ac.DSBIV.enrichment, WT.H4K16ac.DSBV.enrichment,
                   KO.H4K16ac.DSB3.enrichment, KO.H4K16ac.DSB4.enrichment, KO.H4K16ac.DSBIII.enrichment, KO.H4K16ac.DSBIV.enrichment, KO.H4K16ac.DSBV.enrichment),
    sd = c(WT.H4K16ac.DSB3.sd, WT.H4K16ac.DSB4.sd, WT.H4K16ac.DSBIII.sd, WT.H4K16ac.DSBIV.sd, WT.H4K16ac.DSBV.sd,
           KO.H4K16ac.DSB3.sd, KO.H4K16ac.DSB4.sd, KO.H4K16ac.DSBIII.sd, KO.H4K16ac.DSBIV.sd, KO.H4K16ac.DSBV.sd)
)

# specify order of the plot
overview$Genotype <- factor(overview$Genotype, levels = c("WT", "KO"))


# plot the dataframe ---- 

# svg("qRT PCR H4.svg")
ggplot(data = overview, aes(x=primer, y=enrichment, fill=Genotype))+
    geom_errorbar(aes(ymin = enrichment-sd, ymax=enrichment+sd),
                  position=position_dodge(0.45), width=0.3, size=1.23)+
    geom_bar(stat="identity", position=position_dodge(), width=0.45)+
    ylab("Enrichment of input [%]")+
    xlab("")+
    scale_fill_manual(values = c("black", "#919191"))+
    theme_classic(base_size = 14.5)+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 6), breaks = seq(0,6,by=1))
# dev.off()
