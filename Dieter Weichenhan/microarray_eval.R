library(readxl)
library(dplyr)

# read in excel file which contains array results (provided by Supervisor)
setwd("data/")
top5 <- read_excel(path = "IslandArray.xlsx", sheet=1)
min2 <- read_excel(path = "IslandArray.xlsx", sheet=3)
island <- read_excel(path = "IslandArray.xlsx", sheet=4)
promoter <- read_excel(path = "IslandArray.xlsx", sheet=6)

# determine the percentage of regions that are located in genes, promoters or downstream of annotated genes.

inside <- as.numeric(tally(island, "INSIDE" == CpGLocation)/nrow(island))
promoter <- as.numeric(tally(island, "PROMOTER" == CpGLocation)/nrow(island))
downstream <- as.numeric(tally(island, "DOWNSTREAM" == CpGLocation)/nrow(island))
unknown <- as.numeric(tally(island, "Unknown" == CpGLocation)/nrow(island))
DIVERGENT_PROMOTER <- as.numeric(tally(island, "DIVERGENT_PROMOTER" == CpGLocation)/nrow(island))

# load in pathway enrichment analysis results

kegg <- read.table("KEGG_Chart.txt", sep = "\t", header=TRUE)
panther <- read.table("PANTHER_chart.txt", sep = "\t", header=TRUE)
