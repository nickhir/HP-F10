library(readxl)
library(dplyr)

# read in excel file which contains array results (provided by Supervisor)
setwd("data/")
top5 <- read_excel(path = "IslandArray.xlsx", sheet=1)

# determine the percentage of regions that are located in genes, promoters or downstream of annotated genes.

inside <- as.numeric(tally(top5, "INSIDE" == CpGLocation)/nrow(top5))
promoter <- as.numeric(tally(top5, "PROMOTER" == CpGLocation)/nrow(top5))
downstream <- as.numeric(tally(top5, "DOWNSTREAM" == CpGLocation)/nrow(top5))

# determine the number of analyzed regions:

number_DMR <- length(unique(top5$RegNo))

min2 <- read_excel(path = "IslandArray.xlsx", sheet=2)

CGI <- read_excel(path = "IslandArray.xlsx", sheet=3)
