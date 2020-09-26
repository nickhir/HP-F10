library(ggplot2)
library(dplyr)

# data was taken from web of science by searching for topic "Epigenetics"
# read in data
df <- read.table("citation.txt", sep = "\t", header=TRUE)

# we are only interested in years greater than 2004 (pretty random cutoff)
df <- df %>%
    filter(Year>2004 & Year !=2020) %>%
    select(c("Year", "records"))

# plot data

#svg("citations_epigenetics.svg")
ggplot(df, aes(x=Year, y=records))+
    geom_bar(stat="identity", fill="#333333")+
    scale_y_continuous(breaks = seq(0,2500, by=500))+
    scale_x_continuous(breaks= seq(2005,2019, by=2))+
    ylab("Number of citations")+
    xlab("Year")
#dev.off()
