# required packages
library(ggplot2)
library(reshape2)
library(drc)
library(ggpubr)

# the path might be adjusted and should point to the indicated directory (they are included in the github repository).
setwd("./digestion-analysis")

# create a function, which reads in the txt file from the supervisor and calculates the necessary metrics. Afterwards return a data frame with the information
read.data <- function(file) {
    df <- read.table(file, sep="\t", header = TRUE)
    
    # select the desired data. We are only interested in the optical (density-background)/mm^2 which is the last column. Last two rows are also not needed
    od <- df[1:16, 10]
    
    # first 8 values correspond to the upper band. Last 8 to the lower band
    upper <- od[1:8]
    lower <- od[9:16]
    
    # creat a dataframe which annotates the values
    output.dataframe <- data.frame(
        upper_band = upper,
        lower_band = lower,
        row.names = c("positive", "negative", "SAH", "50", "25", "12.5", "6.25", "3.125")
    )
    
    # Per lane, calculate the percentage for each band for the respective lane
    output.dataframe <- cbind(output.dataframe, 
                              upper_band_percentage = upper/(upper+lower), 
                              lower_band_percentage = lower/(upper+lower))
    
    
    # For the negative control we dont expect to see a upper band, because everything should be unmethylated and thus cleaved. 
    # Hence, we assume, that the signal we get for "negative" condition in the upper band is noise/background.
    # This is why in the next step the substract this value from each upper band
    
    output.dataframe <- cbind(output.dataframe, 
                              upper_band_percentage_adjusted = output.dataframe$upper_band_percentage-output.dataframe$upper_band_percentage[2])
    
    
    # For the positive control we dont expect to see a lower band, because everything should be methylated. 
    # Hence, we assume, that the signal we get for "positive" condition in the lower band is noise/background.
    # This is why in the next step the substract this value from each lower band
    
    output.dataframe <- cbind(output.dataframe, 
                              lower_band_percentage_adjusted = output.dataframe$lower_band_percentage-output.dataframe$lower_band_percentage[1])
    
    # Lastly, the upperband is normalized with regard to the adjusted "positive" condition and 
    # the lowerband is normalized with regard to the adjusted "negative" condition
    
    output.dataframe <- cbind(output.dataframe, 
                              normalized_upper = output.dataframe$upper_band_percentage_adjusted/output.dataframe$upper_band_percentage_adjusted[1],
                              normalized_lower = output.dataframe$lower_band_percentage_adjusted/output.dataframe$lower_band_percentage_adjusted[2])
    
    # for some experminets we observe percentages above 100 or below 0. These make biologically no sense and are set to 100 or 0 respectively.
    # identify position of value in normalized upper  which is greater than 1. If one exists, set the value to one. Other conditions are very similar.
    output.dataframe$normalized_upper[which(output.dataframe$normalized_upper > 1)] = 1
    output.dataframe$normalized_upper[which(output.dataframe$normalized_upper < 0)] = 0
    
    output.dataframe$normalized_lower[which(output.dataframe$normalized_lower > 1)] = 1
    output.dataframe$normalized_lower[which(output.dataframe$normalized_lower < 0)] = 0
    
    # At the very end perform a quick sanity check: 
    # Adding normalized_upper and normalized_lower should equal 1.
    ifelse(round(output.dataframe$normalized_upper + output.dataframe$normalized_lower, digits=6) == 1, 
           print("Loading successful"), 
           print("!!!ERROR!!! SOMETHING WENT WRONG WHILE LOADING AND MANIPULATING THE INPUT DATA.\n
                 DO NOT CONTINUE THIS WORKFLOW UNTIL LOADING IS SUCCESSFUL"))
    
    
    return(output.dataframe)
}


# read in the data 
# use these three samlpes because we were working at the same bench
nh <- read.data("nh.txt")
he <- read.data("he.txt")
jh <- read.data("jh.txt")


# for our further analysis we are only interested in the normalized upper band of every person 
# and the different concentrations of the inhibitor we used.

analysis <- data.frame(
    dose = c(50, 25, 12.5, 6.25, 3.125),
    nh = nh$normalized_upper[4:8],
    he = he$normalized_upper[4:8],
    jh = jh$normalized_upper[4:8],
)


# we want values between 0 and 1 for our responds

# calculate average
analysis <- cbind(analysis, mean = rowMeans(analysis[2:4]))

# calculate logistic model function
model <- drm(mean ~ dose, data = analysis, fct = LL.4())

# we want to plot the model with ggplot
# no direct way to isolate function of curve. Workaround is to predict many values that fall
# into our dose range

DummyDoses <- expand.grid(conc = seq(3, 51, length = 400))
pm <- predict(model, newdata = DummyDoses)

fitted_curve <- cbind(DummyDoses, pm)

# for ggplot to work properly dataframe has to be in the _long_ formt
# also generate data frame wiht average and sd
long_df <- melt(analysis[1:4], id.vars = "dose")

average <- data.frame(
    dose = c(50, 25, 12.5, 6.25, 3.125),
    average = rowMeans(analysis[2:4]),
    sd = matrixStats::rowSds(as.matrix(analysis[2:4]))
)


#svg("dose-responds-curve-all.svg")
ggplot() +
    # plot error bars
    #geom_errorbar(data=average, aes(x=dose, ymin=average-sd, ymax=average+sd), width=1.1, size=1, alpha=1)+
    # plot the curve
    geom_line(data = fitted_curve, aes(x = conc, y = pm), size = 1.25, color = "black") +
    # plot the actual data
    geom_point(data = long_df, aes(x = dose, y = value), shape = 18, size = 2.5, color = "black") +
    # plot the average as a bar
    geom_point(
        data = average, aes(x = dose, y = average), shape = 95, size = 10, color = "#C60B0B") +
    xlab("DD880 concentration [µM]") +
    ylab("M.Sss1 activity")+
    scale_x_continuous(breaks = seq(0, 55, by = 5))+
    theme_classic(base_size = 14.5)+
    grids(linetype="dashed")

#dev.off()

sessionInfo()

