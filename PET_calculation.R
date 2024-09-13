# Script to run Penman-Monteith Equation
# Works with functions in PMFunctions.R

# Load Libraries
library(xts)

# Set working Directory
setwd('C:\\Users\\akafle1\\OneDrive - Lamar University\\Aalok_Kushum\\PET')
source('PM_Function.R')


# Read Station Data
lat <- 34.3898 # Station Latitude (decimal degrees)
z <- 1310  # Station Elevation (m)
alpha <- 0.23  # Default Surface albedo

# Read Datasets
a <- read.csv('CompleteData_PET.csv', header=TRUE)
bgn <- as.Date('2024-07-09')
end <- as.Date('2024-08-16')
bgnx <- as.Date('2024-07-09')
datex <- seq.Date(bgn,end,'day')
dates <- seq.Date(bgnx,end,'day')
aa <- data.frame(datex,a) # join date
aa <- subset(aa,datex > as.Date('2024-07-08'))
Tempmax <- aa$Temp_Max
Tempmin <- aa$Temp_Min
Tempdew <- aa$DewPtAvg
Tmax <- (5/9)*(Tempmax-32)    #temperature from F to C
Tmin <- (5/9)*(Tempmin-32)
Tdew <- (5/9)*(Tempdew-32)
u2 <- aa$WS2M  # Windspeed at 2m height
Rs <- aa$ALLSKY_SFC_SW_DWN  # Solar Radiation in Mj/m^2/d


# Calculate ET Solar Radiation MJ/m^2/d
Rax <- Ra(lat,dates)

# Compute Solar Radiation using Hargreaves-Samani Equation
Rsh <- Rshar(Tmax,Tmin,Rax,1)  #1 implies inland location

# Calculate saturated vapor pressure
es <- svp(Tmax,Tmin)

# Calculate actual vapor pressure
ea <- eadew(Tdew) # Using Dewpoint Data

# Calculate vapor pressure deficit
vpdx <- vpd(es,ea)

# Compute Penman-Monteith using
ETHS <- mat.or.vec(length(Tmax),3) #...solar radiation estimation from tmax-tmin
ETPW <- mat.or.vec(length(Tmax),3) #...actual solar radiation

for(i in seq(1,length(Tmax),1))
{
  ETHS[i,1:3] <- pm(Tmax[i],Tmin[i],u2[i],Rsh[i],Rax[i],ea[i],lat,z,alpha) # Based on solar radiation estimation from tmax-tmin
  ETPW[i,1:3] <- pm(Tmax[i],Tmin[i],u2[i],Rs[i],Rax[i],ea[i],lat,z,alpha)  # Using actual solar radiation
}
# Attach Date object
ETHS <- data.frame(dates,ETHS)
ETPW <- data.frame(dates,ETPW)


# Create time-series and aggregate
ETHSRAD.day <- xts(ETHS[,2],dates,1)
ETHSWND.day <- xts(ETHS[,3],dates,1)
ETHSTOT.day <- xts(ETHS[,4],dates,1)


ETHSRAD.mon <- apply.monthly(ETHSRAD.day,'sum')
ETHSWND.mon <- apply.monthly(ETHSWND.day,'sum')
ETHSTOT.mon <- apply.monthly(ETHSTOT.day,'sum')

ETHSRAD.yrs <- apply.yearly(ETHSRAD.day,'sum')
ETHSWND.yrs <- apply.yearly(ETHSWND.day,'sum')
ETHSTOT.yrs <- apply.yearly(ETHSTOT.day,'sum')

#Comparing yearly ET values based on HS and actual solar radiation
ETPWTOT.day <- xts(ETPW[,4],dates,1)
ETPWTOT.yrs <- apply.yearly(ETPWTOT.day,'sum')
compareHSPW <- cbind(ETHSTOT.yrs,ETPWTOT.yrs)

#Combining the evapotranspiration calculated
# Remove the date column from ETPW
ETPW <- ETPW[, -1]
combined_df <- cbind(ETHS, ETPW)

#Rename the columns
colnames(combined_df) <- c("Date", "ETHS_solar_mm/day", "ETHS_wind_mm/day", "ETHS_Total_mm/day", "ETPW_solar_mm/day", "ETPW_wind_mm/day", "ETPW_Total_mm/day")

new_df <- combined_df[, c("Date","ETHS_wind_mm/day", "ETPW_solar_mm/day", "ETPW_Total_mm/day")]
colnames(new_df) <- c("Date", "ETWind", "ETSolar", "ETTotal")

# Load the ggplot2 library
library(ggplot2)

# Create a time series plot with increased Y-axis label size and other modifications
ggplot(new_df, aes(x = as.Date(Date))) +
  geom_line(aes(y = ETWind, color = "ETWind"), size = 1) +
  geom_line(aes(y = ETSolar, color = "ETSolar"), size = 1) +
  geom_line(aes(y = ETTotal, color = "ETTotal"), size = 1) +
  labs(title = "Time Series Plot of ETWind, ETSolar, and ETTotal",
       x = "Date",
       y = "Evapotranspiration (mm/day)",
       color = "Legend") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add box around the plot
        plot.title = element_text(hjust = 0.5),  # Center the title
        axis.text.y = element_text(size = 14),   # Increase Y-axis label size
        axis.title.y = element_text(size = 16))  # Increase Y-axis title size

