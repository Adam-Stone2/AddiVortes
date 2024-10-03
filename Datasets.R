# Load necessary library if you haven't already
if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}
# Load the necessary library
library(readr)

# Import the datasets into R
abalone <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/abalone_Data.csv")
X_abalone <- as.matrix(abalone[,2:9])
Y_abalone <- as.numeric(as.matrix(abalone[,10]))

Baskball <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/Baskball_Data.csv")
X_Baskball <- as.matrix(Baskball[,2:5])
Y_Baskball <- as.numeric(as.matrix(Baskball[,6]))

Boston <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/BostonHousing_Data.csv")
X_Boston <- as.matrix(Boston[,2:14])
Y_Boston <- as.numeric(as.matrix(Boston[,15]))

Enroll <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/Enroll_Data.csv")
X_Enroll <- as.matrix(Enroll[,1:5])
Y_Enroll <- as.numeric(as.matrix(Enroll[,6]))

Fat <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/FAT_Data.csv")
X_Fat <- as.matrix(Fat[,2:14])
Y_Fat <- as.numeric(as.matrix(Fat[,15]))

HatCo <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/HatCo_Data.csv")
X_HatCo <- as.matrix(HatCo[,2:13])
Y_HatCo <- as.numeric(as.matrix(HatCo[,14]))

Labour <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/Labour_Data.csv")
X_Labour <- as.matrix(Labour[,2:20])
Y_Labour <- as.numeric(as.matrix(Labour[,21]))

Mpg <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/MPG_Data.csv")
X_Mpg <- as.matrix(Mpg[,2:8])
Y_Mpg <- as.numeric(as.matrix(Mpg[,9]))

Medicare <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/Medicare_Data.csv")
X_Medicare <- as.matrix(Medicare[,1:21])
Y_Medicare <- as.numeric(as.matrix(Medicare[,22]))

Price <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/Price_Data.csv")
X_Price <- as.matrix(Price[,2:16])
Y_Price <- as.numeric(as.matrix(Price[,17]))

Rate <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/RATE_Data.csv")
X_Rate <- as.matrix(Rate[,2:10])
Y_Rate <- as.numeric(as.matrix(Rate[,11]))

Edu <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/edu_Data.csv")
X_Edu <- as.matrix(Edu[,2:7])
Y_Edu <- as.numeric(as.matrix(Edu[,8]))

Ozone <- read_csv("https://raw.githubusercontent.com/anonymous2738/AddiVortesAlgorithm/DataSets/ozone_Data.csv")
X_Ozone <- as.matrix(Ozone[,2:9])
Y_Ozone <- as.numeric(as.matrix(Ozone[,10]))


