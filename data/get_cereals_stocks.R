### This script uses USDA-PSD data to calculate stocks for cereals 
###  between 1986 and 2011, and outputs a (countries x years) matrix 
# Author: Philippe Marchand (pmarchand@sesync.org)

library(tidyr)
library(dplyr)
library(stringr)

## (1) Load data and get reserves in kcal

# USDA data (download from http://apps.fas.usda.gov/psdonline/psdDownload.aspx)
psd <- read.csv("data/psd_grains_pulses.csv")

# Extract year-end stocks 1986-2011 (units: 1000 metric tons)
psd8611 <- select(psd, Country = Country_Name, Year = Market_Year,
                  Product = Commodity_Description,
                  Variable = Attribute_Description, Value) %>%
    filter(Year %in% 1986:2011) %>%
    spread(key = Variable, value = Value) %>%
    select(Country, Year, Product, R = `Ending Stocks`)

# Read table of kcal/tonnes conversion factors and match commodity names
conv_kcal <- read.delim("data/crop_list.tsv", header = TRUE, as.is = TRUE)
subs_list <- c("Rice Milled" = "Rice, Milled", "Maize" = "Corn", 
               "Mixed grain" = "Mixed Grain")
conv_kcal$cropname <- str_replace_all(conv_kcal$cropname, subs_list)

# Calculate R in kcal, reshape to Country x Year
Rkbyc <- inner_join(psd8611, conv_kcal, by = c("Product" = "cropname")) %>%
         select(Country, Year, Product, R, convf = kcal.ton) %>%
         mutate(R = R * 1000 * as.numeric(convf)) %>%
         group_by(Country, Year) %>%
         summarise(R = sum(R, na.rm = TRUE)) %>%
         spread(key = Year, value = R, fill = 0)


## (2) Various fixes to ensure countries match between FAOSTAT and PSD data

# Fix non-matching country names
ciso <- read.table("ciso3.txt", stringsAsFactors = FALSE)
colnames(ciso) <- c("FAO", "iso3", "name")
ciso$name[56] <- "Cote d'Ivoire" 
ciso$name[185] <- "Reunion"
countries_to_match <- unique(Rkbyc$Country[!(Rkbyc$Country %in% ciso$name)])
country_repl <- c("Bolivia" = "Bolivia (Plurinational State of)", "Brunei" = 
                  "Brunei Darussalam", "Burkina" =  "Burkina Faso", "Burma" = 
                  "Myanmar", "Former Czechoslovakia" = "Czechoslovakia", 
                  "Former Yugoslavia" = "Yugoslav SFR", "Gambia, The" = "Gambia",
                  "Iran" = "Iran (Islamic Republic of)", "Korea, North" = 
                  "Democratic People's Republic of Korea", "Korea, South" = 
                  "Republic of Korea", "Laos" = "Lao People's Democratic Republic",
                  "Macedonia" = "The former Yugoslav Republic of Macedonia",
                  "Moldova" = "Republic of Moldova", "Russia" = 
                  "Russian Federation", "Syria" = "Syrian Arab Republic",
                  "Tanzania" = "United Republic of Tanzania", 
                  "Union of Soviet Socialist Repu" = "USSR", "United States" =
                  "United States of America", "Venezuela" = 
                  "Venezuela (Bolivarian Republic of)", "Vietnam" = "Viet Nam")
Rkbyc$Country <- str_replace_all(Rkbyc$Country, country_repl)
Rkbyc$Country <- str_replace(Rkbyc$Country, fixed("Congo (Brazzaville)"), "Congo") %>%
                 str_replace(fixed("Congo (Kinshasa)"), "Democratic Republic of the Congo") %>%
                 str_replace(fixed("Yemen (Aden)"), "Yemen Dem") %>%
                 str_replace(fixed("Yemen (Sanaa)"), "Yemen Ar Rp")
cnames <- Rkbyc$Country
Rkbyc <- as.matrix(Rkbyc[, -1])
rownames(Rkbyc) <- cnames

# Add Taiwan and HK to China, South Sudan to Sudan, and two Yemens
Rkbyc["China", ] <- Rkbyc["China", ] + Rkbyc["Hong Kong", ] + Rkbyc["Taiwan", ]
Rkbyc["Sudan", ] <- Rkbyc["Sudan", ] + Rkbyc["South Sudan", ]
Rkbyc["Yemen", ] <- Rkbyc["Yemen", ] + Rkbyc["Yemen Ar Rp", ] + Rkbyc["Yemen Dem", ]
Rkbyc <- Rkbyc[!(rownames(Rkbyc) %in% c("Hong Kong", "Taiwan", "South Sudan", 
                                        "Yemen Ar Rp", "Yemen Dem")), ]

# Move pre-1993 Ethiopia to Ethiopia PDR
Rkbyc <- rbind(Rkbyc, "Ethiopia PDR" = rep(0, ncol(Rkbyc)))
Rkbyc["Ethiopia PDR", 1:7] <- Rkbyc["Ethiopia", 1:7]
Rkbyc["Ethiopia", 1:7] <- 0

# Match break points for USSR, Yugoslavia, Czechoslovakia
match_list <- list(
    list(former = "USSR", year = 1992, 
         new_list = c("Armenia", "Azerbaijan", "Belarus", "Estonia", "Georgia", 
                      "Kazakhstan", "Kyrgyzstan", "Latvia", "Lithuania", 
                      "Republic of Moldova", "Russian Federation", "Tajikistan", 
                      "Turkmenistan", "Ukraine", "Uzbekistan")),
    list(former = "Yugoslav SFR", year = 1992,
         new_list = c("Bosnia and Herzegovina", "Croatia", "Serbia and Montenegro", 
                      "Slovenia", "The former Yugoslav Republic of Macedonia")),
    list(former = "Czechoslovakia", year = 1993,
         new_list = c("Czech Republic", "Slovakia"))
)

for (l in match_list) {
    inew <- which(rownames(Rkbyc) %in% l$new_list)
    yrs <- which(as.numeric(colnames(Rkbyc)) < l$year)
    Rkbyc[l$former, yrs] <- Rkbyc[l$former, yrs] + colSums(Rkbyc[inew, yrs])
    Rkbyc[inew, yrs] <- 0
}


## (3) Apportion European Union reserves proportionally to production

eu15_list <- c("Austria", "Belgium", "Denmark", "Finland", "France", "Germany",
               "Greece", "Ireland", "Italy", "Luxembourg", "Netherlands", 
               "Portugal", "Spain", "Sweden", "United Kingdom", 
               "Belgium-Luxembourg") 
eu_list <- c(eu15_list, "Bulgaria", "Croatia", "Cyprus", "Czech Republic",
             "Estonia", "Hungary", "Latvia", "Lithuania", "Malta", "Poland", 
             "Romania", "Slovakia", "Slovenia")
eu15_yrs <- which(as.numeric(colnames(Rkbyc)) < 1998)
eu_yrs <- which(as.numeric(colnames(Rkbyc)) >= 1998)

# Load cereals production data (in object Pkbyc)
load("data/cereals_prod_trade.RData")
rownames(Pkbyc) <- ciso$name

prop_eu15 <- scale(Pkbyc[eu15_list, eu15_yrs], center = FALSE,
                   scale = colSums(Pkbyc[eu15_list, eu15_yrs]))
prop_eu <- scale(Pkbyc[eu_list, eu_yrs], center = FALSE,
                 scale = colSums(Pkbyc[eu_list, eu_yrs]))

# Add individual countries to matrix
Rkbyc <- rbind(Rkbyc, matrix(0, nrow = length(eu15_list), ncol = ncol(Rkbyc), 
                             dimnames = list(eu15_list)))
Rkbyc[eu15_list, eu15_yrs] <- sweep(prop_eu15, 2, Rkbyc["EU-15", eu15_yrs], "*")
Rkbyc[eu_list, eu_yrs] <- sweep(prop_eu, 2, Rkbyc["European Union", eu_yrs], "*")
Rkbyc <- Rkbyc[-which(rownames(Rkbyc) %in% c("EU-15", "European Union")), ]


## (4) Create new matrix with all countries in Pkbyc fill reserves by matching names
Rkbyc_all <- matrix(0, nrow = nrow(Pkbyc), ncol = ncol(Pkbyc))
Rkbyc_all[match(rownames(Rkbyc), ciso$name),] <- Rkbyc

saveRDS(Rkbyc_all, "data/cereals_stocks.RData")

