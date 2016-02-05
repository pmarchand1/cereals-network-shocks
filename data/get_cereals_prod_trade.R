### This script uses FAOSTAT data from 1986 to 2011 data to calculate 
###    the cereals production (countries x years matrix) 
###    and trade (countries x countries x years matrix)
# Author: Philippe Marchand (pmarchand@sesync.org)

library(tidyr)
library(dplyr)
library(FAOSTAT)
library(reshape2)

# Load crop list
group_id <- 1  # Cereals
crop_list <- read.delim("data/crop_list.tsv", header = TRUE, as.is = TRUE) %>%
    filter(Groupnum == group_id)
prod_crops <- filter(crop_list, ProductionL == 1, cropname != "Bulgur") %>%
    select(cropid, cropname, kcal.ton)
trade_crops <- filter(crop_list, TradeL == 1) %>%
    select(cropid, kcal.ton)
country_list <- read.csv("data/country_list.csv")
yr_range <- 1986:2011


## PART 1 : Production data (acquired via FAOSTAT package)

# Get production data from FAOSTAT
#  QC: production-crops domain, 5510: production in tonnes
prod_dat <- getFAOtoSYB(name = prod_crops$cropname, itemCode = prod_crops$cropid, 
                        domainCode = rep("QC", nrow(prod_crops)),
                        elementCode = rep(5510, nrow(prod_crops)),
                        yearRange = yr_range, countrySet = country_list$FAOST_CODE,
                        useCHMT = FALSE, outputFormat = "long")
prod_dat <- prod_dat$entity
prod_dat <- prod_dat[, c("FAOST_CODE", "Year", "Value", "itemCode", "name")]

# Add kcal conversion factor
prod_agg <- inner_join(prod_dat, prod_crops, by = c("itemCode" = "cropid")) %>%
    mutate(value_kcal = Value * as.numeric(kcal.ton)) %>%
    group_by(FAOST_CODE, Year) %>%
    summarise(pkcal = sum(value_kcal, na.rm = TRUE))

# Convert to country x year matrix
prod_mat <- spread(prod_agg, key = Year, value = pkcal, fill = 0)
prod_mat <- prod_mat[match(country_list$FAOST_CODE, prod_mat$FAOST_CODE), ]
prod_mat <- as.matrix(prod_mat[, -1])
prod_mat[is.na(prod_mat)] <- 0
# Only keep data for country with pop > 500k
prod_mat[!country_list$gt.500k, ] <- 0
colnames(prod_mat) <- NULL

# Clear unneded files
rm(prod_dat, prod_agg)


### PART 2: Trade (Detailed trade matrix pre-downloaded trade from FAOSTAT)
trade_dat <- read.csv("data/Trade_DetailedTradeMatrix_E_All_Data_(Norm).csv",
                      stringsAsFactors = FALSE)

trade_dat <- select(trade_dat, reporter = Reporter.Country.Code, 
                    partner = Partner.Country.Code, cropid = Item.Code, 
                    element = Element, year = Year, value = Value) %>%
    filter(year %in% yr_range)

# Put parts of China under same ID
china_ids <- c(41, 96, 128, 214)
trade_dat$reporter[trade_dat$reporter %in% china_ids] <- 351
trade_dat$partner[trade_dat$partner %in% china_ids] <- 351

# Put item 30 (rice total - milled eq.) under 31 (rice milled)
trade_dat$cropid[trade_dat$cropid == 30] <- 31

# Keep countries in country_list and remove trade between a country and itself
trade_dat <- filter(trade_dat, reporter %in% country_list$FAOST_CODE,
                    partner %in% country_list$FAOST_CODE) %>%
    filter(reporter != partner)

# Extract exports
exp_dat <- filter(trade_dat, element == "Export Quantity") %>% select(-element)

# Convert to kcal and aggregate across crops
exp_agg <- inner_join(exp_dat, trade_crops) %>%
    mutate(value_kcal = value * as.numeric(kcal.ton)) %>%
    group_by(reporter, partner, year) %>%
    summarise(ekcal = sum(value_kcal, na.rm = TRUE))

# Clear unneeded large dataframes
rm(trade_dat, exp_dat)


# Convert export datafarme to a country x country x year matrix
clist <- as.character(country_list$FAOST_CODE)
exp_mat <- acast(exp_agg, reporter ~ partner ~ year, fill = 0)
exp_mat <- exp_mat[match(clist, dimnames(exp_mat)[[1]]), , ]
exp_mat <- exp_mat[, match(clist, dimnames(exp_mat)[[2]]), ]
exp_mat[is.na(exp_mat)] <- 0

# Remove countries where max population < 500k over time period
exp_mat[!country_list$gt.500k, , ] <- 0
exp_mat[, !country_list$gt.500k, ] <- 0
# Remove PRI, REU, ESH (incomplete time series)
exp_mat[, c(181, 185, 248), ] <- 0
exp_mat[c(181, 185, 248), , ] <- 0
dimnames(exp_mat) <- NULL


# Save production and trade to single .RData file
Pkbyc <- prod_mat
Tkbyc <- exp_mat
save(Pkbyc, Tkbyc, "cereals_prod_trade.RData")



