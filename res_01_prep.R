library(dplyr)
library(tidyr)
library(ggplot2)

data_dir <- "/nfs/tradefoodsec-data/data/fao_networks"
res_dir <- "/nfs/tradefoodsec-data/analysis/sim_results/cereals"

source("tfs_cascade_funcs.R")
source("tfs_map_funcs.R")
source("res_00_funcs.R")

# Load base trade data
years <- c(1996, 2003, 2009)
prod_trade_file <- file.path(data_dir, "1 - CEREALS AND CEREAL PRODUCTS.RData")
stocks_file <- file.path(data_dir,"cereals_stocks.RData")
trade_dat <- lapply(years, get_trade_data, prod_trade_file, stocks_file,
                    mov_avg = 2)
names(trade_dat) <- years

# Calculate summaries of trade data
trade_st <- bind_rows(lapply(trade_dat, get_trade_stats_sum))
# % change between years
# data.frame(lapply(trade_st, function(x) x / lag(x)))

# Load base map and convert to df
data(wrld_simpl)
world_df <- fortify(wrld_simpl)

# To run and save simulation results...
# res_02_05 <- lapply(trade_dat, function(x) sim_allc(0.2, x$P0, 0.5*x$R0, x$E0, cfrac = 0.01)) 
# saveRDS(res_02_05, file.path(res_dir, "res_02_05.RData"))


# Run sims with 2009 Rtot/Ptot ratio
#Rmult <- trade_st$RP_ratio[3] / trade_st$RP_ratio
#for (i in 1:3) trade_dat[[i]]$R0 <- trade_dat[[i]]$R0 * Rmult[i]
#res_Rscaled <- lapply(trade_dat, function(x) sim_allc(0.2, x$P0, 0.5*x$R0, x$E0, cfrac = 0.01))
#saveRDS(res_Rscaled, file.path(res_dir, "res_Rscaled.RData"))

# Run sims with 2009 R0/S0 for all C
#S0 <- lapply(trade_dat, function(x) x$P0 + colSums(x$E0) - rowSums(x$E0))
#R0_S0 <- lapply(1:3, function(i) trade_dat[[i]]$R0 / S0[[i]])
#for (i in 1:3) trade_dat[[i]]$R0 <- R0_S0[[3]][names(S0[[i]])] * S0[[i]] 
#res_RSscaled <- lapply(trade_dat, function(x) sim_allc(0.2, x$P0, 0.5*x$R0, x$E0, cfrac = 0.01))
#saveRDS(res_RSscaled, file.path(res_dir, "res_RSscaled.RData"))

###

# Get ALL the results
res_list <- c("02_05", "01_025", "03_1", "015_05", 
              "02_05_005", "02_05_01", "Rscaled", "RSscaled")

trade_stats_byc <- bind_rows(lapply(trade_dat, get_trade_stats_by_cty), 
                             .id = "year")

get_stats_from_sim_name <- function(sim_name) {
    res <- readRDS(file.path(res_dir, paste0("res_", sim_name, ".RData")))
    st <- lapply(res, get_stats_allc)
    stats_df <- cbind(params = sim_name, 
                      bind_rows(lapply(st, get_sim_stats_by_cty)))
}

stats_by_cty <- bind_rows(lapply(res_list, function(x) {
    cbind(trade_stats_byc, get_stats_from_sim_name(x))
}))
rm(trade_stats_byc)

# rank_dS_S0 for all sims
get_rank_dS_S0_from_sim_name <- function(sim_name) {
    res <- readRDS(file.path(res_dir, paste0("res_", sim_name, ".RData")))
    st <- lapply(res, get_stats_allc)
    res_df <- cbind(params = sim_name, 
                    bind_rows(lapply(st, get_rank_dS_S0), .id = "year"))
}
rank_dS_S0 <- bind_rows(lapply(res_list, get_rank_dS_S0_from_sim_name))


###

# Load results
#res_02_05 <- readRDS(file.path(res_dir, "res_02_05.RData"))
#st <- lapply(res_02_05, get_stats_allc)
#rm(res_02_05)

# Get all stats by country in a data frame
#stats_by_cty <- cbind
#    bind_rows(lapply(trade_dat, get_trade_stats_by_cty), .id = "year"),
#    bind_rows(lapply(st, get_sim_stats_by_cty))
#)

# Add gdp_percap from World Bank and more derived stats
library(WDI)
gdp_pc <- WDI("all", "NY.GDP.PCAP.PP.KD", start = 1996, end = 2009, extra = TRUE)
gdp_pc <- select(gdp_pc, cty = iso3c, year, gdp_pc = NY.GDP.PCAP.PP.KD,
                 income = income) %>% mutate(year = as.factor(year))
stats_by_cty <- left_join(stats_by_cty, gdp_pc) %>%
    mutate(dC_by_hit = tot_dC_S0 / hitsC, 
           dS_by_hit = tot_dS_S0 / hits,
           tot_dR_S0 = tot_dS_S0 - tot_dC_S0) %>%
    arrange(year, cty)
rm(gdp_pc)


# Overall clusters 1995-2011
cls <- cluster_ma(2003, prod_trade_file, stocks_file, mov_avg = 8)
tmean <- get_trade_data(2003, prod_trade_file, stocks_file, mov_avg = 8)
tdat <- data.frame(country = names(tmean$P0), P = tmean$P0, R = tmean$R0,
                   I = colSums(tmean$E0), E = rowSums(tmean$E0))
tdat[,-1] <- sweep(tdat[,-1], 1, tdat$P + tdat$I - tdat$E, "/")
grps_df <- data.frame(cty = tdat$country, group = cutree(cls, 5))
grps_df <- grps_df[!duplicated(grps_df$cty), ]

# Merge group info with other data 
stats_by_cty <- inner_join(stats_by_cty, grps_df)
stats_by_cty$group <- as.factor(stats_by_cty$group)

rank_dS_S0 <- inner_join(rank_dS_S0, grps_df)
rank_dS_S0$group <- as.factor(rank_dS_S0$group)

rm(cls, tmean, tdat, grps_df)



# Save stats_by_cty and rank_dS_S0
saveRDS(stats_by_cty, file.path(res_dir, "stats_by_cty.RData"))
saveRDS(rank_dS_S0, file.path(res_dir, "rank_dS_S0.RData"))


