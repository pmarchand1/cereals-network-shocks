### Some utility functions for the Dec.2015 results


# Get summary stats from trade data (i.e. output of get_trade_data function)
get_trade_stats_sum <- function(trade_data) {
    Ptot <- sum(trade_data$P0)
    Rtot <- sum(trade_data$R0)
    nlinks <- sum(trade_data$E0 > 0)
    totflow <- sum(trade_data$E0)
    RP_ratio <- Rtot / Ptot
    flowP_ratio <- totflow / Ptot
    data.frame(Ptot = Ptot, Rtot = Rtot, nlinks = nlinks, totflow = totflow,
               RP_ratio = RP_ratio, flowP_ratio = flowP_ratio)
}


# Get derived stats by country from trade data (output of get_trade_data function)
# "year" added in as a parameter to facilitate calling this for multiple years
#  and rbind-ing the output
get_trade_stats_by_cty <- function(trade_data) {
    tgraph <- graph_from_adjacency_matrix(trade_data$E0 != 0, mode = "directed")
    S0 <- trade_data$P0 + colSums(trade_data$E0) - rowSums(trade_data$E0)
    data.frame(cty = names(trade_data$P0), S0 = S0,
               R0_S0 = trade_data$R0 / S0, 
               I0_S0 = colSums(trade_data$E0) / S0,
               E0_S0 = rowSums(trade_data$E0) / S0,
               indeg = degree(tgraph, mode = "in"),
               outdeg = degree(tgraph, mode = "out"))
}


# Utility function to compute Pielou evenness of a vector
evenness <- function(x) {
    x <- x / sum(x)
    - sum(x[x > 0] * log(x[x > 0])) / log(length(x))
}

# Utility function to get self_dC vector from dC matrix
get_dC_self <- function(dC) {
    vapply(rownames(dC), function(cty) {
        if (cty %in% colnames(dC)) dC[cty, cty]
        else 0
    }, 0)
}

# Get a few summary statistics by country from the output of get_stat_allc function
get_sim_stats_by_cty <-  function(stat_list) {
    data_frame(hits = stat_list$hits_by_cty,
               links_hit = stat_list$avg_links_hit_by_cty,
               hitsC = rowSums(stat_list$dCrel != 0),
               tot_dC_S0 = rowSums(stat_list$dC_S0),
               self_dC_S0 = get_dC_self(stat_list$dC_S0),
               tot_dS_S0 = rowSums(stat_list$dS_S0),
               even = apply(stat_list$dS_S0, 1, evenness))
}


# Maps country-level stats; assumes that world_df is in environment,
#  that df has "cty" column with iso3 codes 
#  as well as "year" and "params" columns (unless corresponding arg. is NA)
stat_cty_map <- function(df, stat, year = NA, params = NA,
                         low = "red", high = "white", div = FALSE) {
    select_expr <- rep(TRUE, nrow(df))
    if (!is.na(year)) select_expr <- select_expr & df$year == year
    if (!is.na(params)) select_expr <- select_expr & df$params == params
    stat_select <- df[select_expr, stat]
    names(stat_select) <- df[select_expr, "cty"]
    world_df <- add_country_data(world_df, stat_select)
    if (div) {
        ggplot() + geom_polygon(aes(long, lat, group = group, fill = value),
                                data = world_df, colour = "grey") +
            scale_fill_gradient2(low = low, high = high, na.value = "grey") + 
            coord_fixed() + theme_map()
    } else {
        ggplot() + geom_polygon(aes(long, lat, group = group, fill = value),
                                data = world_df, colour = "grey") +
            scale_fill_gradient(low = low, high = high, na.value = "grey") + 
            coord_fixed() + theme_map()        
    }     
}


# Computes difference in stats by countr b/w two years
#  assumes df has "cty" and "year" columns, as well as params (if arg. not NA)
#  assumes all countries in y1 also in y2 (but not vice versa)
stats_dif <- function(df, y1, y2, params = NA) {
    if (!is.na(params)) df <- df[df$params == params, ]
    df1 <- df[df$year == y1, ]
    df2 <- df[df$year == y2, ]
    df2 <- df2[df2$cty %in% df1$cty, ]
    if (!identical(df1$cty, df2$cty)) stop("countries don't match")
    num_cols <- sapply(df, is.numeric)
    cbind(cty = df1$cty, df2[, num_cols] - df1[, num_cols])
}


# Returns a hierarchical clustering of countries in year (+/- mov_avg)
#  using production, stocks, exports and imports data (scaled by net supply)
#  in addition to specified files, loads ciso3.txt for ISO3 country codes
cluster_ma <- function(year, prod_trade_file, stocks_file, mov_avg) {
    # Get country names
    iso3 <- read.table("ciso3.txt")
    cnames <- iso3[, 2]
    # Year range in our data networks
    years <- 1986:2011
    # Get production and trade data
    load(prod_trade_file)
    dimnames(Pkbyc) <- list(cnames, years)
    dimnames(Tkbyc) <- list(cnames, cnames, years)
    Rkbyc <- readRDS(stocks_file)
    dimnames(Rkbyc) <- list(cnames, years)
    # Calculate moving average for year +/- 2
    yrange <- which(years == year) + -mov_avg:mov_avg
    P <- rowMeans(Pkbyc[, yrange])
    R <- rowMeans(Rkbyc[, yrange])
    Tmat <- apply(Tkbyc[, , yrange], 1:2, mean)
    I <- colSums(Tmat)
    E <- rowSums(Tmat)
    all_dat <- cbind(P, R, I, E)
    # Remove countries with no trade data
    no_dat <- which(I == 0 & E == 0)
    all_dat <- all_dat[-no_dat, ]
    S <- all_dat[, "P"] + all_dat[, "I"] - all_dat[, "E"]
    all_dat <- sweep(all_dat, 1, S, "/")
    # Cluster
    dmat <- dist(all_dat)
    hcl <- hclust(dmat, method = "ward.D")
}


# Produces a data frame of dS_S0 vs. sim.rank (in order of most impact) by country
#   frac: dS_S0 as fraction to total impact to country
#   cumul: also include cumulative dS_S0 (for all ranks <= rank)
get_rank_dS_S0 <- function(stat_list, frac = TRUE, cumul = TRUE) {
    rank_dS_S0 <- -apply(stat_list$dS_S0, 1, sort)
    if (frac) {
        rank_dS_S0 <- scale(rank_dS_S0, center = FALSE, 
                            scale = colSums(rank_dS_S0))
    }
    rank_dS_S0 <- data.frame(rank = 1:nrow(rank_dS_S0), rank_dS_S0) %>%
        gather(key = "cty", value = "dS_S0", -rank)
    if (cumul) {
        rank_dS_S0 <- group_by(rank_dS_S0, cty) %>%
            mutate(cumul = cumsum(dS_S0))
    }
    rank_dS_S0
}
