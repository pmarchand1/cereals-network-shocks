# Functions for Trade and Food Security cascade model

library(igraph) # for network metrics

# Function to load input data for the simulation model
#  based on the structure of .RData files in data/fao_networks
# Also requires "ciso3.txt" in working directory to map ISO3 country codes
# mov_avg is the number of years to average on each side of the focal year
get_trade_data <- function(year, prod_trade_file, stocks_file = NA, mov_avg = 0) {
    # Get country names
    if (!("ciso3.txt" %in% dir())) stop("ciso3.txt not found in working directory")
    iso3 <- read.table("ciso3.txt")
    cnames <- iso3[, 2]
    # Check that year range is in data
    data_years <- 1986:2011
    yrange <- year + -mov_avg:mov_avg
    miss_years <- !(yrange %in% data_years)
    if (any(miss_years)) stop(paste("no data for", yrange[miss_years]))
    # Get production and trade data
    load(prod_trade_file)
    dimnames(Pkbyc) <- list(cnames, data_years)
    dimnames(Tkbyc) <- list(cnames, cnames, data_years)
    yi <- which(data_years %in% yrange)
    if (mov_avg == 0) {
        P0 <- Pkbyc[, yi]
        E0 <- Tkbyc[, , yi]             
    } else {
        P0 <- rowMeans(Pkbyc[, yi])
        E0 <- apply(Tkbyc[, , yi], 1:2, mean)   
    }
    # Get reserves if applicable
    if (is.na(stocks_file)) {
        R0 <- rep(0, length(cnames))
        names(R0) <- cnames
    } else {
        Rkbyc <- readRDS(stocks_file)
        dimnames(Rkbyc) <- list(cnames, data_years)
        if (mov_avg == 0)  {
            R0 <- Rkbyc[, yi]
        } else {
            R0 <- rowMeans(Rkbyc[, yi])   
        }        
    }
    # Remove countries with no trade data
    no_dat <- which(rowSums(E0) == 0 & colSums(E0) == 0)
    P0 <- P0[-no_dat]
    R0 <- R0[-no_dat]
    E0 <- E0[-no_dat, -no_dat]
    
    list(P0 = P0, R0 = R0, E0 = E0)
}


# Function to simulate cascade starting with a change (decline) in production dP
#  P0, R0 and E0 are the initial production vector, reserve vector and trade matrix
#  cfrac is the fraction of shock absorbed by C before changing trade
#  if asym = TRUE, shocked countries can't increase their exports
#  kmax is the max number of iterations
#  any shock below amin (as fraction of net supply) is negligible
#  
#  anim_out produces a different output to be sent to sim_animation
sim_cascade <- function(dP, P0, R0, E0, cfrac = 0, asym = TRUE,
                        kmax = 50, amin = 1E-5, anim_out = FALSE) {
    nc <- length(P0)
    cnames <- names(P0)
    
    # Check input validity
    if (!is.vector(P0)) stop("P0 must be a vector.")
    if (!is.vector(R0)) stop("R0 must be a vector.")
    if (!is.matrix(E0)) stop("E0 must be a matrix.")
    if (length(R0) != length(P0)) stop("dimensions of P0 and R0 don't match.")
    if (any(dim(E0) != length(P0))) stop("dimensions of P0 and E0 don't match.")

    if (cfrac < 0 || cfrac > 1) stop("cfrac must be between 0 and 1")
    
    # Initialize arrays
    R <- matrix(0, nc, kmax, dimnames = list(cnames))
    R[, 1] <- R0
    dR <- matrix(0, nc, kmax, dimnames = list(cnames))
    dC <- matrix(0, nc, kmax, dimnames = list(cnames))
    E <- array(0, dim = c(nc, nc, kmax), dimnames = list(cnames, cnames))
    E[, , 1] <- E0
    dE <- array(0, dim = c(nc, nc, kmax), dimnames = list(cnames, cnames))
    S <- matrix(0, nc, kmax, dimnames = list(cnames))
    # Net supply = production + imports - exports
    S[, 1] <- P0 + colSums(E0) - rowSums(E0)
    shocked <- rep(FALSE, nc)

    # Initial shock (dS = drop in supply)
    dS <- dP
    
    # Main loop
    for (k in 1:(kmax-1)) {
        # Shock first absorbed by reserves, some consumption (cfrac of shock)
        # then trade (up to Tvol = E + I), then consumption. 
        # Countries don't pass a shock smaller than amin * S, just absorb
        # if asym = TRUE, shocked countries with R = 0 can't export more,
        #  i.e. E[i,j] doesn't count in j's imports if shocked[i]
        dR[, k] <- pmax(dS, -R[, k])
        res_shock <- dS - dR[, k]
        if (asym) shocked[res_shock < 0] <- TRUE
        dC[, k] <- cfrac * res_shock
        res_shock <- (1 - cfrac) * res_shock
        Tvol <- rowSums(E[, , k]) + colSums(E[!shocked, , k])
        dTvol <- ifelse(abs(res_shock) < amin * S[, k], 0, pmax(res_shock, -Tvol))
        dC[, k] <- dC[, k] + res_shock - dTvol

        # Countries allocate change in trade volume (+ imports, - exports) 
        #  proportionally to the initial volume of trade on each link
        #  (see dEalloc function)
        prop_dT <- ifelse(Tvol == 0, 0, dTvol / Tvol)
        dE[, , k] <- dEalloc(prop_dT, shocked) * E[, , k]

        # Update values for next iteration
        # Decrease in net supply corresponds to shock absorbed internally (R, C)
        R[, k+1] <- R[, k] + dR[, k]
        E[, , k+1] <- E[, , k] + dE[, , k]
        S[, k+1] <- S[, k] + dR[, k] + dC[, k]  
        
        # Decrease in I-E due to other countries corresponds to next dS
        dS <- colSums(dE[, , k]) - rowSums(dE[, , k]) + dTvol

        # If no more decrease in production, returns a list with all simulation
        #  parameters and the sum (over iterations) of dC, dE and dR
        # 
        # If anim_out is TRUE, return dR, dC and dE for each iteration
        if (isTRUE(all.equal(min(dS), 0))) {
            if (anim_out) {
                return(list(dP = dP, P0 = P0, R0 = R0, E0 = E0, dR = dR[, 1:k], 
                            dC = dC[, 1:k], dE = dE[, , 1:k]))
            } else {
                return(list(dP = dP, P0 = P0, R0 = R0, E0 = E0, dR = rowSums(dR),
                            dC = rowSums(dC), dE = apply(dE, 1:2, sum)))
            }
        }
    }
    warning("No equilibrium reached after maximum number of iterations")
    if (anim_out) {
        list(dP = dP, P0 = P0, R0 = R0, E0 = E0, dR = dR, dC = dC, dE = dE)
    } else {
        list(dP = dP, P0 = P0, R0 = R0, E0 = E0, dR = rowSums(dR), 
             dC = rowSums(dC), dE = apply(dE, 1:2, sum))
    }
}


# Function to calculate the relative change in each trade link
#  from the proportion of trade volume (by country) affected by shock (prop_dT)
#  ana a logical vector indicating shocked countries (that an't export any more)
dEalloc <- function(prop_dT, shocked) {
    nc <- length(shocked)
    # Results in dEprop[i,j] = prop_dT[i] - (!shocked[i]) * prop_dT[j] 
    dEprop <- rep(prop_dT, nc) - 
              rep(!shocked, nc) * rep(prop_dT, each = nc)
    dim(dEprop) <- c(nc, nc)
    diag(dEprop) <- 0
    dEprop
}


# This function takes the output of sim_cascade (or sim_1c)
#  and tests whether it respects equation: S = P + I - E = R + C
sim_diagnostics <- function(sim_res, tol = 1E-5) {
    S0 <- sim_res$P0 + colSums(sim_res$E0) - rowSums(sim_res$E0) 
    dS0 <- sum(sim_res$dP)
    cnames <- names(sim_res$P0)
    
    # Test 1: Total dR + dC must match initial shock
    # (within relative difference of tol)
    discrep <- (sum(sim_res$dR) + sum(sim_res$dC))/dS0 - 1
    if (abs(discrep) > tol)
        return(paste("Total dC + dR does not match initial shock.",
                     "Relative difference:", discrep))
    
    # Test 2: After initial shock, dI - dE = dR + dC by country
    # (within relative difference of tol * final net supply)
    dS_byc <- sim_res$dP + colSums(sim_res$dE) - rowSums(sim_res$dE)
    discrep <- which(abs((dS_byc - sim_res$dR - sim_res$dC)/(dS_byc + S0)) > tol)
    if (length(discrep) > 0)
        return(paste("Net supply change does not match dC + dR for countries:",
                     paste(cnames[discrep], collapse = ", ")))
    
    return("All tests passed.")
}


# Runs a sim_cascade starting with country c_init losing a fraction a_init of P
sim_1c <- function(c_init, a_init, P0, R0, E0, cfrac = 0, asym = TRUE,
                   kmax = 50, amin = 1E-5, anim_out = FALSE) {
    nc <- length(P0)
    cnames <- names(P0)
    if (!(c_init %in% cnames) && !(c_init %in% 1:nc)) {
        stop("invalid value of c_init")
    }
    nc <- length(P0)
    if (a_init < 0 || a_init > 1) stop("a_init must be between 0 and 1")
    
    dP <- rep(0, nc)
    names(dP) <- cnames
    dP[c_init] <- -a_init * P0[c_init]
    sim_cascade(dP, P0, R0, E0, cfrac, asym, kmax, amin, anim_out)
}


# Runs sim_1c for each country represented in input data where P0 > 0
# The list it returns is similar to that of sim_cascade, except that
#  dR, dC and dE have an additional dimension corresponding to the country
#  where the cascade was initiated.
sim_allc <- function(a_init, P0, R0, E0, cfrac = 0, asym = TRUE, 
                     kmax = 50, amin = 1E-5) {
    cnames <- names(P0)
    csim <- cnames[P0 > 0]
    nc <- length(cnames)
    nsim <- length(csim)
    res_multi <- lapply(csim, sim_1c, a_init, P0, R0, E0, 
                        cfrac, asym, kmax, amin, anim_out = FALSE)
    list(a_init = a_init, P0 = P0, R0 = R0, E0 = E0, 
         dR = matrix(vapply(res_multi, `[[`, "dR", FUN.VALUE = rep(0, nc)), 
                     nrow = nc, dimnames = list(cnames, csim)),
         dC = matrix(vapply(res_multi, `[[`, "dC", FUN.VALUE = rep(0, nc)), 
                     nrow = nc, dimnames = list(cnames, csim)),                            
         dE = array(vapply(res_multi, `[[`, "dE", FUN.VALUE = rep(0, nc*nc)), 
                    dim = c(nc, nc, nsim), dimnames = list(cnames, cnames, csim)))
}


# Calculates the depth of a single cascade simulation based on c_init and dE
#  Depth is the maximum graph distance from c_init to any country hit by shock
cascade_depth <- function(c_init, dE) {
    if (all(dE == 0)) return(0) # Depth = 0 if no change in E
    trade_graph <- graph_from_adjacency_matrix(dE != 0, mode = "directed")
    # Remove countries that are not part of graph (i.e. not affected)
    trade_graph <- delete_vertices(trade_graph, which(degree(trade_graph) == 0))
    spaths <- shortest_paths(trade_graph, from = c_init, mode = "all")
    max(sapply(spaths$vpath, length))
}


# This function calculates different summary stats for one sim_res_multi
#  (i.e. outcome of sim_allc) and arranges them in a list
get_stats_allc <- function(sim_res_multi) {
    # Initial net supply by country
    S0 <- with(sim_res_multi, P0 + colSums(E0) - rowSums(E0))
    # Magnitude of initial shock by sim
    csims <- colnames(sim_res_multi$dR)
    nsim <- length(csims)
    dS0 <- with(sim_res_multi, a_init * P0[csims]) 
    
    # Values of dR, dC, dE relative to initial shock
    dRrel <- sweep(sim_res_multi$dR, 2, dS0, "/")
    dCrel <- sweep(sim_res_multi$dC, 2, dS0, "/")
    dErel <- sweep(sim_res_multi$dE, 3, dS0, "/")
    dSrel <- dRrel + dCrel  
    # Values of dC and dS relative to initial net supply
    dC_S0 <- sweep(sim_res_multi$dC, 1, S0, "/")
    dS_S0 <- sweep(sim_res_multi$dR, 1, S0, "/") + dC_S0
    
    # Statistics by simulation
    depth_by_sim <- vapply(csims, 
                       function(x) cascade_depth(x, sim_res_multi$dE[, , x]), 0)
    countries_hit_by_sim <- colSums(dSrel != 0) 
    links_hit_by_sim <- apply(dErel != 0, 3, sum)
    total_dC_by_sim <- colSums(dCrel)
        
    # Statistics by country and link
    hits_by_cty <- rowSums(dSrel != 0) 
    hits_by_link <- apply(dErel !=0, 1:2, sum)
    avg_links_hit_by_cty <- (rowSums(hits_by_link) + colSums(hits_by_link)) / nsim
    
    # Rank-effect plots of shocks over countries (avg. by sim)
    #  and over sims (avg. by country)
    avg_dSrel_by_cty_rank <- rowMeans(-apply(dSrel, 2, sort))
    avg_dS_S0_by_sim_rank <- rowMeans(-apply(dS_S0, 1, sort)) 

    
    list(dRrel = dRrel, dCrel = dCrel, dC_S0 = dC_S0, dS_S0 = dS_S0,
         depth_by_sim = depth_by_sim, countries_hit_by_sim = countries_hit_by_sim, 
         links_hit_by_sim = links_hit_by_sim, total_dC_by_sim = total_dC_by_sim,
         hits_by_cty = hits_by_cty, hits_by_link = hits_by_link, 
         avg_links_hit_by_cty = avg_links_hit_by_cty,
         avg_dSrel_by_cty_rank = avg_dSrel_by_cty_rank, 
         avg_dS_S0_by_sim_rank = avg_dS_S0_by_sim_rank)
}


# Get summary stats from trade data (i.e. output of get_trade_data function)
get_trade_stats_sum <- function(trade_data) {
    nc <- length(trade_data$P0)
    Ptot <- sum(trade_data$P0)
    Rtot <- sum(trade_data$R0)
    nlinks <- sum(trade_data$E0 > 0)
    totflow <- sum(trade_data$E0)
    RP_ratio <- Rtot / Ptot
    flowP_ratio <- totflow / Ptot
    data.frame(nc = nc, Ptot = Ptot, Rtot = Rtot, nlinks = nlinks, totflow = totflow,
               RP_ratio = RP_ratio, flowP_ratio = flowP_ratio)
}


# Get derived stats by country from trade data (output of get_trade_data function)
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
