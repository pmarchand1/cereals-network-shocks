library(ggplot2)
library(grid) # for plot arrows
library(geosphere) # for great circles
library(maptools) # for "wrld_simpl" map
library(tidyr) # for gather function
library(animation)

###
# Functions to help map output of tfs_cascade
###


# Global parameter: points by great circle line
pts_by_gc <- 100


# Theme for mapping
theme_map <- function(base_size = 12, base_family = "Helvetica") {
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
        theme(rect = element_blank(), line = element_blank(),
              axis.text = element_blank(), axis.title = element_blank(),
              axis.ticks.margin = unit(0, "lines"))
}


# This function adds country-level metrics to a map data frame for plotting
#   the metrics vector vect has names corresponding to the country codes
#   under "id" in map_df
add_country_data <- function(map_df, vect) {
    # clear previous values
    if ("value" %in% colnames(map_df)) map_df$value <- NULL 
    metrics <- data.frame(names(vect), vect)
    colnames(metrics) <- c("id", "value")
    map_df <- merge(map_df, metrics, by = "id", all.x = TRUE, sort = FALSE)
    # reorder df to make sure polygons are drawn properly
    map_df[order(map_df$group, map_df$order),] 
}


# This function converts a trade matrix trade_mat
#   (i.e. country_from in rows, country_to in columns, 
#    row and column names must match country ids)
# into a table with 3 columns: "from", "to", "value"
#   keeps a maximum of max_rows in order of decreasing absolute value
trade_tab <- function(trade_mat, max_rows) {
    trade_df <- as.data.frame(trade_mat)
    trade_df <- cbind(from = colnames(trade_df), trade_df)
    trade_df <- gather(trade_df, "to", "value", -from)
    # Return max_lines rows with largest absolute values
    trade_df[order(-abs(trade_df$value))[1:max_rows],]
}


# This function creates great circle lines for any trade matrix 
#   (i.e. country_from in rows, country_to in columns, 
#    row and column names must match country ids)
#   and a data frame of country centroid coordinates centr;
#   keeps a maximum of max_lines in order of descending absolute value;
#   returns a data frame that can be added to a ggplot
get_trade_lines <- function(trade_mat, centers, max_lines) {
    # Convert matrix to "from", "to", "value" table, then add country coordinates
    trade_df <- trade_tab(trade_mat, max_lines)
    trade_df <- merge(trade_df, centers, by.x = "from", by.y = "id")
    trade_df <- merge(trade_df, centers, by.x = "to", by.y = "id",
                      suffixes = c("", "_to"))
    # Crate trade lines data frame 
    routes <- gcIntermediate(trade_df[,c("long", "lat")], 
                             trade_df[,c("long_to", "lat_to")], pts_by_gc - 2, 
                             breakAtDateLine = TRUE, addStartEnd = TRUE, sp = TRUE)
    routes <- fortify(as(routes, "SpatialLinesDataFrame"))
    # Add trade values to lines
    routes$value <- rep(trade_df$value, each = pts_by_gc)
    routes
}


# This function takes a single simulation output and produces an animated map 
#  (nfr frames by iteration in the simulation and tfr seconds by frame)
#  where country colors represent cumulative supply decrease
#  and trade line colors represent change in exports along that route
#  for that iteration
sim_animation <- function(sim_out, nfr = 5, tfr = 0.5, max_lines = 50) {
    # Load world map
    data(wrld_simpl)
    centr <- coordinates(wrld_simpl)
    centr <- data.frame(id = rownames(centr), long = centr[, 1], lat = centr[, 2])
    world_df <- fortify(wrld_simpl)
    
    # Transform data and find bounds
    n_iter <- dim(sim_out$dR)[2]
    dS <- sim_out$dR + sim_out$dC
    dS <- t(apply(dS, 1, cumsum)) / 
                        (sim_out$P0 + colSums(sim_out$E0) - rowSums(sim_out$E0))
    dE <- sign(sim_out$dE) * log1p(abs(sim_out$dE))
    dSlim <- c(min(dS), 0)
    dElim <- c(min(dE), max(dE))
    
    # Great circle cut points for frames
    cut_pts <- seq(0, pts_by_gc, length.out = nfr)

    # Create animation    
    saveHTML({
        ani.options(interval = tfr)
        for (i in 1:(n_iter - 1)) {
            # Create geographical layers for export/import changes
            #   and cumulative decrease of countries' supplies at this iter.
            # Note: cap at 50 trade lines by iteration
            country_data <- add_country_data(world_df, dS[, i])
            trade_data <- get_trade_lines(dE[, , i], centr, max_lines)
            for (j in 1:nfr) {
                sim_map <- ggplot() + 
                    geom_polygon(aes(long, lat, group = group, fill = value),
                                 data = country_data, colour = "grey") +
                    scale_fill_gradient(name = "dS/S0", low = "red", high = "white",
                                        na.value = "grey", limits = dSlim) +
                    geom_path(aes(long, lat, group = group, colour = value), 
                              data = trade_data[trade_data$order <= cut_pts[j], ],
                              arrow = arrow(type = "closed", 
                                            length = unit(0.5, "lines"))) +
                    scale_color_gradient2(name = "log1p(dE)", na.value = "white",
                                          limits = dElim) +
                    scale_x_continuous(limits = c(-180, 180), expand = c(0, 0)) +
                    scale_y_continuous(limits = c(-90, 90), expand = c(0, 0)) +
                    coord_fixed() + theme_map()
            print(sim_map)
            ani.pause()
            }
        }
        # Plot final dS
        country_data <- add_country_data(world_df, dS[, n_iter])
        sim_map <- ggplot() + 
            geom_polygon(aes(long, lat, group = group, fill = value),
                         data = country_data, colour = "grey") +
            scale_fill_gradient(name = "dS/S0", low = "red", high = "white",
                                na.value = "grey", limits = dSlim) +
            scale_x_continuous(limits = c(-180, 180), expand = c(0, 0)) +
            scale_y_continuous(limits = c(-90, 90), expand = c(0, 0)) +
            coord_fixed() + theme_map() 
        print(sim_map)
        ani.pause()
    }, img.name = "sim_iter", single.opts = "utf8: false", autoplay = FALSE,
    interval = tfr, imgdir = "sim_maps", htmlfile = "sim_map.html",
    ani.height = 600, ani.width = 800, 
    title = "Simulation of cascading trade supply shock", description = "")
}
