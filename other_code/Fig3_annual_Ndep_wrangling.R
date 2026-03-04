library(sf)

gdb_path <- "./data/Ndep_data_single_point.gdb"

# list layers
st_layers(dsn = gdb_path)

# read in
my_spatial_data <- st_read(dsn = gdb_path, layer = "T43649424010478")

length(c(1985:2021))
noxi <- c(my_spatial_data[7:43])
noxi_vec <- unlist(noxi)[1:37]
nred <- c(my_spatial_data[44:80])
nred_vec <- unlist(nred)[1:37]

ndep <- data.frame(names = names(noxi_vec),
                   noxi = noxi_vec,
                   nred = nred_vec)
rownames(ndep) <- NULL

ndep2 <- ndep %>%
  separate(names, into = c("junk", "year"), sep = "_") %>%
  select(-junk) %>%
  mutate(datetime = as.Date(paste0(year,"-07-01")),
         ndep = noxi + nred)
write.csv(ndep2, "./data/Fig3_annual_plot_Ndep.csv",row.names = FALSE)  
