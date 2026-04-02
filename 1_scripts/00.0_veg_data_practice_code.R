# Code for clipping data from the backfilled vegetation layer to the polygon grid using PostGIS 

library(sf)
library(terra)
library(tidyverse)
library(ggplot2)
library(RPostgreSQL)


# Note: Once you enter these, they will be visible in your R environment, so I remove them after connecting to the database
username <- rstudioapi::askForPassword("Database username")
password <- rstudioapi::askForPassword("Database password")


pg = dbDriver("PostgreSQL")
# Local Postgres.app database; no password by default
# Of course, you fill in your own database information here.
con = dbConnect(pg, user = username, password = password,
                host = "localhost", port = 5432, dbname = "abc_program_data"); rm(list = c("username", "password"))


## ------------------------------------------------------------------------------------------------------
## Use SQL to clip the backfilled vegetation data to each grid polygon 
## ------------------------------------------------------------------------------------------------------

system.time({
  b <- st_read(con, query = "SELECT 
             c.fid AS clip_id, 
             d.*,
             ST_Intersection(d.geom, c.geometry)
             FROM osr_g1_grid c
             JOIN bfveg_2021_g2 d
             ON
             ST_Intersects(d.geom, c.geometry)
             LIMIT 10;")
})

write.csv(b, "0_data/processed/bfveg_fmu_group_grid/test_veg.csv")





# Takes ~ 2 hours

c <- st_read(con, query = "SELECT DISTINCT combined_chgbycwcs FROM backfill_veghf_2021;")

st_read(con, query = "SELECT * FROM osr_g1_grid LIMIT 10;")
