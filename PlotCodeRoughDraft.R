data = read.csv("ghana_district_covariates.csv")

library(terra)
library(sf)
library(magrittr)

s = read_sf("geoBoundaries-GHA-ADM2.shp")
plot(s)

s %<>% arrange(shapeName)
data %<>% arrange(district_name)


s %<>% mutate(elevation = data$elev_mean)
s %<>% mutate(precipitation = data$precip_mean)
s %<>% mutate(poverty = data$poverty_H)
s %<>% mutate(literacy = data$literacy_rate)
s %<>% mutate(travel_time = data$travel_time_mean)
s %<>% mutate(population_density = data$pop_density)
s %<>% mutate(population = data$pop_total)


plot(s["elevation"])
plot(s["poverty"])
plot(s["precipitation"])
plot(s["literacy"])
plot(s["travel_time"])

xcel = read.csv("gha-rainfall-subnat-full.csv")
xcel %<>% filter(date == "2010-09-21", adm_level == 2) %>% group_by(PCODE) %>% arrange(PCODE)
xcel2 = xcel %>% select(PCODE, r3h)
write.csv(xcel2, "GHA2010Rain.csv")

tmap_mode("view")
tm_shape(s["precipitation"]) + 
  tm_polygons("precipitation", palette = "-RdYlGn") + 
  tm_layout(frame = FALSE, 
            outer.margins = 0,
            inner.margins = c(0.15, 0.0, 0.15, 0.3),
            legend.outside = FALSE,        
            legend.position = c("right", "center"),
            legend.bg.color = "white",)
