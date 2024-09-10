library(sf)
library(ggplot2)
library(ggspatial)
china_map <- sf::st_read("https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json") 
class(china_map)

###
main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=0.5, colour="NA"),
                   axis.line.y=element_line(size=0.5, colour="NA"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=11),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=9),
                   legend.title = element_text(size=10),
                   plot.tag = element_text(size = 14, face = "bold"))
library(ggplot2)
mytheme = theme( panel.background = element_rect(fill='white', colour='black'),
                 legend.position = "none",
                 panel.grid=element_blank(), 
                 legend.title = element_text(size = 9),
                 legend.text = element_text(size = 8),
                 legend.background = element_rect(fill = NA), #axis.ticks.length = unit(0.4,"lines"), 
                 axis.ticks = element_line(color='black'),
                 axis.line = element_line(colour = "black"), 
                 axis.title.x = element_text(colour='black', size=13),
                 axis.title.y = element_text(colour='black', size=13),
                 axis.text = element_text(colour='black',size=11),
                 plot.tag = element_text(size = 14, face = "bold")) 

### 
datasel <- data.frame(province = c("Guangxi","Guangdong","Hunan","Hubei","Henan","Shandong"),
                      lon = c(108.79440,113.42992,111.71165,112.27130,113.61972,118.18776),
                      lat = c(23.83338,23.33464,27.62922,30.98753,33.90265,36.37609))
library(plyr)
province_mapping <- c("广西壮族自治区" = "Guangxi", "广东省" = "Guangdong", "湖南省" = "Hunan", "湖北省" = "Hubei","河南省" = "Henan", "山东省" = "Shandong")
china_map$name <- province_mapping[china_map$name]

ggplot(china_map)+
  geom_sf(color='black',fill=NA,size=0.8)+
  geom_sf(data = china_map[c(15:20),],size=1,aes(fill=factor(name, levels = rev(c("Guangdong","Guangxi","Hunan","Hubei","Henan","Shandong"))))) + 
  #geom_text(data=datasel,aes(x=lon, y=lat-0.6 ,label=province),size=3,colour="black")+
  geom_point(data = final_site_coordinates, aes(x=Longitude,y=Latitude),
             size=2,alpha=1, shape = 16, color = "black") + 
  main_theme+
  scale_fill_manual(values=c("Shandong" = "#376694", "Henan" = "#73BCD5", "Hubei" = "#ABDCE0",
                             "Hunan" = "#EFBA55","Guangxi" = "#E68C51","Guangdong" = "#994240"))+
  scale_color_manual(name = "Site", values = c("black","black","black","black","black","black"))+
  annotation_scale(location = "br", style = "ticks",line_width = 0.1,pad_y = unit(0.5, "cm"),text_cex = 0.9) + 
  annotation_north_arrow(location = "tl", which_north = F, pad_x = unit(0.05, "in"), pad_y = unit(0.05, "in"),style = north_arrow_minimal)+
  #coord_sf(crs = "+proj=laea +lat_0=40 +lon_0=104")+
  guides(fill = guide_legend(title = "Province",ncol = 1, override.aes = list(shape = 22, size=5))) +
  theme(text = element_text(size = 13),
        legend.position = c(0.22,0.26),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(color = "grey",size=0.2),
        axis.line = element_blank())+
  labs(x='', y='', tag = "(a)") -> p1a; p1a

ggplot(china_map)+
  geom_sf(data = china_map,fill="grey95",size=1) + 
  xlim(105,122)+ ylim(20,40)+ 
  geom_sf(data = china_map[c(15:20),],size=1,aes(fill=factor(name, levels = rev(c("Guangdong","Guangxi","Hunan","Hubei","Henan","Shandong"))))) + 
  #geom_text(data=datasel,aes(x=lon, y=lat-0.6 ,label=province),size=3,colour="black")+
  geom_point(data = final_site_coordinates, aes(x=Longitude,y=Latitude),
             size=1,alpha=1, shape = 16, color = "grey20") + 
  main_theme+
  scale_fill_manual(values=c("Shandong" = "#376694", "Henan" = "#73BCD5", "Hubei" = "#ABDCE0",
                             "Hunan" = "#EFBA55","Guangxi" = "#E68C51","Guangdong" = "#994240"))+
  scale_color_manual(name = "Site", values = c("black","black","black","black","black","black"))+
  annotation_scale(location = "br", style = "ticks",line_width = 0.1,pad_y = unit(0.5, "cm"), text_cex = 0.9) + 
  annotation_north_arrow(location = "tl", which_north = "true", height = unit(1, "cm"),width = unit(1, "cm"),pad_x = unit(0.3, "cm"), pad_y = unit(0.3, "cm"),style = north_arrow_fancy_orienteering) +
  #coord_sf(crs = "+proj=laea +lat_0=40 +lon_0=104")+
  guides(fill = guide_legend(title = "Province",ncol = 1, nrow = 6, override.aes = list(shape = 22, size=5))) +
  theme(text = element_text(size = 13),
        panel.background = element_rect(fill='white', colour='black'),
        legend.position = "right",
        panel.grid.major = element_line(color = "grey",size=0.2))+
  labs(x='', y='') -> p1a; p1a

##### 
library(xml2)
library(rvest)
library(dplyr)
library(stringr)
library(rjson)
library(jsonlite)

gGetLocation = function(address)  {
  key = "###########################" ### Amap API-key
  url = paste0("https://restapi.amap.com/v3/geocode/geo?key=",key,"&address=",address)
  data = read_html(url, encoding='utf-8') %>% html_text()
  df = as.data.frame(fromJSON(data))  
  return (df['geocodes.location']) 
}

site_location = rbind(city1,city2,city3,city4,city5,city6)
#### 
final_site_coordinates = NULL
for (i in 1:nrow(site_location)) {
  site <- site_location$city[i]
  tryCatch({
    coordinates <- as.character(gGetLocation(site))
    coords_split <- strsplit(coordinates, ",")[[1]]
    longitude <- as.numeric(coords_split[1])
    latitude <- as.numeric(coords_split[2])
  }, error = function(e) {
    message("Error occurred for site: ", site, ". Setting coordinates to NA.")
    longitude <- NA
    latitude <- NA
  })
  site_coordinates <- data.frame(City = site, Longitude = longitude, Latitude = latitude)
  final_site_coordinates <- rbind(final_site_coordinates, site_coordinates)
  print(paste0("slope", i))
}

head(final_site_coordinates)
