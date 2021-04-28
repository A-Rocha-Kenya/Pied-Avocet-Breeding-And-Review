library(tidyverse)
library(RColorBrewer)
library(mgcv)
library(visreg)
library(gridExtra)
library(lubridate)
library(leaflet)
library(stringr)
library(readxl)
theme_set(theme_bw())
setwd('C:/Users/rnussba1/OneDrive/ARK/ARK - Science/01. Publications/Avocet_breeding_publication/figures')

#####
# Sabaki
d <- read.csv('../data/Sabaki_avocet_count.csv') %>% 
  filter(source=="WBC") %>% 
  filter(count!='X') %>% 
  mutate(
    julian = as.numeric(format(as.Date(date, format = "%d/%m/%Y"), "%j")),
    year = as.numeric(format(as.Date(date, format = "%d/%m/%Y"), "%Y")),
    month = as.numeric(format(as.Date(date, format = "%d/%m/%Y"), "%m")),
    count = as.numeric(paste(count)),
         )
dcyc <- 
  bind_rows(
    d,
    d %>% mutate(julian = julian+365),
    d %>% mutate(julian = julian-365),
  )

m <- gam(count ~ s(julian) + year, family="poisson", data=dcyc)

par(mfrow=c(2,1))
visreg(m)

yp = seq(2000,2020)
jp = seq(0,365)
dfp = data.frame(year=rep(yp,times=length(jp)),julian=rep(jp, each=length(yp)))
dfp <- dfp %>% mutate(
  fit.link = predict.gam(m,dfp, type = "link"),
  se.fit.link = predict.gam(m,dfp, type = "link", se.fit = T) %>% .$se.fit,
  fit = m$family$linkinv(fit.link),
  lwr =  m$family$linkinv(fit.link-3*se.fit.link),
  upr =  m$family$linkinv(fit.link+3*se.fit.link),
)


ggplot() +
  geom_histogram(data = d, aes(x=julian), breaks = c(as.numeric(format(ISOdate(2004,1:12,1),"%j")),366)) +
  geom_point(data = d %>% mutate(countt=ifelse(count>100,100,count)), aes( x=julian , y=countt, size=count )) + 
  #geom_smooth(aes( x=julian , y=count ), method = "gam", formula = y ~ s(x), method.args = list(family = "ziP")) +
  #geom_smooth(data = dcyc, aes( x=julian , y=count ), method = "gam", formula = y ~ s(x), method.args = list(family = "poisson")) +
  geom_ribbon(data = dfp %>% filter(year %in% 2020), aes(x=julian, ymin = lwr, ymax = upr), color = NA, alpha=0.3) +
  geom_line(data = dfp %>% filter(year %in% 2020), aes(x=julian, y = fit), size = 1) +
  xlab('Day of Year') +
  scale_x_continuous(breaks=as.numeric(format(ISOdate(2004,1:13,1),"%j")),
                     labels=format(ISOdate(2004,1:13,1),"%B"),
                     #limits = c(-250,400),
                     expand = c(0,0)) +
  scale_y_continuous(breaks= seq(0,100,20),
                     limits = c(0,100),
                     expand = c(0,0)) +
  theme(panel.grid.minor = element_blank(),aspect.ratio=9/16,legend.position="top") +
  ggsave(file="figure3.pdf")


d %>% ggplot() +
  geom_histogram(aes(x=year), binwidth = 1) +
  geom_point( aes( x=year , y=count, size=count )) + 
  #geom_smooth(aes( x=year , y=count ), method = "gam", formula = y ~ x, method.args = list(family = "poisson")) +
  geom_ribbon(data = dfp %>% filter(julian==which.max(dfp$fit[dfp$year==2020])), aes(x=year, ymin = lwr, ymax = upr), color = NA, alpha=0.3) +
  geom_line(data = dfp %>% filter(julian==which.max(dfp$fit[dfp$year==2020])), aes(x=year, y = fit), size = 1) +
  xlab('Year') +
  scale_x_continuous(
                     expand = c(0,0)) +
  scale_y_continuous(breaks= seq(0,120,20),
                     limits = c(0,120),
                     expand = c(0,0)) +
    theme(panel.grid.minor = element_blank(),aspect.ratio=9/16,legend.position="top") +
  ggsave(file="figure4.pdf")




#####
# eBird

mylist <- list()
u=0
for (f in list.files('../data/',pattern = "*linegraphs*")){
  eb <- readLines(paste0('../data/',f))
  u=u+1
  mylist[[u]] <- data.frame(
    freq = as.numeric(strsplit( eb[5],'\t')[[1]][-1][-1]),
    abun = as.numeric(strsplit( eb[9],'\t')[[1]][-1][-1]),
    hgct = as.numeric(strsplit( eb[17],'\t')[[1]][-1][-1]),
    totals = as.numeric(strsplit( eb[21],'\t')[[1]][-1][-1]),
    avgct = as.numeric(strsplit( eb[5],'\t')[[1]][-1][-1]),
    sample_size =as.numeric(strsplit( eb[6],'\t')[[1]][-1][-1]),
    date = head(seq(from=ymd("2020-01-01"), to = ymd("2021-01-01"),length.out=48+1),-1),
    cat = str_match(f,"_\\s*(.*?)\\s*_")[,2]
  )
}

eb <- do.call("rbind",mylist) %>% 
  mutate(
    cat = factor(cat, levels=c("AE,QA,KW","IL,EG,JO","ET,ER,DJ,SD","TZ,KE,UG","ZA,ZW,NA,MZ,BW"))
  )

eb %>% group_by(cat) %>% summarise(sample_size)

ebcyc = bind_rows(
  eb,
  eb %>% 
    filter(month(date)>8 | month(date)<4) %>% 
    mutate(
    date = if_else(month(date)>6,date-years(1),date+years(1))
  )
)


p<-ebcyc %>% 
  ggplot(aes(x=date, y=abun, color = cat)) +
    geom_point(shape=1)+
  geom_smooth(formula=y~x, method="loess",span=0.4, se = F) +
  scale_x_continuous(
    expand = c(0,0),
    # limits = c(eb$date[1], eb$date[nrow(eb)]),
    breaks=seq(from=ymd("2020-01-01"), to = ymd("2021-01-01"),by="month"),
    labels=format(seq(from=ymd("2020-01-01"), to = ymd("2021-01-01"),by="month"),"%B"),
  ) +
  scale_y_continuous(expand = c(0,0))
  # facet_grid(rows = vars(cat))
  

print(p)

ggsave(filename = 'eBird_abun.pdf', plot=p)



#####
# SABAP

sabap <- read.csv('../data/269.csv') %>%
  mutate(
    present = ifelse(Spp=='-',F,T),
    date = dmy(StartDate),
    week = week(date),
    lat = as.numeric(substr(Pentad, 1, 2)) + as.numeric(substr(Pentad, 3, 4))/60,
    lon = as.numeric(substr(Pentad, 6, 7)) + as.numeric(substr(Pentad, 8, 9))/60,
    letter = substr(Pentad, 5, 5),
    lat = ifelse((letter == '_' | letter == 'a'), -lat, lat),
    lon = ifelse((letter == 'a' | letter == 'b'), -lon, lon),
    cat = 'none',
    cat = ifelse(lat>-5 & lat <5.5 & lon>29& lon<41,'kenya',cat),
  )

lat_r = c(-14,-21,-26,-31,-36)
lon_r = c(10, 25, 37)
p <- list()
u=0
for (j in 1:(length(lon_r)-1)) {
  for (i in 1:(length(lat_r)-1)) {
    u=u+1;
    sabap <- sabap %>% mutate(
      cat = ifelse(lat<lat_r[i] & lat>lat_r[i+1] & lon>lon_r[j] & lon <lon_r[j+1],as.character(u),cat)
    )
  }
}


leaflet(data = sabap %>% filter(present) %>% group_by(lat,lon) %>% summarise(Pentad=Pentad, cat = first(cat)) %>% unique()) %>% 
  addTiles() %>%
  addMarkers(~lon, ~lat, label = ~as.character(cat), clusterOptions = markerClusterOptions(), group=~cat) %>% 
  addLayersControl(
    overlayGroups = sabap$cat %>% unique()
  )


sabap %>% ggplot(aes(x=lat)) +
  geom_histogram(binwidth = 1) +
  scale_x_continuous(
    limits = c(-36,-12),
    expand = c(0,0),
    breaks= seq(-40,0,1)
    ) 

p <- sabap %>% 
  group_by(cat,week) %>% 
  summarise(
    nb_card=n(),
    nb_sightings=sum(present),
  ) %>%
  mutate(week = ifelse(week<26,week+53,week)) %>% 
  ggplot(aes(x=week,y=nb_sightings/nb_card)) +
  geom_point() +
  geom_smooth(formula=y~s(x),method='gam')  +
  facet_grid(rows = vars(cat))

print(p)

ggsave(filename = 'sabap_zone.pdf', plot=p) 




mod <- gam(present ~ s(week) + lat, family = 'binomial', data= sabap %>% 
        filter(lat < -13) %>% 
        mutate(present = ifelse(present,1,0)))


sabap2 = sabap %>% 
  filter( !(cat %in% c('none'))) %>% 
  mutate(cat = ifelse(cat=='kenya','kenya','southernafrica')) %>% 
  group_by(week,cat) %>% 
  summarise(
    nb_card=n(),
    nb_sightings=sum(present), 
  )

sabap2cyc <- bind_rows(
  sabap2,
  sabap2 %>% mutate( week = week-53),
  sabap2 %>% mutate( week = week+53)
)

sabap2cyc = bind_rows(
  sabap2,
  sabap2 %>% 
    #filter(month(date)>8 | month(date)<4) %>% 
    mutate(
      week = if_else(week>26,week-53,week+53)
    )
)

p <- sabap2cyc %>%
  # mutate(week = ifelse(week<26,week+53,week)) %>% 
  ggplot(aes(x=week,y=nb_sightings/nb_card, color=cat)) +
  geom_point() +
  geom_smooth(formula=y~x,method='loess', span=0.4) +
  scale_x_continuous(
    expand = c(0,0),
    # limits = c(eb$date[1], eb$date[nrow(eb)]),
    breaks=c(1,54),
    labels=c('Jan','Jan'),
  ) +
  scale_y_continuous(expand = c(0,0))
 # facet_grid(rows = vars(cat))

print(p)

ggsave(filename = 'sabap_saea.pdf', plot=p) 
  







#####
# Count in South Africa
ctsa <- read.csv('../data/south_africa/CWAC_data_Avocet_selectionOfSItes.csv') %>% 
  mutate(dayofyear = yday(Date))


ctsa %>% 
  ggplot(aes(x=dayofyear,y=Count, color=Region)) +
  geom_point() +
  geom_smooth()








#####
# Count in East Africa

# This first dataset was gathered includeding various source (eBird, scopus, Don Turner pers com.)
ctea <- read.csv('../data/East_africa_avocet_count.csv') %>% 
  mutate(
    count = as.numeric(str_replace(count,',','')),
    date = dmy(date),
    year = year(date),
    month = month(date),
    season = ifelse(month %in% c(6,7,8),'summer','none'),
    season = ifelse(month %in% c(1,2,11,12,10),'winter',season),
    placeCT = ifelse(place %in% c('kenya','tanzania'),1,0)
  ) %>% arrange(date)

p <- ctea %>% 
  filter(placeCT==1) %>% 
  ggplot(aes(x=year, y=count, shape=season, color=place)) +
  geom_point(size=5) + 
  geom_line()+
  theme_minimal() + 
  scale_color_manual(values=c("#BB0000", "#1FB53A"))

ggsave(filename = 'East_africa_counts.pdf', plot=p) 




# We aggregated more data from East African bird count. 
ctea <- read_excel('../data/waterbirdcountNMK.xlsx', sheet = "Avocet January") 
# delete the average column and the row below the table
ctea <- ctea[1:26,1:33] 
# Main place
pla <- c('Lake Elmenteita','Lake Magadi','Lake Nakuru')
p <- ctea %>% 
  pivot_longer(starts_with('19')|starts_with('20'),names_to = "year",values_to = "count") %>% 
  mutate(
    year=as.numeric(year),
    Place = ifelse(Place %in% pla ,Place,'other')
    ) %>% 
  ggplot(aes(x=year, y=count, fill=Place)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  scale_x_continuous(minor_breaks = seq(1990, 2020, 1))

ggsave(filename = 'NMK_East_africa_counts_jan.pdf', plot=p) 


ctea <- read_excel('../data/waterbirdcountNMK.xlsx', sheet = "Avocet July") 
# delete the average column and the row below the table
ctea <- ctea[1:25,1:25] 
# Main place
p <- ctea %>% 
  pivot_longer(starts_with('19')|starts_with('20'),names_to = "year",values_to = "count") %>% 
  mutate(
    year=as.numeric(year),
    Place = ifelse(Place %in% pla ,Place,'other')
  ) %>% 
  ggplot(aes(x=year, y=count, fill=Place)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  scale_x_continuous(minor_breaks = seq(1990, 2020, 1))
  
ggsave(filename = 'NMK_East_africa_counts_jul.pdf', plot=p) 
