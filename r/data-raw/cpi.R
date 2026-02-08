## code to prepare `cpi` dataset goes here

pacman::p_load(tidyverse, lubridate)

cpi = read_csv("https://fred.stlouisfed.org/graph/fredgraph.csv?bgcolor=%23e1e9f0&chart_type=line&drp=0&fo=open%20sans&graph_bgcolor=%23ffffff&height=450&mode=fred&recession_bars=on&txtcolor=%23444444&ts=12&tts=12&width=718&nt=0&thu=0&trc=0&show_legend=yes&show_axis_titles=yes&show_tooltip=yes&id=CPIAUCSL&scale=left&cosd=1947-01-01&coed=2023-11-01&line_color=%234572a7&link_values=false&line_style=solid&mark_type=none&mw=3&lw=2&ost=-99999&oet=99999&mma=0&fml=a&fq=Monthly&fam=avg&fgst=lin&fgsnd=2020-02-01&line_index=1&transformation=lin&vintage_date=2023-12-13&revision_date=2023-12-13&nd=1947-01-01")
cpi = cpi %>% rename(date = DATE, cpi = CPIAUCSL) %>% mutate(date = ymd(date))
cpi = cpi %>% group_by(year = year(date)) %>% summarise(cpi = mean(cpi, na.rm = T)) %>% ungroup()
try({saveRDS(cpi, "/Users/theo/Documents/dlm/data/cpi.RDS")})

usethis::use_data(cpi, overwrite = TRUE)
