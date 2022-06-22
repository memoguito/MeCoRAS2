# Git ----
usethis::use_git() # repository
usethis::use_github() # push
# Libraries ----
library(data.table)
library(readxl)
library(pbapply)
library('tidyverse')
library('devtools')
library(ggplot2)
library("leaflet")
library("magrittr")
library("leaflet.extras")
library(spatstat)
library(RColorBrewer)
library(geog4ga3)
library(pbapply)
library(MASS)
library(ComplexHeatmap)

# script patient 12 ----
data <- lapply(list.files(path="../Data/", pattern = ".xlsx", recursive = F), function(x) openxlsx::read.xlsx(paste0("../Data/", x)))
names(data) <- sapply(strsplit(list.files(path="../Data/", pattern = ".xlsx", recursive = F), split=".xlsx"), function(x) x[1])
jpk <- data.frame(openxlsx::read.xlsx("../JPK_pos.xlsx"))

for ( i in names(data)) {
  data[[i]]$x <- as.numeric(data[[i]]$X.Position)*1000000 + jpk[jpk$file==i,"start.X"] - jpk[jpk$file==i, "origin.X"]
  data[[i]]$y <- as.numeric(data[[i]]$Y.Position)*1000000 + jpk[jpk$file==i,"start.Y"] - jpk[jpk$file==i, "origin.Y"]
  data[[i]]$Ym <- as.numeric(data[[i]]$`Young's.Modulus.[Pa]`)/1000
}

df <- do.call(rbind, data)
df.info <- t(sapply(strsplit(rownames(df), split="_"), 
                  function(x) c(patient=x[1],
                                tissue=x[2],
                                slice=x[3],
                                map=paste(unlist(strsplit(x[4], split=".", fixed=T))[1:3], collapse="."))))

df <- data.frame(df, df.info)

p <- ggplot(data=df) + 
  geom_point(aes(x=x, y=y, color=sqrt(Ym)), shape = 15, size=5) +
  scale_color_scico(palette="batlow") +
  xlim(c(550,850)) +
  ylim(c(-250, 50))

plot(p)

df.sf <- subset(df) %>% st_as_sf(coords=c("x", "y"))
vpolyg <- do.call(c, st_geometry(df.sf)) %>%
    st_voronoi() %>%
    st_collection_extract()
df.v <- df.sf
df.v$geometry <- vpolyg[unlist(st_intersects(df.sf, vpolyg))]
df.v <- df.v %>%
  st_intersection(st_polygon(list(rbind(c(550,-250), c(850,-250), c(850,50), c(550,50), c(550,-250)))))


p <- ggplot(df %>%
              st_as_sf(coords=c("x","y"))) +
  geom_sf(data=df.sf, aes(fill=sqrt(Ym)), size=3) +
  scale_fill_scico(palette="batlow")

plot(p)


samples$pat17.v <- do.call(rbind, lapply(unique(samples$pat17$type), function(s) {
  pat.sf <- subset(samples$pat17, type==s) %>%
    st_as_sf(coords=c("x", "y"))
  vpolyg <- do.call(c, st_geometry(pat.sf)) %>%
    st_voronoi() %>%
    st_collection_extract()
  pat.v <- pat.sf
  pat.v$geometry <- vpolyg[unlist(st_intersects(pat.sf, vpolyg))]
  # W.bbox <- st_polygon(list(rbind(c(-60,-60),c(60,-60),c(60, 60),c(-60, 60),c(0-60,-60))))
  W.bbox <- st_polygon(list(rbind(c(-30,-30),c(30,-30),c(30, 30),c(-30, 30),c(-30,-30))))
  W.owin <- as.owin(W.bbox)
  pat.v <- pat.v %>%
    st_intersection(W.bbox)
  return(pat.v)
}))

samples$pat17.idw <- lapply(unique(samples$pat17$type), function(s) {
  pat.ppp <- as.ppp(X=samples$pat17[samples$pat17$type==s&!is.na(samples$pat17$z),1:4], W=owin(xrange=c(-30,30),yrange=c(-30,30)))
  z_p.idw <- spatstat::idw(pat.ppp, power = 2)
  idw_df <- data.frame(expand.grid(x= z_p.idw$z$xcol, y = z_p.idw$z$yrow),
                       z = as.vector(t(z_p.idw$z$v)), type=s,
                       tissue=sapply(strsplit(s, split=""), function(x) paste0(x[x%in%letters], collapse="")))
  return(idw_df)
})
samples$pat17.idw <- do.call(rbind, samples$pat17.idw)

ggplot(samples$pat17.idw) + 
  geom_tile(aes(x, y, fill=log(z))) +
  geom_sf(data = samples$pat17.v,
          color = adjustcolor("white", alpha.f = .25), fill = NA) +
  scale_fill_distiller(palette="Spectral") +
  # scale_fill_gradientn(colours = rev(brewer.pal(n = 7, "Spectral")), values=rescale(c(0,50,100,200,500,1500,3000))) +
  facet_wrap(~factor(type, levels = c("epithelium", "cancer1", "cancer2", "stroma1")),
             dir="v", nrow=2)