# Modelos espacio-temporales para datos beta-distribuidos de Cropland.

# Cargar las librerías necesarias. ----

library(INLA)
library(fmesher)
library(inlabru)
library(ggplot2)
library(gridExtra)
library(viridis)
library(sf)
library(dplyr)
library(parallel)

# Algunas funciones de carácter general ----

colsc <- function(...) { # función para que las escalas entre distintas gráficas sean la misma
  scale_fill_gradientn(
    colours = turbo(n = 11),
    limits = range(..., na.rm = TRUE)
  )
}

combinatorialFunction <- function(list){ # Función para generar combinaciones
  DF.old <- list
  DF.new <- list
  for(i in 2:length(list)){
    DF.new[[i-1]] <- rep(DF.old[[i-1]], length(DF.old[[i]]))
    DF.new[[i]] <- rep(DF.old[[i]], each=length(DF.old[[i-1]]))
    DF.old <- DF.new
  }
  DF <- as.data.frame(DF.new)
  return(DF)
}

# El primer paso es cargar los datos y hacer una pequeña exploración de los mismos. ----

DF <- readRDS("./Data/NUTS3_sf_complete.RDS")

# Podemos realizar la representación gráfica de la variable para ver si varía mucho espacial y/o temporalmente
ggplot() + geom_sf(data = DF, mapping = aes(fill = Cropland)) + facet_wrap(facets = ~year) + theme_minimal()

logit <- function(x){log(x/(1-x))} # Vamos a evaluar visualmente la evolución temporal en la escala del predictor lineal del modelo
DF %>% select(which(names(DF) %in% c("Cropland"))) %>% st_drop_geometry(.) %>% 
  aggregate(x = ., by = list(DF$NUTS0, DF$year), FUN = "mean") %>% rename(., all_of(c(NUTS0 = "Group.1", year = "Group.2"))) %>% 
  ggplot(data = .) + geom_line(mapping = aes(x = year, y = logit(Cropland))) + facet_wrap(facets = ~NUTS0) + theme_bw()

# rep(DF$NUTS0 %>% table(.) %>% names(.), times = DF$NUTS0 %>% table(.) %>% as.numeric(.)/12)

DF %>% select(which(names(DF) %in% c("Cropland"))) %>% st_drop_geometry(.) %>% 
  aggregate(x = ., by = list(DF$NUTS3, DF$year), FUN = "mean") %>% rename(., all_of(c(NUTS3 = "Group.1", year = "Group.2"))) %>%  
  mutate(., NUTS0 = rep(rep(DF$NUTS0 %>% table(.) %>% names(.), times = DF$NUTS0 %>% table(.) %>% as.numeric(.)/length(unique(DF$year))), length(unique(DF$year)))) %>% 
  ggplot(data = .) + geom_line(mapping = aes(x = year, y = logit(Cropland), group = NUTS3)) + facet_wrap(facets = ~NUTS0) + theme_bw()

# Análisis de la correlación entre las covariables ----
## Vamos a eliminar las variables explicativas de latitude y longitude ya que incluirlas en el modelo supondría evaluar el efecto de las coordenadas x e y 
## conjuntamente a la estructura espacial, por lo que estarían explicando lo mismo
names_corr <- function(DFCovariates, corr_val){
  ind_lwtr <- as.vector(lower.tri(cor(DFCovariates)))
  correlation_values <- as.vector(cor(DFCovariates))[ind_lwtr]
  indx_corr_val <- abs(correlation_values)>corr_val
  DF_corr_cov <- data.frame(x=factor(rep(names(DFCovariates), times=ncol(DFCovariates)))[ind_lwtr],
                            y=factor(rep(names(DFCovariates), each=ncol(DFCovariates)))[ind_lwtr],
                            Corr=round(correlation_values,digits=3))
  CorrDF <- DF_corr_cov[indx_corr_val&DF_corr_cov[,1]!=DF_corr_cov[,2],]
  DFCovNames <- setdiff(names(DFCovariates), CorrDF$y)
  return(DFCovNames)
}
names_cov <- names_corr(DFCovariate = DF %>% select(-c(1:8, which(names(DF) %in% c("latitude", "longitude")))) %>% st_drop_geometry(.), corr_val = 0.75)
DF %>% select(all_of(names_cov)) %>% st_drop_geometry(.) %>% summary(.) # Summary de las covariables seleccionadas

# Transformación de las variables explicativas ----
DF[,names_cov] <- DF %>% select(all_of(names_cov)) %>% st_drop_geometry(.) %>% apply(X = ., MARGIN = 2, FUN = function(i){log((i-min(i))/sd(i) + 1)})

# Modelización de modelos espacio-temporales en áreas (besag, bym, lerroux; ar1). Selección de variables mediante métodos de regresión stepwise (both directions) ----

sp_priors <- c('besagproper', 'bym', 'besagproper2') # distribuciones previas para el efecto espacial: previas de Besag, Besag-York-Mollié, y Lerroux
t_prior <- c('ar1', 'rw1', 'rw2') # distribuciones previas para el efecto temporal: autorregresivo de orden 1, random walk de orden 1 y de orden 2

A <- list(sp = sp_priors, t = t_prior)
sp.t_models <- combinatorialFunction(list = A)

prior_pc.prec <- list(prec=list(prior="pc.prec", param = c(1, 0.01)))
prior_pc.cor <- list(prec=list(prior="pc.cor0", param = c(0.5, 0.9)))

## Emplearemos un método stepwise en ambas direcciones para la selección de variables que han de incluirse en el modelo final, evaluando su mejora según el criterio de WAIC. ----
## 1. Empezaremos con un conjunto vacío de variables explicativas.
## 2. Realizar un proceso stepwise en la dirección forward hasta que no puedan añadirse más variables explicativas.
## 3. Realizar un proceso stepwise en la dirección backward hasta que no puedan extraerse más variables.
## 4. Repetir los pasos 2 y 3 hasta que no se puedan implementar más cambios en el modelo para mejorar según el cirterio empleado.

# spdep::nb2INLA(file="./Data/NUTS3_nb", nb=NUTS3_nb) # escribimos el fichero para la matriz de adyacencia
g <- inla.read.graph(filename = "./Data/NUTS3_nb") # leemos el fichero con la matriz de adyacencia

cov_selected <- readRDS(file = "./Data/covselected.RDS")

for(i in 1:nrow(sp.t_models)){
  list_data_base <- list(y = DF$Cropland, intercept = rep(1, nrow(DF)), sp_idx = rep(1:g$n, times = length(unique(DF$year))))
  formula_inla_base = y ~ -1 + intercept + f(sp_idx, model = sp.t_models[i,1], graph = g, constr = TRUE, hyper = list(prior_pc.prec)) + f(temp_idx, model = sp.t_models[i,2], constr = TRUE, hyper = list(prior_pc.prec))
  
  list_data <- c(list_data_base, DF[,c(cov_selected)] %>% st_drop_geometry(.))
  formula_inla <- update.formula(old=formula_inla_base, new = paste("~ . ", paste(c(cov_selected), collapse = " + "), sep = " + ") %>% as.formula(.))
  
  model_beta <-
      inla(family = "beta", data = list_data, formula = formula_inla,
           control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = FALSE),
           inla.mode = "compact", safe = TRUE, verbose = FALSE)
}


# Modelo de dowscaling ----


Centroids <- DF[DF$year==unique(DF$year)[1],] %>% 
  st_centroid() %>% 
  st_geometry() 

Centroids %>% ggplot(.) + geom_sf()

NUTS3_sp <- as(DF[DF$year==unique(DF$year)[1],] %>% st_geometry(.), 'Spatial')

spde.puregrid <- as.matrix(expand.grid(x=seq(NUTS3_sp@bbox[1,1], NUTS3_sp@bbox[1,2],length.out=1E3),
                                       y=seq(NUTS3_sp@bbox[2,1], NUTS3_sp@bbox[2,2],length.out=1E3)))
spde.grid <- rbind(spde.puregrid, st_coordinates(Centroids))

Sp_spde.grid <- SpatialPoints(coords=spde.grid)
Sp_spde.grid@proj4string <- NUTS3_sp@proj4string
Spde.grid_over_indx <- over(x=Sp_spde.grid, y=NUTS3_sp)
spde.grid <- spde.grid[!is.na(Spde.grid_over_indx),]
Spde.grid_over_indx <- Spde.grid_over_indx[!is.na(Spde.grid_over_indx)]

max_edge <- quantile(as.vector(dist(as.matrix(st_multipoint(st_coordinates(st_convex_hull(st_union(Centroids)))))[,1:2])), p=c(0.03,0.1))
cutoff <- quantile(as.vector(dist(as.matrix(st_multipoint(st_coordinates(st_convex_hull(st_union(Centroids)))))[,1:2])), p=c(0.01))
mesh <- fm_mesh_2d_inla(loc=spde.grid, cutoff=cutoff, max.edge=max_edge,
                        boundary=list(fm_nonconvex_hull(x=Centroids, convex=-0.05), fm_nonconvex_hull(x=Centroids, convex=-0.2)))

Agrid <- inla.spde.make.A(mesh=mesh, loc=as.matrix(spde.grid))
Ablock_uni <- inla.spde.make.block.A(A=Agrid, block=as.vector(Spde.grid_over_indx), rescale="count")

Ablock_inf <- Ablock_uni
for(i in 2:length(unique(DF$year))){
  Ablock_inf <- rbind(Ablock_inf, Ablock_uni)
}

spde <- inla.spde2.pcmatern(mesh=mesh, alpha=2, prior.range=c(max_edge[1],0.1), prior.sigma=c(1,0.5), constr = TRUE)
spde.index <- inla.spde.make.index(name="spatial", n.spde=spde$n.spde)

list_nonsp <- list(intercept = rep(1, nrow(DF)), temp_idx = DF$year)
effects_list <- list(list(sp_idx = 1:mesh$n), list_nonsp)

inf_stk <- inla.stack(data = list(y = DF$Cropland),
                      A = list(Ablock_inf, 1),
                      effects = effects_list,
                      tag = "inf_stk")

for(i in 1:3){
  formula_inla_base = y ~ -1 + intercept + f(sp_idx, model = spde) + f(temp_idx, model = t_prior[l], constr = TRUE, hyper = list(prior_pc.prec))
  list_nonsp <- c(list(intercept = rep(1, nrow(DF)), temp_idx = DF$year), DF[,c(cov_selected)] %>% st_drop_geometry(.))
  
  effects_list <- list(list(sp_idx = 1:mesh$n), list_nonsp)
  
  inf_stk <- inla.stack(data = list(y = DF$Cropland),
                        A = list(Ablock_inf, 1),
                        effects = effects_list,
                        tag = "inf_stk")
  model_beta <-
    inla(family = "beta", data = inla.stack.data(inf_stk), formula = formula_inla,
         control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = FALSE),
         control.predictor = list(A = inla.stack.A(inf_stk)),
         inla.mode = "compact", 
         control.family = list(hyper = list(prior_pc.prec)), verbose = FALSE)
}
