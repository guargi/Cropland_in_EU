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

combinatorialFunction <- function(list){
  DF.old <- list
  DF.new <- list
  for(i in 2:length(A)){
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

# Ahora correspondería hacer un análisis exploratorio ----
names_cov <-  names_cov[!names_cov %in% c("Urban_rent", "Grassland_rent", "Forest_rent", "Cropland_rent", "gva_A", "gva_A_sh")]
DF[,names_cov] <- DF %>% select(all_of(names_cov)) %>% st_drop_geometry(.) %>% apply(X = ., MARGIN = 2, FUN = function(i){log((i-min(i))/sd(i) + 1)})

# Modelización de modelos espacio-temporales en áreas (besag, bym, lerroux; ar1). Selección de variables mediante métodos de regresión stepwise (both directions) ----

sp_priors <- c('besagproper', 'bym', 'besagproper2') # distribución previa de Besag, distribución previa de Besag-York-Mollié, distribución previa de Lerroux
t_prior <- c('rw1', 'rw2','ar1')

# spdep::nb2INLA(file="./Data/NUTS3_nb", nb=NUTS3_nb) # escribimos el fichero para la matriz de adyacencia
g <- inla.read.graph(filename = "./Data/NUTS3_nb") # leemos el fichero con la matriz de adyacencia

# La desviación estándar marginal no devería ser superior a 10: 
# supongamos que comparamos valores de sigma = {5,10}, entonces exp(-5)*(1+exp(10))/(1+exp(5)) \approx 1 
prior_pc.prec <- list(prec=list(prior="pc.prec", param = c(10, 0.01))) 

list_data_base <- list(y = DF$Cropland, intercept = rep(1, nrow(DF)))
formula_inla_base = y ~ -1 + intercept
t1 <- Sys.time()
model_beta_base0 <- inla(family = "beta", data = list_data_base, formula = formula_inla_base,
                         control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = FALSE),
                         inla.mode = "compact", num.threads = 2,
                         verbose = FALSE)
t2 <- Sys.time()
difftime(t2,t1)
# model_beta_base0$waic$waic-model_beta_base_t$waic$waic
waic <- model_beta_base0$waic$waic

# Exploración echaustiva de efectos fijos ----
A <- list()
for(i in 1:length(names_cov)){
  A[[names_cov[i]]] <- c(FALSE, TRUE)
}

DFCovCombinatorialModels <- combinatorialFunction(list=A)

formula_inla_base = y ~ -1 + intercept
list_data_base <- list(y = DF$Cropland, intercept = rep(1, nrow(DF)))
list_res <- list()

t0 <- Sys.time()
ngroup <- 20
res_joint <- data.frame(waic=c(), dic=c(), lcpo=c())
for(j in 1:ceiling(nrow(DFCovCombinatorialModels)/ngroup)){ # round(nrow(DFCovCombinatorialModels)/ngroup)
  cat(sprintf("Group %i of %i.\n", j, ceiling(nrow(DFCovCombinatorialModels)/ngroup)))
  t1 <- Sys.time()
  res <- mclapply(X = (ngroup*(j-1)+1):min(c(ngroup*j, nrow(DFCovCombinatorialModels))), mc.cores = ngroup, FUN = function(i){
    list_data <- c(list_data_base, DF[,c(names_cov[which(DFCovCombinatorialModels[i,] %>% unlist(.) %>% as.vector(.))])] %>% st_drop_geometry(.))
    if(i==1){
      formula_inla <- formula_inla_base
    } else{
      formula_inla <- update.formula(old=formula_inla_base, new = paste("~ . ", paste(names_cov[which(DFCovCombinatorialModels[i,] %>% unlist(.) %>% as.vector(.))], collapse = " + "), sep = " + ") %>% as.formula(.))  
    }
    
    model_beta <- 
      try(
        inla(family = "beta", data = list_data, formula = formula_inla,
             control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = FALSE),
             inla.mode = "compact", num.threads = 20,
             # control.inla = list(strategy = "laplace", int.strategy = "auto"),
             control.inla = list(cmin = 0),
             # control.mode = list(theta = c(internal.hyperpar_config), x = c(0,modal_config[which(DFCovCombinatorialModels[i,] %>% unlist(.) %>% as.vector(.))]), fixed = FALSE, restart = FALSE),
             verbose = FALSE),
        silent = TRUE)
    if(class(model_beta)=="try-error"){
      return(data.frame(waic=NA, dic=NA, lcpo=NA))
    } else{
      return(data.frame(waic=model_beta$waic$waic, dic=model_beta$dic$dic, lcpo=-sum(log(model_beta$cpo$cpo))))
    }
  }) %>% do.call(what = rbind, .)
  res_joint <- rbind(res_joint, res)
  
  # saveRDS(object = res_joint, file = "./Results/res_joint_model_imp.RDS")
  # saveRDS(object = j, file  = "./Results/step_model_imp.RDS")
  t2 <- Sys.time()
  print(difftime(t2,t1))
}
t3 <- Sys.time()
difftime(t3,t0)

