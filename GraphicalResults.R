# Representación gráfica de los resultados 

# Librerías ----

library(INLA)
library(fmesher)
library(inlabru)
library(ggplot2)
library(gridExtra)
library(viridis)
library(sf)
library(dplyr)
library(stringr)

# Algunas funciones ----

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

# Carga de datos ----

DF <- readRDS("./Data/NUTS3_sf_complete.RDS")
mesh <- readRDS(file = "./Results/mesh.RDS")
spde.grid_puregrid <- readRDS(file = "./Results/spde.grid_puregrid.RDS")

besagproper_model <- readRDS(file = "./Results/model_beta_spt1.RDS")
bym_model <- readRDS(file = "./Results/model_beta_spt2.RDS")
leroux_model <- readRDS(file = "./Results/model_beta_spt3.RDS")

besagproper_model_ar1 <- readRDS(file = "./Results/model_beta_spt1.RDS")
besagproper_model_rw1 <- readRDS(file = "./Results/model_beta_spt4.RDS")
besagproper_model_rw2 <- readRDS(file = "./Results/model_beta_spt7.RDS")

model_beta_leroux.ar1 <- readRDS(file = "./Results/model_beta_spt4.RDS")
model_beta_besagproper.ar1 <- readRDS(file = "./Results/model_beta_spt1.RDS")

model_beta_downscaling.ar1 <- readRDS(file = "./Results/model_beta_downscaling.ar1.RDS")
model_beta_downscaling.rw1 <- readRDS(file = "./Results/model_beta_downscaling.rw1.RDS")
model_beta_downscaling.rw2 <- readRDS(file = "./Results/model_beta_downscaling.rw2.RDS")

# Graficando resultados ----

## Gráfica de la malla ----
ggplot() + 
  gg(data = mesh) + 
  theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

## Gráfica de los puntos de integración ----
ggplot() + 
  gg(data = mesh) + 
  geom_point(data = data.frame(id = "Integration points", spde.grid_puregrid), mapping = aes(x = x, y = y)) +
  facet_wrap(facets = ~id) + 
  theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

## Gráficas de los efectos epsaciales de los modelos de áreas ----
ggplot() + 
  geom_sf(data = st_sf(data.frame(besagproper_model$summary.random$sp_idx, id = "Spatial effect (Besag)"), geometry = DF[DF$year==unique(DF$year)[1],] %>% st_geometry(.)), mapping = aes(fill = mean)) +
  geom_sf(data = st_sf(data.frame(bym_model$summary.random$sp_idx[1:g$n,], id = "Spatial effect (BYM)"), geometry = DF[DF$year==unique(DF$year)[1],] %>% st_geometry(.)), mapping = aes(fill = mean)) +
  geom_sf(data = st_sf(data.frame(leroux_model$summary.random$sp_idx, id = "Spatial effect (Leroux)"), geometry = DF[DF$year==unique(DF$year)[1],] %>% st_geometry(.)), mapping = aes(fill = mean)) +
  scale_fill_viridis_c(option = "turbo") + theme_bw() + facet_wrap(facets = ~ id, ncol = 3) + theme(strip.text = element_text(size=14, face = "bold"))

## Gráficas de los efectos temporales ----
ggplot() + 
  geom_ribbon(data = data.frame(besagproper_model_ar1$summary.random$temp_idx, id = "Area model temporal effect (ar1)") %>% rename(., all_of(c(q1='X0.025quant', q3='X0.975quant'))), mapping = aes(x = ID, ymin = q1, ymax = q3), alpha=0.4) +
  geom_line(data = data.frame(besagproper_model_ar1$summary.random$temp_idx, id = "Area model temporal effect (ar1)"), mapping = aes(x=ID, y=mean)) +
  geom_ribbon(data = data.frame(besagproper_model_rw1$summary.random$temp_idx, id = "Area model temporal effect (rw1)") %>% rename(., all_of(c(q1='X0.025quant', q3='X0.975quant'))), mapping = aes(x = ID, ymin = q1, ymax = q3), fill = "red", alpha=0.4) +
  geom_line(data = data.frame(besagproper_model_rw1$summary.random$temp_idx, id = "Area model temporal effect (rw1)"), mapping = aes(x=ID, y=mean), color = "red") +
  geom_ribbon(data = data.frame(besagproper_model_rw2$summary.random$temp_idx, id = "Area model temporal effect (rw2)") %>% rename(., all_of(c(q1='X0.025quant', q3='X0.975quant'))), mapping = aes(x = ID, ymin = q1, ymax = q3), fill = "blue", alpha=0.4) +
  geom_line(data = data.frame(besagproper_model_rw2$summary.random$temp_idx, id = "Area model temporal effect (rw2)"), mapping = aes(x=ID, y=mean), color = "blue") +
  theme_bw() + facet_wrap(facets = ~ id, scales = "free", ncol = 1) +  xlab(label = "Year") + ylab(label = "Values") +
  theme(strip.text = element_text(size=14, face = "bold"))

## Gráficas distribuciones efectos fijos e hiperparámetros ----
ggplot() + # Marginales de los efectos fijos
  geom_line(data = bind_rows(lapply(X = besagproper_model_ar1$marginals.fixed, FUN = function(x){as.data.frame(x)}), .id = "id"), mapping = aes(x = x, y = y)) + 
  theme_bw() + facet_wrap(facets = ~ id, ncol = 4, scales = "free") + theme(strip.text = element_text(size=14, face = "bold"))

ggplot() + # Marginales de los hiperparámetros
  geom_line(data = bind_rows(lapply(X = besagproper_model_rw1$marginals.hyperpar, FUN = function(x){as.data.frame(x)}), .id = "id"), mapping = aes(x = x, y = y)) + 
  theme_bw() + facet_wrap(facets = ~ id, ncol = 2, scales = "free") + theme(strip.text = element_text(size=14, face = "bold"))

# Gráficas residuos ----
residual_areamodel <- DF$Cropland - model_beta_leroux.ar1$summary.fitted.values[model_beta_leroux.ar1$summary.fitted.values %>% rownames(.) %>% str_detect(string = ., pattern = "fitted.Predictor."), "mean"]
residual_downscaling <- DF$Cropland - model_beta_downscaling.ar1$summary.fitted.values[model_beta_downscaling.ar1$summary.fitted.values %>% rownames(.) %>% str_detect(string = ., pattern = "fitted.APredictor."), "mean"]

DFresidual_areamodel <- data.frame(NUTS3 = DF$NUTS3, year = DF$year, res = residual_areamodel) %>% st_sf(., geometry = st_geometry(DF))
DFresidual_downscaling <- data.frame(NUTS3 = DF$NUTS3, year = DF$year, res = residual_downscaling) %>% st_sf(., geometry = st_geometry(DF))

## Distribución de los valores de los reiduos, sin agregación ----
ggplot() + 
  geom_density(data = data.frame(x = residual_arealmodel, id = "Area model"), mapping = aes(x = x), color = "red") +
  geom_density(data = data.frame(x = residual_downscaling, id = "Downscaling model"), mapping = aes(x = x), color = "blue") +
  theme_bw() + facet_wrap(facets = ~ id, nrow = 1, scales = "free") + theme(strip.text = element_text(size=14, face = "bold"))

## Residuos agregados temporalmente ----
ggplot() +
  geom_line(data = data.frame(DFresidual_areamodel$res %>% st_drop_geometry(.) %>% aggregate(x = ., by = list(DFresidual_areamodel$year), FUN = "sum"), id = "Area model residuals (spatial aggreation)"), mapping = aes(y = x, x = Group.1)) +
  geom_line(data = data.frame(DFresidual_downscaling$res %>% st_drop_geometry(.) %>% aggregate(x = ., by = list(DFresidual_downscaling$year), FUN = "sum"), id = "Downscaling model residuals (spatial aggreation)"), mapping = aes(y = x, x = Group.1)) +
  theme_bw() + facet_wrap(facets = ~ id, nrow = 1) + theme(strip.text = element_text(size=14, face = "bold"))

## Residuos agregados espacialmente, según NUTS3
ggplot() +
  geom_sf(data = st_sf(DFresidual_areamodel$res %>% st_drop_geometry(.) %>% aggregate(x = ., by = list(DFresidual_areamodel$NUTS3), FUN = "sum"),
                       geometry = st_geometry(DF[DF$year==unique(DF$year)[1],]), id = "Area model residuals (temporal aggreation)"), 
          mapping = aes(fill = x)) +
  geom_sf(data = st_sf(DFresidual_downscaling$res %>% st_drop_geometry(.) %>% aggregate(x = ., by = list(DFresidual_downscaling$NUTS3), FUN = "sum"),
                       geometry = st_geometry(DF[DF$year==unique(DF$year)[1],]), id = "Downscaling model residuals (temporal aggreation)"), 
          mapping = aes(fill = x)) + labs(fill = "Values") +
  scale_fill_viridis_c(option = "turbo") + theme_bw() + facet_wrap(facets = ~ id, nrow = 1) + theme(strip.text = element_text(size=14, face = "bold"))

## Residuos sin agregación, a lo largo de todos los nodos temporales y espaciales
ggplot() +
  geom_sf(data = st_sf(DFresidual_areamodel %>% st_drop_geometry(.), id = "Area model residuals", geometry = st_geometry(DF)), mapping = aes(fill = res)) +
  geom_sf(data = st_sf(DFresidual_downscaling %>% st_drop_geometry(.), id = "Downscaling model residuals", geometry = st_geometry(DF)), mapping = aes(fill = res)) +
  scale_fill_viridis_c(option = "turbo") + theme_bw() + facet_wrap(facets = ~ id + year, nrow = 4)

## Residuos sin agregación, a lo largo de todos los nodos temporales y espaciales (sólo para el modelo de áreas)
ggplot() +
  geom_sf(data = st_sf(DFresidual_areamodel %>% st_drop_geometry(.), id = "Area model residuals", geometry = st_geometry(DF)), mapping = aes(fill = res)) +
  scale_fill_viridis_c(option = "turbo") + theme_bw() + facet_wrap(facets = ~ id + year, nrow = 3)

## Residuos sin agregación, a lo largo de todos los nodos temporales y espaciales (sólo para el modelo de downscaling)
ggplot() +
  geom_sf(data = st_sf(DFresidual_downscaling %>% st_drop_geometry(.), id = "Downscaling model residuals", geometry = st_geometry(DF)), mapping = aes(fill = res)) +
  scale_fill_viridis_c(option = "turbo") + theme_bw() + facet_wrap(facets = ~ id + year, nrow = 3)

DFresidual_areamodel$res %>% st_drop_geometry(.) %>% aggregate(x = ., by = list(DFresidual_areamodel$NUTS3), FUN = "sum")

# Resultados modelos de downscaling ----

ggplot() + 
  geom_ribbon(data = model_beta_downscaling.rw2$summary.random$temp_idx %>% rename(., all_of(c(q1='0.025quant', q3='0.975quant'))), mapping = aes(x = ID, ymin = q1, ymax = q3), alpha=0.4) +
  geom_line(data = model_beta_downscaling.rw2$summary.random$temp_idx, mapping = aes(x=ID, y=mean)) + 
  geom_ribbon(data = model_beta_downscaling.rw1$summary.random$temp_idx %>% rename(., all_of(c(q1='0.025quant', q3='0.975quant'))), mapping = aes(x = ID, ymin = q1, ymax = q3), fill = "red", alpha=0.4) +
  geom_line(data = model_beta_downscaling.rw1$summary.random$temp_idx, mapping = aes(x=ID, y=mean), color = "red") +
  geom_ribbon(data = model_beta_downscaling.ar1$summary.random$temp_idx %>% rename(., all_of(c(q1='0.025quant', q3='0.975quant'))), mapping = aes(x = ID, ymin = q1, ymax = q3), fill = "blue", alpha=0.4) +
  geom_line(data = model_beta_downscaling.ar1$summary.random$temp_idx, mapping = aes(x=ID, y=mean), color = "blue") + theme_bw()

gg_spde.mean <- ggplot() + geom_tile(data = data.frame(spde.grid_puregrid,
                                                       z=drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.ar1$summary.random$sp_idx$mean)),
                                     mapping = aes(x=x, y=y, fill=z)) + scale_fill_viridis_c(option = "turbo") + theme_minimal()

gg_spde.sd <- ggplot() + geom_tile(data = data.frame(spde.grid_puregrid,
                                                     z=drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.ar1$summary.random$sp_idx$sd)),
                                   mapping = aes(x=x, y=y, fill=z)) + scale_fill_viridis_c(option = "turbo") + theme_minimal()

gg_spde.q1 <- ggplot() + geom_tile(data = data.frame(spde.grid_puregrid,
                                                     z=drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.ar1$summary.random$sp_idx$`0.025quant`)),
                                   mapping = aes(x=x, y=y, fill=z)) + scale_fill_viridis_c(option = "turbo") + theme_minimal()

gg_spde.q3 <- ggplot() + geom_tile(data = data.frame(spde.grid_puregrid, z=drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.ar1$summary.random$sp_idx$`0.975quant`)),
                                   mapping = aes(x=x, y=y, fill=z)) + scale_fill_viridis_c(option = "turbo") + theme_minimal()

csc <- colsc(drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.ar1$summary.random$sp_idx$mean),
             drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.ar1$summary.random$sp_idx$sd),
             drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.ar1$summary.random$sp_idx$`0.025quant`),
             drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.ar1$summary.random$sp_idx$`0.975quant`))

grid.arrange(arrangeGrob(grobs = list(gg_spde.mean, gg_spde.sd, gg_spde.q1, gg_spde.q3), nrow = 1))

ggplot() + 
  geom_tile(data = data.frame(spde.grid_puregrid, id = "A.1 Mean (ar1)", z=drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.ar1$summary.random$sp_idx$mean)), mapping = aes(x=x, y=y, fill=z)) +
  geom_tile(data = data.frame(spde.grid_puregrid, id = "A.2 Stdev. (ar1)", z=drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.ar1$summary.random$sp_idx$sd)), mapping = aes(x=x, y=y, fill=z)) +
  geom_tile(data = data.frame(spde.grid_puregrid, id = "A.3 Quantile 0.025 (ar1)", z=drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.ar1$summary.random$sp_idx$`0.025quant`)), mapping = aes(x=x, y=y, fill=z)) +
  geom_tile(data = data.frame(spde.grid_puregrid, id = "A.4 Quantile 0.975 (ar1)", z=drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.ar1$summary.random$sp_idx$`0.975quant`)), mapping = aes(x=x, y=y, fill=z)) +
  
  geom_tile(data = data.frame(spde.grid_puregrid, id = "B.1 Mean (rw1)", z=drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.rw1$summary.random$sp_idx$mean)), mapping = aes(x=x, y=y, fill=z)) +
  geom_tile(data = data.frame(spde.grid_puregrid, id = "B.2 Stdev. (rw1)", z=drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.rw1$summary.random$sp_idx$sd)), mapping = aes(x=x, y=y, fill=z)) +
  geom_tile(data = data.frame(spde.grid_puregrid, id = "B.3 Quantile 0.025 (rw1)", z=drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.rw1$summary.random$sp_idx$`0.025quant`)), mapping = aes(x=x, y=y, fill=z)) +
  geom_tile(data = data.frame(spde.grid_puregrid, id = "B.4 Quantile 0.975 (rw1)", z=drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.rw1$summary.random$sp_idx$`0.975quant`)), mapping = aes(x=x, y=y, fill=z)) +
  
  geom_tile(data = data.frame(spde.grid_puregrid, id = "C.1 Mean (rw2)", z=drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.rw2$summary.random$sp_idx$mean)), mapping = aes(x=x, y=y, fill=z)) +
  geom_tile(data = data.frame(spde.grid_puregrid, id = "C.2 Stdev. (rw2)", z=drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.rw2$summary.random$sp_idx$sd)), mapping = aes(x=x, y=y, fill=z)) +
  geom_tile(data = data.frame(spde.grid_puregrid, id = "C.3 Quantile 0.025 (rw2)", z=drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.rw2$summary.random$sp_idx$`0.025quant`)), mapping = aes(x=x, y=y, fill=z)) +
  geom_tile(data = data.frame(spde.grid_puregrid, id = "C.4 Quantile 0.975 (rw2)", z=drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.rw2$summary.random$sp_idx$`0.975quant`)), mapping = aes(x=x, y=y, fill=z)) +
  
  scale_fill_viridis_c(option = "turbo") + theme_bw() + facet_wrap(facets = ~ id, nrow = 3) + labs(fill = expression(u(s))) + 
  theme(strip.text = element_text(size=14, face = "bold"), 
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()
  )

ggplot() + 
  geom_tile(data = data.frame(spde.grid_puregrid, id = "Mean", z=drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.ar1$summary.random$sp_idx$mean)), mapping = aes(x=x, y=y, fill=z)) +
  geom_tile(data = data.frame(spde.grid_puregrid, id = "Stdev.", z=drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.ar1$summary.random$sp_idx$sd)), mapping = aes(x=x, y=y, fill=z)) +
  geom_tile(data = data.frame(spde.grid_puregrid, id = "Quantile 0.025", z=drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.ar1$summary.random$sp_idx$`0.025quant`)), mapping = aes(x=x, y=y, fill=z)) +
  geom_tile(data = data.frame(spde.grid_puregrid, id = "Quantile 0.975", z=drop(fm_basis(x=mesh, loc=spde.grid_puregrid) %*% model_beta_downscaling.ar1$summary.random$sp_idx$`0.975quant`)), mapping = aes(x=x, y=y, fill=z)) +
  scale_fill_viridis_c(option = "turbo") + theme_bw() + facet_wrap(facets = ~ factor(id, levels = c("Mean", "Stdev.", "Quantile 0.025", "Quantile 0.975")), nrow = 1) + 
  labs(fill = expression(u(s))) + 
  theme(strip.text = element_text(size=14, face = "bold"), 
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()
  )

ggplot() + 
  geom_ribbon(data = data.frame(model_beta_downscaling.ar1$summary.random$temp_idx, id = "Dowscaling temporal effect (ar1)") %>% rename(., all_of(c(q1='X0.025quant', q3='X0.975quant'))), mapping = aes(x = ID, ymin = q1, ymax = q3), alpha=0.4) +
  geom_line(data = data.frame(model_beta_downscaling.ar1$summary.random$temp_idx, id = "Dowscaling temporal effect (ar1)"), mapping = aes(x=ID, y=mean)) +
  geom_ribbon(data = data.frame(model_beta_downscaling.rw1$summary.random$temp_idx, id = "Dowscaling temporal effect (rw1)") %>% rename(., all_of(c(q1='X0.025quant', q3='X0.975quant'))), mapping = aes(x = ID, ymin = q1, ymax = q3), fill = "red", alpha=0.4) +
  geom_line(data = data.frame(model_beta_downscaling.rw1$summary.random$temp_idx, id = "Dowscaling temporal effect (rw1)"), mapping = aes(x=ID, y=mean), color = "red") +
  geom_ribbon(data = data.frame(model_beta_downscaling.ar1$summary.random$temp_idx, id = "Dowscaling temporal effect (rw2)") %>% rename(., all_of(c(q1='X0.025quant', q3='X0.975quant'))), mapping = aes(x = ID, ymin = q1, ymax = q3), fill = "blue", alpha=0.4) +
  geom_line(data = data.frame(model_beta_downscaling.ar1$summary.random$temp_idx, id = "Dowscaling temporal effect (rw2)"), mapping = aes(x=ID, y=mean), color = "blue") +
  theme_bw() + facet_wrap(facets = ~ id, scales = "free", ncol = 1) +  xlab(label = "Year") + ylab(label = "Values") +
  theme(strip.text = element_text(size=14, face = "bold"))

ggplot() + 
  geom_line(data = bind_rows(lapply(X = model_beta_downscaling.ar1$marginals.fixed, FUN = function(x){as.data.frame(x)}), .id = "id"), mapping = aes(x = x, y = y)) + 
  theme_bw() + facet_wrap(facets = ~ id, ncol = 5, scales = "free") + theme(strip.text = element_text(size=14, face = "bold"))

ggplot() + 
  geom_line(data = bind_rows(lapply(X = model_beta_downscaling.rw1$marginals.fixed, FUN = function(x){as.data.frame(x)}), .id = "id"), mapping = aes(x = x, y = y)) + 
  theme_bw() + facet_wrap(facets = ~ id, ncol = 4, scales = "free")

ggplot() + 
  geom_line(data = bind_rows(lapply(X = model_beta_downscaling.rw2$marginals.fixed, FUN = function(x){as.data.frame(x)}), .id = "id"), mapping = aes(x = x, y = y)) + 
  theme_bw() + facet_wrap(facets = ~ id, ncol = 4, scales = "free")

ggplot() + 
  geom_line(data = bind_rows(lapply(X = model_beta_downscaling.ar1$marginals.hyperpar, FUN = function(x){as.data.frame(x)}), .id = "id"), mapping = aes(x = x, y = y)) + 
  theme_bw() + facet_wrap(facets = ~ id, ncol = 2, scales = "free") + theme(strip.text = element_text(size=14, face = "bold"))

ggplot() + 
  geom_line(data = bind_rows(lapply(X = model_beta_downscaling.rw1$marginals.hyperpar, FUN = function(x){as.data.frame(x)}), .id = "id"), mapping = aes(x = x, y = y)) + 
  theme_bw() + facet_wrap(facets = ~ id, ncol = 2, scales = "free") + theme(strip.text = element_text(size=14, face = "bold"))



