#### Pacotes ####  

library(readxl)
library(sf)
library(terra)
library(landscapemetrics)
library(tidyverse)
library(writexl)
library(vegan)
library(DHARMa)
library(car)
library(pROC)

#### Reclassificando raster ####

ma_2023 <- rast("D:/CETESB/Rasters MapBiomas/ma_2023.tif")
ma_2022 <- rast("D:/CETESB/Rasters MapBiomas/ma_2022.tif")
ma_2021 <- rast("D:/CETESB/Rasters MapBiomas/ma_2021.tif")
ma_2020 <- rast("D:/CETESB/Rasters MapBiomas/ma_2020.tif")
ma_2019 <- rast("D:/CETESB/Rasters MapBiomas/ma_2019.tif")
ma_2018 <- rast("D:/CETESB/Rasters MapBiomas/ma_2018.tif")
ma_2017 <- rast("D:/CETESB/Rasters MapBiomas/ma_2017.tif")
ma_2016 <- rast("D:/CETESB/Rasters MapBiomas/ma_2016.tif")
ma_2015 <- rast("D:/CETESB/Rasters MapBiomas/ma_2015.tif")
ma_2014 <- rast("D:/CETESB/Rasters MapBiomas/ma_2014.tif")
ma_2013 <- rast("D:/CETESB/Rasters MapBiomas/ma_2013.tif")
ma_2012 <- rast("D:/CETESB/Rasters MapBiomas/ma_2012.tif")
ma_2011 <- rast("D:/CETESB/Rasters MapBiomas/ma_2011.tif")
ma_2010 <- rast("D:/CETESB/Rasters MapBiomas/ma_2010.tif")
ma_2009 <- rast("D:/CETESB/Rasters MapBiomas/ma_2009.tif")
ma_2008 <- rast("D:/CETESB/Rasters MapBiomas/ma_2008.tif")
ma_2007 <- rast("D:/CETESB/Rasters MapBiomas/ma_2007.tif")
ma_2006 <- rast("D:/CETESB/Rasters MapBiomas/ma_2006.tif")
ma_2005 <- rast("D:/CETESB/Rasters MapBiomas/br_2005.tif")
ma_2004 <- rast("D:/CETESB/Rasters MapBiomas/br_2004.tif")
ma_2003 <- rast("D:/CETESB/Rasters MapBiomas/ma_2003.tif")
ma_2002 <- rast("D:/CETESB/Rasters MapBiomas/ma_2002.tif")
ma_2001 <- rast("D:/CETESB/Rasters MapBiomas/ma_2001.tif")
ma_2000 <- rast("D:/CETESB/Rasters MapBiomas/ma_2000.tif")
ma_1999 <- rast("D:/CETESB/Rasters MapBiomas/ma_1999.tif")


cetriz_reclas <- c(3, 3, # Formacao florestal
                   6, 3, # Floresta alagavel
                   4, 3, # Savana
                   12, 12, # Formacao campestre
                   15, 15, # Pastagem
                   39, 15, # Soja
                   20, 15, # Cana
                   40, 15, # Arroz
                   62, 15, # Algodao
                   41, 15, # Outras lavouras temporarias
                   46, 15, # Cafe
                   47, 15, # Citrus
                   48, 15, # Outras lavouras perenes
                   21, 15, # Mosaico de usos
                   9, 9, # Silvicultura
                   24, 24, # Area urbanizada
                   0, 0,
                   5, 0,
                   11, 0,
                   23, 0,
                   25, 0, 
                   29, 0,
                   30, 0,
                   31, 0,
                   32, 0,
                   33, 0,
                   49, 0,
                   50, 0)


# Converta para cetriz
cetriz_reclas <- matrix(cetriz_reclas, ncol=2, byrow=TRUE)

# Aplicando a reclassificacao

ma_2023 <- classify(ma_2023, cetriz_reclas)
writeRaster(ma_2023, filename="rec_ma_2023.tif", overwrite=TRUE)

ma_2022 <- classify(ma_2022, cetriz_reclas)
writeRaster(ma_2022, filename="rec_ma_2022.tif", overwrite=TRUE)

ma_2021 <- classify(ma_2021, cetriz_reclas)
writeRaster(ma_2021, filename="rec_ma_2021.tif", overwrite=TRUE)

ma_2020 <- classify(ma_2020, cetriz_reclas)
writeRaster(ma_2020, filename="rec_ma_2020.tif", overwrite=TRUE)

ma_2019 <- classify(ma_2019, cetriz_reclas)
writeRaster(ma_2019, filename="rec_ma_2019.tif", overwrite=TRUE)

ma_2018 <- classify(ma_2018, cetriz_reclas)
writeRaster(ma_2018, filename="rec_ma_2018.tif", overwrite=TRUE)

ma_2017 <- classify(ma_2017, cetriz_reclas)
writeRaster(ma_2017, filename="rec_ma_2017.tif", overwrite=TRUE)

ma_2016 <- classify(ma_2016, cetriz_reclas)
writeRaster(ma_2016, filename="rec_ma_2016.tif", overwrite=TRUE)

ma_2015 <- classify(ma_2015, cetriz_reclas)
writeRaster(ma_2015, filename="rec_ma_2015.tif", overwrite=TRUE)

ma_2014 <- classify(ma_2014, cetriz_reclas)
writeRaster(ma_2014, filename="rec_ma_2014.tif", overwrite=TRUE)

ma_2013 <- classify(ma_2013, cetriz_reclas)
writeRaster(ma_2013, filename="rec_ma_2013.tif", overwrite=TRUE)

ma_2012 <- classify(ma_2012, cetriz_reclas)
writeRaster(ma_2012, filename="rec_ma_2012.tif", overwrite=TRUE)

ma_2011 <- classify(ma_2011, cetriz_reclas)
writeRaster(ma_2011, filename="rec_ma_2011.tif", overwrite=TRUE)

ma_2010 <- classify(ma_2010, cetriz_reclas)
writeRaster(ma_2010, filename="rec_ma_2010.tif", overwrite=TRUE)

ma_2009 <- classify(ma_2009, cetriz_reclas)
writeRaster(ma_2009, filename="rec_ma_2009.tif", overwrite=TRUE)

ma_2008 <- classify(ma_2008, cetriz_reclas)
writeRaster(ma_2008, filename="rec_ma_2008.tif", overwrite=TRUE)

ma_2007 <- classify(ma_2007, cetriz_reclas)
writeRaster(ma_2007, filename="rec_ma_2007.tif", overwrite=TRUE)

ma_2006 <- classify(ma_2006, cetriz_reclas)
writeRaster(ma_2006, filename="rec_ma_2006.tif", overwrite=TRUE)

ma_2005 <- classify(ma_2005, cetriz_reclas)
writeRaster(ma_2005, filename="rec_br_2005.tif", overwrite=TRUE)

ma_2004 <- classify(ma_2004, cetriz_reclas)
writeRaster(ma_2004, filename="rec_br_2004.tif", overwrite=TRUE)

ma_2003 <- classify(ma_2003, cetriz_reclas)
writeRaster(ma_2003, filename="rec_ma_2003.tif", overwrite=TRUE)

ma_2002 <- classify(ma_2002, cetriz_reclas)
writeRaster(ma_2002, filename="rec_ma_2002.tif", overwrite=TRUE)

ma_2001 <- classify(ma_2001, cetriz_reclas)
writeRaster(ma_2001, filename="rec_ma_2001.tif", overwrite=TRUE)

ma_2000 <- classify(ma_2000, cetriz_reclas)
writeRaster(ma_2000, filename="rec_ma_2000.tif", overwrite=TRUE)

ma_1999 <- classify(ma_1999, cetriz_reclas)
writeRaster(ma_1999, filename="rec_ma_1999.tif", overwrite=TRUE)

# Chamando arquivos reclassificados 

ma_2023 <- rast("rec_ma_2023.tif")
ma_2022 <- rast("rec_ma_2022.tif")
ma_2021 <- rast("rec_ma_2021.tif")
ma_2020 <- rast("rec_ma_2020.tif")
ma_2019 <- rast("rec_ma_2019.tif")
ma_2018 <- rast("rec_ma_2018.tif")
ma_2017 <- rast("rec_ma_2017.tif")
ma_2016 <- rast("rec_ma_2016.tif")
ma_2015 <- rast("rec_ma_2015.tif")
ma_2014 <- rast("rec_ma_2014.tif")
ma_2013 <- rast("rec_ma_2013.tif")
ma_2012 <- rast("rec_ma_2012.tif")
ma_2011 <- rast("rec_ma_2011.tif")
ma_2010 <- rast("rec_ma_2010.tif")
ma_2009 <- rast("rec_ma_2009.tif")
ma_2008 <- rast("rec_ma_2008.tif")
ma_2007 <- rast("rec_ma_2007.tif")
ma_2006 <- rast("rec_ma_2006.tif")
ma_2005 <- rast("rec_br_2005.tif")
ma_2004 <- rast("rec_br_2004.tif")
ma_2003 <- rast("rec_ma_2003.tif")
ma_2003 <- rast("rec_ma_2003.tif")
ma_2002 <- rast("rec_ma_2002.tif")
ma_2001 <- rast("rec_ma_2001.tif")
ma_2000 <- rast("rec_ma_2000.tif")

#### Reprojecao  ####

proj_ma_2023 <- project(ma_2023, "EPSG:5880", method = "near")
writeRaster(proj_ma_2023, filename="proj_ma_2023.tif", overwrite=TRUE)

proj_ma_2022 <- project(ma_2022, "EPSG:5880", method = "near")
writeRaster(proj_ma_2022, filename="proj_ma_2022.tif", overwrite=TRUE)

proj_ma_2021 <- project(ma_2021, "EPSG:5880", method = "near")
writeRaster(proj_ma_2021, filename="proj_ma_2021.tif", overwrite=TRUE)

proj_ma_2020 <- project(ma_2020, "EPSG:5880", method = "near")
writeRaster(proj_ma_2020, filename="proj_ma_2020.tif", overwrite=TRUE)

proj_ma_2019 <- project(ma_2019, "EPSG:5880", method = "near")
writeRaster(proj_ma_2019, filename="proj_ma_2019.tif", overwrite=TRUE)

proj_ma_2018 <- project(ma_2018, "EPSG:5880", method = "near")
writeRaster(proj_ma_2018, filename="proj_ma_2018.tif", overwrite=TRUE)

proj_ma_2017 <- project(ma_2017, "EPSG:5880", method = "near")
writeRaster(proj_ma_2017, filename="proj_ma_2017.tif", overwrite=TRUE)

proj_ma_2016 <- project(ma_2016, "EPSG:5880", method = "near")
writeRaster(proj_ma_2016, filename="proj_ma_2016.tif", overwrite=TRUE)

proj_ma_2015 <- project(ma_2015, "EPSG:5880", method = "near")
writeRaster(proj_ma_2015, filename="proj_ma_2015.tif", overwrite=TRUE)

proj_ma_2014 <- project(ma_2014, "EPSG:5880", method = "near")
writeRaster(proj_ma_2014, filename="proj_ma_2014.tif", overwrite=TRUE)

proj_ma_2013 <- project(ma_2013, "EPSG:5880", method = "near")
writeRaster(proj_ma_2013, filename="proj_ma_2013.tif", overwrite=TRUE)

proj_ma_2012 <- project(ma_2012, "EPSG:5880", method = "near")
writeRaster(proj_ma_2012, filename="proj_ma_2012.tif", overwrite=TRUE)

proj_ma_2011 <- project(ma_2011, "EPSG:5880", method = "near")
writeRaster(proj_ma_2011, filename="proj_ma_2011.tif", overwrite=TRUE)

proj_ma_2010 <- project(ma_2010, "EPSG:5880", method = "near") # Continuar daqui
writeRaster(proj_ma_2010, filename="proj_ma_2010.tif", overwrite=TRUE)

proj_ma_2009 <- project(ma_2009, "EPSG:5880", method = "near")
writeRaster(proj_ma_2009, filename="proj_ma_2009.tif", overwrite=TRUE)

proj_ma_2008 <- project(ma_2008, "EPSG:5880", method = "near")
writeRaster(proj_ma_2008, filename="proj_ma_2008.tif", overwrite=TRUE)

proj_ma_2007 <- project(ma_2007, "EPSG:5880", method = "near")
writeRaster(proj_ma_2007, filename="proj_ma_2007.tif", overwrite=TRUE)

proj_ma_2006 <- project(ma_2006, "EPSG:5880", method = "near")
writeRaster(proj_ma_2006, filename="proj_ma_2006.tif", overwrite=TRUE)

proj_ma_2005 <- project(ma_2005, "EPSG:5880", method = "near")
writeRaster(proj_ma_2005, filename="proj_br_2005.tif", overwrite=TRUE)

proj_ma_2004 <- project(ma_2004, "EPSG:5880", method = "near")
writeRaster(proj_ma_2004, filename="proj_br_2004.tif", overwrite=TRUE)

proj_ma_2003 <- project(ma_2003, "EPSG:5880", method = "near")
writeRaster(proj_ma_2003, filename="proj_ma_2003.tif", overwrite=TRUE)

proj_ma_2002 <- project(ma_2002, "EPSG:5880", method = "near")
writeRaster(proj_ma_2002, filename="proj_ma_2002.tif", overwrite=TRUE)

proj_ma_2001 <- project(ma_2001, "EPSG:5880", method = "near")
writeRaster(proj_ma_2001, filename="proj_ma_2001.tif", overwrite=TRUE)

proj_ma_2000 <- project(ma_2000, "EPSG:5880", method = "near")
writeRaster(proj_ma_2000, filename="proj_ma_2000.tif", overwrite=TRUE)

proj_ma_1999 <- project(ma_1999, "EPSG:5880", method = "near")
writeRaster(proj_ma_1999, filename="proj_ma_1999.tif", overwrite=TRUE)

#### Carregando dados ####

tt_ma <- read_excel("D:/CETESB/Rodovias e pontos de ausencia/bin_tt_ma_1999_2023.xlsx", sheet = 1,na = c(" ", "NA")) |> 
  st_as_sf(coords = c("coord_x", "coord_y"), crs=4674)|>
  st_transform(5880)
agua_ma <- st_read("D:/Shapefiles/agua_ma_fei.shp")|>
  st_make_valid()|>
  st_transform(5880)
rod <- st_read("D:/CETESB/Rodovias e pontos de ausencia/rod_final.shp")|>
  st_transform(5880)
UC <- st_read("D:/Shapefiles/UCs 2025/cnuc_2025_03.shp")|>
  st_make_valid()|>
  st_transform(5880)
ma_2023 <- rast("proj_ma_2023.tif")
ma_2022 <- rast("proj_ma_2022.tif")
ma_2021 <- rast("proj_ma_2021.tif")
ma_2020 <- rast("proj_ma_2020.tif")
ma_2019 <- rast("proj_ma_2019.tif")
ma_2018 <- rast("proj_ma_2018.tif")
ma_2017 <- rast("proj_ma_2017.tif")
ma_2016 <- rast("proj_ma_2016.tif")
ma_2015 <- rast("proj_ma_2015.tif")
ma_2014 <- rast("proj_ma_2014.tif")
ma_2013 <- rast("proj_ma_2013.tif")
ma_2012 <- rast("proj_ma_2012.tif")
ma_2011 <- rast("proj_ma_2011.tif")
ma_2010 <- rast("proj_ma_2010.tif")
ma_2009 <- rast("proj_ma_2009.tif")
ma_2008 <- rast("proj_ma_2008.tif")
ma_2007 <- rast("proj_ma_2007.tif")
ma_2006 <- rast("proj_ma_2006.tif")
ma_2005 <- rast("proj_br_2005.tif")
ma_2004 <- rast("proj_br_2004.tif")

#### Verificando duplicatas de pseudo-ausencia ####

dup_tt_ma <- tt_ma |>
  group_by(geometry) |>
  filter(n() > 1) 
dup_tt_ma


tt_ma <- tt_ma|>
  select(-Mes)|>
  dplyr::rename(Bin = id)

#### Inserindo ID nas ocorrencias ####

tt_ma$id_unico <- ave(
  tt_ma$Bin,             
  tt_ma$Bin,             
  FUN = function(x) paste0(x, ".", seq_along(x)))
tt_ma

#### Criando a funcao do indice de proximidade ####

prox <- function(raster_list, target_class = 3, directions = 8, progress = TRUE) {
  
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required")
  }
  
  # Initialize results
  results_df <- data.frame(
    id_unico = names(raster_list),
    prox_median = NA_real_,
    stringsAsFactors = FALSE
  )
  
  if (progress) {
    message("Processing ", length(raster_list), " rasters for class ", target_class)
    pb <- txtProgressBar(min = 0, max = length(raster_list), style = 3)
  }
  
  for (i in seq_along(raster_list)) {
    raster_id <- names(raster_list)[i]
    raster_obj <- raster_list[[i]]
    
    tryCatch({
      # Calculate PROX for current raster using cell centers method
      result <- calculate_prox_single_mall_manters(raster_obj, target_class, directions)
      
      results_df$prox_median[i] <- result$prox_median
      
    }, error = function(e) {
      warning("Error processing raster ", raster_id, ": ", e$message)
      results_df$prox_median[i] <- NA
    })
    
    if (progress) {
      setTxtProgressBar(pb, i)
    }
  }
  
  if (progress) {
    close(pb)
  }
  
  return(results_df)
}

calculate_prox_single_mall_manters <- function(landscape, target_class = 3, directions = 8) {
  
  # Get raster resolution with 4 decimal places
  res_val <- round(terra::res(landscape)[1], 4)
  
  # Create binary raster for target class
  class_raster <- landscape == target_class
  class_raster[class_raster == 0] <- NA
  
  # Check if class exists
  class_malls <- terra::global(class_raster, "sum", na.rm = TRUE)[[1]]
  if (is.na(class_malls) || class_malls == 0) {
    return(list(prox_median = 0, n_patches = 0))
  }
  
  # Get patches
  patches <- terra::patches(class_raster, directions = directions, zeroAsNA = TRUE)
  patch_ids <- terra::unique(patches)[[1]]
  patch_ids <- patch_ids[!is.na(patch_ids)]
  n_patches <- length(patch_ids)
  
  # If only one patch, no neighbors
  if (n_patches <= 1) {
    return(list(prox_median = 0, n_patches = n_patches))
  }
  
  # Pre-calculate patch areas and cell coordinates
  patch_data <- list()
  for (patch_id in patch_ids) {
    patch_mask <- patches == patch_id
    area_malls <- terra::global(patch_mask, "sum", na.rm = TRUE)[[1]]
    
    # Get coordinates of all cells in this patch
    cell_coords <- terra::xyFromCell(patch_mask, which(terra::values(patch_mask) == 1))
    
    patch_data[[as.character(patch_id)]] <- list(
      area_m2 = round(area_malls * (res_val^2), 4),
      cell_coords = cell_coords
    )
  }
  
  # Store PROX values per PATCH
  patch_prox_values <- numeric(0)
  
  # Calculate PROX for each focal patch
  for (i in seq_along(patch_ids)) {
    patch_id_i <- patch_ids[i]
    focal_coords <- patch_data[[as.character(patch_id_i)]]$cell_coords
    
    prox_value <- 0
    
    # Sum PROX contributions from ALL other patches
    for (j in seq_along(patch_ids)) {
      if (i != j) {
        patch_id_j <- patch_ids[j]
        neighbor_coords <- patch_data[[as.character(patch_id_j)]]$cell_coords
        neighbor_area <- patch_data[[as.character(patch_id_j)]]$area_m2
        
        # Calculate minimum distance between cell centers (edge-to-edge)
        min_distance <- calculate_min_distanca_cell_meters(focal_coords, neighbor_coords)
        min_distance <- round(min_distance, 4)
        
        # Add to PROX sum (area of neighbor / distance^2)
        if (min_distance > 0 && is.finite(min_distance)) {
          prox_contrib <- neighbor_area / (min_distance^2)
          prox_value <- prox_value + prox_contrib
        }
      }
    }
    
    # Store the PROX value for this patch
    patch_prox_values <- c(patch_prox_values, prox_value)
  }
  
  # Calculate statistics of PATCH PROX values
  prox_median <- if (length(patch_prox_values) > 0) median(patch_prox_values) else 0
  
  return(list(
    prox_median = prox_median, 
    n_patches = n_patches
  ))
}

# Função otimizada para calcular distância mínima entre células
calculate_min_distanca_cell_meters <- function(coords1, coords2) {
  # Para melhor performance com patches grandes, amostramos se necessário
  if (nrow(coords1) > 1000) {
    coords1 <- coords1[sample(nrow(coords1), 1000), , drop = FALSE]
  }
  if (nrow(coords2) > 1000) {
    coords2 <- coords2[sample(nrow(coords2), 1000), , drop = FALSE]
  }
  
  # Calcula a distância mínima entre quaisquer duas células
  min_dist <- Inf
  for (i in 1:nrow(coords1)) {
    dists <- sqrt((coords1[i, 1] - coords2[, 1])^2 + (coords1[i, 2] - coords2[, 2])^2)
    current_min <- min(dists)
    if (current_min < min_dist) {
      min_dist <- current_min
    }
  }
  return(min_dist)
}

#### Criando a funcao do indice de forma ####

shape <- function(raster_list, class_value = NULL) {
  
  # Verificar se é uma lista
  if(!is.list(raster_list)) {
    stop("O input deve ser uma lista de rasters")
  }
  
  # Função interna para processar um único raster
  process_single_raster <- function(raster_obj, id_unico = NULL) {
    
    # Verificar se o raster é válido
    if(!inherits(raster_obj, "SpatRaster")) {
      warning("Item não é um SpatRaster. Pulando...")
      return(NULL)
    }
    
    # Verificar se o raster é categórico
    if(!is.factor(raster_obj)) {
      warning("O raster '", id_unico, "' não é categórico. Convertendo para fatores...")
      raster_obj <- as.factor(raster_obj)
    }
    
    # Se class_value especificado, verificar se a classe existe no raster
    if(!is.null(class_value)) {
      # Verificar se a classe existe no raster
      unique_vals <- unique(values(raster_obj))
      unique_vals <- unique_vals[!is.na(unique_vals)]
      
      if(!class_value %in% unique_vals) {
        warning("A classe ", class_value, " não existe no raster '", id_unico, "'. Retornando NA.")
        return(tibble(
          id_unico = ifelse(!is.null(id_unico), id_unico, "unknown"),
          median_shape = NA_real_
        ))
      }
      
      # Criar máscara para essa classe
      mask_raster <- raster_obj == class_value
      raster_obj <- mask(raster_obj, mask_raster, maskvalues = 0)
    }
    
    # Verificar se ainda há dados após o masking
    if(all(is.na(values(raster_obj)))) {
      warning("Nenhum dado válido após masking no raster '", id_unico, "'. Retornando NA.")
      return(tibble(
        id_unico = ifelse(!is.null(id_unico), id_unico, "unknown"),
        median_shape = NA_real_
      ))
    }
    
    tryCatch({
      # Converter raster em polígonos (patches)
      patches <- as.polygons(raster_obj) 
      
      # Converter para sf para calcular perímetro
      patches_sf <- st_as_sf(patches) |>
        st_cast("POLYGON")
      
      # Vetor para armazenar os valores de SHAPE
      shape_values <- c()
      
      # Calcular SHAPE para cada patch
      for(i in 1:nrow(patches_sf)) {
        patch <- patches_sf[i, ]
        
        # Calcular área em metros quadrados
        area_m2 <- as.numeric(st_area(patch))
        
        # Calcular perímetro em metros
        perimeter_m <- as.numeric(st_length(st_boundary(patch)))
        
        # Calcular SHAPE index (apenas se área > 0)
        if(area_m2 > 0) {
          shape_index = 0.25 * perimeter_m / sqrt(area_m2)
          shape_values <- c(shape_values, shape_index)
        }
      }
      
      # Calcular estatísticas
      if(length(shape_values) > 0) {
        median_shape <- median(shape_values, na.rm = TRUE)
        
        return(tibble(
          id_unico = ifelse(!is.null(id_unico), id_unico, "unknown"),
          median_shape = median_shape  # Mantém todas as casas decimais
        ))
      } else {
        warning("Nenhum patch válido encontrado no raster '", id_unico, "'.")
        return(tibble(
          id_unico = ifelse(!is.null(id_unico), id_unico, "unknown"),
          median_shape = NA_real_
        ))
      }
    }, error = function(e) {
      warning("Erro ao processar raster '", id_unico, "': ", e$message, ". Retornando NA.")
      return(tibble(
        id_unico = ifelse(!is.null(id_unico), id_unico, "unknown"),
        median_shape = NA_real_
      ))
    })
  }
  
  # Processar cada raster da lista
  results_list <- list()
  
  for(i in 1:length(raster_list)) {
    id_unico <- ifelse(!is.null(names(raster_list)[i]), 
                       names(raster_list)[i], 
                       paste0("raster_", i))
    
    cat("Processando:", id_unico, "\n")
    
    result <- process_single_raster(raster_list[[i]], id_unico)
    
    if(!is.null(result)) {
      results_list[[i]] <- result
    }
  }
  
  # Combinar todos os resultados em um tibble
  if(length(results_list) > 0) {
    final_results <- bind_rows(results_list)
    
    # Estatísticas do processamento
    total_processed <- nrow(final_results)
    successful <- sum(!is.na(final_results$median_shape))
    failed <- total_processed - successful
    
    cat("\n=== RESUMO DO PROCESSAMENTO ===\n")
    cat("Total de rasters processados:", total_processed, "\n")
    cat("Rasters com sucesso:", successful, "\n")
    cat("Rasters com falha:", failed, "\n")
    
    # Configurar opções para mostrar mais casas decimais
    options(pillar.sigfig = 10)
    
    return(final_results)
  } else {
    warning("Nenhum resultado válido obtido.")
    return(tibble(id_unico = character(), median_shape = numeric()))
  }
}

#### Gerando buffer 500m ####

bf_tt_ma_2023_500m <- tt_ma |>
  filter(Ano == "2023") |>
  st_buffer(dist = 500)

bf_tt_ma_2022_500m <- tt_ma |>
  filter(Ano == "2022") |>
  st_buffer(dist = 500)

bf_tt_ma_2021_500m <- tt_ma |>
  filter(Ano == "2021") |>
  st_buffer(dist = 500)

bf_tt_ma_2020_500m <- tt_ma |>
  filter(Ano == "2020") |>
  st_buffer(dist = 500)

bf_tt_ma_2019_500m <- tt_ma |>
  filter(Ano == "2019") |>
  st_buffer(dist = 500)

bf_tt_ma_2018_500m <- tt_ma |>
  filter(Ano == "2018") |>
  st_buffer(dist = 500)

bf_tt_ma_2017_500m <- tt_ma |>
  filter(Ano == "2017") |>
  st_buffer(dist = 500)

bf_tt_ma_2016_500m <- tt_ma |>
  filter(Ano == "2016") |>
  st_buffer(dist = 500)

bf_tt_ma_2015_500m <- tt_ma |>
  filter(Ano == "2015") |>
  st_buffer(dist = 500)

bf_tt_ma_2014_500m <- tt_ma |>
  filter(Ano == "2014") |>
  st_buffer(dist = 500)

bf_tt_ma_2013_500m <- tt_ma |>
  filter(Ano == "2013") |>
  st_buffer(dist = 500)

bf_tt_ma_2012_500m <- tt_ma |>
  filter(Ano == "2012") |>
  st_buffer(dist = 500)

bf_tt_ma_2011_500m <- tt_ma |>
  filter(Ano == "2011") |>
  st_buffer(dist = 500)

bf_tt_ma_2010_500m <- tt_ma |>
  filter(Ano == "2010") |>
  st_buffer(dist = 500)

bf_tt_ma_2009_500m <- tt_ma |>
  filter(Ano == "2009") |>
  st_buffer(dist = 500)

bf_tt_ma_2008_500m <- tt_ma |>
  filter(Ano == "2008") |>
  st_buffer(dist = 500)

bf_tt_ma_2007_500m <- tt_ma |>
  filter(Ano == "2007") |>
  st_buffer(dist = 500)

bf_tt_ma_2006_500m <- tt_ma |>
  filter(Ano == "2006") |>
  st_buffer(dist = 500)

bf_tt_ma_2005_500m <- tt_ma |>
  filter(Ano == "2005") |>
  st_buffer(dist = 500)

bf_tt_ma_2004_500m <- tt_ma |>
  filter(Ano == "2004") |>
  st_buffer(dist = 500)

bf_tt_ma_2003_500m <- tt_ma |>
  filter(Ano == "2003") |>
  st_buffer(dist = 500)

bf_tt_ma_2002_500m <- tt_ma |>
  filter(Ano == "2002") |>
  st_buffer(dist = 500)

bf_tt_ma_2001_500m <- tt_ma |>
  filter(Ano == "2001") |>
  st_buffer(dist = 500)

bf_tt_ma_2000_500m <- tt_ma |>
  filter(Ano == "2000") |>
  st_buffer(dist = 500)

bf_tt_ma_1999_500m <- tt_ma |>
  filter(Ano == "1999") |>
  st_buffer(dist = 500)

#### Cortando raster 500m ####

# 2023

ras_tt_ma_2023_500m <- list()

for (i in 1:nrow(bf_tt_ma_2023_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2023_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2023_500m <- crop(ma_2023, bf_tt_ma_2023_500m[i, ])
  mask_tt_ma_2023_500m <- mask(crop_tt_ma_2023_500m, bf_tt_ma_2023_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2023_500m[[buffer_id]] <- mask_tt_ma_2023_500m
}

output_dir <- "ras_tt_ma_2023_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2023_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2023_500m.tif"))
  writeRaster(
    ras_tt_ma_2023_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2022

ras_tt_ma_2022_500m <- list()

for (i in 1:nrow(bf_tt_ma_2022_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2022_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2022_500m <- crop(ma_2022, bf_tt_ma_2022_500m[i, ])
  mask_tt_ma_2022_500m <- mask(crop_tt_ma_2022_500m, bf_tt_ma_2022_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2022_500m[[buffer_id]] <- mask_tt_ma_2022_500m
}

output_dir <- "ras_tt_ma_2022_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2022_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2022_500m.tif"))
  writeRaster(
    ras_tt_ma_2022_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2021

ras_tt_ma_2021_500m <- list()

for (i in 1:nrow(bf_tt_ma_2021_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2021_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2021_500m <- crop(ma_2021, bf_tt_ma_2021_500m[i, ])
  mask_tt_ma_2021_500m <- mask(crop_tt_ma_2021_500m, bf_tt_ma_2021_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2021_500m[[buffer_id]] <- mask_tt_ma_2021_500m
}

output_dir <- "ras_tt_ma_2021_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2021_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2021_500m.tif"))
  writeRaster(
    ras_tt_ma_2021_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2020

ras_tt_ma_2020_500m <- list()

for (i in 1:nrow(bf_tt_ma_2020_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2020_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2020_500m <- crop(ma_2020, bf_tt_ma_2020_500m[i, ])
  mask_tt_ma_2020_500m <- mask(crop_tt_ma_2020_500m, bf_tt_ma_2020_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2020_500m[[buffer_id]] <- mask_tt_ma_2020_500m
}

output_dir <- "ras_tt_ma_2020_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2020_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2020_500m.tif"))
  writeRaster(
    ras_tt_ma_2020_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2019

ras_tt_ma_2019_500m <- list()

for (i in 1:nrow(bf_tt_ma_2019_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2019_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2019_500m <- crop(ma_2019, bf_tt_ma_2019_500m[i, ])
  mask_tt_ma_2019_500m <- mask(crop_tt_ma_2019_500m, bf_tt_ma_2019_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2019_500m[[buffer_id]] <- mask_tt_ma_2019_500m
}

output_dir <- "ras_tt_ma_2019_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2019_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2019_500m.tif"))
  writeRaster(
    ras_tt_ma_2019_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2018

ras_tt_ma_2018_500m <- list()

for (i in 1:nrow(bf_tt_ma_2018_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2018_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2018_500m <- crop(ma_2018, bf_tt_ma_2018_500m[i, ])
  mask_tt_ma_2018_500m <- mask(crop_tt_ma_2018_500m, bf_tt_ma_2018_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2018_500m[[buffer_id]] <- mask_tt_ma_2018_500m
}

output_dir <- "ras_tt_ma_2018_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2018_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2018_500m.tif"))
  writeRaster(
    ras_tt_ma_2018_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2017

ras_tt_ma_2017_500m <- list()

for (i in 1:nrow(bf_tt_ma_2017_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2017_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2017_500m <- crop(ma_2017, bf_tt_ma_2017_500m[i, ])
  mask_tt_ma_2017_500m <- mask(crop_tt_ma_2017_500m, bf_tt_ma_2017_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2017_500m[[buffer_id]] <- mask_tt_ma_2017_500m
}

output_dir <- "ras_tt_ma_2017_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2017_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2017_500m.tif"))
  writeRaster(
    ras_tt_ma_2017_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2016

ras_tt_ma_2016_500m <- list()

for (i in 1:nrow(bf_tt_ma_2016_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2016_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2016_500m <- crop(ma_2016, bf_tt_ma_2016_500m[i, ])
  mask_tt_ma_2016_500m <- mask(crop_tt_ma_2016_500m, bf_tt_ma_2016_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2016_500m[[buffer_id]] <- mask_tt_ma_2016_500m
}

output_dir <- "ras_tt_ma_2016_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2016_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2016_500m.tif"))
  writeRaster(
    ras_tt_ma_2016_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2015

ras_tt_ma_2015_500m <- list()

for (i in 1:nrow(bf_tt_ma_2015_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2015_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2015_500m <- crop(ma_2015, bf_tt_ma_2015_500m[i, ])
  mask_tt_ma_2015_500m <- mask(crop_tt_ma_2015_500m, bf_tt_ma_2015_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2015_500m[[buffer_id]] <- mask_tt_ma_2015_500m
}

output_dir <- "ras_tt_ma_2015_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2015_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2015_500m.tif"))
  writeRaster(
    ras_tt_ma_2015_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2014

ras_tt_ma_2014_500m <- list()

for (i in 1:nrow(bf_tt_ma_2014_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2014_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2014_500m <- crop(ma_2014, bf_tt_ma_2014_500m[i, ])
  mask_tt_ma_2014_500m <- mask(crop_tt_ma_2014_500m, bf_tt_ma_2014_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2014_500m[[buffer_id]] <- mask_tt_ma_2014_500m
}

output_dir <- "ras_tt_ma_2014_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2014_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2014_500m.tif"))
  writeRaster(
    ras_tt_ma_2014_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2013

ras_tt_ma_2013_500m <- list()

for (i in 1:nrow(bf_tt_ma_2013_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2013_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2013_500m <- crop(ma_2013, bf_tt_ma_2013_500m[i, ])
  mask_tt_ma_2013_500m <- mask(crop_tt_ma_2013_500m, bf_tt_ma_2013_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2013_500m[[buffer_id]] <- mask_tt_ma_2013_500m
}

output_dir <- "ras_tt_ma_2013_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2013_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2013_500m.tif"))
  writeRaster(
    ras_tt_ma_2013_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2012

ras_tt_ma_2012_500m <- list()

for (i in 1:nrow(bf_tt_ma_2012_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2012_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2012_500m <- crop(ma_2012, bf_tt_ma_2012_500m[i, ])
  mask_tt_ma_2012_500m <- mask(crop_tt_ma_2012_500m, bf_tt_ma_2012_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2012_500m[[buffer_id]] <- mask_tt_ma_2012_500m
}

output_dir <- "ras_tt_ma_2012_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2012_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2012_500m.tif"))
  writeRaster(
    ras_tt_ma_2012_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2011

ras_tt_ma_2011_500m <- list()

for (i in 1:nrow(bf_tt_ma_2011_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2011_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2011_500m <- crop(ma_2011, bf_tt_ma_2011_500m[i, ])
  mask_tt_ma_2011_500m <- mask(crop_tt_ma_2011_500m, bf_tt_ma_2011_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2011_500m[[buffer_id]] <- mask_tt_ma_2011_500m
}

output_dir <- "ras_tt_ma_2011_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2011_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2011_500m.tif"))
  writeRaster(
    ras_tt_ma_2011_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2010

ras_tt_ma_2010_500m <- list()

for (i in 1:nrow(bf_tt_ma_2010_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2010_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2010_500m <- crop(ma_2010, bf_tt_ma_2010_500m[i, ])
  mask_tt_ma_2010_500m <- mask(crop_tt_ma_2010_500m, bf_tt_ma_2010_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2010_500m[[buffer_id]] <- mask_tt_ma_2010_500m
}

output_dir <- "ras_tt_ma_2010_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2010_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2010_500m.tif"))
  writeRaster(
    ras_tt_ma_2010_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2009

ras_tt_ma_2009_500m <- list()

for (i in 1:nrow(bf_tt_ma_2009_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2009_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2009_500m <- crop(ma_2009, bf_tt_ma_2009_500m[i, ])
  mask_tt_ma_2009_500m <- mask(crop_tt_ma_2009_500m, bf_tt_ma_2009_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2009_500m[[buffer_id]] <- mask_tt_ma_2009_500m
}

output_dir <- "ras_tt_ma_2009_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2009_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2009_500m.tif"))
  writeRaster(
    ras_tt_ma_2009_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2008

ras_tt_ma_2008_500m <- list()

for (i in 1:nrow(bf_tt_ma_2008_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2008_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2008_500m <- crop(ma_2008, bf_tt_ma_2008_500m[i, ])
  mask_tt_ma_2008_500m <- mask(crop_tt_ma_2008_500m, bf_tt_ma_2008_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2008_500m[[buffer_id]] <- mask_tt_ma_2008_500m
}

output_dir <- "ras_tt_ma_2008_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2008_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2008_500m.tif"))
  writeRaster(
    ras_tt_ma_2008_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2007

ras_tt_ma_2007_500m <- list()

for (i in 1:nrow(bf_tt_ma_2007_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2007_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2007_500m <- crop(ma_2007, bf_tt_ma_2007_500m[i, ])
  mask_tt_ma_2007_500m <- mask(crop_tt_ma_2007_500m, bf_tt_ma_2007_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2007_500m[[buffer_id]] <- mask_tt_ma_2007_500m
}

output_dir <- "ras_tt_ma_2007_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2007_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2007_500m.tif"))
  writeRaster(
    ras_tt_ma_2007_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2006

ras_tt_ma_2006_500m <- list()

for (i in 1:nrow(bf_tt_ma_2006_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2006_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2006_500m <- crop(ma_2006, bf_tt_ma_2006_500m[i, ])
  mask_tt_ma_2006_500m <- mask(crop_tt_ma_2006_500m, bf_tt_ma_2006_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2006_500m[[buffer_id]] <- mask_tt_ma_2006_500m
}

output_dir <- "ras_tt_ma_2006_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2006_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2006_500m.tif"))
  writeRaster(
    ras_tt_ma_2006_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2005

ras_tt_ma_2005_500m <- list()

for (i in 1:nrow(bf_tt_ma_2005_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2005_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2005_500m <- crop(ma_2005, bf_tt_ma_2005_500m[i, ])
  mask_tt_ma_2005_500m <- mask(crop_tt_ma_2005_500m, bf_tt_ma_2005_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2005_500m[[buffer_id]] <- mask_tt_ma_2005_500m
}

output_dir <- "ras_tt_ma_2005_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2005_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2005_500m.tif"))
  writeRaster(
    ras_tt_ma_2005_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2004

ras_tt_ma_2004_500m <- list()

for (i in 1:nrow(bf_tt_ma_2004_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2004_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2004_500m <- crop(ma_2004, bf_tt_ma_2004_500m[i, ])
  mask_tt_ma_2004_500m <- mask(crop_tt_ma_2004_500m, bf_tt_ma_2004_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2004_500m[[buffer_id]] <- mask_tt_ma_2004_500m
}

output_dir <- "ras_tt_ma_2004_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2004_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2004_500m.tif"))
  writeRaster(
    ras_tt_ma_2004_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2003

ras_tt_ma_2003_500m <- list()

for (i in 1:nrow(bf_tt_ma_2003_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2003_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2003_500m <- crop(ma_2003, bf_tt_ma_2003_500m[i, ])
  mask_tt_ma_2003_500m <- mask(crop_tt_ma_2003_500m, bf_tt_ma_2003_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2003_500m[[buffer_id]] <- mask_tt_ma_2003_500m
}

output_dir <- "ras_tt_ma_2003_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2003_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2003_500m.tif"))
  writeRaster(
    ras_tt_ma_2003_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2002

ras_tt_ma_2002_500m <- list()

for (i in 1:nrow(bf_tt_ma_2002_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2002_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2002_500m <- crop(ma_2002, bf_tt_ma_2002_500m[i, ])
  mask_tt_ma_2002_500m <- mask(crop_tt_ma_2002_500m, bf_tt_ma_2002_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2002_500m[[buffer_id]] <- mask_tt_ma_2002_500m
}

output_dir <- "ras_tt_ma_2002_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2002_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2002_500m.tif"))
  writeRaster(
    ras_tt_ma_2002_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2001

ras_tt_ma_2001_500m <- list()

for (i in 1:nrow(bf_tt_ma_2001_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2001_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2001_500m <- crop(ma_2001, bf_tt_ma_2001_500m[i, ])
  mask_tt_ma_2001_500m <- mask(crop_tt_ma_2001_500m, bf_tt_ma_2001_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2001_500m[[buffer_id]] <- mask_tt_ma_2001_500m
}

output_dir <- "ras_tt_ma_2001_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2001_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2001_500m.tif"))
  writeRaster(
    ras_tt_ma_2001_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2000

ras_tt_ma_2000_500m <- list()

for (i in 1:nrow(bf_tt_ma_2000_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2000_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2000_500m <- crop(ma_2000, bf_tt_ma_2000_500m[i, ])
  mask_tt_ma_2000_500m <- mask(crop_tt_ma_2000_500m, bf_tt_ma_2000_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2000_500m[[buffer_id]] <- mask_tt_ma_2000_500m
}

output_dir <- "ras_tt_ma_2000_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2000_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2000_500m.tif"))
  writeRaster(
    ras_tt_ma_2000_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 1999

ras_tt_ma_1999_500m <- list()

for (i in 1:nrow(bf_tt_ma_1999_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_1999_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_1999_500m <- crop(ma_1999, bf_tt_ma_1999_500m[i, ])
  mask_tt_ma_1999_500m <- mask(crop_tt_ma_1999_500m, bf_tt_ma_1999_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_1999_500m[[buffer_id]] <- mask_tt_ma_1999_500m
}

output_dir <- "ras_tt_ma_1999_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_1999_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_1999_500m.tif"))
  writeRaster(
    ras_tt_ma_1999_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# Chamando os recortes

output_dir <- "ras_tt_ma_2023_500m"
ras_tt_ma_2023_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_tt_ma_2023_500m\\.tif$")) |>
  map(rast) 

output_dir <- "ras_tt_ma_2022_500m"
ras_tt_ma_2022_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2022_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2021_500m"
ras_tt_ma_2021_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2021_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2020_500m"
ras_tt_ma_2020_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2020_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2019_500m"
ras_tt_ma_2019_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2019_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2018_500m"
ras_tt_ma_2018_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2018_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2017_500m"
ras_tt_ma_2017_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2017_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2016_500m"
ras_tt_ma_2016_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2016_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2015_500m"
ras_tt_ma_2015_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2015_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2014_500m"
ras_tt_ma_2014_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2014_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2013_500m"
ras_tt_ma_2013_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2013_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2012_500m"
ras_tt_ma_2012_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2012_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2011_500m"
ras_tt_ma_2011_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2011_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2010_500m"
ras_tt_ma_2010_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2010_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2009_500m"
ras_tt_ma_2009_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2009_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2008_500m"
ras_tt_ma_2008_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2008_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2007_500m"
ras_tt_ma_2007_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2007_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2006_500m"
ras_tt_ma_2006_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2006_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2005_500m"
ras_tt_ma_2005_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2005_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2004_500m"
ras_tt_ma_2004_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_tt_ma_2004_500m\\.tif$")) |>
  map(rast) 

output_dir <- "ras_tt_ma_2003_500m"
ras_tt_ma_2003_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_tt_ma_2003_500m\\.tif$")) |>
  map(rast) 

output_dir <- "ras_tt_ma_2002_500m"
ras_tt_ma_2002_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_tt_ma_2002_500m\\.tif$")) |>
  map(rast) 

output_dir <- "ras_tt_ma_2001_500m"
ras_tt_ma_2001_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_tt_ma_2001_500m\\.tif$")) |>
  map(rast) 

output_dir <- "ras_tt_ma_2000_500m"
ras_tt_ma_2000_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_tt_ma_2000_500m\\.tif$")) |>
  map(rast) 

output_dir <- "ras_tt_ma_1999_500m"
ras_tt_ma_1999_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_tt_ma_1999_500m\\.tif$")) |>
  map(rast) 

#### Metricas de paisagem 500m ####

# 2023

id_unico <- names(ras_tt_ma_2023_500m)

met_tt_ma_2023_500m <- map_df(seq_along(ras_tt_ma_2023_500m), function(i) {
  raster <- ras_tt_ma_2023_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2023_500m <- prox(ras_tt_ma_2023_500m, 3)

shape_tt_ma_2023_500m <- shape(ras_tt_ma_2023_500m, class_value = 3)

lsm_tt_ma_2023_500m <- met_tt_ma_2023_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2023")))|>
  inner_join(prox_tt_ma_2023_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2023_500m, by = "id_unico")



# 2022

id_unico <- names(ras_tt_ma_2022_500m)

met_tt_ma_2022_500m <- map_df(seq_along(ras_tt_ma_2022_500m), function(i) {
  raster <- ras_tt_ma_2022_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2022_500m <- prox(ras_tt_ma_2022_500m, 3)

shape_tt_ma_2022_500m <- shape(ras_tt_ma_2022_500m, class_value = 3)

lsm_tt_ma_2022_500m <- met_tt_ma_2022_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2022")))|>
  inner_join(prox_tt_ma_2022_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2022_500m, by = "id_unico")

# 2021

id_unico <- names(ras_tt_ma_2021_500m)

met_tt_ma_2021_500m <- map_df(seq_along(ras_tt_ma_2021_500m), function(i) {
  raster <- ras_tt_ma_2021_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2021_500m <- prox(ras_tt_ma_2021_500m, 3)

shape_tt_ma_2021_500m <- shape(ras_tt_ma_2021_500m, class_value = 3)

lsm_tt_ma_2021_500m <- met_tt_ma_2021_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2021")))|>
  inner_join(prox_tt_ma_2021_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2021_500m, by = "id_unico")

# 2020

id_unico <- names(ras_tt_ma_2020_500m)

met_tt_ma_2020_500m <- map_df(seq_along(ras_tt_ma_2020_500m), function(i) {
  raster <- ras_tt_ma_2020_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2020_500m <- prox(ras_tt_ma_2020_500m, 3)

shape_tt_ma_2020_500m <- shape(ras_tt_ma_2020_500m, class_value = 3)

lsm_tt_ma_2020_500m <- met_tt_ma_2020_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2020")))|>
  inner_join(prox_tt_ma_2020_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2020_500m, by = "id_unico")

# 2019

id_unico <- names(ras_tt_ma_2019_500m)

met_tt_ma_2019_500m <- map_df(seq_along(ras_tt_ma_2019_500m), function(i) {
  raster <- ras_tt_ma_2019_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2019_500m <- prox(ras_tt_ma_2019_500m, 3)

shape_tt_ma_2019_500m <- shape(ras_tt_ma_2019_500m, class_value = 3)

lsm_tt_ma_2019_500m <- met_tt_ma_2019_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2019")))|>
  inner_join(prox_tt_ma_2019_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2019_500m, by = "id_unico")

# 2018

id_unico <- names(ras_tt_ma_2018_500m)

met_tt_ma_2018_500m <- map_df(seq_along(ras_tt_ma_2018_500m), function(i) {
  raster <- ras_tt_ma_2018_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2018_500m <- prox(ras_tt_ma_2018_500m, 3)

shape_tt_ma_2018_500m <- shape(ras_tt_ma_2018_500m, class_value = 3)

lsm_tt_ma_2018_500m <- met_tt_ma_2018_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2018")))|>
  inner_join(prox_tt_ma_2018_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2018_500m, by = "id_unico")

# 2017

id_unico <- names(ras_tt_ma_2017_500m)

met_tt_ma_2017_500m <- map_df(seq_along(ras_tt_ma_2017_500m), function(i) {
  raster <- ras_tt_ma_2017_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2017_500m <- prox(ras_tt_ma_2017_500m, 3)

shape_tt_ma_2017_500m <- shape(ras_tt_ma_2017_500m, class_value = 3)

lsm_tt_ma_2017_500m <- met_tt_ma_2017_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2017")))|>
  inner_join(prox_tt_ma_2017_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2017_500m, by = "id_unico")


# 2016

id_unico <- names(ras_tt_ma_2016_500m)

met_tt_ma_2016_500m <- map_df(seq_along(ras_tt_ma_2016_500m), function(i) {
  raster <- ras_tt_ma_2016_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2016_500m <- prox(ras_tt_ma_2016_500m, 3)

shape_tt_ma_2016_500m <- shape(ras_tt_ma_2016_500m, class_value = 3)

lsm_tt_ma_2016_500m <- met_tt_ma_2016_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2016")))|>
  inner_join(prox_tt_ma_2016_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2016_500m, by = "id_unico")

# 2015

id_unico <- names(ras_tt_ma_2015_500m)

met_tt_ma_2015_500m <- map_df(seq_along(ras_tt_ma_2015_500m), function(i) {
  raster <- ras_tt_ma_2015_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2015_500m <- prox(ras_tt_ma_2015_500m, 3)

shape_tt_ma_2015_500m <- shape(ras_tt_ma_2015_500m, class_value = 3)

lsm_tt_ma_2015_500m <- met_tt_ma_2015_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2015")))|>
  inner_join(prox_tt_ma_2015_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2015_500m, by = "id_unico")

# 2014

id_unico <- names(ras_tt_ma_2014_500m)

met_tt_ma_2014_500m <- map_df(seq_along(ras_tt_ma_2014_500m), function(i) {
  raster <- ras_tt_ma_2014_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2014_500m <- prox(ras_tt_ma_2014_500m, 3)

shape_tt_ma_2014_500m <- shape(ras_tt_ma_2014_500m, class_value = 3)

lsm_tt_ma_2014_500m <- met_tt_ma_2014_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2014")))|>
  inner_join(prox_tt_ma_2014_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2014_500m, by = "id_unico")

# 2013

id_unico <- names(ras_tt_ma_2013_500m)

met_tt_ma_2013_500m <- map_df(seq_along(ras_tt_ma_2013_500m), function(i) {
  raster <- ras_tt_ma_2013_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2013_500m <- prox(ras_tt_ma_2013_500m, 3)

shape_tt_ma_2013_500m <- shape(ras_tt_ma_2013_500m, class_value = 3)

lsm_tt_ma_2013_500m <- met_tt_ma_2013_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2013")))|>
  inner_join(prox_tt_ma_2013_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2013_500m, by = "id_unico")

# 2012
id_unico <- names(ras_tt_ma_2012_500m)

met_tt_ma_2012_500m <- map_df(seq_along(ras_tt_ma_2012_500m), function(i) {
  raster <- ras_tt_ma_2012_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2012_500m <- prox(ras_tt_ma_2012_500m, 3)

shape_tt_ma_2012_500m <- shape(ras_tt_ma_2012_500m, class_value = 3)

lsm_tt_ma_2012_500m <- met_tt_ma_2012_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2012")))|>
  inner_join(prox_tt_ma_2012_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2012_500m, by = "id_unico")

# 2011

id_unico <- names(ras_tt_ma_2011_500m)

met_tt_ma_2011_500m <- map_df(seq_along(ras_tt_ma_2011_500m), function(i) {
  raster <- ras_tt_ma_2011_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2011_500m <- prox(ras_tt_ma_2011_500m, 3)

shape_tt_ma_2011_500m <- shape(ras_tt_ma_2011_500m, class_value = 3)

lsm_tt_ma_2011_500m <- met_tt_ma_2011_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2011")))|>
  inner_join(prox_tt_ma_2011_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2011_500m, by = "id_unico")
# 2010

id_unico <- names(ras_tt_ma_2010_500m)

met_tt_ma_2010_500m <- map_df(seq_along(ras_tt_ma_2010_500m), function(i) {
  raster <- ras_tt_ma_2010_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2010_500m <- prox(ras_tt_ma_2010_500m, 3)

shape_tt_ma_2010_500m <- shape(ras_tt_ma_2010_500m, class_value = 3)

lsm_tt_ma_2010_500m <- met_tt_ma_2010_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2010")))|>
  inner_join(prox_tt_ma_2010_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2010_500m, by = "id_unico")

# 2009
id_unico <- names(ras_tt_ma_2009_500m)

met_tt_ma_2009_500m <- map_df(seq_along(ras_tt_ma_2009_500m), function(i) {
  raster <- ras_tt_ma_2009_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2009_500m <- prox(ras_tt_ma_2009_500m, 3)

shape_tt_ma_2009_500m <- shape(ras_tt_ma_2009_500m, class_value = 3)

lsm_tt_ma_2009_500m <- met_tt_ma_2009_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2009")))|>
  inner_join(prox_tt_ma_2009_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2009_500m, by = "id_unico")

# 2008

id_unico <- names(ras_tt_ma_2008_500m)

met_tt_ma_2008_500m <- map_df(seq_along(ras_tt_ma_2008_500m), function(i) {
  raster <- ras_tt_ma_2008_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2008_500m <- prox(ras_tt_ma_2008_500m, 3)

shape_tt_ma_2008_500m <- shape(ras_tt_ma_2008_500m, class_value = 3)

lsm_tt_ma_2008_500m <- met_tt_ma_2008_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2008")))|>
  inner_join(prox_tt_ma_2008_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2008_500m, by = "id_unico")

# 2007

id_unico <- names(ras_tt_ma_2007_500m)

met_tt_ma_2007_500m <- map_df(seq_along(ras_tt_ma_2007_500m), function(i) {
  raster <- ras_tt_ma_2007_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2007_500m <- prox(ras_tt_ma_2007_500m, 3)

shape_tt_ma_2007_500m <- shape(ras_tt_ma_2007_500m, class_value = 3)

lsm_tt_ma_2007_500m <- met_tt_ma_2007_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2007")))|>
  inner_join(prox_tt_ma_2007_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2007_500m, by = "id_unico")

# 2006

id_unico <- names(ras_tt_ma_2006_500m)

met_tt_ma_2006_500m <- map_df(seq_along(ras_tt_ma_2006_500m), function(i) {
  raster <- ras_tt_ma_2006_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2006_500m <- prox(ras_tt_ma_2006_500m, 3)

shape_tt_ma_2006_500m <- shape(ras_tt_ma_2006_500m, class_value = 3)

lsm_tt_ma_2006_500m <- met_tt_ma_2006_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2006")))|>
  inner_join(prox_tt_ma_2006_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2006_500m, by = "id_unico")

# 2005

id_unico <- names(ras_tt_ma_2005_500m)

met_tt_ma_2005_500m <- map_df(seq_along(ras_tt_ma_2005_500m), function(i) {
  raster <- ras_tt_ma_2005_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2005_500m <- prox(ras_tt_ma_2005_500m, 3)

shape_tt_ma_2005_500m <- shape(ras_tt_ma_2005_500m, class_value = 3)

lsm_tt_ma_2005_500m <- met_tt_ma_2005_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2005")))|>
  inner_join(prox_tt_ma_2005_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2005_500m, by = "id_unico")

# 2004

id_unico <- names(ras_tt_ma_2004_500m)

met_tt_ma_2004_500m <- map_df(seq_along(ras_tt_ma_2004_500m), function(i) {
  raster <- ras_tt_ma_2004_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2004_500m <- prox(ras_tt_ma_2004_500m, 3)

shape_tt_ma_2004_500m <- shape(ras_tt_ma_2004_500m, class_value = 3)

lsm_tt_ma_2004_500m <- met_tt_ma_2004_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2004")))|>
  inner_join(prox_tt_ma_2004_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2004_500m, by = "id_unico")

# 2003

id_unico <- names(ras_tt_ma_2003_500m)

met_tt_ma_2003_500m <- map_df(seq_along(ras_tt_ma_2003_500m), function(i) {
  raster <- ras_tt_ma_2003_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2003_500m <- prox(ras_tt_ma_2003_500m, 3)

shape_tt_ma_2003_500m <- shape(ras_tt_ma_2003_500m, class_value = 3)

lsm_tt_ma_2003_500m <- met_tt_ma_2003_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2003")))|>
  inner_join(prox_tt_ma_2003_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2003_500m, by = "id_unico")

# 2002

id_unico <- names(ras_tt_ma_2002_500m)

met_tt_ma_2002_500m <- map_df(seq_along(ras_tt_ma_2002_500m), function(i) {
  raster <- ras_tt_ma_2002_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2002_500m <- prox(ras_tt_ma_2002_500m, 3)

shape_tt_ma_2002_500m <- shape(ras_tt_ma_2002_500m, class_value = 3)

lsm_tt_ma_2002_500m <- met_tt_ma_2002_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2002")))|>
  inner_join(prox_tt_ma_2002_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2002_500m, by = "id_unico")

# 2001

id_unico <- names(ras_tt_ma_2001_500m)

met_tt_ma_2001_500m <- map_df(seq_along(ras_tt_ma_2001_500m), function(i) {
  raster <- ras_tt_ma_2001_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2001_500m <- prox(ras_tt_ma_2001_500m, 3)

shape_tt_ma_2001_500m <- shape(ras_tt_ma_2001_500m, class_value = 3)

lsm_tt_ma_2001_500m <- met_tt_ma_2001_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2001")))|>
  inner_join(prox_tt_ma_2001_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2001_500m, by = "id_unico")

# 2000

id_unico <- names(ras_tt_ma_2000_500m)

met_tt_ma_2000_500m <- map_df(seq_along(ras_tt_ma_2000_500m), function(i) {
  raster <- ras_tt_ma_2000_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2000_500m <- prox(ras_tt_ma_2000_500m, 3)

shape_tt_ma_2000_500m <- shape(ras_tt_ma_2000_500m, class_value = 3)

lsm_tt_ma_2000_500m <- met_tt_ma_2000_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2000")))|>
  inner_join(prox_tt_ma_2000_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_2000_500m, by = "id_unico")



# 1999

id_unico <- names(ras_tt_ma_1999_500m)

met_tt_ma_1999_500m <- map_df(seq_along(ras_tt_ma_1999_500m), function(i) {
  raster <- ras_tt_ma_1999_500m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_1999_500m <- prox(ras_tt_ma_1999_500m, 3)

shape_tt_ma_1999_500m <- shape(ras_tt_ma_1999_500m, class_value = 3)

lsm_tt_ma_1999_500m <- met_tt_ma_1999_500m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("1999")))|>
  inner_join(prox_tt_ma_1999_500m, by = "id_unico")|>
  inner_join(shape_tt_ma_1999_500m, by = "id_unico")

#### Unificando tabelas 500m ####

lsm_tt_ma_500m <- bind_rows(lsm_tt_ma_2023_500m, lsm_tt_ma_2022_500m, lsm_tt_ma_2021_500m, lsm_tt_ma_2020_500m, lsm_tt_ma_2019_500m, lsm_tt_ma_2018_500m, lsm_tt_ma_2017_500m, lsm_tt_ma_2016_500m, lsm_tt_ma_2015_500m, lsm_tt_ma_2014_500m, lsm_tt_ma_2013_500m, lsm_tt_ma_2012_500m, lsm_tt_ma_2011_500m, lsm_tt_ma_2010_500m, lsm_tt_ma_2009_500m, lsm_tt_ma_2008_500m, lsm_tt_ma_2007_500m, lsm_tt_ma_2006_500m, lsm_tt_ma_2005_500m, lsm_tt_ma_2004_500m, lsm_tt_ma_2003_500m, lsm_tt_ma_2002_500m, lsm_tt_ma_2001_500m, lsm_tt_ma_2000_500m, lsm_tt_ma_1999_500m) |>
  dplyr::select(-pland_0,-lpi_0, -ed_0, -pd_0, -np_0, -lpi_12, -ed_12, -pd_12, -np_12, -lpi_15, -ed_15, -pd_15, -np_15, -lpi_9, -ed_9, -pd_9, -np_9, -np_24, -pd_24, -lpi_24, -ed_24,  -pland_24)|>
  mutate(Bin = str_extract(id_unico, "^[01]"),
         Bin = as.factor(Bin),     
         Ano = as.factor(Ano))|>
  mutate(pland_9 = ifelse(is.na(pland_9), 0, pland_9))|>
  mutate(pland_15 = ifelse(is.na(pland_15), 0, pland_15))|>
  mutate(pland_12 = ifelse(is.na(pland_12), 0, pland_12)) |>
  rename_with(~ paste0(., "_500m"))

lsm_tt_ma_500m_sem_na <- lsm_tt_ma_500m|>
  na.omit() 

View(lsm_tt_ma_500m_sem_na)

# Identifica as linhas com Bin == 0
zeros_500m <- which(lsm_tt_ma_500m_sem_na$Bin_500m == 0)

# Seleciona zeros aleatórios para remover
remove_rows_500m <- sample(zeros_500m, min(613, length(zeros_500m)))

# Filtra o tibble
lsm_tt_ma_500m_total <- lsm_tt_ma_500m_sem_na[-remove_rows_500m, ]

View(lsm_tt_ma_500m_total)

write_xlsx(lsm_tt_ma_500m_total, "lsm_tt_ma_500m.xlsx")

#### Filtando pontos proximos as coberturas florestais ####

tt_ma_flo <- inner_join(lsm_tt_ma_500m_total, tt_ma, by = c("id_unico_500m" = "id_unico"))|>
  select(Bin, id_unico_500m, Ano, geometry)|>
  dplyr::rename('id_unico'='id_unico_500m')|>
  st_as_sf(crs = 5880)

write_xlsx(tt_ma_flo, "tt_ma_flo.xlsx")

st_write(st_as_sf(tt_ma_flo), "tt_ma_flo.shp")

#### Gerando buffer 1000m ####

bf_tt_ma_2023_1000m <-tt_ma_flo |>
  filter(Ano == "2023") |>
  st_buffer(dist = 1000)

bf_tt_ma_2022_1000m <-tt_ma_flo |>
  filter(Ano == "2022") |>
  st_buffer(dist = 1000)

bf_tt_ma_2021_1000m <-tt_ma_flo |>
  filter(Ano == "2021") |>
  st_buffer(dist = 1000)

bf_tt_ma_2020_1000m <-tt_ma_flo |>
  filter(Ano == "2020") |>
  st_buffer(dist = 1000)

bf_tt_ma_2019_1000m <-tt_ma_flo |>
  filter(Ano == "2019") |>
  st_buffer(dist = 1000)

bf_tt_ma_2018_1000m <-tt_ma_flo |>
  filter(Ano == "2018") |>
  st_buffer(dist = 1000)

bf_tt_ma_2017_1000m <-tt_ma_flo |>
  filter(Ano == "2017") |>
  st_buffer(dist = 1000)

bf_tt_ma_2016_1000m <-tt_ma_flo |>
  filter(Ano == "2016") |>
  st_buffer(dist = 1000)

bf_tt_ma_2015_1000m <-tt_ma_flo |>
  filter(Ano == "2015") |>
  st_buffer(dist = 1000)

bf_tt_ma_2014_1000m <-tt_ma_flo |>
  filter(Ano == "2014") |>
  st_buffer(dist = 1000)

bf_tt_ma_2013_1000m <-tt_ma_flo |>
  filter(Ano == "2013") |>
  st_buffer(dist = 1000)

bf_tt_ma_2012_1000m <-tt_ma_flo |>
  filter(Ano == "2012") |>
  st_buffer(dist = 1000)

bf_tt_ma_2011_1000m <-tt_ma_flo |>
  filter(Ano == "2011") |>
  st_buffer(dist = 1000)

bf_tt_ma_2010_1000m <-tt_ma_flo |>
  filter(Ano == "2010") |>
  st_buffer(dist = 1000)

bf_tt_ma_2009_1000m <-tt_ma_flo |>
  filter(Ano == "2009") |>
  st_buffer(dist = 1000)

bf_tt_ma_2008_1000m <-tt_ma_flo |>
  filter(Ano == "2008") |>
  st_buffer(dist = 1000)

bf_tt_ma_2007_1000m <-tt_ma_flo |>
  filter(Ano == "2007") |>
  st_buffer(dist = 1000)

bf_tt_ma_2006_1000m <-tt_ma_flo |>
  filter(Ano == "2006") |>
  st_buffer(dist = 1000)

bf_tt_ma_2005_1000m <-tt_ma_flo |>
  filter(Ano == "2005") |>
  st_buffer(dist = 1000)

bf_tt_ma_2004_1000m <-tt_ma_flo |>
  filter(Ano == "2004") |>
  st_buffer(dist = 1000)

bf_tt_ma_2003_1000m <-tt_ma_flo |>
  filter(Ano == "2003") |>
  st_buffer(dist = 1000)

bf_tt_ma_2002_1000m <-tt_ma_flo |>
  filter(Ano == "2002") |>
  st_buffer(dist = 1000)

bf_tt_ma_2001_1000m <-tt_ma_flo |>
  filter(Ano == "2001") |>
  st_buffer(dist = 1000)

bf_tt_ma_2000_1000m <-tt_ma_flo |>
  filter(Ano == "2000") |>
  st_buffer(dist = 1000)

bf_tt_ma_1999_1000m <-tt_ma_flo |>
  filter(Ano == "1999") |>
  st_buffer(dist = 1000)

#### Cortando raster 1000m ####

# 2023

ras_tt_ma_2023_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2023_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2023_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2023_1000m <- crop(ma_2023, bf_tt_ma_2023_1000m[i, ])
  mask_tt_ma_2023_1000m <- mask(crop_tt_ma_2023_1000m, bf_tt_ma_2023_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2023_1000m[[buffer_id]] <- mask_tt_ma_2023_1000m
}

output_dir <- "ras_tt_ma_2023_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2023_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2023_1000m.tif"))
  writeRaster(
    ras_tt_ma_2023_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2022

ras_tt_ma_2022_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2022_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2022_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2022_1000m <- crop(ma_2022, bf_tt_ma_2022_1000m[i, ])
  mask_tt_ma_2022_1000m <- mask(crop_tt_ma_2022_1000m, bf_tt_ma_2022_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2022_1000m[[buffer_id]] <- mask_tt_ma_2022_1000m
}

output_dir <- "ras_tt_ma_2022_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2022_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2022_1000m.tif"))
  writeRaster(
    ras_tt_ma_2022_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2021

ras_tt_ma_2021_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2021_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2021_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2021_1000m <- crop(ma_2021, bf_tt_ma_2021_1000m[i, ])
  mask_tt_ma_2021_1000m <- mask(crop_tt_ma_2021_1000m, bf_tt_ma_2021_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2021_1000m[[buffer_id]] <- mask_tt_ma_2021_1000m
}

output_dir <- "ras_tt_ma_2021_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2021_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2021_1000m.tif"))
  writeRaster(
    ras_tt_ma_2021_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2020

ras_tt_ma_2020_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2020_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2020_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2020_1000m <- crop(ma_2020, bf_tt_ma_2020_1000m[i, ])
  mask_tt_ma_2020_1000m <- mask(crop_tt_ma_2020_1000m, bf_tt_ma_2020_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2020_1000m[[buffer_id]] <- mask_tt_ma_2020_1000m
}

output_dir <- "ras_tt_ma_2020_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2020_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2020_1000m.tif"))
  writeRaster(
    ras_tt_ma_2020_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2019

ras_tt_ma_2019_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2019_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2019_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2019_1000m <- crop(ma_2019, bf_tt_ma_2019_1000m[i, ])
  mask_tt_ma_2019_1000m <- mask(crop_tt_ma_2019_1000m, bf_tt_ma_2019_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2019_1000m[[buffer_id]] <- mask_tt_ma_2019_1000m
}

output_dir <- "ras_tt_ma_2019_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2019_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2019_1000m.tif"))
  writeRaster(
    ras_tt_ma_2019_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2018

ras_tt_ma_2018_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2018_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2018_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2018_1000m <- crop(ma_2018, bf_tt_ma_2018_1000m[i, ])
  mask_tt_ma_2018_1000m <- mask(crop_tt_ma_2018_1000m, bf_tt_ma_2018_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2018_1000m[[buffer_id]] <- mask_tt_ma_2018_1000m
}

output_dir <- "ras_tt_ma_2018_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2018_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2018_1000m.tif"))
  writeRaster(
    ras_tt_ma_2018_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2017

ras_tt_ma_2017_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2017_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2017_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2017_1000m <- crop(ma_2017, bf_tt_ma_2017_1000m[i, ])
  mask_tt_ma_2017_1000m <- mask(crop_tt_ma_2017_1000m, bf_tt_ma_2017_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2017_1000m[[buffer_id]] <- mask_tt_ma_2017_1000m
}

output_dir <- "ras_tt_ma_2017_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2017_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2017_1000m.tif"))
  writeRaster(
    ras_tt_ma_2017_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2016

ras_tt_ma_2016_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2016_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2016_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2016_1000m <- crop(ma_2016, bf_tt_ma_2016_1000m[i, ])
  mask_tt_ma_2016_1000m <- mask(crop_tt_ma_2016_1000m, bf_tt_ma_2016_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2016_1000m[[buffer_id]] <- mask_tt_ma_2016_1000m
}

output_dir <- "ras_tt_ma_2016_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2016_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2016_1000m.tif"))
  writeRaster(
    ras_tt_ma_2016_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2015

ras_tt_ma_2015_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2015_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2015_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2015_1000m <- crop(ma_2015, bf_tt_ma_2015_1000m[i, ])
  mask_tt_ma_2015_1000m <- mask(crop_tt_ma_2015_1000m, bf_tt_ma_2015_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2015_1000m[[buffer_id]] <- mask_tt_ma_2015_1000m
}

output_dir <- "ras_tt_ma_2015_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2015_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2015_1000m.tif"))
  writeRaster(
    ras_tt_ma_2015_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2014

ras_tt_ma_2014_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2014_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2014_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2014_1000m <- crop(ma_2014, bf_tt_ma_2014_1000m[i, ])
  mask_tt_ma_2014_1000m <- mask(crop_tt_ma_2014_1000m, bf_tt_ma_2014_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2014_1000m[[buffer_id]] <- mask_tt_ma_2014_1000m
}

output_dir <- "ras_tt_ma_2014_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2014_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2014_1000m.tif"))
  writeRaster(
    ras_tt_ma_2014_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2013

ras_tt_ma_2013_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2013_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2013_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2013_1000m <- crop(ma_2013, bf_tt_ma_2013_1000m[i, ])
  mask_tt_ma_2013_1000m <- mask(crop_tt_ma_2013_1000m, bf_tt_ma_2013_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2013_1000m[[buffer_id]] <- mask_tt_ma_2013_1000m
}

output_dir <- "ras_tt_ma_2013_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2013_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2013_1000m.tif"))
  writeRaster(
    ras_tt_ma_2013_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2012

ras_tt_ma_2012_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2012_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2012_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2012_1000m <- crop(ma_2012, bf_tt_ma_2012_1000m[i, ])
  mask_tt_ma_2012_1000m <- mask(crop_tt_ma_2012_1000m, bf_tt_ma_2012_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2012_1000m[[buffer_id]] <- mask_tt_ma_2012_1000m
}

output_dir <- "ras_tt_ma_2012_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2012_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2012_1000m.tif"))
  writeRaster(
    ras_tt_ma_2012_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2011

ras_tt_ma_2011_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2011_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2011_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2011_1000m <- crop(ma_2011, bf_tt_ma_2011_1000m[i, ])
  mask_tt_ma_2011_1000m <- mask(crop_tt_ma_2011_1000m, bf_tt_ma_2011_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2011_1000m[[buffer_id]] <- mask_tt_ma_2011_1000m
}

output_dir <- "ras_tt_ma_2011_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2011_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2011_1000m.tif"))
  writeRaster(
    ras_tt_ma_2011_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2010

ras_tt_ma_2010_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2010_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2010_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2010_1000m <- crop(ma_2010, bf_tt_ma_2010_1000m[i, ])
  mask_tt_ma_2010_1000m <- mask(crop_tt_ma_2010_1000m, bf_tt_ma_2010_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2010_1000m[[buffer_id]] <- mask_tt_ma_2010_1000m
}

output_dir <- "ras_tt_ma_2010_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2010_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2010_1000m.tif"))
  writeRaster(
    ras_tt_ma_2010_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2009

ras_tt_ma_2009_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2009_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2009_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2009_1000m <- crop(ma_2009, bf_tt_ma_2009_1000m[i, ])
  mask_tt_ma_2009_1000m <- mask(crop_tt_ma_2009_1000m, bf_tt_ma_2009_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2009_1000m[[buffer_id]] <- mask_tt_ma_2009_1000m
}

output_dir <- "ras_tt_ma_2009_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2009_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2009_1000m.tif"))
  writeRaster(
    ras_tt_ma_2009_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2008

ras_tt_ma_2008_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2008_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2008_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2008_1000m <- crop(ma_2008, bf_tt_ma_2008_1000m[i, ])
  mask_tt_ma_2008_1000m <- mask(crop_tt_ma_2008_1000m, bf_tt_ma_2008_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2008_1000m[[buffer_id]] <- mask_tt_ma_2008_1000m
}

output_dir <- "ras_tt_ma_2008_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2008_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2008_1000m.tif"))
  writeRaster(
    ras_tt_ma_2008_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2007

ras_tt_ma_2007_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2007_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2007_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2007_1000m <- crop(ma_2007, bf_tt_ma_2007_1000m[i, ])
  mask_tt_ma_2007_1000m <- mask(crop_tt_ma_2007_1000m, bf_tt_ma_2007_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2007_1000m[[buffer_id]] <- mask_tt_ma_2007_1000m
}

output_dir <- "ras_tt_ma_2007_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2007_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2007_1000m.tif"))
  writeRaster(
    ras_tt_ma_2007_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2006

ras_tt_ma_2006_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2006_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2006_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2006_1000m <- crop(ma_2006, bf_tt_ma_2006_1000m[i, ])
  mask_tt_ma_2006_1000m <- mask(crop_tt_ma_2006_1000m, bf_tt_ma_2006_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2006_1000m[[buffer_id]] <- mask_tt_ma_2006_1000m
}

output_dir <- "ras_tt_ma_2006_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2006_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2006_1000m.tif"))
  writeRaster(
    ras_tt_ma_2006_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2005

ras_tt_ma_2005_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2005_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2005_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2005_1000m <- crop(ma_2005, bf_tt_ma_2005_1000m[i, ])
  mask_tt_ma_2005_1000m <- mask(crop_tt_ma_2005_1000m, bf_tt_ma_2005_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2005_1000m[[buffer_id]] <- mask_tt_ma_2005_1000m
}

output_dir <- "ras_tt_ma_2005_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2005_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2005_1000m.tif"))
  writeRaster(
    ras_tt_ma_2005_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2004

ras_tt_ma_2004_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2004_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2004_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2004_1000m <- crop(ma_2004, bf_tt_ma_2004_1000m[i, ])
  mask_tt_ma_2004_1000m <- mask(crop_tt_ma_2004_1000m, bf_tt_ma_2004_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2004_1000m[[buffer_id]] <- mask_tt_ma_2004_1000m
}

output_dir <- "ras_tt_ma_2004_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2004_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2004_1000m.tif"))
  writeRaster(
    ras_tt_ma_2004_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2003

ras_tt_ma_2003_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2003_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2003_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2003_1000m <- crop(ma_2003, bf_tt_ma_2003_1000m[i, ])
  mask_tt_ma_2003_1000m <- mask(crop_tt_ma_2003_1000m, bf_tt_ma_2003_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2003_1000m[[buffer_id]] <- mask_tt_ma_2003_1000m
}

output_dir <- "ras_tt_ma_2003_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2003_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2003_1000m.tif"))
  writeRaster(
    ras_tt_ma_2003_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2002

ras_tt_ma_2002_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2002_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2002_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2002_1000m <- crop(ma_2002, bf_tt_ma_2002_1000m[i, ])
  mask_tt_ma_2002_1000m <- mask(crop_tt_ma_2002_1000m, bf_tt_ma_2002_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2002_1000m[[buffer_id]] <- mask_tt_ma_2002_1000m
}

output_dir <- "ras_tt_ma_2002_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2002_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2002_1000m.tif"))
  writeRaster(
    ras_tt_ma_2002_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2001

ras_tt_ma_2001_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2001_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2001_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2001_1000m <- crop(ma_2001, bf_tt_ma_2001_1000m[i, ])
  mask_tt_ma_2001_1000m <- mask(crop_tt_ma_2001_1000m, bf_tt_ma_2001_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2001_1000m[[buffer_id]] <- mask_tt_ma_2001_1000m
}

output_dir <- "ras_tt_ma_2001_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2001_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2001_1000m.tif"))
  writeRaster(
    ras_tt_ma_2001_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2000

ras_tt_ma_2000_1000m <- list()

for (i in 1:nrow(bf_tt_ma_2000_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2000_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2000_1000m <- crop(ma_2000, bf_tt_ma_2000_1000m[i, ])
  mask_tt_ma_2000_1000m <- mask(crop_tt_ma_2000_1000m, bf_tt_ma_2000_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2000_1000m[[buffer_id]] <- mask_tt_ma_2000_1000m
}

output_dir <- "ras_tt_ma_2000_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2000_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2000_1000m.tif"))
  writeRaster(
    ras_tt_ma_2000_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 1999

ras_tt_ma_1999_1000m <- list()

for (i in 1:nrow(bf_tt_ma_1999_1000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_1999_1000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_1999_1000m <- crop(ma_1999, bf_tt_ma_1999_1000m[i, ])
  mask_tt_ma_1999_1000m <- mask(crop_tt_ma_1999_1000m, bf_tt_ma_1999_1000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_1999_1000m[[buffer_id]] <- mask_tt_ma_1999_1000m
}

output_dir <- "ras_tt_ma_1999_1000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_1999_1000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_1999_1000m.tif"))
  writeRaster(
    ras_tt_ma_1999_1000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# Chamando os recortes

output_dir <- "ras_tt_ma_2023_1000m"
ras_tt_ma_2023_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_tt_ma_2023_1000m\\.tif$")) |>
  map(rast) 

output_dir <- "ras_tt_ma_2022_1000m"
ras_tt_ma_2022_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2022_1000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2021_1000m"
ras_tt_ma_2021_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2021_1000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2020_1000m"
ras_tt_ma_2020_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2020_1000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2019_1000m"
ras_tt_ma_2019_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2019_1000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2018_1000m"
ras_tt_ma_2018_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2018_1000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2017_1000m"
ras_tt_ma_2017_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2017_1000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2016_1000m"
ras_tt_ma_2016_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2016_1000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2015_1000m"
ras_tt_ma_2015_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2015_1000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2014_1000m"
ras_tt_ma_2014_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2014_1000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2013_1000m"
ras_tt_ma_2013_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2013_1000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2012_1000m"
ras_tt_ma_2012_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2012_1000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2011_1000m"
ras_tt_ma_2011_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2011_1000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2010_1000m"
ras_tt_ma_2010_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2010_1000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2009_1000m"
ras_tt_ma_2009_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2009_1000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2008_1000m"
ras_tt_ma_2008_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2008_1000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2007_1000m"
ras_tt_ma_2007_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2007_1000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2006_1000m"
ras_tt_ma_2006_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2006_1000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2005_1000m"
ras_tt_ma_2005_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2005_1000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2004_1000m"
ras_tt_ma_2004_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_tt_ma_2004_1000m\\.tif$")) |>
  map(rast) 

output_dir <- "ras_tt_ma_2003_1000m"
ras_tt_ma_2003_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_tt_ma_2003_1000m\\.tif$")) |>
  map(rast) 

output_dir <- "ras_tt_ma_2002_1000m"
ras_tt_ma_2002_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_tt_ma_2002_1000m\\.tif$")) |>
  map(rast) 

output_dir <- "ras_tt_ma_2001_1000m"
ras_tt_ma_2001_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_tt_ma_2001_1000m\\.tif$")) |>
  map(rast) 

output_dir <- "ras_tt_ma_2000_1000m"
ras_tt_ma_2000_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_tt_ma_2000_1000m\\.tif$")) |>
  map(rast) 

output_dir <- "ras_tt_ma_1999_1000m"
ras_tt_ma_1999_1000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_tt_ma_1999_1000m\\.tif$")) |>
  map(rast) 

#### Metricas de paisagem 1000m ####

# 2023

id_unico <- names(ras_tt_ma_2023_1000m)

met_tt_ma_2023_1000m <- map_df(seq_along(ras_tt_ma_2023_1000m), function(i) {
  raster <- ras_tt_ma_2023_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2023_1000m <- prox(ras_tt_ma_2023_1000m, 3)

shape_tt_ma_2023_1000m <- shape(ras_tt_ma_2023_1000m, class_value = 3)

lsm_tt_ma_2023_1000m <- met_tt_ma_2023_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2023")))|>
  inner_join(prox_tt_ma_2023_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2023_1000m, by = "id_unico")



# 2022

id_unico <- names(ras_tt_ma_2022_1000m)

met_tt_ma_2022_1000m <- map_df(seq_along(ras_tt_ma_2022_1000m), function(i) {
  raster <- ras_tt_ma_2022_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2022_1000m <- prox(ras_tt_ma_2022_1000m, 3)

shape_tt_ma_2022_1000m <- shape(ras_tt_ma_2022_1000m, class_value = 3)

lsm_tt_ma_2022_1000m <- met_tt_ma_2022_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2022")))|>
  inner_join(prox_tt_ma_2022_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2022_1000m, by = "id_unico")

# 2021

id_unico <- names(ras_tt_ma_2021_1000m)

met_tt_ma_2021_1000m <- map_df(seq_along(ras_tt_ma_2021_1000m), function(i) {
  raster <- ras_tt_ma_2021_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2021_1000m <- prox(ras_tt_ma_2021_1000m, 3)

shape_tt_ma_2021_1000m <- shape(ras_tt_ma_2021_1000m, class_value = 3)

lsm_tt_ma_2021_1000m <- met_tt_ma_2021_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2021")))|>
  inner_join(prox_tt_ma_2021_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2021_1000m, by = "id_unico")

# 2020

id_unico <- names(ras_tt_ma_2020_1000m)

met_tt_ma_2020_1000m <- map_df(seq_along(ras_tt_ma_2020_1000m), function(i) {
  raster <- ras_tt_ma_2020_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2020_1000m <- prox(ras_tt_ma_2020_1000m, 3)

shape_tt_ma_2020_1000m <- shape(ras_tt_ma_2020_1000m, class_value = 3)

lsm_tt_ma_2020_1000m <- met_tt_ma_2020_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2020")))|>
  inner_join(prox_tt_ma_2020_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2020_1000m, by = "id_unico")

# 2019

id_unico <- names(ras_tt_ma_2019_1000m)

met_tt_ma_2019_1000m <- map_df(seq_along(ras_tt_ma_2019_1000m), function(i) {
  raster <- ras_tt_ma_2019_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2019_1000m <- prox(ras_tt_ma_2019_1000m, 3)

shape_tt_ma_2019_1000m <- shape(ras_tt_ma_2019_1000m, class_value = 3)

lsm_tt_ma_2019_1000m <- met_tt_ma_2019_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2019")))|>
  inner_join(prox_tt_ma_2019_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2019_1000m, by = "id_unico")

# 2018

id_unico <- names(ras_tt_ma_2018_1000m)

met_tt_ma_2018_1000m <- map_df(seq_along(ras_tt_ma_2018_1000m), function(i) {
  raster <- ras_tt_ma_2018_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2018_1000m <- prox(ras_tt_ma_2018_1000m, 3)

shape_tt_ma_2018_1000m <- shape(ras_tt_ma_2018_1000m, class_value = 3)

lsm_tt_ma_2018_1000m <- met_tt_ma_2018_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2018")))|>
  inner_join(prox_tt_ma_2018_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2018_1000m, by = "id_unico")

# 2017

id_unico <- names(ras_tt_ma_2017_1000m)

met_tt_ma_2017_1000m <- map_df(seq_along(ras_tt_ma_2017_1000m), function(i) {
  raster <- ras_tt_ma_2017_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2017_1000m <- prox(ras_tt_ma_2017_1000m, 3)

shape_tt_ma_2017_1000m <- shape(ras_tt_ma_2017_1000m, class_value = 3)

lsm_tt_ma_2017_1000m <- met_tt_ma_2017_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2017")))|>
  inner_join(prox_tt_ma_2017_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2017_1000m, by = "id_unico")


# 2016

id_unico <- names(ras_tt_ma_2016_1000m)

met_tt_ma_2016_1000m <- map_df(seq_along(ras_tt_ma_2016_1000m), function(i) {
  raster <- ras_tt_ma_2016_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2016_1000m <- prox(ras_tt_ma_2016_1000m, 3)

shape_tt_ma_2016_1000m <- shape(ras_tt_ma_2016_1000m, class_value = 3)

lsm_tt_ma_2016_1000m <- met_tt_ma_2016_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2016")))|>
  inner_join(prox_tt_ma_2016_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2016_1000m, by = "id_unico")

# 2015

id_unico <- names(ras_tt_ma_2015_1000m)

met_tt_ma_2015_1000m <- map_df(seq_along(ras_tt_ma_2015_1000m), function(i) {
  raster <- ras_tt_ma_2015_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2015_1000m <- prox(ras_tt_ma_2015_1000m, 3)

shape_tt_ma_2015_1000m <- shape(ras_tt_ma_2015_1000m, class_value = 3)

lsm_tt_ma_2015_1000m <- met_tt_ma_2015_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2015")))|>
  inner_join(prox_tt_ma_2015_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2015_1000m, by = "id_unico")

# 2014

id_unico <- names(ras_tt_ma_2014_1000m)

met_tt_ma_2014_1000m <- map_df(seq_along(ras_tt_ma_2014_1000m), function(i) {
  raster <- ras_tt_ma_2014_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2014_1000m <- prox(ras_tt_ma_2014_1000m, 3)

shape_tt_ma_2014_1000m <- shape(ras_tt_ma_2014_1000m, class_value = 3)

lsm_tt_ma_2014_1000m <- met_tt_ma_2014_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2014")))|>
  inner_join(prox_tt_ma_2014_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2014_1000m, by = "id_unico")

# 2013

id_unico <- names(ras_tt_ma_2013_1000m)

met_tt_ma_2013_1000m <- map_df(seq_along(ras_tt_ma_2013_1000m), function(i) {
  raster <- ras_tt_ma_2013_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2013_1000m <- prox(ras_tt_ma_2013_1000m, 3)

shape_tt_ma_2013_1000m <- shape(ras_tt_ma_2013_1000m, class_value = 3)

lsm_tt_ma_2013_1000m <- met_tt_ma_2013_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2013")))|>
  inner_join(prox_tt_ma_2013_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2013_1000m, by = "id_unico")

# 2012
id_unico <- names(ras_tt_ma_2012_1000m)

met_tt_ma_2012_1000m <- map_df(seq_along(ras_tt_ma_2012_1000m), function(i) {
  raster <- ras_tt_ma_2012_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2012_1000m <- prox(ras_tt_ma_2012_1000m, 3)

shape_tt_ma_2012_1000m <- shape(ras_tt_ma_2012_1000m, class_value = 3)

lsm_tt_ma_2012_1000m <- met_tt_ma_2012_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2012")))|>
  inner_join(prox_tt_ma_2012_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2012_1000m, by = "id_unico")

# 2011

id_unico <- names(ras_tt_ma_2011_1000m)

met_tt_ma_2011_1000m <- map_df(seq_along(ras_tt_ma_2011_1000m), function(i) {
  raster <- ras_tt_ma_2011_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2011_1000m <- prox(ras_tt_ma_2011_1000m, 3)

shape_tt_ma_2011_1000m <- shape(ras_tt_ma_2011_1000m, class_value = 3)

lsm_tt_ma_2011_1000m <- met_tt_ma_2011_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2011")))|>
  inner_join(prox_tt_ma_2011_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2011_1000m, by = "id_unico")
# 2010

id_unico <- names(ras_tt_ma_2010_1000m)

met_tt_ma_2010_1000m <- map_df(seq_along(ras_tt_ma_2010_1000m), function(i) {
  raster <- ras_tt_ma_2010_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2010_1000m <- prox(ras_tt_ma_2010_1000m, 3)

shape_tt_ma_2010_1000m <- shape(ras_tt_ma_2010_1000m, class_value = 3)

lsm_tt_ma_2010_1000m <- met_tt_ma_2010_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2010")))|>
  inner_join(prox_tt_ma_2010_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2010_1000m, by = "id_unico")

# 2009
id_unico <- names(ras_tt_ma_2009_1000m)

met_tt_ma_2009_1000m <- map_df(seq_along(ras_tt_ma_2009_1000m), function(i) {
  raster <- ras_tt_ma_2009_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2009_1000m <- prox(ras_tt_ma_2009_1000m, 3)

shape_tt_ma_2009_1000m <- shape(ras_tt_ma_2009_1000m, class_value = 3)

lsm_tt_ma_2009_1000m <- met_tt_ma_2009_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2009")))|>
  inner_join(prox_tt_ma_2009_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2009_1000m, by = "id_unico")

# 2008

id_unico <- names(ras_tt_ma_2008_1000m)

met_tt_ma_2008_1000m <- map_df(seq_along(ras_tt_ma_2008_1000m), function(i) {
  raster <- ras_tt_ma_2008_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2008_1000m <- prox(ras_tt_ma_2008_1000m, 3)

shape_tt_ma_2008_1000m <- shape(ras_tt_ma_2008_1000m, class_value = 3)

lsm_tt_ma_2008_1000m <- met_tt_ma_2008_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2008")))|>
  inner_join(prox_tt_ma_2008_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2008_1000m, by = "id_unico")

# 2007

id_unico <- names(ras_tt_ma_2007_1000m)

met_tt_ma_2007_1000m <- map_df(seq_along(ras_tt_ma_2007_1000m), function(i) {
  raster <- ras_tt_ma_2007_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2007_1000m <- prox(ras_tt_ma_2007_1000m, 3)

shape_tt_ma_2007_1000m <- shape(ras_tt_ma_2007_1000m, class_value = 3)

lsm_tt_ma_2007_1000m <- met_tt_ma_2007_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2007")))|>
  inner_join(prox_tt_ma_2007_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2007_1000m, by = "id_unico")

# 2006

id_unico <- names(ras_tt_ma_2006_1000m)

met_tt_ma_2006_1000m <- map_df(seq_along(ras_tt_ma_2006_1000m), function(i) {
  raster <- ras_tt_ma_2006_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2006_1000m <- prox(ras_tt_ma_2006_1000m, 3)

shape_tt_ma_2006_1000m <- shape(ras_tt_ma_2006_1000m, class_value = 3)

lsm_tt_ma_2006_1000m <- met_tt_ma_2006_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2006")))|>
  inner_join(prox_tt_ma_2006_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2006_1000m, by = "id_unico")

# 2005

id_unico <- names(ras_tt_ma_2005_1000m)

met_tt_ma_2005_1000m <- map_df(seq_along(ras_tt_ma_2005_1000m), function(i) {
  raster <- ras_tt_ma_2005_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2005_1000m <- prox(ras_tt_ma_2005_1000m, 3)

shape_tt_ma_2005_1000m <- shape(ras_tt_ma_2005_1000m, class_value = 3)

lsm_tt_ma_2005_1000m <- met_tt_ma_2005_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2005")))|>
  inner_join(prox_tt_ma_2005_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2005_1000m, by = "id_unico")

# 2004

id_unico <- names(ras_tt_ma_2004_1000m)

met_tt_ma_2004_1000m <- map_df(seq_along(ras_tt_ma_2004_1000m), function(i) {
  raster <- ras_tt_ma_2004_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2004_1000m <- prox(ras_tt_ma_2004_1000m, 3)

shape_tt_ma_2004_1000m <- shape(ras_tt_ma_2004_1000m, class_value = 3)

lsm_tt_ma_2004_1000m <- met_tt_ma_2004_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2004")))|>
  inner_join(prox_tt_ma_2004_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2004_1000m, by = "id_unico")

# 2003

id_unico <- names(ras_tt_ma_2003_1000m)

met_tt_ma_2003_1000m <- map_df(seq_along(ras_tt_ma_2003_1000m), function(i) {
  raster <- ras_tt_ma_2003_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2003_1000m <- prox(ras_tt_ma_2003_1000m, 3)

shape_tt_ma_2003_1000m <- shape(ras_tt_ma_2003_1000m, class_value = 3)

lsm_tt_ma_2003_1000m <- met_tt_ma_2003_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2003")))|>
  inner_join(prox_tt_ma_2003_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2003_1000m, by = "id_unico")

# 2002

id_unico <- names(ras_tt_ma_2002_1000m)

met_tt_ma_2002_1000m <- map_df(seq_along(ras_tt_ma_2002_1000m), function(i) {
  raster <- ras_tt_ma_2002_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2002_1000m <- prox(ras_tt_ma_2002_1000m, 3)

shape_tt_ma_2002_1000m <- shape(ras_tt_ma_2002_1000m, class_value = 3)

lsm_tt_ma_2002_1000m <- met_tt_ma_2002_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2002")))|>
  inner_join(prox_tt_ma_2002_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2002_1000m, by = "id_unico")

# 2001

id_unico <- names(ras_tt_ma_2001_1000m)

met_tt_ma_2001_1000m <- map_df(seq_along(ras_tt_ma_2001_1000m), function(i) {
  raster <- ras_tt_ma_2001_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2001_1000m <- prox(ras_tt_ma_2001_1000m, 3)

shape_tt_ma_2001_1000m <- shape(ras_tt_ma_2001_1000m, class_value = 3)

lsm_tt_ma_2001_1000m <- met_tt_ma_2001_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2001")))|>
  inner_join(prox_tt_ma_2001_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2001_1000m, by = "id_unico")

# 2000

id_unico <- names(ras_tt_ma_2000_1000m)

met_tt_ma_2000_1000m <- map_df(seq_along(ras_tt_ma_2000_1000m), function(i) {
  raster <- ras_tt_ma_2000_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2000_1000m <- prox(ras_tt_ma_2000_1000m, 3)

shape_tt_ma_2000_1000m <- shape(ras_tt_ma_2000_1000m, class_value = 3)

lsm_tt_ma_2000_1000m <- met_tt_ma_2000_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2000")))|>
  inner_join(prox_tt_ma_2000_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2000_1000m, by = "id_unico")



# 1999

id_unico <- names(ras_tt_ma_1999_1000m)

met_tt_ma_1999_1000m <- map_df(seq_along(ras_tt_ma_1999_1000m), function(i) {
  raster <- ras_tt_ma_1999_1000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_1999_1000m <- prox(ras_tt_ma_1999_1000m, 3)

shape_tt_ma_1999_1000m <- shape(ras_tt_ma_1999_1000m, class_value = 3)

lsm_tt_ma_1999_1000m <- met_tt_ma_1999_1000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("1999")))|>
  inner_join(prox_tt_ma_1999_1000m, by = "id_unico")|>
  inner_join(shape_tt_ma_1999_1000m, by = "id_unico")

#### Unificando tabelas 1000m ####

lsm_tt_ma_1000m <- bind_rows(lsm_tt_ma_2023_1000m, lsm_tt_ma_2022_1000m, lsm_tt_ma_2021_1000m, lsm_tt_ma_2020_1000m, lsm_tt_ma_2019_1000m, lsm_tt_ma_2018_1000m, lsm_tt_ma_2017_1000m, lsm_tt_ma_2016_1000m, lsm_tt_ma_2015_1000m, lsm_tt_ma_2014_1000m, lsm_tt_ma_2013_1000m, lsm_tt_ma_2012_1000m, lsm_tt_ma_2011_1000m, lsm_tt_ma_2010_1000m, lsm_tt_ma_2009_1000m, lsm_tt_ma_2008_1000m, lsm_tt_ma_2007_1000m, lsm_tt_ma_2006_1000m, lsm_tt_ma_2005_1000m, lsm_tt_ma_2004_1000m, lsm_tt_ma_2003_1000m, lsm_tt_ma_2002_1000m, lsm_tt_ma_2001_1000m, lsm_tt_ma_2000_1000m, lsm_tt_ma_1999_1000m) |>
  dplyr::select(-pland_0,-lpi_0, -ed_0, -pd_0, -np_0, -lpi_12, -ed_12, -pd_12, -np_12, -lpi_15, -ed_15, -pd_15, -np_15, -lpi_9, -ed_9, -pd_9, -np_9, -np_24, -pd_24, -lpi_24, -ed_24,  -pland_24)|>
  mutate(Bin = str_extract(id_unico, "^[01]"),
         Bin = as.factor(Bin),     
         Ano = as.factor(Ano))|>
  mutate(pland_9 = ifelse(is.na(pland_9), 0, pland_9))|>
  mutate(pland_15 = ifelse(is.na(pland_15), 0, pland_15))|>
  mutate(pland_12 = ifelse(is.na(pland_12), 0, pland_12)) |>
  rename_with(~ paste0(., "_1000m"))

lsm_tt_ma_1000m_sem_na <- lsm_tt_ma_1000m|>
  na.omit() 

View(lsm_tt_ma_1000m_sem_na)

write_xlsx(lsm_tt_ma_1000m, "lsm_tt_ma_1000m.xlsx")

lsm_tt_ma_1000m <- read_excel("D:/Analise estatistica mestrado/lsm_tt_ma_1000m.xlsx", sheet = 1,na = c(" ", "NA")) 

#### Gerando buffer 2000m ####

bf_tt_ma_2023_2000m <- tt_ma_flo|>
  filter(Ano == "2023") |>
  st_buffer(dist = 1000)

bf_tt_ma_2022_2000m <- tt_ma_flo|>
  filter(Ano == "2022") |>
  st_buffer(dist = 1000)

bf_tt_ma_2021_2000m <- tt_ma_flo|>
  filter(Ano == "2021") |>
  st_buffer(dist = 1000)

bf_tt_ma_2020_2000m <- tt_ma_flo|>
  filter(Ano == "2020") |>
  st_buffer(dist = 1000)

bf_tt_ma_2019_2000m <- tt_ma_flo|>
  filter(Ano == "2019") |>
  st_buffer(dist = 1000)

bf_tt_ma_2018_2000m <- tt_ma_flo|>
  filter(Ano == "2018") |>
  st_buffer(dist = 1000)

bf_tt_ma_2017_2000m <- tt_ma_flo|>
  filter(Ano == "2017") |>
  st_buffer(dist = 1000)

bf_tt_ma_2016_2000m <- tt_ma_flo|>
  filter(Ano == "2016") |>
  st_buffer(dist = 1000)

bf_tt_ma_2015_2000m <- tt_ma_flo|>
  filter(Ano == "2015") |>
  st_buffer(dist = 1000)

bf_tt_ma_2014_2000m <- tt_ma_flo|>
  filter(Ano == "2014") |>
  st_buffer(dist = 1000)

bf_tt_ma_2013_2000m <- tt_ma_flo|>
  filter(Ano == "2013") |>
  st_buffer(dist = 1000)

bf_tt_ma_2012_2000m <- tt_ma_flo|>
  filter(Ano == "2012") |>
  st_buffer(dist = 1000)

bf_tt_ma_2010_2000m <- tt_ma_flo|>
  filter(Ano == "2010") |>
  st_buffer(dist = 1000)

bf_tt_ma_2011_2000m <- tt_ma_flo|>
  filter(Ano == "2011") |>
  st_buffer(dist = 1000)

bf_tt_ma_2009_2000m <- tt_ma_flo|>
  filter(Ano == "2009") |>
  st_buffer(dist = 1000)

bf_tt_ma_2008_2000m <- tt_ma_flo|>
  filter(Ano == "2008") |>
  st_buffer(dist = 1000)

bf_tt_ma_2007_2000m <- tt_ma_flo|>
  filter(Ano == "2007") |>
  st_buffer(dist = 1000)

bf_tt_ma_2006_2000m <- tt_ma_flo|>
  filter(Ano == "2006") |>
  st_buffer(dist = 1000)

bf_tt_ma_2005_2000m <- tt_ma_flo|>
  filter(Ano == "2005") |>
  st_buffer(dist = 1000)

bf_tt_ma_2004_2000m <- tt_ma_flo|>
  filter(Ano == "2004") |>
  st_buffer(dist = 1000)


#### Cortando raster 2000m ####

# 2023

ras_tt_ma_2023_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2023_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2023_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2023_2000m <- crop(ma_2023, bf_tt_ma_2023_2000m[i, ])
  mask_tt_ma_2023_2000m <- mask(crop_tt_ma_2023_2000m, bf_tt_ma_2023_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2023_2000m[[buffer_id]] <- mask_tt_ma_2023_2000m
}

output_dir <- "ras_tt_ma_2023_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2023_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2023_2000m.tif"))
  writeRaster(
    ras_tt_ma_2023_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2022

ras_tt_ma_2022_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2022_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2022_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2022_2000m <- crop(ma_2022, bf_tt_ma_2022_2000m[i, ])
  mask_tt_ma_2022_2000m <- mask(crop_tt_ma_2022_2000m, bf_tt_ma_2022_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2022_2000m[[buffer_id]] <- mask_tt_ma_2022_2000m
}

output_dir <- "ras_tt_ma_2022_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2022_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2022_2000m.tif"))
  writeRaster(
    ras_tt_ma_2022_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2021

ras_tt_ma_2021_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2021_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2021_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2021_2000m <- crop(ma_2021, bf_tt_ma_2021_2000m[i, ])
  mask_tt_ma_2021_2000m <- mask(crop_tt_ma_2021_2000m, bf_tt_ma_2021_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2021_2000m[[buffer_id]] <- mask_tt_ma_2021_2000m
}

output_dir <- "ras_tt_ma_2021_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2021_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2021_2000m.tif"))
  writeRaster(
    ras_tt_ma_2021_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2020

ras_tt_ma_2020_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2020_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2020_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2020_2000m <- crop(ma_2020, bf_tt_ma_2020_2000m[i, ])
  mask_tt_ma_2020_2000m <- mask(crop_tt_ma_2020_2000m, bf_tt_ma_2020_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2020_2000m[[buffer_id]] <- mask_tt_ma_2020_2000m
}

output_dir <- "ras_tt_ma_2020_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2020_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2020_2000m.tif"))
  writeRaster(
    ras_tt_ma_2020_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2019

ras_tt_ma_2019_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2019_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2019_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2019_2000m <- crop(ma_2019, bf_tt_ma_2019_2000m[i, ])
  mask_tt_ma_2019_2000m <- mask(crop_tt_ma_2019_2000m, bf_tt_ma_2019_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2019_2000m[[buffer_id]] <- mask_tt_ma_2019_2000m
}

output_dir <- "ras_tt_ma_2019_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2019_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2019_2000m.tif"))
  writeRaster(
    ras_tt_ma_2019_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2018

ras_tt_ma_2018_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2018_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2018_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2018_2000m <- crop(ma_2018, bf_tt_ma_2018_2000m[i, ])
  mask_tt_ma_2018_2000m <- mask(crop_tt_ma_2018_2000m, bf_tt_ma_2018_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2018_2000m[[buffer_id]] <- mask_tt_ma_2018_2000m
}

output_dir <- "ras_tt_ma_2018_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2018_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2018_2000m.tif"))
  writeRaster(
    ras_tt_ma_2018_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2017

ras_tt_ma_2017_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2017_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2017_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2017_2000m <- crop(ma_2017, bf_tt_ma_2017_2000m[i, ])
  mask_tt_ma_2017_2000m <- mask(crop_tt_ma_2017_2000m, bf_tt_ma_2017_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2017_2000m[[buffer_id]] <- mask_tt_ma_2017_2000m
}

output_dir <- "ras_tt_ma_2017_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2017_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2017_2000m.tif"))
  writeRaster(
    ras_tt_ma_2017_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2017

ras_tt_ma_2017_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2017_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2017_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2017_2000m <- crop(ma_2017, bf_tt_ma_2017_2000m[i, ])
  mask_tt_ma_2017_2000m <- mask(crop_tt_ma_2017_2000m, bf_tt_ma_2017_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2017_2000m[[buffer_id]] <- mask_tt_ma_2017_2000m
}

output_dir <- "ras_tt_ma_2017_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2017_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2017_2000m.tif"))
  writeRaster(
    ras_tt_ma_2017_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2016

ras_tt_ma_2016_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2016_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2016_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2016_2000m <- crop(ma_2016, bf_tt_ma_2016_2000m[i, ])
  mask_tt_ma_2016_2000m <- mask(crop_tt_ma_2016_2000m, bf_tt_ma_2016_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2016_2000m[[buffer_id]] <- mask_tt_ma_2016_2000m
}

output_dir <- "ras_tt_ma_2016_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2016_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2016_2000m.tif"))
  writeRaster(
    ras_tt_ma_2016_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2015

ras_tt_ma_2015_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2015_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2015_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2015_2000m <- crop(ma_2015, bf_tt_ma_2015_2000m[i, ])
  mask_tt_ma_2015_2000m <- mask(crop_tt_ma_2015_2000m, bf_tt_ma_2015_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2015_2000m[[buffer_id]] <- mask_tt_ma_2015_2000m
}

output_dir <- "ras_tt_ma_2015_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2015_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2015_2000m.tif"))
  writeRaster(
    ras_tt_ma_2015_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2014

ras_tt_ma_2014_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2014_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2014_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2014_2000m <- crop(ma_2014, bf_tt_ma_2014_2000m[i, ])
  mask_tt_ma_2014_2000m <- mask(crop_tt_ma_2014_2000m, bf_tt_ma_2014_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2014_2000m[[buffer_id]] <- mask_tt_ma_2014_2000m
}

output_dir <- "ras_tt_ma_2014_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2014_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2014_2000m.tif"))
  writeRaster(
    ras_tt_ma_2014_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2013

ras_tt_ma_2013_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2013_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2013_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2013_2000m <- crop(ma_2013, bf_tt_ma_2013_2000m[i, ])
  mask_tt_ma_2013_2000m <- mask(crop_tt_ma_2013_2000m, bf_tt_ma_2013_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2013_2000m[[buffer_id]] <- mask_tt_ma_2013_2000m
}

output_dir <- "ras_tt_ma_2013_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2013_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2013_2000m.tif"))
  writeRaster(
    ras_tt_ma_2013_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2012

ras_tt_ma_2012_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2012_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2012_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2012_2000m <- crop(ma_2012, bf_tt_ma_2012_2000m[i, ])
  mask_tt_ma_2012_2000m <- mask(crop_tt_ma_2012_2000m, bf_tt_ma_2012_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2012_2000m[[buffer_id]] <- mask_tt_ma_2012_2000m
}

output_dir <- "ras_tt_ma_2012_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2012_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2012_2000m.tif"))
  writeRaster(
    ras_tt_ma_2012_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2011

ras_tt_ma_2011_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2011_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2011_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2011_2000m <- crop(ma_2011, bf_tt_ma_2011_2000m[i, ])
  mask_tt_ma_2011_2000m <- mask(crop_tt_ma_2011_2000m, bf_tt_ma_2011_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2011_2000m[[buffer_id]] <- mask_tt_ma_2011_2000m
}

output_dir <- "ras_tt_ma_2011_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2011_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2011_2000m.tif"))
  writeRaster(
    ras_tt_ma_2011_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2010

ras_tt_ma_2010_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2010_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2010_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2010_2000m <- crop(ma_2010, bf_tt_ma_2010_2000m[i, ])
  mask_tt_ma_2010_2000m <- mask(crop_tt_ma_2010_2000m, bf_tt_ma_2010_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2010_2000m[[buffer_id]] <- mask_tt_ma_2010_2000m
}

output_dir <- "ras_tt_ma_2010_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2010_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2010_2000m.tif"))
  writeRaster(
    ras_tt_ma_2010_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2009

ras_tt_ma_2009_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2009_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2009_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2009_2000m <- crop(ma_2009, bf_tt_ma_2009_2000m[i, ])
  mask_tt_ma_2009_2000m <- mask(crop_tt_ma_2009_2000m, bf_tt_ma_2009_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2009_2000m[[buffer_id]] <- mask_tt_ma_2009_2000m
}

output_dir <- "ras_tt_ma_2009_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2009_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2009_2000m.tif"))
  writeRaster(
    ras_tt_ma_2009_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2008

ras_tt_ma_2008_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2008_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2008_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2008_2000m <- crop(ma_2008, bf_tt_ma_2008_2000m[i, ])
  mask_tt_ma_2008_2000m <- mask(crop_tt_ma_2008_2000m, bf_tt_ma_2008_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2008_2000m[[buffer_id]] <- mask_tt_ma_2008_2000m
}

output_dir <- "ras_tt_ma_2008_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2008_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2008_2000m.tif"))
  writeRaster(
    ras_tt_ma_2008_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2007

ras_tt_ma_2007_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2007_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2007_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2007_2000m <- crop(ma_2007, bf_tt_ma_2007_2000m[i, ])
  mask_tt_ma_2007_2000m <- mask(crop_tt_ma_2007_2000m, bf_tt_ma_2007_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2007_2000m[[buffer_id]] <- mask_tt_ma_2007_2000m
}

output_dir <- "ras_tt_ma_2007_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2007_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2007_2000m.tif"))
  writeRaster(
    ras_tt_ma_2007_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2006

ras_tt_ma_2006_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2006_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2006_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2006_2000m <- crop(ma_2006, bf_tt_ma_2006_2000m[i, ])
  mask_tt_ma_2006_2000m <- mask(crop_tt_ma_2006_2000m, bf_tt_ma_2006_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2006_2000m[[buffer_id]] <- mask_tt_ma_2006_2000m
}

output_dir <- "ras_tt_ma_2006_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2006_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2006_2000m.tif"))
  writeRaster(
    ras_tt_ma_2006_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2005

ras_tt_ma_2005_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2005_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2005_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2005_2000m <- crop(ma_2005, bf_tt_ma_2005_2000m[i, ])
  mask_tt_ma_2005_2000m <- mask(crop_tt_ma_2005_2000m, bf_tt_ma_2005_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2005_2000m[[buffer_id]] <- mask_tt_ma_2005_2000m
}

output_dir <- "ras_tt_ma_2005_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2005_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2005_2000m.tif"))
  writeRaster(
    ras_tt_ma_2005_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2004

ras_tt_ma_2004_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2004_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2004_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2004_2000m <- crop(ma_2004, bf_tt_ma_2004_2000m[i, ])
  mask_tt_ma_2004_2000m <- mask(crop_tt_ma_2004_2000m, bf_tt_ma_2004_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2004_2000m[[buffer_id]] <- mask_tt_ma_2004_2000m
}

output_dir <- "ras_tt_ma_2004_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2004_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2004_2000m.tif"))
  writeRaster(
    ras_tt_ma_2004_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}



# Chamando os recortes

output_dir <- "ras_tt_ma_2023_2000m"
ras_tt_ma_2023_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2023_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2022_2000m"
ras_tt_ma_2022_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2022_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2021_2000m"
ras_tt_ma_2021_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2021_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2020_2000m"
ras_tt_ma_2020_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2020_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2019_2000m"
ras_tt_ma_2019_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2019_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2018_2000m"
ras_tt_ma_2018_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2018_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2017_2000m"
ras_tt_ma_2017_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2017_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2015_2000m"
ras_tt_ma_2015_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2015_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2014_2000m"
ras_tt_ma_2014_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2014_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2013_2000m"
ras_tt_ma_2013_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2013_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2012_2000m"
ras_tt_ma_2012_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2012_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2011_2000m"
ras_tt_ma_2011_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2011_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2009_2000m"
ras_tt_ma_2009_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2009_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2008_2000m"
ras_tt_ma_2008_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2008_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2007_2000m"
ras_tt_ma_2007_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2007_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2006_2000m"
ras_tt_ma_2006_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2006_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2005_2000m"
ras_tt_ma_2005_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2005_2000m\\.tif$")) %>%
  map(rast) 


#### Metricas de paisagem 2000m ####

# 2023

id_unico <- names(ras_tt_ma_2023_2000m)

met_tt_ma_2023_2000m <- map_df(seq_along(ras_tt_ma_2023_2000m), function(i) {
  raster <- ras_tt_ma_2023_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2023_2000m <- prox(ras_tt_ma_2023_2000m, 3)

shape_tt_ma_2023_2000m <- shape(ras_tt_ma_2023_2000m, class_value = 3)

lsm_tt_ma_2023_2000m <- met_tt_ma_2023_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2023")))|>
  inner_join(prox_tt_ma_2023_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2023_2000m, by = "id_unico")

# 2022

id_unico <- names(ras_tt_ma_2022_2000m)

met_tt_ma_2022_2000m <- map_df(seq_along(ras_tt_ma_2022_2000m), function(i) {
  raster <- ras_tt_ma_2022_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2022_2000m <- prox(ras_tt_ma_2022_2000m, 3)

shape_tt_ma_2022_2000m <- shape(ras_tt_ma_2022_2000m, class_value = 3)

lsm_tt_ma_2022_2000m <- met_tt_ma_2022_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2022")))|>
  inner_join(prox_tt_ma_2022_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2022_2000m, by = "id_unico")

# 2021

id_unico <- names(ras_tt_ma_2021_2000m)

met_tt_ma_2021_2000m <- map_df(seq_along(ras_tt_ma_2021_2000m), function(i) {
  raster <- ras_tt_ma_2021_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2021_2000m <- prox(ras_tt_ma_2021_2000m, 3)

shape_tt_ma_2021_2000m <- shape(ras_tt_ma_2021_2000m, class_value = 3)

lsm_tt_ma_2021_2000m <- met_tt_ma_2021_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2021")))|>
  inner_join(prox_tt_ma_2021_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2021_2000m, by = "id_unico")

# 2020

id_unico <- names(ras_tt_ma_2020_2000m)

met_tt_ma_2020_2000m <- map_df(seq_along(ras_tt_ma_2020_2000m), function(i) {
  raster <- ras_tt_ma_2020_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2020_2000m <- prox(ras_tt_ma_2020_2000m, 3)

shape_tt_ma_2020_2000m <- shape(ras_tt_ma_2020_2000m, class_value = 3)

lsm_tt_ma_2020_2000m <- met_tt_ma_2020_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2020")))|>
  inner_join(prox_tt_ma_2020_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2020_2000m, by = "id_unico")

# 2019

id_unico <- names(ras_tt_ma_2019_2000m)

met_tt_ma_2019_2000m <- map_df(seq_along(ras_tt_ma_2019_2000m), function(i) {
  raster <- ras_tt_ma_2019_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2019_2000m <- prox(ras_tt_ma_2019_2000m, 3)

shape_tt_ma_2019_2000m <- shape(ras_tt_ma_2019_2000m, class_value = 3)

lsm_tt_ma_2019_2000m <- met_tt_ma_2019_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2019")))|>
  inner_join(prox_tt_ma_2019_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2019_2000m, by = "id_unico")

# 2018

id_unico <- names(ras_tt_ma_2018_2000m)

met_tt_ma_2018_2000m <- map_df(seq_along(ras_tt_ma_2018_2000m), function(i) {
  raster <- ras_tt_ma_2018_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2018_2000m <- prox(ras_tt_ma_2018_2000m, 3)

shape_tt_ma_2018_2000m <- shape(ras_tt_ma_2018_2000m, class_value = 3)

lsm_tt_ma_2018_2000m <- met_tt_ma_2018_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2018")))|>
  inner_join(prox_tt_ma_2018_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2018_2000m, by = "id_unico")

# 2017

id_unico <- names(ras_tt_ma_2017_2000m)

met_tt_ma_2017_2000m <- map_df(seq_along(ras_tt_ma_2017_2000m), function(i) {
  raster <- ras_tt_ma_2017_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2017_2000m <- prox(ras_tt_ma_2017_2000m, 3)

shape_tt_ma_2017_2000m <- shape(ras_tt_ma_2017_2000m, class_value = 3)

lsm_tt_ma_2017_2000m <- met_tt_ma_2017_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2017")))|>
  inner_join(prox_tt_ma_2017_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2017_2000m, by = "id_unico")


# 2016

id_unico <- names(ras_tt_ma_2016_2000m)

met_tt_ma_2016_2000m <- map_df(seq_along(ras_tt_ma_2016_2000m), function(i) {
  raster <- ras_tt_ma_2016_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2016_2000m <- prox(ras_tt_ma_2016_2000m, 3)

shape_tt_ma_2016_2000m <- shape(ras_tt_ma_2016_2000m, class_value = 3)

lsm_tt_ma_2016_2000m <- met_tt_ma_2016_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2016")))|>
  inner_join(prox_tt_ma_2016_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2016_2000m, by = "id_unico")

# 2015

id_unico <- names(ras_tt_ma_2015_2000m)

met_tt_ma_2015_2000m <- map_df(seq_along(ras_tt_ma_2015_2000m), function(i) {
  raster <- ras_tt_ma_2015_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2015_2000m <- prox(ras_tt_ma_2015_2000m, 3)

shape_tt_ma_2015_2000m <- shape(ras_tt_ma_2015_2000m, class_value = 3)

lsm_tt_ma_2015_2000m <- met_tt_ma_2015_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2015")))|>
  inner_join(prox_tt_ma_2015_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2015_2000m, by = "id_unico")

# 2014

id_unico <- names(ras_tt_ma_2014_2000m)

met_tt_ma_2014_2000m <- map_df(seq_along(ras_tt_ma_2014_2000m), function(i) {
  raster <- ras_tt_ma_2014_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2014_2000m <- prox(ras_tt_ma_2014_2000m, 3)

shape_tt_ma_2014_2000m <- shape(ras_tt_ma_2014_2000m, class_value = 3)

lsm_tt_ma_2014_2000m <- met_tt_ma_2014_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2014")))|>
  inner_join(prox_tt_ma_2014_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2014_2000m, by = "id_unico")

# 2013

id_unico <- names(ras_tt_ma_2013_2000m)

met_tt_ma_2013_2000m <- map_df(seq_along(ras_tt_ma_2013_2000m), function(i) {
  raster <- ras_tt_ma_2013_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2013_2000m <- prox(ras_tt_ma_2013_2000m, 3)

shape_tt_ma_2013_2000m <- shape(ras_tt_ma_2013_2000m, class_value = 3)

lsm_tt_ma_2013_2000m <- met_tt_ma_2013_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2013")))|>
  inner_join(prox_tt_ma_2013_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2013_2000m, by = "id_unico")

# 2012
id_unico <- names(ras_tt_ma_2012_2000m)

met_tt_ma_2012_2000m <- map_df(seq_along(ras_tt_ma_2012_2000m), function(i) {
  raster <- ras_tt_ma_2012_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2012_2000m <- prox(ras_tt_ma_2012_2000m, 3)

shape_tt_ma_2012_2000m <- shape(ras_tt_ma_2012_2000m, class_value = 3)

lsm_tt_ma_2012_2000m <- met_tt_ma_2012_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2012")))|>
  inner_join(prox_tt_ma_2012_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2012_2000m, by = "id_unico")

# 2011

id_unico <- names(ras_tt_ma_2011_2000m)

met_tt_ma_2011_2000m <- map_df(seq_along(ras_tt_ma_2011_2000m), function(i) {
  raster <- ras_tt_ma_2011_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2011_2000m <- prox(ras_tt_ma_2011_2000m, 3)

shape_tt_ma_2011_2000m <- shape(ras_tt_ma_2011_2000m, class_value = 3)

lsm_tt_ma_2011_2000m <- met_tt_ma_2011_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2011")))|>
  inner_join(prox_tt_ma_2011_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2011_2000m, by = "id_unico")
# 2010

id_unico <- names(ras_tt_ma_2010_2000m)

met_tt_ma_2010_2000m <- map_df(seq_along(ras_tt_ma_2010_2000m), function(i) {
  raster <- ras_tt_ma_2010_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2010_2000m <- prox(ras_tt_ma_2010_2000m, 3)

shape_tt_ma_2010_2000m <- shape(ras_tt_ma_2010_2000m, class_value = 3)

lsm_tt_ma_2010_2000m <- met_tt_ma_2010_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2010")))|>
  inner_join(prox_tt_ma_2010_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2010_2000m, by = "id_unico")

# 2009
id_unico <- names(ras_tt_ma_2009_2000m)

met_tt_ma_2009_2000m <- map_df(seq_along(ras_tt_ma_2009_2000m), function(i) {
  raster <- ras_tt_ma_2009_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2009_2000m <- prox(ras_tt_ma_2009_2000m, 3)

shape_tt_ma_2009_2000m <- shape(ras_tt_ma_2009_2000m, class_value = 3)

lsm_tt_ma_2009_2000m <- met_tt_ma_2009_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2009")))|>
  inner_join(prox_tt_ma_2009_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2009_2000m, by = "id_unico")

# 2008

id_unico <- names(ras_tt_ma_2008_2000m)

met_tt_ma_2008_2000m <- map_df(seq_along(ras_tt_ma_2008_2000m), function(i) {
  raster <- ras_tt_ma_2008_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2008_2000m <- prox(ras_tt_ma_2008_2000m, 3)

shape_tt_ma_2008_2000m <- shape(ras_tt_ma_2008_2000m, class_value = 3)

lsm_tt_ma_2008_2000m <- met_tt_ma_2008_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2008")))|>
  inner_join(prox_tt_ma_2008_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2008_2000m, by = "id_unico")

# 2007

id_unico <- names(ras_tt_ma_2007_2000m)

met_tt_ma_2007_2000m <- map_df(seq_along(ras_tt_ma_2007_2000m), function(i) {
  raster <- ras_tt_ma_2007_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2007_2000m <- prox(ras_tt_ma_2007_2000m, 3)

shape_tt_ma_2007_2000m <- shape(ras_tt_ma_2007_2000m, class_value = 3)

lsm_tt_ma_2007_2000m <- met_tt_ma_2007_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2007")))|>
  inner_join(prox_tt_ma_2007_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2007_2000m, by = "id_unico")

# 2006

id_unico <- names(ras_tt_ma_2006_2000m)

met_tt_ma_2006_2000m <- map_df(seq_along(ras_tt_ma_2006_2000m), function(i) {
  raster <- ras_tt_ma_2006_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2006_2000m <- prox(ras_tt_ma_2006_2000m, 3)

shape_tt_ma_2006_2000m <- shape(ras_tt_ma_2006_2000m, class_value = 3)

lsm_tt_ma_2006_2000m <- met_tt_ma_2006_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2006")))|>
  inner_join(prox_tt_ma_2006_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2006_2000m, by = "id_unico")

# 2005

id_unico <- names(ras_tt_ma_2005_2000m)

met_tt_ma_2005_2000m <- map_df(seq_along(ras_tt_ma_2005_2000m), function(i) {
  raster <- ras_tt_ma_2005_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2005_2000m <- prox(ras_tt_ma_2005_2000m, 3)

shape_tt_ma_2005_2000m <- shape(ras_tt_ma_2005_2000m, class_value = 3)

lsm_tt_ma_2005_2000m <- met_tt_ma_2005_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2005")))|>
  inner_join(prox_tt_ma_2005_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2005_2000m, by = "id_unico")

# 2004

id_unico <- names(ras_tt_ma_2004_2000m)

met_tt_ma_2004_2000m <- map_df(seq_along(ras_tt_ma_2004_2000m), function(i) {
  raster <- ras_tt_ma_2004_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2004_2000m <- prox(ras_tt_ma_2004_2000m, 3)

shape_tt_ma_2004_2000m <- shape(ras_tt_ma_2004_2000m, class_value = 3)

lsm_tt_ma_2004_2000m <- met_tt_ma_2004_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2004")))|>
  inner_join(prox_tt_ma_2004_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2004_2000m, by = "id_unico")


#### Gerando buffer 2000m ####

bf_tt_ma_2023_2000m <-tt_ma_flo |>
  filter(Ano == "2023") |>
  st_buffer(dist = 2000)

bf_tt_ma_2022_2000m <-tt_ma_flo |>
  filter(Ano == "2022") |>
  st_buffer(dist = 2000)

bf_tt_ma_2021_2000m <-tt_ma_flo |>
  filter(Ano == "2021") |>
  st_buffer(dist = 2000)

bf_tt_ma_2020_2000m <-tt_ma_flo |>
  filter(Ano == "2020") |>
  st_buffer(dist = 2000)

bf_tt_ma_2019_2000m <-tt_ma_flo |>
  filter(Ano == "2019") |>
  st_buffer(dist = 2000)

bf_tt_ma_2018_2000m <-tt_ma_flo |>
  filter(Ano == "2018") |>
  st_buffer(dist = 2000)

bf_tt_ma_2017_2000m <-tt_ma_flo |>
  filter(Ano == "2017") |>
  st_buffer(dist = 2000)

bf_tt_ma_2016_2000m <-tt_ma_flo |>
  filter(Ano == "2016") |>
  st_buffer(dist = 2000)

bf_tt_ma_2015_2000m <-tt_ma_flo |>
  filter(Ano == "2015") |>
  st_buffer(dist = 2000)

bf_tt_ma_2014_2000m <-tt_ma_flo |>
  filter(Ano == "2014") |>
  st_buffer(dist = 2000)

bf_tt_ma_2013_2000m <-tt_ma_flo |>
  filter(Ano == "2013") |>
  st_buffer(dist = 2000)

bf_tt_ma_2012_2000m <-tt_ma_flo |>
  filter(Ano == "2012") |>
  st_buffer(dist = 2000)

bf_tt_ma_2011_2000m <-tt_ma_flo |>
  filter(Ano == "2011") |>
  st_buffer(dist = 2000)

bf_tt_ma_2010_2000m <-tt_ma_flo |>
  filter(Ano == "2010") |>
  st_buffer(dist = 2000)

bf_tt_ma_2009_2000m <-tt_ma_flo |>
  filter(Ano == "2009") |>
  st_buffer(dist = 2000)

bf_tt_ma_2008_2000m <-tt_ma_flo |>
  filter(Ano == "2008") |>
  st_buffer(dist = 2000)

bf_tt_ma_2007_2000m <-tt_ma_flo |>
  filter(Ano == "2007") |>
  st_buffer(dist = 2000)

bf_tt_ma_2006_2000m <-tt_ma_flo |>
  filter(Ano == "2006") |>
  st_buffer(dist = 2000)

bf_tt_ma_2005_2000m <-tt_ma_flo |>
  filter(Ano == "2005") |>
  st_buffer(dist = 2000)

bf_tt_ma_2004_2000m <-tt_ma_flo |>
  filter(Ano == "2004") |>
  st_buffer(dist = 2000)

bf_tt_ma_2003_2000m <-tt_ma_flo |>
  filter(Ano == "2003") |>
  st_buffer(dist = 2000)

bf_tt_ma_2002_2000m <-tt_ma_flo |>
  filter(Ano == "2002") |>
  st_buffer(dist = 2000)

bf_tt_ma_2001_2000m <-tt_ma_flo |>
  filter(Ano == "2001") |>
  st_buffer(dist = 2000)

bf_tt_ma_2000_2000m <-tt_ma_flo |>
  filter(Ano == "2000") |>
  st_buffer(dist = 2000)

bf_tt_ma_1999_2000m <-tt_ma_flo |>
  filter(Ano == "1999") |>
  st_buffer(dist = 2000)

#### Cortando raster 2000m ####

# 2023

ras_tt_ma_2023_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2023_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2023_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2023_2000m <- crop(ma_2023, bf_tt_ma_2023_2000m[i, ])
  mask_tt_ma_2023_2000m <- mask(crop_tt_ma_2023_2000m, bf_tt_ma_2023_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2023_2000m[[buffer_id]] <- mask_tt_ma_2023_2000m
}

output_dir <- "ras_tt_ma_2023_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2023_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2023_2000m.tif"))
  writeRaster(
    ras_tt_ma_2023_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2022

ras_tt_ma_2022_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2022_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2022_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2022_2000m <- crop(ma_2022, bf_tt_ma_2022_2000m[i, ])
  mask_tt_ma_2022_2000m <- mask(crop_tt_ma_2022_2000m, bf_tt_ma_2022_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2022_2000m[[buffer_id]] <- mask_tt_ma_2022_2000m
}

output_dir <- "ras_tt_ma_2022_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2022_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2022_2000m.tif"))
  writeRaster(
    ras_tt_ma_2022_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2021

ras_tt_ma_2021_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2021_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2021_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2021_2000m <- crop(ma_2021, bf_tt_ma_2021_2000m[i, ])
  mask_tt_ma_2021_2000m <- mask(crop_tt_ma_2021_2000m, bf_tt_ma_2021_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2021_2000m[[buffer_id]] <- mask_tt_ma_2021_2000m
}

output_dir <- "ras_tt_ma_2021_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2021_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2021_2000m.tif"))
  writeRaster(
    ras_tt_ma_2021_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2020

ras_tt_ma_2020_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2020_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2020_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2020_2000m <- crop(ma_2020, bf_tt_ma_2020_2000m[i, ])
  mask_tt_ma_2020_2000m <- mask(crop_tt_ma_2020_2000m, bf_tt_ma_2020_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2020_2000m[[buffer_id]] <- mask_tt_ma_2020_2000m
}

output_dir <- "ras_tt_ma_2020_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2020_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2020_2000m.tif"))
  writeRaster(
    ras_tt_ma_2020_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2019

ras_tt_ma_2019_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2019_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2019_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2019_2000m <- crop(ma_2019, bf_tt_ma_2019_2000m[i, ])
  mask_tt_ma_2019_2000m <- mask(crop_tt_ma_2019_2000m, bf_tt_ma_2019_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2019_2000m[[buffer_id]] <- mask_tt_ma_2019_2000m
}

output_dir <- "ras_tt_ma_2019_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2019_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2019_2000m.tif"))
  writeRaster(
    ras_tt_ma_2019_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2018

ras_tt_ma_2018_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2018_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2018_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2018_2000m <- crop(ma_2018, bf_tt_ma_2018_2000m[i, ])
  mask_tt_ma_2018_2000m <- mask(crop_tt_ma_2018_2000m, bf_tt_ma_2018_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2018_2000m[[buffer_id]] <- mask_tt_ma_2018_2000m
}

output_dir <- "ras_tt_ma_2018_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2018_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2018_2000m.tif"))
  writeRaster(
    ras_tt_ma_2018_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2017

ras_tt_ma_2017_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2017_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2017_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2017_2000m <- crop(ma_2017, bf_tt_ma_2017_2000m[i, ])
  mask_tt_ma_2017_2000m <- mask(crop_tt_ma_2017_2000m, bf_tt_ma_2017_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2017_2000m[[buffer_id]] <- mask_tt_ma_2017_2000m
}

output_dir <- "ras_tt_ma_2017_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2017_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2017_2000m.tif"))
  writeRaster(
    ras_tt_ma_2017_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2016

ras_tt_ma_2016_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2016_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2016_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2016_2000m <- crop(ma_2016, bf_tt_ma_2016_2000m[i, ])
  mask_tt_ma_2016_2000m <- mask(crop_tt_ma_2016_2000m, bf_tt_ma_2016_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2016_2000m[[buffer_id]] <- mask_tt_ma_2016_2000m
}

output_dir <- "ras_tt_ma_2016_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2016_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2016_2000m.tif"))
  writeRaster(
    ras_tt_ma_2016_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2015

ras_tt_ma_2015_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2015_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2015_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2015_2000m <- crop(ma_2015, bf_tt_ma_2015_2000m[i, ])
  mask_tt_ma_2015_2000m <- mask(crop_tt_ma_2015_2000m, bf_tt_ma_2015_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2015_2000m[[buffer_id]] <- mask_tt_ma_2015_2000m
}

output_dir <- "ras_tt_ma_2015_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2015_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2015_2000m.tif"))
  writeRaster(
    ras_tt_ma_2015_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2014

ras_tt_ma_2014_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2014_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2014_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2014_2000m <- crop(ma_2014, bf_tt_ma_2014_2000m[i, ])
  mask_tt_ma_2014_2000m <- mask(crop_tt_ma_2014_2000m, bf_tt_ma_2014_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2014_2000m[[buffer_id]] <- mask_tt_ma_2014_2000m
}

output_dir <- "ras_tt_ma_2014_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2014_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2014_2000m.tif"))
  writeRaster(
    ras_tt_ma_2014_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2013

ras_tt_ma_2013_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2013_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2013_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2013_2000m <- crop(ma_2013, bf_tt_ma_2013_2000m[i, ])
  mask_tt_ma_2013_2000m <- mask(crop_tt_ma_2013_2000m, bf_tt_ma_2013_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2013_2000m[[buffer_id]] <- mask_tt_ma_2013_2000m
}

output_dir <- "ras_tt_ma_2013_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2013_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2013_2000m.tif"))
  writeRaster(
    ras_tt_ma_2013_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2012

ras_tt_ma_2012_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2012_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2012_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2012_2000m <- crop(ma_2012, bf_tt_ma_2012_2000m[i, ])
  mask_tt_ma_2012_2000m <- mask(crop_tt_ma_2012_2000m, bf_tt_ma_2012_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2012_2000m[[buffer_id]] <- mask_tt_ma_2012_2000m
}

output_dir <- "ras_tt_ma_2012_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2012_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2012_2000m.tif"))
  writeRaster(
    ras_tt_ma_2012_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2011

ras_tt_ma_2011_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2011_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2011_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2011_2000m <- crop(ma_2011, bf_tt_ma_2011_2000m[i, ])
  mask_tt_ma_2011_2000m <- mask(crop_tt_ma_2011_2000m, bf_tt_ma_2011_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2011_2000m[[buffer_id]] <- mask_tt_ma_2011_2000m
}

output_dir <- "ras_tt_ma_2011_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2011_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2011_2000m.tif"))
  writeRaster(
    ras_tt_ma_2011_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2010

ras_tt_ma_2010_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2010_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2010_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2010_2000m <- crop(ma_2010, bf_tt_ma_2010_2000m[i, ])
  mask_tt_ma_2010_2000m <- mask(crop_tt_ma_2010_2000m, bf_tt_ma_2010_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2010_2000m[[buffer_id]] <- mask_tt_ma_2010_2000m
}

output_dir <- "ras_tt_ma_2010_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2010_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2010_2000m.tif"))
  writeRaster(
    ras_tt_ma_2010_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2009

ras_tt_ma_2009_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2009_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2009_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2009_2000m <- crop(ma_2009, bf_tt_ma_2009_2000m[i, ])
  mask_tt_ma_2009_2000m <- mask(crop_tt_ma_2009_2000m, bf_tt_ma_2009_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2009_2000m[[buffer_id]] <- mask_tt_ma_2009_2000m
}

output_dir <- "ras_tt_ma_2009_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2009_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2009_2000m.tif"))
  writeRaster(
    ras_tt_ma_2009_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2008

ras_tt_ma_2008_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2008_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2008_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2008_2000m <- crop(ma_2008, bf_tt_ma_2008_2000m[i, ])
  mask_tt_ma_2008_2000m <- mask(crop_tt_ma_2008_2000m, bf_tt_ma_2008_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2008_2000m[[buffer_id]] <- mask_tt_ma_2008_2000m
}

output_dir <- "ras_tt_ma_2008_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2008_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2008_2000m.tif"))
  writeRaster(
    ras_tt_ma_2008_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2007

ras_tt_ma_2007_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2007_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2007_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2007_2000m <- crop(ma_2007, bf_tt_ma_2007_2000m[i, ])
  mask_tt_ma_2007_2000m <- mask(crop_tt_ma_2007_2000m, bf_tt_ma_2007_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2007_2000m[[buffer_id]] <- mask_tt_ma_2007_2000m
}

output_dir <- "ras_tt_ma_2007_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2007_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2007_2000m.tif"))
  writeRaster(
    ras_tt_ma_2007_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2006

ras_tt_ma_2006_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2006_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2006_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2006_2000m <- crop(ma_2006, bf_tt_ma_2006_2000m[i, ])
  mask_tt_ma_2006_2000m <- mask(crop_tt_ma_2006_2000m, bf_tt_ma_2006_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2006_2000m[[buffer_id]] <- mask_tt_ma_2006_2000m
}

output_dir <- "ras_tt_ma_2006_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2006_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2006_2000m.tif"))
  writeRaster(
    ras_tt_ma_2006_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2005

ras_tt_ma_2005_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2005_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2005_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2005_2000m <- crop(ma_2005, bf_tt_ma_2005_2000m[i, ])
  mask_tt_ma_2005_2000m <- mask(crop_tt_ma_2005_2000m, bf_tt_ma_2005_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2005_2000m[[buffer_id]] <- mask_tt_ma_2005_2000m
}

output_dir <- "ras_tt_ma_2005_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2005_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2005_2000m.tif"))
  writeRaster(
    ras_tt_ma_2005_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2004

ras_tt_ma_2004_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2004_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2004_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2004_2000m <- crop(ma_2004, bf_tt_ma_2004_2000m[i, ])
  mask_tt_ma_2004_2000m <- mask(crop_tt_ma_2004_2000m, bf_tt_ma_2004_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2004_2000m[[buffer_id]] <- mask_tt_ma_2004_2000m
}

output_dir <- "ras_tt_ma_2004_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2004_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2004_2000m.tif"))
  writeRaster(
    ras_tt_ma_2004_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2003

ras_tt_ma_2003_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2003_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2003_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2003_2000m <- crop(ma_2003, bf_tt_ma_2003_2000m[i, ])
  mask_tt_ma_2003_2000m <- mask(crop_tt_ma_2003_2000m, bf_tt_ma_2003_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2003_2000m[[buffer_id]] <- mask_tt_ma_2003_2000m
}

output_dir <- "ras_tt_ma_2003_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2003_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2003_2000m.tif"))
  writeRaster(
    ras_tt_ma_2003_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2002

ras_tt_ma_2002_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2002_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2002_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2002_2000m <- crop(ma_2002, bf_tt_ma_2002_2000m[i, ])
  mask_tt_ma_2002_2000m <- mask(crop_tt_ma_2002_2000m, bf_tt_ma_2002_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2002_2000m[[buffer_id]] <- mask_tt_ma_2002_2000m
}

output_dir <- "ras_tt_ma_2002_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2002_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2002_2000m.tif"))
  writeRaster(
    ras_tt_ma_2002_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2001

ras_tt_ma_2001_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2001_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2001_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2001_2000m <- crop(ma_2001, bf_tt_ma_2001_2000m[i, ])
  mask_tt_ma_2001_2000m <- mask(crop_tt_ma_2001_2000m, bf_tt_ma_2001_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2001_2000m[[buffer_id]] <- mask_tt_ma_2001_2000m
}

output_dir <- "ras_tt_ma_2001_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2001_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2001_2000m.tif"))
  writeRaster(
    ras_tt_ma_2001_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2000

ras_tt_ma_2000_2000m <- list()

for (i in 1:nrow(bf_tt_ma_2000_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_2000_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_2000_2000m <- crop(ma_2000, bf_tt_ma_2000_2000m[i, ])
  mask_tt_ma_2000_2000m <- mask(crop_tt_ma_2000_2000m, bf_tt_ma_2000_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_2000_2000m[[buffer_id]] <- mask_tt_ma_2000_2000m
}

output_dir <- "ras_tt_ma_2000_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_2000_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_2000_2000m.tif"))
  writeRaster(
    ras_tt_ma_2000_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 1999

ras_tt_ma_1999_2000m <- list()

for (i in 1:nrow(bf_tt_ma_1999_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_tt_ma_1999_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_tt_ma_1999_2000m <- crop(ma_1999, bf_tt_ma_1999_2000m[i, ])
  mask_tt_ma_1999_2000m <- mask(crop_tt_ma_1999_2000m, bf_tt_ma_1999_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_tt_ma_1999_2000m[[buffer_id]] <- mask_tt_ma_1999_2000m
}

output_dir <- "ras_tt_ma_1999_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_tt_ma_1999_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_tt_ma_1999_2000m.tif"))
  writeRaster(
    ras_tt_ma_1999_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# Chamando os recortes

output_dir <- "ras_tt_ma_2023_2000m"
ras_tt_ma_2023_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_tt_ma_2023_2000m\\.tif$")) |>
  map(rast) 

output_dir <- "ras_tt_ma_2022_2000m"
ras_tt_ma_2022_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2022_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2021_2000m"
ras_tt_ma_2021_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2021_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2020_2000m"
ras_tt_ma_2020_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2020_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2019_2000m"
ras_tt_ma_2019_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2019_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2018_2000m"
ras_tt_ma_2018_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2018_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2017_2000m"
ras_tt_ma_2017_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2017_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2016_2000m"
ras_tt_ma_2016_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2016_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2015_2000m"
ras_tt_ma_2015_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2015_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2014_2000m"
ras_tt_ma_2014_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2014_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2013_2000m"
ras_tt_ma_2013_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2013_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2012_2000m"
ras_tt_ma_2012_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2012_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2011_2000m"
ras_tt_ma_2011_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2011_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2010_2000m"
ras_tt_ma_2010_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2010_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2009_2000m"
ras_tt_ma_2009_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2009_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2008_2000m"
ras_tt_ma_2008_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2008_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2007_2000m"
ras_tt_ma_2007_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2007_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2006_2000m"
ras_tt_ma_2006_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2006_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2005_2000m"
ras_tt_ma_2005_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_tt_ma_2005_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_tt_ma_2004_2000m"
ras_tt_ma_2004_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_tt_ma_2004_2000m\\.tif$")) |>
  map(rast) 

output_dir <- "ras_tt_ma_2003_2000m"
ras_tt_ma_2003_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_tt_ma_2003_2000m\\.tif$")) |>
  map(rast) 

output_dir <- "ras_tt_ma_2002_2000m"
ras_tt_ma_2002_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_tt_ma_2002_2000m\\.tif$")) |>
  map(rast) 

output_dir <- "ras_tt_ma_2001_2000m"
ras_tt_ma_2001_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_tt_ma_2001_2000m\\.tif$")) |>
  map(rast) 

output_dir <- "ras_tt_ma_2000_2000m"
ras_tt_ma_2000_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_tt_ma_2000_2000m\\.tif$")) |>
  map(rast) 

output_dir <- "ras_tt_ma_1999_2000m"
ras_tt_ma_1999_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_tt_ma_1999_2000m\\.tif$")) |>
  map(rast) 

#### Metricas de paisagem 2000m ####

# 2023

id_unico <- names(ras_tt_ma_2023_2000m)

met_tt_ma_2023_2000m <- map_df(seq_along(ras_tt_ma_2023_2000m), function(i) {
  raster <- ras_tt_ma_2023_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2023_2000m <- prox(ras_tt_ma_2023_2000m, 3)

shape_tt_ma_2023_2000m <- shape(ras_tt_ma_2023_2000m, class_value = 3)

lsm_tt_ma_2023_2000m <- met_tt_ma_2023_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2023")))|>
  inner_join(prox_tt_ma_2023_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2023_2000m, by = "id_unico")



# 2022

id_unico <- names(ras_tt_ma_2022_2000m)

met_tt_ma_2022_2000m <- map_df(seq_along(ras_tt_ma_2022_2000m), function(i) {
  raster <- ras_tt_ma_2022_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2022_2000m <- prox(ras_tt_ma_2022_2000m, 3)

shape_tt_ma_2022_2000m <- shape(ras_tt_ma_2022_2000m, class_value = 3)

lsm_tt_ma_2022_2000m <- met_tt_ma_2022_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2022")))|>
  inner_join(prox_tt_ma_2022_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2022_2000m, by = "id_unico")

# 2021

id_unico <- names(ras_tt_ma_2021_2000m)

met_tt_ma_2021_2000m <- map_df(seq_along(ras_tt_ma_2021_2000m), function(i) {
  raster <- ras_tt_ma_2021_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2021_2000m <- prox(ras_tt_ma_2021_2000m, 3)

shape_tt_ma_2021_2000m <- shape(ras_tt_ma_2021_2000m, class_value = 3)

lsm_tt_ma_2021_2000m <- met_tt_ma_2021_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2021")))|>
  inner_join(prox_tt_ma_2021_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2021_2000m, by = "id_unico")

# 2020

id_unico <- names(ras_tt_ma_2020_2000m)

met_tt_ma_2020_2000m <- map_df(seq_along(ras_tt_ma_2020_2000m), function(i) {
  raster <- ras_tt_ma_2020_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2020_2000m <- prox(ras_tt_ma_2020_2000m, 3)

shape_tt_ma_2020_2000m <- shape(ras_tt_ma_2020_2000m, class_value = 3)

lsm_tt_ma_2020_2000m <- met_tt_ma_2020_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2020")))|>
  inner_join(prox_tt_ma_2020_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2020_2000m, by = "id_unico")

# 2019

id_unico <- names(ras_tt_ma_2019_2000m)

met_tt_ma_2019_2000m <- map_df(seq_along(ras_tt_ma_2019_2000m), function(i) {
  raster <- ras_tt_ma_2019_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2019_2000m <- prox(ras_tt_ma_2019_2000m, 3)

shape_tt_ma_2019_2000m <- shape(ras_tt_ma_2019_2000m, class_value = 3)

lsm_tt_ma_2019_2000m <- met_tt_ma_2019_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2019")))|>
  inner_join(prox_tt_ma_2019_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2019_2000m, by = "id_unico")

# 2018

id_unico <- names(ras_tt_ma_2018_2000m)

met_tt_ma_2018_2000m <- map_df(seq_along(ras_tt_ma_2018_2000m), function(i) {
  raster <- ras_tt_ma_2018_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2018_2000m <- prox(ras_tt_ma_2018_2000m, 3)

shape_tt_ma_2018_2000m <- shape(ras_tt_ma_2018_2000m, class_value = 3)

lsm_tt_ma_2018_2000m <- met_tt_ma_2018_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2018")))|>
  inner_join(prox_tt_ma_2018_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2018_2000m, by = "id_unico")

# 2017

id_unico <- names(ras_tt_ma_2017_2000m)

met_tt_ma_2017_2000m <- map_df(seq_along(ras_tt_ma_2017_2000m), function(i) {
  raster <- ras_tt_ma_2017_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2017_2000m <- prox(ras_tt_ma_2017_2000m, 3)

shape_tt_ma_2017_2000m <- shape(ras_tt_ma_2017_2000m, class_value = 3)

lsm_tt_ma_2017_2000m <- met_tt_ma_2017_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2017")))|>
  inner_join(prox_tt_ma_2017_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2017_2000m, by = "id_unico")


# 2016

id_unico <- names(ras_tt_ma_2016_2000m)

met_tt_ma_2016_2000m <- map_df(seq_along(ras_tt_ma_2016_2000m), function(i) {
  raster <- ras_tt_ma_2016_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2016_2000m <- prox(ras_tt_ma_2016_2000m, 3)

shape_tt_ma_2016_2000m <- shape(ras_tt_ma_2016_2000m, class_value = 3)

lsm_tt_ma_2016_2000m <- met_tt_ma_2016_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2016")))|>
  inner_join(prox_tt_ma_2016_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2016_2000m, by = "id_unico")

# 2015

id_unico <- names(ras_tt_ma_2015_2000m)

met_tt_ma_2015_2000m <- map_df(seq_along(ras_tt_ma_2015_2000m), function(i) {
  raster <- ras_tt_ma_2015_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2015_2000m <- prox(ras_tt_ma_2015_2000m, 3)

shape_tt_ma_2015_2000m <- shape(ras_tt_ma_2015_2000m, class_value = 3)

lsm_tt_ma_2015_2000m <- met_tt_ma_2015_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2015")))|>
  inner_join(prox_tt_ma_2015_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2015_2000m, by = "id_unico")

# 2014

id_unico <- names(ras_tt_ma_2014_2000m)

met_tt_ma_2014_2000m <- map_df(seq_along(ras_tt_ma_2014_2000m), function(i) {
  raster <- ras_tt_ma_2014_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2014_2000m <- prox(ras_tt_ma_2014_2000m, 3)

shape_tt_ma_2014_2000m <- shape(ras_tt_ma_2014_2000m, class_value = 3)

lsm_tt_ma_2014_2000m <- met_tt_ma_2014_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2014")))|>
  inner_join(prox_tt_ma_2014_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2014_2000m, by = "id_unico")

# 2013

id_unico <- names(ras_tt_ma_2013_2000m)

met_tt_ma_2013_2000m <- map_df(seq_along(ras_tt_ma_2013_2000m), function(i) {
  raster <- ras_tt_ma_2013_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2013_2000m <- prox(ras_tt_ma_2013_2000m, 3)

shape_tt_ma_2013_2000m <- shape(ras_tt_ma_2013_2000m, class_value = 3)

lsm_tt_ma_2013_2000m <- met_tt_ma_2013_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2013")))|>
  inner_join(prox_tt_ma_2013_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2013_2000m, by = "id_unico")

# 2012
id_unico <- names(ras_tt_ma_2012_2000m)

met_tt_ma_2012_2000m <- map_df(seq_along(ras_tt_ma_2012_2000m), function(i) {
  raster <- ras_tt_ma_2012_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2012_2000m <- prox(ras_tt_ma_2012_2000m, 3)

shape_tt_ma_2012_2000m <- shape(ras_tt_ma_2012_2000m, class_value = 3)

lsm_tt_ma_2012_2000m <- met_tt_ma_2012_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2012")))|>
  inner_join(prox_tt_ma_2012_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2012_2000m, by = "id_unico")

# 2011

id_unico <- names(ras_tt_ma_2011_2000m)

met_tt_ma_2011_2000m <- map_df(seq_along(ras_tt_ma_2011_2000m), function(i) {
  raster <- ras_tt_ma_2011_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2011_2000m <- prox(ras_tt_ma_2011_2000m, 3)

shape_tt_ma_2011_2000m <- shape(ras_tt_ma_2011_2000m, class_value = 3)

lsm_tt_ma_2011_2000m <- met_tt_ma_2011_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2011")))|>
  inner_join(prox_tt_ma_2011_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2011_2000m, by = "id_unico")
# 2010

id_unico <- names(ras_tt_ma_2010_2000m)

met_tt_ma_2010_2000m <- map_df(seq_along(ras_tt_ma_2010_2000m), function(i) {
  raster <- ras_tt_ma_2010_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2010_2000m <- prox(ras_tt_ma_2010_2000m, 3)

shape_tt_ma_2010_2000m <- shape(ras_tt_ma_2010_2000m, class_value = 3)

lsm_tt_ma_2010_2000m <- met_tt_ma_2010_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2010")))|>
  inner_join(prox_tt_ma_2010_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2010_2000m, by = "id_unico")

# 2009
id_unico <- names(ras_tt_ma_2009_2000m)

met_tt_ma_2009_2000m <- map_df(seq_along(ras_tt_ma_2009_2000m), function(i) {
  raster <- ras_tt_ma_2009_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2009_2000m <- prox(ras_tt_ma_2009_2000m, 3)

shape_tt_ma_2009_2000m <- shape(ras_tt_ma_2009_2000m, class_value = 3)

lsm_tt_ma_2009_2000m <- met_tt_ma_2009_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2009")))|>
  inner_join(prox_tt_ma_2009_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2009_2000m, by = "id_unico")

# 2008

id_unico <- names(ras_tt_ma_2008_2000m)

met_tt_ma_2008_2000m <- map_df(seq_along(ras_tt_ma_2008_2000m), function(i) {
  raster <- ras_tt_ma_2008_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2008_2000m <- prox(ras_tt_ma_2008_2000m, 3)

shape_tt_ma_2008_2000m <- shape(ras_tt_ma_2008_2000m, class_value = 3)

lsm_tt_ma_2008_2000m <- met_tt_ma_2008_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2008")))|>
  inner_join(prox_tt_ma_2008_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2008_2000m, by = "id_unico")

# 2007

id_unico <- names(ras_tt_ma_2007_2000m)

met_tt_ma_2007_2000m <- map_df(seq_along(ras_tt_ma_2007_2000m), function(i) {
  raster <- ras_tt_ma_2007_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2007_2000m <- prox(ras_tt_ma_2007_2000m, 3)

shape_tt_ma_2007_2000m <- shape(ras_tt_ma_2007_2000m, class_value = 3)

lsm_tt_ma_2007_2000m <- met_tt_ma_2007_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2007")))|>
  inner_join(prox_tt_ma_2007_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2007_2000m, by = "id_unico")

# 2006

id_unico <- names(ras_tt_ma_2006_2000m)

met_tt_ma_2006_2000m <- map_df(seq_along(ras_tt_ma_2006_2000m), function(i) {
  raster <- ras_tt_ma_2006_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2006_2000m <- prox(ras_tt_ma_2006_2000m, 3)

shape_tt_ma_2006_2000m <- shape(ras_tt_ma_2006_2000m, class_value = 3)

lsm_tt_ma_2006_2000m <- met_tt_ma_2006_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2006")))|>
  inner_join(prox_tt_ma_2006_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2006_2000m, by = "id_unico")

# 2005

id_unico <- names(ras_tt_ma_2005_2000m)

met_tt_ma_2005_2000m <- map_df(seq_along(ras_tt_ma_2005_2000m), function(i) {
  raster <- ras_tt_ma_2005_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2005_2000m <- prox(ras_tt_ma_2005_2000m, 3)

shape_tt_ma_2005_2000m <- shape(ras_tt_ma_2005_2000m, class_value = 3)

lsm_tt_ma_2005_2000m <- met_tt_ma_2005_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2005")))|>
  inner_join(prox_tt_ma_2005_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2005_2000m, by = "id_unico")

# 2004

id_unico <- names(ras_tt_ma_2004_2000m)

met_tt_ma_2004_2000m <- map_df(seq_along(ras_tt_ma_2004_2000m), function(i) {
  raster <- ras_tt_ma_2004_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2004_2000m <- prox(ras_tt_ma_2004_2000m, 3)

shape_tt_ma_2004_2000m <- shape(ras_tt_ma_2004_2000m, class_value = 3)

lsm_tt_ma_2004_2000m <- met_tt_ma_2004_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2004")))|>
  inner_join(prox_tt_ma_2004_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2004_2000m, by = "id_unico")

# 2003

id_unico <- names(ras_tt_ma_2003_2000m)

met_tt_ma_2003_2000m <- map_df(seq_along(ras_tt_ma_2003_2000m), function(i) {
  raster <- ras_tt_ma_2003_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2003_2000m <- prox(ras_tt_ma_2003_2000m, 3)

shape_tt_ma_2003_2000m <- shape(ras_tt_ma_2003_2000m, class_value = 3)

lsm_tt_ma_2003_2000m <- met_tt_ma_2003_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2003")))|>
  inner_join(prox_tt_ma_2003_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2003_2000m, by = "id_unico")

# 2002

id_unico <- names(ras_tt_ma_2002_2000m)

met_tt_ma_2002_2000m <- map_df(seq_along(ras_tt_ma_2002_2000m), function(i) {
  raster <- ras_tt_ma_2002_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2002_2000m <- prox(ras_tt_ma_2002_2000m, 3)

shape_tt_ma_2002_2000m <- shape(ras_tt_ma_2002_2000m, class_value = 3)

lsm_tt_ma_2002_2000m <- met_tt_ma_2002_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2002")))|>
  inner_join(prox_tt_ma_2002_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2002_2000m, by = "id_unico")

# 2001

id_unico <- names(ras_tt_ma_2001_2000m)

met_tt_ma_2001_2000m <- map_df(seq_along(ras_tt_ma_2001_2000m), function(i) {
  raster <- ras_tt_ma_2001_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2001_2000m <- prox(ras_tt_ma_2001_2000m, 3)

shape_tt_ma_2001_2000m <- shape(ras_tt_ma_2001_2000m, class_value = 3)

lsm_tt_ma_2001_2000m <- met_tt_ma_2001_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2001")))|>
  inner_join(prox_tt_ma_2001_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2001_2000m, by = "id_unico")

# 2000

id_unico <- names(ras_tt_ma_2000_2000m)

met_tt_ma_2000_2000m <- map_df(seq_along(ras_tt_ma_2000_2000m), function(i) {
  raster <- ras_tt_ma_2000_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_2000_2000m <- prox(ras_tt_ma_2000_2000m, 3)

shape_tt_ma_2000_2000m <- shape(ras_tt_ma_2000_2000m, class_value = 3)

lsm_tt_ma_2000_2000m <- met_tt_ma_2000_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("2000")))|>
  inner_join(prox_tt_ma_2000_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_2000_2000m, by = "id_unico")



# 1999

id_unico <- names(ras_tt_ma_1999_2000m)

met_tt_ma_1999_2000m <- map_df(seq_along(ras_tt_ma_1999_2000m), function(i) {
  raster <- ras_tt_ma_1999_2000m[[i]]
  id_unico <- id_unico[i]
  
  df <- calculate_lsm(
    raster, 
    what = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_ed", "lsm_c_pd", "lsm_c_np"),
    directions = 8,
    count_boundary = FALSE, 
    consider_boundary = TRUE,
    neighbourhood = 8,
    progress = TRUE)
  
  df$id_unico <- id_unico # Adicionar coluna com o nome do raster
  return(df)})

prox_tt_ma_1999_2000m <- prox(ras_tt_ma_1999_2000m, 3)

shape_tt_ma_1999_2000m <- shape(ras_tt_ma_1999_2000m, class_value = 3)

lsm_tt_ma_1999_2000m <- met_tt_ma_1999_2000m |>
  unite("metric_class", metric, class, sep = "_") |>
  select(-level, -id) |>
  pivot_wider(
    names_from = metric_class,
    values_from = value) |>
  group_by(id_unico) |>  
  dplyr::summarise(across(everything(), ~ first(na.omit(.)))) |>
  dplyr::rename(Ano = layer) |>
  mutate(Ano = factor(Ano,
                      levels = c("1"),
                      labels = c("1999")))|>
  inner_join(prox_tt_ma_1999_2000m, by = "id_unico")|>
  inner_join(shape_tt_ma_1999_2000m, by = "id_unico")

#### Unificando tabelas 2000m ####

lsm_tt_ma_2000m <- bind_rows(lsm_tt_ma_2023_2000m, lsm_tt_ma_2022_2000m, lsm_tt_ma_2021_2000m, lsm_tt_ma_2020_2000m, lsm_tt_ma_2019_2000m, lsm_tt_ma_2018_2000m, lsm_tt_ma_2017_2000m, lsm_tt_ma_2016_2000m, lsm_tt_ma_2015_2000m, lsm_tt_ma_2014_2000m, lsm_tt_ma_2013_2000m, lsm_tt_ma_2012_2000m, lsm_tt_ma_2011_2000m, lsm_tt_ma_2010_2000m, lsm_tt_ma_2009_2000m, lsm_tt_ma_2008_2000m, lsm_tt_ma_2007_2000m, lsm_tt_ma_2006_2000m, lsm_tt_ma_2005_2000m, lsm_tt_ma_2004_2000m, lsm_tt_ma_2003_2000m, lsm_tt_ma_2002_2000m, lsm_tt_ma_2001_2000m, lsm_tt_ma_2000_2000m, lsm_tt_ma_1999_2000m) |>
  dplyr::select(-pland_0,-lpi_0, -ed_0, -pd_0, -np_0, -lpi_12, -ed_12, -pd_12, -np_12, -lpi_15, -ed_15, -pd_15, -np_15, -lpi_9, -ed_9, -pd_9, -np_9, -np_24, -pd_24, -lpi_24, -ed_24,  -pland_24)|>
  mutate(Bin = str_extract(id_unico, "^[01]"),
         Bin = as.factor(Bin),     
         Ano = as.factor(Ano))|>
  mutate(pland_9 = ifelse(is.na(pland_9), 0, pland_9))|>
  mutate(pland_15 = ifelse(is.na(pland_15), 0, pland_15))|>
  mutate(pland_12 = ifelse(is.na(pland_12), 0, pland_12)) |>
  rename_with(~ paste0(., "_2000m"))


lsm_tt_ma_2000m_sem_na <- lsm_tt_ma_2000m|>
  na.omit() 

write_xlsx(lsm_tt_ma_2000m_sem_na, "lsm_tt_ma_2000m.xlsx")

View(lsm_tt_ma_2000m_sem_na)



#### Distancia massa dagua ####

dist_agua_tt_ma <- st_distance(tt_ma_flo, agua_ma)

tt_ma_flo$dist_min_agua <- as.numeric(apply(dist_agua_tt_ma, 1, min))

dist_min_agua_ma_tt <- st_drop_geometry(tt_ma_flo[, c("id_unico", "dist_min_agua")])



#### Distancia UCs ####

UC$cria_ano <- as.numeric(format(as.Date(UC$cria_ano, format = "%d-%m-%Y"), "%Y"))

# Método por iteração nos pontos
calcular_dist_minima <- function(pontos_sf, polygons_sf) {
  resultados <- map_dfr(1:nrow(pontos_sf), function(i) {
    ponto <- pontos_sf[i, ]
    
    # Filtrar UCs com ano <= ano do ponto
    UCs_validas <- polygons_sf %>% 
      filter(cria_ano <= ponto$Ano)
    
    if (nrow(UCs_validas) == 0) {
      return(tibble(
        id_unico = ponto$id_unico,
        dist_min_UC = NA_real_,
        uc_id = NA_character_,
        nome_uc = NA_character_,
        cria_ano = NA_real_
      ))
    }
    
    # Calcular distâncias para UCs válidas
    distancias <- st_distance(ponto, UCs_validas)
    min_index <- which.min(distancias)
    
    tibble(
      id_unico = ponto$id_unico,
      ano_ocorrencia = ponto$Ano,
      dist_min_UC = as.numeric(distancias[min_index]),
      uc_id = UCs_validas$uc_id[min_index],
      nome_uc = UCs_validas$nome_uc[min_index],
      cria_ano = UCs_validas$cria_ano[min_index]
    )
  })
  
  return(resultados)
}

# Calcular distâncias mínimas
distancias_result <- calcular_dist_minima(tt_ma_flo, UC)

# Juntar com dados originais
tt_ma_flo <- tt_ma_flo %>%
  left_join(distancias_result, by = "id_unico")

#### Selecionando rodovias com pavimentação ####

table(rod$revestimen)
table(rod$revestim_1)

rod <- as_tibble(rod)

rod_pav <- rod  |> 
  dplyr::filter(revestimen == "Pavimentado" | revestim_1 == "Pavimentado")

rod_pav <- st_as_sf(rod_pav)

#### Densidade de rodovias 500m ####

buffer_500m <- rbind(bf_tt_ma_2023_500m, bf_tt_ma_2022_500m, bf_tt_ma_2021_500m, bf_tt_ma_2020_500m, bf_tt_ma_2019_500m, bf_tt_ma_2018_500m, bf_tt_ma_2017_500m, bf_tt_ma_2016_500m, bf_tt_ma_2015_500m, bf_tt_ma_2014_500m, bf_tt_ma_2013_500m, bf_tt_ma_2012_500m, bf_tt_ma_2011_500m, bf_tt_ma_2010_500m, bf_tt_ma_2009_500m, bf_tt_ma_2008_500m, bf_tt_ma_2007_500m, bf_tt_ma_2006_500m, bf_tt_ma_2005_500m,
bf_tt_ma_2004_500m, bf_tt_ma_2003_500m, bf_tt_ma_2002_500m, bf_tt_ma_2001_500m, bf_tt_ma_2000_500m, bf_tt_ma_1999_500m)

rod_buffer_500m <- st_intersection(rod_pav, buffer_500m)

rod_dis_500m <- aggregate(
  x = rod_buffer_500m ["geometry"],
  by = list(id_unico = rod_buffer_500m$id_unico),
  FUN = function(x) st_union(x)) 

comp_rod_500m <- rod_dis_500m |>
  mutate(comp_rod_500m = st_length(geometry))

# Juntar as informações de comprimento com a área do buffer
den_rod_500m <- buffer_500m |> 
  left_join(comp_rod_500m |> 
              st_drop_geometry() |> 
              mutate(comprimento_rod_500m = as.numeric(comp_rod_500m)),
            by = "id_unico")|>
  mutate(comprimento_rod_500m = ifelse(is.na(comprimento_rod_500m), 0, comprimento_rod_500m),
         den_rod_500m = (comprimento_rod_500m/1000) / ((3.14*500^2)/1e6))


den_500m <- data.frame(id_unico = lsm_tt_ma_500m_total$id_unico_500m) |>
  left_join(den_rod_500m, by = "id_unico") |>
  mutate(den_rod_500m = ifelse(is.na(den_rod_500m), 0, as.numeric(den_rod_500m)))

lsm_tt_ma_500m_total$den_rod_500m <- den_500m$den_rod_500m

den_rod_2023 <- st_drop_geometry(lsm_tt_ma_500m_total[, c("id_unico_500m", "den_rod_500m")])

#### Densidade de rodovias 1000m ####

buffer_1000m <- rbind(bf_tt_ma_2023_1000m, bf_tt_ma_2022_1000m, bf_tt_ma_2021_1000m, bf_tt_ma_2020_1000m, bf_tt_ma_2019_1000m, bf_tt_ma_2018_1000m, bf_tt_ma_2017_1000m, bf_tt_ma_2016_1000m, bf_tt_ma_2015_1000m, bf_tt_ma_2014_1000m, bf_tt_ma_2013_1000m, bf_tt_ma_2012_1000m, bf_tt_ma_2011_1000m, bf_tt_ma_2010_1000m, bf_tt_ma_2009_1000m, bf_tt_ma_2008_1000m, bf_tt_ma_2007_1000m, bf_tt_ma_2006_1000m, bf_tt_ma_2005_1000m,
                     bf_tt_ma_2004_1000m, bf_tt_ma_2003_1000m, bf_tt_ma_2002_1000m, bf_tt_ma_2001_1000m, bf_tt_ma_2000_1000m, bf_tt_ma_1999_1000m)

rod_buffer_1000m <- st_intersection(rod_pav, buffer_1000m)

rod_dis_1000m <- aggregate(
  x = rod_buffer_1000m ["geometry"],
  by = list(id_unico = rod_buffer_1000m$id_unico),
  FUN = function(x) st_union(x)) 

comp_rod_1000m <- rod_dis_1000m |>
  mutate(comp_rod_1000m = st_length(geometry))

# Juntar as informações de comprimento com a área do buffer
den_rod_1000m <- buffer_1000m |> 
  left_join(comp_rod_1000m |> 
              st_drop_geometry() |> 
              mutate(comprimento_rod_1000m = as.numeric(comp_rod_1000m)),
            by = "id_unico")|>
  mutate(comprimento_rod_1000m = ifelse(is.na(comprimento_rod_1000m), 0, comprimento_rod_1000m),
         den_rod_1000m = (comprimento_rod_1000m/1000) / ((3.14*1000^2)/1e6))


den_1000m <- data.frame(id_unico = lsm_tt_ma_1000m$id_unico_1000m) |>
  left_join(den_rod_1000m, by = "id_unico") |>
  mutate(den_rod_1000m = ifelse(is.na(den_rod_1000m), 0, as.numeric(den_rod_1000m)))

lsm_tt_ma_1000m$den_rod_1000m <- den_1000m$den_rod_1000m

den_rod_2023 <- st_drop_geometry(lsm_tt_ma_1000m[, c("id_unico_1000m", "den_rod_1000m")])

#### Densidade de rodovias 2000m ####

buffer_2000m <- rbind(bf_tt_ma_2023_2000m, bf_tt_ma_2022_2000m, bf_tt_ma_2021_2000m, bf_tt_ma_2020_2000m, bf_tt_ma_2019_2000m, bf_tt_ma_2018_2000m, bf_tt_ma_2017_2000m, bf_tt_ma_2016_2000m, bf_tt_ma_2015_2000m, bf_tt_ma_2014_2000m, bf_tt_ma_2013_2000m, bf_tt_ma_2012_2000m, bf_tt_ma_2011_2000m, bf_tt_ma_2010_2000m, bf_tt_ma_2009_2000m, bf_tt_ma_2008_2000m, bf_tt_ma_2007_2000m, bf_tt_ma_2006_2000m, bf_tt_ma_2005_2000m,
                     bf_tt_ma_2004_2000m, bf_tt_ma_2003_2000m, bf_tt_ma_2002_2000m, bf_tt_ma_2001_2000m, bf_tt_ma_2000_2000m, bf_tt_ma_1999_2000m)

rod_buffer_2000m <- st_intersection(rod_pav, buffer_2000m)

rod_dis_2000m <- aggregate(
  x = rod_buffer_2000m ["geometry"],
  by = list(id_unico = rod_buffer_2000m$id_unico),
  FUN = function(x) st_union(x)) 

comp_rod_2000m <- rod_dis_2000m |>
  mutate(comp_rod_2000m = st_length(geometry))

# Juntar as informações de comprimento com a área do buffer
den_rod_2000m <- buffer_2000m |> 
  left_join(comp_rod_2000m |> 
              st_drop_geometry() |> 
              mutate(comprimento_rod_2000m = as.numeric(comp_rod_2000m)),
            by = "id_unico")|>
  mutate(comprimento_rod_2000m = ifelse(is.na(comprimento_rod_2000m), 0, comprimento_rod_2000m),
         den_rod_2000m = (comprimento_rod_2000m/1000) / ((3.14*2000^2)/1e6))


den_2000m <- data.frame(id_unico = lsm_tt_ma_2000m$id_unico_2000m) |>
  left_join(den_rod_2000m, by = "id_unico") |>
  mutate(den_rod_2000m = ifelse(is.na(den_rod_2000m), 0, as.numeric(den_rod_2000m)))

lsm_tt_ma_2000m$den_rod_2000m <- den_2000m$den_rod_2000m

den_rod_2023 <- st_drop_geometry(lsm_tt_ma_2000m[, c("id_unico_2000m", "den_rod_2000m")])
#### Unificando todas as tabelas ####

lsm_tt_ma_500m_uni <- lsm_tt_ma_500m_total |>
  dplyr::rename('Bin'='Bin_500m')|>
  dplyr::rename('Ano'='Ano_500m')|>
  dplyr::rename('id_unico'='id_unico_500m')

lsm_tt_ma_1000m_uni <- lsm_tt_ma_1000m |>
  dplyr::rename('id_unico'='id_unico_1000m')|>
  select(-Bin_1000m, -Ano_1000m)

lsm_tt_ma_2000m_uni <- lsm_tt_ma_2000m |>
  dplyr::rename('id_unico'='id_unico_2000m')|>
  select(-Bin_2000m, -Ano_2000m)

glimpse(c(lsm_tt_ma_500m_uni, lsm_tt_ma_1000m_uni, lsm_tt_ma_2000m_uni))

dist_min_uc_tt_ma <- select(distancias_result, c("id_unico", "dist_min_UC"))

lsm_tt_ma_total <- lsm_tt_ma_500m_uni %>%
  full_join(lsm_tt_ma_1000m_uni, by = 'id_unico') %>%
  full_join(lsm_tt_ma_2000m_uni, by = 'id_unico') %>%
  full_join(dist_min_agua_ma_tt, by = 'id_unico') %>%
  full_join(dist_min_uc_tt_ma, by = 'id_unico')

glimpse(lsm_tt_ma_total)

write_xlsx(lsm_tt_ma_total, "lsm_tt_ma_total.xlsx")

#### Estandartizar ####

# 1. Separar as colunas de identificação e as colunas numéricas
id_cols <- lsm_tt_ma_total |> select(id_unico, Ano, Bin)

num_cols <- lsm_tt_ma_total |> select(-id_unico, -Ano, -Bin)


# 2. Aplicar a padronização (exemplo com método "standardize" - scale para média 0 e sd 1)
stan <- vegan::decostand(num_cols, method = "standardize")

# 3. Juntar novamente com as colunas de identificação
lsm_tt_ma_total_std <- bind_cols(id_cols, stan)

# Verificar o resultado
glimpse(lsm_tt_ma_total_std)

write_xlsx(lsm_tt_ma_total_std, "lsm_tt_ma_total_std.xlsx")

lsm_tt_ma_total_std$Bin <- as.factor(lsm_tt_ma_total_std$Bin)


var_tt_ma <- lsm_tt_ma_total_std 

#### Correlacao ####

num_var_tt_ma <- var_tt_ma|>
  dplyr::select(-id_unico, - Bin, -Ano)

cor_tt_ma <- cor(num_var_tt_ma, method = "spearman")

glimpse(cor_tt_ma)


# Remove combinacoes duplicadas e mesmas metricas com escalas diferentes
cor_sem_dup_tt_ma <- as.data.frame(cor_tt_ma) %>%
  mutate(var1 = rownames(cor_tt_ma)) %>%
  tidyr::pivot_longer(
    cols = -var1,
    names_to = "var2",
    values_to = "correlacao"
  ) %>%
  # Filtrar correlações acima de 0.7 (e remover autocorrelações)
  filter(
    correlacao <= 0.7,  
    var1 != var2
  ) %>%
  # Criar uma chave ordenada de forma mais segura
  mutate(
    # Criar vetor ordenado sem usar rowwise()
    key = purrr::map2_chr(var1, var2, ~paste(sort(c(.x, .y)), collapse = "|"))
  ) %>%
  # Remover duplicatas
  distinct(key, .keep_all = TRUE) %>%
  # Remover coluna auxiliar
  select(-key) %>%
  mutate(
    # Remove a escala (última parte após o último underscore)
    base_var1 = str_remove(var1, "_\\d+m$"),
    base_var2 = str_remove(var2, "_\\d+m$")
  ) %>%
  # Remover se a métrica base for a mesma
  filter(base_var1 != base_var2) %>%
  # Remover colunas auxiliares
  select(-base_var1, -base_var2) %>%
  arrange(desc(correlacao))


cor_sem_dup_tt_ma

View(cor_sem_dup_tt_ma)

#### Modelos ####

# Extrair combinações únicas para modelos
# Extrair combinações únicas para modelos
var_glm <- function(pares_df) {
  
  # 1. Variáveis para modelos individuais (apenas da coluna var1)
  var_individual <- unique(c(pares_df$var1, pares_df$var2))
  
  # 2. Combinações para modelos com pares
  combinacoes_pares <- pares_df %>%
    select(var1, var2) %>%
    # Garantir que não haja duplicatas reversas
    mutate(
      chave = mapply(function(x, y) paste(sort(c(x, y)), collapse = "|"), 
                     var1, var2)
    ) %>%
    distinct(chave, .keep_all = TRUE) %>%
    select(-chave)
  
  # 3. Combinações para modelos com interecoes
  combinacoes_interacoes <- pares_df %>%
    select(var1, var2) %>%
    # Garantir que não haja duplicatas reversas
    mutate(
      chave = mapply(function(x, y) paste(sort(c(x, y)), collapse = "|"), 
                     var1, var2)
    ) %>%
    distinct(chave, .keep_all = TRUE) %>%
    select(-chave)
  
  return(list(
    modelos_individuais = var_individual,
    modelos_pares = combinacoes_pares,
    modelos_interecoes = combinacoes_interacoes
  ))
}

# Extrair combinações
comb <- var_glm(cor_sem_dup_tt_ma)

# Ver resultados
cat("Variáveis para modelos individuais (", length(comb$modelos_individuais), "):\n")
print(comb$modelos_individuais)

cat("\n\nPares para modelos combinados (", nrow(comb$modelos_pares), "):\n")
print(comb$modelos_pares, n = 20)

cat("\n\nInterações para modelos combinados (", nrow(comb$modelos_interecoes), "):\n")
print(comb$modelos_interecoes, n = 20)

# Função para criar fórmulas
criar_formulas <- function(var_individuais, pares, resposta = "Bin") {
  
  formulas <- list()
  
  # Fórmulas individuais
  for (var in var_individuais) {
    nome <- paste("ind", var, sep = "_")
    formulas[[nome]] <- as.formula(paste(resposta, "~", var))
  }
  
  # Fórmulas com pares
  for (i in 1:nrow(pares)) {
    var1 <- pares$var1[i]
    var2 <- pares$var2[i]
    nome <- paste("par", var1, var2, sep = "_")
    formulas[[nome]] <- as.formula(paste(resposta, "~", var1, "+", var2))
  }
  
  
  # Fórmulas com interações
  for (i in 1:nrow(pares)) {
    var1 <- pares$var1[i]
    var2 <- pares$var2[i]
    nome <- paste("int", var1, var2, sep = "_")
    formulas[[nome]] <- as.formula(paste(resposta, "~", var1, "*", var2))
  }
  
  return(formulas)
}

# Criar todas as fórmulas
formulas_glm <- criar_formulas(
  comb$modelos_individuais,
  comb$modelos_pares,
  resposta = "Bin"  # Substitua pelo nome da sua variável resposta
)

# Ver algumas fórmulas
cat("\n=== EXEMPLOS DE FÓRMULAS ===\n")
for (i in 1:min(5, length(formulas_glm))) {
  cat(names(formulas_glm)[i], ": ", 
      format(formulas_glm[[i]]), "\n", sep = "")
}

# Agora você pode iterar pelas fórmulas para criar os modelos
# Exemplo:
modelos <- list()
for (nome in names(formulas_glm)) {
  modelos[[nome]] <- glm(formulas_glm[[nome]], 
                         data = var_tt_ma, 
                         family = binomial)
}
modelos


# Selecionar modelos por VIF>4

preparar_mumin <- function(lista_modelos, max_vif = 4) {
  
  # Função auxiliar para verificar VIF
  verifica_vif <- function(modelo) {
    tryCatch({
      vif_val <- car::vif(modelo)
      if (is.matrix(vif_val)) vif_val <- vif_val[, 1]
      max(vif_val) <= max_vif
    }, error = function(e) FALSE)
  }
  
  # Filtrar: individuais (todos) + múltiplos (VIF ok)
  modelos_filtrados <- lista_modelos[
    grepl("^ind_", names(lista_modelos)) | 
      sapply(lista_modelos, verifica_vif)
  ]
  
  cat(length(lista_modelos), "→", length(modelos_filtrados), "modelos (VIF ≤", max_vif, ")\n")
  
  return(modelos_filtrados)
}

modelos_vif<-preparar_mumin(modelos)

options(max.print = 10000)

aic_tabela <- MuMIn::model.sel(modelos_vif, rank = AIC)
aic_tabela



melhor_modelo_1 <- MuMIn::get.models(aic_tabela, subset = 1)[[1]]
summary(melhor_modelo_1)
simulationOutput <- DHARMa::simulateResiduals(fittedModel = melhor_modelo_1, plot = TRUE)


melhor_modelo_2 <- MuMIn::get.models(aic_tabela, subset = 1)[[1]]
summary(melhor_modelo_2)
simulationOutput <- DHARMa::simulateResiduals(fittedModel = melhor_modelo_2, plot = TRUE)

