#### Pacotes ####  

library(readxl)
library(sf)
library(terra)
library(landscapemetrics)
library(tidyverse)
library(writexl)
library(vegan)
library(glmmTMB)
library(DHARMa)
library(Rmisc)
library(pROC)

#### Carregando dados ####

mt_ma <- read_excel("D:/CETESB/Rodovias e pontos de ausencia/bin_mt_ma_2010_2023.xlsx", sheet = 1,na = c(" ", "NA")) |> 
  st_as_sf(coords = c("coord_x", "coord_y"), crs=4674)|>
  st_transform(5880)
agua_ma <- st_read("D:/Shapefiles/agua_ma_fei.shp")|>
  st_make_valid()|>
  st_transform(5880) 
rod <- st_read("D:/CETESB/Rodovias e pontos de ausencia/rod_final.shp")|>
  st_transform(5880) 
UC <- st_read("D:/Shapefiles/UCs 2025/cnuc_2025_03.shp")|>
  st_make_valid() |>
  st_transform(5880) 

ma_2023 <- rast("D:/Analise estatistica mestrado/proj_ma_2023.tif") # EPSG 5880
ma_2022 <- rast("D:/Analise estatistica mestrado/proj_ma_2022.tif")
ma_2021 <- rast("D:/Analise estatistica mestrado/proj_ma_2021.tif")
ma_2020 <- rast("D:/Analise estatistica mestrado/proj_ma_2020.tif")
ma_2019 <- rast("D:/Analise estatistica mestrado/proj_ma_2019.tif")
ma_2018 <- rast("D:/Analise estatistica mestrado/proj_ma_2018.tif")
ma_2017 <- rast("D:/Analise estatistica mestrado/proj_ma_2017.tif")
ma_2016 <- rast("D:/Analise estatistica mestrado/proj_ma_2016.tif")
ma_2015 <- rast("D:/Analise estatistica mestrado/proj_ma_2015.tif")
ma_2014 <- rast("D:/Analise estatistica mestrado/proj_ma_2014.tif")
ma_2013 <- rast("D:/Analise estatistica mestrado/proj_ma_2013.tif")
ma_2012 <- rast("D:/Analise estatistica mestrado/proj_ma_2012.tif")
ma_2011 <- rast("D:/Analise estatistica mestrado/proj_ma_2011.tif")
ma_2010 <- rast("D:/Analise estatistica mestrado/proj_ma_2010.tif")

#### Verificando duplicatas de pseudo-ausencia ####

dup_mt_ma <- mt_ma |>
  group_by(geometry) |>
  filter(n() > 1) 
dup_mt_ma


mt_ma <- mt_ma|>
  select(-Mes)|>
  dplyr::rename(Bin=id)

#### Inserindo ID nas ocorrencias ####

mt_ma$id_unico <- ave(
  mt_ma$Bin,             
  mt_ma$Bin,             
  FUN = function(x) paste0(x, ".", seq_along(x)))
mt_ma

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
      result <- calculate_prox_single_cell_centers(raster_obj, target_class, directions)
      
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

calculate_prox_single_cell_centers <- function(landscape, target_class = 3, directions = 8) {
  
  # Get raster resolution with 4 decimal places
  res_val <- round(terra::res(landscape)[1], 4)
  
  # Create binary raster for target class
  class_raster <- landscape == target_class
  class_raster[class_raster == 0] <- NA
  
  # Check if class exists
  class_cells <- terra::global(class_raster, "sum", na.rm = TRUE)[[1]]
  if (is.na(class_cells) || class_cells == 0) {
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
    area_cells <- terra::global(patch_mask, "sum", na.rm = TRUE)[[1]]
    
    # Get coordinates of all cells in this patch
    cell_coords <- terra::xyFromCell(patch_mask, which(terra::values(patch_mask) == 1))
    
    patch_data[[as.character(patch_id)]] <- list(
      area_m2 = round(area_cells * (res_val^2), 4),
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
        min_distance <- calculate_min_distance_cell_centers(focal_coords, neighbor_coords)
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
calculate_min_distance_cell_centers <- function(coords1, coords2) {
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

bf_mt_ma_2023_500m <- mt_ma |>
  filter(Ano == "2023") |>
  st_buffer(dist = 500)

bf_mt_ma_2022_500m <- mt_ma |>
  filter(Ano == "2022") |>
  st_buffer(dist = 500)

bf_mt_ma_2021_500m <- mt_ma |>
  filter(Ano == "2021") |>
  st_buffer(dist = 500)

bf_mt_ma_2020_500m <- mt_ma |>
  filter(Ano == "2020") |>
  st_buffer(dist = 500)

bf_mt_ma_2019_500m <- mt_ma |>
  filter(Ano == "2019") |>
  st_buffer(dist = 500)

bf_mt_ma_2018_500m <- mt_ma |>
  filter(Ano == "2018") |>
  st_buffer(dist = 500)

bf_mt_ma_2017_500m <- mt_ma |>
  filter(Ano == "2017") |>
  st_buffer(dist = 500)

bf_mt_ma_2016_500m <- mt_ma |>
  filter(Ano == "2016") |>
  st_buffer(dist = 500)

bf_mt_ma_2015_500m <- mt_ma |>
  filter(Ano == "2015") |>
  st_buffer(dist = 500)

bf_mt_ma_2014_500m <- mt_ma |>
  filter(Ano == "2014") |>
  st_buffer(dist = 500)

bf_mt_ma_2013_500m <- mt_ma |>
  filter(Ano == "2013") |>
  st_buffer(dist = 500)

bf_mt_ma_2012_500m <- mt_ma |>
  filter(Ano == "2012") |>
  st_buffer(dist = 500)

bf_mt_ma_2011_500m <- mt_ma |>
  filter(Ano == "2011") |>
  st_buffer(dist = 500)

bf_mt_ma_2010_500m <- mt_ma |>
  filter(Ano == "2010") |>
  st_buffer(dist = 500)

#### Cortando raster 500m ####

# 2023

ras_mt_ma_2023_500m <- list()

for (i in 1:nrow(bf_mt_ma_2023_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2023_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2023_500m <- crop(ma_2023, bf_mt_ma_2023_500m[i, ])
  mask_mt_ma_2023_500m <- mask(crop_mt_ma_2023_500m, bf_mt_ma_2023_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2023_500m[[buffer_id]] <- mask_mt_ma_2023_500m
}

output_dir <- "ras_mt_ma_2023_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2023_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2023_500m.tif"))
  writeRaster(
    ras_mt_ma_2023_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2022

ras_mt_ma_2022_500m <- list()

for (i in 1:nrow(bf_mt_ma_2022_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2022_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2022_500m <- crop(ma_2022, bf_mt_ma_2022_500m[i, ])
  mask_mt_ma_2022_500m <- mask(crop_mt_ma_2022_500m, bf_mt_ma_2022_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2022_500m[[buffer_id]] <- mask_mt_ma_2022_500m
}

output_dir <- "ras_mt_ma_2022_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2022_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2022_500m.tif"))
  writeRaster(
    ras_mt_ma_2022_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2021

ras_mt_ma_2021_500m <- list()

for (i in 1:nrow(bf_mt_ma_2021_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2021_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2021_500m <- crop(ma_2021, bf_mt_ma_2021_500m[i, ])
  mask_mt_ma_2021_500m <- mask(crop_mt_ma_2021_500m, bf_mt_ma_2021_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2021_500m[[buffer_id]] <- mask_mt_ma_2021_500m
}

output_dir <- "ras_mt_ma_2021_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2021_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2021_500m.tif"))
  writeRaster(
    ras_mt_ma_2021_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2020

ras_mt_ma_2020_500m <- list()

for (i in 1:nrow(bf_mt_ma_2020_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2020_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2020_500m <- crop(ma_2020, bf_mt_ma_2020_500m[i, ])
  mask_mt_ma_2020_500m <- mask(crop_mt_ma_2020_500m, bf_mt_ma_2020_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2020_500m[[buffer_id]] <- mask_mt_ma_2020_500m
}

output_dir <- "ras_mt_ma_2020_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2020_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2020_500m.tif"))
  writeRaster(
    ras_mt_ma_2020_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2019

ras_mt_ma_2019_500m <- list()

for (i in 1:nrow(bf_mt_ma_2019_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2019_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2019_500m <- crop(ma_2019, bf_mt_ma_2019_500m[i, ])
  mask_mt_ma_2019_500m <- mask(crop_mt_ma_2019_500m, bf_mt_ma_2019_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2019_500m[[buffer_id]] <- mask_mt_ma_2019_500m
}

output_dir <- "ras_mt_ma_2019_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2019_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2019_500m.tif"))
  writeRaster(
    ras_mt_ma_2019_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2018

ras_mt_ma_2018_500m <- list()

for (i in 1:nrow(bf_mt_ma_2018_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2018_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2018_500m <- crop(ma_2018, bf_mt_ma_2018_500m[i, ])
  mask_mt_ma_2018_500m <- mask(crop_mt_ma_2018_500m, bf_mt_ma_2018_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2018_500m[[buffer_id]] <- mask_mt_ma_2018_500m
}

output_dir <- "ras_mt_ma_2018_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2018_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2018_500m.tif"))
  writeRaster(
    ras_mt_ma_2018_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2017

ras_mt_ma_2017_500m <- list()

for (i in 1:nrow(bf_mt_ma_2017_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2017_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2017_500m <- crop(ma_2017, bf_mt_ma_2017_500m[i, ])
  mask_mt_ma_2017_500m <- mask(crop_mt_ma_2017_500m, bf_mt_ma_2017_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2017_500m[[buffer_id]] <- mask_mt_ma_2017_500m
}

output_dir <- "ras_mt_ma_2017_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2017_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2017_500m.tif"))
  writeRaster(
    ras_mt_ma_2017_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2016

ras_mt_ma_2016_500m <- list()

for (i in 1:nrow(bf_mt_ma_2016_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2016_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2016_500m <- crop(ma_2016, bf_mt_ma_2016_500m[i, ])
  mask_mt_ma_2016_500m <- mask(crop_mt_ma_2016_500m, bf_mt_ma_2016_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2016_500m[[buffer_id]] <- mask_mt_ma_2016_500m
}

output_dir <- "ras_mt_ma_2016_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2016_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2016_500m.tif"))
  writeRaster(
    ras_mt_ma_2016_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2015

ras_mt_ma_2015_500m <- list()

for (i in 1:nrow(bf_mt_ma_2015_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2015_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2015_500m <- crop(ma_2015, bf_mt_ma_2015_500m[i, ])
  mask_mt_ma_2015_500m <- mask(crop_mt_ma_2015_500m, bf_mt_ma_2015_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2015_500m[[buffer_id]] <- mask_mt_ma_2015_500m
}

output_dir <- "ras_mt_ma_2015_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2015_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2015_500m.tif"))
  writeRaster(
    ras_mt_ma_2015_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2014

ras_mt_ma_2014_500m <- list()

for (i in 1:nrow(bf_mt_ma_2014_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2014_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2014_500m <- crop(ma_2014, bf_mt_ma_2014_500m[i, ])
  mask_mt_ma_2014_500m <- mask(crop_mt_ma_2014_500m, bf_mt_ma_2014_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2014_500m[[buffer_id]] <- mask_mt_ma_2014_500m
}

output_dir <- "ras_mt_ma_2014_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2014_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2014_500m.tif"))
  writeRaster(
    ras_mt_ma_2014_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2013

ras_mt_ma_2013_500m <- list()

for (i in 1:nrow(bf_mt_ma_2013_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2013_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2013_500m <- crop(ma_2013, bf_mt_ma_2013_500m[i, ])
  mask_mt_ma_2013_500m <- mask(crop_mt_ma_2013_500m, bf_mt_ma_2013_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2013_500m[[buffer_id]] <- mask_mt_ma_2013_500m
}

output_dir <- "ras_mt_ma_2013_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2013_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2013_500m.tif"))
  writeRaster(
    ras_mt_ma_2013_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2012

ras_mt_ma_2012_500m <- list()

for (i in 1:nrow(bf_mt_ma_2012_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2012_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2012_500m <- crop(ma_2012, bf_mt_ma_2012_500m[i, ])
  mask_mt_ma_2012_500m <- mask(crop_mt_ma_2012_500m, bf_mt_ma_2012_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2012_500m[[buffer_id]] <- mask_mt_ma_2012_500m
}

output_dir <- "ras_mt_ma_2012_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2012_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2012_500m.tif"))
  writeRaster(
    ras_mt_ma_2012_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2011

ras_mt_ma_2011_500m <- list()

for (i in 1:nrow(bf_mt_ma_2011_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2011_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2011_500m <- crop(ma_2011, bf_mt_ma_2011_500m[i, ])
  mask_mt_ma_2011_500m <- mask(crop_mt_ma_2011_500m, bf_mt_ma_2011_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2011_500m[[buffer_id]] <- mask_mt_ma_2011_500m
}

output_dir <- "ras_mt_ma_2011_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2011_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2011_500m.tif"))
  writeRaster(
    ras_mt_ma_2011_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2010

ras_mt_ma_2010_500m <- list()

for (i in 1:nrow(bf_mt_ma_2010_500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2010_500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2010_500m <- crop(ma_2010, bf_mt_ma_2010_500m[i, ])
  mask_mt_ma_2010_500m <- mask(crop_mt_ma_2010_500m, bf_mt_ma_2010_500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2010_500m[[buffer_id]] <- mask_mt_ma_2010_500m
}

output_dir <- "ras_mt_ma_2010_500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2010_500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2010_500m.tif"))
  writeRaster(
    ras_mt_ma_2010_500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}




# Chamando os recortes

output_dir <- "ras_mt_ma_2023_500m"
ras_mt_ma_2023_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) |>
  set_names(~ str_remove(basename(.), "ras_mt_ma_2023_500m\\.tif$")) |>
  map(rast) 

output_dir <- "ras_mt_ma_2022_500m"
ras_mt_ma_2022_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2022_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2021_500m"
ras_mt_ma_2021_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2021_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2020_500m"
ras_mt_ma_2020_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2020_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2019_500m"
ras_mt_ma_2019_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2019_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2018_500m"
ras_mt_ma_2018_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2018_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2017_500m"
ras_mt_ma_2017_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2017_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2016_500m"
ras_mt_ma_2016_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2016_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2015_500m"
ras_mt_ma_2015_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2015_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2014_500m"
ras_mt_ma_2014_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2014_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2013_500m"
ras_mt_ma_2013_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2013_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2012_500m"
ras_mt_ma_2012_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2012_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2011_500m"
ras_mt_ma_2011_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2011_500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2010_500m"
ras_mt_ma_2010_500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2010_500m\\.tif$")) %>%
  map(rast) 

#### Metricas de paisagem 500m ####

# 2023

id_unico <- names(ras_mt_ma_2023_500m)

met_mt_ma_2023_500m <- map_df(seq_along(ras_mt_ma_2023_500m), function(i) {
  raster <- ras_mt_ma_2023_500m[[i]]
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

prox_mt_ma_2023_500m <- prox(ras_mt_ma_2023_500m, 3)

shape_mt_ma_2023_500m <- shape(ras_mt_ma_2023_500m, class_value = 3)

lsm_mt_ma_2023_500m <- met_mt_ma_2023_500m |>
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
  inner_join(prox_mt_ma_2023_500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2023_500m, by = "id_unico")



# 2022

id_unico <- names(ras_mt_ma_2022_500m)

met_mt_ma_2022_500m <- map_df(seq_along(ras_mt_ma_2022_500m), function(i) {
  raster <- ras_mt_ma_2022_500m[[i]]
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

prox_mt_ma_2022_500m <- prox(ras_mt_ma_2022_500m, 3)

shape_mt_ma_2022_500m <- shape(ras_mt_ma_2022_500m, class_value = 3)

lsm_mt_ma_2022_500m <- met_mt_ma_2022_500m |>
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
  inner_join(prox_mt_ma_2022_500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2022_500m, by = "id_unico")

# 2021

id_unico <- names(ras_mt_ma_2021_500m)

met_mt_ma_2021_500m <- map_df(seq_along(ras_mt_ma_2021_500m), function(i) {
  raster <- ras_mt_ma_2021_500m[[i]]
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

prox_mt_ma_2021_500m <- prox(ras_mt_ma_2021_500m, 3)

shape_mt_ma_2021_500m <- shape(ras_mt_ma_2021_500m, class_value = 3)

lsm_mt_ma_2021_500m <- met_mt_ma_2021_500m |>
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
  inner_join(prox_mt_ma_2021_500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2021_500m, by = "id_unico")

# 2020

id_unico <- names(ras_mt_ma_2020_500m)

met_mt_ma_2020_500m <- map_df(seq_along(ras_mt_ma_2020_500m), function(i) {
  raster <- ras_mt_ma_2020_500m[[i]]
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

prox_mt_ma_2020_500m <- prox(ras_mt_ma_2020_500m, 3)

shape_mt_ma_2020_500m <- shape(ras_mt_ma_2020_500m, class_value = 3)

lsm_mt_ma_2020_500m <- met_mt_ma_2020_500m |>
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
  inner_join(prox_mt_ma_2020_500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2020_500m, by = "id_unico")

# 2019

id_unico <- names(ras_mt_ma_2019_500m)

met_mt_ma_2019_500m <- map_df(seq_along(ras_mt_ma_2019_500m), function(i) {
  raster <- ras_mt_ma_2019_500m[[i]]
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

prox_mt_ma_2019_500m <- prox(ras_mt_ma_2019_500m, 3)

shape_mt_ma_2019_500m <- shape(ras_mt_ma_2019_500m, class_value = 3)

lsm_mt_ma_2019_500m <- met_mt_ma_2019_500m |>
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
  inner_join(prox_mt_ma_2019_500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2019_500m, by = "id_unico")

# 2018

id_unico <- names(ras_mt_ma_2018_500m)

met_mt_ma_2018_500m <- map_df(seq_along(ras_mt_ma_2018_500m), function(i) {
  raster <- ras_mt_ma_2018_500m[[i]]
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

prox_mt_ma_2018_500m <- prox(ras_mt_ma_2018_500m, 3)

shape_mt_ma_2018_500m <- shape(ras_mt_ma_2018_500m, class_value = 3)

lsm_mt_ma_2018_500m <- met_mt_ma_2018_500m |>
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
  inner_join(prox_mt_ma_2018_500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2018_500m, by = "id_unico")

# 2017

id_unico <- names(ras_mt_ma_2017_500m)

met_mt_ma_2017_500m <- map_df(seq_along(ras_mt_ma_2017_500m), function(i) {
  raster <- ras_mt_ma_2017_500m[[i]]
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

prox_mt_ma_2017_500m <- prox(ras_mt_ma_2017_500m, 3)

shape_mt_ma_2017_500m <- shape(ras_mt_ma_2017_500m, class_value = 3)

lsm_mt_ma_2017_500m <- met_mt_ma_2017_500m |>
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
  inner_join(prox_mt_ma_2017_500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2017_500m, by = "id_unico")


# 2016

id_unico <- names(ras_mt_ma_2016_500m)

met_mt_ma_2016_500m <- map_df(seq_along(ras_mt_ma_2016_500m), function(i) {
  raster <- ras_mt_ma_2016_500m[[i]]
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

prox_mt_ma_2016_500m <- prox(ras_mt_ma_2016_500m, 3)

shape_mt_ma_2016_500m <- shape(ras_mt_ma_2016_500m, class_value = 3)

lsm_mt_ma_2016_500m <- met_mt_ma_2016_500m |>
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
  inner_join(prox_mt_ma_2016_500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2016_500m, by = "id_unico")

# 2015

id_unico <- names(ras_mt_ma_2015_500m)

met_mt_ma_2015_500m <- map_df(seq_along(ras_mt_ma_2015_500m), function(i) {
  raster <- ras_mt_ma_2015_500m[[i]]
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

prox_mt_ma_2015_500m <- prox(ras_mt_ma_2015_500m, 3)

shape_mt_ma_2015_500m <- shape(ras_mt_ma_2015_500m, class_value = 3)

lsm_mt_ma_2015_500m <- met_mt_ma_2015_500m |>
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
  inner_join(prox_mt_ma_2015_500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2015_500m, by = "id_unico")

# 2014

id_unico <- names(ras_mt_ma_2014_500m)

met_mt_ma_2014_500m <- map_df(seq_along(ras_mt_ma_2014_500m), function(i) {
  raster <- ras_mt_ma_2014_500m[[i]]
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

prox_mt_ma_2014_500m <- prox(ras_mt_ma_2014_500m, 3)

shape_mt_ma_2014_500m <- shape(ras_mt_ma_2014_500m, class_value = 3)

lsm_mt_ma_2014_500m <- met_mt_ma_2014_500m |>
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
  inner_join(prox_mt_ma_2014_500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2014_500m, by = "id_unico")

# 2013

id_unico <- names(ras_mt_ma_2013_500m)

met_mt_ma_2013_500m <- map_df(seq_along(ras_mt_ma_2013_500m), function(i) {
  raster <- ras_mt_ma_2013_500m[[i]]
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

prox_mt_ma_2013_500m <- prox(ras_mt_ma_2013_500m, 3)

shape_mt_ma_2013_500m <- shape(ras_mt_ma_2013_500m, class_value = 3)

lsm_mt_ma_2013_500m <- met_mt_ma_2013_500m |>
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
  inner_join(prox_mt_ma_2013_500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2013_500m, by = "id_unico")

# 2012
id_unico <- names(ras_mt_ma_2012_500m)

met_mt_ma_2012_500m <- map_df(seq_along(ras_mt_ma_2012_500m), function(i) {
  raster <- ras_mt_ma_2012_500m[[i]]
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

prox_mt_ma_2012_500m <- prox(ras_mt_ma_2012_500m, 3)

shape_mt_ma_2012_500m <- shape(ras_mt_ma_2012_500m, class_value = 3)

lsm_mt_ma_2012_500m <- met_mt_ma_2012_500m |>
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
  inner_join(prox_mt_ma_2012_500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2012_500m, by = "id_unico")

# 2011

id_unico <- names(ras_mt_ma_2011_500m)

met_mt_ma_2011_500m <- map_df(seq_along(ras_mt_ma_2011_500m), function(i) {
  raster <- ras_mt_ma_2011_500m[[i]]
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

prox_mt_ma_2011_500m <- prox(ras_mt_ma_2011_500m, 3)

shape_mt_ma_2011_500m <- shape(ras_mt_ma_2011_500m, class_value = 3)

lsm_mt_ma_2011_500m <- met_mt_ma_2011_500m |>
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
  inner_join(prox_mt_ma_2011_500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2011_500m, by = "id_unico")
# 2010

id_unico <- names(ras_mt_ma_2010_500m)

met_mt_ma_2010_500m <- map_df(seq_along(ras_mt_ma_2010_500m), function(i) {
  raster <- ras_mt_ma_2010_500m[[i]]
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

prox_mt_ma_2010_500m <- prox(ras_mt_ma_2010_500m, 3)

shape_mt_ma_2010_500m <- shape(ras_mt_ma_2010_500m, class_value = 3)

lsm_mt_ma_2010_500m <- met_mt_ma_2010_500m |>
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
  inner_join(prox_mt_ma_2010_500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2010_500m, by = "id_unico")


#### Unificando tabelas 500m ####

lsm_mt_ma_500m <- bind_rows(lsm_mt_ma_2023_500m, lsm_mt_ma_2022_500m, lsm_mt_ma_2021_500m, lsm_mt_ma_2020_500m, lsm_mt_ma_2019_500m, lsm_mt_ma_2018_500m, lsm_mt_ma_2017_500m, lsm_mt_ma_2016_500m, lsm_mt_ma_2015_500m, lsm_mt_ma_2014_500m, lsm_mt_ma_2013_500m, lsm_mt_ma_2012_500m, lsm_mt_ma_2011_500m, lsm_mt_ma_2010_500m) |>
  dplyr::select(-pland_0,-lpi_0, -ed_0, -pd_0, -np_0, -lpi_12, -ed_12, -pd_12, -np_12, -lpi_15, -ed_15, -pd_15, -np_15, -lpi_9, -ed_9, -pd_9, -np_9, -np_24, -pd_24, -lpi_24, -ed_24,  -pland_24)|>
  mutate(Bin = str_extract(id_unico, "^[01]"),
         Bin = as.factor(Bin),     
         Ano = as.factor(Ano))|>
  mutate(pland_9 = ifelse(is.na(pland_9), 0, pland_9))|>
  mutate(pland_15 = ifelse(is.na(pland_15), 0, pland_15))|>
  mutate(pland_12 = ifelse(is.na(pland_12), 0, pland_12)) |>
  rename_with(~ paste0(., "_500m"))

lsm_mt_ma_500m_sem_na <- lsm_mt_ma_500m|>
  na.omit() 

View(lsm_mt_ma_500m_sem_na)

# Identifica as linhas com Bin == 0
zeros_500m <- which(lsm_mt_ma_500m_sem_na$Bin_500m == 0)

# Seleciona zeros aleatórios para remover
remove_rows_500m <- sample(zeros_500m, min(587, length(zeros_500m)))

# Filtra o tibble
lsm_mt_ma_500m_total <- lsm_mt_ma_500m_sem_na[-remove_rows_500m, ]
View(lsm_mt_ma_500m_total)

write_xlsx(lsm_mt_ma_500m_total, "lsm_mt_ma_500m.xlsx")

lsm_mt_ma_500m_total <- read_excel("D:/Analise estatistica M. tridactyla - MA/Analise estatistica M. tridactyla - MA/lsm_mt_ma_500m.xlsx", sheet = 1,na = c(" ", "NA"))

#### Filtando pontos proximos as coberturas florestais ####

mt_ma_flo <- inner_join(lsm_mt_ma_500m_total, mt_ma, by = c("id_unico_500m" = "id_unico"))|>
  select(Bin, id_unico_500m, Ano, geometry)|>
  dplyr::rename("id_unico" = "id_unico_500m")|>
  st_as_sf(crs = 5880)

write_xlsx(mt_ma_flo, "mt_ma_flo.xlsx")

st_write(mt_ma_flo, "mt_ma_flo.shp")

#### Gerando buffer 2000m ####

bf_mt_ma_2023_2000m <- mt_ma_flo|>
  filter(Ano == "2023") |>
  st_buffer(dist = 1000)

bf_mt_ma_2022_2000m <- mt_ma_flo|>
  filter(Ano == "2022") |>
  st_buffer(dist = 1000)

bf_mt_ma_2021_2000m <- mt_ma_flo|>
  filter(Ano == "2021") |>
  st_buffer(dist = 1000)

bf_mt_ma_2020_2000m <- mt_ma_flo|>
  filter(Ano == "2020") |>
  st_buffer(dist = 1000)

bf_mt_ma_2019_2000m <- mt_ma_flo|>
  filter(Ano == "2019") |>
  st_buffer(dist = 1000)

bf_mt_ma_2018_2000m <- mt_ma_flo|>
  filter(Ano == "2018") |>
  st_buffer(dist = 1000)

bf_mt_ma_2017_2000m <- mt_ma_flo|>
  filter(Ano == "2017") |>
  st_buffer(dist = 1000)

bf_mt_ma_2016_2000m <- mt_ma_flo|>
  filter(Ano == "2016") |>
  st_buffer(dist = 1000)

bf_mt_ma_2015_2000m <- mt_ma_flo|>
  filter(Ano == "2015") |>
  st_buffer(dist = 1000)

bf_mt_ma_2014_2000m <- mt_ma_flo|>
  filter(Ano == "2014") |>
  st_buffer(dist = 1000)

bf_mt_ma_2013_2000m <- mt_ma_flo|>
  filter(Ano == "2013") |>
  st_buffer(dist = 1000)

bf_mt_ma_2012_2000m <- mt_ma_flo|>
  filter(Ano == "2012") |>
  st_buffer(dist = 1000)

bf_mt_ma_2011_2000m <- mt_ma_flo|>
  filter(Ano == "2011") |>
  st_buffer(dist = 1000)

bf_mt_ma_2010_2000m <- mt_ma_flo|>
  filter(Ano == "2010") |>
  st_buffer(dist = 1000)



#### Cortando raster 2000m ####

# 2023

ras_mt_ma_2023_2000m <- list()

for (i in 1:nrow(bf_mt_ma_2023_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2023_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2023_2000m <- crop(ma_2023, bf_mt_ma_2023_2000m[i, ])
  mask_mt_ma_2023_2000m <- mask(crop_mt_ma_2023_2000m, bf_mt_ma_2023_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2023_2000m[[buffer_id]] <- mask_mt_ma_2023_2000m
}

output_dir <- "ras_mt_ma_2023_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2023_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2023_2000m.tif"))
  writeRaster(
    ras_mt_ma_2023_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2022

ras_mt_ma_2022_2000m <- list()

for (i in 1:nrow(bf_mt_ma_2022_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2022_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2022_2000m <- crop(ma_2022, bf_mt_ma_2022_2000m[i, ])
  mask_mt_ma_2022_2000m <- mask(crop_mt_ma_2022_2000m, bf_mt_ma_2022_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2022_2000m[[buffer_id]] <- mask_mt_ma_2022_2000m
}

output_dir <- "ras_mt_ma_2022_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2022_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2022_2000m.tif"))
  writeRaster(
    ras_mt_ma_2022_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2021

ras_mt_ma_2021_2000m <- list()

for (i in 1:nrow(bf_mt_ma_2021_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2021_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2021_2000m <- crop(ma_2021, bf_mt_ma_2021_2000m[i, ])
  mask_mt_ma_2021_2000m <- mask(crop_mt_ma_2021_2000m, bf_mt_ma_2021_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2021_2000m[[buffer_id]] <- mask_mt_ma_2021_2000m
}

output_dir <- "ras_mt_ma_2021_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2021_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2021_2000m.tif"))
  writeRaster(
    ras_mt_ma_2021_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2020

ras_mt_ma_2020_2000m <- list()

for (i in 1:nrow(bf_mt_ma_2020_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2020_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2020_2000m <- crop(ma_2020, bf_mt_ma_2020_2000m[i, ])
  mask_mt_ma_2020_2000m <- mask(crop_mt_ma_2020_2000m, bf_mt_ma_2020_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2020_2000m[[buffer_id]] <- mask_mt_ma_2020_2000m
}

output_dir <- "ras_mt_ma_2020_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2020_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2020_2000m.tif"))
  writeRaster(
    ras_mt_ma_2020_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2019

ras_mt_ma_2019_2000m <- list()

for (i in 1:nrow(bf_mt_ma_2019_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2019_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2019_2000m <- crop(ma_2019, bf_mt_ma_2019_2000m[i, ])
  mask_mt_ma_2019_2000m <- mask(crop_mt_ma_2019_2000m, bf_mt_ma_2019_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2019_2000m[[buffer_id]] <- mask_mt_ma_2019_2000m
}

output_dir <- "ras_mt_ma_2019_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2019_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2019_2000m.tif"))
  writeRaster(
    ras_mt_ma_2019_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2018

ras_mt_ma_2018_2000m <- list()

for (i in 1:nrow(bf_mt_ma_2018_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2018_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2018_2000m <- crop(ma_2018, bf_mt_ma_2018_2000m[i, ])
  mask_mt_ma_2018_2000m <- mask(crop_mt_ma_2018_2000m, bf_mt_ma_2018_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2018_2000m[[buffer_id]] <- mask_mt_ma_2018_2000m
}

output_dir <- "ras_mt_ma_2018_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2018_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2018_2000m.tif"))
  writeRaster(
    ras_mt_ma_2018_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2017

ras_mt_ma_2017_2000m <- list()

for (i in 1:nrow(bf_mt_ma_2017_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2017_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2017_2000m <- crop(ma_2017, bf_mt_ma_2017_2000m[i, ])
  mask_mt_ma_2017_2000m <- mask(crop_mt_ma_2017_2000m, bf_mt_ma_2017_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2017_2000m[[buffer_id]] <- mask_mt_ma_2017_2000m
}

output_dir <- "ras_mt_ma_2017_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2017_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2017_2000m.tif"))
  writeRaster(
    ras_mt_ma_2017_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2017

ras_mt_ma_2017_2000m <- list()

for (i in 1:nrow(bf_mt_ma_2017_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2017_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2017_2000m <- crop(ma_2017, bf_mt_ma_2017_2000m[i, ])
  mask_mt_ma_2017_2000m <- mask(crop_mt_ma_2017_2000m, bf_mt_ma_2017_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2017_2000m[[buffer_id]] <- mask_mt_ma_2017_2000m
}

output_dir <- "ras_mt_ma_2017_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2017_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2017_2000m.tif"))
  writeRaster(
    ras_mt_ma_2017_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2016

ras_mt_ma_2016_2000m <- list()

for (i in 1:nrow(bf_mt_ma_2016_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2016_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2016_2000m <- crop(ma_2016, bf_mt_ma_2016_2000m[i, ])
  mask_mt_ma_2016_2000m <- mask(crop_mt_ma_2016_2000m, bf_mt_ma_2016_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2016_2000m[[buffer_id]] <- mask_mt_ma_2016_2000m
}

output_dir <- "ras_mt_ma_2016_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2016_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2016_2000m.tif"))
  writeRaster(
    ras_mt_ma_2016_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2015

ras_mt_ma_2015_2000m <- list()

for (i in 1:nrow(bf_mt_ma_2015_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2015_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2015_2000m <- crop(ma_2015, bf_mt_ma_2015_2000m[i, ])
  mask_mt_ma_2015_2000m <- mask(crop_mt_ma_2015_2000m, bf_mt_ma_2015_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2015_2000m[[buffer_id]] <- mask_mt_ma_2015_2000m
}

output_dir <- "ras_mt_ma_2015_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2015_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2015_2000m.tif"))
  writeRaster(
    ras_mt_ma_2015_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2014

ras_mt_ma_2014_2000m <- list()

for (i in 1:nrow(bf_mt_ma_2014_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2014_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2014_2000m <- crop(ma_2014, bf_mt_ma_2014_2000m[i, ])
  mask_mt_ma_2014_2000m <- mask(crop_mt_ma_2014_2000m, bf_mt_ma_2014_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2014_2000m[[buffer_id]] <- mask_mt_ma_2014_2000m
}

output_dir <- "ras_mt_ma_2014_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2014_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2014_2000m.tif"))
  writeRaster(
    ras_mt_ma_2014_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2013

ras_mt_ma_2013_2000m <- list()

for (i in 1:nrow(bf_mt_ma_2013_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2013_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2013_2000m <- crop(ma_2013, bf_mt_ma_2013_2000m[i, ])
  mask_mt_ma_2013_2000m <- mask(crop_mt_ma_2013_2000m, bf_mt_ma_2013_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2013_2000m[[buffer_id]] <- mask_mt_ma_2013_2000m
}

output_dir <- "ras_mt_ma_2013_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2013_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2013_2000m.tif"))
  writeRaster(
    ras_mt_ma_2013_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2012

ras_mt_ma_2012_2000m <- list()

for (i in 1:nrow(bf_mt_ma_2012_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2012_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2012_2000m <- crop(ma_2012, bf_mt_ma_2012_2000m[i, ])
  mask_mt_ma_2012_2000m <- mask(crop_mt_ma_2012_2000m, bf_mt_ma_2012_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2012_2000m[[buffer_id]] <- mask_mt_ma_2012_2000m
}

output_dir <- "ras_mt_ma_2012_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2012_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2012_2000m.tif"))
  writeRaster(
    ras_mt_ma_2012_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2011

ras_mt_ma_2011_2000m <- list()

for (i in 1:nrow(bf_mt_ma_2011_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2011_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2011_2000m <- crop(ma_2011, bf_mt_ma_2011_2000m[i, ])
  mask_mt_ma_2011_2000m <- mask(crop_mt_ma_2011_2000m, bf_mt_ma_2011_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2011_2000m[[buffer_id]] <- mask_mt_ma_2011_2000m
}

output_dir <- "ras_mt_ma_2011_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2011_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2011_2000m.tif"))
  writeRaster(
    ras_mt_ma_2011_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2010

ras_mt_ma_2010_2000m <- list()

for (i in 1:nrow(bf_mt_ma_2010_2000m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2010_2000m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2010_2000m <- crop(ma_2010, bf_mt_ma_2010_2000m[i, ])
  mask_mt_ma_2010_2000m <- mask(crop_mt_ma_2010_2000m, bf_mt_ma_2010_2000m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2010_2000m[[buffer_id]] <- mask_mt_ma_2010_2000m
}

output_dir <- "ras_mt_ma_2010_2000m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2010_2000m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2010_2000m.tif"))
  writeRaster(
    ras_mt_ma_2010_2000m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}




# Chamando os recortes

output_dir <- "ras_mt_ma_2023_2000m"
ras_mt_ma_2023_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2023_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2022_2000m"
ras_mt_ma_2022_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2022_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2021_2000m"
ras_mt_ma_2021_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2021_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2020_2000m"
ras_mt_ma_2020_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2020_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2019_2000m"
ras_mt_ma_2019_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2019_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2018_2000m"
ras_mt_ma_2018_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2018_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2017_2000m"
ras_mt_ma_2017_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2017_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2015_2000m"
ras_mt_ma_2015_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2015_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2014_2000m"
ras_mt_ma_2014_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2014_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2013_2000m"
ras_mt_ma_2013_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2013_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2012_2000m"
ras_mt_ma_2012_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2012_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2011_2000m"
ras_mt_ma_2011_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2011_2000m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2009_2000m"
ras_mt_ma_2010_2000m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2010_2000m\\.tif$")) %>%
  map(rast) 


#### Metricas de paisagem 2000m ####

# 2023

id_unico <- names(ras_mt_ma_2023_2000m)

met_mt_ma_2023_2000m <- map_df(seq_along(ras_mt_ma_2023_2000m), function(i) {
  raster <- ras_mt_ma_2023_2000m[[i]]
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

prox_mt_ma_2023_2000m <- prox(ras_mt_ma_2023_2000m, 3)

shape_mt_ma_2023_2000m <- shape(ras_mt_ma_2023_2000m, class_value = 3)

lsm_mt_ma_2023_2000m <- met_mt_ma_2023_2000m |>
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
  inner_join(prox_mt_ma_2023_2000m, by = "id_unico")|>
  inner_join(shape_mt_ma_2023_2000m, by = "id_unico")

# 2022

id_unico <- names(ras_mt_ma_2022_2000m)

met_mt_ma_2022_2000m <- map_df(seq_along(ras_mt_ma_2022_2000m), function(i) {
  raster <- ras_mt_ma_2022_2000m[[i]]
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

prox_mt_ma_2022_2000m <- prox(ras_mt_ma_2022_2000m, 3)

shape_mt_ma_2022_2000m <- shape(ras_mt_ma_2022_2000m, class_value = 3)

lsm_mt_ma_2022_2000m <- met_mt_ma_2022_2000m |>
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
  inner_join(prox_mt_ma_2022_2000m, by = "id_unico")|>
  inner_join(shape_mt_ma_2022_2000m, by = "id_unico")

# 2021

id_unico <- names(ras_mt_ma_2021_2000m)

met_mt_ma_2021_2000m <- map_df(seq_along(ras_mt_ma_2021_2000m), function(i) {
  raster <- ras_mt_ma_2021_2000m[[i]]
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

prox_mt_ma_2021_2000m <- prox(ras_mt_ma_2021_2000m, 3)

shape_mt_ma_2021_2000m <- shape(ras_mt_ma_2021_2000m, class_value = 3)

lsm_mt_ma_2021_2000m <- met_mt_ma_2021_2000m |>
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
  inner_join(prox_mt_ma_2021_2000m, by = "id_unico")|>
  inner_join(shape_mt_ma_2021_2000m, by = "id_unico")

# 2020

id_unico <- names(ras_mt_ma_2020_2000m)

met_mt_ma_2020_2000m <- map_df(seq_along(ras_mt_ma_2020_2000m), function(i) {
  raster <- ras_mt_ma_2020_2000m[[i]]
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

prox_mt_ma_2020_2000m <- prox(ras_mt_ma_2020_2000m, 3)

shape_mt_ma_2020_2000m <- shape(ras_mt_ma_2020_2000m, class_value = 3)

lsm_mt_ma_2020_2000m <- met_mt_ma_2020_2000m |>
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
  inner_join(prox_mt_ma_2020_2000m, by = "id_unico")|>
  inner_join(shape_mt_ma_2020_2000m, by = "id_unico")

# 2019

id_unico <- names(ras_mt_ma_2019_2000m)

met_mt_ma_2019_2000m <- map_df(seq_along(ras_mt_ma_2019_2000m), function(i) {
  raster <- ras_mt_ma_2019_2000m[[i]]
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

prox_mt_ma_2019_2000m <- prox(ras_mt_ma_2019_2000m, 3)

shape_mt_ma_2019_2000m <- shape(ras_mt_ma_2019_2000m, class_value = 3)

lsm_mt_ma_2019_2000m <- met_mt_ma_2019_2000m |>
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
  inner_join(prox_mt_ma_2019_2000m, by = "id_unico")|>
  inner_join(shape_mt_ma_2019_2000m, by = "id_unico")

# 2018

id_unico <- names(ras_mt_ma_2018_2000m)

met_mt_ma_2018_2000m <- map_df(seq_along(ras_mt_ma_2018_2000m), function(i) {
  raster <- ras_mt_ma_2018_2000m[[i]]
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

prox_mt_ma_2018_2000m <- prox(ras_mt_ma_2018_2000m, 3)

shape_mt_ma_2018_2000m <- shape(ras_mt_ma_2018_2000m, class_value = 3)

lsm_mt_ma_2018_2000m <- met_mt_ma_2018_2000m |>
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
  inner_join(prox_mt_ma_2018_2000m, by = "id_unico")|>
  inner_join(shape_mt_ma_2018_2000m, by = "id_unico")

# 2017

id_unico <- names(ras_mt_ma_2017_2000m)

met_mt_ma_2017_2000m <- map_df(seq_along(ras_mt_ma_2017_2000m), function(i) {
  raster <- ras_mt_ma_2017_2000m[[i]]
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

prox_mt_ma_2017_2000m <- prox(ras_mt_ma_2017_2000m, 3)

shape_mt_ma_2017_2000m <- shape(ras_mt_ma_2017_2000m, class_value = 3)

lsm_mt_ma_2017_2000m <- met_mt_ma_2017_2000m |>
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
  inner_join(prox_mt_ma_2017_2000m, by = "id_unico")|>
  inner_join(shape_mt_ma_2017_2000m, by = "id_unico")


# 2016

id_unico <- names(ras_mt_ma_2016_2000m)

met_mt_ma_2016_2000m <- map_df(seq_along(ras_mt_ma_2016_2000m), function(i) {
  raster <- ras_mt_ma_2016_2000m[[i]]
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

prox_mt_ma_2016_2000m <- prox(ras_mt_ma_2016_2000m, 3)

shape_mt_ma_2016_2000m <- shape(ras_mt_ma_2016_2000m, class_value = 3)

lsm_mt_ma_2016_2000m <- met_mt_ma_2016_2000m |>
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
  inner_join(prox_mt_ma_2016_2000m, by = "id_unico")|>
  inner_join(shape_mt_ma_2016_2000m, by = "id_unico")

# 2015

id_unico <- names(ras_mt_ma_2015_2000m)

met_mt_ma_2015_2000m <- map_df(seq_along(ras_mt_ma_2015_2000m), function(i) {
  raster <- ras_mt_ma_2015_2000m[[i]]
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

prox_mt_ma_2015_2000m <- prox(ras_mt_ma_2015_2000m, 3)

shape_mt_ma_2015_2000m <- shape(ras_mt_ma_2015_2000m, class_value = 3)

lsm_mt_ma_2015_2000m <- met_mt_ma_2015_2000m |>
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
  inner_join(prox_mt_ma_2015_2000m, by = "id_unico")|>
  inner_join(shape_mt_ma_2015_2000m, by = "id_unico")

# 2014

id_unico <- names(ras_mt_ma_2014_2000m)

met_mt_ma_2014_2000m <- map_df(seq_along(ras_mt_ma_2014_2000m), function(i) {
  raster <- ras_mt_ma_2014_2000m[[i]]
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

prox_mt_ma_2014_2000m <- prox(ras_mt_ma_2014_2000m, 3)

shape_mt_ma_2014_2000m <- shape(ras_mt_ma_2014_2000m, class_value = 3)

lsm_mt_ma_2014_2000m <- met_mt_ma_2014_2000m |>
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
  inner_join(prox_mt_ma_2014_2000m, by = "id_unico")|>
  inner_join(shape_mt_ma_2014_2000m, by = "id_unico")

# 2013

id_unico <- names(ras_mt_ma_2013_2000m)

met_mt_ma_2013_2000m <- map_df(seq_along(ras_mt_ma_2013_2000m), function(i) {
  raster <- ras_mt_ma_2013_2000m[[i]]
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

prox_mt_ma_2013_2000m <- prox(ras_mt_ma_2013_2000m, 3)

shape_mt_ma_2013_2000m <- shape(ras_mt_ma_2013_2000m, class_value = 3)

lsm_mt_ma_2013_2000m <- met_mt_ma_2013_2000m |>
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
  inner_join(prox_mt_ma_2013_2000m, by = "id_unico")|>
  inner_join(shape_mt_ma_2013_2000m, by = "id_unico")

# 2012
id_unico <- names(ras_mt_ma_2012_2000m)

met_mt_ma_2012_2000m <- map_df(seq_along(ras_mt_ma_2012_2000m), function(i) {
  raster <- ras_mt_ma_2012_2000m[[i]]
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

prox_mt_ma_2012_2000m <- prox(ras_mt_ma_2012_2000m, 3)

shape_mt_ma_2012_2000m <- shape(ras_mt_ma_2012_2000m, class_value = 3)

lsm_mt_ma_2012_2000m <- met_mt_ma_2012_2000m |>
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
  inner_join(prox_mt_ma_2012_2000m, by = "id_unico")|>
  inner_join(shape_mt_ma_2012_2000m, by = "id_unico")

# 2011

id_unico <- names(ras_mt_ma_2011_2000m)

met_mt_ma_2011_2000m <- map_df(seq_along(ras_mt_ma_2011_2000m), function(i) {
  raster <- ras_mt_ma_2011_2000m[[i]]
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

prox_mt_ma_2011_2000m <- prox(ras_mt_ma_2011_2000m, 3)

shape_mt_ma_2011_2000m <- shape(ras_mt_ma_2011_2000m, class_value = 3)

lsm_mt_ma_2011_2000m <- met_mt_ma_2011_2000m |>
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
  inner_join(prox_mt_ma_2011_2000m, by = "id_unico")|>
  inner_join(shape_mt_ma_2011_2000m, by = "id_unico")
# 2010

id_unico <- names(ras_mt_ma_2010_2000m)

met_mt_ma_2010_2000m <- map_df(seq_along(ras_mt_ma_2010_2000m), function(i) {
  raster <- ras_mt_ma_2010_2000m[[i]]
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

prox_mt_ma_2010_2000m <- prox(ras_mt_ma_2010_2000m, 3)

shape_mt_ma_2010_2000m <- shape(ras_mt_ma_2010_2000m, class_value = 3)

lsm_mt_ma_2010_2000m <- met_mt_ma_2010_2000m |>
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
  inner_join(prox_mt_ma_2010_2000m, by = "id_unico")|>
  inner_join(shape_mt_ma_2010_2000m, by = "id_unico")



#### Unificando tabelas 2000m ####

lsm_mt_ma_2000m <- bind_rows(lsm_mt_ma_2023_2000m, lsm_mt_ma_2022_2000m, lsm_mt_ma_2021_2000m, lsm_mt_ma_2020_2000m, lsm_mt_ma_2019_2000m, lsm_mt_ma_2018_2000m, lsm_mt_ma_2017_2000m, lsm_mt_ma_2016_2000m, lsm_mt_ma_2015_2000m, lsm_mt_ma_2014_2000m, lsm_mt_ma_2013_2000m, lsm_mt_ma_2012_2000m, lsm_mt_ma_2011_2000m, lsm_mt_ma_2010_2000m) |>
  dplyr::select(-pland_0,-lpi_0, -ed_0, -pd_0, -np_0, -lpi_12, -ed_12, -pd_12, -np_12, -lpi_15, -ed_15, -pd_15, -np_15, -lpi_9, -ed_9, -pd_9, -np_9, -np_24, -pd_24, -lpi_24, -ed_24,  -pland_24)|>
  mutate(Bin = str_extract(id_unico, "^[01]"),
         Bin = as.factor(Bin),     
         Ano = as.factor(Ano))|>
  mutate(pland_9 = ifelse(is.na(pland_9), 0, pland_9))|>
  mutate(pland_15 = ifelse(is.na(pland_15), 0, pland_15))|>
  mutate(pland_12 = ifelse(is.na(pland_12), 0, pland_12)) |>
  rename_with(~ paste0(., "_2000m"))


lsm_mt_ma_2000m_sem_na <- lsm_mt_ma_2000m|>
  na.omit() 

View(lsm_mt_ma_2000m_sem_na)

write_xlsx(lsm_mt_ma_2000m_sem_na, "lsm_mt_ma_2000m.xlsx")



#### Gerando buffer 3500m ####

bf_mt_ma_2023_3500m <- mt_ma_flo|>
  filter(Ano == "2023") |>
  st_buffer(dist = 3500)

bf_mt_ma_2022_3500m <- mt_ma_flo|>
  filter(Ano == "2022") |>
  st_buffer(dist = 3500)

bf_mt_ma_2021_3500m <- mt_ma_flo|>
  filter(Ano == "2021") |>
  st_buffer(dist = 3500)

bf_mt_ma_2020_3500m <- mt_ma_flo|>
  filter(Ano == "2020") |>
  st_buffer(dist = 3500)

bf_mt_ma_2019_3500m <- mt_ma_flo|>
  filter(Ano == "2019") |>
  st_buffer(dist = 3500)

bf_mt_ma_2018_3500m <- mt_ma_flo|>
  filter(Ano == "2018") |>
  st_buffer(dist = 3500)

bf_mt_ma_2017_3500m <- mt_ma_flo|>
  filter(Ano == "2017") |>
  st_buffer(dist = 3500)

bf_mt_ma_2016_3500m <- mt_ma_flo|>
  filter(Ano == "2016") |>
  st_buffer(dist = 3500)

bf_mt_ma_2015_3500m <- mt_ma_flo|>
  filter(Ano == "2015") |>
  st_buffer(dist = 3500)

bf_mt_ma_2014_3500m <- mt_ma_flo|>
  filter(Ano == "2014") |>
  st_buffer(dist = 3500)

bf_mt_ma_2013_3500m <- mt_ma_flo|>
  filter(Ano == "2013") |>
  st_buffer(dist = 3500)

bf_mt_ma_2012_3500m <- mt_ma_flo|>
  filter(Ano == "2012") |>
  st_buffer(dist = 3500)

bf_mt_ma_2011_3500m <- mt_ma_flo|>
  filter(Ano == "2011") |>
  st_buffer(dist = 3500)

bf_mt_ma_2010_3500m <- mt_ma_flo|>
  filter(Ano == "2010") |>
  st_buffer(dist = 3500)


#### Cortando raster 3500m ####

# 2023

ras_mt_ma_2023_3500m <- list()

for (i in 1:nrow(bf_mt_ma_2023_3500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2023_3500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2023_3500m <- crop(ma_2023, bf_mt_ma_2023_3500m[i, ])
  mask_mt_ma_2023_3500m <- mask(crop_mt_ma_2023_3500m, bf_mt_ma_2023_3500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2023_3500m[[buffer_id]] <- mask_mt_ma_2023_3500m
}

output_dir <- "ras_mt_ma_2023_3500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2023_3500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2023_3500m.tif"))
  writeRaster(
    ras_mt_ma_2023_3500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2022

ras_mt_ma_2022_3500m <- list()

for (i in 1:nrow(bf_mt_ma_2022_3500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2022_3500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2022_3500m <- crop(ma_2022, bf_mt_ma_2022_3500m[i, ])
  mask_mt_ma_2022_3500m <- mask(crop_mt_ma_2022_3500m, bf_mt_ma_2022_3500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2022_3500m[[buffer_id]] <- mask_mt_ma_2022_3500m
}

output_dir <- "ras_mt_ma_2022_3500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2022_3500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2022_3500m.tif"))
  writeRaster(
    ras_mt_ma_2022_3500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2021

ras_mt_ma_2021_3500m <- list()

for (i in 1:nrow(bf_mt_ma_2021_3500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2021_3500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2021_3500m <- crop(ma_2021, bf_mt_ma_2021_3500m[i, ])
  mask_mt_ma_2021_3500m <- mask(crop_mt_ma_2021_3500m, bf_mt_ma_2021_3500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2021_3500m[[buffer_id]] <- mask_mt_ma_2021_3500m
}

output_dir <- "ras_mt_ma_2021_3500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2021_3500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2021_3500m.tif"))
  writeRaster(
    ras_mt_ma_2021_3500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2020

ras_mt_ma_2020_3500m <- list()

for (i in 1:nrow(bf_mt_ma_2020_3500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2020_3500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2020_3500m <- crop(ma_2020, bf_mt_ma_2020_3500m[i, ])
  mask_mt_ma_2020_3500m <- mask(crop_mt_ma_2020_3500m, bf_mt_ma_2020_3500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2020_3500m[[buffer_id]] <- mask_mt_ma_2020_3500m
}

output_dir <- "ras_mt_ma_2020_3500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2020_3500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2020_3500m.tif"))
  writeRaster(
    ras_mt_ma_2020_3500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2019

ras_mt_ma_2019_3500m <- list()

for (i in 1:nrow(bf_mt_ma_2019_3500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2019_3500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2019_3500m <- crop(ma_2019, bf_mt_ma_2019_3500m[i, ])
  mask_mt_ma_2019_3500m <- mask(crop_mt_ma_2019_3500m, bf_mt_ma_2019_3500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2019_3500m[[buffer_id]] <- mask_mt_ma_2019_3500m
}

output_dir <- "ras_mt_ma_2019_3500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2019_3500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2019_3500m.tif"))
  writeRaster(
    ras_mt_ma_2019_3500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2018

ras_mt_ma_2018_3500m <- list()

for (i in 1:nrow(bf_mt_ma_2018_3500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2018_3500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2018_3500m <- crop(ma_2018, bf_mt_ma_2018_3500m[i, ])
  mask_mt_ma_2018_3500m <- mask(crop_mt_ma_2018_3500m, bf_mt_ma_2018_3500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2018_3500m[[buffer_id]] <- mask_mt_ma_2018_3500m
}

output_dir <- "ras_mt_ma_2018_3500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2018_3500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2018_3500m.tif"))
  writeRaster(
    ras_mt_ma_2018_3500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2017

ras_mt_ma_2017_3500m <- list()

for (i in 1:nrow(bf_mt_ma_2017_3500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2017_3500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2017_3500m <- crop(ma_2017, bf_mt_ma_2017_3500m[i, ])
  mask_mt_ma_2017_3500m <- mask(crop_mt_ma_2017_3500m, bf_mt_ma_2017_3500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2017_3500m[[buffer_id]] <- mask_mt_ma_2017_3500m
}

output_dir <- "ras_mt_ma_2017_3500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2017_3500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2017_3500m.tif"))
  writeRaster(
    ras_mt_ma_2017_3500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2017

ras_mt_ma_2017_3500m <- list()

for (i in 1:nrow(bf_mt_ma_2017_3500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2017_3500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2017_3500m <- crop(ma_2017, bf_mt_ma_2017_3500m[i, ])
  mask_mt_ma_2017_3500m <- mask(crop_mt_ma_2017_3500m, bf_mt_ma_2017_3500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2017_3500m[[buffer_id]] <- mask_mt_ma_2017_3500m
}

output_dir <- "ras_mt_ma_2017_3500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2017_3500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2017_3500m.tif"))
  writeRaster(
    ras_mt_ma_2017_3500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2016

ras_mt_ma_2016_3500m <- list()

for (i in 1:nrow(bf_mt_ma_2016_3500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2016_3500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2016_3500m <- crop(ma_2016, bf_mt_ma_2016_3500m[i, ])
  mask_mt_ma_2016_3500m <- mask(crop_mt_ma_2016_3500m, bf_mt_ma_2016_3500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2016_3500m[[buffer_id]] <- mask_mt_ma_2016_3500m
}

output_dir <- "ras_mt_ma_2016_3500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2016_3500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2016_3500m.tif"))
  writeRaster(
    ras_mt_ma_2016_3500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2015

ras_mt_ma_2015_3500m <- list()

for (i in 1:nrow(bf_mt_ma_2015_3500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2015_3500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2015_3500m <- crop(ma_2015, bf_mt_ma_2015_3500m[i, ])
  mask_mt_ma_2015_3500m <- mask(crop_mt_ma_2015_3500m, bf_mt_ma_2015_3500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2015_3500m[[buffer_id]] <- mask_mt_ma_2015_3500m
}

output_dir <- "ras_mt_ma_2015_3500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2015_3500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2015_3500m.tif"))
  writeRaster(
    ras_mt_ma_2015_3500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2014

ras_mt_ma_2014_3500m <- list()

for (i in 1:nrow(bf_mt_ma_2014_3500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2014_3500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2014_3500m <- crop(ma_2014, bf_mt_ma_2014_3500m[i, ])
  mask_mt_ma_2014_3500m <- mask(crop_mt_ma_2014_3500m, bf_mt_ma_2014_3500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2014_3500m[[buffer_id]] <- mask_mt_ma_2014_3500m
}

output_dir <- "ras_mt_ma_2014_3500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2014_3500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2014_3500m.tif"))
  writeRaster(
    ras_mt_ma_2014_3500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2013

ras_mt_ma_2013_3500m <- list()

for (i in 1:nrow(bf_mt_ma_2013_3500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2013_3500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2013_3500m <- crop(ma_2013, bf_mt_ma_2013_3500m[i, ])
  mask_mt_ma_2013_3500m <- mask(crop_mt_ma_2013_3500m, bf_mt_ma_2013_3500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2013_3500m[[buffer_id]] <- mask_mt_ma_2013_3500m
}

output_dir <- "ras_mt_ma_2013_3500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2013_3500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2013_3500m.tif"))
  writeRaster(
    ras_mt_ma_2013_3500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2012

ras_mt_ma_2012_3500m <- list()

for (i in 1:nrow(bf_mt_ma_2012_3500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2012_3500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2012_3500m <- crop(ma_2012, bf_mt_ma_2012_3500m[i, ])
  mask_mt_ma_2012_3500m <- mask(crop_mt_ma_2012_3500m, bf_mt_ma_2012_3500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2012_3500m[[buffer_id]] <- mask_mt_ma_2012_3500m
}

output_dir <- "ras_mt_ma_2012_3500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2012_3500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2012_3500m.tif"))
  writeRaster(
    ras_mt_ma_2012_3500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2011

ras_mt_ma_2011_3500m <- list()

for (i in 1:nrow(bf_mt_ma_2011_3500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2011_3500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2011_3500m <- crop(ma_2011, bf_mt_ma_2011_3500m[i, ])
  mask_mt_ma_2011_3500m <- mask(crop_mt_ma_2011_3500m, bf_mt_ma_2011_3500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2011_3500m[[buffer_id]] <- mask_mt_ma_2011_3500m
}

output_dir <- "ras_mt_ma_2011_3500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2011_3500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2011_3500m.tif"))
  writeRaster(
    ras_mt_ma_2011_3500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}

# 2010

ras_mt_ma_2010_3500m <- list()

for (i in 1:nrow(bf_mt_ma_2010_3500m)) {
  # Extrai o ID único do buffer
  buffer_id <- bf_mt_ma_2010_3500m$id_unico[i]
  
  # Recorta e mascara o raster
  crop_mt_ma_2010_3500m <- crop(ma_2010, bf_mt_ma_2010_3500m[i, ])
  mask_mt_ma_2010_3500m <- mask(crop_mt_ma_2010_3500m, bf_mt_ma_2010_3500m[i, ])
  
  # Armazena o raster na lista usando o ID como nome
  ras_mt_ma_2010_3500m[[buffer_id]] <- mask_mt_ma_2010_3500m
}

output_dir <- "ras_mt_ma_2010_3500m"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Salvando os recortes

for (nome in names(ras_mt_ma_2010_3500m)) {
  arquivo_saida <- file.path(output_dir, paste0(nome, "ras_mt_ma_2010_3500m.tif"))
  writeRaster(
    ras_mt_ma_2010_3500m[[nome]],
    filename = arquivo_saida,
    filetype = "GTiff",
    overwrite = TRUE)}


# Chamando os recortes

output_dir <- "ras_mt_ma_2023_3500m"
ras_mt_ma_2023_3500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2023_3500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2022_3500m"
ras_mt_ma_2022_3500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2022_3500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2021_3500m"
ras_mt_ma_2021_3500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2021_3500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2020_3500m"
ras_mt_ma_2020_3500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2020_3500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2019_3500m"
ras_mt_ma_2019_3500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2019_3500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2018_3500m"
ras_mt_ma_2018_3500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2018_3500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2017_3500m"
ras_mt_ma_2017_3500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2017_3500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2015_3500m"
ras_mt_ma_2015_3500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2015_3500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2014_3500m"
ras_mt_ma_2014_3500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2014_3500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2013_3500m"
ras_mt_ma_2013_3500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2013_3500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2012_3500m"
ras_mt_ma_2012_3500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2012_3500m\\.tif$")) %>%
  map(rast) 

output_dir <- "ras_mt_ma_2011_3500m"
ras_mt_ma_2011_3500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2011_3500m\\.tif$")) %>%
  map(rast) 



#### Metricas de paisagem 3500m ####

# 2023

id_unico <- names(ras_mt_ma_2023_3500m)

met_mt_ma_2023_3500m <- map_df(seq_along(ras_mt_ma_2023_3500m), function(i) {
  raster <- ras_mt_ma_2023_3500m[[i]]
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

prox_mt_ma_2023_3500m <- prox(ras_mt_ma_2023_3500m, 3)

shape_mt_ma_2023_3500m <- shape(ras_mt_ma_2023_3500m, class_value = 3)

lsm_mt_ma_2023_3500m <- met_mt_ma_2023_3500m |>
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
  inner_join(prox_mt_ma_2023_3500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2023_3500m, by = "id_unico")

# 2022

id_unico <- names(ras_mt_ma_2022_3500m)

met_mt_ma_2022_3500m <- map_df(seq_along(ras_mt_ma_2022_3500m), function(i) {
  raster <- ras_mt_ma_2022_3500m[[i]]
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

prox_mt_ma_2022_3500m <- prox(ras_mt_ma_2022_3500m, 3)

shape_mt_ma_2022_3500m <- shape(ras_mt_ma_2022_3500m, class_value = 3)

lsm_mt_ma_2022_3500m <- met_mt_ma_2022_3500m |>
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
  inner_join(prox_mt_ma_2022_3500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2022_3500m, by = "id_unico")

# 2021

id_unico <- names(ras_mt_ma_2021_3500m)

met_mt_ma_2021_3500m <- map_df(seq_along(ras_mt_ma_2021_3500m), function(i) {
  raster <- ras_mt_ma_2021_3500m[[i]]
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

prox_mt_ma_2021_3500m <- prox(ras_mt_ma_2021_3500m, 3)

shape_mt_ma_2021_3500m <- shape(ras_mt_ma_2021_3500m, class_value = 3)

lsm_mt_ma_2021_3500m <- met_mt_ma_2021_3500m |>
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
  inner_join(prox_mt_ma_2021_3500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2021_3500m, by = "id_unico")

# 2020

id_unico <- names(ras_mt_ma_2020_3500m)

met_mt_ma_2020_3500m <- map_df(seq_along(ras_mt_ma_2020_3500m), function(i) {
  raster <- ras_mt_ma_2020_3500m[[i]]
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

prox_mt_ma_2020_3500m <- prox(ras_mt_ma_2020_3500m, 3)

shape_mt_ma_2020_3500m <- shape(ras_mt_ma_2020_3500m, class_value = 3)

lsm_mt_ma_2020_3500m <- met_mt_ma_2020_3500m |>
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
  inner_join(prox_mt_ma_2020_3500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2020_3500m, by = "id_unico")

# 2019

id_unico <- names(ras_mt_ma_2019_3500m)

met_mt_ma_2019_3500m <- map_df(seq_along(ras_mt_ma_2019_3500m), function(i) {
  raster <- ras_mt_ma_2019_3500m[[i]]
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

prox_mt_ma_2019_3500m <- prox(ras_mt_ma_2019_3500m, 3)

shape_mt_ma_2019_3500m <- shape(ras_mt_ma_2019_3500m, class_value = 3)

lsm_mt_ma_2019_3500m <- met_mt_ma_2019_3500m |>
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
  inner_join(prox_mt_ma_2019_3500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2019_3500m, by = "id_unico")

# 2018

id_unico <- names(ras_mt_ma_2018_3500m)

met_mt_ma_2018_3500m <- map_df(seq_along(ras_mt_ma_2018_3500m), function(i) {
  raster <- ras_mt_ma_2018_3500m[[i]]
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

prox_mt_ma_2018_3500m <- prox(ras_mt_ma_2018_3500m, 3)

shape_mt_ma_2018_3500m <- shape(ras_mt_ma_2018_3500m, class_value = 3)

lsm_mt_ma_2018_3500m <- met_mt_ma_2018_3500m |>
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
  inner_join(prox_mt_ma_2018_3500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2018_3500m, by = "id_unico")

# 2017

id_unico <- names(ras_mt_ma_2017_3500m)

met_mt_ma_2017_3500m <- map_df(seq_along(ras_mt_ma_2017_3500m), function(i) {
  raster <- ras_mt_ma_2017_3500m[[i]]
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

prox_mt_ma_2017_3500m <- prox(ras_mt_ma_2017_3500m, 3)

shape_mt_ma_2017_3500m <- shape(ras_mt_ma_2017_3500m, class_value = 3)

lsm_mt_ma_2017_3500m <- met_mt_ma_2017_3500m |>
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
  inner_join(prox_mt_ma_2017_3500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2017_3500m, by = "id_unico")


# 2016

id_unico <- names(ras_mt_ma_2016_3500m)

met_mt_ma_2016_3500m <- map_df(seq_along(ras_mt_ma_2016_3500m), function(i) {
  raster <- ras_mt_ma_2016_3500m[[i]]
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

prox_mt_ma_2016_3500m <- prox(ras_mt_ma_2016_3500m, 3)

shape_mt_ma_2016_3500m <- shape(ras_mt_ma_2016_3500m, class_value = 3)

lsm_mt_ma_2016_3500m <- met_mt_ma_2016_3500m |>
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
  inner_join(prox_mt_ma_2016_3500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2016_3500m, by = "id_unico")

# 2015

id_unico <- names(ras_mt_ma_2015_3500m)

met_mt_ma_2015_3500m <- map_df(seq_along(ras_mt_ma_2015_3500m), function(i) {
  raster <- ras_mt_ma_2015_3500m[[i]]
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

prox_mt_ma_2015_3500m <- prox(ras_mt_ma_2015_3500m, 3)

shape_mt_ma_2015_3500m <- shape(ras_mt_ma_2015_3500m, class_value = 3)

lsm_mt_ma_2015_3500m <- met_mt_ma_2015_3500m |>
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
  inner_join(prox_mt_ma_2015_3500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2015_3500m, by = "id_unico")

# 2014

id_unico <- names(ras_mt_ma_2014_3500m)

met_mt_ma_2014_3500m <- map_df(seq_along(ras_mt_ma_2014_3500m), function(i) {
  raster <- ras_mt_ma_2014_3500m[[i]]
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

prox_mt_ma_2014_3500m <- prox(ras_mt_ma_2014_3500m, 3)

shape_mt_ma_2014_3500m <- shape(ras_mt_ma_2014_3500m, class_value = 3)

lsm_mt_ma_2014_3500m <- met_mt_ma_2014_3500m |>
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
  inner_join(prox_mt_ma_2014_3500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2014_3500m, by = "id_unico")

# 2013

id_unico <- names(ras_mt_ma_2013_3500m)

met_mt_ma_2013_3500m <- map_df(seq_along(ras_mt_ma_2013_3500m), function(i) {
  raster <- ras_mt_ma_2013_3500m[[i]]
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

prox_mt_ma_2013_3500m <- prox(ras_mt_ma_2013_3500m, 3)

shape_mt_ma_2013_3500m <- shape(ras_mt_ma_2013_3500m, class_value = 3)

lsm_mt_ma_2013_3500m <- met_mt_ma_2013_3500m |>
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
  inner_join(prox_mt_ma_2013_3500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2013_3500m, by = "id_unico")

# 2012
id_unico <- names(ras_mt_ma_2012_3500m)

met_mt_ma_2012_3500m <- map_df(seq_along(ras_mt_ma_2012_3500m), function(i) {
  raster <- ras_mt_ma_2012_3500m[[i]]
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

prox_mt_ma_2012_3500m <- prox(ras_mt_ma_2012_3500m, 3)

shape_mt_ma_2012_3500m <- shape(ras_mt_ma_2012_3500m, class_value = 3)

lsm_mt_ma_2012_3500m <- met_mt_ma_2012_3500m |>
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
  inner_join(prox_mt_ma_2012_3500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2012_3500m, by = "id_unico")

# 2011

id_unico <- names(ras_mt_ma_2011_3500m)

met_mt_ma_2011_3500m <- map_df(seq_along(ras_mt_ma_2011_3500m), function(i) {
  raster <- ras_mt_ma_2011_3500m[[i]]
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

prox_mt_ma_2011_3500m <- prox(ras_mt_ma_2011_3500m, 3)

shape_mt_ma_2011_3500m <- shape(ras_mt_ma_2011_3500m, class_value = 3)

lsm_mt_ma_2011_3500m <- met_mt_ma_2011_3500m |>
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
  inner_join(prox_mt_ma_2011_3500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2011_3500m, by = "id_unico")
# 2010

id_unico <- names(ras_mt_ma_2010_3500m)

met_mt_ma_2010_3500m <- map_df(seq_along(ras_mt_ma_2010_3500m), function(i) {
  raster <- ras_mt_ma_2010_3500m[[i]]
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

prox_mt_ma_2010_3500m <- prox(ras_mt_ma_2010_3500m, 3)

shape_mt_ma_2010_3500m <- shape(ras_mt_ma_2010_3500m, class_value = 3)

lsm_mt_ma_2010_3500m <- met_mt_ma_2010_3500m |>
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
  inner_join(prox_mt_ma_2010_3500m, by = "id_unico")|>
  inner_join(shape_mt_ma_2010_3500m, by = "id_unico")



#### Unificando tabelas 3500m ####

lsm_mt_ma_3500m <- bind_rows(lsm_mt_ma_2023_3500m, lsm_mt_ma_2022_3500m, lsm_mt_ma_2021_3500m, lsm_mt_ma_2020_3500m, lsm_mt_ma_2019_3500m, lsm_mt_ma_2018_3500m, lsm_mt_ma_2017_3500m, lsm_mt_ma_2016_3500m, lsm_mt_ma_2015_3500m, lsm_mt_ma_2014_3500m, lsm_mt_ma_2013_3500m, lsm_mt_ma_2012_3500m, lsm_mt_ma_2011_3500m, lsm_mt_ma_2010_3500m) |>
  dplyr::select(-pland_0,-lpi_0, -ed_0, -pd_0, -np_0, -lpi_12, -ed_12, -pd_12, -np_12, -lpi_15, -ed_15, -pd_15, -np_15, -lpi_9, -ed_9, -pd_9, -np_9, -np_24, -pd_24, -lpi_24, -ed_24, -pland_24)|>
  mutate(Bin = str_extract(id_unico, "^[01]"),
         Bin = as.factor(Bin),     
         Ano = as.factor(Ano))|>
  mutate(pland_9 = ifelse(is.na(pland_9), 0, pland_9))|>
  mutate(pland_15 = ifelse(is.na(pland_15), 0, pland_15))|>
  mutate(pland_12 = ifelse(is.na(pland_12), 0, pland_12)) |>
  rename_with(~ paste0(., "_3500m"))


lsm_mt_ma_3500m_sem_na <- lsm_mt_ma_3500m|>
  na.omit() 

write_xlsx(lsm_mt_ma_3500m_sem_na, "lsm_mt_ma_3500m.xlsx")


lsm_mt_ma_3500m <- read_excel("D:/Analise estatistica M. tridactyla - MA/Analise estatistica M. tridactyla - MA/lsm_mt_ma_3500m.xlsx", sheet = 1,na = c(" ", "NA"))


#### Distancia massa dagua ####

dist_agua_mt_ma <- st_distance(mt_ma_flo, agua_ma)

mt_ma_flo$dist_min_agua <- as.numeric(apply(dist_agua_mt_ma, 1, min))

dist_min_agua_mt_ma <- st_drop_geometry(mt_ma_flo[, c("id_unico", "dist_min_agua")])



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
distancias_result <- calcular_dist_minima(mt_ma_flo, UC)

# Juntar com dados originais
mt_ma_flo <- mt_ma_flo %>%
  left_join(distancias_result, by = "id_unico")

#### Selecionando rodovias com pavimentação ####

table(rod$revestimen)
table(rod$revestim_1)

rod <- as_tibble(rod)

rod_pav <- rod  |> 
  dplyr::filter(revestimen == "Pavimentado" | revestim_1 == "Pavimentado")

rod_pav <- st_as_sf(rod_pav)

#### Densidade de rodovias 500m ####

buffer_500m <- rbind(bf_mt_ma_2023_500m, bf_mt_ma_2022_500m, bf_mt_ma_2021_500m, bf_mt_ma_2020_500m, bf_mt_ma_2019_500m, bf_mt_ma_2018_500m, bf_mt_ma_2017_500m, bf_mt_ma_2016_500m, bf_mt_ma_2015_500m, bf_mt_ma_2014_500m, bf_mt_ma_2013_500m, bf_mt_ma_2012_500m, bf_mt_ma_2011_500m, bf_mt_ma_2010_500m)

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


den_500m <- data.frame(id_unico = lsm_mt_ma_500m_total$id_unico_500m) |>
  left_join(den_rod_500m, by = "id_unico") |>
  mutate(den_rod_500m = ifelse(is.na(den_rod_500m), 0, as.numeric(den_rod_500m)))

lsm_mt_ma_500m_total$den_rod_500m <- den_500m$den_rod_500m

den_rod_2023 <- st_drop_geometry(lsm_mt_ma_500m_total[, c("id_unico_500m", "den_rod_500m")])

#### Densidade de rodovias 1000m ####

buffer_1000m <- rbind(bf_mt_ma_2023_1000m, bf_mt_ma_2022_1000m, bf_mt_ma_2021_1000m, bf_mt_ma_2020_1000m, bf_mt_ma_2019_1000m, bf_mt_ma_2018_1000m, bf_mt_ma_2017_1000m, bf_mt_ma_2016_1000m, bf_mt_ma_2015_1000m, bf_mt_ma_2014_1000m, bf_mt_ma_2013_1000m, bf_mt_ma_2012_1000m, bf_mt_ma_2011_1000m, bf_mt_ma_2010_1000m)

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
         den_rod_1000m = (comprimento_rod_1000m/1000) / ((3.14*500^2)/1e6))


den_1000m <- data.frame(id_unico = lsm_mt_ma_1000m$id_unico_1000m) |>
  left_join(den_rod_1000m, by = "id_unico") |>
  mutate(den_rod_1000m = ifelse(is.na(den_rod_1000m), 0, as.numeric(den_rod_1000m)))

lsm_mt_ma_1000m$den_rod_1000m <- den_1000m$den_rod_1000m

den_rod_2023 <- st_drop_geometry(lsm_mt_ma_1000m[, c("id_unico_1000m", "den_rod_1000m")])

#### Densidade de rodovias 2000m ####

buffer_2000m <- rbind(bf_mt_ma_2023_2000m, bf_mt_ma_2022_2000m, bf_mt_ma_2021_2000m, bf_mt_ma_2020_2000m, bf_mt_ma_2019_2000m, bf_mt_ma_2018_2000m, bf_mt_ma_2017_2000m, bf_mt_ma_2016_2000m, bf_mt_ma_2015_2000m, bf_mt_ma_2014_2000m, bf_mt_ma_2013_2000m, bf_mt_ma_2012_2000m, bf_mt_ma_2011_2000m, bf_mt_ma_2010_2000m)

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
         den_rod_2000m = (comprimento_rod_2000m/1000) / ((3.14*500^2)/1e6))


den_2000m <- data.frame(id_unico = lsm_mt_ma_2000m$id_unico_2000m) |>
  left_join(den_rod_2000m, by = "id_unico") |>
  mutate(den_rod_2000m = ifelse(is.na(den_rod_2000m), 0, as.numeric(den_rod_2000m)))

lsm_mt_ma_2000m$den_rod_2000m <- den_2000m$den_rod_2000m

den_rod_2023 <- st_drop_geometry(lsm_mt_ma_2000m[, c("id_unico_2000m", "den_rod_2000m")])
#### Selecionando linhas ####

lsm_mt_ma_total <- read_excel("lsm_mt_ma_total.xlsx", sheet = 1,na = c(" ", "NA")) 

lsm_mt_ma_total <- lsm_mt_ma_total|>
  na.omit() 

View(lsm_mt_ma_total)

# Identifica as linhas com Bin == 0
zeros <- which(lsm_mt_ma_total$Bin == 0)

# Seleciona 7 zeros aleatórios para remover
remove_rows <- sample(zeros, min(15, length(zeros)))

# Filtra o tibble
lsm_mt_ma_total <- lsm_mt_ma_total[-remove_rows, ]

#### Estandartizar ####

# 1. Separar as colunas de identificação e as colunas numéricas
id_cols <- lsm_mt_ma_total |> select(id_unico, Ano, Bin)

num_cols <- lsm_mt_ma_total |> select(-id_unico, -Ano, -Bin)


# 2. Aplicar a padronização (exemplo com método "standardize" - scale para média 0 e sd 1)
stan <- vegan::decostand(num_cols, method = "standardize")

# 3. Juntar novamente com as colunas de identificação
lsm_mt_ma_total_std <- bind_cols(id_cols, stan)

# Verificar o resultado
glimpse(lsm_mt_ma_total_std)

lsm_mt_ma_total_std$Bin <- as.factor(lsm_mt_ma_total_std$Bin)


var_mt_ma <- lsm_mt_ma_total_std 

#### Correlacao ####

num_var_mt_ma <- var_mt_ma|>
  dplyr::select(-id_unico, - Bin, -Ano)

cor_mt_ma <- cor(num_var_mt_ma, method = "spearman", use = "pairwise.complete.obs")

cor_mt_ma

#### GLMM ####

var_mt_ma_total <- var_mt_ma|>
  dplyr::select(-id_unico)

var_resposta <- "Bin"
var_aleatoria <- "Ano"
vars_preditoras <- setdiff(names(var_mt_ma_total), 
                           c(var_resposta, var_aleatoria, "Ano"))

# 2. Definir pares de variáveis correlacionadas para excluir (exemplo com |r| > 0.7)
#Correlações positivas > 0.7
pares_positivos <- which(cor_mt_ma > 0.7 & upper.tri(cor_mt_ma), arr.ind = TRUE)

# Correlações negativas < -0.7
pares_negativos <- which(cor_mt_ma < -0.7 & upper.tri(cor_mt_ma), arr.ind = TRUE)

# Combinar ambos (equivalente ao seu código atual)
pares_excluir <- rbind(pares_positivos, pares_negativos)
pares_excluir <- data.frame(
  var1 = rownames(cor_mt_ma)[pares_excluir[, 1]],
  var2 = colnames(cor_mt_ma)[pares_excluir[, 2]],
  correlacao = cor_mt_ma[cbind(pares_excluir[, 1], pares_excluir[, 2])]  # Opcional: adiciona o valor da correlação
)

gerar_combinacoes_validas <- function(vars, pares_excluir, max_vars = 2) {
  # Garante que max_vars não seja maior que o número de variáveis disponíveis
  max_vars <- min(max_vars, length(vars))
  
  # Gera combinações apenas se houver variáveis suficientes
  if(length(vars) == 0) return(list())
  
  # Gerar todas combinações e aplicar filtro
  todas_combinacoes <- map(1:max_vars, ~combn(vars, .x, simplify = FALSE)) |> 
    flatten() |>
    keep(~{
      if(length(.x) < 2) TRUE else {
        todos_pares <- combn(.x, 2, simplify = FALSE)
        !any(map_lgl(todos_pares, ~{
          any((pares_excluir$var1 %in% .x & pares_excluir$var2 %in% .x) |
                (pares_excluir$var2 %in% .x & pares_excluir$var1 %in% .x))
        }))
      }
    })
  
  return(todas_combinacoes)
}


# 4. Gerar todas as combinações válidas
combinacoes_validas <- gerar_combinacoes_validas(vars_preditoras, pares_excluir)

# 5. Função para ajustar modelos
ajustar_modelo_modificado <- function(vars, idx = NULL) {
  if(length(vars) == 0) return(NULL)
  
  # Mensagem informando qual modelo está sendo ajustado
  if(!is.null(idx)) {
    message("\nAjustando modelo ", idx, " de ", length(combinacoes_validas), 
            " com variáveis: ", paste(vars, collapse = ", "))
  } else {
    message("\nAjustando modelo com variáveis: ", paste(vars, collapse = ", "))
  }
  
  formula <- as.formula(paste(var_resposta, "~", paste(vars, collapse = "+")))
  
  tryCatch({
    modelo <- glmmTMB(formula, 
                      data = var_mt_ma_total, 
                      family = binomial())
    
    # Mostra um resumo do modelo ajustado
    message("\nModelo ajustado com sucesso!")
    print(summary(modelo))  # Mostra um resumo do modelo
    
    return(modelo)
  }, error = function(e) {
    message("\nErro ao ajustar modelo com variáveis: ", paste(vars, collapse = ", "))
    message("Mensagem de erro: ", e$message)
    return(NULL)
  })
}


# 6. Ajustar todos os modelos (pode demorar)
modelos <- map(combinacoes_validas, ajustar_modelo_modificado)


#### Selecao dos modelos ####

# Calcular AICc e pesos
aicc_tabela <- MuMIn::model.sel(modelos)
weights <- MuMIn::Weights(aicc_tabela)

# Melhor modelo
melhor_modelo_1 <- MuMIn::get.models(aicc_tabela, , subset = cumsum(weights) <= 0.95)[[1]]
melhor_modelo_1
summary(melhor_modelo_1)
simulationOutput <- DHARMa::simulateResiduals(fittedModel = melhor_modelo_1, plot = TRUE)

melhor_modelo_2 <- MuMIn::get.models(aicc_tabela, , subset = cumsum(weights) <= 0.95)[[2]]
melhor_modelo_2
summary(melhor_modelo_2)
simulationOutput <- DHARMa::simulateResiduals(fittedModel = melhor_modelo_2, plot = TRUE)


calculate_auc_for_selected <- function(model) {
  # Obter as predições (ignorando efeitos aleatórios para avaliação geral)
  preds <- predict(model, type = "response", re.form = NA)
  
  # Calcula AUC
  roc_obj <- roc(response = var_mt_ma_total[[var_resposta]], 
                 predictor = preds)
  
  # Retornar AUC e objeto ROC para plotagem
  return(list(AUC = round(auc(roc_obj), 3),
              ROC = roc_obj))
}

## Calcular AUC para o melhor modelo ##
auc_1 <- calculate_auc_for_selected(melhor_modelo_1) # 0.711
auc_2 <- calculate_auc_for_selected(melhor_modelo_2) # 0.624

#### Resultados ####


#> summary(melhor_modelo_1)
#Family: binomial  ( logit )
#Formula:          Bin ~ np_4_3500
#Data: var_mt_ma_total

#AIC       BIC    logLik -2*log(L)  df.resid 
#170.2     175.8     -83.1     166.2       122 


#Conditional model:
#  Estimate Std. Error z value Pr(>|z|)  
#(Intercept) -0.01264    0.18433  -0.069   0.9453  
#np_4_3500   -0.46727    0.21018  -2.223   0.0262 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#> summary(melhor_modelo_2)
#Family: binomial  ( logit )
#Formula:          Bin ~ pd_4_3500
#Data: var_mt_ma_total

#AIC       BIC    logLik -2*log(L)  df.resid 
#170.3     175.9     -83.1     166.3       122 


#Conditional model:
#  Estimate Std. Error z value Pr(>|z|)  
#(Intercept) -0.01232    0.18424  -0.067   0.9467  
#pd_4_3500   -0.46270    0.20966  -2.207   0.0273 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

resultados <- as_tibble(aicc_tabela)
resultados
write_xlsx(resultados, "resultados_mt_ma.xlsx")
