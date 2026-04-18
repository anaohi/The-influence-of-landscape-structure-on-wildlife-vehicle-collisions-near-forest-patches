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
library(MuMIn)

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
  st_make_valid()|>
  st_transform(5880)

ma_2023 <- rast("D:/Analise estatistica mestrado/proj_ma_2023.tif")
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
ma_2009 <- rast("D:/Analise estatistica mestrado/proj_ma_2009.tif")
ma_2008 <- rast("D:/Analise estatistica mestrado/proj_ma_2008.tif")
ma_2007 <- rast("D:/Analise estatistica mestrado/proj_ma_2007.tif")
ma_2006 <- rast("D:/Analise estatistica mestrado/proj_ma_2006.tif")
ma_2005 <- rast("D:/Analise estatistica T. tetradactyla - Cerrado/Analise estatistica T. tetradactyla - Cerrado/proj_br_2005.tif")
ma_2004 <- rast("D:/Analise estatistica T. tetradactyla - Cerrado/Analise estatistica T. tetradactyla - Cerrado/proj_br_2004.tif")
ma_2003 <- rast("D:/Analise estatistica mestrado/proj_ma_2003.tif")
ma_2002 <- rast("D:/Analise estatistica mestrado/proj_ma_2002.tif")
ma_2001 <- rast("D:/Analise estatistica mestrado/proj_ma_2001.tif")
ma_2000 <- rast("D:/Analise estatistica mestrado/proj_ma_2000.tif")
ma_1999 <- rast("D:/Analise estatistica mestrado/proj_ma_1999.tif")


#### Verificando duplicatas de pseudo-ausencia ####

dup_mt_ma <- mt_ma |>
  group_by(geometry) |>
  filter(n() > 1) 
dup_mt_ma


mt_ma <- select(mt_ma, c(-Mes))

#### Inserindo ID nas ocorrencias ####

mt_ma$id_unico <- ave(
  mt_ma$id,             
  mt_ma$id,             
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

#### Criando a funcao do indide de forma ####

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

#### Filtando pontos proximos as coberturas florestais ####

mt_ma_flo <- inner_join(lsm_mt_ma_500m_total, mt_ma, by = c("id_unico_500m" = "id_unico"))|>
  select(id, id_unico_500m, Ano, geometry)|>
  rename("id_unico" = "id_unico_500m")|>
  st_as_sf(crs = 5880)

write_xlsx(mt_ma_flo, "mt_ma_flo.xlsx")

st_write(st_as_sf(mt_ma_flo), "mt_ma_flo_2.shp", append=FALSE)

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

output_dir <- "ras_mt_ma_2010_2000m"
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

write_xlsx(lsm_mt_ma_2000m_sem_na, "lsm_mt_ma_2000m.xlsx")

View(lsm_mt_ma_2000m_sem_na)



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

output_dir <- "ras_mt_ma_2010_3500m"
ras_mt_ma_2010_3500m <- list.files(output_dir, pattern = "\\.tif$", full.names = TRUE) %>%
  set_names(~ str_remove(basename(.), "ras_mt_ma_2010_3500m\\.tif$")) %>%
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

write_xlsx(lsm_mt_ma_3500m, "lsm_mt_ma_3500m.xlsx")

lsm_mt_ma_3500m_sem_na <- lsm_mt_ma_3500m|>
  na.omit() 

View(lsm_mt_ma_3500m_sem_na)


#### Distancia massa dagua ####

dist_agua_mt_ma_2023 <- st_distance(mt_ma_flo, agua_ma)

mt_ma_flo$dist_min_agua <- as.numeric(apply(dist_agua_mt_ma_2023, 1, min))

dist_min_agua_ma_mt_2023 <- st_drop_geometry(mt_ma_flo[, c("id_unico", "dist_min_agua")])



#### Distancia UCs ####

UC$cria_ano <- as.numeric(format(as.Date(UC$cria_ano, format = "%d-%m-%Y"), "%Y"))

# 2023

uc_2023 <- UC|>
  dplyr::filter(cria_ano <= 2023)
mt_ma_2023 <- mt_ma_flo|>
  dplyr::filter(Ano <= 2023)
dist_uc_mt_ma_2023 <- st_distance(mt_ma_2023, uc_2023)

mt_ma_2023$dist_min_uc <- as.numeric(apply(dist_uc_mt_ma_2023, 1, min))

dist_min_uc_2023 <- st_drop_geometry(mt_ma_2023[, c("id_unico", "dist_min_uc")])



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

#### Densidade de rodovias 3500m ####

buffer_3500m <- rbind(bf_mt_ma_2023_3500m, bf_mt_ma_2022_3500m, bf_mt_ma_2021_3500m, bf_mt_ma_2020_3500m, bf_mt_ma_2019_3500m, bf_mt_ma_2018_3500m, bf_mt_ma_2017_3500m, bf_mt_ma_2016_3500m, bf_mt_ma_2015_3500m, bf_mt_ma_2014_3500m, bf_mt_ma_2013_3500m, bf_mt_ma_2012_3500m, bf_mt_ma_2011_3500m, bf_mt_ma_2010_3500m)

rod_buffer_3500m <- st_intersection(rod_pav, buffer_3500m)

rod_dis_3500m <- aggregate(
  x = rod_buffer_3500m ["geometry"],
  by = list(id_unico = rod_buffer_3500m$id_unico),
  FUN = function(x) st_union(x)) 

comp_rod_3500m <- rod_dis_3500m |>
  mutate(comp_rod_3500m = st_length(geometry))

# Juntar as informações de comprimento com a área do buffer
den_rod_3500m <- buffer_3500m |> 
  left_join(comp_rod_3500m |> 
              st_drop_geometry() |> 
              mutate(comprimento_rod_3500m = as.numeric(comp_rod_3500m)),
            by = "id_unico")|>
  mutate(comprimento_rod_3500m = ifelse(is.na(comprimento_rod_3500m), 0, comprimento_rod_3500m),
         den_rod_3500m = (comprimento_rod_3500m/3500) / ((3.14*500^2)/1e6))


den_3500m <- data.frame(id_unico = lsm_mt_ma_3500m$id_unico_3500m) |>
  left_join(den_rod_3500m, by = "id_unico") |>
  mutate(den_rod_3500m = ifelse(is.na(den_rod_3500m), 0, as.numeric(den_rod_3500m)))

lsm_mt_ma_3500m$den_rod_3500m <- den_3500m$den_rod_3500m

den_rod_2023 <- st_drop_geometry(lsm_mt_ma_3500m[, c("id_unico_3500m", "den_rod_3500m")])

#### Unificando todas as tabelas ####

lsm_mt_ma_500m_uni <- lsm_mt_ma_500m_total |>
  dplyr::rename('Bin'='Bin_500m')|>
  dplyr::rename('Ano'='Ano_500m')|>
  dplyr::rename('id_unico'='id_unico_500m')

lsm_mt_ma_2000m_uni <- lsm_mt_ma_2000m |>
  dplyr::rename('id_unico'='id_unico_2000m')|>
  select(-Bin_2000m, -Ano_2000m)

lsm_mt_ma_3500m_uni <- lsm_mt_ma_3500m |>
  dplyr::rename('id_unico'='id_unico_3500m')|>
  select(-Bin_3500m, -Ano_3500m)

glimpse(c(lsm_mt_ma_500m_uni, lsm_mt_ma_2000m_uni, lsm_mt_ma_3500m_uni))

dist_min_uc_mt_ma <- select(distancias_result, c("id_unico", "dist_min_UC"))

lsm_mt_ma_total <- lsm_mt_ma_500m_uni %>%
  full_join(lsm_mt_ma_2000m_uni, by = 'id_unico') %>%
  full_join(lsm_mt_ma_3500m_uni, by = 'id_unico') %>%
  full_join(dist_min_agua_ma_mt_2023, by = 'id_unico') %>%
  full_join(dist_min_uc_mt_ma, by = 'id_unico')|>
  na.omit()

glimpse(lsm_mt_ma_total)

write_xlsx(lsm_mt_ma_total, "lsm_mt_ma_total.xlsx")

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

write_xlsx(lsm_mt_ma_total_std, "lsm_mt_ma_total_std.xlsx")

lsm_mt_ma_total_std$Bin <- as.factor(lsm_mt_ma_total_std$Bin)


var_mt_ma <- lsm_mt_ma_total_std 

#### Correlacao ####

num_var_mt_ma <- var_mt_ma|>
  dplyr::select(-id_unico, - Bin, -Ano)

cor_mt_ma <- cor(num_var_mt_ma, method = "spearman")

glimpse(cor_mt_ma)


# Remove combinacoes duplicadas e mesmas metricas com escalas diferentes
cor_sem_dup_mt_ma <- as.data.frame(cor_mt_ma) %>%
  mutate(var1 = rownames(cor_mt_ma)) %>%
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


cor_sem_dup_mt_ma


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
comb <- var_glm(cor_sem_dup_mt_ma)

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
                         data = var_mt_ma, 
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


aic_tabela <- MuMIn::model.sel(modelos_vif, rank = AIC)
options(max.print = 10056)
aic_tabela


melhor_modelo_1 <- MuMIn::get.models(aic_tabela, subset = 1)[[1]]
summary(melhor_modelo_1)
simulationOutput <- DHARMa::simulateResiduals(fittedModel = melhor_modelo_1, plot = TRUE)

install.packages("car")

library(car)

vif(melhor_modelo_1)

melhor_modelo_2 <- MuMIn::get.models(aic_tabela, subset = 2)[[1]]
summary(melhor_modelo_2)
simulationOutput <- DHARMa::simulateResiduals(fittedModel = melhor_modelo_2, plot = TRUE)

vif(melhor_modelo_2)


melhor_modelo_3 <- MuMIn::get.models(aic_tabela, subset = 3)[[1]]
summary(melhor_modelo_3)
simulationOutput <- DHARMa::simulateResiduals(fittedModel = melhor_modelo_3, plot = TRUE)

melhor_modelo_4 <- MuMIn::get.models(aic_tabela, subset = 4)[[1]]
summary(melhor_modelo_4)
simulationOutput <- DHARMa::simulateResiduals(fittedModel = melhor_modelo_4, plot = TRUE)
