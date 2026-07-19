#### Pacote ####  

library(terra)


#### Reclassificando raster ####

ma_2023 <- rast("D:/CETESB/Rasters MapBiomas/ma_2023.tif") # WGS 84
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
ma_2012 <- rast("D:/CETESB/Rasters MapBiomas/ma_12.tif")
ma_2011 <- rast("D:/CETESB/Rasters MapBiomas/ma_11.tif")
ma_2010 <- rast("D:/CETESB/Rasters MapBiomas/ma_10.tif")
ma_2009 <- rast("D:/CETESB/Rasters MapBiomas/ma_09.tif")
ma_2008 <- rast("D:/CETESB/Rasters MapBiomas/ma_08.tif")
ma_2007 <- rast("D:/CETESB/Rasters MapBiomas/ma_07.tif")
ma_2006 <- rast("D:/CETESB/Rasters MapBiomas/ma_06.tif")
ma_2005 <- rast("ma_final_2005.tif")
ma_2004 <- rast("proj_rec_ma_2004.tif")
ma_2003 <- rast("D:/CETESB/Rasters MapBiomas/ma_03.tif")
ma_2002 <- rast("D:/CETESB/Rasters MapBiomas/ma_02.tif")
ma_2001 <- rast("D:/CETESB/Rasters MapBiomas/ma_01.tif")
ma_2000 <- rast("D:/CETESB/Rasters MapBiomas/ma_00.tif")
ma_1999 <- rast("D:/CETESB/Rasters MapBiomas/ma_99.tif")

matriz_reclas <- c(3, 3, # Formacao florestal
                   4, 3, # Formacao savanica (floresta)
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
                   50, 0,
                   35, 0)


# Converta para matriz
matriz_reclas <- matrix(matriz_reclas, ncol=2, byrow=TRUE)

# Aplicando a reclassificacao

ma_2023 <- classify(ma_2023, matriz_reclas)
writeRaster(ma_2023, filename="rec_ma_2023.tif", overwrite=TRUE)

ma_2022 <- classify(ma_2022, matriz_reclas)
writeRaster(ma_2022, filename="rec_ma_2022.tif", overwrite=TRUE)

ma_2021 <- classify(ma_2021, matriz_reclas)
writeRaster(ma_2021, filename="rec_ma_2021.tif", overwrite=TRUE)

ma_2020 <- classify(ma_2020, matriz_reclas)
writeRaster(ma_2020, filename="rec_ma_2020.tif", overwrite=TRUE)

ma_2019 <- classify(ma_2019, matriz_reclas)
writeRaster(ma_2019, filename="rec_ma_2019.tif", overwrite=TRUE)

ma_2018 <- classify(ma_2018, matriz_reclas)
writeRaster(ma_2018, filename="rec_ma_2018.tif", overwrite=TRUE)

ma_2017 <- classify(ma_2017, matriz_reclas)
writeRaster(ma_2017, filename="rec_ma_2017.tif", overwrite=TRUE)

ma_2016 <- classify(ma_2016, matriz_reclas)
writeRaster(ma_2016, filename="rec_ma_2016.tif", overwrite=TRUE)

ma_2015 <- classify(ma_2015, matriz_reclas)
writeRaster(ma_2015, filename="rec_ma_2015.tif", overwrite=TRUE)

ma_2014 <- classify(ma_2014, matriz_reclas)
writeRaster(ma_2014, filename="rec_ma_2014.tif", overwrite=TRUE)

ma_2013 <- classify(ma_2013, matriz_reclas)
writeRaster(ma_2013, filename="rec_ma_2013.tif", overwrite=TRUE)

ma_2012 <- classify(ma_2012, matriz_reclas)
writeRaster(ma_2012, filename="rec_ma_2012.tif", overwrite=TRUE)

ma_2011 <- classify(ma_2011, matriz_reclas)
writeRaster(ma_2011, filename="rec_ma_2011.tif", overwrite=TRUE)

ma_2010 <- classify(ma_2010, matriz_reclas)
writeRaster(ma_2010, filename="rec_ma_2010.tif", overwrite=TRUE)

ma_2009 <- classify(ma_2009, matriz_reclas)
writeRaster(ma_2009, filename="rec_ma_2009.tif", overwrite=TRUE)

ma_2008 <- classify(ma_2008, matriz_reclas)
writeRaster(ma_2008, filename="rec_ma_2008.tif", overwrite=TRUE)

ma_2007 <- classify(ma_2007, matriz_reclas)
writeRaster(ma_2007, filename="rec_ma_2007.tif", overwrite=TRUE)

ma_2006 <- classify(ma_2006, matriz_reclas)
writeRaster(ma_2006, filename="rec_ma_2006.tif", overwrite=TRUE)

ma_2005 <- classify(ma_2005, matriz_reclas)
writeRaster(ma_2005, filename="ma_final_2005.tif", overwrite=TRUE)

ma_2004 <- classify(ma_2004, matriz_reclas) 
writeRaster(ma_2004, filename="ma_final_2004.tif", overwrite=TRUE)

ma_2003 <- classify(ma_2003, matriz_reclas)
writeRaster(ma_2003, filename="rec_ma_2003.tif", overwrite=TRUE)

ma_2002 <- classify(ma_2002, matriz_reclas)
writeRaster(ma_2002, filename="rec_ma_2002.tif", overwrite=TRUE)

ma_2001 <- classify(ma_2001, matriz_reclas)
writeRaster(ma_2001, filename="rec_ma_2001.tif", overwrite=TRUE)

ma_2000 <- classify(ma_2000, matriz_reclas)
writeRaster(ma_2000, filename="rec_ma_2000.tif", overwrite=TRUE)

ma_1999 <- classify(ma_1999, matriz_reclas)
writeRaster(ma_1999, filename="rec_ma_1999.tif", overwrite=TRUE)