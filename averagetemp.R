library(ncdf)
library(lubridate)
library(fields)
library(animation)
lats = seq(89.5,-89.5,-1) 
lons = seq(0.5,359.5,1) 

lon1 = 5
lon2 = 355
lat1 = 89
lat2 = -89

lon1a = which.min(abs(lon1 - lons)) 
lon2a = which.min(abs(lon2 - lons)) 
lat1a = which.min(abs(lat1 - lats)) 
lat2a = which.min(abs(lat2 - lats)) 

nlon = (lon2a - lon1a) + 1 
nlat = (lat2a - lat1a) + 1

nc = open.ncdf("sst.wkmean.1990-present.nc")
ncdates = nc$dim$time$vals
ncdates = as.Date(ncdates,origin = '1800-1-1') 

saveGIF({
        for (i in 1:23){
                date1 = '1990-1-1'
                date2 = '1990-12-31'
                date1 = as.Date(date1, format = "%Y-%m-%d")+years(i) 
                date2 = as.Date(date2, format = "%Y-%m-%d")+years(i)
                date1a = which.min(abs(date1 - ncdates))
                date2a = which.min(abs(date2 - ncdates))
                ndates = (date2a - date1a) + 1
                sstout = get.var.ncdf(nc, varid = 'sst', start = c(lon1a,lat1a,date1a),count = c(nlon,nlat,ndates))
                datesout = ncdates[date1a:date2a]
                nc2 = open.ncdf('lsmask.nc') #open land-sea mask
                mask = get.var.ncdf(nc2, varid = "mask",start = c(lon1a,lat1a,1),count = c(nlon,nlat,1))                 
                mask = ifelse(mask == 0,NA,1) 
                dims = dim(sstout)
                        for (i in 1:dims[3]){
                                sstout[,,i] = sstout[,,i] * mask
                        }
                # get the mean
                temp <- apply(sstout, c(1,2), mean)
                goodlons = lons[lon1a:lon2a]
                goodlats = lats[lat1a:lat2a]
                par(mar = c(5,5,4,6)) #widen margins
                image(goodlons,rev(goodlats),temp,
                      ylim = rev(range(goodlats)), #reverse y-axis
                      col = tim.colors(32),
                      yaxt = "n", xaxt = "n",
                      ylab = "Latitude (degrees north)",
                      xlab = "Longitude (degrees east)",
                      main = paste("Year of ",year(date1),sep = ""))
                axis(1,at = pretty(goodlons),labels = pretty(goodlons), las = 1)
                axis(2,at = rev(pretty(goodlats)),labels = pretty(goodlats), las = 1)
                image.plot(zlim = range(temp,na.rm = TRUE),nlevel = 32,
                           legend.only = TRUE, horizontal = FALSE, col = tim.colors(32), legend.args = list(text = "Temperature,\u00B0C", cex = 1.2, side = 4, line = 2))
        }              
}, movie.name = "oceantemp.gif", ani.width = 800, ani.height = 500)
