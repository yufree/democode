library(tidyverse)
library(readxl)
library(ggplot2)
library(lubridate)
library(showtext)
library(Cairo)
pm2.5 <- read_xls('data/pm.xls')
colnames(pm2.5) <- c('time','PM2.5')
# extract data

dataset <- read_csv("data/pm2.5_13-17.csv",
                    col_types = cols(DateTime = col_datetime(format = "%Y-%m-%d %H:%M")))



# Make a dataframe
dat<-dataset

# We will facet by year ~ month, and each subgraph will
# show week-of-month versus weekday
# the year is simple
dat$year<-as.numeric(as.POSIXlt(dat$DateTime)$year+1900)
# the month too 
dat$month<-as.numeric(as.POSIXlt(dat$DateTime)$mon+1)
# but turn months into ordered facors to control the appearance/ordering in the presentation
dat$monthf<-factor(dat$month,levels=as.character(1:12),labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),ordered=TRUE)
# the day of week is again easily found
dat$weekday = as.POSIXlt(dat$DateTime)$wday
# again turn into factors to control appearance/abbreviation and ordering
# I use the reverse function rev here to order the week top down in the graph
# you can cut it out to reverse week order
dat$weekdayf<-factor(dat$weekday,levels=rev(1:7),labels=rev(c("Mon","Tue","Wed","Thu","Fri","Sat","Sun")),ordered=TRUE)
# the monthweek part is a bit trickier 
# first a factor which cuts the data into month chunks
dat$yearmonth<-as.yearmon(dat$DateTime)
dat$yearmonthf<-factor(dat$yearmonth)
# then find the "week of year" for each day
dat$week <- as.numeric(format(dat$DateTime,"%W"))
# and now for each monthblock we normalize the week to start at 1 
dat<-ddply(dat,.(yearmonthf),transform,monthweek=1+week-min(week))

# Now for the plot
P<- ggplot(dat, aes(monthweek, weekdayf, fill = `ug.m3`)) + 
        geom_tile(colour = "white") + facet_grid(year~monthf) + scale_fill_gradient(low="red", high="yellow") +
        xlab("Week of Month") + ylab("")


dataset <- dataset %>%
        mutate(hour = hour(DateTime),day = day(DateTime),month = month(DateTime),year = year(DateTime),week = weekdays.Date(DateTime), weeks = week(DateTime)) 

dataset$week2 <- factor(dataset$week, levels = c("Monday","Tuesday",'Wednesday','Thursday','Friday','Saturday','Sunday'))

library(ggridges)
CairoPDF("year.pdf", width = 9, height = 4)
showtext.auto()
font.add("test","Kaiti.ttf")
dataset %>%
        ggplot(aes(`ug/m3`, year, group = year)) + geom_density_ridges2(scale = 0.95)+xlab(expression(paste('PM'[2.5],'(',mu,"g/m"^{3},')')))+scale_y_continuous('年',breaks = c(2013:2017), labels = c('2013','2014','2015','2016','2017'))
dev.off()

CairoPDF("month.pdf", width = 9, height = 4)
showtext.auto()
font.add("test","Kaiti.ttf")
dataset %>%
        ggplot(aes(`ug/m3`, month, group = month)) + geom_density_ridges2(scale = 0.95)+facet_wrap(~year,nrow = 1)+xlab(expression(paste('PM'[2.5],'(',mu,"g/m"^{3},')')))+scale_y_continuous('月',breaks = c(1:12), labels = c('一月','二月','三月','四月','五月','六月','七月','八月','九月','十月','十一月','十二月'))
dev.off()

CairoPDF("week.pdf", width = 9, height = 4)
showtext.auto()
font.add("test","Kaiti.ttf")
dataset %>%
        ggplot(aes(`ug/m3`, weeks, group = weeks)) + geom_density_ridges2(scale = 0.95)+facet_wrap(~year,nrow = 1)+xlab(expression(paste('PM'[2.5],'(',mu,"g/m"^{3},')')))+scale_y_continuous('周')
dev.off()

CairoPDF("monthday.pdf", width = 9, height = 4)
showtext.auto()
font.add("test","Kaiti.ttf")
dataset %>%
        ggplot(aes(`ug/m3`, day, group = day)) + geom_density_ridges2(scale = 0.95)+facet_wrap(~year,nrow = 1)+xlab(expression(paste('PM'[2.5],'(',mu,"g/m"^{3},')')))+scale_y_continuous('月内变化',breaks = c(1:31), labels = c(1:31))
dev.off()

CairoPDF("weekday.pdf", width = 9, height = 4)
showtext.auto()
font.add("test","Kaiti.ttf")
dataset %>%
        ggplot(aes(`ug/m3`, week2, group = week2)) + geom_density_ridges2(scale = 0.95)+facet_wrap(~year,nrow = 1)+xlab(expression(paste('PM'[2.5],'(',mu,"g/m"^{3},')')))+scale_y_discrete('周内变化',labels = c('星期一','星期二','星期三','星期四','星期五','星期六','星期日'))
dev.off()

CairoPDF("daily.pdf", width = 9, height = 4)
showtext.auto()
font.add("test","Kaiti.ttf")
dataset %>%
        ggplot(aes(`ug/m3`, hour, group = hour)) + geom_density_ridges2(scale = 0.95)+facet_wrap(~year,nrow = 1) + xlab(expression(paste('PM'[2.5],'(',mu,"g/m"^{3},')')))+scale_y_continuous('日内变化',breaks = c(0:23), labels = c('0:00','1:00','2:00','3:00','4:00','5:00','6:00','7:00','8:00','9:00','10:00','11:00','12:00','13:00','14:00','15:00','16:00','17:00','18:00','19:00','20:00','21:00','22:00','23:00'))
dev.off()

CairoPDF("dailychange.pdf", width = 9, height = 4)
showtext.auto()
font.add("test","Kaiti.ttf")
dataset %>%
        ggplot(aes(hour,week2, height =`ug/m3`, group = week2)) + geom_density_ridges(stat = 'identity',scale = 0.95)+facet_wrap(~year,nrow = 1) + xlab(expression(paste('PM'[2.5],'(',mu,"g/m"^{3},')')))+scale_x_continuous('日内变化')+scale_y_discrete('周内变化',labels = c('星期一','星期二','星期三','星期四','星期五','星期六','星期日'))
dev.off()

CairoPDF("monthchange.pdf", width = 9, height = 4)
showtext.auto()
font.add("test","Kaiti.ttf")
dataset %>%
        ggplot(aes(hour,month, height =`ug/m3`, group = month)) + geom_density_ridges(stat = 'identity',scale = 0.95)+facet_wrap(~year,nrow = 1) + xlab(expression(paste('PM'[2.5],'(',mu,"g/m"^{3},')')))+scale_x_continuous('日内变化')+scale_y_continuous('月',breaks = c(1:12), labels = c('一月','二月','三月','四月','五月','六月','七月','八月','九月','十月','十一月','十二月'))
dev.off()

CairoPDF("pm2.pdf", width = 8, height = 5);
showtext.auto()
font.add("test","Kaiti.ttf")
ggplot(df.tmp[-c(1:24,193),], aes(hour, day, fill = PM2.5)) + 
        geom_tile(colour = "white") + 
        scale_fill_gradient(name = expression(paste("PM"[2.5], "(ppm)"))) +
        theme(panel.background = element_blank(),text = element_text(family = 'test',size=8), axis.text.y = element_text(size = 20))+
        scale_x_continuous('',breaks = c(0:23), labels = c('0:00','1:00','2:00','3:00','4:00','5:00','6:00','7:00','8:00','9:00','10:00','11:00','12:00','13:00','14:00','15:00','16:00','17:00','18:00','19:00','20:00','21:00','22:00','23:00'))+
scale_y_continuous('',breaks = c(13:19), labels=c('星期一','星期二','星期三','星期四','星期五','星期六','星期日'))
dev.off()

ggplot(df.tmp[-c(1:24,193),], aes(x = factor(hour), y = day)) + 
        geom_bar(aes(fill = day), stat = "identity", position = "stack") + 
        scale_fill_manual(values = c(blank = "white", dark = "black")) + 
        theme_minimal() + 
        guides(fill = FALSE)

randomData = expand.grid(1:5, 1:100)    # The two discrete variables
colnames(randomData) = c("index1", "index2")

randomData$z = runif(nrow(randomData)) 
d <- rbind(
        mutate(randomData, col = "dark"), 
        mutate(randomData, col = "blank", z = 1 - z)
) %>% arrange(d, index2, index1, col)


calendarHeat <- function(dates, 
                         values,
                         colors,
                         ncolors=99,
                         title,
                         date.form = "%Y-%m-%d", ...) {
        require(lattice, quietly = TRUE)
        require(grid, quietly = TRUE)
        
        if (class(dates) == "character" | class(dates) == "factor" ) {
                dates <- strptime(dates, date.form)
        }
        caldat <- data.frame(value = values, dates = dates)
        min.date <- as.Date(paste(format(min(dates), "%Y"),
                                  "-1-1",sep = ""))
        max.date <- as.Date(paste(format(max(dates), "%Y"),
                                  "-12-31", sep = ""))
        dates.f <- data.frame(date.seq = seq(min.date, max.date, by="days"))
        
        # Merge moves data by one day, avoid
        caldat <- data.frame(date.seq = seq(min.date, max.date, by="days"), value = NA)
        dates <- as.Date(dates) 
        caldat$value[match(dates, caldat$date.seq)] <- values
        
        caldat$dotw <- as.numeric(format(caldat$date.seq, "%w"))
        caldat$woty <- as.numeric(format(caldat$date.seq, "%U")) + 1
        caldat$yr <- as.factor(format(caldat$date.seq, "%Y"))
        caldat$month <- as.numeric(format(caldat$date.seq, "%m"))
        yrs <- as.character(unique(caldat$yr))
        d.loc <- as.numeric()                        
        for (m in min(yrs):max(yrs)) {
                d.subset <- which(caldat$yr == m)  
                sub.seq <- seq(1,length(d.subset))
                d.loc <- c(d.loc, sub.seq)
        }  
        caldat <- cbind(caldat, seq=d.loc)
        
        #color styles
        if (missing(colors)) colors <- c("#FFFFBD","#D61818","#711603")                          
        
        calendar.pal <- colorRampPalette(colors, space = "Lab")
        def.theme <- lattice.getOption("default.theme")
        cal.theme <-
                function() {  
                        theme <-
                                list(
                                        strip.background = list(col = "transparent"),
                                        strip.border = list(col = "transparent"),
                                        axis.line = list(col="transparent"),
                                        par.strip.text=list(cex=0.8))
                }
        lattice.options(default.theme = cal.theme)
        yrs <- (unique(caldat$yr))
        nyr <- length(yrs)
        print(cal.plot <- levelplot(value~woty*dotw | yr, data=caldat,
                                    as.table=TRUE,
                                    aspect=.12,
                                    layout = c(1, nyr%%7),
                                    between = list(x=0, y=c(0.5,0.5)),
                                    strip=TRUE,
                                    main = ifelse(missing(title), "", title),
                                    scales = list(
                                            x = list(
                                                    at= c(seq(2.9, 52, by=4.42)),
                                                    labels = c('一月','二月','三月','四月','五月','六月','七月','八月','九月','十月','十一月','十二月'),
                                                    alternating = c(1, rep(0, (nyr-1))),
                                                    tck=0,
                                                    cex = 0.7),
                                            y=list(
                                                    at = c(0, 1, 2, 3, 4, 5, 6),
                                                    labels = c("星期天", "星期一", "星期二", "星期三", "星期四",
                                                               "星期五", "星期六"),
                                                    alternating = 1,
                                                    cex = 0.6,
                                                    tck=0)),
                                    xlim =c(0.4, 54.6),
                                    ylim=c(6.6,-0.6),
                                    cuts= ncolors - 1,
                                    col.regions = (calendar.pal(ncolors)),
                                    xlab="" ,
                                    ylab="",
                                    colorkey= list(col = calendar.pal(ncolors), width = 0.6, height = 0.5),
                                    subscripts=TRUE
        ) )
        panel.locs <- trellis.currentLayout()
        for (row in 1:nrow(panel.locs)) {
                for (column in 1:ncol(panel.locs))  {
                        if (panel.locs[row, column] > 0)
                        {
                                trellis.focus("panel", row = row, column = column,
                                              highlight = FALSE)
                                xyetc <- trellis.panelArgs()
                                subs <- caldat[xyetc$subscripts,]
                                dates.fsubs <- caldat[caldat$yr == unique(subs$yr),]
                                y.start <- dates.fsubs$dotw[1]
                                y.end   <- dates.fsubs$dotw[nrow(dates.fsubs)]
                                dates.len <- nrow(dates.fsubs)
                                adj.start <- dates.fsubs$woty[1]
                                
                                for (k in 0:6) {
                                        if (k < y.start) {
                                                x.start <- adj.start + 0.5
                                        } else {
                                                x.start <- adj.start - 0.5
                                        }
                                        if (k > y.end) {
                                                x.finis <- dates.fsubs$woty[nrow(dates.fsubs)] - 0.5
                                        } else {
                                                x.finis <- dates.fsubs$woty[nrow(dates.fsubs)] + 0.5
                                        }
                                        grid.lines(x = c(x.start, x.finis), y = c(k -0.5, k - 0.5), 
                                                   default.units = "native", gp=gpar(col = "grey", lwd = 1))
                                }
                                if (adj.start <  2) {
                                        grid.lines(x = c( 0.5,  0.5), y = c(6.5, y.start-0.5), 
                                                   default.units = "native", gp=gpar(col = "grey", lwd = 1))
                                        grid.lines(x = c(1.5, 1.5), y = c(6.5, -0.5), default.units = "native",
                                                   gp=gpar(col = "grey", lwd = 1))
                                        grid.lines(x = c(x.finis, x.finis), 
                                                   y = c(dates.fsubs$dotw[dates.len] -0.5, -0.5), default.units = "native",
                                                   gp=gpar(col = "grey", lwd = 1))
                                        if (dates.fsubs$dotw[dates.len] != 6) {
                                                grid.lines(x = c(x.finis + 1, x.finis + 1), 
                                                           y = c(dates.fsubs$dotw[dates.len] -0.5, -0.5), default.units = "native",
                                                           gp=gpar(col = "grey", lwd = 1))
                                        }
                                        grid.lines(x = c(x.finis, x.finis), 
                                                   y = c(dates.fsubs$dotw[dates.len] -0.5, -0.5), default.units = "native",
                                                   gp=gpar(col = "grey", lwd = 1))
                                }
                                for (n in 1:51) {
                                        grid.lines(x = c(n + 1.5, n + 1.5), 
                                                   y = c(-0.5, 6.5), default.units = "native", gp=gpar(col = "grey", lwd = 1))
                                }
                                x.start <- adj.start - 0.5
                                
                                if (y.start > 0) {
                                        grid.lines(x = c(x.start, x.start + 1),
                                                   y = c(y.start - 0.5, y.start -  0.5), default.units = "native",
                                                   gp=gpar(col = "black", lwd = 1.75))
                                        grid.lines(x = c(x.start + 1, x.start + 1),
                                                   y = c(y.start - 0.5 , -0.5), default.units = "native",
                                                   gp=gpar(col = "black", lwd = 1.75))
                                        grid.lines(x = c(x.start, x.start),
                                                   y = c(y.start - 0.5, 6.5), default.units = "native",
                                                   gp=gpar(col = "black", lwd = 1.75))
                                        if (y.end < 6  ) {
                                                grid.lines(x = c(x.start + 1, x.finis + 1),
                                                           y = c(-0.5, -0.5), default.units = "native",
                                                           gp=gpar(col = "black", lwd = 1.75))
                                                grid.lines(x = c(x.start, x.finis),
                                                           y = c(6.5, 6.5), default.units = "native",
                                                           gp=gpar(col = "black", lwd = 1.75))
                                        } else {
                                                grid.lines(x = c(x.start + 1, x.finis),
                                                           y = c(-0.5, -0.5), default.units = "native",
                                                           gp=gpar(col = "black", lwd = 1.75))
                                                grid.lines(x = c(x.start, x.finis),
                                                           y = c(6.5, 6.5), default.units = "native",
                                                           gp=gpar(col = "black", lwd = 1.75))
                                        }
                                } else {
                                        grid.lines(x = c(x.start, x.start),
                                                   y = c( - 0.5, 6.5), default.units = "native",
                                                   gp=gpar(col = "black", lwd = 1.75))
                                }
                                
                                if (y.start == 0 ) {
                                        if (y.end < 6  ) {
                                                grid.lines(x = c(x.start, x.finis + 1),
                                                           y = c(-0.5, -0.5), default.units = "native",
                                                           gp=gpar(col = "black", lwd = 1.75))
                                                grid.lines(x = c(x.start, x.finis),
                                                           y = c(6.5, 6.5), default.units = "native",
                                                           gp=gpar(col = "black", lwd = 1.75))
                                        } else {
                                                grid.lines(x = c(x.start + 1, x.finis),
                                                           y = c(-0.5, -0.5), default.units = "native",
                                                           gp=gpar(col = "black", lwd = 1.75))
                                                grid.lines(x = c(x.start, x.finis),
                                                           y = c(6.5, 6.5), default.units = "native",
                                                           gp=gpar(col = "black", lwd = 1.75))
                                        }
                                }
                                for (j in 1:12)  {
                                        last.month <- max(dates.fsubs$seq[dates.fsubs$month == j])
                                        x.last.m <- dates.fsubs$woty[last.month] + 0.5
                                        y.last.m <- dates.fsubs$dotw[last.month] + 0.5
                                        grid.lines(x = c(x.last.m, x.last.m), y = c(-0.5, y.last.m),
                                                   default.units = "native", gp=gpar(col = "black", lwd = 1.75))
                                        if ((y.last.m) < 6) {
                                                grid.lines(x = c(x.last.m, x.last.m - 1), y = c(y.last.m, y.last.m),
                                                           default.units = "native", gp=gpar(col = "black", lwd = 1.75))
                                                grid.lines(x = c(x.last.m - 1, x.last.m - 1), y = c(y.last.m, 6.5),
                                                           default.units = "native", gp=gpar(col = "black", lwd = 1.75))
                                        } else {
                                                grid.lines(x = c(x.last.m, x.last.m), y = c(- 0.5, 6.5),
                                                           default.units = "native", gp=gpar(col = "black", lwd = 1.75))
                                        }
                                }
                        }
                }
                trellis.unfocus()
        } 
        lattice.options(default.theme = def.theme)
}

dataset <- dataset %>%
        mutate(hour = hour(DateTime),day = day(DateTime),month = month(DateTime),year = year(DateTime),week = weekdays.Date(DateTime), weeks = week(DateTime)) %>%
        mutate(ymd = as.Date(DateTime,format = "%m/%d/%y"))

dat2 <- dataset %>%
        select(ymd,`ug/m3`)

dat3 <- aggregate(dat2[, 2], list(dat2$ymd), mean,na.rm = T)

CairoPDF("pmall.pdf", width = 10, height = 6.18);
showtext.auto()
font.add("test","Kaiti.ttf")
calendarHeat(dates=dat3$Group.1,
             values=dat3$`ug/m3`)
dev.off()
