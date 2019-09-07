library(enviGCMS)
library(pmd)
# 将csv文件放到你的工作目录下
water <- getmzrtcsv('xcms.csv')
# 保留时间换成秒
water$rt <- water$rt*60
# pmd分析
# 先找到母离子 346.92809
(mass <- water$mz[water$mz>346.5&water$mz<347])
# 应该是346.9274
# 寻找高频pmd
pmd <- globalstd(water,ng=NULL,digits = 3,top = 20)
pmdlist <- unique(pmd$sda$diff2)
# 构建高频pmd网络，diff里指定pmd，mass 里指定母离子，digits 指定质量精确度
pmdnet <- getchain(pmd,diff = pmdlist,mass = 346.9274, digits = 3)
# 作图
library(igraph)
library(RColorBrewer)
pal <- (grDevices::colorRampPalette(rev(
                        RColorBrewer::brewer.pal(11,
                                                 "RdYlBu")
                )))(20)
net <- graph_from_data_frame(pmdnet$sdac,directed = F)
plot(net,vertex.size =5,edge.width = 5,edge.color = pal[as.numeric(as.factor(E(net)$diff2))],vertex.label=ifelse(round(as.numeric(V(net)$name),4) == 346.9274,'Compounds',NA),vertex.label.dist=1,vertex.color=ifelse(round(as.numeric(V(net)$name),4) == 346.9274,'red','black'),main = 'PMD network')
legend("topright",bty = "n",
       legend=unique(E(net)$diff2),
       fill=unique(pal[as.numeric(as.factor(E(net)$diff2))]), border=NA,horiz = F)
