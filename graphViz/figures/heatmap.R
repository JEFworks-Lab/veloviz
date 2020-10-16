### composite distance heatmap 
setwd("/Users/lylaatta/Documents/GitHub/veloviz/graphViz/figures")
library(ggplot2)
library(reshape2)
library(viridis)


#### Calculate composite distance for range of angles and euclidean distances
logd = seq(-3,0,0.01)
d = 10^logd #euclidean distance bw neighbor and projected state
a = as.numeric(formatC(seq(0,2*pi,((2*pi)/500)),digits = 3)) #angle bw velocity vector and 
negsim = -1*cos(a)
inveuc = 1/d
pwiseCDs = sapply(inveuc,function(x) x*negsim) #composite distance col=dist, row=angle 

#### Make heatmap of composite dist 
#melt for ggplot
colnames(pwiseCDs) = as.character(d)
rownames(pwiseCDs) = as.character(a)
cdLong = melt(pwiseCDs)
colnames(cdLong) = c("Angle","Distance","CD")
cdLong$Angle = as.character(cdLong$Angle)
cdLong$Distance = as.character(cdLong$Distance)
#plot heatmap
hm = ggplot(cdLong, aes(x=Angle, y=Distance, fill=CD)) + geom_tile() + labs(x="Angle",y="Euclidean Distance",fill = "Composite\nDistance") + 
  scale_fill_viridis(limits = c(-1000,1000), breaks = c(-1000,0,1000)) +
  scale_y_discrete(expand=c(0,0), breaks = c(0.001,0.01,0.1,1,10,100)) + 
  scale_x_discrete(expand=c(0,0), breaks = c(0,3.14,6.28), labels = c("0","180","360")) +
  theme(text = element_text(size = 20))
hm
ggsave(file = "heatmap.svg", plot = hm, width = 8, height = 6)
ggsave(file = "heatmap.png", plot = hm, width = 8, height = 6)


#### Make heatmap of edge weights 
ewLong = cdLong
ewLong$CD = max(ewLong$CD) - ewLong$CD
colnames(ewLong) = c("Angle","Distance","EW")
hm2 = ggplot(ewLong, aes(x=Angle, y=Distance, fill=EW)) + geom_tile() + labs(x="Angle",y="Euclidean Distance",fill = "Edge\nWeight") + 
  scale_fill_viridis(limits = c(0,2000), breaks = c(0,1000,2000)) +
  scale_y_discrete(expand=c(0,0), breaks = c(0.001,0.01,0.1,1,10,100)) + 
  scale_x_discrete(expand=c(0,0), breaks = c(0,3.14,6.28), labels = c("0","180","360")) +
  theme(text = element_text(size = 20))
hm2
ggsave(file = "heatmap_ew.svg", plot = hm2, width = 8, height = 6)
ggsave(file = "heatmap_ew.png", plot = hm2, width = 8, height = 6)




#plot(c(1:length(cdLong$Distance)),as.numeric(cdLong$Distance),log="y")


