### composite distance heatmap 
setwd("/Users/lylaatta/Documents/GitHub/veloviz/graphViz/figures")
#using velo_merfish environment 
library(ggplot2)
library(gridExtra)
library(reshape2)
library(viridis)
library(cowplot)
library(svglite)



#### Calculate composite distance for range of angles and euclidean distances
#### 
#EucDis weight = 1 
logd = seq(-3,1,0.01)
d = 10^logd #euclidean distance bw neighbor and projected state
a = as.numeric(formatC(seq(0,2*pi,((2*pi)/500)),digits = 3)) #angle bw velocity vector and 
negsim = -1*cos(a)

inveuc = 1/(1+(d))
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
h <- ggplot(cdLong, aes(x=Angle, y=Distance, fill=CD)) + geom_tile() + 
                        labs(x="Angle",y="Euclidean Distance",fill = "Composite\nDistance") + 
                        scale_fill_viridis(limits = c(-1,1), breaks = c(-1,0,1)) +
                        scale_y_discrete(breaks = c(0.001,0.01,0.1,1,10)) + 
                        scale_x_discrete(breaks = c(0,3.14,6.28), labels = c("0","180","360")) +
                        theme(text = element_text(size = 18))
h
ggsave(file = "heatmap_cd.png", plot = h, width = 8, height = 6)

#corresponding edge weights 
#### Make heatmap of edge weights 
ewLong = cdLong
ewLong$CD = 1 - ewLong$CD
colnames(ewLong) = c("Angle","Distance","EW")
hm2 = ggplot(ewLong, aes(x=Angle, y=Distance, fill=EW)) + geom_tile() + labs(x="Angle",y="Euclidean Distance",fill = "Edge\nWeight") + 
  scale_fill_viridis(limits = c(0,2), breaks = c(0,1,2)) +
  scale_y_discrete(expand=c(0,0), breaks = c(0.001,0.01,0.1,1,10,100)) + 
  scale_x_discrete(expand=c(0,0), breaks = c(0,3.14,6.28), labels = c("0","180","360")) +
  theme(text = element_text(size = 18))
hm2
ggsave(file = "heatmap_ew.png", plot = hm2, width = 8, height = 6)


#Change EucDis weight
ws = c(0.01,0.1,0.5,1,5,10)
hs = vector("list", length = 6)
for (w in c(1:length(ws))){
  currW = ws[w]
  inveuc = 1/(1+(currW*d))
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
  hs[[w]] <- ggplotGrob(ggplot(cdLong, aes(x=Angle, y=Distance, fill=CD)) + geom_tile() + 
    labs(title = paste("Distance Weight:",currW) ,x="Angle",y="Euclidean Distance",fill = "Composite\nDistance") + 
    scale_fill_viridis(limits = c(-1,1), breaks = c(-1,0,1)) +
    scale_y_discrete(breaks = c(0.001,0.01,0.1,1,10)) + 
    scale_x_discrete(breaks = c(0,3.14,6.28), labels = c("0","180","360")) +
    theme(text = element_text(size = 14), plot.title = element_text(size = 14, hjust = 0.5)))
}

g = grid.arrange(hs[[1]],hs[[2]],hs[[3]],hs[[4]],hs[[5]],hs[[6]],ncol = 3)
ggsave("weightedCD.png",plot = g,width = 12, height = 6)

# ggsave(file = "heatmap.svg", plot = hm, width = 8, height = 6)
# ggsave(file = "heatmap.png", plot = hm, width = 8, height = 6)



