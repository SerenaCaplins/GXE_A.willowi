
require(ggplot2)


#get the data from map_pop
map_pop<-read.csv("C:/Users/SAPCaps/Documents/Projects/Alderia_papers_pics//pics/IndivPlast_MappingPops/Date_forR.csv", header=TRUE, sep=",")

#look at mean for individual egg masses

map_pop_eggmean<-aggregate(map_pop$Capsule, list(map_pop$date, map_pop$ID, map_pop$eggmass), mean)

colnames(map_pop_eggmean)<- c("date", "ID", "eggmass", "mean_cap")

map_pop_eggmean_tog<-merge(map_pop_eggmean, map_pop, by=c("date", "ID", "eggmass"))

#remove the F2 in low salinity data
map_pop_eggmean_test<-map_pop_eggmean_tog[!(map_pop_eggmean_tog$gen == "F2" & map_pop_eggmean_tog$original == "Low Salinity 16 ppt"),]

#set the factor levels for gen to get the plot to plot p first
map_pop_eggmean_tog$gen<-factor(map_pop_eggmean_tog$gen, levels=c("p", "F1", "F2", "F3"))

#change the salinity labels
levels(map_pop_eggmean_tog$original)

levels(map_pop_eggmean_tog$original)<-c("Low Salinity 16 ppt", "High Salinity 32 ppt")

#remove F2 from the low salinity treatment
map_pop_eggmean_tog_fin<-map_pop_eggmean_tog[c(map_pop_eggmean_tog$gen == "F2" & map_pop_eggmean_tog$original == "Low Salinity 16 ppt")==F,]


#change the levels of gen
levels(map_pop_eggmean_tog_fin$gen)<- c("p","S1","S2","S3")

levels(map_pop$gen)<- c("S1","S2","S3","P")

ggplot(data=subset(map_pop, gen %in% c("P", "S1", "S2")))+
  #geom_boxplot(aes(y=mean_cap,x=gen, color=gen))+
  geom_histogram(position="identity", bins=40, aes(x=Capsule*1000, y=..density..,  group=gen, color=gen, fill=gen), alpha=0.2, size=1)+
  geom_density(aes(x=Capsule*1000, y=..density.., group=gen,  color=gen, fill=gen), alpha=0.1, size=1)+
  scale_color_manual(values=c("aquamarine3", "gold","blue"))+
  scale_fill_manual(values=c("aquamarine3", "gold", "blue"))+
  geom_vline(xintercept=150, linetype=2)+
  facet_wrap(~original)+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), 
                     plot.title = element_text(size = 16, face = "bold"),
                     text = element_text(size = 16))+ 
  guides(fill=guide_legend(title="Generation"), 
         color=guide_legend(title="Generation"))+
  labs(y="Density", x=expression(paste("Egg capsule size (", mu,"m)")))

