# Egg diameter figure for panel 1E

#load the data

map_pop<-read.csv("Map_pop_data.csv", header=TRUE, sep=",")
SvN<-read.csv("SoCalvNorCal_Data.csv", header=TRUE, sep=",")
egg_di<-read.csv("EggDiameterMeasurements.csv", header=TRUE, sep=",")

#make svn do the dates right
SvN$date<-as.Date(SvN$date, "%m/%d/%Y")
map_pop$date<-as.Date(map_pop$date, "%m/%d/%Y")

#rbind SvN and maping pops to make a master high versus low (hvl) data set

HvL<-rbind(SvN[,c("Capsule", "original", "date", "Family", "ID", "eggmass")],
           map_pop[, c("Capsule", "original", "date", "Family", "ID", "eggmass")])
HvL$date<-as.Date(HvL$date, "%m/%d/%Y")


#look at mean of egg di and mean of egg capsule

#calculate the mean per egg mass for egg di
egg_di_mean<-aggregate(egg_di$EggDi, list(egg_di$date, egg_di$ID, egg_di$eggmass), mean)
colnames(egg_di_mean)<-c("date", "ID", "eggmass", "eggdi")
egg_di_mean$eggmass<-as.numeric(egg_di_mean$eggmass)
egg_di_mean$date<-as.Date(egg_di_mean$date, "%m/%d/%Y")

#again for hvl
HvL_mean<-aggregate(HvL$Capsule, list(HvL$date, HvL$ID, HvL$eggmass), mean)
colnames(HvL_mean)<-c("date", "ID", "eggmass", "Capsule")
HvL$eggmass<-as.numeric(HvL$eggmass)

#merge the two datasets 
egg_di_map<-merge(egg_di_mean, HvL_mean, by=c("date", "ID", "eggmass"))
head(egg_di_map)

#change to um
egg_di_map$egg_di_um<-egg_di_map$eggdi*1000
egg_di_map$egg_cp_um<-egg_di_map$Capsule*1000

#generate simple plot in basic R
plot(data=egg_di_map, egg_di_um~egg_cp_um, col=eggcol, pch=16, cex=1.5, xlab="Egg Capsule Diameter (um)", ylab="Egg Diameter (um)")
abline(v=150, lty=2)


#some of these points are close together. Are these mixed egg masses? Let's see
#write a function to define eggs by the egg capsule diameter
map_pop_mlp<-aggregate(HvL$Capsule, list(HvL$date, HvL$eggmass, HvL$ID), function(x) ifelse(x > 0.15, "l", "p"))

colnames(map_pop_mlp)<-c("date", "eggmass", "ID", "x")
map_pop_mlp$is_l<-grepl("l", map_pop_mlp$x)
map_pop_mlp$is_p<-grepl("p", map_pop_mlp$x)
map_pop_mlp$is_m<- map_pop_mlp$is_l == T & map_pop_mlp$is_p == T
map_pop_mlp$eggmass_type<-ifelse(map_pop_mlp$is_m == T, "m", ifelse(map_pop_mlp$is_l == T, "l", "p"))
table(map_pop_mlp$eggmass_type)

head(map_pop_mlp)

egg_id<-merge(map_pop_mlp, egg_di_map, by=c("date", "ID", "eggmass"))

table(egg_id$eggmass_type)

ggplot(data=egg_id, aes(x=egg_cp_um, y=egg_di_um, col=eggmass_type))+
  geom_point(size=3, cex=2)+
  scale_color_manual(breaks=c("l", "m",  "p"),
                    labels=c("Lecithotrophic", "Mixed", "Planktotrophic"),
                    name="Developmental mode",
                    values=c("goldenrod", "aquamarine3",  "blue"))+
  labs(x=expression(paste("Mean Egg Capsule Diameter (", mu, "m)")), 
       y=expression(paste("Mean Egg Diameter (", mu,"m)")))+
  theme_classic()+theme(legend.position = "right", 
                        axis.text=element_text(size=12),
                        axis.title=element_text(size=14,face="bold"), 
                        legend.title = element_text(size=12), 
                        legend.text=element_text(size=12))+
  geom_vline(xintercept = 150, lty=2)
  



