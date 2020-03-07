require(ggplot2)

#get the data

map_pop<-read.csv("Map_pop_data.csv", header=TRUE, sep=",")
SvN<-read.csv("SoCalvNorCal_Data.csv", header=TRUE, sep=",")

#pull-out first generation
p_map_pop<-map_pop[map_pop$gen == "p",]
p_map_pop$date<-as.Date(p_map_pop$date, "%m/%d/%Y")

#make svn do the dates right

SvN$date<-as.Date(SvN$date, "%m/%d/%Y")

#add a generation column to svn
SvN$gen<-"p"

#rbind SvN and maping pops to make a master high versus low (hvl) data set

HvL<-rbind(SvN[,c("Capsule", "original", "date", "Family", "ID", "gen")],p_map_pop[, c("Capsule", "original", "date", "Family", "ID", "gen")])

#calculate the mean and mode for HvL

mean(HvL$Capsule)
#[1] 0.1303721

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

getmode(HvL$Capsule)
#[1] 0.11


find_modes<- function(x) {
  modes <- NULL
  for ( i in 2:(length(x)-1) ){
    if ( (x[i] > x[i-1]) & (x[i] > x[i+1]) ) {
      modes <- c(modes,i)
    }
  }
  if ( length(modes) == 0 ) {
    modes = 'This is a monotonic distribution'
  }
  return(modes)
}
plot(density(HvL$Capsule))

mymodes<-find_modes(density(HvL$Capsule)$y)

density(HvL$Capsule)$x[mymodes]
#[1] 0.06824804 0.11277401 0.16583922 0.18230773 0.19572652 0.23049337 0.25306133 0.27257957 0.30246687 0.34699284 0.35919174

density(HvL$Capsule)$y[mymodes]
#[1]  0.12200067 23.53159239  1.74180354  2.59124926  2.58617408  1.36721120  0.41914932  0.17099219
#[9]  0.06258092  0.02288615  0.02293676

#thus the two highest modes are 0.1127 and 0.1823

#how many individuals per salinity treatment per generation?

table(HvL[unique(HvL$ID),]$gen, HvL[unique(HvL$ID),]$original)
#20C16ppt 20C32ppt
#p      189      244
table(map_pop[unique(map_pop$ID),]$gen, map_pop[unique(map_pop$ID),]$original)

### LETS LOOK AT HvL by individual mean capsule size. this is important
### use the function aggregate, it outputs a dataframe instead of tapply which outputs a vector


ID_mean<-aggregate(HvL$Capsule, list(HvL$ID) , mean)
colnames(ID_mean)<-c("ID", "mean_cap")

#now merge them together with the original data so we can get the Family etc info
ID_mean_tog<-merge(ID_mean, HvL, by=c("ID"))
ID_mean_tog<-subset(ID_mean_tog, select = -c(Capsule, date))
ID_mean_tog_uni<-unique.data.frame(ID_mean_tog)

#now add on the type for different sizes

ID_mean_tog_uni$type_150 <- ifelse(ID_mean_tog_uni$mean_cap>0.150,"l","p")


freq_150<-as.data.frame(table(ID_mean_tog_uni$type_150, ID_mean_tog_uni$Family:ID_mean_tog_uni$original))
head(freq_150)

#get the dataframe ready to calulate the proportion
freq_150$Var2<-as.vector(freq_150$Var2)
splity<-unlist(strsplit(freq_150$Var2, ":"))

Family<-splity[seq(1, length(splity), by=2)]
ENv<-splity[seq(2, length(splity), by=2)]

freq_150<-cbind(freq_150, Family, ENv)
freq_150<-within(freq_150, rm(Var2))
head(freq_150)

#this has values for all of the families that we have in the map_pop dataframe. 
#we want to remove the 2nd and 3rd gen families - we will tell it which fams to keep

the_fams<-c("SCW46", "SCW47", "SCW48", "SCW49", "SCW50",
            "W128", "W129", "W130", "W131", "W132", "W133", "W134",
            "W143", "W144", "W148","W149", "W150")

freq_l_150_fams<-freq_150[is.element(freq_150$Family, the_fams), ]


head(freq_l_150_fams)
str(freq_l_150_fams)
length(unique(freq_l_150_fams$Family))

#write a for loop to calulate the proportion
#create an empty vector to hold the values of prop_l

prop_l_150<-seq(1:length(freq_l_150_fams$Var1))

for (i in seq(1,length(freq_l_150_fams$Var1),2)){
  calc_prop<-freq_l_150_fams$Freq[i]/(freq_l_150_fams$Freq[i+1]+freq_l_150_fams$Freq[i])
  prop_l_150[i]<-calc_prop
  
}

head(prop_l_150)

#bind it to the data frame with the family and env info

prop_l_150<-cbind(prop_l_150, freq_l_150_fams)

head(prop_l_150)

#pull out just the lecithotrophic props(since those are the only ones that make sense)

prop_l_150_l<-prop_l_150[prop_l_150$Var1 == "l",]

#there are repeating columns in the dataframe pullout the ones we want

prop_l_final<-prop_l_150_l[,c(1,4,5)]

colors <- c("black", "black", "black", "black", "black", #socal black
            "#74C476", "#74C476", "#74C476", "#74C476", "#74C476", #tomales = the greens
            "lightblue", "lightblue", "lightblue", "lightblue", "lightblue", "lightblue", "lightblue" ) #mill valley = the blues


prop_l_final$Family<-sort(prop_l_final$Family)

levels(prop_l_final$ENv)<- c("Low", "High")

#now plot that shit

fig2a<-ggplot(data=prop_l_final, aes(x=ENv, y=prop_l_150, group=Family, color=Family))+
  geom_point(size=4)+ geom_line(size=2)+
  #scale_shape_manual(values=shapes)+
  scale_color_manual(values=colors)+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), 
                     plot.title = element_text(size = 16, face = "bold"),
                     text = element_text(size = 16),
                     legend.position="none")+
  labs(x="Salinity", y="Proportion lecithotrophic", title="A")


#now for figure 2b ploting the proportion of low salinity against high salinity

#make a table of the freq of l and p

type_150 <- ifelse(ID_mean_tog_uni$mean_cap > 0.150,"l","p")
ID_mean_tog_uni<-cbind(ID_mean_tog_uni, type_150)
head(ID_mean_tog_uni)

freq_150<-as.data.frame(table(ID_mean_tog_uni$type_150, ID_mean_tog_uni$Family:ID_mean_tog_uni$original))
head(freq_150)
#get the dataframe ready to calulate the proportion
freq_150$Var2<-as.vector(freq_150$Var2)
splity<-unlist(strsplit(freq_150$Var2, ":"))

Family<-splity[seq(1, length(splity), by=2)]
ENv<-splity[seq(2, length(splity), by=2)]

freq_150<-cbind(freq_150, Family, ENv)
freq_150<-within(freq_150, rm(Var2))
head(freq_150)

#this has values for all of the families that we have in the map_pop dataframe. 
#we want to remove the 2nd and 3rd gen families - we will tell it which fams to keep

the_fams<-c("SCW46", "SCW47", "SCW48", "SCW49", "SCW50",
            "W128", "W129", "W130", "W131", "W132", "W133", "W134",
            "W143", "W144", "W148","W149", "W150")

freq_l_150_fams<-freq_150[is.element(freq_150$Family, the_fams), ]


head(freq_l_150_fams)
str(freq_l_150_fams)
length(unique(freq_l_150_fams$Family))

#write a for loop to calulate the proportion
#create an empty vector to hold the values of prop_l

prop_l_150<-seq(1:length(freq_l_150_fams$Var1))

for (i in seq(1,length(freq_l_150_fams$Var1),2)){
  calc_prop<-freq_l_150_fams$Freq[i]/(freq_l_150_fams$Freq[i+1]+freq_l_150_fams$Freq[i])
  prop_l_150[i]<-calc_prop
  
}

head(prop_l_150)


#bind it to the data frame with the family and env info

prop_l_150<-cbind(prop_l_150, freq_l_150_fams)

head(prop_l_150)

#pull out just the lecithotrophic props(since those are the only ones that make sense)

prop_l_150_l<-prop_l_150[prop_l_150$Var1 == "l",]

names(colors)<-the_fams

coef(lm(subset(prop_l_150_l, ENv %in% "20C16ppt")[,1]~subset(prop_l_150_l, ENv %in% "20C32ppt")[,1]))
# (Intercept) 
# 0.009200238 
# subset(prop_l_150_l, ENv %in% "20C32ppt")[, 1] 
# 0.628004209  

#we need to sort the families to get colors to work properly

prop_l_150_l$Family<-sort(prop_l_150_l$Family)

#now plot that shit

fig2b<-ggplot()+
  geom_point(aes(y=subset(prop_l_150_l, ENv %in% "20C16ppt")[,1],x=subset(prop_l_150_l, ENv %in% "20C32ppt")[,1]), 
             group=colors, color=colors, size=4)+
  geom_abline(intercept = 0.009200238 , slope = 0.628004209 , size=2)+
  #scale_color_manual(values=colors)+
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"), 
                     plot.title = element_text(size = 16, face = "bold"),
                     text = element_text(size = 16),
                     legend.position="right")+
  labs(x="Proportion lecithotrophic in high salinity", 
       y="Proportion lecithotrophic in low salinity", title = "B")


#now that we've done that we can put them on a single plot

require(gridExtra)

lay <- rbind(c(1,2),
             c(1,2))


grid.arrange(fig2a, fig2b, layout_matrix=lay)

#do I want legends on this plot?? 
  # YES by region

#divide sites by region

head(prop_l_final)

unique(prop_l_final$Family)
# SCW46 SCW47 SCW48 SCW49 SCW50 W128  W129  W130  W131  W132  W133  W134  W143  W144  W148  W149  W150 


prop_l_final$region[prop_l_final$Family == "SCW46"] = "SoCal" 
prop_l_final$region[prop_l_final$Family == "SCW47"] = "SoCal" 
prop_l_final$region[prop_l_final$Family == "SCW48"] = "SoCal" 
prop_l_final$region[prop_l_final$Family == "SCW49"] = "SoCal" 
prop_l_final$region[prop_l_final$Family == "SCW50"] = "SoCal" 
prop_l_final$region[prop_l_final$Family == "W128"] = "TomBay" 
prop_l_final$region[prop_l_final$Family == "W129"] = "TomBay" 
prop_l_final$region[prop_l_final$Family == "W130"] = "TomBay"      
prop_l_final$region[prop_l_final$Family == "W131"] = "TomBay" 
prop_l_final$region[prop_l_final$Family == "W132"] = "TomBay" 
prop_l_final$region[prop_l_final$Family == "W133"] = "TomBay" 
prop_l_final$region[prop_l_final$Family == "W134"] = "TomBay" 
prop_l_final$region[prop_l_final$Family == "W143"] = "MillVal" 
prop_l_final$region[prop_l_final$Family == "W144"] = "MillVal"
prop_l_final$region[prop_l_final$Family == "W148"] = "MillVal"
prop_l_final$region[prop_l_final$Family == "W149"] = "MillVal"
prop_l_final$region[prop_l_final$Family == "W150"] = "MillVal"
        
head(prop_l_final)

#calculate the mean by region
tapply(prop_l_final$prop_l_150, prop_l_final$region, mean)


