## Figure 5 selection in high and low salinity

library(ggplot2)

#get the data

F1_sal<-read.table("F1_Capsule_data.csv", header=TRUE, sep=",")
SvN<-read.csv("SoCalvNorCal_Data.csv", header=TRUE, sep=",")



#make a plot for before and after selection by family
#data for before selection is in SvN
head(SvN)
#add a column for before selection
SvN$Generation<-as.vector(strrep("Parental", 1))
#adjust the column names
colnames(SvN)<- c("date","order","Capsule","Family","ID","eggmass","Salinity","current","region","starved","Generation")

#add a column to F1_sal for after selection
F1_sal$Generation<-ifelse(F1_sal$mom_env == "20C16ppt","low_selection","high_selection")

#rename off_id to ID and offspring_env to Salinity
colnames(F1_sal)<-c("date", "order", "Capsule", "mom_id", "mom_env", "ID", "eggmass", "Salinity", "Family", "Generation")

#put svn and f1_sal together
fig_5<-rbind(SvN[,c("date", "Capsule", "ID", "Salinity", "Generation", "Family", "eggmass")], 
             F1_sal[,c("date", "Capsule", "ID", "Salinity", "Generation", "Family", "eggmass")])

#by proportion

#get the mean egg capsule size by individual
fig5_ID_mean<-aggregate(fig_5$Capsule, list(fig_5$date, fig_5$ID, fig_5$eggmass), mean)
head(fig5_ID_mean)
colnames(fig5_ID_mean)<-c("date", "ID", "eggmass", "mean_cap")

#now merge them together with the original data so we can get the Family etc info
fig5_ID_mean_tog<-merge(fig5_ID_mean, fig_5, by=c("ID"))
fig5_ID_mean_tog<-subset(fig5_ID_mean_tog, select = -c(Capsule))
fig5_mean<-unique.data.frame(fig5_ID_mean_tog)

head(fig5_mean)


#now define L and P

fig5_mean$type_150 <- ifelse(fig5_mean$mean_cap > 0.150,"l","p")

#then define the proportions that we want using table, we must do this for before and after selection separately
#before
fig5_data_sel_bs<-fig5_mean[fig5_mean$Generation == "Parental",]
freq_150_bs<-as.data.frame(table(fig5_data_sel_bs$type_150, fig5_data_sel_bs$Family:fig5_data_sel_bs$Salinity))

#after low
fig5_data_sel_low<-fig5_mean[fig5_mean$Generation == "low_selection",]
freq_150_low<-as.data.frame(table(fig5_data_sel_low$type_150, fig5_data_sel_low$Family:fig5_data_sel_low$Salinity))

#after high
fig5_data_sel_high<-fig5_mean[fig5_mean$Generation == "high_selection",]
freq_150_high<-as.data.frame(table(fig5_data_sel_high$type_150, fig5_data_sel_high$Family:fig5_data_sel_high$Salinity))
#now stitch them back together

freq_150_bs$Generation<-as.vector(strrep("Parental", 1))
freq_150_low$Generation<-as.vector(strrep("Low Salinity", 1))
freq_150_high$Generation<-as.vector(strrep("High Salinity", 1))

freq_150<-rbind(freq_150_bs, freq_150_high, freq_150_low)

#we need to get it back to being a dataframe
freq_150$Var2<-as.vector(freq_150$Var2)
splity<-unlist(strsplit(freq_150$Var2, ":"))

Family<-splity[seq(1, length(splity), by=2)]
ENv<-splity[seq(2, length(splity), by=2)]

freq_150<-cbind(freq_150, Family, ENv)
freq_150<-within(freq_150, rm(Var2))
head(freq_150)

#now we can set up a for loop to calculate the proportion

fig5_prop_l_150<-seq(1:length(freq_150$Var1))

for (i in seq(1,length(freq_150$Var1),2)){
  calc_prop<-freq_150$Freq[i]/(freq_150$Freq[i+1]+freq_150$Freq[i])
  fig5_prop_l_150[i]<-calc_prop
  
}
#this just returns the proportions, see
head(fig5_prop_l_150)

fig5_prop_l_150<-cbind(fig5_prop_l_150, freq_150)

head(fig5_prop_l_150)

#yay, remove the p because they don't make sense

fig5_prop_l_150<-fig5_prop_l_150[fig5_prop_l_150$Var1 == "l",]

head(fig5_prop_l_150)#looks like several families did not have data for after selection,   
                      #they are lowering the prop (having zero values)

#to make the plot pretty change env labels
levels(fig5_prop_l_150$ENv)[levels(fig5_prop_l_150$ENv) == "20C16ppt"]<-"Low 16 ppt"
levels(fig5_prop_l_150$ENv)[levels(fig5_prop_l_150$ENv) == "20C32ppt"]<-"High 32 ppt"


#removing families SCW46, SCW47, W143, W149
fig5_prop_l_150<-fig5_prop_l_150[fig5_prop_l_150$Freq != 0,]

#now plot it up

ggplot(data=fig5_prop_l_150, aes(x=ENv, y=fig5_prop_l_150, group=Generation, color=Generation))+
  stat_summary(fun.y=mean, geom="point", size=1)+
  stat_summary(fun.y=mean, geom="line", size=1)+
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=0.1), size=0.5)+
  #facet_grid(~Generation)+
  scale_color_manual(values=c("black", "grey", "blue"))+
  theme_classic()+
  labs(x="Salinity", y="Proportion Lecithotrophic")


sel_mod_prop<-lm(data=fig5_prop_l_150, fig5_prop_l_150~Generation*ENv-1)
summary(sel_mod_prop)
