setwd("C:\\Users\\nt78066\\OneDrive - University of Georgia\\Documents\\Abrasion paper\\")

#load packages
library(ggplot2)
library(Rmisc)
library(lme4)
library(MuMIn)
library(wesanderson)
library(lmerTest)
library(dplyr)

#load data 
d<-read.csv("respiration_data.csv")
head(d)

volume<-32.29/1000 # L of water in incubation.

d$DO.change.mg<-(d$start.do.mg.l-d$end.do.mg.l)*volume
d$DO.change.rate.mg.hr<-d$DO.change.mg/d$tim.epapsed.decimal

blanks<-aggregate(DO.change.rate.mg.hr~sample.date,data=d[d$leaf.species=="blank",],mean)

names(blanks)[names(blanks) == "DO.change.rate.mg.hr"] <- "DO.change.rate.mg.hr.blank"

d<-d[d$leaf.species!="blank",]

d.blanks<-merge(d,blanks,by="sample.date")

d.blanks$DO.change.rate.mg.hr.corrected<-d.blanks$DO.change.rate.mg.hr-d.blanks$DO.change.rate.mg.hr.blank
d.blanks$do.change.rate.mg.hr.g<-d.blanks$DO.change.rate.mg.hr.corrected/d.blanks$leaf.mass.g

#split data into fine and coarse mesh bags
fine<-d.blanks[d.blanks$bag.type=="fine",]
coarse<-d.blanks[d.blanks$bag.type=="coarse",]

#new data frames so that columns have correct names
fine2<-data.frame(ID=fine$sampleID,bag.type="fine",leaf.type=fine$leaf.species,fine.resp=fine$do.change.rate.mg.hr.g)

coarse2<-data.frame(ID=coarse$sampleID,bag.type="coarse",leaf.type=coarse$leaf.species,coarse.resp=coarse$do.change.rate.mg.hr.g)

#merge coarse and fine
d2<-merge(fine2,coarse2,by="ID")
head(d2)

#difference of logged respiration rates
d2$difference<-log(d2$coarse.resp)-log(d2$fine.resp)

# quick look 
hist(d2$difference)


#testing significance
t.test(d2[d2$leaf.type.x=="acer","difference"])
t.test(d2[d2$leaf.type.x=="rhodo","difference"])

#estimating actual differences for acer
d2$difference<-d2$fine.resp-d2$coarse.resp
d2$proportional.diff<-d2$fine.resp/d2$coarse.resp
mean(d2[d2$leaf.type.x=="acer","difference"],na.rm=TRUE)
mean(d2[d2$leaf.type.x=="acer","proportional.diff"],na.rm=TRUE)

#estimating difference for rhodo
d2$difference<-d2$coarse.resp-d2$fine.resp
d2$proportional.diff<-d2$coarse.resp/d2$fine.resp
mean(d2[d2$leaf.type.x=="rhodo","difference"],na.rm=TRUE)
mean(d2[d2$leaf.type.x=="rhodo","proportional.diff"],na.rm=TRUE)




# summarizing resp data for plotting
d.s<-summarySE(d.blanks,measurevar = "do.change.rate.mg.hr.g",groupvars = c("bag.type","leaf.species"),na.rm=TRUE)



respiration_raw<-ggplot(d.s[d.s$leaf.species!="blank",],aes(x=leaf.species,y=do.change.rate.mg.hr.g,fill=bag.type))+
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=do.change.rate.mg.hr.g-ci,ymax=do.change.rate.mg.hr.g+ci),width=0.01,position=position_dodge(0.9))+
  theme_classic()+xlab("Leaf Species")+scale_x_discrete(labels=c("acer"="Acer","rhodo"="Rhododendron"))+
  theme(text = element_text(size=20))+ylab(expression("Respiration Rate ("*"mg"~O[2]~"g"^-1~"hr"^-1*")"))+
  guides(fill=guide_legend(title=""))+scale_fill_manual(values=wes_palette("Darjeeling1")[2:3])

d.blanks$x.position<-1  
d.blanks[d.blanks$leaf.species=="acer","x.position"]<-d.blanks[d.blanks$leaf.species=="acer","x.position"]+1.5
d.blanks[d.blanks$bag.type=="fine","x.position"]<-d.blanks[d.blanks$bag.type=="fine","x.position"]+1

# plotting paired samples
respiration_paired<-ggplot(d.blanks,aes(x=x.position,y=do.change.rate.mg.hr.g,color=leaf.species))+
  geom_point(size=2)+geom_line(aes(group = sampleID),size=1)+
    scale_x_continuous(breaks = c(1,2,2.5,3.5), labels = paste0(c("Coarse","Fine","Coarse","Fine")))+
  theme_classic()+xlab("Bag Type")+
  theme(text = element_text(size=20))+ylab(expression("Respiration Rate ("*"mg"~O[2]~"g"^-1~"hr"^-1*")"))+
  guides(fill=guide_legend(title=""))+scale_color_manual(values=wes_palette("Darjeeling1")[4:5], labels = c("Acer","Rhodo"),name="")

tiff(filename="paired_respiration_5aug2021.tiff",units="in",res=800,width=6,height=10,compression="lzw")
plot_grid(respiration_raw,respiration_paired,labels="AUTO",ncol=1,label_x=0.9,label_y=0.9,label_size = 20)
dev.off()

#### centering leaf mass so that the intercepts are a clearer estimate of effect size.
acer.model.data<-d.blanks[d.blanks$leaf.species=="acer",]
acer.model.data$leaf.mass.centered<-acer.model.data$leaf.mass.g-mean(acer.model.data$leaf.mass.g,na.rm=TRUE)

# linear mixed effects models 
model2<-lmer(do.change.rate.mg.hr.g~bag.type*leaf.mass.centered+(1|sampleID),data=acer.model.data)
plot(model2)
anova(model2)
summary(model2)
r.squaredGLMM(model2)

rhodo.model.data<-d.blanks[d.blanks$leaf.species=="rhodo",]
rhodo.model.data$leaf.mass.centered<-rhodo.model.data$leaf.mass.g-mean(rhodo.model.data$leaf.mass.g,na.rm=TRUE)


model3<-lmer(do.change.rate.mg.hr.g~bag.type+leaf.mass.centered+(1|sampleID),data=rhodo.model.data)
plot(model3)
anova(model3)
summary(model3)
r.squaredGLMM(model3)

ggplot(acer.model.data,aes(x=leaf.mass.g,y=do.change.rate.mg.hr.g,color=bag.type))+geom_point()+geom_smooth(method="lm",se=FALSE)
