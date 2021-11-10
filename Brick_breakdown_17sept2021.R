
## Load packages
library(lme4)
library(ggplot2)
library(lmerTest)
library(Rmisc)
library(reshape2)
library(MuMIn)
library(cowplot)
library(tidyverse)
library(wesanderson)

# Load data
setwd("C:\\Users\\nt78066\\OneDrive - University of Georgia\\Documents\\Abrasion paper\\blocks\\")

flow.sum<-read.csv("summarized_flow_data.csv")

d<-read.csv("landscape_brick_breakdown_master_2.csv")

#format dates
d$date.deployed<-as.Date(strptime(d$date.deployed,format="%Y-%m-%d"))
d$date.collected<-as.Date(strptime(d$date.collected,format="%Y-%m-%d"))


#merge data
d<-merge(d,flow.sum,by="stream")


## percent mass remaining and breakdown rates
d$percent.remain<-d$final.wt.g/d$start.wt.g
d$first.order.decay<-(-log(d$percent.remain))/d$days.deployed

#remove observations where the blocks were physically broken
d2<-d[d$notes!="crumbled",]

######## plot precent remaining and the ratio of breakdown rates

d.s<-summarySE(d2,measurevar = "percent.remain",groupvars = c("fine.coarse"),na.rm=TRUE)


mass_loss_raw<-ggplot(d.s,aes(y=percent.remain,x=fine.coarse))+
  geom_bar(stat="identity",fill=wes_palette("Darjeeling1")[2])+
  geom_errorbar(aes(ymin=percent.remain-ci,ymax=percent.remain+ci),width=0.01)+
  theme_classic()+xlab("")+scale_x_discrete(labels=c("coarse"="Coarse","fine"="Fine"))+
  theme(text = element_text(size=20))+ylab(expression("Percent Mass Remaining"))

mass_remaining_paired<-ggplot(d2,aes(y=percent.remain,x=fine.coarse))+
  geom_point(size=2,color=wes_palette("Darjeeling1")[2])+
  geom_line(aes(group = bag.number),size=1,color=wes_palette("Darjeeling1")[2])+
  theme_classic()+xlab("Bag Type")+
  theme(text = element_text(size=20))+ylab(expression("Percent Mass Remaining"))

tiff(filename="paired_bricks_29oct2021.tiff",units="in",res=800,width=6,height=10,compression="lzw")
plot_grid(mass_loss_raw,mass_remaining_paired,ncol=1,labels="auto",label_x=0.9,label_y=0.9,label_size = 20)
dev.off()


# pairing fine and coarse measurements
fine<-d2[d2$fine.coarse=="fine",]
coarse<-d2[d2$fine.coarse=="coarse",]

fine2<-data.frame(fine.decay=fine$first.order.decay,bag.number=fine$bag.number)
coarse2<-data.frame(stream=coarse$stream,coarse.decay=coarse$first.order.decay,bag.number=coarse$bag.number,mean.flow=coarse$mean.flow)

d3<-merge(coarse2,fine2,by="bag.number")

################ How much faster are rates in the coarse mesh bags

d3$ratio.of.rates<-d3$coarse.decay/d3$fine.decay
d3$fragmentation<-d3$coarse.decay-(d3$fine.decay-d3$coarse.decay)/(log(d3$fine.decay)-log(d3$coarse.decay))


mean.line<-data.frame(y=c(0,9),x=c(mean(d3$ratio.of.rates,na.rm=TRUE),mean(d3$ratio.of.rates,na.rm=TRUE)))

density_percent_faster<-ggplot(d3, aes(x=ratio.of.rates))+theme_classic() + geom_histogram(fill=wes_palette("Darjeeling1")[2],bins=10)+
  theme(text = element_text(size=20))+ylab(expression("Frequency"))+
  labs(x=expression(paste("Ratio of block breakdown ",italic(k[c]*":"*k[f]))))+
  geom_line(data=mean.line,aes(x=x,y=y),size=2,linetype="dashed")
density_percent_faster

tiff(filename="alt_paired_bricks_29oct2021.tiff",units="in",res=800,width=6,height=10,compression="lzw")
#multiplot(mass_remaining_paired,density_percent_faster,cols=1)
plot_grid(mass_remaining_paired,density_percent_faster,labels="auto",ncol=1,label_x=0.9,label_y=0.9,label_size = 20)
dev.off()


#### Looking at drivers of the differences in breakdown rates


breakdown_summary<-summarySE(d3,measurevar = "fragmentation" , groupvars = "stream",na.rm=TRUE)

mean_coarse<-aggregate(coarse.decay~stream,data=coarse2,mean,na.rm=TRUE)


breakdown_summary2<-merge(breakdown_summary,flow.sum,by="stream")
breakdown_summary3<-merge(breakdown_summary2,mean_coarse,by="stream")

decay.fit<-data.frame(coarse.decay=seq(from=range(d3$coarse.decay,na.rm=TRUE)[1],to=range(d3$coarse.decay,na.rm=TRUE)[2],by=0.00001))
decay.fit$fit<--1.681e-04+3.464e-01*decay.fit$coarse.decay ## values from lmer below

discharge.fit<-data.frame(q=seq(from=range(d3$mean.flow)[1],to=range(d3$mean.flow)[2],by=0.1))
discharge.fit$fit<--1.897e-04+log(discharge.fit$q)*1.542e-04 ## values from lmer below

dischage_percent<-ggplot(breakdown_summary3,aes(x=mean.flow,y=fragmentation))+geom_point()+theme_classic()+
  #geom_smooth(se=FALSE,method="lm")+
  geom_line(data=discharge.fit,aes(x=q,y=fit),size=2,col="#00A08A")+
  ylab(expression("Block fragmentation rate ("*italic(lambda*"F"[B])~day^-1*")"))+
  xlab(expression("Mean discharge (L s"^-1*")"))+
  scale_x_log10()+
  ylim(0,0.001005759)+
  geom_errorbar(aes(ymin=fragmentation-se,ymax=fragmentation+se),width=0.001)+
  theme(text = element_text(size=20))
dischage_percent

rate_percent<-ggplot(d3,aes(x=coarse.decay,y=fragmentation))+geom_point()+theme_classic()+
  #geom_smooth(se=FALSE,method="lm")+
  geom_line(data=decay.fit,aes(x=coarse.decay,y=fit),size=2,col="#00A08A")+
  ylab(expression("Block fragmentation rate ("*italic(lambda*"F"[B])~day^-1*")"))+
  xlab(expression("Block abrasion  ("*italic(k[c])~day^-1*")"))+
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))+
  ylim(0,0.001005759)+
  #scale_y_log10()+scale_x_log10()+
  theme(text = element_text(size=20))
rate_percent

tiff(filename="brick_loss_predictions_29oct2021.tiff",units="in",res=800,width=6,height=11,compression="lzw")
plot_grid(rate_percent,dischage_percent,labels="auto",ncol=1,label_x=0.9,label_y=0.9,label_size = 20)
dev.off()

### Statistics 

mean.flow<-lmer(fragmentation~log(mean.flow)+(1|stream),data=d3)
anova(mean.flow)
r.squaredGLMM(mean.flow) ## r squared is 0.23
summary(mean.flow)

k.percent<-lmer(fragmentation~coarse.decay+(1|stream),data=d3)
anova(k.percent)
summary(k.percent)
r.squaredGLMM(k.percent)

t.test(d3$ratio.of.rates,mu=1)
t.test(d3$fragmentation)

k.percent<-lmer(fragmentation~coarse.decay+(1|stream),data=d3)


################# Litter bag data ####################################################################

leaves2<-read.csv("jan_2018_litterbags.csv")



mean.abrasion<-summarySE(measurevar="first.order.decay",groupvars="stream",data=d[d$fine.coarse=="coarse",],na.rm=TRUE)

leaves2$fragmentation<-leaves2$abs_val_k_coarse-(leaves2$abs_val_k_fine-leaves2$abs_val_k_coarse)/(log(leaves2$abs_val_k_fine)-log(leaves2$abs_val_k_coarse))




leaves3<-summarySE(measurevar="fragmentation",groupvars=c("stream","rhodo_acer"),data=leaves2,na.rm=TRUE)

plot.data<-merge(leaves3,mean.abrasion,by="stream")



abrasion.breakdown.plot<-ggplot(plot.data,aes(x=first.order.decay,y=fragmentation,color=rhodo_acer))+
  theme_classic()+
  geom_smooth(se=FALSE,method="lm")+
  #geom_line(data=decay.fit,aes(x=coarse.decay,y=fit),size=2,col="#00A08A")+
  ylab(expression("Leaf fragmentation rate ("*italic(lambda*"F"[L])~day^-1*")"))+
  xlab(expression("Block abrasion  ("*italic(k[c])~day^-1*")"))+
  #scale_y_log10()+scale_x_log10()+
  theme(text = element_text(size=20))+
  geom_errorbar(aes(ymin=fragmentation-se.x,ymax=fragmentation+se.x),width=0.00001)+
  geom_errorbarh(aes(xmin=first.order.decay-se.y,xmax=first.order.decay+se.y),height=0.00001)+
  theme(legend.position = c(0.15, 0.85))+
  scale_color_manual(name = "", labels = c("Acer", "Rhodo"),values=wes_palette("Darjeeling1")[4:5])+
  geom_point(size=2)
abrasion.breakdown.plot

tiff(filename="leaves.and.abrasion.tiff",units="in",res=800,width=6,height=6,compression="lzw")
abrasion.breakdown.plot
dev.off()

model.leaves.bricks<-lm(fragmentation~first.order.decay+rhodo_acer,data=plot.data)
summary(model.leaves.bricks)
anova(model.leaves.bricks)
