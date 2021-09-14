
## Load packages
library(lme4)
library(ggplot2)
library(lmerTest)
library(Rmisc)
library(reshape2)
library(MuMIn)
library(cowplot)

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

d.s<-summarySE(d,measurevar = "percent.remain",groupvars = c("fine.coarse"),na.rm=TRUE)


mass_loss_raw<-ggplot(d.s,aes(y=percent.remain,x=fine.coarse))+
  geom_bar(stat="identity",fill=wes_palette("Darjeeling1")[2])+
  geom_errorbar(aes(ymin=percent.remain-ci,ymax=percent.remain+ci),width=0.01)+
  theme_classic()+xlab("")+scale_x_discrete(labels=c("coarse"="Coarse","fine"="Fine"))+
  theme(text = element_text(size=20))+ylab(expression("Percent Mass Remaining"))

mass_remaining_paired<-ggplot(d,aes(y=percent.remain,x=fine.coarse))+
  geom_point(size=2,color=wes_palette("Darjeeling1")[2])+
  geom_line(aes(group = bag.number),size=1,color=wes_palette("Darjeeling1")[2])+
  theme_classic()+xlab("Bag Type")+
  theme(text = element_text(size=20))+ylab(expression("Percent Mass Remaining"))

tiff(filename="paired_bricks_5aug2021.tiff",units="in",res=800,width=6,height=10,compression="lzw")
multiplot(mass_loss_raw,mass_remaining_paired,cols=1)
dev.off()


# pairing fine and coarse measurements
fine<-d2[d2$fine.coarse=="fine",]
coarse<-d2[d2$fine.coarse=="coarse",]

fine2<-data.frame(fine.decay=fine$first.order.decay,bag.number=fine$bag.number)
coarse2<-data.frame(stream=coarse$stream,coarse.decay=coarse$first.order.decay,bag.number=coarse$bag.number,mean.flow=coarse$mean.flow)

d3<-merge(coarse2,fine2,by="bag.number")

################ How much faster are rates in the coarse mesh bags

d3$percent_faster<-(d3$coarse.decay-d3$fine.decay)/d3$fine.decay

mean.line<-data.frame(y=c(0,9),x=c(mean(d3$percent_faster,na.rm=TRUE),mean(d3$percent_faster,na.rm=TRUE)))

density_percent_faster<-ggplot(d3, aes(x=percent_faster))+theme_classic() + geom_histogram(fill=wes_palette("Darjeeling1")[2],bins=10)+
  theme(text = element_text(size=20))+ylab(expression("Frequency"))+xlab("Proportion faster in Coarse Mesh")+
  geom_line(data=mean.line,aes(x=x,y=y),size=2,linetype="dashed")
density_percent_faster

tiff(filename="alt_paired_bricks_5aug2021.tiff",units="in",res=800,width=6,height=10,compression="lzw")
#multiplot(mass_remaining_paired,density_percent_faster,cols=1)
plot_grid(mass_remaining_paired,density_percent_faster,labels="AUTO",cols=1,label_x=0.9,label_y=0.9,label_size = 20)
dev.off()


#### Looking at drivers of the differences in breakdown rates


breakdown_summary<-summarySE(d3,measurevar = "percent_faster" , groupvars = "stream",na.rm=TRUE)

mean_coarse<-aggregate(coarse.decay~stream,data=coarse2,mean,na.rm=TRUE)


breakdown_summary2<-merge(breakdown_summary,flow.sum,by="stream")
breakdown_summary3<-merge(breakdown_summary2,mean_coarse,by="stream")

decay.fit<-data.frame(coarse.decay=seq(from=range(d3$coarse.decay,na.rm=TRUE)[1],to=range(d3$coarse.decay,na.rm=TRUE)[2],by=0.00001))
decay.fit$fit<-0.1343+362.3953*decay.fit$coarse.decay ## values from lmer below

discharge.fit<-data.frame(q=seq(from=range(d3$mean.flow)[1],to=range(d3$mean.flow)[2],by=0.1))
discharge.fit$fit<-0.14488+log(discharge.fit$q)*0.14982 ## values from lmer below

dischage_percent<-ggplot(breakdown_summary3,aes(x=mean.flow,y=percent_faster))+geom_point()+theme_classic()+
  #geom_smooth(se=FALSE,method="lm")+
  geom_line(data=discharge.fit,aes(x=q,y=fit),size=2,col="#00A08A")+
  ylab("Proportion faster in Coarse Mesh")+xlab("Mean Discharge (L/s)")+
  scale_x_log10()+geom_errorbar(aes(ymin=percent_faster-se,ymax=percent_faster+se),width=0.001)+
  theme(text = element_text(size=20))+scale_y_log10()
dischage_percent

rate_percent<-ggplot(d3,aes(x=coarse.decay,y=percent_faster))+geom_point()+theme_classic()+
  #geom_smooth(se=FALSE,method="lm")+
  geom_line(data=decay.fit,aes(x=coarse.decay,y=fit),size=2,col="#00A08A")+
  ylab("Proportion faster in Coarse Mesh")+xlab(expression("Coarse Mesh Decay Rate ("*day^-1*")"))+
  scale_y_log10()+scale_x_log10()+
  theme(text = element_text(size=20))
rate_percent

tiff(filename="brick_loss_predictions_5aug2021.tiff",units="in",res=800,width=6,height=10,compression="lzw")
plot_grid(rate_percent,dischage_percent,labels="AUTO",cols=1,label_x=0.9,label_y=0.9,label_size = 20)
dev.off()

### Statistics 

mean.flow<-lmer(percent_faster~log(mean.flow)+(1|stream),data=d3)
anova(mean.flow)
r.squaredGLMM(mean.flow) ## r squared is 0.23
summary(mean.flow)

k.percent<-lmer(percent_faster~coarse.decay+(1|stream),data=d3)
anova(k.percent)
summary(k.percent)
r.squaredGLMM(k.percent)

t.test(d3$percent_faster)
