# packages needed
library(readr)
library(ggplot2)
library(ggpubr)
library(survRM2)
library(survival)
library(survminer)
library(gridExtra)
library(grid)
library(lattice)
library(RColorBrewer)
library(dplyr)
library(pvaluefunctions)

# 2 colors in 3 gradients
blue1<-"#1A5276"
blue2<-"#2980B9"
blue3<-"#A9CCE3"
red1<-"#7B241C"
red2<-"#C0392B"
red3<-"#E6B0AA"


##################################################
#### DATA PREPARATION
##################################################

# Read in files
path1<-"/Users/ankerichters/Documents/Werk/IPTW RMST/RCT re-analysis/Data Dogliotti et al/"
path2<-"/Users/ankerichters/Documents/Werk/IPTW RMST/RCT re-analysis/Data Dreicer et al/"
path3<-"/Users/ankerichters/Documents/Werk/IPTW RMST/RCT re-analysis"
dogliotti.arm1<-"Dogliotti Arm 1 GP def.txt"      # Output file for IPD
dogliotti.arm2<-"Dogliotti Arm 2 GC def.txt"      # Output file for IPD
dreicer.arm2<-"Dreicer arm1 def.txt"      # Output file for IPD
dreicer.arm1<-"Dreicer arm2 def.txt"      # Output file for IPD

do.arm1 <- read.table(paste(path1,dogliotti.arm1,sep=""),header=TRUE)
do.arm2 <- read.table(paste(path1,dogliotti.arm2,sep=""),header=TRUE)
dr.arm1 <- read.table(paste(path2,dreicer.arm1,sep=""),header=TRUE)
dr.arm2 <- read.table(paste(path2,dreicer.arm2,sep=""),header=TRUE)

# Merge arms
do.arm1$arm<-1
do.arm2$arm<-2
dr.arm1$arm<-1
dr.arm2$arm<-2
do.total <- rbind(do.arm1,do.arm2)
do.total$study<-1
dr.total <- rbind(dr.arm1,dr.arm2)
dr.total$study<-2
dr.total$Time2<-NULL

# Merge studies
df<-rbind(do.total,dr.total)
df$arm2<-0
df$arm2<-ifelse(df$arm==1 & df$study==1, 1, df$arm2)
df$arm2<-ifelse(df$arm==2 & df$study==1, 2, df$arm2)
df$arm2<-ifelse(df$arm==1 & df$study==2, 3, df$arm2)
df$arm2<-ifelse(df$arm==2 & df$study==2, 4, df$arm2)
df$study <- factor(df$study,
                   labels = c("Dogliotti", "Dreicer"))


##################################################
#### ALTERNATIVE SCENARIO FOR DOGLIOTTI ET AL
##################################################

# Adjusting time and censoring per scenario
#  - Scenario: all patients censored would have died after 4 months
do.total$Event_sx<-ifelse(do.total$Event==0, 1, do.total$Event)        # Turn censored events into death events
do.total$Time_sx<-ifelse(do.total$Event==0, do.total$Time+6, do.total$Time)  # Turn censor times into death times + 6 months

dr.total2<-dr.total
dr.total2$Time_sx<-0
dr.total2$Event_sx<-0
df2<-rbind(do.total,dr.total2)
df2$Event_sx<-ifelse(df2$study==2,df2$Event,df2$Event_sx)
df2$Time_sx<-ifelse(df2$study==2,df2$Time,df2$Time_sx)
#df2$study <- factor(df2$study, labels = c("Dogliotti", "Dreicer"))

dog<-df2[which(df2$study==1),]
dog$Event_sx<-NULL
dog$Time_sx<-NULL
dog$arm<-dog$arm+2

dox<-df2[which(df2$study==1),]
dox$Time<-dox$Time_sx
dox$Event<-dox$Event_sx
dox$Event_sx<-NULL
dox$Time_sx<-NULL
dox$arm<-dox$arm+4

dre<-df2[which(df2$study==2),]
dre$Event_sx<-NULL
dre$Time_sx<-NULL

dflong<-rbind(dre,dog)
dflong<-rbind(dflong,dox)

##################################################
#### MAKE KAPLAN-MEIER GRAPHS
##################################################

KM.dr<-      survfit(Surv(Time, Event)~arm, data=dr.total,type="kaplan-meier")              #1 Dreicer (normal)
KM.do<-      survfit(Surv(Time, Event)~arm, data=do.total,type="kaplan-meier")              #2 Dogliotti (normal)
KM.do.sx<-   survfit(Surv(Time_sx, Event_sx)~arm, data=do.total,type="kaplan-meier")        #3 Dogliotti (alternative)
KM.drdo<-    survfit(Surv(Time, Event)~arm, data=df,type="kaplan-meier")                    #4 Total (studies combined)
KM.drdo.sx<- survfit(Surv(Time_sx, Event_sx)~arm, data=df2,type="kaplan-meier")             #5 Total, altnernative (studies combined)
KM.drdo2<-   survfit(Surv(Time, Event)~arm2, data=df,type="kaplan-meier")                   #6 Dreicer vs Dogliotti 
KM.all<-     survfit(Surv(Time, Event)~arm, data=dflong,type="kaplan-meier")                #7 Total, all studies combined

plot1<-ggsurvplot(KM.dr, title="A. Dreicer et al. (2004)",legend=c(0.2,0.2), font.title=c(10,"bold"),font.x=c(10), font.y=c(10),font.tickslab=c(8),pval=F,
                  legend.labs=c("Cisplatin","Carboplatin"),xlim=c(0,39.5),xlab="Months",ylab="Overall survival",linetype=c("solid","solid"),
                  legend.title="Treatment arm", break.x.by=5,break.y.by=0.2,censor.shape=124, palette=c(blue2, red2), size=0.5, censor.size=3)
plot2<-ggsurvplot(KM.do, title="B. Dogliotti et al. (2007)",font.title=c(10,"bold"),font.x=c(10), font.y=c(10),font.tickslab=c(8),linetype=c("solid","solid"),pval=F,
                  legend=c(0.2,0.2), legend.labs=c("Cisplatin","Carboplatin"),xlim=c(0,15.5),xlab="Months",ylab=NULL,
                  legend.title="Treatment arm", break.x.by=2,break.y.by=0.2,censor.shape=124, palette=c(blue2, red2), size=0.5, censor.size=3)
plot3<-ggsurvplot(KM.do.sx, title="C. Dogliotti - alternative censoring",font.title=c(10,"bold"),font.x=c(10), font.y=c(10),font.tickslab=c(8),pval=F,
                  legend=c(0.2,0.2), legend.labs=c("Cisplatin","Carboplatin"),xlim=c(0,15.5),xlab="Months",ylab=NULL,linetype=c("solid","solid"),
                  legend.title="Treatment arm", break.x.by=2,break.y.by=0.2,censor.shape=124, palette=c(blue2, red2), size=0.5, censor.size=3)
plot4<-ggsurvplot(KM.drdo2, title="D. Overlay of A & B",font.title=c(10,"bold"),font.x=c(10), font.y=c(10),font.tickslab=c(8),
                  legend=c(0.2,0.35), legend.labs=c("Cisplatin (Dogliotti)","Carboplatin (Dogliotti)",
                                                    "Cisplatin (Dreicer)","Carboplatin (Dreicer)"),
                  xlim=c(0,15.5),xlab="Months",ylab="Overall survival",
                  legend.title="Treatment arm", break.x.by=2,break.y.by=0.2,censor.shape=124, 
                  palette=c(blue1,red1,blue3,red3),
                  linetype=c("solid","solid","solid","solid"), 
                  size=0.5, censor.size=3)


plot5<-ggsurvplot(KM.drdo, title="E. Dreicer & Dogliotti, combined",font.title=c(10,"bold"),font.x=c(10), font.y=c(10),font.tickslab=c(8),pval=F,
                  legend=c(0.2,0.2), legend.labs=c("Cisplatin","Carboplatin"),xlim=c(0,15.5),xlab="Months",ylab="Overall survival",
                  legend.title="Treatment arm", break.x.by=2,break.y.by=0.2,censor.shape=124, palette=c(blue2, red2), size=0.5, censor.size=3)
plot6<-ggsurvplot(KM.drdo.sx, title="F. Dreicer & Dogliotti (alt. cens.), combined",font.title=c(10,"bold"),font.x=c(10), font.y=c(10),font.tickslab=c(8),pval=F,
                  legend=c(0.2,0.2), legend.labs=c("Cisplatin","Carboplatin"),xlim=c(0,15.5),xlab="Months",ylab=NULL,
                  legend.title="Treatment arm", break.x.by=2,break.y.by=0.2,censor.shape=124, palette=c(blue2, red2), size=0.5, censor.size=3)

plot1$plot<-plot1$plot+theme(legend.title = element_text(size=8),
  legend.text = element_text(size = 8),
  legend.key = element_rect(colour = "transparent"),
  legend.key.height = unit(0.15, "cm"),
  legend.justification = c(0,0), legend.position = c(0,0),
  legend.background=element_blank())
plot2$plot<-plot2$plot+theme(legend.title = element_text(size=8),
  legend.text = element_text(size = 8),
  legend.key = element_rect(colour = "transparent"),
  legend.key.height = unit(0.15, "cm"),
  legend.justification = c(0,0), legend.position = c(0,0),
  legend.background=element_blank())
plot3$plot<-plot3$plot+theme(legend.title = element_text(size=8),
  legend.text = element_text(size = 8),
  legend.key = element_rect(colour = "transparent"),
  legend.key.height = unit(0.2, "cm"),
  legend.justification = c(0,0), legend.position = c(0,0),
  legend.background=element_blank())
plot4$plot<-plot4$plot+theme(legend.title = element_text(size=8),
  legend.text = element_text(size = 8),
  legend.key = element_rect(colour = "transparent"),
  legend.key.height = unit(0.4, "cm"),
  legend.justification = c(0,0), legend.position = c(0,0),
  legend.background=element_blank())

plot5$plot<-plot5$plot+
  theme(legend.title = element_text(size=8),
    legend.text = element_text(size = 8),
    legend.key = element_rect(colour = "transparent"),
    legend.key.height = unit(0.15, "cm"),
    legend.justification = c(0,0), legend.position = c(0,0),
    legend.background=element_blank())+
    geom_vline(xintercept=6,lty=2,color="gray")+
    geom_vline(xintercept=12,lty=2,color="gray")
plot6$plot<-plot6$plot+
          theme(legend.title = element_text(size=8),
              legend.text = element_text(size = 8),
              legend.key = element_rect(colour = "transparent"),
              legend.key.height = unit(0.15, "cm"),
              legend.justification = c(0,0), legend.position = c(0,0),
              legend.background=element_blank())+
          geom_vline(xintercept=6,lty=2,color="gray")+
          geom_vline(xintercept=12,lty=2,color="gray")

#PLOT 1
lay3<-rbind(c(1,1,2,2,3,3),
            c(1,1,2,2,3,3),
            c(4,4,4,4,4,4),
            c(4,4,4,4,4,4),
            c(4,4,4,4,4,4),
            c(4,4,4,4,4,4),
            c(4,4,4,4,4,4),
            c(5,5,5,6,6,6),
            c(5,5,5,6,6,6))
t <- textGrob(" ")
plotA3<-grid.arrange(plot1$plot,plot2$plot,plot3$plot,plot4$plot,plot5$plot,plot6$plot, layout_matrix=lay3)
ggsave(filename="Plot1_rev.jpg", plot=plotA3, path=path3, width=20, height=26, units=c("cm"), dpi=300, limitsize=F)






#### RESTRICTED MEAN SURVIVAL TIMES 
df2$armx<-df2$arm-1
df$armx<-df$arm-1
do.total$arm<-do.total$arm-1
dr.total$arm<-dr.total$arm-1


#DREICER
rmst2(dr.total$Time, dr.total$Event, dr.total$arm, tau=6)
rmst2(dr.total$Time, dr.total$Event, dr.total$arm, tau=12)

#DOGLIOTTI
rmst2(do.total$Time, do.total$Event, do.total$arm, tau=6)
rmst2(do.total$Time, do.total$Event, do.total$arm, tau=12)

#Actual data
rmst2(df$Time, df$Event, df$armx, tau=6)
rmst2(df$Time, df$Event, df$armx, tau=12)
# Alt cens data
rmst2(df2$Time_sx, df2$Event_sx, df2$armx, tau=6)
rmst2(df2$Time_sx, df2$Event_sx, df2$armx, tasu=12)


surv_median(KM.dr)

