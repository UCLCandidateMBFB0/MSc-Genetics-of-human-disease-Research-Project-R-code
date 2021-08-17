### This script describes the simulation model used in the analysis##


##****generating result t vs. number of mutations.
library(reshape2)
library(ggplot2)
library(ggstream)
library(ggthemes)

##The illustration of parameters and variables is shown in the dissertation
divide=function(m,s,md,mn){
  d1<-m
  d2<-m
  a=runif(1,0,1)
  #replicate the cell into two daughter cells
  if (m$status==0){
    #0 is stem cell, 1 is differentiated cell
    if (a>s){
      k<-c(1,0)
      l=sample(0:1,1)
      u=k[k!=l]
      d1$status=l
      d2$status=u
    }
  }
  d1$nd<-d1$nd+rpois(1,md/2) #introduce kd ∼ Pois(md /2) driver mutations
  d1$nn<-d1$nn+rpois(1,mn/2) #introduce kn ∼ Pois(mn/2) neutral mutations
  d2$nd<-d2$nd+rpois(1,md/2)
  d2$nn<-d2$nn+rpois(1,mn/2)
  A=rbind(d1,d2)
  return(A)
}

#x: a variable, decrease of the cell death probability per driver mutation
test<-function(md,mn,f,x,g0,d0,dd0,s,pc,tmax){
  t=0;
  p=1;
  md=md;mn=mn;f=f;x=x;g0=g0;d0=d0;dd0=dd0;s=s;pc=pc;tmax=tmax;
  stem_cell_i<-as.data.frame(matrix(NA,1,5))
  colnames(stem_cell_i)<-c("nd","nn","status","pc","t")
  stem_cell_i[1,1]=as.numeric(0)
  stem_cell_i[1,2]=as.numeric(0)
  stem_cell_i[1,3]=as.numeric(0)
  stem_cell_i[1,4]=Inf
  stem_cell_i[1,5]=t
  history_cell<-stem_cell_i
  stem_cell<-stem_cell_i
  while ((t<tmax)){
    print(t);
    #print(p);
    divide_x<-NULL
    temporary_store_gen<-stem_cell_i[-1,];
    
    for (i in 1:nrow(stem_cell)) {
      #print(i);
      g=(1-p/pc)*g0*f^(stem_cell[i,1]);
      d=d0*x^(-(stem_cell[i,1]));
      a=runif(1,0,1);
      if (stem_cell[i,]$status==1){
        d=(dd0)*(x^(-(stem_cell[i,1])))}
      if(a<0.5){
        if(a<g){
          divide_x<-divide(stem_cell[i,],s,md,mn);
        }
        if((p>1) && (a<d)){
          divide_x<-divide_x[-(sample(1:2,1)),];
        }
      }
      else{
        if((p>1) && (a<d)){
          stem_cell<-stem_cell[-i,]
          #print(happened)
        }
        if (a<g){
          divide_x<-divide(stem_cell[i,],s,md,mn);
        }
      }
      print(divide_x)
      temporary_store_gen<-rbind(temporary_store_gen,divide_x) #*********
      row.names(temporary_store_gen)<-NULL;
    } 
    t=t+1
    if(nrow(temporary_store_gen)==0){
      stem_cell$t<-t
      p<-nrow(stem_cell)
      history_cell<-rbind(history_cell,stem_cell)
    }
    else{  
      p<-nrow(temporary_store_gen);
      temporary_store_gen$t<-t
      stem_cell<-temporary_store_gen
      history_cell<-rbind(history_cell,temporary_store_gen);
      #print(history_cell)
    }
  }
  return(list(stem_cell,history_cell,p,t)) 
  #stem cell: the final result of the last division
  #history_cell: the result for the whole period
  #p: number of cells after the final iteration
  #t: total generation
}

## g0 is changed for trail 2 (dormant condition), g0 remains unchanged for trail 1 (cycling condition)
trail1<-test(md=3,mn=10,f=10,x=1,g0=10,d0=0,dd0=0,s= 10^(-1),pc=1/0,tmax=15)
trail2<-test(md=3,mn=10,f=10,x=1,g0=10^(0),d0=0,dd0=0,s= 10^(-1),pc=1/0,tmax=15)

trail2<-test(md=3,mn=10,f=10,x=1,g0=10^(-0.5),d0=0,dd0=0,s= 10^(-1),pc=1/0,tmax=15)

trail2<-test(md=3,mn=10,f=10,x=1,g0=10^(-1),d0=0,dd0=0,s= 10^(-1),pc=1/0,tmax=15)

trail2<-test(md=3,mn=10,f=10,x=1,g0=10^(-1.3),d0=0,dd0=0,s= 10^(-1),pc=1/0,tmax=15)

trail2<-test(md=3,mn=10,f=10,x=1,g0=10^(-1.6),d0=0,dd0=0,s= 10^(-1),pc=1/0,tmax=15)

data1<-trail1[[2]]
data2<-trail2[[2]]
data1$key<-"fastcycling"
data2$key<-"highlyquiescent"

a<-length(unique(data1$t))
data_1_1<-as.data.frame(matrix(NA,a,4))
colnames(data_1_1)<-c("nd","nn","sum","t")
data_1_1$t<-unique(data1$t)
for (i in 1:(nrow(data_1_1))) {
  data_1_1$nd[i]<-sum(data1[data1$t==(i-1),]$nd);
  data_1_1$nn[i]<-sum(data1[data1$t==(i-1),]$nn);
  data_1_1$sum[i]<-sum(data1[data1$t==(i-1),]$nn)+sum(data1[data1$t==(i-1),]$nd);
}
data_1_1$type<-"cycling"

b<-length(unique(data2$t))
data_2_1<-as.data.frame(matrix(NA,b,4))
colnames(data_2_1)<-c("nd","nn","sum","t")
data_2_1$t<-unique(data2$t)
for (i in 1:nrow(data_2_1)) {
  data_2_1$nd[i]<-sum(data2[data2$t==(i-1),]$nd);
  data_2_1$nn[i]<-sum(data2[data2$t==(i-1),]$nn);
  data_2_1$sum[i]<-sum(data2[data2$t==(i-1),]$nn)+sum(data2[data2$t==(i-1),]$nd);
}
data_2_1$type<-"quiescent"
data<-rbind(data_1_1,data_2_1)

data2<-data
t<-max(data2$t)
datamelt<-melt(data2[,-c(3,4)])
datamelt$t<-rep(0:t,4)

# Give to the neutral mutations negative values.
datamelt[datamelt$variable%in%"nn","value"]<--1*datamelt[datamelt$variable%in%"nn","value"]

datamelt$new_status<-paste(datamelt$type,datamelt$variable,sep="_")


##visualisation
A=ggplot(datamelt, aes(x = t, y = value, fill = new_status)) +
  geom_area()+labs(y= "Number of  mutations", x = "generation") + theme_economist() + 
  scale_color_economist()+ggtitle("number of mutation accumulated, 
  g0_dormant=10^(-0.5),g0_cycling=10;tmax=15")

A+ylim(-6e+06,1e+06)
