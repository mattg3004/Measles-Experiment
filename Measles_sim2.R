library("prob")
library("gplots")
library("ggplot2")
library("RColorBrewer")

source("Measles_func2.R")
#source("Measles_initial_conds.R")


# make transition matrix for all ages
for(j in 1:num.steps){
  N                           =       sum(disease.state)
  average.births              =       birth.rate*N*time.step
  total.births                =       rpois(1,average.births)
  disease.state[1]            =       disease.state[1] + total.births
  migrant.infecteds           =       rpois(  length(demographic.ages[,1])  ,  av.migrants.per.age.per.day*time.step)
  #print(paste("infected.migrants =",sum(migrant.infecteds)))
  updated.state               =       matrix( 0  ,  num.comps*length(demographic.ages[,1])  ,  1)
  
  disease.state[seq(inf.comp,length(disease.state),num.comps)]    =     disease.state[seq( inf.comp  ,  length(disease.state)  ,  num.comps )] + migrant.infecteds
  # print(paste("infecteds =",sum(disease.state[seq(inf.comp,length(updated.state),num.comps)])))
  pre.infected                =       sum(disease.state[seq( inf.comp  ,  length(disease.state)  ,  num.comps )])
  
  new.infected                =       0
  
  for ( i in 1:length(demographic.ages[,1]) ){
    
    age                         =      demographic.ages[i,1]
    
    foi.age                     =      force.of.infection.by.age(demographic.ages[i,1]  ,  mixing.matrix  ,  disease.state  ,  beta  ,  gamma  ,  time.step  , inf.comp  ,  num.comps)
    
    change.matrix.by.age        =      stochastic.transmission.matrix(age , ceiling(disease.state*prob.survival) , v , foi.age , maternal.immunity.loss , demographic.ages , time.step  ,sigma  ,  rho)
    new.infected                =      new.infected   +    change.matrix.by.age[13]
    #print(paste("change matrix output =",change.matrix.by.age))
    if(age == max(demographic.ages[,1])){
      
      updated.state[seq((((i-1)*num.comps)+1),num.comps*(i))]        <-      change.matrix.by.age[1:num.comps]
    }
    
    else{
      
      updated.state[seq((((i-1)*num.comps)+1),num.comps*(i+1))]      <-      updated.state[seq((((i-1)*num.comps)+1),num.comps*(i+1))] + change.matrix.by.age[seq(1,2*num.comps)]
    }
    #print(paste("updated state =",updated.state))
  }
  # print(paste("pop size post =",sum(updated.state)))
  total.infecteds               =       sum(updated.state[seq(inf.comp,length(updated.state),num.comps)])
  #print(paste("total.infecteds =",total.infecteds ))
  #print(paste("total.exposed =",sum(updated.state[seq(inf.comp-1,length(updated.state),num.comps)])))
  #print(paste("new.infecteds =",new.infected))
  average.infection.age[j]      =       sum(updated.state[seq(inf.comp,length(updated.state),num.comps)]*demographic.ages[,1])
  infecteds.by.time[j]          =       new.infected
  susceptibles.by.time[j]       =       sum(updated.state[seq(2,length(updated.state),num.comps)])
  disease.state                 =       updated.state
  t                             =       t   +    time.step
  if (j %% 100 == 0)
  {
    print(paste("j =",j))
   # Sys.sleep(0.01)
  #  qplot(seq(1,num.steps),infecteds.by.time,geom="line")
  #  Sys.sleep(1)
  }
}

