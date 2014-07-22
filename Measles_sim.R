library("prob")
library("gplots")
library("ggplot2")
library("RColorBrewer")
source("Measles_func.R")

demographic.ages              =       read.csv("Ages.csv")
contacts                      <-      read.csv("Contacts.csv")

num.comps                     =       6          # number of compartments in the model
inf.comp                      =       4          # what compartment keeps track of the number of infectious for each age group
t                             =       0
R_0                           =       15
time.step                     =       1          # number of days to step forward
average.age.at.vaccination    =       0.5
vacc.prop                     =       1
vacc.success                  =       0.9
vacc.immune                   =       vacc.prop*vacc.success
v                             =       vacc.success*min(1,time.step/60)*min(11/12,(1+time.step/60)/6)
maternal.immunity.loss        =       min(1,time.step/250)
birth.rate                    =       44.3/(1000*365)
death.rate                    =       14/(1000*365)
prob.survival                 =       1-(death.rate)*(time.step)
initial.prop.susceptible      =       0.01
years.per.extra.vaccination   =       3
gamma                         =       0.97
exposed.days                  =       12          # number of days spent in the exposed class on average
infectious.days               =       8           # number of days spent in the infected class on average
sigma                         =       min(1,time.step/exposed.days)
rho                           =       min(1,time.step/infectious.days)
max.age                       =       10          # assume that when R_0 was calculated previously, this was the maximum age of the people spreading measles

mixing.matrix                 <-      full.mixing.matrix(contacts,demographic.ages)      # average number of people met of each age grooup, stratified by age
mean.number.contacts          <-      average.contacts.by.age(mixing.matrix,max.age)
beta_0                        =       R_0/(mean.number.contacts*infectious.days)         # average number of people infected during infectiousness = R_0, 
                                                                                         # therefore mean transmission rate is given by this expression
beta_1                        =       0.8
beta                          =       beta_0 * ( 1 + beta_1 * cos( 2 * pi * t) )
disease.state                 <-      initial.disease.state(demographic.ages  ,  vacc.immune  ,  maternal.immunity.loss  ,  average.age.at.vaccination  ,  initial.prop.susceptible  ,  num.comps)

av.migrants.per.age.per.day   =       1/(365)
prob.survival                 =       1 - death.rate*time.step
updated.state                 =       matrix(0,num.comps*length(demographic.ages[,1]),1)
num.steps                     =       1000
infecteds.by.time             =       matrix(0,num.steps,1)
susceptibles.by.time          =       matrix(0,num.steps,1)
average.infection.age         =       matrix(0,num.steps,1)

# make transition matrix for all ages
for(j in 1:num.steps){
  N                           =       sum(disease.state)
  average.births              =       birth.rate*N*time.step
  total.births                =       rpois(1,average.births)
  disease.state[1]            =       disease.state[1] + total.births
  migrant.infecteds           =       rpois(  length(demographic.ages[,1])  ,  av.migrants.per.age.per.day*time.step)
  print(paste("infected.migrants =",sum(migrant.infecteds)))
  updated.state               =       matrix( 0  ,  num.comps*length(demographic.ages[,1])  ,  1)
  
  disease.state[seq(inf.comp,length(disease.state),num.comps)]    =     disease.state[seq( inf.comp  ,  length(disease.state)  ,  num.comps )] + migrant.infecteds
  print(paste("infecteds =",sum(disease.state[seq(inf.comp,length(updated.state),num.comps)])))
  pre.infected                =       sum(disease.state[seq( inf.comp  ,  length(disease.state)  ,  num.comps )])
 # print(paste("pop size pre =",sum(disease.state)))
#  print(paste("pop size post should be =",ceiling(sum(disease.state)*prob.survival)))
    for ( i in 1:length(demographic.ages[,1]) ){
      
      age                         =      demographic.ages[i,1]

      age.population              =      sum(disease.state[((num.comps*(age))+1):(num.comps*(age+1))])
      foi.age                     =      force.of.infection.by.age(demographic.ages[i,1]  ,  mixing.matrix  ,  disease.state  ,  beta  ,  gamma  ,  time.step  , inf.comp  ,  num.comps)
      #print(foi.age)
      change.matrix.by.age        =      stochastic.transmission.matrix(age , ceiling(disease.state*prob.survival) , v , foi.age , maternal.immunity.loss , demographic.ages , time.step  ,sigma  ,  rho)
      
      
      if(age == max(demographic.ages[,1])){
        
        updated.state[seq((((i-1)*num.comps)+1),num.comps*(i))]        <-      change.matrix.by.age[1:num.comps]
      }
  
      else{
        
        updated.state[seq((((i-1)*num.comps)+1),num.comps*(i+1))]      <-      updated.state[seq((((i-1)*num.comps)+1),num.comps*(i+1))] + change.matrix.by.age
      }
    
    }
 # print(paste("pop size post =",sum(updated.state)))
  total.infecteds               =       sum(updated.state[seq(inf.comp,length(updated.state),num.comps)])
  print(paste("total.infecteds =",total.infecteds))
  print(paste("total.exposed =",sum(updated.state[seq(inf.comp-1,length(updated.state),num.comps)])))
  average.infection.age[j]      =       sum(updated.state[seq(inf.comp,length(updated.state),num.comps)]*demographic.ages[,1])
  infecteds.by.time[j]          =       total.infecteds
  susceptibles.by.time[j]       =       sum(updated.state[seq(2,length(updated.state),num.comps)])
  disease.state                 =       updated.state
        t                       =       t   +    time.step
}

