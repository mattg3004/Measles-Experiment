library("prob")
library("gplots")
library("ggplot2")
library("RColorBrewer")
source("Simpler_measles.R")

demographic.ages              =       read.csv("Ages.csv")
contacts                      <-      read.csv("Contacts.csv")

num.comps                     =       5          # number of compartments in the model
inf.comp                      =       3          # what compartment keeps track of the number of infectious for each age group
t                             =       0
R_0                           =       15
time.step                     =       7          # number of days to step forward
vacc.prop                     =       1
vacc.success                  =       0.9
v                             =       vacc.prop*vacc.success

birth.rate                    =       44.3/(1000*365)
death.rate                    =       14/(1000*365)
prob.survival                 =       1-(death.rate)*(time.step)
initial.prop.susceptible      =       0
days.per.extra.vaccination    =       3 * 365     # do an additonal vaccination campaign every this many days
number.sup.vacs               =       0           # number of supplementary vaccination campaigns so far
max.age.for.supp.vac          =       5           # max age for supplementary vaccinations
supp.vac                      =       0.7         # proportion of people who are up to the age of the max age for supplementary vaccination who will move to vaccinated class
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
av.migrants.per.age.per.day   =       100/(365 * length(demographic.ages[,1]))
prob.survival                 =       1 - death.rate*time.step
updated.state                 =       matrix(0,num.comps*length(demographic.ages[,1]),1)
num.steps                     =       100000
infecteds.by.time             =       matrix(0,num.steps,1)
susceptibles.by.time          =       matrix(0,num.steps,1)
average.infection.age         =       matrix(0,num.steps,1)
foi.by.time                   =       matrix(0,num.steps,1)
all.infectious.by.time        =       matrix(0,num.steps,1)
all.exposed.by.time           =       matrix(0,num.steps,1)
# make transition matrix for all ages
for(j in 1:num.steps){
  N                           =       sum(disease.state)
  births.average              =       birth.rate*N*time.step
  births.total                =       rpois(1,births.average)
  births.vac                  =       rbinom(1,births.total,v)              # number of childrn born who go straight to immune class
  disease.state[1]            =       disease.state[1] + (births.total - births.vac)
  disease.state[5]            =       disease.state[5] + births.vac
  migrant.infecteds           =       rpois(  length(demographic.ages[,1])  ,  av.migrants.per.age.per.day*time.step)
  #print(paste("infected.migrants =",sum(migrant.infecteds)))
  updated.state               =       matrix( 0  ,  num.comps*length(demographic.ages[,1])  ,  1)
  
  disease.state[seq(inf.comp,length(disease.state),num.comps)]    =     disease.state[seq( inf.comp  ,  length(disease.state)  ,  num.comps )] + migrant.infecteds
  # print(paste("infecteds =",sum(disease.state[seq(inf.comp,length(updated.state),num.comps)])))
  pre.infected                =       sum(disease.state[seq( inf.comp  ,  length(disease.state)  ,  num.comps )])
  
  new.infected                =       0
  num.spreader                =       0
  all.exposed                 =       0
  
  for ( i in 1:length(demographic.ages[,1]) ){
    
    age                         =      demographic.ages[i,1]
    
    foi.age                     =      force.of.infection.by.age(age  ,  mixing.matrix  ,  disease.state  ,  beta  ,  gamma  ,  time.step  , inf.comp  ,  num.comps)
   # print(paste("age =",age,",foi =",foi.age))
   if(age == 0)
   {
     foi.by.time[j]      =      foi.age
   }
    change.matrix.by.age        =      stochastic.transmission.matrix2(age , ceiling(disease.state*prob.survival) , v , foi.age , demographic.ages , time.step  ,sigma  ,  rho)
    
    #change.matrix.by.age        =      deterministic.transmission.matrix(age , ceiling(disease.state*prob.survival) , v , foi.age , demographic.ages , time.step  ,sigma  ,  rho)
    new.infected                =      new.infected   +    change.matrix.by.age[11]
    num.spreader                =      num.spreader   +    change.matrix.by.age[12]
    all.exposed                 =      all.exposed    +    change.matrix.by.age[13]
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
  number.by.age                 =       number.of.each.age(demographic.ages,updated.state)
  prop.infected.by.age          =       updated.state[seq(inf.comp,length(updated.state),num.comps)]/number.by.age
  prop.susceptible              =       updated.state[seq(1,length(updated.state),num.comps)]/number.by.age
  total.infecteds               =       sum(updated.state[seq(inf.comp,length(updated.state),num.comps)])
  
 # print(paste("total.infecteds =",total.infecteds ))
 # print(paste("total.exposed =",sum(updated.state[seq(inf.comp-1,length(updated.state),num.comps)])))
 # print(paste("new.infecteds =",new.infected))
  average.infection.age[j]      =       sum(updated.state[seq(inf.comp,length(updated.state),num.comps)]*demographic.ages[,1])
  infecteds.by.time[j]          =       new.infected
  all.infectious.by.time[j]     =       num.spreader
 all.exposed.by.time[j]         =       all.exposed
  susceptibles.by.time[j]       =       sum(updated.state[seq(1,length(updated.state),num.comps)])
  disease.state                 =       updated.state
           t                    =       t   +    time.step
 
 if( t > days.per.extra.vaccination * (number.sup.vacs +1)){
   number.sup.vacs              =       number.sup.vacs   +   1
   disease.state                =       supplementary.vacc(disease.state  ,  supp.vac  ,  max.age.for.supp.vac)
 }
   
   
  if (j %% 50 == 0)
  {
    print(paste("j =",j))
    # Sys.sleep(0.01)
    par(mfrow=c(3,3))
    plot(1:j,infecteds.by.time[1:j],type="l",ylab = "new inf")
    plot(1:j,cumsum(infecteds.by.time[1:j]),type="l",ylab = "cum inf")
    plot(1:j,foi.by.time[1:j],type="l",ylab = "foi.for.0")
    
    plot(demographic.ages[,1],prop.infected.by.age,type="l",ylab = "inf prop",xlab = "age")
    plot(1:j,susceptibles.by.time[1:j],type="b",ylab = "num sus",xlab = "step")
    plot(demographic.ages[,1],updated.state[seq(1,length(updated.state),num.comps)],type="b",ylab = "num sus",xlab = "age")
    
    plot(1:j,all.infectious.by.time[1:j],type="l",ylab = "all inf")
    plot(1:j,all.exposed.by.time[1:j],type="l",ylab = "all exp")
    years   =   ceiling(t/365)
    if (years > 0)
    {
      infections.per.year        =    matrix(0,years,1)
      for(pp in 1:years)
      {
        infections.per.year[pp]   =   sum(infecteds.by.time[((floor(pp-1)*(365/time.step)) +1):floor(pp*365/time.step)])
      }
      plot(1:years,infections.per.year,type="b",ylab = "num infs",xlab = "year")
    }
    
  }
}

View(disease.state)
