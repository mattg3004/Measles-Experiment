library("prob")
library("gplots")
library("ggplot2")
library("RColorBrewer")
source("Measles_func2.R")

demographic.ages              =       read.csv("Ages.csv")
contacts                      <-      read.csv("Contacts.csv")

num.comps                     =       6          # number of compartments in the model
inf.comp                      =       4          # what compartment keeps track of the number of infectious for each age group
t                             =       0
R_0                           =       15
time.step                     =       14          # number of days to step forward
average.age.at.vaccination    =       0.5
vacc.prop                     =       1
vacc.success                  =       0.9
vacc.immune                   =       vacc.prop*vacc.success
v                             =       vacc.immune
maternal.immunity.loss        =       min(1,time.step/250)
birth.rate                    =       44.3/(1000*365)
death.rate                    =       14/(1000*365)
prob.survival                 =       1-(death.rate)*(time.step)
initial.prop.susceptible      =       0.02
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

av.migrants.per.age.per.day   =       1/(365 * length(demographic.ages[,1])) 
prob.survival                 =       1 - death.rate*time.step
updated.state                 =       matrix(0,num.comps*length(demographic.ages[,1]),1)
num.steps                     =       100000
infecteds.by.time             =       matrix(0,num.steps,1)
susceptibles.by.time          =       matrix(0,num.steps,1)
average.infection.age         =       matrix(0,num.steps,1)
