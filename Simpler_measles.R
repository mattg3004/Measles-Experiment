#

full.mixing.matrix<- function(contacts,demographic.ages)
{
  polymod.ages                =    contacts[,1]
  contacts                    <-   contacts[,2:16]
  max.age                     <-   max(demographic.ages[,1])
  contact.matrix              =    matrix(0,max.age+1,max.age+1)
  colnames(contact.matrix)    =    c(0:max.age)
  rownames(contact.matrix)    =    c(0:max.age)
  num.polymod.ages            <-   length(polymod.ages)
  row.start.position          =    1
  
  for(i in 1:(num.polymod.ages)){
    col.start.position        =    1
    col.end.position          =    1
    
    if (i == num.polymod.ages){
      row.end.position        =    length(contact.matrix[1,])
    }
    else
    {
      row.end.position        =    row.start.position + polymod.ages[i+1] - polymod.ages[i] -1
    }
    
    for (k in 1:(num.polymod.ages)){
      
      if (k < num.polymod.ages){
        num.repeats           =    polymod.ages[k+1] - polymod.ages[k]
        col.end.position      =    col.start.position + num.repeats - 1
      }
      else
      {
        num.repeats           =    length(contact.matrix[1,]) - polymod.ages[k]
        col.end.position      =    length(contact.matrix[1,])
      }
      num.repeats             =    max((row.end.position-row.start.position+1),(col.end.position-col.start.position+1))
      contact.block           =    matrix((contacts[i,k]/num.repeats),(row.end.position-row.start.position+1),(col.end.position-col.start.position+1))
      
      
      contact.matrix[row.start.position  :  row.end.position  ,  col.start.position  :  col.end.position] = contact.block
      
      col.start.position      =    col.end.position + 1
    }
    row.start.position        =    row.end.position + 1   
  }
  
  return(contact.matrix)
}

average.contacts.by.age<- function(mixing.matrix,max.age){
  
  contacts.by.age       =     colSums(mixing.matrix[,seq(1,(max.age+1))])
  av.num.contacts       =     mean(contacts.by.age)
  return(av.num.contacts)
}


deterministic.transmission.matrix<-function(age , disease.state , vacc.immune , foi  , demographic.ages , time.step ,  sigma  ,  rho){
  
  age.disease.state         =    disease.state[((age*5)+1):((age+1)*5)]
  number.of.age             =    sum(age.disease.state)
  u                         =    time.step/365
  change.matrix             =    matrix(0,11,1)
  new.infecteds             =    0
  
  A                         =    matrix(0,11,1)
  A[1]                      =    round((1-u)*(1-foi) * age.disease.state[1])
  A[2]                      =    round((1-u)*(foi) * age.disease.state[1])   +   round((1-u)*(1-sigma) * age.disease.state[2])
  A[3]                      =    round((1-u)*(sigma) * age.disease.state[2])   +   round((1-u)*(1-rho) * age.disease.state[3])
  A[4]                      =    round((1-u)*(rho) * age.disease.state[3])   +   round((1-u) * age.disease.state[4])
  A[5]                      =    round((1-u)   *  age.disease.state[5])
  A[6]                      =    round((u)*(1-foi) * age.disease.state[1])
  A[7]                      =    round((u)*(foi) * age.disease.state[1])   +   round((u)*(1-sigma) * age.disease.state[2])
  A[8]                      =    round((u)*(sigma) * age.disease.state[2])   +   round((u)*(1-rho) * age.disease.state[3])
  A[9]                      =    round((u)*(rho) * age.disease.state[3])   +   round((u) * age.disease.state[4])
  A[10]                     =    round((u)   *  age.disease.state[5])
  
  A[11]                     =    round((1-u)*(foi) * age.disease.state[1])   +   round((u)*(foi) * age.disease.state[1])
  
  return(A)
}


number.of.each.age <- function(demographic.ages,disease.state){
  
  number.of.age = matrix(0,length(demographic.ages[,1]),1)
  for (i in 1:length(demographic.ages[,1]))
  {
    age = demographic.ages[i,1]
    number.of.age[i]          =    sum(disease.state[((age*5)+1):((age+1)*5)])
  }
  return(number.of.age)
}


stochastic.transmission.matrix2<-function(age , disease.state , vacc.immune , foi  , demographic.ages , time.step ,  sigma  ,  rho){
  
  age.disease.state         =    disease.state[((age*5)+1):((age+1)*5)]
  number.of.age             =    sum(age.disease.state)
  prob.age.change           =    time.step/365
  change.matrix             =    matrix(0,13,1)
  new.infecteds             =    0
  num.spreaders             =    0
  all.exposed               =    0
  # Assume that no one 1 or older gets vaccinnated
  if (age > 0)
  {
    vacc.immune      =    0
  }
  u          =    prob.age.change
  v          =    vacc.immune
  
  susceptibles     =     infecteds    =     recovered     =     vac    =   exposed   =    matrix(0,12,1)
  
  if (age.disease.state[1] > 0){
      susceptibles     =    rmultinom(1, age.disease.state[1] , c( (1-u)*(1-foi) , (1-u)*foi , 0 , 0 , 0 , u*(1-foi) , u*foi , 0 , 0 , 0))
      new.infecteds    =    susceptibles[2]    +    susceptibles[7]
  }
  
  if (age.disease.state[2] > 0){
    exposed          =    rmultinom(1, age.disease.state[2] , c(0 , (1-u)*(1-sigma), (1-u)*sigma , 0 , 0 , 0 , u*(1-sigma) , u*sigma , 0 , 0))
    num.spreaders    =    exposed[3]     +     exposed[8]
    all.exposed      =    new.infecteds   +  exposed[2]   +   exposed[7]
  }
  
  if (age.disease.state[3] > 0){
    #number.recover   =    rbinom(age.disease.state[4],1,min(1,365*time.step/15))
    #prop.recover     =    number.recover/age.disease.state[4]
    infecteds        =    rmultinom(1,age.disease.state[3],c( 0 , 0 , (1-u)*(1-rho) , (1-u)*rho , 0 , 0 , 0 , u*(1-rho) , u*rho,0))
    num.spreaders     =    infecteds[3]   +   infecteds[8]
  }
  
  if (age.disease.state[4] > 0){
    recovered        =    rmultinom(1,age.disease.state[4],c( 0 , 0 , 0 , (1-u) , 0 , 0 , 0 , 0 , u , 0))
  }
  
  if (age.disease.state[5] > 0){
    vac              =    rmultinom(1,age.disease.state[5],c(0 , 0 , 0 , 0 , (1-u) , 0 , 0 , 0 , 0 , u))
  }
  
  change.matrix[seq(1,10)]      =    c( susceptibles[1:5] + exposed[1:5] + infecteds[1:5] + recovered[1:5] + vac[1:5],
                                        susceptibles[6:10] + exposed[6:10] + infecteds[6:10] + recovered[6:10] + vac[6:10])
  change.matrix[11] = new.infecteds
  change.matrix[12] = num.spreaders
  change.matrix[13] = all.exposed
  #print(new.infecteds)
  return(change.matrix)
}

contacts.per.age.group <- function(mixing.matrix  ,  demographic.ages){
  mixing.times.population <- mixing.matrix
  contacts.per.age <- mixing.matrix
  for (i in 1:length(demographic.ages[,1])){
    k<-matrix(demographic.ages[i,2],length(demographic.ages[,1]),1)
    mixing.times.population[,i] <- k*mixing.matrix[,i]
  }
  for (i in 1:length(demographic.ages[,1])){
    contacts.per.age[i,] <- mixing.times.population[i,]/demographic.ages[i,2]
  }
  
  return(rowSums(contacts.per.age))
}

grouped.ages <-function(contacts,demographic.ages)
{
  group.ages <- matrix(0,length(contacts[,1]),2)
  group.ages[,1] <- contacts[,1]
  count = 1
  for (i in 1:length(demographic.ages[,1]))
  {
    if(demographic.ages[i,1] > tail(group.ages[,1],1)){
      group.ages[length(group.ages[,1]),2] <- group.ages[length(group.ages[,1]),2] + demographic.ages[i,2]
    }
    else if(group.ages[count+1,1] > demographic.ages[i,1]){
      group.ages[count,2] <- group.ages[count,2] + demographic.ages[i,2] 
    }
    else {
      count <- count + 1
      group.ages[count,2] <- group.ages[count,2] + demographic.ages[i,2] 
    }
  }
  return(group.ages)
}


force.of.infection.by.age <- function(age  ,  mixing.matrix  ,  disease.state  ,  beta  ,  gamma  ,  time.step  ,  inf.comp  ,  num.comps){
  number.age.brackets        =      length(disease.state)/num.comps
  infectious.indices         =      seq(inf.comp,length(disease.state),num.comps)
  number.infectious.by.age   =      disease.state[infectious.indices]
  mixing.by.age              =      mixing.matrix[,age+1]*time.step
  population.by.age          =      matrix(0,number.age.brackets,1)
  for(q in 1:number.age.brackets){
    population.by.age[q]     =      max(1,sum(disease.state[seq(((num.comps*(q-1))+1),num.comps*q)]))
  }
  foi.by.age                 =      1 - exp(-sum  ( beta* ( (number.infectious.by.age)^gamma )*mixing.by.age/population.by.age  )  )
  
  return(foi.by.age)
}

initial.disease.state <- function(demographic.ages  ,  vacc.immune  ,  maternal.immunity.loss  ,  average.age.at.vaccination  ,  initial.prop.susceptible  ,  num.comps){
  disease.state                 =      matrix(0,num.comps*length(demographic.ages[,1]),1)
  
  for (i in 1:length(demographic.ages[,1])){
    disease.state[(((i-1)*num.comps)+1):((i)*num.comps)]      =    ceiling(c( demographic.ages[i,2]*initial.prop.susceptible , 0 , 0 , 0, demographic.ages[i,2]*(1-initial.prop.susceptible)))
  }
  
  return(disease.state)
}

supplementary.vacc <- function(disease.state,supp.vac,max.age){
  disease.state
  for (age in 0:max.age)
  {
    age.disease.state                 =    disease.state[((age*5)+1):((age+1)*5)]
    p = disease.state[((age*5)+1)]
    l                                 =    round(disease.state[((age*5)+1)]*(1-supp.vac))
    disease.state[((age*5)+1)]        =    l
    disease.state[(age+1)*5]          =    round(disease.state[((age*5)+1)]*(supp.vac)) + disease.state[((age+1)*5)]
  }
  return(disease.state)
}
  