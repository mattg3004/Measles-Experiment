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


deterministic.transmission.matrix<-function(age,vacc.prop,vacc.success,force.of.infection,maternal.immunity.loss){
  vacc.prop.success          =    vacc.prop*vacc.success
  Trans.matrix               =    matrix(0,5,5)
  Trans.matrix[1,1]          =    1
  
  if (age == 0){
    Trans.matrix[1,1]        =    1-maternal.immunity.loss
    Trans.matrix[2,1]        =    maternal.immunity.loss
  }
  Trans.matrix[2,2]          =    (1-force.of.infection)*(1-vacc.prop.success)
  Trans.matrix[3,2]          =    force.of.infection*(1-vacc.prop.success)
  Trans.matrix[5,2]          =    vacc.prop.success
  Trans.matrix[4,3]          =    Trans.matrix[4,4]    =   Trans.matrix[5,5] = 1
  
  return(Trans.matrix)
}



stochastic.transmission.matrix<-function(age , disease.state , vacc.immune , foi , maternal.immunity.loss , demographic.ages , time.step ,  sigma  ,  rho){
  
  age.disease.state         =    disease.state[((age*6)+1):((age+1)*6)]
  number.of.age             =    sum(age.disease.state)
  prob.age.change           =    time.step/365
  change.matrix             =    matrix(0,13,1)
  new.infecteds             =    0
  # Assume that no one 1 or older gets vaccinnated
  if (age > 0)
  {
    vacc.immune      =    0
  }
  u          =    prob.age.change
  d          =    maternal.immunity.loss
  v          =    vacc.immune
  
  mat.immune         =    susceptibles     =     infecteds    =     recovered     =     vac    =   exposed   =    matrix(0,12,1)
  if((age == 0)&(age.disease.state[1] > 0))  
  {
    mat.immune       =    rmultinom(1 , age.disease.state[1] , c((1-u)*(1-d) , (1-u)*d*(1-v) , 0 , 0 , 0 , (1-u)*d*v , 0 , u*(1-v) , 0 , 0 , 0 ,u*v))
  }else{
    mat.immune       =    matrix(0,12,1)
  }
  
  if (age.disease.state[2] > 0){
    if(age == 0 ){
      susceptibles     =    rmultinom(1, age.disease.state[2] , c(0, (1-u)*(1-foi)*(1-v) , (1-u)*foi*(1-v) , 0 , 0 , v*(1-u) , 0 , u*(1-foi)*(1-v) , u*foi*(1-v) , 0 , 0 , u*v))   
    }
    else{
      susceptibles     =    rmultinom(1, age.disease.state[2] , c(0, (1-u)*(1-foi) , (1-u)*foi , 0 , 0 , 0 , 0 , u*(1-foi) , u*foi , 0 , 0 , 0))
    }
    new.infecteds    =    susceptibles[3]    +    susceptibles[9]
  }
  
  if (age.disease.state[3] > 0){
    exposed          =    rmultinom(1, age.disease.state[3] , c(0, 0, (1-u)*(1-sigma) , (1-u)*sigma , 0 , 0 , 0 , 0 , u*(1-sigma) , u*sigma , 0 , 0))
  }
  
  if (age.disease.state[4] > 0){
    #number.recover   =    rbinom(age.disease.state[4],1,min(1,365*time.step/15))
    #prop.recover     =    number.recover/age.disease.state[4]
    infecteds        =    rmultinom(1,age.disease.state[4],c(0 , 0 , 0 , (1-u)*(1-rho) , (1-u)*rho , 0 , 0 , 0 , 0 , u*(1-rho) , u*rho,0))
  }
  
  if (age.disease.state[5] > 0){
    recovered        =    rmultinom(1,age.disease.state[5],c(0 , 0 , 0 , 0 , (1-u) , 0 , 0 , 0 , 0 , 0 , u , 0))
  }
  
  if (age.disease.state[6] > 0){
    vac              =    rmultinom(1,age.disease.state[6],c(0 , 0 , 0 , 0 , 0 , (1-u) , 0 , 0 , 0 , 0 , 0 , u))
  }
  
  change.matrix[seq(1,12)]      =    c(mat.immune[1:6] + susceptibles[1:6] + exposed[1:6] + infecteds[1:6] + recovered[1:6] + vac[1:6],
                            mat.immune[7:12] + susceptibles[7:12] + exposed[7:12] + infecteds[7:12] + recovered[7:12] + vac[7:12])
  change.matrix[13] = new.infecteds
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
  disease.state[1]              =      ceiling(demographic.ages[1,2]*(1-maternal.immunity.loss))
  disease.state[2]              =      ceiling(demographic.ages[1,2]*(maternal.immunity.loss)*((average.age.at.vaccination)*(1-vacc.immune) + (1-average.age.at.vaccination)))
  disease.state[num.comps]      =      ceiling(demographic.ages[1,2]*(maternal.immunity.loss)*average.age.at.vaccination*(vacc.immune))
  
  for (i in 2:length(demographic.ages[,1])){
    disease.state[(((i-1)*num.comps)+1):((i)*num.comps)]      =    ceiling(c(0 , demographic.ages[i,2]*initial.prop.susceptible , 0 , 0 , 0 , demographic.ages[i,2]*(1-initial.prop.susceptible)))
  }
  
  return(disease.state)
}

