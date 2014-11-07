# Script to accompany manuscript titled 'Pathway-level models to predict non-indigenous species establishment using propagule pressure, environmental tolerance and trait data'

# Published in Journal of Applied Ecology
# Authors: Johanna Bradie and Brian Leung
# Updated May 5, 2014

rm(list=ls())     # Clear R memory

setwd('HERE')  ## Put working directory HERE in quotations-- folder must contain appropriate data files

pathway_data<-read.csv("FILENAME.csv") #csv of species names, PP data, species traits/tolerances

trait_columns<-c(HERE) #Put column numbers containing numeric trait data separated by a comma HERE i.e. "trait_columns<-c(5,6,7,8,9,10)"
trait_data<-pathway_data[,trait_columns]

cor(trait_data)  ## View correlation structure between variables.  
trait_data<-trait_data[,-c(HERE)]  # Remove highly correlated variables (>0.5) by indicating column numbers HERE in brackets separated by a comma

## This replaces NA trait values with mean values ##
for (i in 1:length(trait_data[1,]))
{
  NA_vals<-which(is.na(trait_data[,i]==TRUE))
  col_mean<-mean(trait_data[,i],na.rm=TRUE)
  trait_data[NA_vals,i]<-col_mean
}

PP_Est_columns<-c(HERE) # Put column numbers HERE that contain the PP and Est data separated by a comma # i.e. "PP_Est_columns<-c(3,4)" #1ST column must be PP, 2nd column must be Establishment data
PP_Est<-pathway_data[,PP_Est_columns]
colnames(PP_Est)<-c("PP","Est")
  
PP_Est<-as.matrix(PP_Est) # Runs faster as matrix in optimization function
trait_data<-as.matrix(trait_data)

### Likelihood function 2nd order terms for q parameters #############

q_fun2<-function(base_par,params,params2,traits)
  {
    ## Calculate q based on traits and parameters
    z<-base_par+sum(params*traits,na.rm=TRUE)+sum(params2*traits^2,na.rm=TRUE)
    q<-1/(1+exp(-z))
    return(q)
  }

lifish2<-function (pars,sp,tr,give_pos=T)  
  {
    if(pars[1]<=0)
    {
      return(NaN)
    }
    
    par_c<-pars[1]
    ntraits<-length(tr[1,])
    li=0
    
    for (i in 1:length(sp[,1]))
    { 
      if(sp[i,"Est"]==1) {
        li=li+log(1-(q_fun2(base_par=pars[2],params=pars[3:(ntraits+2)],params2=pars[(ntraits+3):(2*ntraits+2)],traits=tr[i,])^((sp[i,"PP"])^par_c))) 
        
      } else {
        li=li+log(q_fun2(base_par=pars[2],params=pars[3:(ntraits+2)],params2=pars[(ntraits+3):(2*ntraits+2)],traits=tr[i,])^((sp[i,"PP"])^par_c)) 

      }
    }
    
    if(give_pos) {
      return(li)
    } else {
      return(-li)
    }
  }

#############################################################################

ntraits<-length(trait_data[1,])
initial_par_values<-c(.5,1,rep(0,1,2*ntraits))  # These can be edited if the optimization will not run.

opt_val<-optim(par=initial_par_values,lifish2,sp=PP_Est,tr=trait_data,give_pos=F,control=list(maxit=100000))
## Check opt_val output to ensure convergence = 0.  If convergence is not 0, you may need to remove collinear variables and/or increase the number of iterations (maxit).

good_parameters<-opt_val$par
c_value=good_parameters[1]
base_par=good_parameters[2]
params=good_parameters[3:(ntraits+2)]
params2=good_parameters[(ntraits+3):((2*ntraits)+2)]

#########################################################
## Generate list of q values.  q is the complement of p; the probability of a single propagule of a given species establishing
#########################################################

q_list<-NULL           
for (i in 1:length(trait_data[,1]))
  {
    q_list<-c(q_list,q_fun2(base_par=base_par,params=params,params2=params2,traits=trait_data[i,]))
  }

pathway_data<-cbind(PP_Est,trait_data,q_list)
pathway_data<-as.array(pathway_data)

########################################################
## Generate a list of probabilities of establishment
########################################################

PE_list<-NULL
PE_list<-1-pathway_data[,which(colnames(pathway_data)=="q_list")]^pathway_data[,which(colnames(pathway_data)=="PP")]^good_parameters[1]
pathway_data<-cbind(pathway_data,PE_list)


### Display top 20 species by q (i.e. per propagule establishment probability)

temp<-pathway_data[order(pathway_data[,which(colnames(pathway_data)=="q_list")],decreasing=FALSE),]
temp[1:20,]

## Display top 20 species by PE

temp<-species_data[order(species_data$PE_list,decreasing=TRUE),]
temp[1:20,]

## Display a histogram of species' q values

hist(q_list,xlab="Species' q value",main="Histogram of species' q values",col=2)

###### Plot graph with all species curves and historical establishment data ###########

plot(log(pathway_data[,which(colnames(pathway_data)=="PP")]),pathway_data[,which(colnames(pathway_data)=="Est")],xlab="log(Imports)",ylab="Probability of Establishment",col="red")
for (i in 1:length(trait_data[,1]))
{
  curve((1-(q_list[i]^exp(x)^good_parameters[1])),add=T, lty=2,col="blue") 
}

## To calculate change in establishment risk for a species if propagule pressure is reduced, use formula Probability of Establishment = 1 - q^N^c   Simply substitute species' q value and predicted c parameter into the model and change the propagule pressure as required.
