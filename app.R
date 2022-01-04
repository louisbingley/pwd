
options(scipen = 999)
options(warn = -1)
bfun=function(x,n,p){
  output=choose(n,x)*(p^x)*((1-p)^(n-x))
  return(output)
};
Bfun=function(xx,n,p){
  
  B0=choose(n,0)*(p^0)*(1-p)^n
  B=rep(NA,xx)
  for (i in 1:xx){
    B[i]=choose(n,i)*(p^i)*((1-p)^(n-i)) 
  }
  output=B0+sum(B)
  return(output)
};
beta1fun=function(theta0,k,N0,N1,d1,d2){theta1=theta0+d1;theta2=theta0+d2;

innersum=rep(NA,N1-N0+1) #innersum is a vector storing each individual innersum value

for(x in N0:N1){
  
  f0=bfun(x,N1,theta2)*1*(Bfun(x-1,N1,theta1))^(k-1)
  f=rep(NA,k-1)
  for (j in 1:(k-1)){
    f[j]=choose(k-1,j)*
      bfun(x,N1,theta2)*
      (bfun(x,N1,theta1))^j*
      (Bfun(x-1,N1,theta1))^(k-1-j)*
      1/(j+1) }
  
  innersum[x-(N0-1)]=f0+sum(f) }

output=sum(innersum)
return(output)

}
n2funn=function(theta0,a,BETA2,xi1,xi2,d1,d2){ #my function n2
  if(BETA2<1){
    BETA2<-BETA2;
    if(1-a<1){a<- a} else{ a<-0};
    N2=(qnorm(1-a)*sqrt((xi1*theta0+xi2*(theta0+d2))*(1-xi1*theta0-xi2*(theta0+d2))/(xi1*xi2))+qnorm(BETA2)*sqrt(theta0*(1-theta0)/xi1+(theta0+d2)*(1-theta0-d2)/xi2)  )^2/d2^2
    #  N2=ceiling(N2)
  }
  else{BETA2<-1;
  if(1-a<1){ a<-a} else{ a<-0};
  N2=(qnorm(1-a)*sqrt((xi1*theta0+xi2*(theta0+d2))*(1-xi1*theta0-xi2*(theta0+d2))/(xi1*xi2))+qnorm(BETA2)*sqrt(theta0*(1-theta0)/xi1+(theta0+d2)*(1-theta0-d2)/xi2)  )^2/d2^2
  # N2=ceiling(N2)
  }
  return(N2)
}


# for fixed n1, run n0=1:n1, locate the min(EN/Nmax) for the given n1 
# pet0_low=.7;pet0_up=.9;a2_low=0.08;a2_up=0.15;
f_inner=function(theta0,d1,d2,k,n1,alphalv,beta,xi1,xi2,dz,om,ao_low,ao_up,a2_low,a2_up){ 
  ## define empty vec
  {  pi0=rep(NA,n1)
  pi1=rep(NA,n1)
  oneminuspi0=rep(NA,n1)
  beta1=rep(NA,n1)
  beta2min=rep(NA,n1)
  lambda2=rep(NA,n1)
  N2=rep(NA,n1); alpha2max=rep(NA,n1); overalla=rep(NA,n1); EN=rep(NA,n1); E0N=rep(NA,n1);Nmax=rep(NA,n1);Nmax_show=rep(NA,n1)
  }  
  ## for fixed n1, run n0=1:n1, locate the min(EN/Nmax) for the given n1
  for(n0 in 1:n1){
    beta1[n0]=beta1fun(theta0,k,n0,n1,d1,d2)                                        ;beta1[n0]=round(beta1[n0],4)
    beta2min[n0]=beta/beta1[n0]                                                     ;beta2min[n0]=round( beta2min[n0],4) 
    pi0[n0]=1-(Bfun(n0-1,n1,theta0))^k                                              ;pi0[n0]=round(pi0[n0],4)
    theta1=theta0+d1;theta2=theta0+d2 
    pi1[n0]=1-(Bfun(n0-1,n1,theta1))^(k-1)*Bfun(n0-1,n1,theta2)                     ; pi1[n0]=round(pi1[n0],4)
    oneminuspi0[n0]=1-pi0[n0]
    
    if(dz=='s2'){ 
      
      if (ao_low<ao_up ){
        pet0_low=1- ao_up/alphalv   
        pet0_up=1- ao_low/alphalv
        alpha2max[n0]=alphalv; overalla[n0]=alpha2max[n0]*pi0[n0]
        N2[n0]=n2funn(theta0,alpha2max[n0],beta2min[n0],xi1,xi2,d1,d2)           
        EN[n0]=k*n1+N2[n0]*(pi0[n0]+pi1[n0])/2
        E0N[n0]=k*n1+N2[n0]*pi0[n0]
        Nmax[n0]=k*n1+N2[n0]
        Nmax_show[n0]=k*n1+ceiling(N2[n0]*xi1)+ceiling(N2[n0]*xi2)
        overallbeta=beta
        n0=seq(1,n1,1) 
        tb=cbind.data.frame(k,theta0,overalla,overallbeta,n1,n0,beta1,beta2min, alpha2max,
                            ceiling(N2*xi1),ceiling(N2*xi2),   EN,   E0N,  Nmax, Nmax_show,oneminuspi0,pi1);  
        if(  om=='o'){ tb$Nb =tb$EN } else 
           if(om=='m'){ tb $Nb =tb$Nmax }else
             {stop("Please enter valid om or dz.")}
        therows=tb[which(tb$Nb>0 & tb$oneminuspi0<=pet0_up & tb$oneminuspi0>=pet0_low),]
      }else {stop("The lower bound of Overall Type I error should be less than its upper bound.") }
      
      
    }else if(dz=='overall'){ 
      
      if(a2_low<a2_up){
        alpha2max[n0]=alphalv/pi0[n0]; overalla[n0]=alpha2max[n0]*pi0[n0]
        N2[n0]=n2funn(theta0,alpha2max[n0],beta2min[n0],xi1,xi2,d1,d2)           
        EN[n0]=k*n1+N2[n0]*(pi0[n0]+pi1[n0])/2
        E0N[n0]=k*n1+N2[n0]*pi0[n0]
        Nmax[n0]=k*n1+N2[n0]
        Nmax_show[n0]=k*n1+ceiling(N2[n0]*xi1)+ceiling(N2[n0]*xi2)
        overallbeta=beta
        n0=seq(1,n1,1) 
        tb=cbind.data.frame(k,theta0,overalla,   overallbeta,n1,n0,beta1,beta2min, alpha2max,   
                            ceiling(N2*xi1),ceiling(N2*xi2),   EN,   E0N,  Nmax, Nmax_show,oneminuspi0,pi1);  
        if(  om=='o'){ tb$Nb =tb$EN } else 
          if(om=='m'){ tb $Nb =tb$Nmax }else
            {stop("Please enter valid om or dz.")}
        therows=tb[which(tb$Nb>0 & tb$alpha2max<=a2_up &tb$alpha2max>=a2_low),]
      }else {stop("The lower bound of StageII TypeI error should be less than its upper bound.") }
     
    }
  
  }#n0 loop
  
 
  theNb=min(therows$Nb ) 
  therow=tb[which(tb$Nb ==theNb ),];
  return(therow)
  
} 

# define function f: compute each single row in the tables (locate the n1 yielding the min Nbenchmark)
f=function(theta0,d1,d2,k,alphalv,beta,xi1,xi2,dz,om,ao_low,ao_up,a2_low,a2_up){    
#  if(dz=='s2'){
#    if(pet0_low>0 & pet0_low<pet0_up & pet0_up<1){a2_low=0;a2_up=1}else{stop("Please enter valid boundries for PET0!")}
#  }else if(dz=='overall'){
#    if(a2_low>0 & a2_low<a2_up & a2_up<1){pet0_low=0;pet0_up=1}else{stop("Please enter valid boundries for a2!")}
#  }else{pet0_low=NA;pet0_up=NA;a2_low=NA;a2_up=NA}
  
  ## call function step1_locateEN_flex for all n1
  stacked=NULL; length(n1range)
  for (i in 1:length(n1range)){
    coming=f_inner(theta0,d1,d2,k,n1range[i],alphalv,beta,xi1,xi2,dz,om,ao_low,ao_up,a2_low,a2_up) 
    stacked=rbind(stacked,coming) 
  } 
  stacked$y2=qnorm(1-stacked$alpha2max)
  stacked=as.data.frame(stacked)
  #   colnames(stacked)[ncol(stacked)-1]<-"Nb";stacked
  row=stacked[which(stacked$Nb==min(stacked$Nb)),]
  row=row[,c(3,9,4,7,8,5,6,10,11,12,13,15,16,19)]
  return(row)
}



n1range=seq(10,100,1)

############################################
# shiny R                                  #
############################################

library(shiny)
library(shinythemes)
library(data.table)
library(randomForest)

####################################
# User Interface                   #
####################################
ui <- fluidPage(theme = shinytheme("united"),
                navbarPage("Flexible Pick the Winner Design (PWD) Sample Size Calculator",
                           
                           tabPanel("Home",
                                    # Input values
                                    sidebarPanel(HTML("<h3>Input parameters</h3>"),
                                                 selectInput("dz", label = "Type1Error Control Level", 
                                                             choices = list(  "Stage2 level"="s2",   "overall level"="overall"),
                                                             selected = "s2"),
                                                 selectInput("om", label = "Optimization Method", 
                                                             choices = list(  "Optimal"="o", "Minimax"="m"),
                                                             selected = "o"),

                                                 sliderInput("d1", label = HTML("&delta;1"), value = 0.05, min = 0, max = 0.3,step = 0.05),
                                                 sliderInput("d2", label = HTML("&delta;2"), value = 0.2, min = 0,max = 0.6,step = 0.1), 
                                                 sliderInput("theta0", label = HTML("&theta;0"), value = 0.2, min = 0, max = 1,step = 0.05),
                                                 sliderInput("k", label = "K", value = 2, min = 2,max = 5,step = 1),
                                                 sliderInput("alphalv",label=HTML("&alpha; level"),value=0.05,min =0.05 ,max =0.05 ,step =0 ),
                                                 sliderInput("beta",label="Overall Power",value=0.7,min =0.7 ,max =1 ,step =0.05 ),
                                                 sliderInput("xi1",label="Allocation Ratio (Ctl Arm)",value=0.5,min =0.1 ,max =1 ,step =0.1 ),
                                                 sliderInput("xi2",label="Allocation Ratio (Trmt Arm)",value=0.5,min =0.1 ,max =1 ,step =0.1 ),
                                                
                                                 sliderInput("ao_low",label="Please set a lower bound for Overall type I error if Type1Error Control Level = 'Stage2 level' ",value=0.005,min =0.001 ,max =.01 ,step =0.001 ),
                                                 sliderInput("ao_up",label="Please set an upper bound for Overall type I error if Type1Error Control Level = 'Stage2 level' ",value=0.015,min =0.005,max =.015 ,step =0.001 ),
                                                 sliderInput("a2_low",label="Please set a lower bound for Stage2 type I error if Type1Error Control Level = 'Overall level' ",value=0.08,min =0.05 ,max =0.1 ,step =0.01 ),
                                                 sliderInput("a2_up",label="Please set an upper bound for Stage2 type I error if Type1Error Control Level = 'Overall level' ",value=0.15,min =0.1 ,max =.15 ,step =0.01 ),
                                                 
                                                 actionButton("submitbutton", "Submit",class = "btn btn-primary")),
                                    mainPanel(
                                      tags$label(h3('Status/Output')), # Status/Output Text Box
                                      verbatimTextOutput('contents'),
                                      tableOutput('tabledata1'),  # Results table
                                      tableOutput('tabledata2')  # Results table
                                   
                                    ) # mainPanel()
                                    
                           ),#tabpanel 
                           
                           tabPanel("Result Interpretation", 
                                    titlePanel("Result Interpretation"), 
                                    div(includeMarkdown("aboutt.html"), 
                                        align="justify")
                           ) #tabPanel(), About
                           
                           
                ) # navbarPage()
)# fluidPage()

####################################
# Server                           #
####################################
server <- function(input, output, session) {
  
  # Input Data
  datasetInput <- reactive({  
    d1=input$d1  
    d2=input$d2  
    theta0=input$theta0
    k=input$k
    alphalv=input$alphalv
    beta=input$beta
    xi1=input$xi1
    xi2=input$xi2
    ao_low=input$ao_low
    ao_up=input$ao_up
    a2_low=input$a2_low
    a2_up=input$a2_up
    dz=as.character(input$dz)
    om=as.character(input$om)
    
    c <- f(theta0,d1,d2,k,alphalv,beta,xi1,xi2,dz,om,ao_low,ao_up,a2_low,a2_up)
    c <- data.frame(c)
    names(c) <- c("alpha_overall","alpha2","power_overall","beta1","beta2","n1","y1","n2.Ctl","n2.Trmt","EN","EN0","Nmax","PET0","y2")
    print(c)
   # print(c)
  })
  

  
  # Status/Output Text Box
  output$contents <- renderPrint({
    if (input$submitbutton>0) { 
      isolate("Calculation is complete.") 
    } else {
      return("Please enter parameter values for calculation.")
    }
  })
  
  # Prediction results table
  output$tabledata1 <- renderTable({
    if (input$submitbutton>0) {  
      c1 <-isolate(datasetInput()[,1:5])
      } 
  })
  output$tabledata2 <- renderTable({
    if (input$submitbutton>0) {  
      c2 <-isolate(datasetInput()[,6:14])
    } 
  })
}


####################################
# Create Shiny App                 #
####################################
shinyApp(ui = ui, server = server)
