#Article: "Spatio-Temporal Clustering Analysis Using Generalized Lasso with an Application to Reveal the Spread of Covid-19 Cases in Japan"
#Authors: Septian RAHARDIANTORO* & Wataru SAKAMOTO
#Corresponding Email*: septianrahardiantoro@apps.ipb.ac.id 
#Last Modified: 2023-01-14

#####Code of Section 5.	Real case data application: Weekly Covid-19 cases in Japan

#####Data preparation
dt <- read.csv("20230114 Weekly Covid-19 Cases in Japan.csv",header=T)
dim(dt) #47 81

##Weekly Covid-19 Cases
dt1 <- dt[,4:81]
dim(dt1) #47 78

##Response variable y_it = ln(Weekly Covid-19 Cases/Population)
dt3 <- log((dt1+1)/dt$Pop)

##Preparing data matrix
period <- rep(names(dt1),each=47)
pref <- rep(dt$Prefecture,78)

cases.y2 <- dt3[,1]
for(i in 2:78)
	{
		cases.y2 <- c(cases.y2,dt3[,i])
	} 
length(cases.y2) #3666
data.ok <- data.frame(pref,period,cases.y2)


#####D penalty matrix preparation
library(genlasso);library(sp);library(spdep);library(igraph)
library(rgdal);library(pracma)

##Importing D for fused lasso (D1) based on National Statistics Center (2016)
D <- read.csv("20230114 D land bto.csv",header=T)
D <- as.matrix(D[,-1]) #D1
dim(D) #93 47

##Creating D for temporal effect (D2 for trend filtering)
row <- 76
D1 <- diag(rep(1,row))
dim(D1) #76 76
D2 <- (-2)*D1
D3 <- D1
nol <- matrix(0,row,1)
D11 <- cbind(D1,nol,nol)
D21 <- cbind(nol,D2,nol)
D31 <- cbind(nol,nol,D3)
D.period <- D11+D21+D31 #D2
dim(D.period) #76 78

#####Generalized Lasso Preparation
#---------------------------------
#####Temporal effect detection
##Response variable preparation
y2.tbar <- c()
for(m in 1:78){
	id <- ((m-1)*47+1):(m*47)
	m.y2 <- mean(data.ok$cases.y2[id])
	y2.tbar <- c(y2.tbar,m.y2)
	}
#y2.tbar

##Predictor matrix preparation
Xa <- diag(rep(1,length(y2.tbar))) 
#---------------------------------

#---------------------------------
#####Spatial effect detection
##Response variable preparation
y2.it <- data.ok$cases.y2
##Predictor matrix preparation
Xb1 <- diag(rep(1,47))
#---------------------------------

##Define y
y.tbar <- y2.tbar
yit <- y2.it

#####Generalized Lasso: Temporal Effects Detection
##Setting maximum degree of freedom
df1.max <- 78*3/4

##Trend filtering
beta.tbar0 <- 0

  zt <- y.tbar-beta.tbar0
  mod1 <- genlasso(zt,Xa,D=D.period)
  
  ### Start of selecting lambda_T ###
  df.T <- mod1$df[mod1$df<=df1.max]
  lam1 <- mod1$lambda[mod1$df<=df1.max]
  lam1 <- exp(seq(log(min(lam1)), log(max(lam1)), length=101))  
  alpha <- coef(mod1, lambda=lam1, type="both") 
  beta1 <- alpha$beta 
  u1 <- alpha$u
  df1 <- alpha$df
  
  ####ALOCV & GCV
  alot <- c()
  gt <- c()
  for(j in 1:length(lam1)){
    E <- which(abs( abs(u1[,j])-lam1[j] )<0.00001)
    print(E)
    if (length(E) == 0) {Dbaru <- D.period
	  } else {Dbaru <- D.period[-E,]}
    if (length(E) == nrow(D.period)-1) Dbaru <- matrix(Dbaru,1,ncol(D.period)) 
    A <- Xa %*% nullspace(Dbaru)
    Aplus <- MASS::ginv(A)
    H <- A%*%Aplus
    hj <- diag(H)
    trH <- sum(hj)/78
    hj <- ifelse(hj >=0.999,0.999,hj)
    alo <- mean(((zt-beta1[,j])/(1-hj))^2)
    alot <- c(alot,alo)
    go <- mean(((zt-beta1[,j])/(1-trH))^2)
    gt <- c(gt,go)
  }
  #ALOCV
  best.lam1 <- lam1[which.min(alot)]
  best.df1 <- df1[which.min(alot)]
  #GCV
  best.lam2 <- lam1[which.min(gt)]
  best.df2 <- df1[which.min(gt)]

##Table 3
Method <- c("GenLasso with ALOCV","GenLasso with GCV")
Lambda_T <- c(best.lam1,best.lam2)
df <- c(best.df1,best.df2)
table3 <- data.frame(Method,Lambda_T,df)
table3
#               Method   Lambda_T df
#1 GenLasso with ALOCV 0.27593490 37
#2   GenLasso with GCV 0.04109102 51

##Temporal effects estimates: aplha_t hat
#ALOCV  
alpha.hat <- coef(mod1, lambda=best.lam1)$beta
#GCV
alpha.hat2 <- coef(mod1, lambda=best.lam2)$beta

##Figure 6
plot(alpha.hat,xlab="Week",type="o", col="red") #ALOCV
lines(alpha.hat2,type="o",col="blue") #GCV
#######################################################################################

#####Generalized Lasso: Spatial Effects Detection
##Setting maximum degree of freedom
df2.max <- 47*3/4 #for spatial effects detection

#####Fused Lasso for a Common lambda_S
##Variable preparation
rit <- yit-rep(y.tbar,each=47)
l.time <- length(alpha.hat)

##Checking lambda_S and DF for each t
c.lam <- c()
c.df <- c()
for(z in 1:l.time){
	ix <- ((z-1)*47+1):(z*47)
	md <- genlasso(rit[ix],Xb1,D=D)
  	mdl <- md$lambda[md$df<=df2.max]
	mdd <- md$df[md$df<=df2.max]
	c.lam <- c(c.lam,mdl)
	c.df	<- c(c.df,mdd)
}
range(c.lam)
#[1] 0.03470326 9.58184508
range(c.df)
#[1]  1 35

##Start of selecting lambda_S using ALOCV & GCV
lam2 <- seq(min(c.lam),max(c.lam),length=101)

beta.k.all <- NULL
alot2.all <- NULL
gt2.all <- NULL
df2 <- NULL

for(t in 1:l.time){
  index <- ((t-1)*47+1):(t*47)
  mod2 <- genlasso(rit[index],Xb1,D=D)
  beta <- coef(mod2, lambda=lam2, type="both") 
  beta2 <- beta$beta
  u2 <- beta$u
  df2 <- rbind(df2,beta$df)

  alot2 <- c()
  gt2 <- c()

  for(k in 1:length(lam2)){
    E2 <- which(abs( abs(u2[,k])-lam2[k] )<0.00001)
    if (length(E2) == 0) {Dbaru2 <- D
	} else {Dbaru2 <- D[-E2,]}
    if (length(E2) == nrow(D)-1) Dbaru2 <- matrix(Dbaru2,1,ncol(D)) 
    A2 <- Xb1 %*% nullspace(Dbaru2)
    Aplus2 <- MASS::ginv(A2)
    H2 <- A2%*%Aplus2
    hj2 <- diag(H2)
    trH2 <- sum(hj2)/47
    hj2 <- ifelse(hj2 >=0.999,0.999,hj2)
    alo2 <- mean(((rit[index]-beta2[,k])/(1-hj2))^2)
    alot2 <- c(alot2,alo2)
    go2 <- mean(((rit[index]-beta2[,k])/(1-trH2))^2)
    gt2 <- c(gt2,go2)
  }
  beta.k.all <- rbind(beta.k.all,beta2)
  alot2.all <- rbind(alot2.all,alot2)
  gt2.all <- rbind(gt2.all,gt2)
}
dfS <- apply(df2,2,mean)

alocv.error <- apply(alot2.all,2,mean) #being checked
ind.lamb <- which.min(alocv.error)
best.lamS <- lam2[ind.lamb]
beta.hatS <- beta.k.all[,ind.lamb]
best.dfS <- dfS[ind.lamb]
  
gcv.error <- apply(gt2.all,2,mean)
ind.lamb2 <- which.min(gcv.error)
best.lamS2 <- lam2[ind.lamb2]
beta.hatS2 <- beta.k.all[,ind.lamb2]
best.dfS2 <- dfS[ind.lamb2]
## End of selecting lambda_S ###

##Table 4
Lambda_S <- c(best.lamS,best.lamS2)
df_Spatial <- c(best.dfS,best.dfS2)
table4 <- data.frame(Method,Lambda_S,df_Spatial)
table4
#               Method Lambda_S df_Spatial
#1 GenLasso with ALOCV 4.139974   1.294872
#2   GenLasso with GCV 0.893946   6.641026

##Checking unstable ALOCV
dt.alocvS <- cbind(lam2,alocv.error)
range(dt.alocvS[(alocv.error>100),1]) #0.03470326 3.47167432

#####Figure 7, 8, 9
##Figure 7 Unpenalized beta
beta.hat.unp <- matrix(yit-rep(y.tbar,each=47),47,78)
pref.names <- dt$Prefecture
rownames(beta.hat.unp) <- pref.names
dd<-seq(as.Date("2020/3/21"), as.Date("2021/9/11"),length.out=78)
colnames(beta.hat.unp)<- dd

library(RColorBrewer);library(gplots);library(heatmaply)
mypalette<-brewer.pal(9,"OrRd")
pheat <- heatmaply(beta.hat.unp,dendrogram = "none",colors=mypalette,fontsize_row = 7, fontsize_col = 7,plot_method="ggplot",
                   labCol=as.character(dd)) #as.character(seq(as.Date("2020/3/21"), as.Date("2021/9/11"),length.out=78)))
pheat 

##Figure 8 ALOCV
beta.hat.alo <- matrix(beta.hatS,47,78)
rownames(beta.hat.alo) <- pref.names
colnames(beta.hat.alo)<- dd

pheat2 <- heatmaply(beta.hat.alo,dendrogram = "none",colors=mypalette,fontsize_row = 7, fontsize_col = 7,plot_method="ggplot",
                   labCol=as.character(dd)) #as.character(seq(as.Date("2020/3/21"), as.Date("2021/9/11"),length.out=78)))
pheat2

##Figure 9 GCV
beta.hat.g <- matrix(beta.hatS2,47,78)
rownames(beta.hat.g) <- pref.names
colnames(beta.hat.g)<- dd

pheat3 <- heatmaply(beta.hat.g,dendrogram = "none",colors=mypalette,fontsize_row = 7, fontsize_col = 7,plot_method="ggplot",
                    labCol=as.character(dd)) #as.character(seq(as.Date("2020/3/21"), as.Date("2021/9/11"),length.out=78)))
pheat3

#####Figure 10
##selected GCV result
pheat4 <- heatmaply(beta.hat.g,dendrogram = "row",colors=mypalette,fontsize_row = 7, fontsize_col = 7,plot_method="ggplot",
                    labCol=as.character(dd)) #as.character(seq(as.Date("2020/3/21"), as.Date("2021/9/11"),length.out=78)))
pheat4

#######################################################################################
###End of analysis