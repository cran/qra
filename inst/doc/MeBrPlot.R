## ----linkfun, message=FALSE---------------------------------------------------
library(qra, quietly=TRUE)

## ----echo=FALSE---------------------------------------------------------------
pkg <- "glmmTMB"
pcheck <- suppressWarnings(requireNamespace(pkg, quietly = TRUE))
if(pcheck) pcheck & packageVersion("glmmTMB") >= "1.1.2"
noglmmTMB <- !pcheck

## ----y1988, fig.width=8, fig.height=2.5,  out.width='100%'--------------------
qra::graphSum(df=qra::codling1988, link="cloglog", logScale=FALSE,
                     dead="dead", tot="total", dosevar="ct", Rep="rep",
                     fitRep=NULL, fitPanel=NULL,
                     byFacet=~Cultivar,
                     maint="1988: Codling moth, MeBr",
                     xlab=expression(bold("CT ")*"(gm.h."*m^{-3}*")"))

## ----y1989, fig.width=6.0, fig.height=2.25, out.width="80%"-------------------
qra::graphSum(df=qra::codling1989, link="cloglog", logScale=FALSE,
                     dead="dead", tot="total", dosevar="ct", Rep="rep",
                     fitRep=NULL, fitPanel=NULL,
                     byFacet=~Cultivar,
                     maint="1989: Codling moth, MeBr",
                     xlab=expression(bold("CT ")*"(gm.h."*m^{-3}*")"))

## ----yr88, fig.width=8.5, fig.height=3.0, out.width="98%"---------------------
library(lattice)
cloglog <- make.link("cloglog")$linkfun
cod88 <- subset(qra::codling1988, dose>0)
cm88 <- subset(qra::codling1988, dose==0)
cmMatch <- match(cod88$cultRep,cm88$cultRep)
cod88$cm <- cm88[cmMatch,'PropDead']
cod88$apobs <- with(cod88, (PropDead-cm)/(1-cm))
xyplot(cloglog(apobs)~ct|Cultivar, groups=rep, 
       data=subset(cod88, apobs>0), layout=c(5,1),
       maint="1988: Codling moth, MeBr")

## ----yr89, fig.width=6.5, fig.height=3.0, out.width="70%"---------------------
cod89 <- subset(qra::codling1989, dose>0)
cm89 <- subset(qra::codling1989, dose==0)
cmMatch <- match(cod89$cultRep,cm89$cultRep)
cod89$cm <- cm89[cmMatch,'PropDead']
cod89$apobs <- with(cod89, (PropDead-cm)/(1-cm))
xyplot(cloglog(apobs)~ct|Cultivar, groups=rep, data=cod89, layout=c(3,1),
       maint="1989: Codling moth, MeBr")

## ----noglmmTMB, echo=FALSE----------------------------------------------------
if (noglmmTMB) {
  message(" This vignette requires a version of `glmmTMB` >= 1.1.2")
  message(" Earlier versions of `glmmTMB` may have issues\n for matching versions of `TMB` and `Matrix`")
knit_exit()
}

## ----cm-----------------------------------------------------------------------
ctl <- glmmTMB::glmmTMBControl(optimizer=optim,
                      optArgs=list(method="BFGS"))
cm88.TMB <- glmmTMB::glmmTMB(cbind(dead,total-dead)~Cultivar,
               family=glmmTMB::betabinomial(link='logit'), 
               data=cm88)
cm88.glm <- glm(cbind(dead,total-dead)~Cultivar, 
               family=quasibinomial(link='logit'), 
               data=cm88)
cm89.TMB <- glmmTMB::glmmTMB(cbind(dead,total-dead)~Cultivar, 
                   family=glmmTMB::betabinomial(link='logit'), 
                   data=cm89)
cm89.glm <- glm(cbind(dead,total-dead)~Cultivar, 
                family=quasibinomial(link='logit'), 
                data=cm89)

## ----varmult------------------------------------------------------------------
mults <- matrix(nrow=2, ncol=6)
dimnames(mults) <- list(c('88','89'),c('phi','nmin','nmax',
                                       'Phimin','Phimax','PhiGLM'))
phi88 <- glmmTMB::sigma(cm88.TMB)
nrange <- range(cm88$total)
Phirange <- 1+phi88^{-1}*(nrange-1)
mults['88',] <- c(phi88, nrange, Phirange, 
                  summary(cm88.glm)$dispersion)
phi89 <- sigma(cm89.TMB)
nrange <- range(cm89$total)
Phirange <- 1+phi89^{-1}*(nrange-1)
mults['89',] <- c(phi89, nrange, Phirange, 
                  summary(cm89.glm)$dispersion)
round(mults,2)

## ----more-on-cm89-------------------------------------------------------------
cm89d.TMB <- update(cm89.TMB, dispformula=~0+Cultivar,    
                   control=ctl)
nranges <- with(cm89,
                sapply(split(total,Cultivar),range))
phi89d <- exp(coef(summary(cm89d.TMB))$disp[,1])
Phiranges <- 1+t(nranges-1)*phi89d^-1
colnames(Phiranges) <- c("Min","Max")
round(Phiranges,2)

## ----cmDF88-cults-------------------------------------------------------------
cults <- unique(cm88$Cultivar)
mat <- matrix(nrow=length(cults),ncol=5)
dimnames(mat) <- list(cults,
                      c("phi","nmin","nmax","PhiMin","PhiMax"))
for(cult in cults){
df <- subset(cm88, cult==Cultivar)
obj <- glmmTMB::glmmTMB(cbind(dead,total-dead)~1, control=ctl,      
               family=glmmTMB::betabinomial(link='logit'), data=df)
phi <- glmmTMB::sigma(obj)
nrange <- range(df$total)
mat[cult,] <- c(phi, nrange, 1+(nrange-1)*phi^-1)
}
round(mat,1)

