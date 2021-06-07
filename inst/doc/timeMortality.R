## ----setup, include=FALSE-------------------------------------------
knitr::opts_chunk$set(echo = FALSE, comment=NA, width=70)
options(show.signif.stars=FALSE)

## ----prepareData, echo=TRUE-----------------------------------------
suppressPackageStartupMessages(library(qra))
HawCon <- qra::HawCon
## Change name "CommonName" to "CN", for more compact output.
CCnum <- match("CommonName", names(HawCon))
names(HawCon)[CCnum] <- "CN"
## trtGp will identify species & lifestage combination
## trtGpRep will identify species, lifestage, and rep
## cTime is centered version of TrtTime
## scTime is centered and scaled version of TrtTime,
## needed to get some mixed model fits to converge
HawCon <- within(HawCon, {
  trtGp <- factor(paste0(CN,LifestageTrt, sep=":"))
  trtGpRep <- paste0(CN,LifestageTrt,":",RepNumber)
  scTime <- scale(TrtTime)
  obs <- factor(1:nrow(HawCon))
})

## ----cap1, echo=FALSE-----------------------------------------------
cap1 <- " Graphs are designed to give an indication of the pattern, 
when mortalities are shown on a complementary log-log scale, of 
mortality response with days in coolstorage."

## ----plots, fig.width=7, fig.height=7.5, fig.align='center', out.width="75%", fig.cap=cap1----
library(ggplot2)
qra::graphSum(df=HawCon, link="cloglog", logScale=FALSE,
         dead="Dead", tot="Total", dosevar="TrtTime",
         Rep="RepNumber", fitRep=NULL, fitPanel=NULL,
         byFacet=~trtGp, layout=LifestageTrt~Species,
         xlab="Days", maint="Hawaian contemporary data")

## ----glmmTMB-altFits, message=FALSE, warning=FALSE, echo=FALSE------
if("VGAM" %in% (.packages()))
  detach("package:VGAM", unload=TRUE)
# Load packages that will be used
suppressMessages(
  {library(lme4); library(glmmTMB); library(splines); 
    library(DHARMa)})
form <- cbind(Dead,Live)~0+trtGp/TrtTime+(1|trtGpRep)
form2s <- cbind(Dead,Live)~0+trtGp/TrtTime+ns(scTime,2)+(1|trtGpRep)
HCbb.cll <- glmmTMB(form, dispformula=~trtGp+ns(scTime,2), 
                    family=betabinomial(link="cloglog"), data=HawCon)
HCbb2s.cll <- update(HCbb.cll, formula=form2s)
HCbb.logit <- glmmTMB(form, dispformula=~trtGp+ns(TrtTime,2), 
                      family=betabinomial(link="logit"),data=HawCon)
HCbb2s.logit <- update(HCbb.logit, formula=form2s)

## ----glmmTMB-altFits, eval=FALSE, echo=-(1:2)-----------------------
#  if("VGAM" %in% (.packages()))
#    detach("package:VGAM", unload=TRUE)
#  # Load packages that will be used
#  suppressMessages(
#    {library(lme4); library(glmmTMB); library(splines);
#      library(DHARMa)})
#  form <- cbind(Dead,Live)~0+trtGp/TrtTime+(1|trtGpRep)
#  form2s <- cbind(Dead,Live)~0+trtGp/TrtTime+ns(scTime,2)+(1|trtGpRep)
#  HCbb.cll <- glmmTMB(form, dispformula=~trtGp+ns(scTime,2),
#                      family=betabinomial(link="cloglog"), data=HawCon)
#  HCbb2s.cll <- update(HCbb.cll, formula=form2s)
#  HCbb.logit <- glmmTMB(form, dispformula=~trtGp+ns(TrtTime,2),
#                        family=betabinomial(link="logit"),data=HawCon)
#  HCbb2s.logit <- update(HCbb.logit, formula=form2s)

## ----cap7, echo=FALSE-----------------------------------------------
cap7 <- "Differences are shown, between fitted degree 2 normal spline 
curves and fittes lines.
Panel A is for the models that use a complementary log-log
(cloglog) link, while Panel B is for a logit link."

## ----glmmTMB-altFits-gph1, fig.width=9, fig.height=4.0, fig.show='hold', top=2, out.width="100%",  fig.align='center', message=FALSE, warning=0, fig.cap=cap7, echo=FALSE----
library(lattice)
clog <- make.link('cloglog')$linkfun
logit <- make.link('logit')$linkfun
cloginv <- make.link('cloglog')$linkinv
logitinv <- make.link('logit')$linkinv
my.panel.bands <- function(x, y, upper, lower, fill, col,
                           subscripts, ..., font, fontface)
{
  upper <- upper[subscripts]
  lower <- lower[subscripts]
  panel.lines(x,lower, ...)
  panel.lines(x,upper, ...)
}
panel2 <- function(x, y, ...){
  panel.superpose(x, y, 
                  ## panel.groups = my.panel.bands, 
                  type='l', ...)
  panel.xyplot(x, y, type='l', lwd=2, cex=0.6, ...)
}
parset <- simpleTheme(col=rep(1:4,rep(2,4)), lty=rep(1:2,4), lwd=2)
## c('solid','1141')
dat <- expand.grid(trtGp=factor(levels(HawCon$trtGp), levels=levels(HawCon$trtGp)),
                   TrtTime=pretty(range(HawCon$TrtTime),15), link=c('cloglog','logit'))
dat$scTime <- scale(dat$TrtTime)
dat$trtGpRep <- rep('new',nrow(dat))
hatClog <- predict(HCbb.cll, newdata=subset(dat, link=='cloglog'), se=TRUE, allow.new.levels=TRUE)
hatClog2 <- predict(HCbb2s.cll, newdata=subset(dat, link=='cloglog'), se=TRUE,  allow.new.levels=TRUE)
diffClog <- cloginv(hatClog2$fit)-cloginv(hatClog$fit)
hatLgt <- predict(HCbb.logit, newdata=subset(dat, link=='logit'), se=TRUE, allow.new.levels=TRUE)
hatLgt2 <- predict(HCbb2s.logit, newdata=subset(dat, link=='logit'), se=TRUE,  allow.new.levels=TRUE)
diffLgt <- logitinv(hatLgt2$fit)-logitinv(hatLgt$fit)
dat <- within(dat, {diff<-c(diffClog,diffLgt)
})
## dat <- within(dat, {lwr<-fit-2*se.fit; upr<-fit+2*se.fit})

gph <- xyplot(diff~TrtTime|link, outer=TRUE, data=dat, groups=trtGp,
              # upper = dat$upr, lower = dat$lwr, 
              panel = panel2,  
              xlab="Days in coolstorage", ylab="Difference in fitted value", 
              auto.key=list(text=levels(HawCon$trtGp), columns=4, points=FALSE, lines=TRUE), 
              par.settings=parset, layout=c(2,1), 
              scales=list(x=list(at=c(2,6,10,14)), 
              y=list(relation='free'), alternating=c(1,1)))
gph2 <- update(gph, strip=strip.custom(factor.levels=
            c("A: Complementary log-log link", "B: Logit link")))
# parset3 <- simpleTheme(col=rep(1:4,rep(2,4)), lty=rep(2, 8), lwd=rep(2:1,4), alpha=0.35)
# gph3 <- xyplot(fit2~TrtTime|link, outer=TRUE, data=dat, groups=trtGp,
#               ## upper = dat$upr, lower = dat$lwr, 
#              panel = panel2, par.settings=parset3, layout=c(2,1))
gph2 

## ----cap3, echo=FALSE-----------------------------------------------
cap3 <- "Panel A shows the quantile-quantile plot,
for the linear model with a complementary log-log link.
Panel B plots estimated quantiles against mortality, 
while Panel C plots estimated quantiles against total 
number, on a logarithmic scale."

## ----DHARMa, fig.width=3.75, fig.asp=0.95, bot=-1, top=1.5, out.width="32%", warning=0, fig.align='center', fig.show="hold", message=FALSE, echo=FALSE, mar=c(3.1,3.1,2.6,1.1), fig.cap=cap3----
set.seed(29)
simRef <- simulateResiduals(HCbb.cll, n=250, seed=FALSE)
plotQQunif(simRef)
plotResiduals(simRef)
plotResiduals(simRef, form=log(HawCon$Total), xlab="log(Total)")

## ----cap4, echo=FALSE-----------------------------------------------
cap4 <- "Diagnostic plots, for the model with a logit link."

## ----residBYgp, fig.width=3.75, fig.asp=0.95, bot=-1, top=1.5, out.width="32%", warning=0, fig.align='center', fig.show="hold", message=FALSE, echo=FALSE, mar=c(3.1,3.1,2.6,1.1), fig.cap=cap4----
simResBB <- suppressWarnings(simulateResiduals(HCbb.cll, n=250) )
plotQQunif(simResBB)
plotResiduals(simResBB, xlab= "Predictions, Complementary log-log model")
simResLGT <- suppressWarnings(simulateResiduals(HCbb.logit, n=250) )
plotResiduals(simResLGT, xlab= "Predictions, logit model")

## Alternative -- group names are not shown:
## plotResiduals(simRes, form=HawCon$trtGp)

## ----cap4.5, echo=FALSE---------------------------------------------
cap4.5 <- "Quantile residuals, by treatment group, for the 
betabinomial model"

## ----trtGp, echo=FALSE, fig.width=5, fig.asp=0.7, bot=-2.5, warning=FALSE, fig.align='center', fig.show="hold", message=FALSE, warning=FALSE, out.width="49%", echo=FALSE, fig.cap=cap4.5----
dotplot(scaledResiduals~HawCon$trtGp, data=simResBB, 
        scales=list(x=list(rot=30)), ylab="Quantile Residuals",
       main=list(expression(plain("A: Residuals, by treatment group")),
                      x=0.05, y=-0.2, just="left"))
bwplot(scaledResiduals~HawCon$trtGp, data=simResBB,  ylab="",
       scales=list(x=list(rot=30)),
       main=list(expression(plain("B: Boxplot comparison of residuals")),
                 x=0.05, y=-0.2, just="left"))

## ----cap4.75--------------------------------------------------------
cap4.75 <- paste0("Diagnostics for model fitted to strongly overdispersed
binomial type data. Notice that the overdispersion results in an
S-shaped distribution of the residuals around the line $y=x$. The
boxplot is, in this case, uninformative.")

## ----overdispSim, echo=FALSE, fig.width=5, fig.asp=0.85, bot=-2.5, warning=FALSE, fig.align='center', fig.show="hold", out.width="49%", echo=FALSE, fig.cap=cap4.75----
yes <- rbinom(n=100, size=50, prob=0.5)
sim <- cbind(yes=yes*4, no=200-yes*4)
sim.TMB <- glmmTMB(sim~1, family=binomial)
sim250 <- simulateResiduals(sim.TMB)
plotQQunif(sim250)
plotResiduals(sim250)

## ----obslevel, echo=FALSE, eval=TRUE--------------------------------
ctl <- glmmTMBControl(optimizer=optim,
                               optArgs=list(method="BFGS"))
HCbiObs.cll <- 
  glmmTMB(cbind(Dead, Live) ~ 0 + trtGp/scTime +
                       (1|obs) + (scTime|trtGpRep),
                       family=binomial(link='cloglog'),
                       control=ctl, data=HawCon)
HCbiObs.logit <- 
  glmmTMB(cbind(Dead, Live) ~ 0 + trtGp/scTime +
                          (1|obs) + (scTime|trtGpRep),
                      family=binomial(link='logit'),
                      control=ctl, data=HawCon)

## ----glmmTMB-altFits-AIC, echo=FALSE--------------------------------
aicStats <- 
  AIC(HCbb.cll,HCbb2s.cll,HCbb.logit, HCbb2s.logit, 
              HCbiObs.cll, HCbiObs.logit)
rownames(aicStats) <- substring(rownames(aicStats),3)
round(t(aicStats),2)

## ----cap8, echo=FALSE-----------------------------------------------
cap8 <- "Panels A and B show intra-class correlation estimates 
for, respectively, a complementary log-log link and a logit link. 
Both models assume a betabinomial error. Panel C shows, for 
the complementary log-log model, the dispersion factors that 
result."

## ----glmmTMB-altFits-gph-disp, fig.width=9, fig.height=3.5, top=2, out.width="100%",  fig.align='center', fig.show="hold", fig.cap=cap8, echo=FALSE, R.options=par(oma=c(0,0,2,0))----
parset <- simpleTheme(col=rep(1:4,rep(2,4)),pch=rep(c(1,2), 4), lwd=2)
HawCon$rho2clog <- qra::getRho(HCbb.cll)
HawCon$dispClog <- with(HawCon, 1+(Total-1)*rho2clog)
titles=c(expression("A: "*rho*", cloglog link"),expression("B: "*rho*", logit link"),
         "C: Dispersion, cloglog link")
library(lattice)
HawCon$rho2logit <- qra::getRho(HCbb.logit)
xyplot(rho2clog+rho2logit+dispClog ~ TrtTime, groups=trtGp, data=HawCon,
       outer=TRUE, between=list(x=0.25),
       par.settings=parset,
       scales=list(x=list(alternating=FALSE), y=list(relation='free',tick.number=4)),
       auto.key=list(columns=4, between.columns=2, between=1),
       xlab="Days in coolstorage", ylab="Parameter Estimate",
       strip=strip.custom(factor.levels=titles))

## ----obslevel, echo=FALSE, eval=TRUE--------------------------------
ctl <- glmmTMBControl(optimizer=optim,
                               optArgs=list(method="BFGS"))
HCbiObs.cll <- 
  glmmTMB(cbind(Dead, Live) ~ 0 + trtGp/scTime +
                       (1|obs) + (scTime|trtGpRep),
                       family=binomial(link='cloglog'),
                       control=ctl, data=HawCon)
HCbiObs.logit <- 
  glmmTMB(cbind(Dead, Live) ~ 0 + trtGp/scTime +
                          (1|obs) + (scTime|trtGpRep),
                      family=binomial(link='logit'),
                      control=ctl, data=HawCon)

## ----cap2-----------------------------------------------------------
cap2 <- "AICs are compared between models with binomial family
errors and observation level random effects.  The two sets of
points make the comparison for, respectively, data that have 
been simulated from a model with logit link and a model with 
complementary log-log link.  The large triangle makes the
comparison for the models fitted to the 'HawCon' data."

## ----only-obslevel, fig.width=3.75, fig.asp=0.85, bot=-1, top=1.5, out.width="55%",  fig.align='center', fig.show="hold", message=FALSE, echo=FALSE, mar=c(3.1,3.1,2.6,1.1), fig.cap = cap2----
## Notice the use of the 'BFGS' optimization method in place 
## of the default.  
aicData <- setNames(AIC(HCbiObs.logit,HCbiObs.cll)[,2], c('logit','cll'))
set.seed(17)
sim <- simulate(HCbiObs.logit, nsim=10)
simcll <- simulate(HCbiObs.cll, nsim=10)
aic.cll <- aic2.cll <- aic.logit <- aic2.logit <- numeric(10)
for (i in 1:10){
  zlogit <- 
    glmmTMB(sim[[i]] ~ 0 + trtGp/scTime +
                 (1|obs) + (1|trtGpRep),
             family=binomial(link='logit'),
             data=HawCon) 
  aic.logit[i] <- AIC(zlogit)
  zcll <- 
    glmmTMB(sim[[i]] ~ 0 + trtGp/scTime +
                 (1|obs) + (1|trtGpRep),
             family=binomial(link='cloglog'),
             data=HawCon) 
  aic.cll[i] <- AIC(zcll)
  zlogit2 <- 
    glmmTMB(simcll[[i]] ~ 0 + trtGp/scTime +
                 (1|obs) + (1|trtGpRep),
             family=binomial(link='logit'),
             data=HawCon) 
  aic2.logit[i] <- AIC(zlogit2)
  zcll2 <- 
    glmmTMB(simcll[[i]] ~ 0 + trtGp/scTime +
                 (1|obs) + (1|trtGpRep),
             family=binomial(link='cloglog'),
             data=HawCon) 
  aic2.cll[i] <- AIC(zcll2)  
}
gph1 <- xyplot(aic.cll~aic.logit,
               key=list(text=list(c("logit model","cll model", "Data")), 
               points=list(pch=c(1,1,2),
                           col=c('blue','red','black')),columns=3))
gph2 <- xyplot(aic2.cll~aic2.logit, col='red')
gph3 <- gph1+latticeExtra::as.layer(gph2)+
  latticeExtra::layer(panel.abline(0,1), 
          panel.points(aicData[['logit']],
                       aicData[['cll']], pch=2, cex=2, col=1))
update(gph3, xlim=range(c(aic.logit, aic2.logit,
                          aicData['logit']))*c(.98,1.02), 
                        ylim=range(c(aic.cll, aic2.cll, 
                                     aicData['cll']))*c(.98,1.02))

## ----cap5, echo=FALSE-----------------------------------------------
cap5 <- "Diagnostics for model with binomial errors and observation level
random effects."

## ----biObsCLL, fig.width=3.75, fig.asp=0.95, bot=-1, top=2.5, out.width="49%",  fig.align='center', fig.show="hold", message=FALSE, warning=0, echo=FALSE, mar=c(3.1,3.1,2.6,1.1), fig.cap=cap5----
set.seed(29)
simRefcll <- suppressWarnings(
  simulateResiduals(HCbiObs.cll, n=250, seed=FALSE) )
plotResiduals(simRefcll, xlab='cll: model prediction')
plotResiduals(simRefcll, form=log(HawCon$Total), 
              xlab="cll: log(Total)")
simReflogit <- suppressWarnings(
  simulateResiduals(HCbiObs.logit, n=250, seed=FALSE) )
plotResiduals(simReflogit, xlab='logit: model prediction')
plotResiduals(simReflogit, form=log(HawCon$Total), 
              xlab="logit: log(Total)")

## ----biObsCLL, eval=FALSE, echo=TRUE--------------------------------
#  set.seed(29)
#  simRefcll <- suppressWarnings(
#    simulateResiduals(HCbiObs.cll, n=250, seed=FALSE) )
#  plotResiduals(simRefcll, xlab='cll: model prediction')
#  plotResiduals(simRefcll, form=log(HawCon$Total),
#                xlab="cll: log(Total)")
#  simReflogit <- suppressWarnings(
#    simulateResiduals(HCbiObs.logit, n=250, seed=FALSE) )
#  plotResiduals(simReflogit, xlab='logit: model prediction')
#  plotResiduals(simReflogit, form=log(HawCon$Total),
#                xlab="logit: log(Total)")

## ----shorten, echo=FALSE--------------------------------------------
shorten <- function(nam, leaveout=c('trtGp','Fly',':')){
  for(txt in leaveout){
    nam <- gsub(txt,'', nam, fixed=TRUE)
  }
  nam
}

## ----crude-cll, echo=FALSE, warning=F-------------------------------
## Fit two simplistic and unsatisfactory models.
HCbbDisp1.cll <- update(HCbb.cll, dispformula=~1)
HCbin.cll <- update(HCbb.cll, family=binomial(link="cloglog"))

## ----extract-BB-LTcll, echo=FALSE-----------------------------------
LTbb.cll <- qra::extractLT(p=0.99, obj=HCbb.cll, link="cloglog",
                              a=1:8, b=9:16, eps=0, df.t=NULL)[,-2]
rownames(LTbb.cll) <- shorten(rownames(LTbb.cll))

## ----extract-BI-obsRE, echo=FALSE-----------------------------------
## NB: The formula has used the scaled value of time.
## Below, `offset` is used to retrieve the scaling parameters 
## `center` ## and `scale` in `(x-center)/scale`.
offset <- qra::getScaleCoef(HawCon$scTime)
LTbiObs.cll <- qra::extractLT(p=0.99, obj=HCbiObs.cll,
                          a=1:8, b=9:16, eps=0, offset=offset,
                          df.t=NULL)[,-2]
rownames(LTbiObs.cll) <- shorten(rownames(LTbiObs.cll))

## ----extract-BB-LTlogit, echo=FALSE---------------------------------
LTbb.logit <- qra::extractLT(p=0.99, obj=HCbb.logit, link="logit",
                          a=1:8, b=9:16, eps=0, offset=0,
                          df.t=NULL)[,-2]
rownames(LTbb.logit) <- shorten(rownames(LTbb.logit))

## ----extract-crude-LTcll, echo=FALSE--------------------------------
LTbbDisp1.cll <- 
  qra::extractLT(p=0.99, obj=HCbbDisp1.cll, 
                 a=1:8, b=9:16, eps=0, df.t=NULL)[,-2]
rownames(LTbbDisp1.cll) <- shorten(rownames(LTbbDisp1.cll))
LTbin.cll <- 
  qra::extractLT(p=0.99, obj=HCbin.cll, 
                 a=1:8, b=9:16, eps=0, df.t=NULL)[,-2]
rownames(LTbin.cll) <- shorten(rownames(LTbin.cll))

## ----cap9, echo=FALSE-----------------------------------------------
cap9 <- "LT99 $95\\%$ confidence intervals are compared between
five different models."

## ----plotCI, echo=FALSE, fig.width=7.0, bot=1.0, fig.asp=0.625, warning=FALSE, fig.align='center', message=FALSE, out.width="75%", echo=FALSE, fig.cap=cap9----
gpNam <- rownames(LTbb.cll)
ordEst <- order(LTbb.cll[,1])
library(plotrix)
col5 <- c("blue","lightslateblue","brown",'gray','gray50')
plotCI(1:8-0.34, y=LTbb.cll[ordEst,1], ui=LTbb.cll[ordEst,3],
       li=LTbb.cll[ordEst,2], lwd=2, col=col5[1], xaxt="n", 
       fg="gray", xlab="", ylab="LT99 Estimate (days)", 
       xlim=c(0.8,8.2), ylim=c(0,29), sfrac=0.008)
plotCI(1:8-0.17, y=LTbiObs.cll[ordEst,1], ui=LTbiObs.cll[ordEst,3],
       li=LTbiObs.cll[ordEst,2], lwd=2, col=col5[2], xaxt="n", 
       fg="gray", xlab="", ylab="LT99 Estimate (days)", 
       xlim=c(0.8,8.2), ylim=c(0,29), add=TRUE, sfrac=0.008)
plotCI(1:8, y=LTbb.logit[ordEst,1], ui=LTbb.logit[ordEst,3],
       li=LTbb.logit[ordEst,2], lwd=2, col=col5[3], xaxt="n", 
       add=TRUE, sfrac=0.008)
plotCI(1:8+0.17, y=LTbbDisp1.cll[ordEst,1], ui=LTbbDisp1.cll[ordEst,3],
       li=LTbbDisp1.cll[ordEst,2], lwd=2, col=col5[4], xaxt="n", 
       add=TRUE, sfrac=0.008)
plotCI(1:8+0.34, y=LTbin.cll[ordEst,1], ui=LTbin.cll[ordEst,3],
       li=LTbin.cll[ordEst,2], lwd=2, col=col5[5], xaxt="n", 
       add=TRUE, sfrac=0.008)
axis(1, at=1:8, labels=gpNam[ordEst], las=2, lwd=0, 
     lwd.ticks=0.5, col.ticks="gray")
legend("topleft", legend=c("BB-cll (cll=cloglog)",  "BB-cll-ObsRE", "BB-logit",
                           "BB-cll, const Disp factor", 
                           "Binomial-cll"),
       inset=c(0.01,0.01), lty=c(1,1), col=col5[1:5],
       text.col=col5[1:5], bty="n",y.intersp=0.85)

## ----extract-BB-LTcll, eval=FALSE, echo=TRUE------------------------
#  LTbb.cll <- qra::extractLT(p=0.99, obj=HCbb.cll, link="cloglog",
#                                a=1:8, b=9:16, eps=0, df.t=NULL)[,-2]
#  rownames(LTbb.cll) <- shorten(rownames(LTbb.cll))

## ----obsLevel1, echo=FALSE, eval=FALSE------------------------------
#  dMedEgg <- with(HawCon, dummy(trtGp,"MedFlyEgg:"))
#  dMedL1 <- with(HawCon, dummy(trtGp,"MedFlyL1:"))
#  dMedL2 <- with(HawCon, dummy(trtGp,"MedFlyL2:"))
#  dMedL3 <- with(HawCon, dummy(trtGp,"MedFlyL3:"))
#  dMelEgg <- with(HawCon, dummy(trtGp,"MelonFlyEgg:"))
#  dMelL1 <- with(HawCon, dummy(trtGp,"MelonFlyL1:"))
#  dMelL2 <- with(HawCon, dummy(trtGp,"MelonFlyL2:"))
#  dMelL3 <- with(HawCon, dummy(trtGp,"MelonFlyL3:"))
#  formXre <- cbind(Dead, Live) ~ 0 + trtGp/scTime +
#    (1|obs) + (1|trtGpRep) +
#    (0 + dMedEgg| trtGpRep) + (0 + dMedL1 | trtGpRep) +
#    (0 + dMedL2 | trtGpRep) + (0 + dMedL3 | trtGpRep) +
#    (0 + dMelEgg| trtGpRep) + (0 + dMelL1 | trtGpRep) +
#    (0 + dMelL2 | trtGpRep) + (0 + dMelL3 | trtGpRep)
#  ## NB: The formula has used the scaled value of time.
#  ## Below, `offset` is used to record the scaling parameters
#  ## `center` ## and `scale` in `(x-center)/scale`.
#  offset <- with(attributes(HawCon$scTime),
#                 c(`scaled:center`, `scaled:scale`))
#  HCXre.biObs <- glmmTMB(formXre, family=binomial(link='cloglog'),
#                         control=ctl, data=HawCon)
#  ## Could not get the following to converge
#  ## formXreS <- update(formXre, ~ .+ trtGp/splines::ns(scTime,2))

## ----glmerFits------------------------------------------------------
## Comparisons using `glmer()`
HCglmerBIobs.cll <- 
  glmer(cbind(Dead, Live) ~ 0 + trtGp/scTime + (1|obs) + (1|trtGpRep),
        family=binomial(link='cloglog'), nAGQ=1, data=HawCon,
        control=glmerControl(optimizer='bobyqa'))
HCglmerBIobs2s.cll <- suppressWarnings(
  glmer(cbind(Dead, Live) ~ 0 + trtGp/scTime + splines::ns(scTime,2) +
        (1|obs) + (1|trtGpRep), family=binomial(link='cloglog'),
        nAGQ=0, data=HawCon))
HCglmerBIobs.logit <- glmer(cbind(Dead, Live) ~ 0 + trtGp/scTime +
                           (1|obs) + (1|trtGpRep),
        family=binomial(link='logit'), nAGQ=0, data=HawCon, 
        control=glmerControl(optimizer='bobyqa'))
HCglmerBIobs2s.logit <- 
  glmer(cbind(Dead, Live) ~ 0 + trtGp/scTime + splines::ns(scTime,2) +
        (1|obs) + (1|trtGpRep), family=binomial(link='logit'), nAGQ=0,
                       data=HawCon, control=glmerControl(optimizer='bobyqa'))

## ----cfAICs---------------------------------------------------------
cfAIC <- 
  AIC(HCbiObs.cll,HCglmerBIobs.cll,HCbiObs.logit,HCglmerBIobs.logit,
             HCglmerBIobs2s.cll,HCglmerBIobs2s.logit)
rownames(cfAIC) <- c('TMB:cll','mer:cll','TMB:lgt','mer:lgt',
                     'mer:cllCurve', 'mer:lgtCurve')
(tab <- t(round(cfAIC,2)))

## ----allFit, results='hide', echo=TRUE------------------------------
check <- (requireNamespace('dfoptim',quietly=TRUE)&
    requireNamespace('optimx',quietly=TRUE))
if(check)
ss<-suppressWarnings(summary(allFit(HCglmerBIobs.cll)))

## ----names-ss, hold=FALSE, echo=-2----------------------------------
stopifnot(check)
options(width=70)
names(ss)
ss$msgs
ss$llik

## ----LT99gauss.LTcll------------------------------------------------
cloglog <- make.link('cloglog')$linkfun
HCgauss.cll <- glmmTMB(cloglog((PropDead+0.002)/1.004)~0+
                         trtGp/TrtTime+(TrtTime|trtGpRep), 
                       family=gaussian(), data=HawCon)
LTgauss.cll <- qra::extractLT(p=0.99, obj=HCgauss.cll, link="cloglog",
                                 a=1:8, b=9:16, eps=0.002, offset=c(0,1),                                 df.t=NULL)[,-2]
rownames(LTgauss.cll) <- shorten(rownames(LTgauss.cll))

## ----cap11, echo=FALSE----------------------------------------------
cap11 <- "Residuals versus predicted quantile deviations, for
the linear mixed model, 
with \\(log(1-log((p+0.002)/(1+0.004)))\\) as the dependent 
variable, complementary log-log link, and gaussian error."

## ----Gauss-simRes, echo=FALSE, w=4.5, fig.asp=0.75, left=-0.5, bot=-1, mgp=c(3,1,0), crop=TRUE, fig.align='center', out.width="50%", fig.cap=cap11----
sim <- simulateResiduals(HCgauss.cll)
plotResiduals(sim)

## ----cap12, echo=FALSE----------------------------------------------
cap12 <- "Comparison of estimates and
upper and lower $95\\%$ confidence limits, between the 
preferred betabinomial complementary log-log model and this
model."

## ----BBvsGauss, out.width="100%",  message=FALSE, echo=FALSE--------
library(kableExtra)
cfLTs <- cbind(LTbb.cll, LTgauss.cll)
colnames(cfLTs) <- c(rep('bb',3),rep('gauss',3))
tab <- round(cfLTs[,c(c(1,4),c(1,4)+1,c(1,4)+2)],1)
library(magrittr)
linesep <- c('', '\\addlinespace') 
knitr::kable(tab, booktabs=TRUE, linesep=linesep, format='latex', caption=cap12, 
             format.args=list(justify="right", width=9)) %>%
  kable_styling(latex_options = "striped", stripe_index = 5:8, font_size=9, full_width=FALSE) %>%
  column_spec(6:7, bold=TRUE) %>%
  add_header_above(c(' '=1,'Estimate'=2,'Lower'=2, 'Upper'=2), align='r') 

