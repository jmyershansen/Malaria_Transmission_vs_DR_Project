
install.packages("dlnm")
library(dlnm)
data(chicagoNMMAPS)
head(chicagoNMMAPS,3)
View(chicagoNMMAPS)

cb1.pm <- crossbasis(chicagoNMMAPS$pm10, lag=15, argvar=list(fun="lin"),
                       arglag=list(fun="poly",degree=4))
cb1.temp <- crossbasis(chicagoNMMAPS$temp, lag=3, argvar=list(df=5),
                         arglag=list(fun="strata",breaks=1))

summary(cb1.pm)
library(splines)
model1 <- glm(death ~ cb1.pm + cb1.temp + ns(time, 7*14) + dow,
                family=quasipoisson(), chicagoNMMAPS)
summary(model1)

pred1.pm <- crosspred(cb1.pm, model1, at=0:20, bylag=0.2, cumul=TRUE)

plot(pred1.pm, "slices", var=10, col=3, ylab="RR", ci.arg=list(density=15,lwd=2),
     main="Association with a 10-unit increase in PM10")

plot(pred1.pm, "slices", var=10, col=2, cumul=TRUE, ylab="Cumulative RR",
       main="Cumulative association with a 10-unit increase in PM10")
pred1.pm$allRRfit["10"]
cbind(pred1.pm$allRRlow, pred1.pm$allRRhigh)["10",]

#EXAMPLE2

chicagoNMMAPSseas <- subset(chicagoNMMAPS, month %in% 6:9)
cb2.o3 <- crossbasis(chicagoNMMAPSseas$o3, lag=5,
                     argvar=list(fun="thr",thr=40.3), arglag=list(fun="integer"),
                     group=chicagoNMMAPSseas$year)
cb2.temp <- crossbasis(chicagoNMMAPSseas$temp, lag=10,
                         argvar=list(fun="thr",thr=c(15,25)), arglag=list(fun="strata",breaks=c(2,6)),
                         group=chicagoNMMAPSseas$year)
model2 <- glm(death ~ cb2.o3 + cb2.temp + ns(doy, 4) + ns(time,3) + dow,
              family=quasipoisson(), chicagoNMMAPSseas)
pred2.o3 <- crosspred(cb2.o3, model2, at=c(0:65,40.3,50.3))

plot(pred2.o3, "slices", var=50.3, ci="bars", type="p", col=2, pch=19,
     ci.level=0.80, main="Lag-response a 10-unit increase above threshold (80CI)")
plot(pred2.o3,"overall",xlab="Ozone", ci="l", col=3, ylim=c(0.9,1.3), lwd=2,
       ci.arg=list(col=1,lty=3), main="Overall cumulative association for 5 lags")
pred2.o3$allRRfit["50.3"]


#EXAMPLE3
cb3.pm <- crossbasis(chicagoNMMAPS$pm10, lag=1, argvar=list(fun="lin"),
                     arglag=list(fun="strata"))
varknots <- equalknots(chicagoNMMAPS$temp,fun="bs",df=5,degree=2)
lagknots <- logknots(30, 3)
cb3.temp <- crossbasis(chicagoNMMAPS$temp, lag=30, argvar=list(fun="bs",
                                                                 knots=varknots), arglag=list(knots=lagknots))
model3 <- glm(death ~ cb3.pm + cb3.temp + ns(time, 7*14) + dow,
              family=quasipoisson(), chicagoNMMAPS)

pred3.temp <- crosspred(cb3.temp, model3, cen=21, by=1)
plot(pred3.temp, xlab="Temperature", zlab="RR", theta=200, phi=40, lphi=30,
       main="3D graph of temperature effect")
plot(pred3.temp, "contour", xlab="Temperature", key.title=title("RR"),
       plot.title=title("Contour plot",xlab="Temperature",ylab="Lag"))

plot(pred3.temp, "slices", var=-20, ci="n", col=1, ylim=c(0.95,1.25), lwd=1.5,
     main="Lag-response curves for different temperatures, ref. 21C")
for(i in 1:3) lines(pred3.temp, "slices", var=c(0,27,33)[i], col=i+1, lwd=1.5)
legend("topright",paste("Temperature =",c(-20,0,27,33)), col=1:4, lwd=1.5)
plot(pred3.temp, "slices", var=c(-20,33), lag=c(0,5), col=4,
       ci.arg=list(density=40,col=grey(0.7)))

###
###

df= data(drug)
head(drug,3)
View(drug)

df2 = df %>% pivot_longer(cols = 4:7)