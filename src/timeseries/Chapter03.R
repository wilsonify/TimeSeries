# Title     : Chapter 3
# Objective : ARIMA Models
# Created by: thom
# Created on: 12/22/20



    #Example 3.2

    par(mfrow=c(2,1))
    # in the expressions below, ~ is a space and == is equal
    tsplot(sarima.sim(ar= .9, n=100), ylab="x", main=(expression(AR(1)~~~phi==+.9)))
    tsplot(sarima.sim(ar=-.9, n=100), ylab="x", main=(expression(AR(1)~~~phi==-.9)))

    #Example 3.5

    par(mfrow=c(2,1))
    tsplot(sarima.sim(ma= .9, n=100), ylab="x", main=(expression(MA(1)~~~theta==+.9)))
    tsplot(sarima.sim(ma=-.9, n=100), ylab="x", main=(expression(MA(1)~~~theta==-.9)))

    #Example 3.7

    set.seed(8675309)         # Jenny, I got your number
    x = rnorm(150, mean=5)    # Jenerate iid N(5,1)s
    arima(x, order=c(1,0,1))  # Jenstimation

    #Example 3.8

    ARMAtoMA(ar = .9,  ma = .5,  10)   # first 10 psi-weights
    ARMAtoMA(ar = -.5, ma = -.9, 10)   # first 10 pi-weights

    #Example 3.11

    z = c(1,-1.5,.75)    # coefficients of the polynomial
    (a = polyroot(z)[1]) # = 1+0.57735i,  print one root which is 1 + i 1/sqrt(3)
    arg = Arg(a)/(2*pi)  # arg in cycles/pt
    1/arg                # = 12,  the period

    set.seed(8675309)    # Jenny, it's me again
    ar2 = sarima.sim(ar=c(1.5,-.75), n = 144)
    plot(ar2, axes=FALSE, xlab="Time")
    axis(2); axis(1, at=seq(0,144,by=12)); box()  # work the plot machine
    abline(v=seq(0,144,by=12), lty=2, col=8)

    ACF = ARMAacf(ar=c(1.5,-.75), ma=0, 50)
    plot(ACF, type="h", xlab="lag")
    abline(h=0)

    #Example 3.12

    ARMAtoMA(ar=.9, ma=.5, 50)       #  for a list
    plot(ARMAtoMA(ar=.9, ma=.5, 50)) #  for a graph

    #Example 3.16

    ar2.acf = ARMAacf(ar=c(1.5,-.75), ma=0, 24)[-1]
    ar2.pacf = ARMAacf(ar=c(1.5,-.75), ma=0, 24, pacf=TRUE)
    par(mfrow=c(1,2))
    plot(ar2.acf, type="h", xlab="lag")
    abline(h=0)
    plot(ar2.pacf, type="h", xlab="lag")
    abline(h=0)

    #Example 3.18

    acf2(rec, 48)     # will produce values and a graphic
    (regr = ar.ols(rec, order=2, demean=F, intercept=TRUE))  # regression
    regr$asy.se.coef  # standard errors

    #Example 3.25

    regr = ar.ols(rec, order=2, demean=FALSE, intercept=TRUE)
    fore = predict(regr, n.ahead=24)
    ts.plot(rec, fore$pred, col=1:2, xlim=c(1980,1990), ylab="Recruitment")
    lines(fore$pred, type="p", col=2)
    lines(fore$pred+fore$se, lty="dashed", col=4)
    lines(fore$pred-fore$se, lty="dashed", col=4)

    #Example 3.26

    set.seed(666)
    x = sarima.sim(ar =.9, ma=.5, n = 100)
    xr = rev(x) # xr is the reversed data
    pxr = predict(arima(xr, order=c(1,0,1)), 10) # predict the reversed data
    pxrp = rev(pxr$pred) # reorder the predictors (for plotting)
    pxrse = rev(pxr$se) # reorder the SEs
    nx = ts(c(pxrp, x), start=-9) # attach the backcasts to the data
    tsplot(nx, ylab=expression(X[~t]), main='Backcasting', ylim=c(-7,4))
    U = nx[1:10] + pxrse; L = nx[1:10] - pxrse
    xx = c(-9:0, 0:-9); yy = c(L, rev(U))
    polygon(xx, yy, border = 8, col = gray(0.6, alpha = 0.2))
    lines(-9:0, nx[1:10], col=2, type='o')

    #Example 3.28

    rec.yw = ar.yw(rec, order=2)
    rec.yw$x.mean  # = 62.26278 (mean estimate)
    rec.yw$ar      # = 1.3315874, -.4445447  (parameter estimates)
    sqrt(diag(rec.yw$asy.var.coef))  # = .04222637, .04222637  (standard errors)
    rec.yw$var.pred  # = 94.79912 (error variance estimate)

    rec.pr = predict(rec.yw, n.ahead=24)
    U = rec.pr$pred + rec.pr$se
    L = rec.pr$pred - rec.pr$se
    minx = min(rec,L); maxx = max(rec,U)
    ts.plot(rec, rec.pr$pred, xlim=c(1980,1990), ylim=c(minx,maxx))
    lines(rec.pr$pred, col="red", type="o")
    lines(U, col="blue", lty="dashed")
    lines(L, col="blue", lty="dashed")

    #Example 3.29

    set.seed(1)
    ma1 = sarima.sim(ma=.9, n=50)
    acf1(ma1, 1, plot=FALSE)  # [1] .536 (lag 1 sample ACF)

    #Example 3.31

    # Note: I'm not convinced this is really the MLE...
    #  ... but eventually 'sarima()' will be used
    rec.mle = ar.mle(rec, order=2)
    rec.mle$x.mean
    rec.mle$ar
    sqrt(diag(rec.mle$asy.var.coef))
    rec.mle$var.pred

    #Example 3.33

    x = diff(log(varve))
    # Evaluate Sc on a Grid
    c(0) -> w -> z -> Sc -> Sz -> Szw
    num = length(x)
    th  = seq(-.3,-.94,-.01)
    for (p in 1:length(th)){
    for (i in 2:num){ w[i] = x[i]-th[p]*w[i-1] }
    Sc[p] = sum(w^2) }
    tsplot(th, Sc, type="l", ylab=expression(S[c](theta)), xlab=expression(theta), lwd=2)
    # Gauss-Newton Estimation
    r = acf(x, lag=1, plot=FALSE)$acf[-1]
    rstart = (1-sqrt(1-4*(r^2)))/(2*r) # from (3.105)
    c(0) -> w -> z
    c() -> Sc -> Sz -> Szw -> para
    niter = 12
    para[1] = rstart
    for (p in 1:niter){
      for (i in 2:num){ w[i] = x[i]-para[p]*w[i-1]
      z[i] = w[i-1]-para[p]*z[i-1]
      }
     Sc[p]     = sum(w^2)
     Sz[p]     = sum(z^2)
     Szw[p]    = sum(z*w)
     para[p+1] = para[p] + Szw[p]/Sz[p]
    }
    round(cbind(iteration=0:(niter-1), thetahat=para[1:niter] , Sc , Sz ), 3)
    abline(v = para[1:12], lty=2)
    points(para[1:12], Sc[1:12], pch=16)

    #Example 3.36

    # generate data
    set.seed(101010)
    e   = rexp(150, rate=.5)
    u   = runif(150,-1,1)
    de  = e*sign(u)
    dex = 50 + sarima.sim(n=100, ar=.95, innov=de, burnin=0)
    #-- as in text - same result
    # dex = 50 + arima.sim(n=100, list(ar=.95), innov=de, n.start=50)
    tsplot(dex, ylab=expression(X[~t]))

    # small sample and asymptotic distn
    set.seed(111)
    phi.yw = rep(NA, 1000)
    for (i in 1:1000){
     e = rexp(150, rate=.5); u = runif(150,-1,1); de = e*sign(u)
     x = 50 + arima.sim(n=100,list(ar=.95), innov=de, n.start=50)
     phi.yw[i] = ar.yw(x, order=1)$ar
    }
    hist(phi.yw, prob=TRUE, main="", ylim=c(0,14), xlim=c(.70,1.05))
    lines(density(phi.yw, bw=.015))
    u = seq(.75, 1.1, by=.001)
    lines(u, dnorm(u, mean=.96, sd=.03), lty=2, lwd=2)

    # Bootstrap
    set.seed(666)                 # not that 666
    fit     = ar.yw(dex, order=1) # assumes the data were retained
    m       = fit$x.mean          # estimate of mean
    phi     = fit$ar              # estimate of phi
    nboot   = 250                 # number of bootstrap replicates
    resids  = fit$resid[-1]       # the 99 innovations
    x.star  = dex                 # initialize x*
    phi.star.yw = rep(NA, nboot)
    #- start it up
    for (i in 1:nboot) {
      resid.star = sample(resids, replace=TRUE)
      for (t in 1:99){ x.star[t+1] = m + phi*(x.star[t]-m) + resid.star[t]
      }
     phi.star.yw[i] = ar.yw(x.star, order=1)$ar
    }
    # Picture
    culer = rgb(.5,.7,1,.5)
    hist(phi.star.yw, 15, main="", prob=TRUE, xlim=c(.65,1.05), ylim=c(0,14),
    col=culer, xlab=expression(hat(phi)))
    lines(density(phi.yw, bw=.02), lwd=2) # from previous simulation
    u = seq(.75, 1.1, by=.001) # normal approximation
    lines(u, dnorm(u, mean=.96, sd=.03), lty=2, lwd=2)
    legend(.65, 14, legend=c('true distribution', 'bootstrap distribution',
    'normal approximation'), bty='n', lty=c(1,0,2), lwd=c(2,0,2),
    col=1, pch=c(NA,22,NA), pt.bg=c(NA,culer,NA), pt.cex=2.5)

    #Example 3.38

    set.seed(666)
    x = sarima.sim(ma = -0.8, d=1, n = 100)
    (x.ima = HoltWinters(x, beta=FALSE, gamma=FALSE))  # α is 1-λ here
    plot(x.ima)

    #Example 3.39, 3.40, and 3.43

    tsplot(gnp)
    acf2(gnp, 50)
    gnpgr = diff(log(gnp))      # growth rate
    tsplot(gnpgr)
    acf2(gnpgr, 24)
    sarima(gnpgr, 1, 0, 0)      # AR(1)
    sarima(gnpgr, 0, 0, 2)      # MA(2)
    ARMAtoMA(ar=.35, ma=0, 10)  # prints psi-weights

    #Example 3.41

    sarima(log(varve), 0, 1, 1, no.constant=TRUE)   # ARIMA(0,1,1)
    dev.new()
    sarima(log(varve), 1, 1, 1, no.constant=TRUE)   # ARIMA(1,1,1)

    #Example 3.44

    trend  = time(cmort)
    temp   = tempr - mean(tempr)
    temp2  = temp^2
    summary(fit <- lm(cmort~trend + temp + temp2 + part, na.action=NULL))
    acf2(resid(fit), 52) # implies AR2
    sarima(cmort, 2,0,0, xreg=cbind(trend,temp,temp2,part) )

    #Example 3.45

    # Note: this could benefit from a seasonal model fit, but it hasn't
    #  been talked about yet - you could come back to this after the next section
    dummy = ifelse(soi<0, 0, 1)
    fish = ts.intersect(rec, soiL6=lag(soi,-6), dL6=lag(dummy,-6), dframe=TRUE)
    summary(fit <- lm(rec ~soiL6*dL6, data=fish, na.action=NULL))
    attach(fish)
    plot(resid(fit))
    acf2(resid(fit))     # indicates AR(2)
    intract = soiL6*dL6  # interaction term
    sarima(rec,2,0,0, xreg = cbind(soiL6, dL6, intract))
    # not in text, but this works better
    # sarima(rec,2,0,0,0,1,1,12, xreg = cbind(soiL6, dL6, intract))


    #Example 3.46

    set.seed(666)
    phi  = c(rep(0,11),.9)
    sAR  = arima.sim(list(order=c(12,0,0), ar=phi), n=37)
    sAR  = ts(sAR, freq=12)
    layout(matrix(c(1,1,2, 1,1,3), nc=2))
    par(mar=c(3,3,2,1), mgp=c(1.6,.6,0))
    plot(sAR, axes=FALSE, main='seasonal AR(1)', xlab="year", type='c')
    Months = c("J","F","M","A","M","J","J","A","S","O","N","D")
    points(sAR, pch=Months, cex=1.25, font=4, col=1:4)
    axis(1, 1:4)
    abline(v=1:4, lty=2, col=gray(.7))
    axis(2)
    box()
    ACF  = ARMAacf(ar=phi, ma=0, 100)
    PACF = ARMAacf(ar=phi, ma=0, 100, pacf=TRUE)
    plot(ACF,type="h", xlab="LAG", ylim=c(-.1,1))
    abline(h=0)
    plot(PACF, type="h", xlab="LAG", ylim=c(-.1,1))
    abline(h=0)

    #Example 3.47

    phi  = c(rep(0,11),.8)
    ACF  = ARMAacf(ar=phi, ma=-.5, 50)[-1] # [-1] removes 0 lag
    PACF = ARMAacf(ar=phi, ma=-.5, 50, pacf=TRUE)
    par(mfrow=c(1,2))
    plot(ACF, type="h", xlab="LAG", ylim=c(-.4,.8)); abline(h=0)
    plot(PACF, type="h", xlab="LAG", ylim=c(-.4,.8)); abline(h=0)

    #Example 3.49

    x     = AirPassengers
    lx    = log(x)
    dlx   = diff(lx)
    ddlx  = diff(dlx, 12)
    plot.ts(cbind(x, lx, dlx, ddlx), main="")
    # below of interest for showing seasonal RW (not shown here):
    par(mfrow=c(2,1))
    monthplot(dlx)
    monthplot(ddlx)

    sarima(lx, 1,1,1, 0,1,1, 12)   # model 1
    sarima(lx, 0,1,1, 0,1,1, 12)   # model 2 (the winner)
    sarima(lx, 1,1,0, 0,1,1, 12)   # model 3

    sarima.for(lx, 12, 0,1,1, 0,1,1,12)  # forecasts