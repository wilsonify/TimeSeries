# Title     : Chapter 5
# Objective : Additional Time Domain Topics
# Created by: thom
# Created on: 12/22/20


    #Example 5.1

    # NOTE: I think 'fracdiff' is a dinosaur and I should have changed
    # this in the new edition... so just below, I'll do it using 'arfima',
    # which seems to work well
    library(fracdiff)
    lvarve = log(varve) - mean(log(varve))
    varve.fd = fracdiff(lvarve, nar=0, nma=0, M=30)
    varve.fd$d  # = 0.3841688
    varve.fd$stderror.dpq  # = 4.589514e-06 (If you believe this, I have a bridge for sale.)

    p = rep(1,31)
    for (k in 1:30){ p[k+1] = (k-varve.fd$d)*p[k]/(k+1) }
    plot(1:30, p[-1], ylab=expression(pi(d)), xlab="Index", type="h", lwd=2)
    res.fd = diffseries(log(varve), varve.fd$d)       # frac diff resids
    res.arima = resid(arima(log(varve), order=c(1,1,1))) # arima resids

    dev.new()
    par(mfrow=c(2,1))
    acf(res.arima, 100, xlim=c(4,97), ylim=c(-.2,.2), main="arima resids")
    acf(res.fd, 100, xlim=c(4,97), ylim=c(-.2,.2), main="frac diff resids")

    #Example 5.1 redux

    library(arfima)
    summary(varve.fd <- arfima(log(varve)))  # d.hat = 0.3728, se(d,hat) = 0.0273
    # residual stuff
    innov = resid(varve.fd)
    plot.ts(innov[[1]])
    acf(innov[[1]])

    ## ... much better ...  sorry I didn't ...
    ## ... get it in for the newest edition ..
    ## ... once in awhile, they slip on by ...

    #Example 5.2

    series = log(varve)  # specify series to be analyzed
    d0 = .1              # initial value of d
    n.per = nextn(length(series))
    m = (n.per)/2  - 1
    per = abs(fft(series-mean(series))[-1])^2  # remove 0 freq
    per = per/n.per      # R doesn't scale fft by sqrt(n)
    g = 4*(sin(pi*((1:m)/n.per))^2)

    # Function to calculate -log.likelihood
    whit.like = function(d){
     g.d=g^d
     sig2 = (sum(g.d*per[1:m])/m)
     log.like = m*log(sig2) - d*sum(log(g)) + m
     return(log.like)
    }

    # Estimation (?optim for details - output not shown)
    (est = optim(d0, whit.like, gr=NULL, method="L-BFGS-B", hessian=TRUE, lower=-.5, upper=.5,
                 control=list(trace=1,REPORT=1)))

    # Results  [d.hat = .380, se(dhat) = .028]
    cat("d.hat =", est$par, "se(dhat) = ",1/sqrt(est$hessian),"\n")
    g.dhat = g^est$par
    sig2 = sum(g.dhat*per[1:m])/m
    cat("sig2hat =",sig2,"\n")  # sig2hat = .229

    u = spec.ar(log(varve), plot=FALSE)  # produces AR(8)
    g = 4*(sin(pi*((1:500)/2000))^2)
    fhat = sig2*g^{-est$par}             # long memory spectral estimate
    plot(1:500/2000, log(fhat), type="l", ylab="log(spectrum)", xlab="frequency")
    lines(u$freq[1:250], log(u$spec[1:250]), lty="dashed")
    ar.mle(log(varve))                   # to get AR(8) estimates

    # GPH estimate with big bandwidth or else estimate sucks
    fdGPH(log(varve), bandw=.9)   # m = n^bandw

    #Example 5.3

    library(tseries)
    adf.test(log(varve), k=0)  # DF test
    adf.test(log(varve))       # ADF test
    pp.test(log(varve))        # PP test

    #Example 5.4

    gnpgr = diff(log(gnp))          # get the returns
    u     = sarima(gnpgr, 1, 0, 0)  # fit an AR(1)
    acf2(resid(u$fit), 20)          # get (p)acf of the squared residuals

    library(fGarch)
    summary(garchFit(~arma(1,0)+garch(1,0), gnpgr))

    #Example 5.5 and 5.6

    library(xts)   # needed to handle djia
    djiar = diff(log(djia$Close))[-1]
    acf2(djiar)    # exhibits some autocorrelation (not shown)
    acf2(djiar^2)  # oozes autocorrelation (not shown)
    library(fGarch)
    # GARCH fit
    summary(djia.g <- garchFit(~arma(1,0)+garch(1,1), data=djiar, cond.dist='std'))
    plot(djia.g)    # to see all plot options
    # APARCH fit
    summary(djia.ap <- garchFit(~arma(1,0)+aparch(1,1), data=djiar, cond.dist='std'))
    plot(djia.ap)

    #Example 5.7

    tsplot(flu, type="c")
    Months = c("J","F","M","A","M","J","J","A","S","O","N","D")
    points(flu, pch=Months, cex=.8, font=2)
    # Start analysis
    dflu = diff(flu)
    lag1.plot(dflu, corr=FALSE) # scatterplot with lowess fit
    thrsh = .05 # threshold
    Z = ts.intersect(dflu, lag(dflu,-1), lag(dflu,-2), lag(dflu,-3),
    lag(dflu,-4) )
    ind1 = ifelse(Z[,2] < thrsh, 1, NA) # indicator < thrsh
    ind2 = ifelse(Z[,2] < thrsh, NA, 1) # indicator >= thrsh
    X1 = Z[,1]*ind1
    X2 = Z[,1]*ind2
    summary(fit1 <- lm(X1~ Z[,2:5]) ) # case 1
    summary(fit2 <- lm(X2~ Z[,2:5]) ) # case 2
    D = cbind(rep(1, nrow(Z)), Z[,2:5]) # design matrix
    p1 = D %*% coef(fit1) # get predictions
    p2 = D %*% coef(fit2)
    prd = ifelse(Z[,2] < thrsh, p1, p2)
    plot(dflu, ylim=c(-.5,.5), type='p', pch=3)
    lines(prd)
    prde1 = sqrt(sum(resid(fit1)^2)/df.residual(fit1) )
    prde2 = sqrt(sum(resid(fit2)^2)/df.residual(fit2) )
    prde = ifelse(Z[,2] < thrsh, prde1, prde2)
    tx = time(dflu)[-(1:4)]
    xx = c(tx, rev(tx))
    yy = c(prd-2*prde, rev(prd+2*prde))
    polygon(xx, yy, border=8, col=gray(.6, alpha=.25) )
    abline(h=.05, col=4, lty=6)

    # Using tsDyn (not in text)
    library(tsDyn)
    # vignette("tsDyn")   # for package details (it's quirky, so you'll need this)
    dflu = diff(flu)
    (u = setar(dflu, m=4, thDelay=0))  # fit model and view results (thDelay=0 is lag 1 delay)
    BIC(u); AIC(u)                     # if you want to try other models ... m=3 works well too
    plot(u)                            # graphics -  ?plot.setar for information

    #Example 5.8 and 5.9

    soi.d   = resid(lm(soi~time(soi), na.action=NULL)) # detrended SOI
    acf2(soi.d)
    fit     = arima(soi.d, order=c(1,0,0))
    ar1     = as.numeric(coef(fit)[1]) # = 0.5875
    soi.pw  = resid(fit)
    rec.fil = filter(rec, filter=c(1, -ar1), sides=1)
    ccf2(soi.pw, rec.fil, na.action=na.omit)

    fish  = ts.intersect(rec, RL1=lag(rec,-1), SL5=lag(soi.d,-5))
    (u    = lm(fish[,1]~fish[,2:3], na.action=NULL))
    acf2(resid(u)) # suggests ar1
    (arx  = sarima(fish[,1], 1, 0, 0, xreg=fish[,2:3])) # final model
    pred  = rec + resid(arx$fit) # 1-step-ahead predictions
    ts.plot(pred, rec, col=c('gray90',1), lwd=c(7,1))

    #Example 5.10 and 5.11

    library(vars)
    x = cbind(cmort, tempr, part)
    summary(VAR(x, p=1, type="both"))  # "both" fits constant + trend

    VARselect(x, lag.max=10, type="both")
    summary(fit <- VAR(x, p=2, type="both"))
    acf(resid(fit), 52)
    serial.test(fit, lags.pt=12, type="PT.adjusted")

    (fit.pr = predict(fit, n.ahead = 24, ci = 0.95))  # 4 weeks ahead
    dev.new()
    fanchart(fit.pr)  # plot prediction + error

    #Example 5.12

    library(marima)
    model   = define.model(kvar=3, ar=c(1,2), ma=c(1))
    arp     = model$ar.pattern
    map     = model$ma.pattern
    cmort.d = resid(detr <- lm(cmort~ time(cmort), na.action=NULL))
    xdata   = matrix(cbind(cmort.d, tempr, part), ncol=3)  # strip ts attributes
    fit     = marima(xdata, ar.pattern=arp, ma.pattern=map, means=c(0,1,1), penalty=1)
    # resid analysis (not displayed)
    innov   = t(resid(fit))
    plot.ts(innov)
    acf(innov)
    # fitted values for cmort
    pred    = ts(t(fitted(fit))[,1], start=start(cmort), freq=frequency(cmort)) +
    detr$coef[1] + detr$coef[2]*time(cmort)
    plot(pred, ylab="Cardiovascular Mortality", lwd=2, col=4)
    points(cmort)
    # print estimates and corresponding t^2-statistic
    short.form(fit$ar.estimates, leading=FALSE)
    short.form(fit$ar.fvalues,   leading=FALSE)
    short.form(fit$ma.estimates, leading=FALSE)
    short.form(fit$ma.fvalues,   leading=FALSE)
    fit$resid.cov # estimate of noise cov matrix
