# ARIMA+GARCH Strategy - fit & predict
ag.forecast <- function(ret, foreLength, windowLength, lookahead) {

    # Create vectors to store the predictions
    fc.pred <- vector(mode="numeric", length=foreLength)
    fc.sigma <- vector(mode="numeric", length=foreLength)
    fc.dates <- vector(mode="character", length=foreLength)
    
    # step to print output
    print.step <- round((foreLength+1)/5)

    start_time <- Sys.time()

    for (d in 1:(foreLength+1)) {

        # Obtain the S&P500 rolling window for this day
        retOffset <- ret[d:(windowLength+d)]

        # Fit the ARIMA model
        fit <- tryCatch(auto.arima(retOffset, seasonal=FALSE, 
                                   ic='aicc', 
                                   #ic='aic', 
                                   d=0, 
                                   trace=FALSE),
                        error=function(err) {FALSE},
                        warning=function(err) {FALSE} )

        if( !is.logical(fit) ) {
            final.order <- as.numeric(arimaorder(fit))
        } else {
            final.order <- c(0,0,0)
        }

        # Specify and fit the GARCH model
        spec <- ugarchspec(
            variance.model=list(garchOrder=c(1,1)),
            mean.model=list(armaOrder=c(final.order[1], final.order[3]), include.mean=T),
            distribution.model="sged"
        )
        fit <- tryCatch(ugarchfit(spec, retOffset, solver='hybrid'), 
                        error=function(e) {e}, 
                        warning=function(w) {w})

        # If the GARCH model does not converge, set the direction to "long" else
        # choose the correct forecast direction based on the returns prediction
        # Output the results to the screen and the forecasts vector
        if(is(fit, "warning")) {
            fc.pred[d] <- 0 # assume zero return
            fc.sigma[d] <- 0
            fc.dates[d] <- as.character(index(tail(retOffset, 1)[1]))
        } else {
            fore <- ugarchforecast(fit, n.ahead=1)
            ind <- fore@forecast$seriesFor
            sig <- fore@forecast$sigmaFor
            # save d + lookahead prediction at d date
            fc.pred[d] <- ind[1]
            fc.sigma[d] <- sig[1]
            fc.dates[d] <- colnames(ind)
        }

        if ((d %% print.step == 0) | (d == foreLength+1)) {
            a <- Sys.time() - start_time
            u <- attr(a, "units")
            v <- as.numeric(a)
            print(sprintf("%0.0f %% done (%0.0f %s).", 100 * d/(foreLength+1), v, u))
        }
    }

    pred <- xts(cbind(fc.pred, fc.sigma), order.by=as.Date(fc.dates, "%Y-%m-%d"))
    colnames(pred) <- c(paste('series T+', lookahead, sep=''), paste('sigma T+', lookahead, sep=''))
    return(pred)
}


# ARIMA+GARCH Strategy - plot the result of ag.forecast
# pred: return value of ag.forecast
ag.plotforecast <- function(pred, trueCurve, lookahead, 
                            label='S&P 500', n_sigma=2) 
{
    trueCurve <- trueCurve[index(pred)]
    predSer <- pred[,1]
    predSig <- pred[,2]

    getPrice <- function(x) exp(x) * trueCurve # calc returns
    z <- getPrice(predSer)
    zup <- getPrice(predSer + n_sigma*predSig)
    zdn <- getPrice(predSer - n_sigma*predSig)

    # merge forecast & true values
    result <- merge(z, zup)
    result <- merge(result, zdn)
    result <- lag(result, lookahead) # shift prediction by lookahead
    result <- merge(result, trueCurve, all=T)
    result <- result[!is.na(result[,1])]
    colnames(result) <- c('Z', 'Zup', 'Zdn', 'price')

    # plot
    labelPr <- label
    labelFc <- paste("Forecast w/ conditional ", n_sigma, "-sigma bands", sep="")
    num_lines <- 2
    gg <- (ggplot(result, aes(Index, Z))
             + geom_line(aes(color=labelFc), size = 1)
             + geom_ribbon(aes(ymin = Zdn, ymax = Zup), fill=hue_pal()(num_lines)[1], alpha=0.1)
             + geom_line(aes(y=price, color=labelPr), size = 1)
             + ylab('Index') + xlab("")
             + theme(legend.position="bottom", legend.title=element_blank())
             + guides(color=guide_legend(nrow=num_lines))
    )
    return(gg)
}

# ARIMA+GARCH Strategy - trading strategy
# pred: return value of ag.forecast
ag.backtest <- function(pred, ret, lookahead) {
    if (dim(pred)[2]>1) stop("dim of pred is invalid.")
    
    # shift prediction to its realized day before merge with ret
    position <- lag(pred, lookahead)
    position <- pred[!is.na(position)]
    # define position
    position <- ifelse(position < 0, -1, ifelse(position > 0, 1, 0))
    # merge position with ret
    position <- merge(position, ret, all=F )
    # calc result of position
    result <- position[,1] * position[,2]

    # Create the backtests for ARIMA+GARCH and Buy & Hold
    # assuming reinvestment of returns at every day
    # remember return is logarithmic one
    resArimaGarchCurve = cumsum(result)
    resBuyHoldCurve = cumsum(position[,2])
    
    resCombinedCurve = merge(resArimaGarchCurve, resBuyHoldCurve, all=F )

    # Plot the equity curves
    xyp <- ag.plot(resCombinedCurve, labels=c("ARIMA+GARCH", "Buy & Hold"))
    return(xyp)
}

# see ag.backtest
ag.plot <- function(data, 
                    labels=c("ARIMA+GARCH", "Buy & Hold"),
                    colors=c("darkred", "darkblue")) {
    xyp <- xyplot( 
        data,
        superpose=T,
        col=colors,
        lwd=2,
        key=list(
            text=list(labels),
            lines=list(lwd=2, col=colors)
        )
    )
    return(xyp)
}



# ARIMA+GARCH Forecasts - fit & predict
ag2.forecast <- function(ret, foreLength, out.sample=0)
{
    # get train data for ARIMA model by dropping out-of-sample
    train <- ret[1:(nrow(as.xts(ret))-foreLength)]
    # Fit the ARIMA model
    fit <- tryCatch(auto.arima(train, seasonal=FALSE, 
                               ic='aicc', 
                               #ic='aic', 
                               d=0, 
                               trace=FALSE),
                    error=function(err) {FALSE},
                    warning=function(err) {FALSE} )

    if( !is.logical(fit) ) {
        final.order <- as.numeric(arimaorder(fit))
    } else {
        final.order <- c(0,0,0)
    }

    # Specify and fit the GARCH model
    spec <- ugarchspec(
        variance.model=list(garchOrder=c(1,1)),
        mean.model=list(armaOrder=c(final.order[1], final.order[3]), include.mean=T),
        distribution.model="sged"
    )
    fit <- tryCatch(ugarchfit(spec, ret, solver='hybrid', out.sample=out.sample), 
                    error=function(e) {e}, 
                    warning=function(w) {w})

    # If the GARCH model does not converge, set the direction to "long" else
    # choose the correct forecast direction based on the returns prediction
    # Output the results to the screen and the forecasts vector
    if(is(fit, "warning")) {
        print('GARCH model does not converge')
        fore <- NA
    } else {
        fore <- ugarchforecast(fit, n.ahead=foreLength, n.roll=out.sample)
    }
    #plot(fore, which=2)
    return(fore)
}

# ARIMA+GARCH Forecasts - plot
# fore: return value of ag2.forecast
ag2.plot <- function(fore, var.mode='unconditional', plot.mode='return',
                     price=NULL, lookahead=NULL,
                     figsize=c(10,6),
                     plotFUN = paste("ag2.plot.garchforecast", 1:2, sep = ".")
                    ) 
{
    my.figsize(figsize[1],figsize[2])
    vmode <- ifelse(var.mode=='unconditional', 1, 2)
    if (vmode == 2) {
        if (forc@forecast$n.roll < 5) {
            return('ERROR: Run ag2.forecast with out.sample >= 5')
        }
    }
    
    if (plot.mode=='return') {
        plot(fore, which=vmode)
    } else {
        if ((is.null(price) | (is.null(lookahead)))) {
            return('ERROR: price data & lookahead required.')
        } else {
            func <- match.fun(plotFUN[vmode])
            func(fore, price, lookahead)
        }
    }
}


# ARIMA+GARCH Forecasts - plotFUN of ag2.plot in the case of unconditional variance
ag2.plot.garchforecast.1 <- function(x, price, lookahead,
                                    size.title=1.0, color.x.sig='grey', 
                                    ylab='Index',
                                    #history.length=42,
                                    history.length=25,
                                    ...)
{
    vmodel = x@model$modeldesc$vmodel
    # 1. Time Series:
    #nr = x@forecast$n.roll
    #if(n.roll > nr) stop("plot-->error: n.roll choice is invalid", call. = FALSE)
    #n.roll <- nr
    n.roll <- 0
    n = x@forecast$n.ahead
    N = x@forecast$N - x@forecast$n.start
    forseries = x@forecast$seriesFor[,n.roll+1]
    forsigma = x@forecast$sigmaFor[,n.roll+1]
    xdates = x@model$modeldata$index[(N+n.roll-min(N,history.length)):(N+n.roll)]
    fdates = seq(tail(xdates,1), by = x@model$modeldata$period, length.out = n+1)[-1]
    series = x@model$modeldata$data[(N+n.roll-min(N,history.length)):(N+n.roll)]

    xforseries = c(series, forseries)
    series = c(series, rep(NA, n))
    Zup = forseries+1*forsigma
    Zdn = forseries-1*forsigma
    
    ### added
    getPrice <- function(x, d) {
        logret <- as.xts(x, order.by=d)
        price.init <- tail(price[index(price) < d[1]], lookahead)
        price.rec <- my.recprice(logret, lookahead, price.init, message=FALSE)
        return(as.numeric(price.rec))
    }
    xforseries <- getPrice(xforseries, c(xdates, fdates))
    series <- getPrice(series, c(xdates, fdates))
    forseries <- getPrice(forseries, fdates)
    Zup <- getPrice(Zup, fdates)
    Zdn <- getPrice(Zdn, fdates)
    
    #ylim=c(0.95*min(xforseries,na.rm=TRUE), 1.05*max(xforseries,na.rm=TRUE))
    ylim=c(0.99*min(xforseries, Zdn, na.rm=TRUE), 1.01*max(xforseries, Zup, na.rm=TRUE))
    plot(c(xdates, fdates), as.numeric(xforseries), type="l", col="steelblue",
            main = paste("Forecast Series\n w/th unconditional 1-Sigma bands", sep = "") ,
            ylab="Series",xlab="Time/Horizon", 
            ylim = ylim, 
            cex.main = 0.7, cex.axis = 0.8, cex.lab = 0.9)
    abline(h = 0, col = "grey", lty = 3)
    for(i in 2:n) rect(fdates[i-1], Zdn[i-1], fdates[i], Zup[i], col = colors()[142], border=NA)
    lines(c(xdates, fdates), series, col="steelblue")
    lines(fdates, forseries, col="tomato1")
    mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
    if(vmodel=="fGARCH"){
        mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
        mtext(paste("Horizon: ",n,sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
    } else{
        mtext(paste("Horizon: ",n,sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
    }
    lg.txt = c("Actual","Forecast")
    legend("topleft", legend = lg.txt, col = c("steelblue", "tomato1"), y.intersp = 1.5, pch = 21, cex = 0.7, bty="n")
    box()
    grid()
}


# ARIMA+GARCH Forecasts - plotFUN of ag2.plot in the case of conditional variance
ag2.plot.garchforecast.2 <- function(x, price, lookahead,
                                    size.title=1.0, color.x.sig='grey', 
                                    ylab='Index',
                                    #history.length=42,
                                    history.length=25,
                                    ...)
{
    vmodel = x@model$modeldesc$vmodel
    # 1. Time Series:
    nr = x@forecast$n.roll
    if(nr<5) stop("\nn.roll less than 5!...does not make sense to provide this plot.")
    N = x@forecast$N - x@forecast$n.start
    fdata  = x@forecast$seriesFor[1,]
    fsigma = x@forecast$sigmaFor[1,]
    xdata = x@model$modeldata$data[((N+1)-min(N, history.length)):(N+nr)]
    xsigma = x@model$modeldata$sigma[((N+1)-min(N, history.length)):(N+nr)]
    xdates = x@model$modeldata$index[((N+1)-min(N, history.length)):(N+nr)]
    ns = length(xdata)
    xplus =  xdata + 2*xsigma
    xminus = xdata - 2*xsigma
    fplus = c(rep(NA, (ns - nr)), fdata[-length(fdata)] +  2*fsigma[-length(fsigma)])
    fminus = c(rep(NA, (ns - nr)), fdata[-length(fdata)] - 2*fsigma[-length(fsigma)])
    fdata = c(rep(NA, (ns - nr)), fdata[-length(fdata)])
    
    ### added
    price.shitfed <- lag(price, lookahead)[xdates]
    getPrice <- function(x) exp(x) * price.shitfed # calc returns
    xdata <- getPrice(xdata)
    xplus <- getPrice(xplus)
    xminus <- getPrice(xminus)
    fdata <- getPrice(fdata)
    fplus <- getPrice(fplus)
    fminus <- getPrice(fminus)

    #ylim=c(0.95*min(xminus,na.rm=TRUE), 1.2*max(xplus,na.rm=TRUE))
    ylim=c(0.99*min(xminus, fminus, na.rm=TRUE), 1.01*max(xplus, fplus, na.rm=TRUE))
    plot(xdates, xdata, type="l", col="black",
            main = paste("Rolling Forecast vs Actual Series\n w/th conditional 2-Sigma bands", sep = "") ,
            ylab=ylab,xlab="Time/Horizon", 
            ylim = ylim, 
            cex.main = size.title,
            cex.axis = 0.8, cex.lab = 0.9)
    for(i in 26:(nr+25)) rect(xdates[i-1], fminus[i-1], xdates[i], fplus[i], col = colors()[142], border=NA)
    lines(xdates, xdata, col = "steelblue")
    lines(xdates, xplus, col = color.x.sig, lwd = 0.5)
    lines(xdates, xminus, col = color.x.sig, lwd = 0.5)
    abline(h = 0, col = "grey", lty = 3)
    lines(xdates, fdata, col = "tomato1", lwd = 2.5)
    #lines(xts(fplus,xdates),  col = "brown", lwd = 0.5)
    #lines(xts(fminus,xdates),  col = "brown", lwd = 0.5)

    mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
    if(vmodel=="fGARCH"){
        mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
        mtext(paste("Horizon: ", nr, sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
    } else{
        mtext(paste("Horizon: ", nr, sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
    }
    lg.txt = c("Actual","Forecast")
    legend("topleft", legend = lg.txt, col = c("steelblue", "tomato1"), y.intersp = 1.5, pch = 21, cex = 0.7, bty="n")
    box()
    grid()
}



