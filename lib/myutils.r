library(repr)
library(forecast)
library(ggplot2)


my.figsize <- function(w, h) {
    if(missing(w)) {
        x <- getOption("repr.plot.width")
        y <- getOption("repr.plot.height")
        sprintf("width: %s, height %s", x, y)
    } else {
        options(repr.plot.width=w, repr.plot.height=h)
    }
}


# fit: (result) value of Arima in forecast pkg
# future: number of times to forecast
# past: number of past data to plot with forecast
# test: actual data to compare with forecasts
my.plot_forecast <- function(fit, future=NULL, past=NULL, test=NULL, xreg=NULL) {
    
    series <- eval(parse(text=fit$series))
    l <- length(series)
    
    if (is.null(future)) {
        future <- as.integer(0.1 * l)
    }
    
    if (is.null(past)) {
        past <- l
    }
    
    y <- forecast(fit, h=future, xreg=xreg)
    gg <- autoplot(y, include=past)
    
    if (!is.null(test)) {
        if (!is.ts(test)) {
            test <- ts(as.numeric(test), start=start(y$mean)[1], end=end(y$mean)[1])
        }
        gg <- gg +
               autolayer(test, series="Actuals") +
               autolayer(y$mean, series="Forecasts")
    }
    return(gg)
}


# updated from following source
# https://stats.stackexchange.com/questions/431545/why-isnt-the-tscv-function-allowing-for-step-size-other-than-1
my.tsCV <- function (y, forecastfunction, h = 1, window = NULL, xreg = NULL, 
          initial = 0, step = 1, silent=TRUE, count.freq=0.1, ...) 
{
    y <- as.ts(y)
    n <- length(y)
    step <- round(step)
    step_ind <- seq(step, n - 1L, by = step) ### Added line

    if (initial >= n) 
        stop("initial period too long")

    if (!is.null(xreg)) {
        xreg <- ts(as.matrix(xreg))
        if (NROW(xreg) != length(y)) 
            stop("xreg must be of the same size as y")
        tsp(xreg) <- tsp(y)
    }
    
    if (is.null(window)) {
        indx <- seq(1 + initial, n - 1L)
    } else {
        indx <- seq(window + initial, n - 1L, by = 1L)
    }
    indx <- intersect(indx, step_ind) ### Added line
    
    e.cols <- c('forecast_start', 'forecast_end', 'rmse', 'mape')
    e <- ts(matrix(NA_real_, nrow = floor(n/step), ncol = length(e.cols)))
    colnames(e) <- e.cols
    
    cnt <- 0
    print.when <- seq(0, length(indx), by=round(count.freq*length(indx)))
    for (i in indx) {
		# get new start of subset of y & xreg
    	if (is.null(window)) {
		    start <- 1L
		} else {
		    if (i - window >= 0L) {
		        start <- i - window + 1L
		    } else {
		        stop("small window")
		    }
		}
		
		y_subset <- subset(y, start=start, end = i)
		if (is.null(xreg)) {
		    fc <- try(suppressWarnings(forecastfunction(y_subset, 
		                                          h = h, ...)), silent = silent)
		    #                                      h = h)), silent = silent)
		} else {
		    xreg_subset <- as.matrix(subset(xreg, start=start, end=i))
		    fc <- try(suppressWarnings(forecastfunction(y_subset, 
		                                          h = h, xreg = xreg_subset, ...)), silent = silent)
		    #                                      h = h, xreg = xreg_subset)), silent = silent)
		}
        if (!is.element("try-error", class(fc))) {
            rmse <- sqrt(mean((fc$mean - y[i + (1:h)])^2))
            mape <- mean(abs(1 - fc$mean / y[i + (1:h)]))
            e[i/step, ] <- c(i+1, i+h, rmse, mape)
        }
        
        cnt <- cnt + 1
        if (cnt %in% print.when) {
            print(sprintf("%0.0f %% done.", 100*cnt/length(indx)))
        }
    }
    #return(na.omit(e)) # times of NA kept in e as attr(na.action)
    return(e)
}


# type 2 of time series cross validation 
# forecastfunction.2 returns avg and std of errors
my.tsCV.2 <- function (y, forecastfunction.2, h = 1, window = NULL, xreg = NULL, 
          initial = 0, step = 1, silent=TRUE, count.freq=0.1, ...) 
{
    y <- ts(y, start=1, frequency=1) # converted to n * 1 ts
    n <- length(y)
    step <- round(step)
    step_ind <- seq(step, n - 1L, by = step) ### Added line

    if (initial >= n) 
        stop("initial period too long")

    if (!is.null(xreg)) {
        xreg <- ts(xreg, start=1, frequency=1) # work even if n * 1?
        if (NROW(xreg) != length(y)) 
            stop("xreg must be of the same size as y")
        tsp(xreg) <- tsp(y)
    }
    
    if (is.null(window)) {
        indx <- seq(1 + initial, n - 1L)
    } else {
        indx <- seq(window + initial, n - 1L, by = 1L)
    }
    indx <- intersect(indx, step_ind) ### Added line
    
    e.cols <- c('forecast_start', 'forecast_end', 
                'rmse.mean', 'rmse.sigma', 'mape.mean', 'mape.sigma')
    e <- ts(matrix(NA_real_, nrow = floor(n/step), ncol = length(e.cols)))
    colnames(e) <- e.cols
    
    cnt <- 0
    print.when <- seq(0, length(indx), by=round(count.freq*length(indx)))
    for (i in indx) {
        # get new start of subset of y & xreg
        if (is.null(window)) {
            start <- 1L
        } else {
            if (i - window >= 0L) {
                start <- i - window + 1L
            } else {
                stop("small window")
            }
        }
        if (i+h > length(y)) {
            break
        } else {
            y_subset <- window(y, start, i+h)
        }
        if (is.null(xreg)) {
            fc <- try(suppressWarnings(forecastfunction.2(y_subset, 
                                                  h = h, ...)), silent = silent)
            #                                      h = h)), silent = silent) # uncomment when debugging this function
        } else {
            xreg_subset <- window(xreg, start, i+h)
            fc <- try(suppressWarnings(forecastfunction.2(y_subset, 
                                                  h = h, xreg = xreg_subset, ...)), silent = silent)
            #                                      h = h, xreg = xreg_subset)), silent = silent) # uncomment when debugging this function
        }
        if (!is.element("try-error", class(fc))) {
            e[i/step, ] <- c(i+1, i+h, fc$rmse.mean, fc$rmse.sigma, fc$mape.mean, fc$mape.sigma)
        }
        
        cnt <- cnt + 1
        if (cnt %in% print.when) {
            print(sprintf("%0.0f %% done.", 100*cnt/length(indx)))
        }
    }
    return(e)
}



# x: xts/ts obj
my.minmaxscale <- function(x, center=NULL, scale=NULL) {

	x <- na.omit(x)
	
	if (!is.null(center)) { # inv.transf if the args given
		sc <- center
		ss <- scale	
	} else if (!is.null(attr(x, 'scaled:scale'))) { # inv.transf if the attr given
		sc <- attr(x, 'scaled:center')
		ss <- attr(x, 'scaled:scale')	
	} else { # scale
	    sc <- NULL
	}

	if ('xts' %in% class(x)) {
		# set time zone same as that of x
		# as the apply function will lose the time zone
		tz <- tzone(x)
		Sys.setenv(TZ = tz)
		as.time <- function(x) {as.xts(x, tz=tz)}
	} else {
		x <- as.xts(x)
		as.time <- as.ts
	}

	# check if single time series
	if (is.null(dim(x))) {
		TO.TS <- TRUE
	} else if (dim(x)[2] == 1) {
		TO.TS <- TRUE
	} else {
		TO.TS <- FALSE
	}

	# save col names
	cols <- colnames(x)


	# transform or inv.
	if (is.null(sc)) { # transform
		sc <- apply(x, 2, mean)
		ss <- (apply(x, 2, max) - apply(x, 2, min)) / 2
		x <- apply(x, 1, function(r) {(r - sc) / ss})
	} else { # inv.transform
		# arrts reset
		x <- apply(x, 1, function(r) {r * ss + sc})
		sc <- NULL
	}

	# convert to xts with original time zone
	if (TO.TS) {
		x <- as.time(x) 
	} else {
		x <- as.time(t(x))
	}
	colnames(x) <- cols
	if (!is.null(sc)) {
		attr(x, 'scaled:center') <- sc
		attr(x, 'scaled:scale') <- ss
	}
	return(x)
}


# get price from log. return
# logret: logret of period
# period: period for logret
# prices.init: init. prices for period before first logret
my.recprice <- function(logret, period, prices.init, message=TRUE) {
    N <- nrow(logret)
    price <- NULL
    for (i in 1:period) {
        logret.sum <- tryCatch(logret[seq(i, N, by=period)],
         					   error=function(e) {
                                         if (message) {print(e)}
                                         # use NULL instead of NA to avoid error (the condition has length > 1)
                                         # if logret.sum is series
                                         return(NULL)
                                     }
                              )
        if (is.null(logret.sum)) next 
        price.rec <- rep(prices.init[i], nrow(logret.sum)) * exp(cumsum(logret.sum))
        if (is.null(price)) {
            price <- price.rec
        } else {
            price <- rbind(price, price.rec)
        }
    }
    return(price)
}


my.get_result <- function(x, group, errors=c('rmse','mape'), group.col='cs') {
    x <- x[,errors]
    x <- as.data.frame(x)
    x <- na.omit(x)
    x[group.col] <- group
    return(x)
}


my.plot_errors <- function(results, metrics=c('rmse','mape'), group.col='cs', loc='right',
                           ylog=FALSE) 
{
    func <- function(x, q, l){
        y <- (quantile(x, probs=q) + median(x))*.5
        return(data.frame(y = y, label = l)) 
    }
    if (ylog) {
        func.m <- function(x) {func(x, 0.75, signif(mean(10^x), 3))}
    } else {
        func.m <- function(x) {func(x, 0.75, signif(mean(x), 3))}
    }
    func.n <- function(x) {func(x, 0.25, paste('n', length(x), sep='='))}

    p <- NULL
    for (m in metrics) {
        p.m <- (ggplot(results, aes_string(group.col, m, fill=group.col)) 
                + geom_boxplot()
                + theme(legend.position=loc)
                + stat_summary(fun.data=func.m, fun=median, geom="text")
                + stat_summary(fun.data=func.n, fun=median, geom="text")
        )
        if (is.null(p)) {
            p <- p.m
        } else {
            p <- p + p.m
        }
    }
    if (ylog) {
            p <- (p + scale_y_continuous(trans='log10') 
                    + labs(y = paste(p$labels$y, " (log scale)")))
    }
    return(p)
}

