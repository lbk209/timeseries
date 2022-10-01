library(xgboost)


xgb.forecast <- function(
    train.label, train.features, test.label, test.features,
    nrounds = 1000, early_stopping_rounds = 3,
    max_depth = 6, 
    eta = 0.3, # learning rate
    # In linear regression task, this simply corresponds to minimum number of instances needed to be in each node. 
    # The larger min_child_weight is, the more conservative the algorithm will be. range: [0,?]
    #min_child_weight = 1 ,
    # Minimum loss reduction required to make a further partition on a leaf node of the tree. 
    # The larger gamma is, the more conservative the algorithm will be. range: [0,?]
    #gamma = 0,
    verbose=0,
    sample.n=0,
    result.error=FALSE
) 
{
    if (is.null(dim(test.features))) {
        # conver to matrix for the case of 1 day prediction
        test.features <- t(test.features)
    }
    model <- xgboost(data = train.features,
                     label = train.label,
                     nrounds = nrounds,
                     objective = "reg:squarederror",
                     early_stopping_rounds = early_stopping_rounds,
                     max_depth = max_depth,
                     eta = eta,
                     verbose=verbose)

    pred <- predict(model, test.features)
    
    if (result.error) {
        h = length(pred)
        idx = 1:h
        if ((sample.n>0) & (sample.n<h)) {
            idx <- sort(sample(idx, sample.n))
        }

        # calc errors
        rmse <- sqrt((pred[idx] - test.label[idx])^2)
        mape <- abs(1 - pred[idx] / test.label[idx])
        result <- list(rmse.mean=mean(rmse), rmse.sigma=sd(rmse), 
                       mape.mean=mean(mape), mape.sigma=sd(mape))
        return(result)
    } else {
        res <- list(model=model, pred=pred)
        return(res)
    }
}


xgb.gridsearch <- function(label, features, test.h, hyper_grid, verbose=1) {
    len <- length(label)
    idx.train <- 1:(len-test.h)
    idx.test <- (len-test.h+1):(len)
    train.label <- label[idx.train]
    train.features <- features[idx.train,]
    test.label <- label[idx.test]
    test.features <- features[idx.test,]
    
    xgb_test_rmse <- NULL
    xgb_test_mape <- NULL

    for (j in 1:nrow(hyper_grid)) {
        #set.seed(123)
        
        errors <- xgb.forecast(train.label, train.features, test.label, test.features,
                            nrounds = 1000, early_stopping_rounds = 3,
                            max_depth = hyper_grid$max_depth[j], 
                            eta = hyper_grid$eta[j],
                            verbose=0, result.error=TRUE)
        # calc errors
        xgb_test_rmse[j] <- errors$rmse.mean
        xgb_test_mape[j] <- errors$mape.mean
    }

    #ideal hyperparamters
    r <- hyper_grid[which.min(xgb_test_rmse), ]
    if (verbose>0) {
        print(r)
    }
    return(r)
}

# run gridsearch if comb of max_depth & eta g.t. 1
# train length in each fold for gridseatch is max value of h & step
xgb.tsCV <- function (label, features, max_depth = 6, eta = .25,
                      h = 1, window = NULL, initial = 0, step = 1, 
                      count.freq=0.1, ...) 
{
    y <- as.ts(label)
    n <- length(y)
    step <- round(step)
    step_ind <- seq(step, n - 1L, by = step)

    if (initial >= n) 
        stop("initial period too long")

    xreg <- ts(as.matrix(features))
    if (NROW(xreg) != length(y)) 
        stop("features must be of the same size as label")
    tsp(xreg) <- tsp(y)
 
    if (is.null(window)) {
        indx <- seq(1 + initial, n - 1L)
    } else {
        indx <- seq(window + initial, n - 1L, by = 1L)
    }
    indx <- intersect(indx, step_ind)

    e.cols <- c('forecast_start', 'forecast_end', 
                'rmse.mean', 'rmse.sigma', 'mape.mean', 'mape.sigma')
    e <- ts(matrix(NA_real_, nrow = floor(n/step), ncol = length(e.cols)))
    colnames(e) <- e.cols
    
    ###
    hyper_grid <- expand.grid(max_depth = max_depth, eta = eta)
    if (nrow(hyper_grid)>1) {
        hyper_grid.flag <- TRUE
    } else {
        hyper_grid.flag <- FALSE
    }
    

    indx.len <- length(indx)
    by <- round(count.freq*indx.len)
    by <- max(1, by)
    print.when <- seq(0, indx.len, by=by)
    
    cnt <- 0
    
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
        train.label <- subset(y, start=start, end = i)
        train.features <- as.matrix(subset(xreg, start=start, end=i))
        
        # get test data
        start <- i+1
        end <- i+h
        if (end <= nrow(xreg)) {
            test.label <- subset(y, start=start, end=end)
            test.features <- as.matrix(subset(xreg, start=start, end=end))
        } else {
            next
        }
        
        # tune hyperparams
        if (hyper_grid.flag) {
            res <- xgb.gridsearch(train.label, train.features, 
                                  max(h, step), 
                                  hyper_grid)
            max_depth.best <- res$max_depth
            eta.best <- res$eta
        } else {
            max_depth.best <- max_depth
            eta.best <- eta
        }
        
        # train model
        errors <- xgb.forecast(train.label, train.features, test.label, test.features,
                                nrounds = 1000, early_stopping_rounds = 3,
                                max_depth = max_depth.best,
                                eta = eta.best, result.error=TRUE, ...)
        # calc errors
        e[i/step, ] <- c(start, end, errors$rmse.mean, errors$rmse.sigma, 
                                     errors$mape.mean, errors$mape.sigma)

        cnt <- cnt + 1
        if (cnt %in% print.when) {
            message(sprintf("%0.0f %% done.", 100*cnt/length(indx)))
        }
    }
    #return(na.omit(e)) # times of NA kept in e as attr(na.action)
    return(e)
}


xgb.tsCV.mean <- function(label, features, cols=c(1,2,3,5), ...) {
    e <- xgb.tsCV(label, features, ...)
    result <- e[,cols]
    colnames(result) <- c('forecast_start', 'forecast_end', 'rmse', 'mape')
    return(result)
}