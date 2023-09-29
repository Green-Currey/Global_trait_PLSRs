
print.message <- function(text, more.text=NULL) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message(paste("[", timestamp, "]", text, more.text))
}

brightness.norm <- function(x) {
    x / sqrt(sum(x^2))
}

jackknife.test <- function(plsr.dataset, data.var, n.comps = 15, iterations = 20, prop = 0.05, plots = F) {
    require(pls)
    require(tidyr)
    pls.options(plsralg = "oscorespls")
    pls.options("plsralg")
    
    #lets start a timer to see how long this takes to run
    start.time <- Sys.time()
    
    dims <- dim(plsr.dataset)
    
    #lets create an empty matrix to store our results in
    jk.out <- matrix(data = NA, nrow = iterations, ncol = n.comps) 
    
    
    
    #lets run through all the different iterations of this
    for (i in 1:iterations) {
        
        #remind ourselves what iteration we are on
        # print(paste("Iteration: ", i, sep = ""))
        
        #lets take a sample from our dataset to test this on
        rows <- sample(1:nrow(plsr.dataset), floor(prop*nrow(plsr.dataset)))
        sub.data <- plsr.dataset[rows,]
        
        #lets run our PLSR model now
        plsr.out <- plsr(as.formula(paste0(data.var,'~spectra')), scale = FALSE, 
                         ncomp = n.comps, validation = "LOO",
                         trace = F, data = sub.data)
        
        #lets save our press statistic in our empty matrix
        resPRESS <- as.vector(plsr.out$validation$PRESS)
        jk.out[i,seq(plsr.out$validation$ncomp)] = resPRESS
    }
    
    
    #lets change our output matrix to a dataframe for easier manipulation
    pressDF <- as.data.frame(jk.out)
    
    #lets name the columns
    names(pressDF) <- as.character(seq(n.comps))
    
    #lets write this as a csv for later use
    # write.csv(pressDF, file = paste0(in.var, "_", all, 
    #                                  "_Jackkife_PLSR_Coefficients",
    #                                  today,
    #                                  ".csv"), 
    #           row.names = FALSE)
    
    #lets melt the data for easier plotting
    pressDFres <- pivot_longer(cbind.data.frame(ID = seq(nrow(pressDF)), pressDF), cols = -1, names_to = 'Comp', values_to = 'PRESS')
    pressDFres <- pressDFres %>% mutate(Comp = as.numeric(Comp))
    
    #lets see what our press statistics look like. small is better for this.
    if (plots == T) {
        boxplot(pressDFres$PRESS ~ pressDFres$Comp, 
                xlab = "n Components",
                ylab = "PRESS",
                main = 'PLSR_LMA')
        
        print(Sys.time() - start.time)
    }
    
    return(pressDFres)
    
}

runPLSR <- function (plsr.df, data.var, train.size, band.prefix, jk.test = TRUE,
                     jk.comps = 5, jk.iterations = 30, jk.prop = 0.10, 
                     plots = F, wl = seq(400, 2500, 10)) {
    
    # Debugging
    # plsr.df <- plsr.data
    # data.var <- lma.name
    # train.size <- 2000
    # band.prefix <- 'wave'
    # jk.test <- FALSE
    # jk.comps <- 5
    # jk.iterations <- 30
    # jk.prop <- 0.10
    # plots <- TRUE
    # wl <- seq(400, 2500, 10)
    
    
    require(pls)
    pls.options(plsralg = "oscorespls")
    pls.options("plsralg")
    require(tidyverse)
    
    n <- seq(nrow(plsr.df))
    split <- sample(n, replace = F, size = train.size) 
    
    train <- plsr.df[split, ]
    test <- plsr.df[-split, ]
    
    plsr.spectra <- as.matrix(train[,grep(band.prefix, names(plsr.df))])
    plsr.dataset <- data.frame(data = train[data.var], 
                               spectra = I(plsr.spectra))
    
    test.spectra <- as.matrix(test[,grep(band.prefix, names(plsr.df))])
    test.dataset <- data.frame(data = test[data.var],
                               spectra = I(test.spectra))
    
    
    #lets take a look at the correlations between the spectra and the biochemical data
    
    spectra.cor <- cor(plsr.dataset[data.var], 
                       plsr.dataset$spectra, 
                       use = "complete.obs")
    if (plots == T) {
        #lets take a look at the spectra - due to the naming of the bands in the csv there is not an x axis range
        plot(wl, 
             seq(from = -1, to = 1, length = length(wl)),
             type = "n",
             xlab = "Wavelength (nm)", 
             ylab = "Correlation / Refl (scaled)")
        
        matplot(wl, t(plsr.spectra)*2, 
                type = "l", 
                lty = 1, 
                ylab = "Reflectance (%)", 
                xlab = "Wavelength (nm)",
                add = TRUE)
        
        #lets take a look at the correlation between the spectra and biochemical data
        points(wl,
               spectra.cor, 
               lwd = 4)
        abline(h = 0,lty = 2, lwd = 1.5, col = "grey80")
    }
    
    # now let's figure out which spectra to keep
    spec.cors <- as.data.frame(cbind(wl, t(spectra.cor)))
    names(spec.cors) <- c(band.prefix,"cor")
    
    # # now let's throw out some wavelengths!
    # keep.cors <- which(spec.cors$cor > 0.2 | spec.cors$cor < -0.1)
    # 
    # plsr.dataset$spectra <- plsr.dataset$spectra[,keep.cors]
    
    
    #lets do a jackknife test to find the number of components to include in our PLSR model
    
    #first lets find the dimensions of our dataset and set some parameters
    if (jk.test) {
        print.message('Running jackknife test.')
        jk.df <- jackknife.test(plsr.dataset, data.var, jk.comps, jk.iterations, jk.prop, plots)
        
        # How many components? Can use this to determine if next largest is sig different than lower.  Then lower is best. 
        # We can do this with a simple T-Test
        # a smaller PRESS statistic is better. so lets see where this starts to vary. we want the lowest number of components so that
        # we don't over predict our model.
        loc <- 2
        pval <- 0
        while (pval < 0.1 && loc <= jk.comps) {
            ttest <- t.test(jk.df$PRESS[which(jk.df$Comp == loc-1)], 
                            jk.df$PRESS[which(jk.df$Comp == loc)]); 
            pval <- ttest$p.value;
            loc <- loc+1
        }
        
        print.message(loc-1,'components determined.')
        # since we see a low p-value we can see that there is no difference between the two variables now. so lets go with the smaller value.
        #Now that we know the number of test components lets run our PLSR model again with that number of components.
        nComps <- loc-1
    } else {
        print.message('No jackknife test.')
        print.message('Number of comps = ', paste(jk.comps))
        nComps <- jk.comps
    }
    
    print.message('Running PLSR')
    plsr.out <- plsr(as.formula(paste0(data.var,'~spectra')), scale = FALSE,
                     ncomp = nComps, validation = "LOO",
                     trace = F, data = plsr.dataset)
    if (plots == T) { plot(plsr.out) }
    
    
    if (plots == T) { validationplot(plsr.out) }
    
    #lets save our fitted values
    fit1 <- (plsr.out$fitted.values[, , nComps])
    
    #lets plot them to see what they look like
    
    if (plots == T) {
        par(mar = c(4,4,3,3))
        plot(unlist(train[data.var]), unlist(train[data.var]), 
             type = "n",
             xlab = "Modelled",
             ylab = "Measured",
             main = "PLSR - TRAINING DATA")
        
        abline(lm(unlist(train[data.var]) ~ fit1), lwd = 2)
        abline(0, 1, col = "red", lwd = 2, lty = 2)
        graphics::text(fit1, 
                       unlist(train[data.var]), fit1, 
                       labels = paste0('[',round(train$x),', ', round(train$y),']'),
                       cex = 0.6)
    }
    
    # dev.off()
    
    actual <- unlist(train[data.var])
    predicted <- fit1
    rmse <- sqrt(mean((predicted - actual)^2)) %>% round(3)
    mae <- mean(abs(predicted - actual)) %>% round(3)
    print.message('Training RMSE: ', rmse)
    print.message('Training %RMSE: ', round(rmse/mean(actual)*100,3))
    print.message('Training MAE: ', mae)
    print.message('Training %MAE: ', round(mae/mean(actual)*100,3))
    print.message('Training R2: ', round(summary(lm(unlist(train[data.var]) ~ fit1))$r.squared, 3))
    
    # predict LMA on the test data then plot against measured values
    plsr.predicted.test <- predict(plsr.out, 
                                   ncomp = nComps, 
                                   newdata = test.spectra)
    
    if (plots == T) {
        plot(unlist(train[data.var]), unlist(train[data.var]),
             type = "n",
             xlab = "Modelled",
             ylab = "Measured",
             main = "PLSR - TESTING DATA")
        points(plsr.predicted.test, unlist(test[data.var]))
        abline(lm(unlist(test[data.var]) ~ plsr.predicted.test), lwd = 2)
        abline(0, 1, col = "red", lwd = 2, lty = 2)
    }
    
    actual <- unlist(test[data.var])
    predicted <- plsr.predicted.test
    rmse <- sqrt(mean((predicted - actual)^2)) %>% round(3)
    mae <- mean(abs(predicted - actual)) %>% round(3)
    print.message('Testing RMSE: ', rmse)
    print.message('Testing %RMSE: ', round(rmse/mean(actual)*100,3))
    print.message('Testing MAE: ', mae)
    print.message('Testing %MAE: ', round(mae/mean(actual)*100,3))
    print.message('Testing R2: ',round(summary(lm(unlist(test[data.var]) ~ plsr.predicted.test))$r.squared, 3))


    return(coef(plsr.out, intercept=TRUE) %>% as.vector)
}


wl.interp <- function(y, wavelength) {
    out <- approx(x = seq(400,2500,10), y = y, xout = wavelength, 
                  method = "linear", rule = 2)[[2]]
    
    return(out)
}

trait.map <- function(raster, coeffs, intercept = 0, coeffs_wl, na.rm = F) {
    require(terra)
    
    # raster[raster==0] <- NA
    if (dim(raster)[3] != length(coeffs_wl)) {
        # resample wavelengths 
        raster <- app(raster, function (y, w) wl.interp(y, coeffs_wl))
        raster[raster==0] <- NA
    }
    
    # appling coeffs
    traitmap <- app(raster, function (x) sum(x * coeffs, na.rm) + as.numeric(intercept))
    # traitmap[traitmap<0] <- NA
    
    return(traitmap)    
    
}
