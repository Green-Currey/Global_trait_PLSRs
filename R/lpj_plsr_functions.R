
print.message <- function(text, more.text=NULL) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message(paste("[", timestamp, "]", text, more.text))
}

brightness.norm <- function(x) {
    x / sqrt(sum(x^2))
}

jackknife.test <- function(plsr.dataset, data.var, n.comps = 15, iterations = 20, prop = 0.05, plots = T) {
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

runPLSR <- function (plsr.df, data.var, train.size,
                     jk.comps = 10, jk.iterations = 30, jk.prop = 0.10, 
                     plots = T, wl = seq(400, 2500, 10)) {
    
    require(pls)
    pls.options(plsralg = "oscorespls")
    pls.options("plsralg")
    require(tidyverse)
    
    n <- seq(nrow(plsr.df))
    split <- sample(n, replace = F, size = train.size) 
    
    train <- plsr.df[split,]
    test <- plsr.df[-split,]
    
    plsr.spectra <- as.matrix(train[,7:ncol(train)])
    plsr.dataset <- data.frame(data = train[data.var], 
                               spectra = I(plsr.spectra))
    
    test.spectra <- as.matrix(test[,7:ncol(test)])
    test.dataset <- data.frame(data = test[data.var],
                               spectra = I(test.spectra))
    
    
    #lets take a look at the correlations between the spectra and the biochemical data
    
    #lets take a quick look at the correlations between the spectra and biochemical data
    # Kyla note: I don't understand why this grep fcn is so complicated...
    # spectra.cor <- cor(plsr.spectra, 
    #                    plsr.dataset[grep(in.var, names(plsr.dataset), fixed = TRUE)], 
    #                    use = "complete.obs")
    
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
    names(spec.cors) <- c("wave","cor")
    
    # # now let's throw out some wavelengths!
    # keep.cors <- which(spec.cors$cor > 0.2 | spec.cors$cor < -0.1)
    # 
    # plsr.dataset$spectra <- plsr.dataset$spectra[,keep.cors]
    
    
    #lets do a jackknife test to find the number of components to include in our PLSR model
    
    #first lets find the dimensions of our dataset and set some parameters
    print.message('Running jackknife test.')
    jk.df <- jackknife.test(plsr.dataset, data.var, jk.comps, jk.iterations, jk.prop)
    dim(jk.df)
    
    # I have no idea how to determine what works here so I'm just guessing. Need to ask Shawn or reread his paper.
    
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
    # By examining the out put we can determine what the best number of components are to avoid overfitting. Need to ask Shawn about this.
    
    
    # since we see a low p-value we can see that there is no difference between the two variables now. so lets go with the smaller value.
    #Now that we know the number of test components lets run our PLSR model again with that number of components.
    
    nComps <- loc-1
    
    print.message('Running PLSR')
    plsr.out <- plsr(as.formula(paste0(data.var,'~spectra')), scale = FALSE,
                     ncomp = nComps, validation = "LOO",
                     trace = F, data = plsr.dataset)
    if (plots == T) { plot(plsr.out) }
    
    # plot(plsr.dataset$PLSR_N, plsr.out$fitted.values[,,1], type = "n")
    # graphics::text(plsr.dataset$PLSR_N, plsr.out$fitted.values[,,1], 
    #                labels = upper.data$pointIDs, 
    #                cex = 0.8)
    
    if (plots == T) { validationplot(plsr.out) }
    
    #lets save our fitted values
    fit1 <- (plsr.out$fitted.values[, , nComps])
    
    #lets plot them to see what they look like
    # png(filename = paste0(fig.dir, "toc_LMA_obs_pred_", today, ".png"),
    #     width = 7, height = 7, units = "in", res = 300)
    
    if (plots == T) {
        par(mar = c(4,4,3,3))
        plot(unlist(train[data.var]), unlist(train[data.var]), 
             type = "n",
             xlab = "Modelled LMA (g/m2)",
             ylab = "Measured LMA (g/m2)",
             main = "PLSR Modelled Top of Canopy LMA - TRAINING DATA")
        
        abline(lm(unlist(train[data.var]) ~ fit1), lwd = 2)
        abline(0, 1, col = "red", lwd = 2, lty = 2)
        graphics::text(fit1, 
                       unlist(train[data.var]), fit1, 
                       labels = paste0('[',round(train$x),', ', round(train$y),']'),
                       cex = 0.6)
    }
    
    # dev.off()
    
    
    print.message('Training RMSE: ', round(sqrt(mean((unlist(train[data.var]) - fit1)^2)),3))
    print.message('Training R2: ',round(summary(lm(unlist(train[data.var]) ~ fit1))$r.squared, 3))
    
    # predict LMA on the test data then plot against measured values
    plsr.predicted.test <- predict(plsr.out, 
                                   ncomp = nComps, 
                                   newdata = test.spectra)
    
    if (plots == T) {
        plot(unlist(train[data.var]), unlist(train[data.var]),
             type = "n",
             xlab = "Modelled LMA (g/m2)",
             ylab = "Measured LMA (g/m2)",
             main = "PLSR Modelled Top of Canopy LMA - TEST DATA")
        points(plsr.predicted.test, unlist(test[data.var]))
        abline(lm(unlist(test[data.var]) ~ plsr.predicted.test), lwd = 2)
        abline(0, 1, col = "red", lwd = 2, lty = 2)
    }
    
    print.message('Testing RMSE: ', round(sqrt(mean((unlist(test[data.var]) - plsr.predicted.test)^2)),3))
    print.message('Testing R2: ',round(summary(lm(unlist(test[data.var]) ~ plsr.predicted.test))$r.squared, 3))
    
    return(coef(plsr.out, intercept=TRUE) %>% as.vector)
}

applyCoeff <- function(spectra, coeffs, intercept = 0, scale = 1) {
    temp <- sum(spectra * as.vector(coeffs) * scale)
    trait <- temp + as.numeric(intercept)
    return(as.numeric(trait))
}

wl.interp <- function(y, wavelength) {
    out <- approx(x = seq(400,2500,10), y = y, xout = wavelength, 
                  method = "linear", rule = 2)[[2]]
    
    return(out)
}

trait.map <- function(raster, coeffs, coeffs_wl) {
    require(terra)
    raster[raster==0] <- NA
    if (dim(raster)[3] != length(coeffs_wl)) {
        
        # resample wavelengths 
        raster <- app(raster, function (y, w) wl.interp(y, coeffs_wl))
        raster[raster==0] <- NA
        
    }
    
    # appling coeffs
    traitmap <- app(raster, function (x) applyCoeff(x, coeffs = coeffs[-1], intercept = coeffs[1]))
    traitmap[traitmap<0] <- NA
    
    return(traitmap)    
    
}
