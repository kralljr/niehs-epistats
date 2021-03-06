

#' Function to perform all analyses
#' 
#' @param dat dataframe of outcome y, confounders z, and covariates x
#' @param y name of outcome variable
#' @param confound equation for model of confounding variables
#' @param covar vector of names of covariates of interest
#' @param nfac number of PCA factors
#' @param plot whether to plot results
niehs_outer <- function(dat, y, confound, covar, cp1, nfac = NULL, plot = T, 
  std1 = T, log1 = F, labels1 = NULL) {
    # Get C&RT
    crt0 <- run_crt(dat, y, confound, covar, 0)$tree
    crt1 <- run_crt(dat, y, confound, covar, cp1)

    residxy <- crt1$residxy
    crt1 <- crt1$tree

    if(log1 == T) {
      dat[, covar] <- log(dat[, covar])
    }
    
    # Run pca
    pca1b <- pca1(dat, covar, nfac = nfac)
    sc1 <- pca1b$scores
    nc <- ncol(sc1)
    pcan <- paste0("rotPC", seq(1, nc))
    colnames(sc1) <- pcan
   
    # Get variance explained
    varexp <- sum(pca1b$values[1 : nc]) / length(pca1b$values)
    varexp <- round(varexp * 100, 1)



    # Get variability in y explained by confoundersi
    eqn1 <- paste(y, "~", confound)
    varconf <- summary(lm(eval(eqn1), data = dat))$r.squared
    varconf <- round(varconf * 100, 1)


    # Plot loadings
    pload <- plot_load(pca1b)
  
    # Scale data
    sd1 <- apply(dat[, covar], 2, sd)
    dat[, covar] <- sweep(dat[, covar], 2, sd1, "/")  
    
    # Get full dataset
    full_dat <- data.frame(dat, sc1)
    
    # Do standard regression
    reg1 <- mix_regress(full_dat, y, confound, covar, std = std1)
    
    # Do regression on PCA
    regPCA <- mix_regress(full_dat, y, confound, pcan, std = std1)

    # get groupings
    load1 <- pca1b$load[, 1 : nc]
    pca2 <- pcan[apply(load1, 1, which.max)]
    groupings <- data.frame(rownames(load1), pca2)
    colnames(groupings) <- c("covar", "rotPC")
    
    # Plot regression and PCA results
    preg <- plot_reg(reg1, regPCA, groupings, labels1)
    
    # Get output 
    out <- list(crt0 = crt0, crt1 = crt1, pca1 = pca1b, pload = pload,
        dat = full_dat, reg1 = reg1, regPCA = regPCA,
        groupings = groupings, preg = preg, varexp = varexp,
	varconf = varconf, residxy = residxy, labels1 = labels1, groupings = groupings)
    
    if(plot) {
        fit.control <- rpart.control(xval = 100, cp = 0, 
          minbucket = 5, maxcompete = 4)
        plotcp(crt0)
        fit.control <- rpart.control(xval = 100, cp = cp1, 
          minbucket = 5, maxcompete = 4)
        #par(oma = c(1, 1, 1, 1)) 
        #plot(crt1)
        #text(crt1)
        #fancyRpartPlot(crt1)
	bcol <- brewer.pal(5, "Blues")
        prp(crt1, branch = 1, extra = 1, box.col = bcol[1], 
	    split.box.col = bcol[2], fallen.leaves = T)


	print(pload)
        print(preg$g1)
        fit.control <- rpart.control(xval = 100, cp = 0, 
          minbucket = 5, maxcompete = 4)
    }
    
    return(out)
}


#' Function to perform PCA on exposures data
#' 
#' @param dat dataframe of outcome y, confounders z, and covariates x
#' @param nfac number of latent factors
pca1 <- function(dat, covar, nfac = NULL) {
    # Limit to covariates
    dat <- dat[, covar]
    
    # Find number of factors if not specified
    if(is.null(nfac)) {
        
        # Number of eigenvalues > 1
        nfac <- length(which(prcomp(dat, scale = T)$sdev > 1))
        
    }
    
    # Apply PCA with varimax rotation
    pr1 <- principal(dat, nfactors = nfac)
    
    return(pr1)
}



#' Function to perform C&RT on mixtures data
#' 
#' @param dat dataframe of outcome y, confounders z, and covariates x
#' @param y name of outcome variable
#' @param confound equation for model of confounding variables
#' @param covar vector of names of covariates of interest
run_crt <- function(dat, y, confound, covar, cp1) {
    # Set output
    fit.control <- rpart.control(xval = 100, cp = cp1, minbucket = 5, maxcompete = 4)
    
    # Get matrix of y and covar
    residxy <- dat[, c(y, covar)]

    # First regress out effects of confounders
    confound1 <- paste("~", confound)
    for(i in 1 : ncol(residxy)) {
        eqn1 <- paste(colnames(residxy)[i], confound1)
        resid1 <- lm(eval(eqn1), data  = dat)$resid 
        residxy[, i] <- resid1
    }
    residxy <- data.frame(residxy)


    colnames(residxy) <- paste0("r", colnames(residxy))
    y <- paste0("r", y)
    covar <- paste0("r", covar)


    # Get equation for tree
    eqn2 <- paste(y, "~", paste(covar, collapse = " + "))
    # Find tree
    tree2 <- rpart(eval(eqn2) , data = residxy, control = fit.control)
    
    return(list(tree = tree2, residxy = residxy))
}

# Can use plotcp(run_crt(dat)), plot(run_crt(dat))/ text(run_crt(dat))



#' Standard regression of mixtures
#' 
#' @param dat dataframe of outcome y, confounders z, and covariates x
#' @param y name of outcome variable
#' @param confound equation for model of confounding variables
#' @param covar vector of names of covariates of interest
mix_regress <- function(dat, y, confound, covar, std = T) {
    
    # Specify outcome and confounding
    #eqn1 <- paste(y, "~", confound)
    eqn1 <- paste(y, "~") 


    # Get univariate results
    lmout <- matrix(nrow = length(covar), ncol = 4)
    for(i in 1 : length(covar)) {
        eqn2 <- paste(eqn1, "+", covar[i])
        lm1 <- lm(eqn2, data = dat)
        temp <- lm1 %>% summary
        lmout[i, ] <- temp$coef[covar[i], ]
    }
    colnames(lmout) <- colnames(summary(lm1)$coef)
    rownames(lmout) <- covar
    
    # Get multivariate results
    eqn1 <- paste(y, "~", confound)
    eqn2 <- paste(eqn1, "+", paste(covar, collapse = "+"))
    temp <- lm(eqn2, data = dat) %>% summary
    mlmout <- temp$coef[covar, ]
        
    # Add type of regression and merge
    lmout <- data.frame(lmout)
    mlmout <- data.frame(mlmout)
    lmout <- mutate(lmout, Type = "lm.univar", Variable = rownames(lmout))
    mlmout <- mutate(mlmout, Type = "lm.multivar", Variable = rownames(mlmout))
    lmout <- full_join(lmout, mlmout)
    lmout$Type <- factor(lmout$Type, levels = c("lm.univar", "lm.multivar"))

# Add confidence intervals
    tstat1 <- qt(0.976, nrow(dat) - 1)
    LB <- lmout[, 1] - tstat1 * lmout[, 2]
    UB <- lmout[, 1] + tstat1 * lmout[, 2]
    lmout <- mutate(lmout, LB, UB)
    
    
    # Standarize results (for SD increase)
    if(std) {
        lmout[, 2] <- NA
        # Get std for each variable
        sd <- apply(dat, 2, sd)
        names(sd) <- colnames(dat)
        sd1 <- sd[lmout$Variable]
        
        # Rescale estimate, LB, UB
        lmout[, c("Estimate", "LB", "UB")] <- lmout[, c("Estimate", "LB", "UB")] * sd1
        
    }
    
    return(lmout)
}



#' Function to plot loadings from PCA
#' 
#' @param pca1b results from pca1 function
plot_load <- function(pca1b) {
    # Find number of factors
    nc <- ncol(pca1b$scores)
    
    # Get loadings
    load <- pca1b$load[, 1 : nc]
    
    # Find direction of maximum
    dir1 <- sign(colSums(load))
    load <- sweep(load, 2, dir1, "*")
    
    # Get in format to plot
    colnames(load) <- paste0("rotPC", seq(1 : nc))
    mload <- melt(load)
    colnames(mload)[1 : 2] <- c("variable", "PC")
    
    # Specify colors
    colsload <- c(brewer.pal(8, "Dark2"), brewer.pal(9, "Set1"), 
        brewer.pal(8, "Accent"))
    # Plot
    gload <- ggplot(mload, aes(x = variable, y = value, fill = variable)) + 
        geom_bar(stat = "identity") + theme_bw() + theme(legend.position = "none") +
        scale_fill_manual(guide = "none", values = colsload) + xlab("") +
        ylab("Rotated PCA loadings") +
        facet_wrap(~ PC, nrow = 1) + theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1)) 

    gload
}

#' Wrapper for plotting results from niehs_outer
#' 
#' @param simout results from niehs_outer
#' @param size1 size of lines
#' @param size2 size of points
plot_reg_outer <- function(simout, size1, size2, line1 = T, rmun = F) {
    lmout1 <- simout$reg1
    lmoutPCA <-  simout$regPCA
    groupings <- simout$groupings
    labels1 <- simout$labels1

    plot_reg(lmout1, lmoutPCA, groupings, labels1, size1, size2, line1, rmun)$g1


}


#' Function to plot regression results
#' 
#' @param lmout1 results from mix_regress for standard regression
#' @param lmoutPCA results from mix_regress for PCA
#' @param groupings match lmout1 results to lmoutPCA results
plot_reg <- function(lmout1, lmoutPCA, groupings, labels1 = NULL, size1 = 1.1, size2 = 1.7, line1 = F, rmun = F) {
    
    # number of pcs
    nc <- length(unique(groupings[, 2]))
    
    # Add grouping and model type
    rownames(groupings) <- groupings[, 1]
    gr1 <- groupings[lmout1$Variable, 2]
    lmout1 <- mutate(lmout1, Group = gr1, Model = "regression")
    
    lmoutPCA <- mutate(lmoutPCA, Group = Variable, Model = "PCA")
    
    # Merge data together and add things to plot
    regall <- full_join(lmout1, lmoutPCA)
    regall <- mutate(regall, corsx = 1.5)
    regall <- mutate(regall, coltype = paste0(substr(Type, 4, 4), 
                                              substr(Model, 1, 1)))
    regall$coltype <- factor(regall$coltype, levels = c( "uP", "mP", "ur", "mr"))
    
    
    # Reorder x axis
    regall$Variable <- factor(regall$Variable)
    lev1 <- levels(regall$Variable)
    wh1 <- which(substr(lev1, 1, 2) == "ro")
    wh2 <- which(substr(lev1, 1, 2) != "ro")
    lev1 <- lev1[c(wh1, wh2)]
    regall$Variable <- factor(regall$Variable, levels = lev1)
    
    if(!is.null(labels1)) {
        regall$Variable <- factor(regall$Variable, levels = labels1[, 1], labels = labels1[, 2])
    }

    if(rmun) {

       regall <- filter(regall, Type == "lm.multivar")
    }


    # Plot output
    pd <- position_dodge(.4)
    cols <- brewer.pal(3, "Dark2")[2 : 3]
    g1 <- ggplot(regall, aes(x = Variable, y = Estimate, colour = Model, alpha = Type, shape = Type)) +
        geom_errorbar(aes(ymin = LB, ymax = UB, alpha = Type), size = size1, width = 0,
                      position = pd) +
        geom_point(size = size2, width = 0,
                   position = pd) +
        theme_bw() +     geom_hline(aes(yintercept = 0), size = size1, colour = "grey50", 
                                    linetype = "dashed") +

        scale_color_manual(values = cols, name = "Approach", labels = c("PCA", "Exposures")) +

        scale_shape_manual(values = c(16, 16), guide = "none") + xlab("") +
        ylab("Change per SD increase") + 
        theme(text = element_text(size = 10)) + 
        theme(legend.position = "right") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))


        if(rmun) {
           g1 <- g1 + scale_alpha_manual(guide = F, values = c(1), labels = c( "Adjusted")) 
        }else {
           g1 <- g1 + scale_alpha_manual(values = c(0.5, 1), labels = c("Unadjusted", "Adjusted")) 
}


 
        if(line1) {
          g1 <- g1 + geom_vline(aes(xintercept = corsx), size = size1)
        }

    
    g1 <- g1 + facet_wrap(~ Group, scales = "free_x", nrow = 1)
    

 
        
    list(g1 = g1, data = regall)
    
}




#' Function to get 
#' 
#' @param dat1 results from niehs_outer
#' @param fills matrix of how to fill in, variable, direction, value, order
#' @param cols vector of colors for variables
#' @param lim1 x-axis limits for density plots
#' @param size1 size of text
getden <- function(dat1, fills, cols, lim1, size1 = 20) { 
	# Get residuals
	datres <- dat1$resid
	mres <- melt(datres)

	for(i in 1 : nrow(fills)) { 
		# restrict to one variable
		mres2 <- filter(mres, variable == fills[i, 1])

		# get density
		den1 <- density(mres2$value)


		den2 <- den1
		if(fills[i, 2] == ">=") {

			den2$y <- ifelse(den1$x >= fills[i, 3], den1$y, NA)
		} else {

			den2$y <- ifelse(den1$x < fills[i, 3], den1$y,  NA)
		}

		den2 <- data.frame(den1$x, den1$y, den2$x, den2$y)
		den2 <- mutate(den2, variable = fills[i, 1], type = fills[i, 4]) 
		colnames(den2) <- c("xo", "yo", "x1", "y1", "variable", "type")

		if(i == 1) { 
			datout <- den2
		}else{
			datout <- full_join(datout, den2)
		}
	}


	datout$variable <- factor(datout$variable, levels = fills[, 1])
	datout$type <- factor(datout$type, levels = fills[, 4])


	dat2 <- datout[complete.cases(datout), ]
	dat2 <- dat2[order(dat2$x1), ]

	g1 <- ggplot(datout, aes(x = xo, y = yo))  +

		scale_x_continuous(limits = lim1) +
		geom_line() +
		theme_bw() + xlab("Residuals") + ylab("Density") +
		theme(text = element_text(size = size1)) +

		geom_ribbon(data = dat2, aes(ymin = 0, x = x1, ymax = y1, fill = variable))  +
		scale_fill_manual(values = cols, guide = F) +
		facet_wrap(~ variable, ncol = 1)



	return(list(datout = datout, g1 = g1))
}



