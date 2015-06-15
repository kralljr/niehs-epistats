

#' Function to perform all analyses
#' 
#' @param dat dataframe of outcome y, confounders z, and covariates x
#' @param y name of outcome variable
#' @param confound vector of names of confounding variables
#' @param covar vector of names of covariates of interest
#' @param nfac number of PCA factors
#' @param plot whether to plot results
niehs_outer <- function(dat, y, confound, covar, cp1, nfac = NULL, plot = T) {
    # Get C&RT
    crt1 <- run_crt(dat, y, confound, covar, cp1)
    
    # Run pca
    pca1b <- pca1(dat, covar)
    sc1 <- pca1b$scores
    nc <- ncol(sc1)
    pcan <- paste0("rPC", seq(1, nc))
    colnames(sc1) <- pcan
    
    # Plot loadings
    pload <- plot_load(pca1b)
    
    
    # Get full dataset
    full_dat <- data.frame(dat, sc1)
    
    # Do standard regression
    reg1 <- mix_regress(full_dat, y, confound, covar, std = T)
    
    # Do regression on PCA
    regPCA <- mix_regress(full_dat, y, confound, pcan, std = T)
    
    # get groupings
    load1 <- pca1b$load[, 1 : nc]
    pca2 <- pcan[apply(load1, 1, which.max)]
    groupings <- data.frame(rownames(load1), pca2)
    colnames(groupings) <- c("covar", "rPC")
    
    # Plot regression and PCA results
    preg <- plot_reg(reg1, regPCA, groupings)
    
    # Get output 
    out <- list(crt = crt1, pca1 = pca1b, pload = pload,
        dat = full_dat, reg1 = reg1, regPCA = regPCA,
        groupings = groupings, preg = preg)
    
    if(plot) {
        plotcp(crt1)
        plot(crt1)
        text(crt1)
        print(pload)
        print(preg)
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
#' @param confound vector of names of confounding variables
#' @param covar vector of names of covariates of interest
run_crt <- function(dat, y, confound, covar, cp1) {
    
    fit.control <- rpart.control(xval = 100, cp = cp1, minbucket = 5, maxcompete = 4)
    
    # Get matrix of y and covar
    residxy <- dat[, c(y, covar)]

    # First regress out effects of confounders
    confound1 <- paste("~", paste(confound, collapse = "+"))
    for(i in 1 : ncol(residxy)) {
        eqn1 <- paste(colnames(residxy)[i], confound1)
        resid1 <- lm(eval(eqn1), data  = dat)$resid 
        residxy[, i] <- resid1
    }
    residxy <- data.frame(residxy)

    
    # Get equation for tree
    eqn2 <- paste(y, "~", paste(covar, collapse = " + "))
    # Find tree
    tree2 <- rpart(eval(eqn2) , data = residxy, control = fit.control)
    
    return(tree2)
}

# Can use plotcp(run_crt(dat)), plot(run_crt(dat))/ text(run_crt(dat))



#' Standard regression of mixtures
#' 
#' @param dat dataframe of outcome y, confounders z, and covariates x
#' @param y name of outcome variable
#' @param confound vector of names of confounding variables
#' @param covar vector of names of covariates of interest
mix_regress <- function(dat, y, confound, covar, std = T) {
    
    # Specify outcome and confounding
    eqn1 <- paste(y, "~", paste(confound, collapse = "+"))

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
    eqn2 <- paste(eqn1, "+", paste(covar, collapse = "+"))
    temp <- lm(eqn2, data = dat) %>% summary
    mlmout <- temp$coef[covar, ]
        
    # Add type of regression and merge
    lmout <- data.frame(lmout)
    mlmout <- data.frame(mlmout)
    lmout <- mutate(lmout, Type = "lm.univar", Variable = rownames(lmout))
    mlmout <- mutate(mlmout, Type = "lm.multivar", Variable = rownames(mlmout))
    lmout <- full_join(lmout, mlmout)
    
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
    colnames(load) <- paste0("rPC", seq(1 : nc))
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
        facet_wrap(~ PC, nrow = 1) + theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, hjust = 1)) 

    gload
}


#' Function to plot regression results
#' 
#' @param lmout1 results from mix_regress for standard regression
#' @param lmoutPCA results from mix_regress for PCA
#' @param groupings match lmout1 results to lmoutPCA results
plot_reg <- function(lmout1, lmoutPCA, groupings) {
    
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
    regall$coltype <- factor(regall$coltype, levels = c("ur", "mr", "uP", "mP"))
    
    
    # Reorder x axis
    regall$Variable <- factor(regall$Variable)
    lev1 <- levels(regall$Variable)
    wh1 <- which(substr(lev1, 1, 2) == "rP")
    wh2 <- which(substr(lev1, 1, 2) != "rP")
    lev1 <- lev1[c(wh1, wh2)]
    regall$Variable <- factor(regall$Variable, levels = lev1)
    
    
    # Plot output
    pd <- position_dodge(.4)
    cols <- brewer.pal(3, "Dark2")[2 : 3]
    g1 <- ggplot(regall, aes(x = Variable, y = Estimate, colour = Model, alpha = Type, shape = Type)) +
        geom_errorbar(aes(ymin = LB, ymax = UB, alpha = Type), size = 1.2, width = 0,
                      position = pd) +
        geom_point(size = 2.2, width = 0,
                   position = pd) +
        theme_bw() +     geom_hline(aes(yintercept = 0), colour = "grey50", 
                                    linetype = "dashed") +
        geom_vline(aes(xintercept = corsx)) +
        scale_color_manual(values = cols, name = "Approach", labels = c("PCA", "Exposures")) +
        scale_alpha_manual(values = c(0.5, 1), labels = c("Unadjusted", "Adjusted")) +
        scale_shape_manual(values = c(16, 16), guide = "none") + xlab("") +
        ylab("Change per SD increase") + 
        theme(text = element_text(size = 14)) + 
        theme(legend.position = "bottom") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    nc1 <- 
    g1 <- g1 + facet_wrap(~ Group, scales = "free_x", ncol = 3)
    

 
        
    g1
    
}

