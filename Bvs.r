my.Bvs <- function (formula, fixed.cov = c("Intercept"), data, prior.betas = "Robust", 
    prior.models = "Constant", n.keep = 10, time.test = TRUE, 
    priorprobs = NULL) 
{
    result <- list()
    wd <- tempdir()
    unlink(paste(wd, "*", sep = "/"))
    if (!is.null(fixed.cov)) {
        lmfull = lm(formula, data = data, y = TRUE, x = TRUE)
        X.full <- lmfull$x
        namesx <- dimnames(X.full)[[2]]
        if (namesx[1] == "(Intercept)") {
            namesx[1] <- "Intercept"
        }
        response <- strsplit(formula, "~")[[1]][1]
        if (length(fixed.cov) == 1) {
            if ("Intercept" %in% fixed.cov) {
                formulanull = paste(response, "~1", sep = "")
            }
            if (!"Intercept" %in% fixed.cov) {
                formulanull = paste(response, "~-1+", fixed.cov, 
                  sep = "")
            }
        }
        if (length(fixed.cov) > 1) {
            if ("Intercept" %in% fixed.cov) {
                formulanull <- paste(response, "~", paste(fixed.cov[-which(fixed.cov == 
                  "Intercept")], collapse = "+"), sep = "")
            }
            if (!"Intercept" %in% fixed.cov) {
                formulanull <- paste(response, "~", paste(fixed.cov, 
                  collapse = "+"), sep = "")
                formulanull = paste(formulanull, "-1", sep = "")
            }
        }
        lmnull <- lm(formula = formulanull, data = data, y = TRUE, 
            x = TRUE)
        namesnull <- dimnames(lmnull$x)[[2]]
        if (namesnull[1] == "(Intercept)") {
            namesnull[1] <- "Intercept"
        }
        if (length(fixed.cov) < length(namesnull)) {
            if (!"Intercept" %in% namesnull) {
                if ("Intercept" %in% namesx) {
                  stop("When using a factor, Intercept should be included either in both or in non of the models")
                }
            }
            if ("Intercept" %in% namesnull) {
                if (!"Intercept" %in% namesx) {
                  stop("When using a factor, Intercept should be included either in both or in non of the models")
                }
            }
        }
        ausent <- NULL
        j <- 0
        for (i in 1:length(namesnull)) {
            if (!namesnull[i] %in% namesx) {
                ausent <- c(ausent, namesnull[i])
                j <- j + 1
            }
        }
        if (j > 0) {
            stop(paste("object '", ausent, "' not found in the full model\n", 
                sep = ""))
        }
        if (length(namesx) == length(namesnull)) {
            stop("The number of fixed covariates is equal to the number of covariates in the full model. No model selection can be done\n")
        }
        fixed.pos <- which(namesx %in% namesnull)
        n <- dim(data)[1]
        Y <- lmnull$residuals
        X0 <- lmnull$x
        P0 <- X0 %*% (solve(t(X0) %*% X0)) %*% t(X0)
        knull <- dim(X0)[2]
        X1 <- lmfull$x[, -fixed.pos]
        if (dim(X1)[1] < n) {
            stop("NA values found for some of the competing variables")
        }
        X <- (diag(n) - P0) %*% X1
        namesx <- dimnames(X)[[2]]
        if (namesx[1] == "(Intercept)") {
            namesx[1] <- "Intercept"
        }
        p <- dim(X)[2]
        if (n.keep > 2^(p)) {
            warning(paste("The number of models to keep (", n.keep, 
                ") is larger than the total number of models (", 
                2^(p), ") and it has been set to ", 2^(p), sep = ""))
            n.keep <- 2^p
        }
    }
    if (is.null(fixed.cov)) {
        lmfull = lm(formula, data, y = TRUE, x = TRUE)
        X.full <- lmfull$x
        namesx <- dimnames(X.full)[[2]]
        if (namesx[1] == "(Intercept)") {
            namesx[1] <- "Intercept"
        }
        X <- lmfull$x
        knull <- 0
        Y <- lmfull$y
        p <- dim(X)[2]
        n <- dim(X)[1]
        if (n.keep > 2^(p)) {
            warning(paste("The number of models to keep (", n.keep, 
                ") is larger than the total number of models (", 
                2^(p), ") and it has been set to ", 2^(p), sep = ""))
            n.keep <- 2^p
        }
    }
    write(Y, ncolumns = 1, file = paste(wd, "/Dependent.txt", 
        sep = ""))
    write(t(X), ncolumns = p, file = paste(wd, "/Design.txt", 
        sep = ""))
    pfb <- substr(tolower(prior.betas), 1, 1)
    if (pfb != "g" && pfb != "r" && pfb != "z" && pfb != "l" && 
        pfb != "f") 
        stop("I am very sorry: prior for betas no valid\n")
    pfms <- substr(tolower(prior.models), 1, 1)
    if (pfms != "c" && pfms != "s" && pfms != "u") 
        stop("I am very sorry: prior for model space not valid\n")
    if (pfms == "u" && is.null(priorprobs)) {
        stop("A valid vector of prior probabilities must be provided\n")
    }
    if (pfms == "u" && length(priorprobs) != (p + 1)) {
        stop("Vector of prior probabilities with incorrect length\n")
    }
    if (pfms == "u" && sum(priorprobs < 0) > 0) {
        stop("Prior probabilities must be positive\n")
    }
    if (pfms == "u" && priorprobs[1] == 0) {
        stop("Vector of prior probabilities not valid: All the theory here implemented works with the implicit assumption that the null model could be the true model\n")
    }
    else {
        write(priorprobs, ncolumns = 1, file = paste(wd, "/priorprobs.txt", 
            sep = ""))
    }
    if (pfms == "c" | pfms == "s") {
        priorprobs <- rep(0, p + 1)
        write(priorprobs, ncolumns = 1, file = paste(wd, "/priorprobs.txt", 
            sep = ""))
    }
    method <- paste(pfb, pfms, sep = "")
    cat("Info. . . .\n")
    cat("Most complex model has", p + knull, "covariates\n")
    if (!is.null(fixed.cov)) {
        if (knull > 1) {
            cat("From those", knull, "are fixed and we should select from the remaining", 
                p, "\n")
        }
        if (knull == 1) {
            cat("From those", knull, "is fixed and we should select from the remaining", 
                p, "\n")
        }
        cat(paste(paste(namesx, collapse = ", ", sep = ""), "\n", 
            sep = ""))
    }
    cat("The problem has a total of", 2^(p), "competing models\n")
    cat("Of these, the ", n.keep, "most probable (a posteriori) are kept\n")
    if (p > 30) {
        stop("Number of covariates too big. . . consider using GibbsBvs\n")
    }
    estim.time <- 0
    if (time.test && p >= 18) {
        cat("Time test. . . .\n")
        aux <- system.time(result <- switch(method, gc = .C("gConst", 
            as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
            as.integer(2^(p - 1) - 1999), as.integer(2^(p - 1) + 
                2000), as.character(wd), as.double(estim.time), 
            as.integer(knull),PACKAGE="BayesVarSel"), gs = .C("gSB", as.character(""), 
            as.integer(n), as.integer(p), as.integer(4000), as.integer(2^(p - 
                1) - 1999), as.integer(2^(p - 1) + 2000), as.character(wd), 
            as.double(estim.time), as.integer(knull),PACKAGE="BayesVarSel"), gu = .C("gUser", 
            as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
            as.integer(2^(p - 1) - 1999), as.integer(2^(p - 1) + 
                2000), as.character(wd), as.double(estim.time), 
            as.integer(knull),PACKAGE="BayesVarSel"), rc = .C("RobustConst", as.character(""), 
            as.integer(n), as.integer(p), as.integer(4000), as.integer(2^(p - 
                1) - 1999), as.integer(2^(p - 1) + 2000), as.character(wd), 
            as.double(estim.time), as.integer(knull),PACKAGE="BayesVarSel"), rs = .C("RobustSB", 
            as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
            as.integer(2^(p - 1) - 1999), as.integer(2^(p - 1) + 
                2000), as.character(wd), as.double(estim.time), 
            as.integer(knull),PACKAGE="BayesVarSel"), ru = .C("RobustUser", as.character(""), 
            as.integer(n), as.integer(p), as.integer(4000), as.integer(2^(p - 
                1) - 1999), as.integer(2^(p - 1) + 2000), as.character(wd), 
            as.double(estim.time), as.integer(knull),PACKAGE="BayesVarSel"), lc = .C("LiangConst", 
            as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
            as.integer(2^(p - 1) - 1999), as.integer(2^(p - 1) + 
                2000), as.character(wd), as.double(estim.time), 
            as.integer(knull),PACKAGE="BayesVarSel"), ls = .C("LiangSB", as.character(""), 
            as.integer(n), as.integer(p), as.integer(4000), as.integer(2^(p - 
                1) - 1999), as.integer(2^(p - 1) + 2000), as.character(wd), 
            as.double(estim.time), as.integer(knull),PACKAGE="BayesVarSel"), lu = .C("LiangUser", 
            as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
            as.integer(2^(p - 1) - 1999), as.integer(2^(p - 1) + 
                2000), as.character(wd), as.double(estim.time), 
            as.integer(knull),PACKAGE="BayesVarSel"), zc = .C("ZSConst", as.character(""), 
            as.integer(n), as.integer(p), as.integer(4000), as.integer(2^(p - 
                1) - 1999), as.integer(2^(p - 1) + 2000), as.character(wd), 
            as.double(estim.time), as.integer(knull),PACKAGE="BayesVarSel"), zs = .C("ZSSB", 
            as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
            as.integer(2^(p - 1) - 1999), as.integer(2^(p - 1) + 
                2000), as.character(wd), as.double(estim.time), 
            as.integer(knull),PACKAGE="BayesVarSel"), zu = .C("ZSUser", as.character(""), 
            as.integer(n), as.integer(p), as.integer(4000), as.integer(2^(p - 
                1) - 1999), as.integer(2^(p - 1) + 2000), as.character(wd), 
            as.double(estim.time), as.integer(knull),PACKAGE="BayesVarSel"), fc = .C("flsConst", 
            as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
            as.integer(2^(p - 1) - 1999), as.integer(2^(p - 1) + 
                2000), as.character(wd), as.double(estim.time), 
            as.integer(knull),PACKAGE="BayesVarSel"), fs = .C("flsSB", as.character(""), 
            as.integer(n), as.integer(p), as.integer(4000), as.integer(2^(p - 
                1) - 1999), as.integer(2^(p - 1) + 2000), as.character(wd), 
            as.double(estim.time), as.integer(knull)), fu = .C("flsUser", 
            as.character(""), as.integer(n), as.integer(p), as.integer(4000), 
            as.integer(2^(p - 1) - 1999), as.integer(2^(p - 1) + 
                2000), as.character(wd), as.double(estim.time), 
            as.integer(knull),PACKAGE="BayesVarSel")))
        estim.time <- result[[8]] * 2^(p)/(60 * 4000)
        cat("The problem would take ", estim.time, "minutes (approx.) to run\n")
        ANSWER <- readline("Do you want to continue?(y/n) then press enter.\n")
        while (substr(ANSWER, 1, 1) != "n" & substr(ANSWER, 1, 
            1) != "y") {
            ANSWER <- readline("")
        }
        if (substr(ANSWER, 1, 1) == "n") {
            return(NULL)
        }
    }
    cat("Working on the problem...please wait.\n")
    result <- switch(method, gc = .C("gConst", as.character(""), 
        as.integer(n), as.integer(p), as.integer(n.keep), as.integer(1), 
        as.integer(2^(p) - 1), as.character(wd), as.double(estim.time), 
        as.integer(knull),PACKAGE="BayesVarSel"), gs = .C("gSB", as.character(""), 
        as.integer(n), as.integer(p), as.integer(n.keep), as.integer(1), 
        as.integer(2^(p) - 1), as.character(wd), as.double(estim.time), 
        as.integer(knull),PACKAGE="BayesVarSel"), gu = .C("gUser", as.character(""), 
        as.integer(n), as.integer(p), as.integer(n.keep), as.integer(1), 
        as.integer(2^(p) - 1), as.character(wd), as.double(estim.time), 
        as.integer(knull),PACKAGE="BayesVarSel"), rc = .C("RobustConst", as.character(""), 
        as.integer(n), as.integer(p), as.integer(n.keep), as.integer(1), 
        as.integer(2^(p) - 1), as.character(wd), as.double(estim.time), 
        as.integer(knull),PACKAGE="BayesVarSel"), rs = .C("RobustSB", as.character(""), 
        as.integer(n), as.integer(p), as.integer(n.keep), as.integer(1), 
        as.integer(2^(p) - 1), as.character(wd), as.double(estim.time), 
        as.integer(knull),PACKAGE="BayesVarSel"), ru = .C("RobustUser", as.character(""), 
        as.integer(n), as.integer(p), as.integer(n.keep), as.integer(1), 
        as.integer(2^(p) - 1), as.character(wd), as.double(estim.time), 
        as.integer(knull),PACKAGE="BayesVarSel"), lc = .C("LiangConst", as.character(""), 
        as.integer(n), as.integer(p), as.integer(n.keep), as.integer(1), 
        as.integer(2^(p) - 1), as.character(wd), as.double(estim.time), 
        as.integer(knull),PACKAGE="BayesVarSel"), ls = .C("LiangSB", as.character(""), 
        as.integer(n), as.integer(p), as.integer(n.keep), as.integer(1), 
        as.integer(2^(p) - 1), as.character(wd), as.double(estim.time), 
        as.integer(knull),PACKAGE="BayesVarSel"), lu = .C("LiangUser", as.character(""), 
        as.integer(n), as.integer(p), as.integer(n.keep), as.integer(1), 
        as.integer(2^(p) - 1), as.character(wd), as.double(estim.time), 
        as.integer(knull),PACKAGE="BayesVarSel"), zc = .C("ZSConst", as.character(""), 
        as.integer(n), as.integer(p), as.integer(n.keep), as.integer(1), 
        as.integer(2^(p) - 1), as.character(wd), as.double(estim.time), 
        as.integer(knull),PACKAGE="BayesVarSel"), zs = .C("ZSSB", as.character(""), 
        as.integer(n), as.integer(p), as.integer(n.keep), as.integer(1), 
        as.integer(2^(p) - 1), as.character(wd), as.double(estim.time), 
        as.integer(knull),PACKAGE="BayesVarSel"), zu = .C("ZSUser", as.character(""), 
        as.integer(n), as.integer(p), as.integer(n.keep), as.integer(1), 
        as.integer(2^(p) - 1), as.character(wd), as.double(estim.time), 
        as.integer(knull),PACKAGE="BayesVarSel"), fc = .C("flsConst", as.character(""), 
        as.integer(n), as.integer(p), as.integer(n.keep), as.integer(1), 
        as.integer(2^(p) - 1), as.character(wd), as.double(estim.time), 
        as.integer(knull),PACKAGE="BayesVarSel"), fs = .C("flsSB", as.character(""), 
        as.integer(n), as.integer(p), as.integer(n.keep), as.integer(1), 
        as.integer(2^(p) - 1), as.character(wd), as.double(estim.time), 
        as.integer(knull),PACKAGE="BayesVarSel"), fu = .C("flsUser", as.character(""), 
        as.integer(n), as.integer(p), as.integer(n.keep), as.integer(1), 
        as.integer(2^(p) - 1), as.character(wd), as.double(estim.time), 
        as.integer(knull),PACKAGE="BayesVarSel"))
    time <- result[[8]]
    integer.base.b_C <- function(x, k) {
        if (x == 0) 
            return(rep(0, k))
        else {
            ndigits <- (floor(logb(x, base = 2)) + 1)
            res <- rep(0, ndigits)
            for (i in 1:ndigits) {
                res[i] <- (x%%2)
                x <- (x%/%2)
            }
            return(c(res, rep(0, k - ndigits)))
        }
    }
    models <- as.vector(t(read.table(paste(wd, "/MostProbModels", 
        sep = ""), colClasses = "numeric")))
    prob <- as.vector(t(read.table(paste(wd, "/PostProb", sep = ""), 
        colClasses = "numeric")))
    incl <- as.vector(t(read.table(paste(wd, "/InclusionProb", 
        sep = ""), colClasses = "numeric")))
    joint <- as.matrix(read.table(paste(wd, "/JointInclusionProb", 
        sep = ""), colClasses = "numeric"))
    dimen <- as.vector(t(read.table(paste(wd, "/ProbDimension", 
        sep = ""), colClasses = "numeric")))
    betahat <- as.vector(t(read.table(paste(wd, "/betahat", sep = ""), 
        colClasses = "numeric")))
    mod.mat <- as.data.frame(cbind(t(rep(0, (p + 1)))))
    names(mod.mat) <- c(namesx, "prob")
    N <- n.keep
    for (i in 1:N) {
        mod.mat[i, 1:p] <- integer.base.b_C(models[i], p)
        varnames.aux <- rep("", p)
        varnames.aux[mod.mat[i, 1:p] == 1] <- "*"
        mod.mat[i, 1:p] <- varnames.aux
    }
    mod.mat[, (p + 1)] <- prob[]
    inclusion <- incl
    result <- list()
    result$time <- time
    result$lmfull <- lmfull
    if (!is.null(fixed.cov)) {
        result$lmnull <- lmnull
    }
    result$variables <- namesx
    result$n <- n
    result$p <- p
    result$k <- knull
    result$HPMbin <- integer.base.b_C(models[1], (p))
    names(result$HPMbin) <- namesx
    result$modelsprob <- mod.mat
    result$inclprob <- inclusion
    names(result$inclprob) <- namesx
    result$jointinclprob <- data.frame(joint[1:p, 1:p], row.names = namesx)
    names(result$jointinclprob) <- namesx
    result$postprobdim <- dimen
    names(result$postprobdim) <- (0:p) + knull
    result$call <- match.call()
    result$method <- "full"
    class(result) <- "Bvs"
    result
}
#<environment: namespace:BayesVarSel>
