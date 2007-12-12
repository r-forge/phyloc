"test.linear" <-
function (phylogeny, dataset, variable.list, check.robust = FALSE, 
    cutoff = 3, use.robust = FALSE, kappa = NULL, filter = NULL) 
{
    all.variables <- c(variable.list, filter)
    nc <- length(dataset[, 1]) - 1
    contrasts <- array(data = NA, dim = c(nc, 4))
    variable.list[3] <- paste(variable.list[2], ".squared", sep = "")
    variable.list[4] <- paste(variable.list[2], ".cubed", sep = "")
    dataset[, length(names(dataset)) + 1] <- dataset[, names(dataset) == 
        variable.list[2]]^2
    names(dataset)[length(names(dataset))] <- paste(variable.list[2], 
        ".squared", sep = "")
    dataset[, length(names(dataset)) + 1] <- dataset[, names(dataset) == 
        variable.list[2]]^3
    names(dataset)[length(names(dataset))] <- paste(variable.list[2], 
        ".cubed", sep = "")
    model1 <- pic.lm(phylogeny, dataset, variable.list[1], variable.list[2:4], 
        check.robust = check.robust, cutoff = cutoff, use.robust = use.robust, 
        kappa = kappa)
    model2 <- pic.lm(phylogeny, dataset, variable.list[1], variable.list[2:3], 
        check.robust = check.robust, cutoff = cutoff, use.robust = use.robust, 
        kappa = kappa)
    model3 <- pic.lm(phylogeny, dataset, variable.list[1], variable.list[2], 
        check.robust = check.robust, cutoff = cutoff, use.robust = use.robust, 
        kappa = kappa)
    print(summary.pic.lm(model1))
    print(summary.pic.lm(model2))
    print(summary.pic.lm(model3))
}

