
# wrapper of cepa.ora.all, cepa.univariate.all
# choose corresponding fucntions according to the arguments
# arguments for cepa.ora.all:
#   dif, bk, pc, cen, cen.name, iter
# arguments for cepa.univariate.all:
#   mat, label, pc, cen, cen.name, nlevel, plevel, iter
cepa.all = function(dif = NULL, bk = NULL, mat = NULL, label = NULL, pc, cen = default.centralities,
    cen.name = sapply(cen, function(x) ifelse(mode(x) == "name", deparse(x), x)), 
    nlevel = "tvalue_abs", plevel = "mean", iter = 1000 ) {
    
    # for those who are lazy to specify argument names
    # if the first argument is a vector, then it is ora method
    if(is.vector(dif)) {
        # do nothing
        if(is.null(bk)) {
            dir = system.file(package = "CePa")
            bk = read.table(paste(dir, "/extdata/bk.genome", sep=""), quote = "", stringsAsFactors = FALSE)[[1]]
        }
    } else if(is.matrix(dif)) {
        mat = dif
        dif = NULL
    } else if(is.data.frame(dif)) {
        mat = as.matrix(dif)
        dif = NULL
    }
    
    if(! is.null(dif)) {     # if dif is specified
        res = cepa.ora.all(dif = dif, bk = bk, pc = pc, cen = cen, cen.name = cen.name, iter = iter)
    } else {
        res = cepa.univariate.all(mat = mat, label = label, pc = pc, cen = cen, cen.name = cen.name, nlevel = nlevel, plevel = plevel, iter = iter)
    }
    
    return(res)
}


# wrapper of cepa.ora, cepa.univariate, cepa.multivariate
# choose corresponding fucntions according to the arguments
# arguments for cepa.ora:
#   dif, bk, pc, pathway, id, cen, cen.name, iter
# arguments for cepa.univariate:
#   mat, label, pc, pathway, id, cen, cen.name, nlevel, nlevel, plevel, iter, gene.level, r.gene.level
cepa = function(dif = NULL, bk = NULL, mat = NULL, label = NULL, pc, pathway = NULL, id = NULL, cen = "equal.weight",
    cen.name = if(is.function(cen)) deparse(substitute(cen)) else if(mode(cen) == "name") deparse(cen) else cen,
    nlevel = "tvalue_abs", plevel = "mean", iter = 1000) {
    
    # if the first argument is a vector, then it is ora method
    if(is.vector(dif)) {
        if(is.null(bk)) {
            dir = system.file(package = "CePa")
            bk = read.table(paste(dir, "/extdata/bk.genome", sep=""), quote = "", stringsAsFactors = FALSE)[[1]]
        }
    } else if(is.matrix(dif)) {
        mat = dif
        dif = NULL
    } else if(is.data.frame(dif)) {
        mat = as.matrix(dif)
        dif = NULL
    }
    
    if(! is.null(dif)) {     # if dif is specified
        cat("  Applying Over-representative analysis.\n")
        res = cepa.ora(dif = dif, bk = bk, pc = pc, pathway = pathway, id = id, cen = cen, cen.name = cen.name, iter = iter)
    } else {
        cat("  Applying Gene-set analysis (univariate procedure).\n")
        res = cepa.univariate(mat = mat, label = label, pathway = pathway, pc = pc, id = id, cen = cen,
                   cen.name = cen.name, iter = iter, nlevel = nlevel, plevel = plevel)
    }
    
    return(res)
    
}

is.ora = function(x) {
    return(class(x) == "cepa.all" && x[[1]][[1]]$framework == "ora")
}

is.gsa = function(x) {
    return(class(x) == "cepa.all" && x[[1]][[1]]$framework == "gsa.univariate")
}



divide = function(x, k) {
    if(length(x) ==1 && is.numeric(x)) {
        x = 1:x
    }
    if(length(x) < k) {
        k = length(x)
    }
    w = floor(length(x)/k)
    q = length(x) - k*w
    d = matrix(0, nrow=k, ncol=2)
    n = 1
    for(i in 1:k) {
        d[i, 1] = n
        d[i, 2] = n+w-1+ifelse(q>0, 1, 0)
        n = d[i,2]+1
        q = ifelse(q > 0, q-1, 0)
    }
    d[k,2] = length(x)
    return(d)
}

combine.cepa.all = function(res) {        
    obj = list()
    for(i in 1:length(res)) {
         obj = c(obj, res[[i]])
    }
    class(obj) = "cepa.all"
    return(obj)
}

cepa.all.parallel = function(dif = NULL, bk = NULL, mat = NULL, label = NULL, pc, cen = default.centralities,
    cen.name = sapply(cen, function(x) ifelse(mode(x) == "name", deparse(x), x)), 
    nlevel = "tvalue_abs", plevel = "mean", iter = 1000, ncores = 2) {
    
    if(length(pc$pathList) < ncores) {
        stop("Number of cores should not be larger than the number of pathways.")
    }
    
    cat("Use snow package to do parallel computing (", ncores, "cores ) ...\n")
    cepa.all.parallel.by.snow(dif = dif, bk = bk, mat = mat, label = label, pc = pc, cen = cen,
        cen.name = cen.name, nlevel = nlevel, plevel = plevel, iter = iter, ncores = ncores)

}

cepa.all.parallel.by.snow = function(dif = NULL, bk = NULL, mat = NULL, label = NULL, pc, cen = default.centralities,
    cen.name = sapply(cen, function(x) ifelse(mode(x) == "name", deparse(x), x)), 
    nlevel = "tvalue_abs", plevel = "mean", iter = 1000, ncores = 2) {
    
    
    # to avoid lazy evaluation
    mode(dif)
    mode(bk)
    mode(mat)
    mode(label)
    mode(pc)
    mode(cen)
    mode(cen.name)
    mode(nlevel)
    mode(plevel)
    mode(iter)

    d = divide(1:length(pc$pathList), ncores)
    if(dim(d)[1] < ncores) {
        ncores = dim(d)[1]
        cat("Since your task can only be divided into", ncores, "parts, modify the number of cores to", ncores, "\n")
    }
    
    cl = makeCluster(ncores, type="SOCK")
    ignore = clusterCall(cl, function() {library(CePa); NULL})
    
    res = clusterApply(cl, 1:ncores, function(i) {
                pc = set.pathway.catalogue(pathList = pc$pathList[d[i, 1]:d[i, 2]],
                                           interactionList = pc$interactionList,
                                           mapping = pc$mapping)
                cepa.all(dif = dif, bk = bk, mat = mat, label = label, pc = pc, cen = cen,
                         cen.name = cen.name, nlevel = nlevel, plevel = plevel, iter = iter)
            })
    stopCluster(cl)
    res.p = combine.cepa.all(res)
    return(res.p)
}
