# p values of all pathways
p.table = function(x, adj.method = NA, cutoff = ifelse(adj.method == "none", 0.01, 0.05)) {

    if(class(x) != "cepa.all") {
        stop("x should be cepa.all object.\n")
    }
    
    n.pathway = length(x)
    p.value = matrix(0, nrow=length(x), ncol= length(x[[1]]))
    for(i in 1:length(x)) {
        p.value[i, ] = sapply(x[[i]], function(x) x$p.value)
    }

    rownames(p.value) = names(x)
    colnames(p.value) = names(x[[1]])
    
    if(!is.na(adj.method)) {
        p.value = apply(p.value, 2, p.adjust, adj.method)
        l = apply(p.value, 1, function(x) { sum(x <= cutoff) > 0})
        p.value = p.value[l, , drop=FALSE]
    }
    
    return(p.value)
}

# return cepa object
get.cepa = function(x, id = NULL, cen = 1) {
    
    if(class(x) != "cepa.all") {
        stop("x should be cepa.all object.\n")
    }
    
    if(is.null(id)) {
        stop("id cannot be null")
    }
    
    if(length(cen) > 1) {
        stop("Length of cen must be equal to 1.\n")
    }
    
    if(is.function(cen)) {
        cen = deparse(substitute(cen))
    }
    else if(mode(cen) == "name") {
        cen = deparse(cen)
    }
    
    return(x[[id]][[cen]])
}
    
