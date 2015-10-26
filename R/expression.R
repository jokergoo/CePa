read.cls <-
function (file, treatment, control) {
    label = readLines(file, n = -1)
    label = label[3]
    label = unlist(strsplit(label, " |\t"));
    return(sampleLabel(label = label, treatment = treatment, control = control))
}

read.gct <-
function (file) {
    expr = read.table(file, skip = 2, header = TRUE, sep = "\t", quote = "")
    rownames(expr) = expr[,1]

    checkName = table(expr[,1])
    if(max(checkName) > 1) {
        stop(paste("Genes in gct file should be unique: ", names(which.max(checkName)), sep = " "))
    }
    expr = expr[,-c(1,2)]
    expr = as.matrix(expr)
    
    return(expr)
}


sampleLabel = function (label, treatment, control) {
    if (sum(label == treatment) == 0 | sum(label == control) == 0) {
        stop("Can not find treatment label or control label.")
    }
    res = list(label = label, treatment = treatment, control = control)
    class(res) = "sampleLabel"
    return(res)
}

.treatment = function(sl) {
    return(which(sl$label == sl$treatment))
}

.control = function(sl) {
    return(which(sl$label == sl$control))
}

.permutate = function(sl) {
    sl$label = sample(sl$label, length(sl$label), replace = FALSE)
    return(sl)
}

.factor = function(sl) {
    return(factor(sl$label))
}