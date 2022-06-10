#' Return the index of the vertex that best approximates each
#' row of data in `newdata`
predict.prtree = function(object, newdata, ...)
{
    V = object$v
    x = data.matrix(newdata)
    stopifnot(ncol(x) == ncol(V))
    #t(apply(x, 1, function(y) {
    #    M = sweep(V, 2, y)
    #    r = -diag(M %*% t(M)) / object$sigma
    #    r = exp(r - max(r))
    #    r / sum(r)
    #}))
    apply(x, 1, function(y) {
        M = sweep(V, 2, y)
        which.min(diag(M %*% t(M)))
    })
}