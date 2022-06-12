#' Layout the principle tree in the plane
#'
#' Uses stress majorization to compute an (x,y) coordinate
#' for each vertex.
#'
#' @references 
#' Graph Drawing by Stress Majorization.
#' Emden R. Gansner, Yehuda Koren and Stephen North.
#' In Proceedings 12th Symposium on Graph Drawing (GD), 
#' pages 239–250, 2004.
prtree_layout = function(object, maxit=10000L, verbose=TRUE)
{
    stopifnot(inherits(object, "prtree"))

    # Initial layout
    if (nrow(object$v) > 2)
        x = cmdscale(dist(object$v))
    else if (nrow(object$v) == 2)
        x = rbind(rnorm(2),c(0,0))
    else
        return (matrix(0, 1, 2))

    # Fix the last vertex at (0,0)
    x = t(sweep(x, 2, x[nrow(x),]))

    .Call(
        C_prtree_layout
        , x
        , object$b
        , t(object$v)
        , as.integer(maxit)
        , as.integer(verbose)
    )

    t(x)
}
