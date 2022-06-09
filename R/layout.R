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
prtree_layout = function(object)
{
    stopifnot(inherits(object, "prtree"))

    # Initial layout
    x = cmdscale(dist(object$v))

    # Fix the last vertex at (0,0)
    x = t(sweep(x, 2, x[nrow(x),]))

    .Call(C_prtree_layout, x, object$b, t(object$v))

    t(x)
}
