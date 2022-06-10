piepoints = function(x, rad, pie, piecol, piebg, ...)
{
    w = par("pin")[1]/diff(par("usr")[1:2])
    h = par("pin")[2]/diff(par("usr")[3:4])
    asp = w/h

    theta = apply(
        pie
        , 1
        , function(p) {
            ang = cumsum((360 * p / sum(p)) * (pi / 180))
            ang = cbind(c(0, ang[-length(ang)]), c(ang[-length(ang)], 2*pi))
            ang
        }
        , simplify=FALSE
    )

    for (i in seq_along(theta))
    {
        xx = x[i, 1]
        yy = x[i, 2]
        th = theta[[i]]
        for (j in 1:nrow(th))
        {
            if ((th[j,2] - th[j,1]) > 0)
            {
                wedges = seq(th[j, 1], th[j, 2], length.out=30)
                xvec = rad[i] * cos(wedges) + xx
                yvec = rad[i] * asp * sin(wedges) + yy
                if (isTRUE(all.equal(unname(th[j,2] - th[j,1]), 2*pi)))
                {
                    polygon(
                        xvec, yvec
                        , border=piecol[j]
                        , col=piebg[j]
                        , ...
                    )
                }
                else
                {
                    polygon(
                        c(xx, xvec), c(yy, yvec)
                        , border=piecol[j]
                        , col=piebg[j]
                        , ...
                    )
                }
            }
            else
                next
        }
    }

}

#' @param x Object of class `prtree`.
#' @param by A grouping factor the same length as the number of
#' rows in the data matrix component of `x`. Alternatively, a numeric
#' matrix with the proportions of each pie.
#' @param piecol A vector of border colors the same length as the 
#' number of levels in the grouping factor `by` (or number of
#' columns if `by` is a matrix).
#' @param piebg A vector of interior colors the same length as the
#' number of levels in the grouping factor `by`.
points.prtree = function(x, by, piecol, piebg, ...)
{
    obj = x
    if (missing(by) || missing(piebg))
    {
        args = list(...)
        cex = if (is.null(args$cex)) 1 else args$cex
        args$cex = cex * colMeans(obj$r)
        args$x = obj$layout
        do.call(points, args)
    }
    else
    {
        if (is.matrix(by))
        {
            stopifnot(nrow(by) == nrow(obj$v))
            stopifnot(length(piebg) == ncol(by))
            pie = sweep(by, 1, rowSums(by), "/")
            with_data = which(!is.nan(rowSums(pie)) & !is.na(rowSums(pie)))
            pie = pie[with_data, ]
        }
        else
        {
            by = as.factor(by)
            stopifnot(length(by) == nrow(obj$x))
            stopifnot(length(piebg) == nlevels(by))
            # map each datum to its nearest graph vertex
            g = predict(obj, obj$x)
            # determine pie proportions of each graph vertex
            pie = table(g, by)
            with_data = as.integer(rownames(pie))
        }

        args = list(...)
        cex = if (is.null(args$cex)) 1 else args$cex
        args$cex = NULL
        cex.vertex = colMeans(obj$r)
        rad = cex * cex.vertex * par("cxy")[2] / pi

        args$x = obj$layout[with_data,]
        args$rad = rad[with_data]
        args$pie = pie
        if (missing(piecol))
            args$piecol = rep(1, length(piebg))
        else
            args$piecol = piecol
        args$piebg = piebg
        do.call(piepoints, args)
    }
}


plot.prtree = function(x, by, piecol, piebg, draw.edges=TRUE, ...) 
{
    obj = x
    
    if (is.null(obj$layout))
        obj$layout = prtree_layout(obj)

    x = obj$layout

    plot(x, type='n', xaxt='n', yaxt='n', bty='n', xlab='', ylab='', ...)

    edge = which(obj$b > 0, arr.ind = TRUE)
    edge = edge[edge[, 1] < edge[, 2], ]

    x0 = x[edge[, 1], 1]
    y0 = x[edge[, 1], 2]
    x1 = x[edge[, 2], 1]
    y1 = x[edge[, 2], 2]

    if (draw.edges)
        segments(x0, y0, x1, y1, ...)
    
    points(obj, by, piecol, piebg, ...)

    invisible(obj)
}
