beta.score <- function(r) {
    r <- sort(r, na.last=TRUE)
    rr <- r[!is.na(r)]
    k <- length(rr)
    return(pbeta(rr, 1:k, k:1))
}

rho.score <- function(r) {
    k <- sum(!is.na(r))
    return(
        min(
            pmax(
                0,
                pmin(beta.score(r) * k, 1, na.rm=TRUE)
            )
        )
    )
}

proportion.norm.factor <- function(...,
                                   selection.ratio)
{
    Call <- match.call(expand.dots=TRUE)
    if (!('selection.ratio' %in% names(Call))) {
        stop('selection.ratio should be set as a vector the same length as previous vectors.')
    }

    selection.ratio <- eval(
        Call[['selection.ratio']],
        envir=parent.frame()
    )

    if (length(Call) - 2 != length(selection.ratio)) {
        stop('The length of selection.ratio should be  the same as previous vectors.')
    }

    Call[['selection.ratio']] <- NULL

    data <- list()

    Call[[1]] <- NULL
    for (i in seq(length(Call))) {
        data[[i]] <- eval(
            Call[[i]],
            envir=parent.frame()
        )
    }

    count.sums <- unlist(Map(sum, data))
    reads.ratio <- count.sums / max(count.sums)
    norm.factor <- selection.ratio / reads.ratio
    norm.factor
}
