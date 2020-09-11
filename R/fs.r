## -*- coding: utf-8 -*-
##' @name fluorescence.selection
##'
##' @title Ranking genes from fluorescence selection data
##'
##' @description Function for calculation of gene ranks from
##'     fluorescence enriched screening data.
##'
##' @param gene vector of character, gene names.
##' @param sgrna vector of character, sgrna sequences or names.
##' @param barcode vector of character, barcode sequences or names.
##' @param reference vector of integer, raw counts of barcodes in
##'     reference group.
##' @param fluorescence vector of integer, raw counts of barcodes in
##'     fluorescence group.
##' @param selection.ratio numeric, value should be in (0, 1), the
##'     ratio of cells in fluorescence group which were selected from
##'     reference group.
##' @param multiplier float, default to be NA. multiplier for
##'     normalization, left NA for automatic calculation.
##' @param n.to.consider integer, default to be NA. the number of
##'     barcodes to be used for ranking of genes. Genes with more
##'     barcodes than n.to.consider will be ranked by randomly sampled
##'     barcodes (with bootstrapping to make stable results).
##'
##' @return returns the fs data type, which is a list containing gene
##'     ranks, barcode data, and analysis parameters.
##'
##' @examples
##' result <- fluorescence.selection(
##'     gene, sgrna, barcode, reference, fluorescence, selection.ratio=0.01
##' )
##'
##' @export

fluorescence.selection <- function(gene,
                                   sgrna,
                                   barcode,
                                   reference,
                                   fluorescence,
                                   selection.ratio,
                                   multiplier=NA,
                                   n.to.consider=NA)
{
    ## checks
    ## length of input should be the same
    input.length <- c(
        length(gene), length(sgrna),
        length(barcode),
        length(reference), length(fluorescence)
    )
    if (length(unique(input.length)) != 1) {
        stop('The length of gene, sgrna, barcode, reference, and fluorescence vector should be the same.')
    }

    ## selection.ratio should be in (0, 1)

    if (selection.ratio <= 0 | selection.ratio >= 1) {
        stop('selection.ratio should be in (0, 1).')
    }

    ## n.to.consider should be NA or an integer

    if (!(is.na(n.to.consider) | is.numeric(n.to.consider))) {
        stop('n.to.consider should be an integer or default NA value.')
    }

    ## make data.frame
    fs.barcode <- data.frame(
        gene=as.character(gene),
        sgrna=as.character(sgrna),
        barcode=as.character(barcode),
        reference.raw=reference,
        fluorescence.raw=fluorescence,
        stringsAsFactors=FALSE
    )

    ## normalization of data
    norm.factor <- proportion.norm.factor(
        fs.barcode$reference.raw,
        fs.barcode$fluorescence.raw,
        selection.ratio=c(1, selection.ratio)
    )

    if (is.na(multiplier)) {
        multiplier <- max(1/norm.factor)
    }

    fs.barcode$reference.norm <- fs.barcode$reference.raw *
        norm.factor[1] * multiplier

    fs.barcode$fluorescence.norm <- fs.barcode$fluorescence.raw *
        norm.factor[2] * multiplier

    ## binomial distibution parameters
    fs.barcode$binomial.mean <- fs.barcode$reference.norm *
        selection.ratio
    fs.barcode$binomial.sd <-  sqrt(
        fs.barcode$reference.norm *
        selection.ratio * (1 - selection.ratio)
    )

    ## statistic test
    ## using normal test instead of binomial test

    fs.barcode$p <- unlist(
        apply(
            fs.barcode[
                c(
                    'fluorescence.norm',
                    'binomial.mean',
                    'binomial.sd'
                )
            ],
            MARGIN=1,
            function(a) {
                if (a['binomial.mean'] != 0) {
                    pnorm(
                        a['fluorescence.norm'],
                        mean=a['binomial.mean'],
                        sd=a['binomial.sd'],
                        lower.tail=FALSE
                    )
                } else {
                    NA
                }
            }
        )
    )

    fs.barcode$p.rank <- rank(
        fs.barcode$p, ties.method='average', na.last=TRUE
    ) / sum(!is.na(fs.barcode$p))

    ## RRA to gene

    fs.ranks <- split(
        fs.barcode[!is.na(fs.barcode$p), c('gene', 'p.rank')],
        fs.barcode$gene[!is.na(fs.barcode$p)]
    )

    ## n.to.consider

    gene.n <- unlist(Map(nrow, fs.ranks))

    if (is.na(n.to.consider)) {
        n.to.consider <- median(gene.n)
    }

    gene.n.le <- gene.n[gene.n <= n.to.consider]

    gene.n.gt <- gene.n[gene.n > n.to.consider]

    ## for gene with barcodes less than or equal to n.to.consider

    fs.ranks.spread <- as.data.frame(t(as.data.frame(lapply(
        fs.ranks[names(gene.n.le)],
        function(a) {
            c(
                a$p.rank,
                rep(NA, times=n.to.consider)
            )[1:n.to.consider]
        }
    ))))

    fs.gene.rra.le <- apply(
        fs.ranks.spread,
        MARGIN=1,
        rho.score
    )

    ## for gene with barcodes greater than n.to.consider

    fs.gene.rra.samples <- list()

    for (i in seq(100)) {
        fs.ranks.spread <- as.data.frame(t(as.data.frame(lapply(
            fs.ranks[names(gene.n.gt)],
            function(a) {
                sample(
                    a$p.rank,
                    n.to.consider,
                    replace=FALSE
                )
            }
        ))))

        fs.gene.rra.samples[[i]] <- apply(
            fs.ranks.spread,
            MARGIN=1,
            rho.score
        )
    }

    fs.gene.rra.gt <- apply(
        as.data.frame(fs.gene.rra.samples),
        MARGIN=1,
        median
    )

    ## merge rra from two groups
    fs.gene.rra <- sort(
        c(fs.gene.rra.le, fs.gene.rra.gt)
    )

    ## gather up results
    result <- list()
    result$barcode <- fs.barcode
    result$gene <- fs.gene.rra
    result$selection.ratio <- selection.ratio
    result$norm.factor <- norm.factor
    result$multiplier <- multiplier
    result$n.to.consider <- n.to.consider
    class(result) <- 'fs'
    return(result)
}


print.fs <- function(x) {
    cat(
        "Fluorescence selection\n",
        "\n",
        "parameters:\n",
        "    selection ratio: ",
        x$selection.ratio, "\n",
        "    total raw counts: ",
        paste(
            sum(x$barcode$reference.raw),
            sum(x$barcode$fluorescence.raw),
            collapse=', '
        ), "\n",
        "    normalization factor: ",
        paste(x$norm.factor, collapse=', '), "\n",
        "    multiplier: ",
        x$multiplier, "\n",
        "    total normalized counts: ",
        paste(
            sum(x$barcode$reference.norm),
            sum(x$barcode$fluorescence.norm),
            collapse=', '
        ), "\n",
        "    total number of barcodes: ",
        nrow(x$barcode), "\n",
        "    number of barcodes to rank: ",
        sum(x$barcode$binomial.mean != 0), "\n",
        "    total number of genes: ",
        length(unique(x$barcode$gene)), "\n",
        "    number of genes to rank: ",
        length(x$gene), "\n",
        "    number of barcodes to consider per gene: ",
        x$n.to.consider, "\n",
        "\n",
        "top 5 genes:\n",
        sep=""
    )
    print(
        data.frame(
            Rank=seq(5),
            Score=head(x$gene, n=5)
        )
    )
    invisible(x)
}
