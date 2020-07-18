# ------------------------------------------------------------------------------
# gene categories for simulation of differential distributions
# ------------------------------------------------------------------------------
#   ee = equivalently expressed
#   ep = equivalent proportions
#   de = 'classical' differential expression (shift in means)
#   dp = differential proportions
#   dm = differential modality
#   db = shift in mean & dm
# ------------------------------------------------------------------------------
names(cats) <- cats <- c("ee", "ep", "de", "dp", "dm", "db")
cats <- factor(cats, levels = cats)

# ------------------------------------------------------------------------------
# for ea. cluster-sample, samples the number of cells for 2 groups
# ------------------------------------------------------------------------------
#   n: single numeric or range to sample from
#   k: character vector of cluster IDs
#   s: sample IDs
# > array of dim. #s x #k x 2 with row names = sample IDs & 
#   column names = cluster IDs. Ea. entry is a list of length 2 
#   with numberic values equal to or in the range of n.
# ------------------------------------------------------------------------------
.sample_n_cells <- function(n, k, s) {
    nk <- length(k)
    ns <- length(s)
    if (length(n) == 1) {
        n <- list(rep(n, 2))
    } else {
        n <- replicate(nk * ns, 
            list(sample(seq(n[1], n[2]), 2)))
    }
    matrix(n, 
        nrow = ns, ncol = nk, 
        dimnames = list(s, k))
}

# ------------------------------------------------------------------------------
# generate a randomized data.frame/colData of cell metadata
# (cluster IDs, sample IDs, and group IDs)
# ------------------------------------------------------------------------------
#   n     = nb. of cells
#   ids   = list of IDs to sample from
#   probs = list of probabilities for ea. set of IDs
#           (in order of cluster, sample, group)
# ------------------------------------------------------------------------------
.sample_cell_md <- function(n, ids, probs = NULL) {
    ns <- vapply(ids, length, numeric(1))
    if (is.null(probs)) 
        probs <- vector("list", 3)
    probs <- lapply(seq_along(probs), function(i) {
        if (!is.null(probs[[i]])) {
            return(probs[[i]])
        } else {
            rep(1 / ns[i], ns[i])
        }
    })
    cd <- vapply(seq_along(probs), function(i) 
        sample(ids[[i]], n, TRUE, probs[[i]]), 
        character(n))
    cd <- data.frame(cd, row.names = NULL)
    colnames(cd) <- c("cluster_id", "sample_id", "group_id")
    cd$group_id <- factor(cd$group_id, levels = ids[[3]])
    return(cd)
}

# ------------------------------------------------------------------------------
# helper to sample from a NB across a grid 
# of dispersions ('size') and means ('mu')
# ------------------------------------------------------------------------------
#' @importFrom stats rnbinom
.nb <- function(cs, ds, m, lfcs = NULL, f = 1) {
    ng <- length(ds)
    nc <- length(cs)
    if (is.null(lfcs)) {
        lfcs <- rep(0, ng)
    } else {
        lfcs[lfcs < 0] <- 0
    }
    fcs <- f * (2 ^ lfcs)
    fcs <- rep(fcs, each = nc)
    ss <- rep(1/ds, each = nc)
    ms <- c(t(m[, cs])) * fcs 
    y <- rnbinom(ng * nc, size = ss, mu = ms)
    y <- matrix(y, byrow = TRUE, 
        nrow = ng, ncol = nc, 
        dimnames = list(names(ds), cs))
    ms <- split(ms, rep(seq_len(nrow(m)), each = nc))
    list(counts = y, means = ms)
}

# ------------------------------------------------------------------------------
# helper to simulate differential distributions 
# ------------------------------------------------------------------------------
#   gs:  character vector of genes to simulate from
#   cs:  character vector of cells to simulate from
#   ng1: nb. of cells in 1st group
#   ng2: nb. of cells in 2nd group
#   m:   G x C matrix of NB means
#   d:   numeric vector of dispersions
#   lfc: numeric vector of logFCs
# ------------------------------------------------------------------------------
.sim <- function(cat = unfactor(cats), cs, ms, ds, lfcs, ep, dp, dm) {
    
    cat <- match.arg(cat)
    ncs <- lapply(cs, length)
    
    res <- switch(cat,
        "ee" = {
            list(
                .nb(cs$group1, ds, ms$group1),
                .nb(cs$group2, ds, ms$group2))
        },
        "ep" = {
            g1_hi <- sample(ncs$group1, round(ncs$group1 * ep))
            g2_hi <- sample(ncs$group2, round(ncs$group2 * ep))
            list(
                .nb(cs$group1[-g1_hi], ds, ms$group1),
                .nb(cs$group1[ g1_hi], ds, ms$group1, lfcs), # 50% g1 hi
                .nb(cs$group2[-g2_hi], ds, ms$group2),
                .nb(cs$group2[ g2_hi], ds, ms$group2, lfcs)) # 50% g2 hi
        },
        "de" = {
            list(
                .nb(cs$group1, ds, ms$group1, -lfcs), # lfcs < 0 => all g1 hi
                .nb(cs$group2, ds, ms$group2,  lfcs)) # lfcs > 0 => all g2 hi
        },
        "dp" = {
            props <- sample(c(dp, 1 - dp), 2)
            g1_hi <- sample(ncs$group1, round(ncs$group1 * props[1]))
            g2_hi <- sample(ncs$group2, round(ncs$group2 * props[2]))
            list(                           
                .nb(cs$group1[-g1_hi], ds, ms$group1), 
                .nb(cs$group1[ g1_hi], ds, ms$group1,  lfcs), # lfcs > 0 => dp/(1-dp)% up
                .nb(cs$group2[-g2_hi], ds, ms$group2), 
                .nb(cs$group2[ g2_hi], ds, ms$group2, -lfcs)) # lfcs < 0 => (1-dp)/dp% up
        },
        "dm" = {
            g1_hi <- sample(ncs$group1, round(ncs$group1 * dm))
            g2_hi <- sample(ncs$group2, round(ncs$group2 * dm))
            list(
                .nb(cs$group1[-g1_hi], ds, ms$group1),
                .nb(cs$group1[ g1_hi], ds, ms$group1, -lfcs), # lfcs < 0 => 50% g1 hi
                .nb(cs$group2[-g2_hi], ds, ms$group2),
                .nb(cs$group2[ g2_hi], ds, ms$group2,  lfcs)) # lfcs > 0 => 50% g2 hi
        }, 
        "db" = {
            if (sample(c(TRUE, FALSE), 1)) {
                # all g1 mi, 50% g2 hi
                g2_hi <- sample(ncs$group2, round(ncs$group2 * 0.5))
                list(
                    .nb(cs$group1, ds, ms$group1, abs(lfcs), 0.5),
                    .nb(cs$group2[-g2_hi], ds, ms$group2, -lfcs), 
                    .nb(cs$group2[ g2_hi], ds, ms$group2,  lfcs)) 
            } else {
                # all g2 mi, 50% g1 hi
                g1_hi <- sample(ncs$group1, round(ncs$group1 * 0.5))
                list(
                    .nb(cs$group2, ds, ms$group2, abs(lfcs), 0.5), 
                    .nb(cs$group1[-g1_hi], ds, ms$group1, -lfcs),  
                    .nb(cs$group1[ g1_hi], ds, ms$group1,  lfcs))  
            }
        }
    )
    cs <- map(res, "counts")
    cs <- do.call("cbind", cs)
    ms <- map(res, "means")
    rmv <- vapply(ms, is.null, logical(1))
    ms <- ms[!rmv] %>% 
        map_depth(2, mean) %>% 
        map_depth(1, unlist) %>% 
        do.call(what = cbind)
    ms <- switch(cat, 
        ee = ms,
        de = ms,
        db = if (ncs$group2 == 0) {
            as.matrix(
                ms[, 1])
        } else {
            cbind(
                ms[, 1],
                rowMeans(ms[, c(2, 3)]))
        }, if (ncs$group2 == 0) {
            as.matrix(
                rowMeans(ms[, c(1, 2)]))
        } else {
            cbind(
                rowMeans(ms[, c(1, 2)]),
                rowMeans(ms[, c(3, 4)]))
        })
    rownames(ms) <- rownames(cs)
    colnames(ms) <- names(which(ncs != 0))
    list(cs = cs, ms = ms)
}
