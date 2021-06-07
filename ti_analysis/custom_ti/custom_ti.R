#!/scratch/opt/R/3.6.0/bin/Rscript

library(dplyr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)
library(dynwrap, warn.conflicts = FALSE)

requireNamespace("princurve", quietly = TRUE)
requireNamespace("cluster", quietly = TRUE)
requireNamespace("irlba", quietly = TRUE)

suppressWarnings(library(slingshot, warn.conflicts = FALSE))

message(sprintf('slingshot version: %s', packageVersion('slingshot')))

run_fun <- function(expression, parameters, priors, verbose) {
  start_id <- priors$start_id
  end_id <- priors$end_id
  dimred <- priors$dimred
  groups_id <- priors$groups_id
  
  start_cell <- if (!is.null(start_id)) { sample(start_id, 1) } else { NULL }
  
  checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))
  
  if (is.null(dimred)) {
    ndim <- parameters$ndim
    if (ncol(expression) <= ndim) {
      message(paste0(
        "ndim is ", ndim, " but number of dimensions is ", ncol(expression),
        ". Won't do dimensionality reduction."
      ))
      rd <- as.matrix(expression)
    } else {
      pca <- irlba::prcomp_irlba(expression, n = ndim)
      if (ndim > 3) {
        x <- 1:ndim
        optpoint1 <- which.min(sapply(2:10, function(i) {
          x2 <- pmax(0, x - i)
          sum(lm(pca$sdev[1:ndim] ~ x + x2)$residuals^2 * rep(1:2,each = 10))
        }))
        x <- cbind(1:ndim, pca$sdev[1:ndim])
        line <- x[c(1, nrow(x)),]
        proj <- princurve::project_to_curve(x, line)
        optpoint2 <- which.max(proj$dist_ind)-1
        optpoint <- max(c(min(c(optpoint1, optpoint2)), 3))
      } else {
        optpoint <- ndim
      }
      
      rd <- pca$x[, seq_len(optpoint)]
      rownames(rd) <- rownames(expression)
    }
  } else {
    message("Using given dimred")
    rd <- dimred
  }
  if (is.null(groups_id)) {
    max_clusters <- min(nrow(expression)-1, 10)
    if (parameters$cluster_method == "pam") {
      if (nrow(rd) > 10000) {
        warning("PAM (the default clustering method) does not scale well to a lot of cells. You might encounter memory issues. This can be resolved by using the CLARA clustering method, i.e. cluster_method = 'clara'.")
      }
      clusterings <- lapply(3:max_clusters, function(K){
        cluster::pam(rd, K) # we generally prefer PAM as a more robust alternative to k-means
      })
    } else if (parameters$cluster_method == "clara") {
      clusterings <- lapply(3:max_clusters, function(K){
        cluster::clara(rd, K) # we generally prefer PAM as a more robust alternative to k-means
      })
    }
    wh.cl <- which.max(sapply(clusterings, function(x){ x$silinfo$avg.width })) + 1
    labels <- clusterings[[min(c(wh.cl, 8))]]$clustering
  } else {
    message("Using given groups/clustering")
    labels <- groups_id %>% deframe()
  }
  
  start.clus <-
    if(!is.null(start_cell)) {
      labels[[start_cell]]
    } else {
      NULL
    }
  end.clus <-
    if(!is.null(end_id)) {
      unique(labels[end_id])
    } else {
      NULL
    }
  sds <- slingshot::slingshot(
    rd,
    labels,
    start.clus = start.clus,
    end.clus = end.clus,
    shrink = parameters$shrink,
    reweight = parameters$reweight,
    reassign = parameters$reassign,
    thresh = parameters$thresh,
    maxit = parameters$maxit,
    stretch = parameters$stretch,
    smoother = parameters$smoother,
    shrink.method = parameters$shrink.method
  )
  
  start_cell <- apply(slingshot::slingPseudotime(sds), 1, min) %>% sort() %>% head(1) %>% names()
  start.clus <- labels[[start_cell]]
  
  checkpoints$method_aftermethod <- as.numeric(Sys.time())
  from <- to <- NULL
  lineages <- slingLineages(sds)
  lineage_ctrl <- slingParams(sds)
  cluster_network <- lineages %>%
    map_df(~ tibble(from = .[-length(.)], to = .[-1])) %>%
    unique() %>%
    mutate(
      length = lineage_ctrl$dist[cbind(from, to)],
      directed = TRUE
    )
  dimred <- reducedDim(sds)
  cluster <- slingClusterLabels(sds)
  adj <- slingAdjacency(sds)
  lin_assign <- apply(slingCurveWeights(sds), 1, which.max)
  progressions <- map_df(seq_along(lineages), function(l) {
    ind <- lin_assign == l
    lin <- lineages[[l]]
    pst.full <- slingPseudotime(sds, na = FALSE)[,l]
    pst <- pst.full[ind]
    means <- sapply(lin, function(clID){
      stats::weighted.mean(pst.full, cluster[,clID])
    })
    non_ends <- means[-c(1,length(means))]
    edgeID.l <- as.numeric(cut(pst, breaks = c(-Inf, non_ends, Inf)))
    from.l <- lineages[[l]][edgeID.l]
    to.l <- lineages[[l]][edgeID.l + 1]
    m.from <- means[from.l]
    m.to <- means[to.l]
    
    pct <- (pst - m.from) / (m.to - m.from)
    pct[pct < 0] <- 0
    pct[pct > 1] <- 1
    
    tibble(cell_id = names(which(ind)), from = from.l, to = to.l, percentage = pct)
  })
  output <-
    dynwrap::wrap_data(
      cell_ids = rownames(expression) # wrap_data$expression
    ) %>%
    dynwrap::add_trajectory(
      milestone_network = cluster_network,
      progressions = progressions
    ) %>%
    dynwrap::add_dimred(
      dimred = dimred
    ) %>%
    dynwrap::add_timings(checkpoints)
  output$sling_out <- list()
  output$sling_out$pseudotime <- slingshot::slingPseudotime(sds) %>% as.data.frame()
  output$sling_out$cell_weights <- slingshot::slingCurveWeights(sds) %>% as.data.frame()
  save(sds, file='sds.RData')
  output
}