#' AFSD order relation
#'
#' @description
#' Compares all rows of a numeric matrix or data frame \code{X}
#' with respect to the AFSD order with parameter eps.
#'
#' @usage compOrd2(X,eps)
#'
#' @param X a numeric matrix or a data frame containing numeric or ordered
#'   factor variables with at least two columns and for which there are no distinct
#'   rows that are exactly equal.
#'
#' @param eps a parameter in the open interval (0,0.5)
#'
#' @return A list containing
#'
#' \item{\code{paths}}{a two-column matrix giving all pairs of indices
#'   \code{(i,j)} which satisfy \code{(X[i, ] <=_{a(eps)} X[j, ])}.}
#'
#' \item{\code{colOrder}}{a matrix of the columnwise orders of \code{X}. Used to
#'   compute paths and required for other function calls.}
#'
#' @keywords internal
# C++ function:


# Classical R function

compOrd_cpp <- function(X, eps) {
  if (eps == 'sd'){
    M = new_func_sd(as.matrix(X))
  } else {
    M = new_func(as.matrix(X), eps)
  }

  paths <- which(M==1,arr.ind = TRUE)
  dimnames(paths) <- list(NULL, c("smaller", "greater"))
  list(M=M, paths = paths,eps=eps)
}

compOrd_cpp_eps <- function(X, perc, epsilon) {
  n <- nrow(X)
  list_d_mat = new_func_eps(as.matrix(X))
  d <- list_d_mat$d
  Intvals <- list_d_mat$mat
  ksel <- list_d_mat$k
  Int <- Intvals[1:ksel, 1]
  Int1 <- Intvals[1:ksel, 2]
  Int2 <- Intvals[1:ksel, 3]
  k <- 1
  indx_list <- list_d_mat$indx[1:ksel, ]
  ##########
  sd_perc <- 2*d / (n*(n-1))
  # (1-sd_perc) <= 1e-9
  if (sd_perc == 1){
    M <- list_d_mat$M
    eps <- 'sd'
  } else {
    print('epsilon grid hard coded')

    percentage = rep(0, length(epsilon))
    for (eps in epsilon){
      bool_val <-((Int1 <= eps*Int) + (Int2 <= eps*Int))
      #vec_check[rowcol[bool_val,]] <- TRUE
      percentage[k] <-2*(sum(bool_val)+d) / (n*(n-1))
      k <- k+1
    }
    percentvals <- c(2*d/(n*(n-1)), percentage)
    print(percentvals)
    indx <- length(percentvals) - findInterval(1-perc, sort(1-percentvals))

    # if indx == 0 or indx == 1
    if (indx == 0){
      #std
      M <- list_d_mat$M
      eps <- 'sd'
    } else {
      if (indx == length(percentvals)){
        eps <- epsilon[length(epsilon)]
      } else {
        eps <- epsilon[indx]
      }
      # not if percentage is not reached maximum epsilon is selected
      M <- list_d_mat$M
      eps <- epsilon[indx]
      indx_rows <- indx_list[, 1] + 1
      indx_cols <- indx_list[, 2] + 1
      bool_val <- as.numeric(Int1 - eps*Int <= 1e-9)
      bool_val2 <- as.numeric(Int2 - eps*Int <= 1e-9)
      for (l in 1:length(indx_rows)){
        M[indx_rows[l], indx_cols[l]] <- bool_val[l]
        M[indx_cols[l], indx_rows[l]] <- bool_val2[l]
      }
      # use vec_check to get position of elements which still needs to be checked
    }

  }
  #########

  paths <- which(M==1,arr.ind = TRUE)
  dimnames(paths) <- list(NULL, c("smaller", "greater"))
  list(M=M, paths = paths,eps=eps)

  #paths <- which(M==1,arr.ind = TRUE)
  #dimnames(paths) <- list(NULL, c("smaller", "greater"))
  #list(M=M, paths = paths,eps=eps)
}

compOrd_ecdf <- function(X, grid, eps) {
  if (eps == 'sd'){
    M <- ecdf_comp_sd(X, grid)
  } else {
    M <- ecdf_comp(X, grid, eps)
  }
  paths <- which(M==1,arr.ind = TRUE)
  dimnames(paths) <- list(NULL, c("smaller", "greater"))
  list(M=M, paths = paths,eps=eps)
}

compOrd_ecdf_eps <- function(X, grid, perc, epsilon) {
  n <- nrow(X)
  list_d_mat =  ecdf_comp_eps(as.matrix(X), grid)
  d <- list_d_mat$d
  Intvals <- list_d_mat$mat
  ksel <- list_d_mat$k
  sum_all <- Intvals[1:ksel, 1]
  left_sum <- Intvals[1:ksel, 2]
  right_sum <- Intvals[1:ksel, 3]
  indx_list <- list_d_mat$indx[1:ksel, ]
  ##########
  sd_perc <- 2*d / (n*(n-1))
  if (sd_perc == 1){
    M <- list_d_mat$M
    eps <- 'sd'
  } else {

    print('epsilon grid hard coded')
    #epsilon <- seq(0.01, 0.49, 0.1)
    #epsilon <- c(0.00001, 0.0001, 0.001, epsilon)
    k <- 1
    percentage = rep(0, length(epsilon))
    for (eps in epsilon){
      bool_val <-((left_sum <= eps*sum_all) + (right_sum <= eps*sum_all))
      #vec_check[rowcol[bool_val,]] <- TRUE
      percentage[k] <-2*(sum(bool_val)+d) / (n*(n-1))
      k <- k+1
    }
    percentvals <- c(2*d/(n*(n-1)), percentage)
    print(percentvals)
    indx <- length(percentvals) - findInterval(1-perc, sort(1-percentvals))

    if (indx == 0){
      #std
      M <- list_d_mat$M
      eps <- 'sd'
    } else {
      if (indx == length(percentvals)){
        eps <- epsilon[length(epsilon)]
      } else {
        eps <- epsilon[indx]
      }
      # not if percentage is not reached maximum epsilon is selected
      M <- list_d_mat$M
      eps <- epsilon[indx]
      indx_rows <- indx_list[, 1] + 1
      indx_cols <- indx_list[, 2] + 1
      bool_val <- as.numeric(left_sum - eps*sum_all <= 1e-9)
      bool_val2 <- as.numeric(right_sum - eps*sum_all <= 1e-9)
      for (l in 1:length(indx_rows)){
        M[indx_rows[l], indx_cols[l]] <- bool_val[l]
        M[indx_cols[l], indx_rows[l]] <- bool_val2[l]
      }
      # use vec_check to get position of elements which still needs to be checked
    }
  }
  #########


  paths <- which(M==1,arr.ind = TRUE)
  dimnames(paths) <- list(NULL, c("smaller", "greater"))
  list(M=M, paths = paths,eps=eps)
}

compOrd_normal <- function(X, eps) {
  if (eps == 'sd'){
    M = normal_comp_sd(as.matrix(X))
  } else {
    M = normal_comp(as.matrix(X), eps)
  }

  paths <- which(M==1,arr.ind = TRUE)
  dimnames(paths) <- list(NULL, c("smaller", "greater"))
  list(M=M, paths = paths,eps=eps)
}

compOrd_normal_eps <- function(X, perc, epsilon) {
  n <- nrow(X)
  list_d_mat = normal_comp_eps(as.matrix(X))
  d <- list_d_mat$d
  Intvals <- list_d_mat$mat
  ksel <- list_d_mat$k
  left <- Intvals[1:ksel, 1]
  diff_mean <- Intvals[1:ksel, 2]
  F1 <- Intvals[1:ksel, 3]
  G2 <- Intvals[1:ksel, 4]
  case_ix <- Intvals[1:ksel, 5]
  k <- 1
  indx_list <- list_d_mat$indx[1:ksel, ]
  ##########
  sd_perc <- 2*d / (n*(n-1))
  if (1 == sd_perc){
    M <- list_d_mat$M
    eps <- 'sd'
  } else {
    print('epsilon grid hard coded')
    #epsilon <- seq(0.01, 0.49, 0.1)
    #epsilon <- c(0.00001, 0.0001, 0.001, epsilon)
    percentage = rep(0, length(epsilon))
    bool_val = rep(0, length(diff_mean))
    for (eps in epsilon){
      bool_val1 <-(((1-2*eps)*left <= diff_mean*(eps * (2*F1-1) + (1-F1))) + ((1-2*eps)*left <= -1*diff_mean*(eps + G2*(1-2*eps))))
      bool_val2 <-(((1-2*eps)*left <= diff_mean*(eps + (1-2*eps)*F1)) + ((1-2*eps)*left <= -1*diff_mean*(eps*(2*G2-1) + (1-G2))))
      #vec_check[rowcol[bool_val,]] <- TRUE
      ix1 = case_ix == 0
      ix2 = case_ix == 1
      bool_val[ix1] = bool_val1[ix1]
      bool_val[ix2] = bool_val2[ix2]
      percentage[k] <-2*(sum(bool_val)+d) / (n*(n-1))
      k <- k+1
    }
    percentvals <- c(2*d/(n*(n-1)), percentage)
    print(percentvals)
    indx <- length(percentvals) - findInterval(1-perc, sort(1-percentvals))

    # if indx == 0 or indx == 1
    if (indx == 0){
      #std
      M <- list_d_mat$M
      eps <- 'sd'
    } else {
      if (indx == length(percentvals)){
        eps <- epsilon[length(epsilon)]
      } else {
        eps <- epsilon[indx]
      }
      # not if percentage is not reached maximum epsilon is selected
      M <- list_d_mat$M
      eps <- epsilon[indx]
      indx_rows <- indx_list[, 1] + 1
      indx_cols <- indx_list[, 2] + 1
      # update M with
      bool_val1 <- (1-2*eps)*left <= -1*diff_mean * (eps + (1-2*eps)*G2)#as.numeric(Int1 - eps*Int <= 1e-9)
      bool_val2 <- (1-2*eps)*left <= diff_mean * (eps*(2*F1-1) + (1-F1))#as.numeric(Int2 - eps*Int <= 1e-9)
      bool_val3 <- (1-2*eps)*left <= -1*diff_mean * (eps + (1-2*eps)*F1)#as.numeric(Int1 - eps*Int <= 1e-9)
      bool_val4 <- (1-2*eps)*left <= diff_mean * (eps*(2*G2-1) + (1-G2))

      for (l in 1:length(indx_rows)){
        if (case_ix[l] == 0){
          M[indx_rows[l], indx_cols[l]] <- bool_val1[l]
          M[indx_cols[l], indx_rows[l]] <- bool_val2[l]
        } else {
          M[indx_rows[l], indx_cols[l]] <- bool_val3[l]
          M[indx_cols[l], indx_rows[l]] <- bool_val4[l]
        }
      }
      # use vec_check to get position of elements which still needs to be checked
    }
  }
  #########

  paths <- which(M==1,arr.ind = TRUE)
  dimnames(paths) <- list(NULL, c("smaller", "greater"))
  list(M=M, paths = paths,eps=eps)

}

#' Transitive Reduction of path matrix
#'
#' @description
#' Computes the transitive reduction of the path matrix of a directed acyclic
#' graph.
#'
#' @usage
#' trReduc(paths, n)
#'
#' @param paths a two column integer matrix containing all pairs of nodes
#'   which are linked by a path of edges.
#' @param n the number of nodes.
#'
#' @details
#'
#' For each \code{k}, all indirect paths going through node \code{k} are
#' removed (if there are any).
#'
#' @return
#' A two column matrix giving all pairs of edges \code{(i,j)} such that there
#' exists a directed edge from node \code{i} to node \code{j}.
#'
#' @keywords internal
trReduc <- function(paths, n) {
  # Construct path matrix from node pairs
  edges <- matrix(nrow = n, ncol = n, FALSE)
  edges[paths] <- TRUE
  diag(edges) <- FALSE
  for (k in 1:n) {
    # Remove all indirect paths which go through node k
    edges[edges[, k], edges[k, ]] <- FALSE
  }
  edges <- which(edges, arr.ind = TRUE)
  colnames(edges) <- c("from", "to")
  edges
}

#' Neighbor points with respect to AFSD order
#'
#' @description Find the neighbor points of the rows of a matrix \code{x} within
#' the rows of \code{X} in the AFSD partial order. That is, for each
#' row \code{x[i, ]}, find all indices \code{k} such that either \code{x[i, ] >=
#' X[k, ]} in AFSD and \code{x[i, ] >= X[j, ] >= X[k, ]} holds for no
#' \code{j} different from \code{k}, or \code{x[i, ] <= X[k, ]} in AFSD
#'  and \code{x[i ,] <= X[j, ] <= X[k, ]} holds for no \code{j}
#' different from \code{k}.
#'
#' @usage
#' neighborPoints(x, X, orderX)
#'
#' @param x numeric matrix with at least two columns.
#' @param X numeric matrix with same number of columns as \code{x}.
#' @param eps parameter for the AFSD order, is put the same as in the fit
#' where X is put
#'
#' @return
#' Lists of length \code{nrow(x)} giving for each \code{x[i, ]} the indices
#' of the smaller and the greater neighbor points within the rows of \code{X}.
#'
#'
#'



neighborPoints <- function(x, X,eps, Ord) {

  x <- t(apply(x, 1, FUN=function(x) sort(x)))
  X <- data.matrix(X, rownames.force = NA)
  #d <- length(thresholds)
  mX <- nrow(X)
  mx <- nrow(x)

  if (eps == 'sd'){
    list_indx <- new_func2_sd(X,x)
  } else {
    list_indx <- new_func2(X,x,eps)
  }

  smaller_indx <- list_indx$smaller
  greater_indx <- list_indx$greater

  # compute smaller_indx and greater_indx: lists with mx elements
  # smaller_indx: for each obs find all indx in X that are smaller
  # greater_indx: for each obs find all indx in X that are greater

  predecessors <- successors <- list(0)
  smaller <- greater <- list(0)

  for (i in 1:mx){
    predecessors[[i]] <- smaller[[i]] <- which(smaller_indx[i,] == 1)
    # all predecessors in X, with indices wrt to positions in y
    successors[[i]] <- greater[[i]] <- which(greater_indx[i,] == 1)
    # all successors in X, with indices wrt to positions in y
    Msmallerreduced <- Ord$M[predecessors[[i]],predecessors[[i]],drop=FALSE]
    # what are the order-relations among the predecessors?
    smaller[[i]] <- predecessors[[i]][which(rowSums(Msmallerreduced)==1)]
    # Indices of direct smaller neighbors in y
    Mgreaterreduced <- Ord$M[successors[[i]],successors[[i]], drop=FALSE]
    # what are the order-relations among the successors?
    greater[[i]] <- successors[[i]][which(colSums(Mgreaterreduced)==1)]
  }

  list(smaller = unname(smaller), greater = unname(greater))
  #list(smaller = smaller_indx, greater = greater_indx)
}


neighborPoints3 <- function(x, X,gridx, gridX, eps, Ord) {
  # train object: X is list or X is matrix
  # test object: x
  if (is.vector(gridx, mode = "numeric") && (is.matrix(x))){

    #d <- length(thresholds)
    #mX <- nrow(X)
    mx <- nrow(x)
    if (is.matrix(X)){

      X <- data.matrix(X, rownames.force = NA)
      if (eps == 'sd'){
        list_indx <- new_func_mat_sd(X,x,gridx, gridX)
      } else {
        list_indx <- new_func_mat(X,x,gridx, gridX, eps)
      }

    } else {
      #stop("'data' must be list")
      if (eps == 'sd'){

        list_indx <- new_func_mat_list_sd(X,x,gridx, gridX)
      } else {
        list_indx <- new_func_mat_list(X,x,gridx, gridX, eps)
      }

    }

  } else if (is.list(x) && is.list(gridx)){

      mx <- length(x)
      if (is.list(X)){
        if (eps == 'sd'){
          list_indx <- new_func_list_sd(X,x,gridx, gridX)
        } else {
          list_indx <- new_func_list(X,x,gridx, gridX, eps)
        }
      } else {
        #stop("'data' must be matrix")
        if (eps == 'sd'){
          list_indx <- new_func_list_mat_sd(X,x,gridx, gridX)
        } else {
          list_indx <- new_func_list_mat(X,x,gridx, gridX, eps)
        }
      }
  } else
      stop("wrong input type for 'data' and / or 'grid'")


  #x <- t(apply(x, 1, FUN=function(x) sort(x)))
  #X <- data.matrix(X, rownames.force = NA)
  #d <- length(thresholds)
  #mX <- nrow(X)
  #mx <- nrow(x)


  #list_indx <- new_func_mat(X,x,gridx, gridX, eps)
  smaller_indx <- list_indx$smaller
  greater_indx <- list_indx$greater

  #print(smaller_indx[48,])
  #print(greater_indx[48,])
  # compute smaller_indx and greater_indx: lists with mx elements
  # smaller_indx: for each obs find all indx in X that are smaller
  # greater_indx: for each obs find all indx in X that are greater

  predecessors <- successors <- list(0)
  smaller <- greater <- list(0)

  for (i in 1:mx){

    predecessors[[i]] <- smaller[[i]] <- which(smaller_indx[i,] == 1)
    # all predecessors in X, with indices wrt to positions in y
    successors[[i]] <- greater[[i]] <- which(greater_indx[i,] == 1)
    # all successors in X, with indices wrt to positions in y
    Msmallerreduced <- Ord$M[predecessors[[i]],predecessors[[i]],drop=FALSE]
    # what are the order-relations among the predecessors?
    smaller[[i]] <- predecessors[[i]][which(rowSums(Msmallerreduced)==1)]
    # Indices of direct smaller neighbors in y
    Mgreaterreduced <- Ord$M[successors[[i]],successors[[i]], drop=FALSE]
    # what are the order-relations among the successors?
    greater[[i]] <- successors[[i]][which(colSums(Mgreaterreduced)==1)]
  }

  list(smaller = unname(smaller), greater = unname(greater))
}

neighborPoints_samegrid <- function(x, X,grid, eps, Ord) {
  mx <- nrow(x)
  if (eps == 'sd'){
    list_indx <- new_func_single_grid_sd(X,x,grid)
  } else {
    list_indx <- new_func_single_grid(X,x,grid, eps)
  }


  smaller_indx <- list_indx$smaller
  greater_indx <- list_indx$greater

  predecessors <- successors <- list(0)
  smaller <- greater <- list(0)

  for (i in 1:mx){

    predecessors[[i]] <- smaller[[i]] <- which(smaller_indx[i,] == 1)
    # all predecessors in X, with indices wrt to positions in y
    successors[[i]] <- greater[[i]] <- which(greater_indx[i,] == 1)
    # all successors in X, with indices wrt to positions in y
    Msmallerreduced <- Ord$M[predecessors[[i]],predecessors[[i]],drop=FALSE]
    # what are the order-relations among the predecessors?
    smaller[[i]] <- predecessors[[i]][which(rowSums(Msmallerreduced)==1)]
    # Indices of direct smaller neighbors in y
    Mgreaterreduced <- Ord$M[successors[[i]],successors[[i]], drop=FALSE]
    # what are the order-relations among the successors?
    greater[[i]] <- successors[[i]][which(colSums(Mgreaterreduced)==1)]
  }

  list(smaller = unname(smaller), greater = unname(greater))
  #list(smaller = smaller_indx, greater = greater_indx)
}




# neighborPoints_norm(x = data, X = X, eps=object$eps, Ord = object$Ord)
neighborPoints_norm <- function(x, X, eps, Ord){
  mx <- nrow(x)
  if (eps == 'sd'){
    list_indx <- indx_norm_sd(X, x)
  } else {
    list_indx <- indx_norm(X, x, eps)
  }


  smaller_indx <- list_indx$smaller
  greater_indx <- list_indx$greater


  predecessors <- successors <- list(0)
  smaller <- greater <- list(0)

  for (i in 1:mx){
    predecessors[[i]] <- smaller[[i]] <- which(smaller_indx[i,] == 1)
    # all predecessors in X, with indices wrt to positions in y
    successors[[i]] <- greater[[i]] <- which(greater_indx[i,] == 1)
    # all successors in X, with indices wrt to positions in y
    Msmallerreduced <- Ord$M[predecessors[[i]],predecessors[[i]],drop=FALSE]
    # what are the order-relations among the predecessors?
    smaller[[i]] <- predecessors[[i]][which(rowSums(Msmallerreduced)==1)]
    # Indices of direct smaller neighbors in y
    Mgreaterreduced <- Ord$M[successors[[i]],successors[[i]], drop=FALSE]
    # what are the order-relations among the successors?
    greater[[i]] <- successors[[i]][which(colSums(Mgreaterreduced)==1)]

  }


  list(smaller = unname(smaller), greater = unname(greater))
  #list(smaller = smaller_indx, greater = greater_indx)
}


neighborPoints_old <- function(x, X,eps) {

  nX <- nrow(X)  # earlier   N <- nrow(X)
  nx <- nrow(x)
  oldNames <- names(X)
  x <- data.frame(t(apply(x, 1, FUN=function(x) sort(x))) )
  names(x) <- oldNames
  XX <- rbind(X,x)
  # x can be sorted row wise, since the ensemble members are exchangeable.
  y <- stats::aggregate(x = 1:(nX+nx), by = XX
                        , FUN = identity, simplify = FALSE)
  names(y)[names(y)=="x"] <- "ind"

  Ord <- compOrd_cpp(y[,-which(names(y)=="ind")],eps)

  ind <- y$ind
  unlist_ind <- unlist(ind)
  ls <- rep(seq_along(ind),times=lengths(ind))
  pos <- ls[order(unlist_ind)] # positions in y of XX[i,]
  train_ind <- sort(pos[1:nX]) # position of training data in y

  predecessors <- successors <- list(0)
  smaller <- smaller2 <-greater2 <- greater <- list(0)

  for (i in 1:nx){
    k <- pos[nX+i] # position of x[i,] in y
    tmp1 <- which(Ord$M[,k]==1) # all predecessors in XX
    tmp2 <- which(Ord$M[k,]==1) # all successors in XX
    predecessors[[i]] <- smaller[[i]] <- tmp1[tmp1 %in% train_ind]
    # all predecessors in X, with indices wrt to positions in y
    successors[[i]] <- greater[[i]] <- tmp2[tmp2 %in% train_ind]
    # all successors in X, with indices wrt to positions in y
    Msmallerreduced <- Ord$M[predecessors[[i]],predecessors[[i]],drop=FALSE]
    # what are the order-relations among the predecessors?
    smaller[[i]] <- predecessors[[i]][which(rowSums(Msmallerreduced)==1)]
    # Indices of direct smaller neighbors in y
    Mgreaterreduced <- Ord$M[successors[[i]],successors[[i]], drop=FALSE]
    # what are the order-relations among the successors?
    greater[[i]] <- successors[[i]][which(colSums(Mgreaterreduced)==1)]
    # Indices of direct greater neighbors in y
    if (identical(smaller[[i]],integer(0))==FALSE){
      tmp1 <- unlist(y[smaller[[i]],]$ind)
      # translate in the original indexation in X
      smaller[[i]] <- tmp1[tmp1 <= nX]
    }
    if (identical(greater[[i]],integer(0))==FALSE){
      tmp2 <- unlist(y[greater[[i]],]$ind)
      # translate in the original indexation in X
      greater[[i]] <- tmp2[tmp2 <= nX]
    }
  }

  list(smaller = unname(smaller), greater = unname(greater))
}



find_eps_ecdf <-function(y, X, grid, perc){
  if (is.list(grid)){
    list_d_mat <- ecdf_list_comp_class_eps(X, grid)
  } else {
    list_d_mat <- ecdf_comp_class_eps(X, grid)
  }
  #n <- nrow(X)
  #list_d_mat = new_func_eps(as.matrix(X))

  class_X <- list_d_mat$ind
  x <- data.frame(y = y, ind = seq_along(y))

  split_classes <- stats::aggregate(x = x, by = data.frame(class_X), FUN = identity, simplify = FALSE)

  cpY <- split_classes[["y"]]
  indices <- split_classes[["ind"]]
  class_X_sorted <- split_classes[, 1]
  vec_indx <- sapply(indices, FUN = function(x) x[1])


  d <- list_d_mat$d

  Intvals <- list_d_mat$mat

  ksel <- list_d_mat$k
  Int <- Intvals[1:ksel, 1]
  Int1 <- Intvals[1:ksel, 2]
  Int2 <- Intvals[1:ksel, 3]
  k <- 1
  indx_list <- list_d_mat$indx[1:ksel, ]
  epsilon <- seq(0.01, 0.49, 0.01)
  percentage = rep(0, length(epsilon))
  n <- length(vec_indx)

  for (eps in epsilon){
    bool_val <-((Int1 <= eps*Int) + (Int2 <= eps*Int))
    #vec_check[rowcol[bool_val,]] <- TRUE
    percentage[k] <-2*(sum(bool_val)+d) / (n*(n-1))
    k <- k+1
  }
  percentvals <- c(2*d/(n*(n-1)), percentage)
  indx <- length(percentvals) - findInterval(1-perc, sort(1-percentvals))
  #print(percentvals)

  # if indx == 0 or indx == 1
  if (indx == 0){
    #std
    M <- list_d_mat$M[vec_indx, vec_indx]
    eps <- 'sd'
  } else {
    if (indx == length(percentvals)){
      eps <- epsilon[length(epsilon)]
    } else {
      eps <- epsilon[indx]
    }
    # not if percentage is not reached maximum epsilon is selected
    M <- list_d_mat$M
    eps <- epsilon[indx]
    indx_rows <- indx_list[, 1] + 1
    indx_cols <- indx_list[, 2] + 1
    bool_val <- as.numeric(Int1 - eps*Int <= 1e-9)
    bool_val2 <- as.numeric(Int2 - eps*Int <= 1e-9)
    for (l in 1:length(indx_rows)){
      M[indx_rows[l], indx_cols[l]] <- bool_val[l]
      M[indx_cols[l], indx_rows[l]] <- bool_val2[l]
    }
    # use vec_check to get position of elements which still needs to be checked
    M <- M[vec_indx, vec_indx]
  }
  #paths <- which(M==1,arr.ind = TRUE)
  #dimnames(paths) <- list(NULL, c("smaller", "greater"))
  list(M=M, indices = indices ,eps=eps, vec_indx = vec_indx, cpY = cpY)



}


find_eps_ecdf2 <-function(X, grid, perc){
  if (is.list(grid)){
    list_d_mat <- ecdf_list_comp_class_eps(X, grid)
    n <- length(X)
  } else {
    list_d_mat <- ecdf_comp_class_eps(X, grid)
    n <- nrow(X)
  }

  #list_d_mat = new_func_eps(as.matrix(X))

  #class_X <- list_d_mat$ind
  #x <- data.frame(y = y, ind = seq_along(y))

  #split_classes <- stats::aggregate(x = x, by = data.frame(class_X), FUN = identity, simplify = FALSE)

  #cpY <- split_classes[["y"]]
  #indices <- split_classes[["ind"]]
  #class_X_sorted <- split_classes[, 1]
  #vec_indx <- sapply(indices, FUN = function(x) x[1])


  d <- list_d_mat$d

  Intvals <- list_d_mat$mat

  ksel <- list_d_mat$k
  Int <- Intvals[1:ksel, 1]
  Int1 <- Intvals[1:ksel, 2]
  Int2 <- Intvals[1:ksel, 3]
  k <- 1
  indx_list <- list_d_mat$indx[1:ksel, ]

  epsilon <- seq(0.01, 0.49, 0.01)
  percentage = rep(0, length(epsilon))
  #n <- length(vec_indx)

  for (eps in epsilon){
    bool_val <-((Int1 <= eps*Int) + (Int2 <= eps*Int))
    #vec_check[rowcol[bool_val,]] <- TRUE
    percentage[k] <-2*(sum(bool_val)+d) / (n*(n-1))
    k <- k+1
  }
  percentvals <- c(2*d/(n*(n-1)), percentage)
  indx <- length(percentvals) - findInterval(1-perc, sort(1-percentvals))
  #print(percentvals)

  # if indx == 0 or indx == 1
  if (indx == 0){
    #std
    M <- list_d_mat$M
    eps <- 'sd'
  } else {
    if (indx == length(percentvals)){
      eps <- epsilon[length(epsilon)]
    } else {
      eps <- epsilon[indx]
    }
    # not if percentage is not reached maximum epsilon is selected
    M <- list_d_mat$M
    #eps <- epsilon[indx]
    indx_rows <- indx_list[, 1] + 1
    indx_cols <- indx_list[, 2] + 1
    bool_val <- as.numeric(Int1 - eps*Int <= 1e-9)
    bool_val2 <- as.numeric(Int2 - eps*Int <= 1e-9)
    for (l in 1:length(indx_rows)){
      M[indx_rows[l], indx_cols[l]] <- bool_val[l]
      M[indx_cols[l], indx_rows[l]] <- bool_val2[l]
    }
    # use vec_check to get position of elements which still needs to be checked
  }
  #paths <- which(M==1,arr.ind = TRUE)
  #dimnames(paths) <- list(NULL, c("smaller", "greater"))
  list(M=M ,eps=eps)

}



testing_ind <- function(X){
  m = nrow(X)
  d = ncol(X)
  classY = rep(0, m)
  class_count <- 1
  for (i in 1:(m-1)){
    check_false = 0
    if (classY[i] == 0){
      classY[i] = class_count
      class_count = class_count + 1
      for (j in (i+1):m) {
        F1 <- X[i,]
        F2 <- X[j,]
        if (sum(!F1 == F2) == 0){
          classY[i] = class_count -1
        }
      }
    }
  }
  if (classY[m] == 0){
    classY[m] = class_count
  }
  return(classY)
}



