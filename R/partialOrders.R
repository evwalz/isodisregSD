#' SD order relation
#'
#' @description
#' Compares all rows of a numeric matrix or data frame \code{X}
#' with respect to the SD order.
#'
#' @usage compOrd_cpp(X)
#'
#' @param X a numeric matrix or a data frame containing numeric or ordered
#'   factor variables with at least two columns and for which there are no distinct
#'   rows that are exactly equal.
#'
#' @return A list containing
#'
#' \item{\code{paths}}{a two-column matrix giving all pairs of indices
#'   \code{(i,j)} which satisfy \code{(X[i, ] <= X[j, ])}.}
#'
#' \item{\code{colOrder}}{a matrix of the columnwise orders of \code{X}. Used to
#'   compute paths and required for other function calls.}
#'
#' @keywords internal
# C++ function:


# Classical R function

compOrd_cpp <- function(X) {
  M = new_func_sd(as.matrix(X))
  paths <- which(M==1,arr.ind = TRUE)
  dimnames(paths) <- list(NULL, c("smaller", "greater"))
  list(M=M, paths = paths)
}


compOrd_ecdf <- function(X, grid) {
  M <- ecdf_comp_sd(X, grid)
  paths <- which(M==1,arr.ind = TRUE)
  dimnames(paths) <- list(NULL, c("smaller", "greater"))
  list(M=M, paths = paths)
}

compOrd_normal <- function(X) {
  M = normal_comp_sd(as.matrix(X))
  paths <- which(M==1,arr.ind = TRUE)
  dimnames(paths) <- list(NULL, c("smaller", "greater"))
  list(M=M, paths = paths)
}


compOrd_normal_ab <- function(X, a, b) {
  M = normal_comp_sd_ab(as.matrix(X), a, b)
  paths <- which(M==1,arr.ind = TRUE)
  dimnames(paths) <- list(NULL, c("smaller", "greater"))
  list(M=M, paths = paths)
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

#' Neighbor points with respect to SD order
#'
#' @description Find the neighbor points of the rows of a matrix \code{x} within
#' the rows of \code{X} in the SD order. That is, for each
#' row \code{x[i, ]}, find all indices \code{k} such that either \code{x[i, ] >=
#' X[k, ]} in SD and \code{x[i, ] >= X[j, ] >= X[k, ]} holds for no
#' \code{j} different from \code{k}, or \code{x[i, ] <= X[k, ]} in SD
#'  and \code{x[i ,] <= X[j, ] <= X[k, ]} holds for no \code{j}
#' different from \code{k}.
#'
#' @usage
#' neighborPoints(x, X, orderX)
#'
#' @param x numeric matrix with at least two columns.
#' @param X numeric matrix with same number of columns as \code{x}.
#'
#' @return
#' Lists of length \code{nrow(x)} giving for each \code{x[i, ]} the indices
#' of the smaller and the greater neighbor points within the rows of \code{X}.
#'
#'
#'
#' @keywords internal

neighborPoints <- function(x, X, Ord) {

  x <- t(apply(x, 1, FUN=function(x) sort(x)))
  X <- data.matrix(X, rownames.force = NA)
  #d <- length(thresholds)
  mX <- nrow(X)
  mx <- nrow(x)

  list_indx <- new_func2_sd(X,x)

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


neighborPoints3 <- function(x, X,gridx, gridX, Ord) {
  # train object: X is list or X is matrix
  # test object: x
  if (is.vector(gridx, mode = "numeric") && (is.matrix(x))){

    #d <- length(thresholds)
    #mX <- nrow(X)
    mx <- nrow(x)
    if (is.matrix(X)){

      X <- data.matrix(X, rownames.force = NA)

      list_indx <- new_func_mat_sd(X,x,gridx, gridX)

    } else {
      #stop("'data' must be list")
      list_indx <- new_func_mat_list_sd(X,x,gridx, gridX)
    }

  } else if (is.list(x) && is.list(gridx)){

      mx <- length(x)
      if (is.list(X)){
          list_indx <- new_func_list_sd(X,x,gridx, gridX)

      } else {
        #stop("'data' must be matrix")
        list_indx <- new_func_list_mat_sd(X,x,gridx, gridX)
      }
  } else
      stop("wrong input type for 'data' and / or 'grid'")


  #x <- t(apply(x, 1, FUN=function(x) sort(x)))
  #X <- data.matrix(X, rownames.force = NA)
  #d <- length(thresholds)
  #mX <- nrow(X)
  #mx <- nrow(x)

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

neighborPoints_samegrid <- function(x, X,grid,  Ord) {
  mx <- nrow(x)
  list_indx <- new_func_single_grid_sd(X,x,grid)


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




neighborPoints_norm <- function(x, X, Ord){
  mx <- nrow(x)
  list_indx <- indx_norm_sd(X, x)

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





