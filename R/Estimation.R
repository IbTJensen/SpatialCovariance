library(EstimationTools)
library(data.table)
library(spatstat)
library(Rcpp)

# Kernels
k <- function(r) ifelse(abs(r) < 1, 3/4*(1-r^2), 0)
k_b <- function(r, b) ifelse(abs(r) < b, k(r/b)/b, 0)

table_construct <- function(X, Z){
  N_tau_W <- Z$n
  Neighbours <- crosspairs(X, Z, rmax = Inf)
  Z_dt <- data.table(v_x = Z$x, v_y = Z$y, Z_v = Z$marks)
  info_dt <- data.table(u_x = Neighbours$xi, u_y = Neighbours$yi,
                        v_x = Neighbours$xj, v_y = Neighbours$yj,
                        dist = Neighbours$d)
  info_dt <- merge(info_dt, Z_dt, by = c("v_x", "v_y"))
  setcolorder(info_dt, c("u_x", "u_y", "v_x", "v_y", "Z_v", "dist"))
  info_dt[
    ,e:=edge.Ripley(X = ppp(x = v_x, y = v_y, window = X$window), r = dist)
  ]
  return(info_dt)
}

hatc0 <- function(info_dt, X, Z, r, b){
  N_tau <- Z$n
  lambda <- X$n/area(X$window)

  info_dt <- info_dt[dist > r-b & dist < r+b]
  info_dt[,k:=k_b(r-dist, b)]

  # Note: Markus' kode og besrkivelse bruger 2*pi*r, mens Rasmus'
  # besrkivelse bruger ||u-v||. Sidstenævnte svarer til at skifte
  # 2*pi*r ud med dist nedenfor. Førstenævnte giver dog et estimat i tråd med
  # forventningen, mens sidstenævnte ikke gør.
  # info_dt[,X_sum_terms:=(k*e)/(lambda*(2*pi*r))]
  info_dt[,X_sum_terms:=Z_v*(k*e)/(lambda*(2*pi*r))]
  # info_sum <- info_dt[,.(sum_terms = sum(Z_v*X_sum_terms)), list(v_x, v_y)]
  # c0 <- sum(info_sum$sum_terms)/N_tau
  c0 <- sum(info_dt$X_sum_terms)/N_tau
  return(c0)
}

# sourceCpp("index_cpp.cpp")
Mise_est <- function(info_dt, X, Z, b, R){
  info_dt_R <- info_dt[dist < R]
  info_dt_R <- info_dt_R[order(dist)]
  N_tau <- Z$n
  lambda <- X$n/area(X$window)

  # pair_ppp <- ppp(x = info_dt_R$u_x - info_dt_R$v_x,
  #                 y = info_dt_R$u_y - info_dt_R$v_y,
  #                 window = owin(xrange = c(-R, R), yrange = c(-R, R)))
  #
  # quartets <- closepairs(pair_ppp, rmax = b)
  # qmat <- matrix(NA, nrow = length(quartets$i), ncol = 2)
  # qmat[,1] <- quartets$i
  # qmat[,2] <- quartets$j
  #
  # a <- Sys.time()
  # test <- dist(info_dt_R$dist)
  # Sys.time()-a

  # a <- Sys.time()
  # n <- nrow(info_dt_R)
  # L <- rep(list(NA), n)
  # i <- 1
  # while(i<=n){
  #   # cat(i, "\r")
  #   d <- 0
  #   j <- i
  #   uvi <- info_dt_R$dist[i]
  #   # forwards
  #   while(d < b & j < n){
  #     j <- j+1
  #     d <- info_dt_R$dist[j]-uvi
  #   }
  #   idx_f <- (i+1):(j-1)
  #   # backwards
  #   j <- i
  #   while(d > -b & j > 1){
  #     j <- j-1
  #     d <- info_dt_R$dist[j]-uvi
  #   }
  #   j <- ifelse(j == 1, 0, j)
  #   idx_b <- (j+1):(i-1)
  #   idx_b <- setdiff(idx_b, 0)
  #   L[[i]] <- c(idx_b, idx_f[-1])
  #   i <- i+1
  # }
  # Sys.time()-a

  # idx_dist_pairs <- neighbors_within_band(info_dt_R$dist, b)
  idx_dist_pairs <- neighbors_matrix_condition(as.matrix(info_dt_R[,c(1:4,6)]), b)
  idx_dist_pairs <- lapply(1:length(idx_dist_pairs),
                           function(i) list(dist = info_dt_R$dist[i],
                                            idx = idx_dist_pairs[[i]])
                           )

  # Find all pairs of indicies (i,j), such that dist[i]-dist[j] in (-b,b)
  # a <- Sys.time()
  # idx_dist_pairs <- lapply(
  #   info_dt_R$dist,
  #   function(x){
  #     list(dist = x, idx = which(between(x - info_dt_R$dist, -b, b)))
  #   }
  # )
  # Sys.time()-a

  # a <- Sys.time()
  # idx_dist_pairs <- lapply(
  #   info_dt_R$dist, function(x) which(abs(x-info_dt_R$dist)<b & x!=info_dt_R$dist)
  # )
  # Sys.time()-a

  # a <- Sys.time()
  # idx_dist_pairs2 <- lapply(
  #   1:nrow(info_dt_R),
  #   function(i){
  #     list(dist = info_dt_R$dist[i], idx = which(
  #       # between(info_dt_R$dist[i] - info_dt_R$dist, -b, b) &
  #       abs(info_dt_R$dist[i] - info_dt_R$dist) < b &
  #         ((info_dt_R$u_x[i] != info_dt_R$u_x | info_dt_R$u_y[i] != info_dt_R$u_y) &
  #         (info_dt_R$v_x[i] != info_dt_R$v_x | info_dt_R$v_y[i] != info_dt_R$v_y))
  #     ))
  #   }
  # )
  # idx_dist_pairs3 <- lapply(idx_dist_pairs2, function(x) x$idx)
  # Sys.time()-a
  # Caluclate the terms in the sum to estimate c0 for each dist
  c0_sum_terms <- lapply(
    idx_dist_pairs,
    function(x){
      info_dt_R$Z_v[x$idx]*k_b(x$dist-info_dt_R$dist[x$idx], b)*info_dt_R$e[x$idx]/
        (lambda*2*pi*x$dist)
    }
  )

  # c0_sum_terms <- lapply(
  #   idx_dist_pairs,
  #   function(x){
  #     info_dt_R$Z_v[x]*k_b(x$dist-info_dt_R$dist[x], b)*info_dt_R$e[x]/
  #       (lambda*2*pi*x$dist)
  #   }
  # )
  rm(idx_dist_pairs); gc()

  # sum over all terms to esimate c0^-(u,v)(||u-v||) for every point-pair (u,v)
  c0_cv_list <- lapply(c0_sum_terms, sum)
  c0_cv <- unlist(c0_cv_list)/N_tau
  rm(c0_sum_terms); gc()

  info_dt_R[,c0_cv:=c0_cv]
  Mise_term2 <- sum(info_dt_R$Z_v*info_dt_R$e*info_dt_R$c0_cv/lambda)/N_tau

  # Quadrature for term 1
  hatc <- function(info_dt, X, Z, r, b){
    sapply(r, function(x) hatc0(info_dt, X, Z, x, b))
  }

  hatcsq <- function(r){
    return(hatc(info_dt, X, Z, r, b)^2*r)
  }

  Mise_term1 <- gauss_quad(
    fun = hatcsq,
    lower = 0, upper = R
  )

  return(Mise_term1 - 2*Mise_term2)
}

bandwidth_selection <- function(X, Z, R){
  info_dt <- table_construct(X, Z)
  MISE_est_fct <- function(b) Mise_est(info_dt, X, Z, b, R)

  O <- optim(par = 0.5, fn = MISE_est_fct, lower = 0, upper = 1,
             method = "L-BFGS-B")

  b <- O$par
}
