my_clhs=function (x, size, include = NULL, cost = NULL, iter = 10000, 
          temp = 1, tdecrease = 0.95, weights = list(numeric = 1, 
                                                     factor = 1, correlation = 1), eta = 1, obj.limit = -Inf, 
          length.cycle = 10, simple = TRUE, progress = TRUE, track = NULL) 
{
  if (tdecrease >= 1) 
    stop("tdecrease should be < 1")
  if (is.null(cost)) {
    cost_mode <- FALSE
    if (!is.null(track)) {
      if (is.numeric(track)) 
        i_cost <- track
      else i_cost <- which(names(x) == track)
      cost <- x[, i_cost, drop = FALSE]
      x <- x[, -1 * i_cost, drop = FALSE]
      cost_mode <- TRUE
      track_mode <- TRUE
    }
    else {
      track_mode <- FALSE
    }
  }
  else {
    if (is.numeric(cost)) 
      i_cost <- cost
    else i_cost <- which(names(x) == cost)
    if (!length(i_cost)) 
      stop("Could not find the cost attribute.")
    cost <- x[, i_cost, drop = FALSE]
    x <- x[, -1 * i_cost, drop = FALSE]
    if (!is.null(include)) {
      cost[include, ] <- 0
    }
    cost_mode <- TRUE
    track_mode <- FALSE
  }
  if (!is.null(include)) {
    if (size <= length(include)) {
      stop(paste0("size (", size, ") should be larger than length of include (", 
                  length(include), ")"))
    }
  }
  i_factor <- which(!sapply(x, is.numeric))
  n_factor <- length(i_factor)
  if (n_factor > 0) {
    data_continuous <- x[, -1 * i_factor, drop = FALSE]
    data_factor <- x[, i_factor, drop = FALSE]
    factor_levels <- apply(data_factor, 2, function(x) {
      ifelse(is.factor(x), res <- levels(x), res <- levels(factor(x)))
      res
    })
  }
  else {
    data_continuous <- x
  }
  metropolis <- exp(-1 * 0/temp)
  n_data <- nrow(data_continuous)
  continuous_strata <- apply(data_continuous, 2, function(x) {
    seq(min(x, na.rm=T), max(x, na.rm = T), length.out = size+1)
  })
  
  cor_mat <- cor(data_continuous, use = "complete.obs")
  if (n_factor == 0) 
    data_factor_sampled <- data.frame(stringsAsFactors = TRUE)
  else factor_obj <- alply(data_factor, 2, function(x) table(x)/n_data)
  sampled_size <- size - length(include)
  not_included <- setdiff(1:n_data, include)
  n_remainings <- n_data - size
  i_sampled <- c(sample(not_included, size = sampled_size, 
                        replace = FALSE), include)
  i_unsampled <- setdiff(1:n_data, i_sampled)
  data_continuous_sampled <- data_continuous[i_sampled, , 
                                             drop = FALSE]
  if (n_factor > 0) 
    data_factor_sampled <- data_factor[i_sampled, , drop = FALSE]
  res <- .lhs_obj(size = size, data_continuous_sampled = data_continuous_sampled, 
                  data_factor_sampled = data_factor_sampled, continuous_strata = continuous_strata, 
                  cor_mat = cor_mat, factor_obj = factor_obj, weights = weights, 
                  eta = eta)
  obj <- res$obj
  delta_obj_continuous <- res$delta_obj_continuous
  if (cost_mode) {
    op_cost <- sum(cost[i_sampled, ])
    op_cost_values <- vector(mode = "numeric", length = iter)
  }
  else op_cost_values <- NULL
  obj_values <- vector(mode = "numeric", length = iter)
  if (progress) 
    pb <- txtProgressBar(min = 1, max = iter, style = 3)
  if (any(duplicated(i_sampled))) 
    browser()
  if (any(i_sampled %in% i_unsampled)) 
    browser()
  for (i in 1:iter) {
    previous <- list()
    previous$obj <- obj
    previous$i_sampled <- i_sampled
    previous$i_unsampled <- i_unsampled
    previous$delta_obj_continuous <- delta_obj_continuous
    if (cost_mode) 
      previous$op_cost <- op_cost
    if (runif(1) < 0.5) {
      idx_removed <- sample(1:length(setdiff(i_sampled, 
                                             include)), size = 1, replace = FALSE)
      spl_removed <- setdiff(i_sampled, include)[idx_removed]
      idx_added <- sample(1:length(i_unsampled), size = 1, 
                          replace = FALSE)
      i_sampled <- setdiff(i_sampled, include)[-idx_removed]
      i_sampled <- c(i_sampled, i_unsampled[idx_added], 
                     include)
      i_unsampled <- i_unsampled[-idx_added]
      i_unsampled <- c(i_unsampled, spl_removed)
      if (any(duplicated(i_sampled))) 
        browser()
      if (any(i_sampled %in% i_unsampled)) 
        browser()
      data_continuous_sampled <- data_continuous[i_sampled, 
                                                 , drop = FALSE]
      if (n_factor > 0) 
        data_factor_sampled <- data_factor[i_sampled, 
                                           , drop = FALSE]
    }
    else {
      worse <- max(delta_obj_continuous[!i_sampled %in% 
                                          include])
      i_worse <- which(delta_obj_continuous[!i_sampled %in% 
                                              include] == worse)
      if (length(i_worse) > 1) 
        i_worse <- sample(i_worse, size = 1)
      spl_removed <- setdiff(i_sampled, include)[i_worse]
      idx_added <- sample(1:n_remainings, size = 1, replace = FALSE)
      i_sampled <- setdiff(i_sampled, include)[-i_worse]
      i_sampled <- c(i_sampled, i_unsampled[idx_added], 
                     include)
      i_unsampled <- i_unsampled[-idx_added]
      i_unsampled <- c(i_unsampled, spl_removed)
      if (any(duplicated(i_sampled))) 
        browser()
      if (any(i_sampled %in% i_unsampled)) 
        browser()
      data_continuous_sampled <- data_continuous[i_sampled, 
                                                 , drop = FALSE]
      if (n_factor > 0) 
        data_factor_sampled <- data_factor[i_sampled, 
                                           , drop = FALSE]
    }
    res <- .lhs_obj(size = size, data_continuous_sampled = data_continuous_sampled, 
                    data_factor_sampled = data_factor_sampled, continuous_strata = continuous_strata, 
                    cor_mat = cor_mat, factor_obj = factor_obj, weights = weights, 
                    eta = eta)
    obj <- res$obj
    delta_obj_continuous <- res$delta_obj_continuous
    delta_obj <- obj - previous$obj
    metropolis <- exp(-1 * delta_obj/temp)
    if (cost_mode) {
      op_cost <- sum(cost[i_sampled, ])
      delta_cost <- op_cost - previous$op_cost
      if (track_mode) 
        metropolis_cost <- Inf
      else metropolis_cost <- exp(-1 * delta_cost/temp)
    }
    else metropolis_cost <- Inf
    if (obj <= obj.limit) {
      warning("\nThe objective function has reached its minimum value, as specified by the obj.limit option.")
      if (progress) {
        setTxtProgressBar(pb, i)
        close(pb)
      }
      obj_values[i] <- obj
      if (cost_mode) 
        op_cost_values[i] <- op_cost
      break
    }
    if (delta_obj > 0 & runif(1) >= metropolis | runif(1) >= 
        metropolis_cost) {
      i_sampled <- previous$i_sampled
      i_unsampled <- previous$i_unsampled
      data_continuous_sampled <- data_continuous[i_sampled, 
                                                 , drop = FALSE]
      if (n_factor > 0) 
        data_factor_sampled <- data_factor[i_sampled, 
                                           , drop = FALSE]
      obj <- previous$obj
      delta_obj_continuous <- previous$delta_obj_continuous
      if (cost_mode) 
        op_cost <- previous$op_cost
    }
    obj_values[i] <- obj
    if (cost_mode) 
      op_cost_values[i] <- op_cost
    if ((i%%length.cycle) == 0) 
      temp <- temp * tdecrease
    if (progress) 
      setTxtProgressBar(pb, i)
  }
  if (progress) 
    close(pb)
  if (n_factor > 0) {
    sampled_data <- data.frame(data_continuous_sampled, 
                               data_factor_sampled, stringsAsFactors = TRUE)
    sampled_data <- sampled_data[, names(x)]
  }
  else sampled_data <- data_continuous_sampled
  if (simple) 
    res <- i_sampled
  else {
    res <- list(initial_object = x, index_samples = i_sampled, 
                sampled_data = sampled_data, obj = obj_values, cost = op_cost_values)
    class(res) = c("cLHS_result", "list")
  }
  res
}