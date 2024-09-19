# User-defined function for mml project


# Generate simulation data
generate_mml <- function(J_k_vec, theta_tilde, X, R) {
  
  #General setting
  K <- length(J_k_vec); n <- nrow(X); p <- ncol(X);
  out_cat_data <- data.frame(out = as.character(1:K), num_cat = J_k_vec)
  combi <- combn(1:K, 2)
  combi_set <- paste("(", combi[1, ], ",", combi[2, ], ")")
  combi_list <- combi_set %>% str_split(",") %>% `names<-`(value = combi_set)
  
  #Calculate joint probabiliy for every possible cases
  ##Calculate UQ
  poss_UQ_mat_list <- Map(list, J_k_vec) %>% lapply(FUN = function(x) `names<-`(x, "J_k"))
  
  for(k in 1:K) {
    tmp_list <- vector("list", K)
    for(idx in 1:K) {
      if(k == idx) {
        tmp_list[[idx]] <- diag(J_k_vec[idx])
      } else {
        tmp_list[[idx]] <- matrix(rep(1, J_k_vec[idx]), ncol = 1)
      }
    }
    tmp_U_k <- reduce(tmp_list, .f = kronecker)
    tmp_Q_k <- rbind(matrix(0, nrow = 1, ncol = J_k_vec[k] - 1), diag(J_k_vec[k] - 1))
    poss_UQ_mat_list[[k]] <- tmp_U_k %*% tmp_Q_k
  }
  ##Calculate W
  poss_W_mat_list <- combi_set %>% str_split(",") %>% `names<-`(value = combi_set) %>%
    lapply(FUN = function(x) str_extract_all(x, one_or_more(DGT), simplify = TRUE) %>% as.numeric)
  for(i in 1:length(poss_W_mat_list)) {
    tmp_combi <- poss_W_mat_list[[i]]
    tmp_list <- vector("list", K) 
    for(idx in 1:K) {
      if (idx %in% tmp_combi) {
        tmp_list[[idx]] <- diag(J_k_vec[idx])
      } else {
        tmp_list[[idx]] <- matrix(rep(1, J_k_vec[idx]), ncol = 1)
      }
    }
    poss_W_mat_list[[i]] <- reduce(tmp_list, .f = kronecker)
  }
  poss_UQ_mat <- reduce(poss_UQ_mat_list, cbind)
  poss_W_mat <- reduce(poss_W_mat_list, cbind)
  
  ##Calculate Z
  J_0 <- sum(J_k_vec - 1); S_0 <- combi_data %>% filter(!str_detect(category, pattern = "1")) %>% nrow();
  X_list <- lapply(seq_len(nrow(X)), FUN = function(i) X[i, ])
  
  alpha_tilde <- theta_tilde[1:J_0]
  poss_alpha <- poss_UQ_mat %*% alpha_tilde
  psi_tilde <- theta_tilde[(J_0+1):(J_0+S_0)]
  poss_psi <- poss_W_mat %*% R %*% psi_tilde
  beta_tilde <- theta_tilde[(J_0+S_0+1):(J_0+S_0+ J_0*p)]
  poss_beta_list <- lapply(X_list, FUN = function(x) kronecker(poss_UQ_mat, t(x)) %*% beta_tilde)
  delta_tilde <- theta_tilde %>% tail(S_0*p)
  poss_delta_list <- lapply(X_list, FUN = function(x) kronecker(poss_W_mat %*% R, t(x)) %*% delta_tilde)
  
  Z_list <- mapply(FUN = function(beta, delta) beta + delta, beta = poss_beta_list, delta = poss_delta_list, SIMPLIFY = FALSE) %>% 
    lapply(FUN = function(x) (x + poss_alpha + poss_psi)) #%>%
  #mapply(FUN = function(Z, eps) Z + eps, Z = ., eps = eps, SIMPLIFY = FALSE)
  
  prob_list <- Z_list %>% lapply(FUN = function(x) exp(x)/sum(exp(x)))
  assign_list <- prob_list %>% lapply(FUN = function(x) sample(1:prod(J_k_vec), size = 1, replace = FALSE, prob = as.vector(x)))
  
  poss <- sapply(J_k_vec, FUN = function(x) 1:x) %>% expand.grid %>% arrange(across(everything()))
  Y <- assign_list %>% lapply(FUN = function(x) poss[x, ]) %>% bind_rows() %>%
    mutate_all(as.factor) %>%
    `names<-`(paste0("y", 1:K))
  
  out <- list(Y = Y, joint_Z = Z_list, assign_list = assign_list, poss_W_mat = poss_W_mat, poss_UQ_mat = poss_UQ_mat)
  return(out);
}


# Calculate the design matrix
calc_design_mat <- function(Data, U.k_ref, W.k_ref) {
  K <- ncol(Data);
  Combi <- Data %>% as.data.frame() %>% unite(col = "combi", sep = ",") %>% mutate(combi = paste("(", combi, ")"))
  U.k_list <- vector(mode = "list", length = K)
  W.k_list <- vector(mode = "list", length = K)
  for(k in 1:K) {
    sub_U.k <- Combi %>% left_join(poss_U_mat_list[[k]], by = "combi") %>% dplyr::select(-combi)
    sub_W.k <- Combi %>% left_join(poss_W_mat_list[[k]], by = "combi") %>% dplyr::select(-combi)
    
    U.k_list[[k]] <- sub_U.k %>% as.matrix()
    W.k_list[[k]] <- sub_W.k %>% as.matrix()
  }
  out <- list(U.k = U.k_list, W.k = W.k_list)
  return(out)
}


# Estimate the coefficient vector using mm algorithm (with group lasso)
est_mml_MM <- function(max_iter, Z_list_0, Z_tilde_list_0, X_tilde_stack, Y, lambda, p, adaptive = FALSE, weight) {
  # Setting
  n <- nrow(Y); K <- ncol(Y);
  J_k_vec <- Y %>% apply(2, unique) %>% lapply(FUN = function(x) x %>% length) %>% unlist
  J_0 <- sum(J_k_vec - 1); S_0 <- (ncol(X_tilde_stack) - J_0*(p+1))/(p+1);
  
  colidx_list <- vector(mode = "list", length = K)
  colidx_sum <- 1
  for(k in 1:K) {
    colidx_list[[k]] <- (colidx_sum):(colidx_sum + (J_k_vec[k] - 2))
    colidx_sum <- colidx_sum + J_k_vec[k] - 1
  }
  o.ik_list <- lapply(seq_len(ncol(Y)), FUN = function(i) Y[, i] %>% as.numeric %>% data.frame("cat" = .)) 
  y.ik_list <- vector(mode = "list", length = K)
  for(k in 1:K) {
    y.ik <-  matrix(0, nrow = n, ncol = J_k_vec[k])
    for(i in 1:n) {
      y.ik[i, o.ik_list[[k]][i,]] <- 1 
    }
    y.ik_list[[k]] <- y.ik
  }
  obs_idx_list <- lapply(seq(1, K), FUN = function(k) Y[, k])
  
  # Initial setting for iteration
  coef_mat <- vector(mode = "numeric", length = ncol(X_tilde_stack))
  err_set <- c()
  cond_likel_set <- c()
  group <- c(rep(1, J_0 + S_0), rep(1:p, J_0 + S_0) + 1)
  
  # Initial value
  Z_list_new <- Z_list_0; Z_tilde_list_new <- Z_tilde_list_0;
  
  for(iter in 1:max_iter) {
    # Replace updated values
    Z_list_prev <- Z_list_new
    Z_tilde_list_prev <- Z_tilde_list_new
    
    # Maximization
    Z_tilde_prev_stack <- Z_tilde_list_prev %>% unlist() 
    data_tmp <- cbind(Z_tilde_prev_stack, X_tilde_stack) %>%
      `colnames<-`(value = c("Z_tilde", paste0("param", 1:ncol(X_tilde_stack)))) %>% 
      data.frame()
    if(lambda == 0) {
      fit <- lm(Z_tilde ~ . - 1, data = data_tmp)
      est_coef <- coef(fit)
    } else {
      if(adaptive == TRUE & missing(weight)) {
        stop("weight argument is necessary for adaptive group lasso")
      }
      if(missing(weight)) {
        weight <- table(group) %>% as.vector() %>% sqrt
        weight[1] <- 0
      }
      # Arrange group argument to make it consecutive
      group_org <- sort(group)
      X_unorg <- data_tmp[, -1] %>% as.matrix()
      X_org <- c()
      for(g in 1:max(group)) {
        X_org <- cbind(X_org, X_unorg[, (group == g)])
      }
      fit <- gglasso(x = X_org, y = data_tmp$Z_tilde, group = group_org, pf = weight, loss = "ls", 
                     lambda = lambda, intercept = FALSE)
      est_coef_prev <- fit$beta[, 1]
      est_coef <- c()
      for(idx in 1:length(group)) {
        est_coef <- c(est_coef, est_coef_prev[order(group) == idx])
      }
    }
    coef_mat <- rbind(coef_mat, est_coef)
    
    #Update Z
    Z_mat_new <- X_tilde_stack %*% t(tail(coef_mat, 1)) %>% matrix(byrow = FALSE, nrow = n)
    Z_list_new <- lapply(colidx_list, FUN = function(idx) Z_mat_new[, idx] %>% matrix(byrow = FALSE, nrow = n))
    
    #Update Z tilde
    denom_Z_new_list <- lapply(Z_list_new, FUN = function(x) apply(x, MARGIN = 1, FUN = function(y) sum(exp(y), 1)))
    
    for(k in 1:K) {
      if(J_k_vec[k] == 2) {
        L_k <- -1/4
      } else if(J_k_vec[k] == 3) {
        L_k <- -1/2
      } else {L_k <- -1}
      for(j in 1:(J_k_vec[k] - 1)) {
        for(i in 1:n) {
          Z_tilde_list_new[[k]][i, j] <- Z_list_new[[k]][i, j] - 
            (1/L_k)*(y.ik_list[[k]][i,(j+1)] - 
                       exp(Z_list_new[[k]][i, j])/denom_Z_new_list[[k]][i])
        }
      }
    }
    
    cond_likel <- calc_cond_likel(Z_list = Z_list_new, obs_idx_list = obs_idx_list)
    cond_likel_set <- c(cond_likel_set, cond_likel)
    err <- tail(coef_mat, 2) %>% apply(MARGIN = 2, diff) %>% abs %>% max
    err_set <- c(err_set, err)
    
    if(abs(err) < 0.1^3) break;
  }
  
  # Output
  out <- list(coef = coef_mat, cond_likel = cond_likel_set, err = err_set, lambda = lambda)
  return(out)
}

# Calculate conditional likelihood
calc_cond_likel <- function(Z_list, obs_idx_list) {
  n <- Z_list[[1]] %>% nrow()
  denom_list <- lapply(Z_list, FUN = function(x) apply(x, MARGIN = 1, FUN = function(y) sum(exp(y), 1)))
  Z_list_complete <- lapply(Z_list, FUN = function(x) cbind(0, x))
  Z_list_obs <- mapply(FUN = function(Z, idx) Z[cbind(seq_len(n), idx)], Z = Z_list_complete, idx = obs_idx_list, SIMPLIFY = FALSE) %>% unlist
  log_sum_exp <- unlist(denom_list) %>% log
  cond_likel <- sum(Z_list_obs - log_sum_exp)
  
  return(cond_likel)
}

# Calculate degree of freedom in group lasso
cal_df_gglasso <- function(beta, beta_ls, group) {
  p_j <- table(group) %>% as.vector()
  
  beta_l2 <- c()
  beta_ls_l2 <- c()
  for(g in 1:max(group)) {
    beta_j <- beta[group == g]
    beta_j_l2 <- beta_j^2 %>% sum %>% sqrt
    beta_l2 <- c(beta_l2, beta_j_l2)
    
    beta_ls_j <- beta_ls[group == g]
    beta_ls_j_l2 <- beta_ls_j^2 %>% sum %>% sqrt
    beta_ls_l2 <- c(beta_ls_l2, beta_ls_j_l2)
  }
  
  out <- sum(beta_l2 > 0) + sum(beta_l2/beta_ls_l2*(p_j - 1))
  
  return(out)
}

# Estimate the coefficient vector using mm algorithm (with gBridge)
est_mml_gBridge <- function(max_iter, Z_list_0, Z_tilde_list_0, X_tilde_stack, Y, lamb_seq, p) {
  #Setting
  n <- nrow(Y); K <- ncol(Y);
  J_k_vec <- Y %>% apply(2, unique) %>% lapply(FUN = function(x) x %>% length) %>% unlist
  J_0 <- sum(J_k_vec - 1); S_0 <- (ncol(X_tilde_stack) - J_0*(p+1))/(p+1);
  
  colidx_list <- vector(mode = "list", length = K)
  colidx_sum <- 1
  for(k in 1:K) {
    colidx_list[[k]] <- (colidx_sum):(colidx_sum + (J_k_vec[k] - 2))
    colidx_sum <- colidx_sum + J_k_vec[k] - 1
  }
  o.ik_list <- lapply(seq_len(ncol(Y)), FUN = function(i) Y[, i] %>% as.numeric %>% data.frame("cat" = .)) 
  y.ik_list <- vector(mode = "list", length = K)
  for(k in 1:K) {
    y.ik <-  matrix(0, nrow = n, ncol = J_k_vec[k])
    for(i in 1:n) {
      y.ik[i, o.ik_list[[k]][i,]] <- 1 
    }
    y.ik_list[[k]] <- y.ik
  }
  obs_idx_list <- lapply(seq(1, K), FUN = function(k) Y[, k])
  
  #Initial setting
  coef_list <- vector(mode = "list", length = length(lamb_seq))
  BIC_set <- vector(mode = "numeric", length = length(lamb_seq))
  err_list <- vector(mode = "list", length = length(lamb_seq))
  cond_likel_list <- vector(mode = "list", length = length(lamb_seq))
  group <- c(rep(1, J_0 + S_0), rep(1:p, J_0 + S_0) + 1)
  group_mult <- table(group) %>% as.vector()
  group_mult[1] <- 0
  
  for(lamb_idx in 0:length(lamb_seq)) {
    #Initial value
    Z_list_new <- Z_list_0; Z_tilde_list_new <- Z_tilde_list_0;
    coef_set <- vector(mode = "numeric", length = ncol(X_tilde_stack))
    cond_likel_set <- NULL
    err_set <- c()
    for(iter in 1:max_iter) {
      #replace updated values
      Z_list_prev <- Z_list_new
      Z_tilde_list_prev <- Z_tilde_list_new
      
      # Minimization
      Z_tilde_prev_stack <- Z_tilde_list_prev %>% unlist() 
      data_tmp <- cbind(Z_tilde_prev_stack, X_tilde_stack) %>%
        `colnames<-`(value = c("Z_tilde", paste0("param", 1:ncol(X_tilde_stack)))) %>% 
        data.frame()
      
      # Fit group lasso or linear regression
      if(lamb_idx == 0) {
        fit <- lm(Z_tilde ~ . - 1, data = data_tmp)
        est_coef <- coef(fit)
      } else{
        group_org <- sort(group)
        X_unorg <- data_tmp[, -1] %>% as.matrix()
        X_org <- c()
        for(g in 1:max(group)) {
          X_org <- cbind(X_org, X_unorg[, (group == g)])
        }
        centered_response <- data_tmp$Z_tilde %>% scale(center = TRUE, scale = FALSE)
        colwise_centered_X_org <- apply(X_org, MARGIN = 2, FUN = scale, center = TRUE, scale = FALSE)
        fit <- gBridge(X = colwise_centered_X_org, y = centered_response, lambda = lamb_seq[lamb_idx], 
                       group = group_org, group.multiplier = group_mult, family = "gaussian") 
        est_coef_prev <- fit$beta[-1, 1]
        est_coef <- c()
        for(idx in 1:length(group)) {
          est_coef <- c(est_coef, est_coef_prev[order(group) == idx])
        }
      }
      coef_set <- rbind(coef_set, est_coef)
      
      #Update Z
      Z_mat_new <- X_tilde_stack %*% t(tail(coef_set, 1)) %>% matrix(byrow = FALSE, nrow = n)
      Z_list_new <- lapply(colidx_list, FUN = function(idx) Z_mat_new[, idx] %>% matrix(byrow = FALSE, nrow = n))
      
      #Update Z tilde
      denom_Z_new_list <- lapply(Z_list_new, FUN = function(x) apply(x, MARGIN = 1, FUN = function(y) sum(exp(y), 1)))
      
      for(k in 1:K) {
        if(J_k_vec[k] == 2) {
          L_k <- -1/4
        } else if(J_k_vec[k] == 3) {
          L_k <- -1/2
        } else {L_k <- -1}
        for(j in 1:(J_k_vec[k] - 1)) {
          for(i in 1:n) {
            Z_tilde_list_new[[k]][i, j] <- Z_list_new[[k]][i, j] - 
              (1/L_k)*(y.ik_list[[k]][i,(j+1)] - 
                         exp(Z_list_new[[k]][i, j])/denom_Z_new_list[[k]][i])
          }
        }
      }
      
      cond_likel <- calc_cond_likel(Z_list = Z_list_new, obs_idx_list = obs_idx_list)
      cond_likel_set <- c(cond_likel_set, cond_likel)
      err <- tail(coef_set, 2) %>% apply(MARGIN = 2, diff) %>% abs %>% max
      err_set <- c(err_set, err)
      
      if(abs(err) < 0.1^3) break;
    }
    if(lamb_idx == 0) {
      coef_ls <- coef_set %>% `rownames<-`(paste0("iter", 0:iter))
      ls_cond_likel <- calc_cond_likel(Z_list = Z_list_new, obs_idx_list = obs_idx_list)
    } else {
      coef_list[[lamb_idx]] <- coef_set %>% `rownames<-`(paste0("iter", 0:iter))
      err_list[[lamb_idx]] <- err_set
      cond_likel_list[[lamb_idx]] <- cond_likel_set
    }
  }
  
  # Output
  out <- list(coef_ls = coef_ls, ls_cond_likel = ls_cond_likel, coef = coef_list, err = err_list, cond_likel = cond_likel_list)
  return(out)
}

