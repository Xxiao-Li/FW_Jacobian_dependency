


#=================Functions used in 3_model_food_webs.R=========================
get_Re_lmax <- function(J) max(Re(eigen(J, only.values = TRUE)$values))

rebuild_web_objects <- function(df) {
  n <- unique(df$richness_eq)
  
  C <- e_link <- J_emp <- matrix(0, n, n)
  B <- rep(NA_real_, n)
  
  C[cbind(df$prey, df$predator)] <- df$cij
  e_link[cbind(df$prey, df$predator)] <- df$e
  B[df$prey] <- df$B_prey
  B[df$predator] <- df$B_predator
  
  J_emp[cbind(df$prey, df$predator)] <- df$J_topdown
  J_emp[cbind(df$predator, df$prey)] <- df$J_bottomup
  diag(J_emp) <- 0
  
  out <- list(C = C, B = B, e_link = e_link, J_emp = J_emp)
  
  if (all(c("J_ode_diag_prey", "J_ode_diag_predator") %in% names(df))) {
    J_ode_diag <- rep(NA_real_, n)
    J_ode_diag[df$prey] <- df$J_ode_diag_prey
    J_ode_diag[df$predator] <- df$J_ode_diag_predator
    
    if (any(is.na(J_ode_diag))) {
      warning("Some species are missing J_ode diagonal values.")
    }
    
    out$J_ode_diag <- J_ode_diag
  }
  
  out
}

build_jacobian <- function(C, B, e_link, J_diag = NULL) {
  n <- length(B)
  J <- matrix(0, n, n)
  links <- which(C > 0, arr.ind = TRUE)
  
  for (k in seq_len(nrow(links))) {
    i <- links[k, 1]
    j <- links[k, 2]
    J[i, j] <- -C[i, j] * B[i]
    J[j, i] <-  e_link[i, j] * C[i, j] * B[j]
  }
  
  if (is.null(J_diag)) {
    diag(J) <- 0
  } else {
    stopifnot(length(J_diag) == n)
    diag(J) <- J_diag
  }
  
  J
}

get_parameter_pools <- function(C, B, e_link) {
  list(
    c_pool = as.numeric(C[C > 0]),
    B_pool = as.numeric(B[B > 0]),
    e_pool = as.numeric(e_link[e_link > 0])
  )
}

get_links <- function(C) {
  links <- which(C > 0, arr.ind = TRUE)
  tibble(prey = links[, 1], predator = links[, 2])
}

randomise_params <- function(C, B, e_link, scenario) {
  n <- length(B)
  links <- get_links(C)
  pools <- get_parameter_pools(C, B, e_link)
  
  c_pool <- pools$c_pool
  B_pool <- pools$B_pool
  e_pool <- pools$e_pool
  
  if (scenario == "full_dependency") {
    C_new <- e_new <- matrix(0, n, n)
    
    for (k in seq_len(nrow(links))) {
      i <- links$prey[k]
      j <- links$predator[k]
      C_new[i, j] <- sample(c_pool, 1, replace = TRUE)
      e_new[i, j] <- sample(e_pool, 1, replace = TRUE)
    }
    
    return(list(C = C_new, B = sample(B_pool, n, replace = TRUE), e_link = e_new))
  }
  
  J <- matrix(0, n, n)
  
  if (scenario == "c_dependency") {
    for (k in seq_len(nrow(links))) {
      i <- links$prey[k]
      j <- links$predator[k]
      c_ij <- sample(c_pool, 1, replace = TRUE)
      J[i, j] <- -c_ij * sample(B_pool, 1, replace = TRUE)
      J[j, i] <-  sample(e_pool, 1, replace = TRUE) * c_ij * sample(B_pool, 1, replace = TRUE)
    }
  }
  
  if (scenario == "B_dependency") {
    B_new <- sample(B_pool, n, replace = TRUE)
    
    for (k in seq_len(nrow(links))) {
      i <- links$prey[k]
      j <- links$predator[k]
      J[i, j] <- -sample(c_pool, 1, replace = TRUE) * B_new[i]
      J[j, i] <-  sample(e_pool, 1, replace = TRUE) * sample(c_pool, 1, replace = TRUE) * B_new[j]
    }
  }
  
  if (scenario == "e_dependency") {
    predators <- sort(unique(links$predator))
    e_by_predator <- rep(0, n)
    e_by_predator[predators] <- sample(e_pool, length(predators), replace = TRUE)
    
    for (k in seq_len(nrow(links))) {
      i <- links$prey[k]
      j <- links$predator[k]
      J[i, j] <- -sample(c_pool, 1, replace = TRUE) * sample(B_pool, 1, replace = TRUE)
      J[j, i] <-  e_by_predator[j] * sample(c_pool, 1, replace = TRUE) * sample(B_pool, 1, replace = TRUE)
    }
  }
  
  if (scenario == "no_dependency") {
    for (k in seq_len(nrow(links))) {
      i <- links$prey[k]
      j <- links$predator[k]
      J[i, j] <- -sample(c_pool, 1, replace = TRUE) * sample(B_pool, 1, replace = TRUE)
      J[j, i] <-  sample(e_pool, 1, replace = TRUE) * sample(c_pool, 1, replace = TRUE) * sample(B_pool, 1, replace = TRUE)
    }
  }
  
  diag(J) <- 0
  list(J = J)
}

compute_randomised_Re_lmax <- function(C, B, e_link, scenario) {
  rand <- randomise_params(C, B, e_link, scenario)
  
  if (!is.null(rand$J)) {
    diag(rand$J) <- 0
    return(get_Re_lmax(rand$J))
  }
  
  get_Re_lmax(build_jacobian(rand$C, rand$B, rand$e_link))
}

get_interaction_diag_ode <- function(C, B, e_link) {
  n <- length(B)
  diag_vec <- rep(0, n)
  links <- which(C > 0, arr.ind = TRUE)
  
  for (k in seq_len(nrow(links))) {
    i <- links[k, 1]  # resource/prey
    j <- links[k, 2]  # consumer/predator
    
    c_ij <- C[i, j]
    e_ij <- e_link[i, j]
    
    diag_vec[i] <- diag_vec[i] - c_ij * B[j]
    diag_vec[j] <- diag_vec[j] + e_ij * c_ij * B[i]
  }
  
  diag_vec
}

build_jacobian_ode_interaction <- function(C, B, e_link, background_diag) {
  n <- length(B)
  J <- matrix(0, n, n)
  diag(J) <- background_diag
  
  links <- which(C > 0, arr.ind = TRUE)
  
  for (k in seq_len(nrow(links))) {
    i <- links[k, 1]  # resource/prey
    j <- links[k, 2]  # consumer/predator
    
    c_ij <- C[i, j]
    e_ij <- e_link[i, j]
    
    # Off-diagonal terms
    J[i, j] <- J[i, j] - c_ij * B[i]
    J[j, i] <- J[j, i] + e_ij * c_ij * B[j]
    
    # Diagonal terms from the same interaction
    J[i, i] <- J[i, i] - c_ij * B[j]
    J[j, j] <- J[j, j] + e_ij * c_ij * B[i]
  }
  
  J
}

randomise_jacobian_ode <- function(C, B, e_link, scenario, background_diag) {
  n <- length(B)
  links <- get_links(C)
  pools <- get_parameter_pools(C, B, e_link)
  
  c_pool <- pools$c_pool
  B_pool <- pools$B_pool
  e_pool <- pools$e_pool
  
  J <- matrix(0, n, n)
  diag(J) <- background_diag
  
  if (scenario == "full_dependency") {
    C_new <- e_new <- matrix(0, n, n)
    B_new <- sample(B_pool, n, replace = TRUE)
    
    for (k in seq_len(nrow(links))) {
      i <- links$prey[k]
      j <- links$predator[k]
      C_new[i, j] <- sample(c_pool, 1, replace = TRUE)
      e_new[i, j] <- sample(e_pool, 1, replace = TRUE)
    }
    
    return(build_jacobian_ode_interaction(C_new, B_new, e_new, background_diag))
  }
  
  if (scenario == "c_dependency") {
    for (k in seq_len(nrow(links))) {
      i <- links$prey[k]
      j <- links$predator[k]
      c_ij <- sample(c_pool, 1, replace = TRUE)
      
      B_i_off <- sample(B_pool, 1, replace = TRUE)
      B_j_off <- sample(B_pool, 1, replace = TRUE)
      B_i_diag <- sample(B_pool, 1, replace = TRUE)
      B_j_diag <- sample(B_pool, 1, replace = TRUE)
      e_off <- sample(e_pool, 1, replace = TRUE)
      e_diag <- sample(e_pool, 1, replace = TRUE)
      
      J[i, j] <- -c_ij * B_i_off
      J[j, i] <-  e_off * c_ij * B_j_off
      J[i, i] <- J[i, i] - c_ij * B_j_diag
      J[j, j] <- J[j, j] + e_diag * c_ij * B_i_diag
    }
  }
  
  if (scenario == "B_dependency") {
    B_new <- sample(B_pool, n, replace = TRUE)
    
    for (k in seq_len(nrow(links))) {
      i <- links$prey[k]
      j <- links$predator[k]
      
      c1 <- sample(c_pool, 1, replace = TRUE)
      c2 <- sample(c_pool, 1, replace = TRUE)
      c3 <- sample(c_pool, 1, replace = TRUE)
      c4 <- sample(c_pool, 1, replace = TRUE)
      e1 <- sample(e_pool, 1, replace = TRUE)
      e2 <- sample(e_pool, 1, replace = TRUE)
      
      J[i, j] <- -c1 * B_new[i]
      J[j, i] <-  e1 * c2 * B_new[j]
      J[i, i] <- J[i, i] - c3 * B_new[j]
      J[j, j] <- J[j, j] + e2 * c4 * B_new[i]
    }
  }
  
  if (scenario == "e_dependency") {
    predators <- sort(unique(links$predator))
    e_by_predator <- rep(0, n)
    e_by_predator[predators] <- sample(e_pool, length(predators), replace = TRUE)
    
    for (k in seq_len(nrow(links))) {
      i <- links$prey[k]
      j <- links$predator[k]
      
      J[i, j] <- -sample(c_pool, 1, replace = TRUE) * sample(B_pool, 1, replace = TRUE)
      J[j, i] <-  e_by_predator[j] * sample(c_pool, 1, replace = TRUE) * sample(B_pool, 1, replace = TRUE)
      J[i, i] <- J[i, i] - sample(c_pool, 1, replace = TRUE) * sample(B_pool, 1, replace = TRUE)
      J[j, j] <- J[j, j] + e_by_predator[j] * sample(c_pool, 1, replace = TRUE) * sample(B_pool, 1, replace = TRUE)
    }
  }
  
  if (scenario == "no_dependency") {
    for (k in seq_len(nrow(links))) {
      i <- links$prey[k]
      j <- links$predator[k]
      
      J[i, j] <- -sample(c_pool, 1, replace = TRUE) * sample(B_pool, 1, replace = TRUE)
      J[j, i] <-  sample(e_pool, 1, replace = TRUE) * sample(c_pool, 1, replace = TRUE) * sample(B_pool, 1, replace = TRUE)
      J[i, i] <- J[i, i] - sample(c_pool, 1, replace = TRUE) * sample(B_pool, 1, replace = TRUE)
      J[j, j] <- J[j, j] + sample(e_pool, 1, replace = TRUE) * sample(c_pool, 1, replace = TRUE) * sample(B_pool, 1, replace = TRUE)
    }
  }
  
  J
}

compute_randomised_Re_lmax_ode <- function(C, B, e_link, scenario, background_diag) {
  J <- randomise_jacobian_ode(C, B, e_link, scenario, background_diag)
  get_Re_lmax(J)
}




safe_sd <- function(x) {
  s <- sd(x, na.rm = TRUE)
  ifelse(is.na(s) || s == 0, NA_real_, s)
}

safe_se <- function(x) {
  n <- sum(!is.na(x))
  s <- sd(x, na.rm = TRUE)
  ifelse(n <= 1 || is.na(s), NA_real_, s / sqrt(n))
}

make_z <- function(x) as.numeric(scale(x))