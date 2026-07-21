


#=================Functions used in model_food_webs.R=========================
#---Basic----
get_Re_lmax <- function(J) {
  max(Re(eigen(J, only.values = TRUE)$values))
}

make_z <- function(x) {
  as.numeric(scale(x))
}


#----Rebuild food-web objects----
rebuild_web_objects <- function(df) {
  n <- unique(df$richness_eq)
  
  C <- matrix(0, n, n)
  e_link <- matrix(0, n, n)
  J_emp <- matrix(0, n, n)
  B <- numeric(n)
  
  link_idx <- cbind(df$prey, df$predator)
  
  C[link_idx] <- df$cij
  e_link[link_idx] <- df$e
  
  B[df$prey] <- df$B_prey
  B[df$predator] <- df$B_predator
  
  J_emp[link_idx] <- df$J_topdown
  J_emp[cbind(df$predator, df$prey)] <- df$J_bottomup
  
  out <- list(
    C = C,
    B = B,
    e_link = e_link,
    J_emp = J_emp
  )
  
  J_ode_diag <- numeric(n)
  
  J_ode_diag[df$prey] <- df$J_ode_diag_prey
  J_ode_diag[df$predator] <- df$J_ode_diag_predator
  
  out$J_ode_diag <- J_ode_diag
  
  out
}


#----Build Jacobian----
build_jacobian <- function(C, B, e_link, J_diag = NULL) {
  n <- length(B)
  J <- matrix(0, n, n)
  links <- which(C > 0, arr.ind = TRUE)
  
  for (k in seq_len(nrow(links))) {
    i <- links[k, 1]
    j <- links[k, 2]
    
    J[i, j] <- -C[i, j] * B[i]
    J[j, i] <- e_link[i, j] * C[i, j] * B[j]
  }
  
  if (!is.null(J_diag)) {
    diag(J) <- J_diag
  }
  
  J
}


#---Interaction-derived diagonal----
get_interaction_diag_ode <- function(C, B, e_link) {
  J_diag <- numeric(length(B))
  links <- which(C > 0, arr.ind = TRUE)
  
  for (k in seq_len(nrow(links))) {
    i <- links[k, 1]
    j <- links[k, 2]
    
    J_diag[i] <- J_diag[i] -
      C[i, j] * B[j]
    
    J_diag[j] <- J_diag[j] +
      e_link[i, j] * C[i, j] * B[i]
  }
  
  J_diag
}

#----Randomise Jacobian---
randomise_jacobian <- function(
    C,
    B,
    e_link,
    scenario,
    background_diag = NULL
) {
  n <- length(B)
  links <- which(C > 0, arr.ind = TRUE)
  
  c_pool <- C[C > 0]
  B_pool <- B[B > 0]
  e_pool <- e_link[e_link > 0]
  
  J <- matrix(0, n, n)
  include_diag <- !is.null(background_diag)
  
  if (include_diag) {
    diag(J) <- background_diag
  }
  
  # Shared species biomass
  if (scenario %in% c("full_dependency", "B_dependency")) {
    B_species <- sample(B_pool, n, replace = TRUE)
  }
  
  # Shared consumer efficiency
  if (scenario %in% c("full_dependency", "e_dependency")) {
    predators <- unique(links[, 2])
    e_predator <- numeric(n)
    
    e_predator[predators] <- sample(
      e_pool,
      length(predators),
      replace = TRUE
    )
  }
  
  for (k in seq_len(nrow(links))) {
    i <- links[k, 1]
    j <- links[k, 2]
    
    # Interaction coefficients
    if (scenario %in% c("full_dependency", "c_dependency")) {
      c_ij <- sample(c_pool, 1, replace = TRUE)
      
      c_top <- c_ij
      c_bottom <- c_ij
      c_prey_diag <- c_ij
      c_pred_diag <- c_ij
    } else {
      c_top <- sample(c_pool, 1, replace = TRUE)
      c_bottom <- sample(c_pool, 1, replace = TRUE)
      c_prey_diag <- sample(c_pool, 1, replace = TRUE)
      c_pred_diag <- sample(c_pool, 1, replace = TRUE)
    }
    
    # Biomass terms
    if (scenario %in% c("full_dependency", "B_dependency")) {
      B_i <- B_species[i]
      B_j <- B_species[j]
    } else {
      B_i <- sample(B_pool, 1, replace = TRUE)
      B_j <- sample(B_pool, 1, replace = TRUE)
    }
    
    # Efficiency terms
    if (scenario %in% c("full_dependency", "e_dependency")) {
      e_bottom <- e_pred_diag <- e_predator[j]
    } else {
      e_bottom <- sample(e_pool, 1, replace = TRUE)
      e_pred_diag <- sample(e_pool, 1, replace = TRUE)
    }
    
    # Off-diagonal terms
    J[i, j] <- -c_top * B_i
    J[j, i] <- e_bottom * c_bottom * B_j
    
    # Interaction-derived diagonal terms
    if (include_diag) {
      if (scenario %in% c("full_dependency", "B_dependency")) {
        B_i_diag <- B_i
        B_j_diag <- B_j
      } else {
        B_i_diag <- sample(B_pool, 1, replace = TRUE)
        B_j_diag <- sample(B_pool, 1, replace = TRUE)
      }
      
      J[i, i] <- J[i, i] -
        c_prey_diag * B_j_diag
      
      J[j, j] <- J[j, j] +
        e_pred_diag * c_pred_diag * B_i_diag
    }
  }
  
  J
}
