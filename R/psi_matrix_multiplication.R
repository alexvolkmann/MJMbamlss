psi_mat_crossprod <- function(Psi, R = NULL, y = NULL) {

  X <- Psi$X[Psi$cp_info[[1]], ]
  ni_obs <- Psi$cp_info[[2]]

  if (!is.null(R)) {

    # Hesse
    R <- R[Psi$cp_info[[1]]]
    .Call("psi_mat_multiplication", X, ni_obs, R)

  } else {

    # Score longitudinal part
    y <- y[Psi$cp_info[[1]]]
    .Call("psi_vec_multiplication", X, ni_obs, y)

  }
  # Score survival part does not seem to profit from C

}
