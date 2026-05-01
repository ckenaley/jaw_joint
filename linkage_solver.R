# ============================================================
# Symmetric linkage solver
#
# This file builds and solves the symmetric bass linkage model. It
# includes low-level vector geometry, model construction, global-pose
# handling, D/C/B/E path solving, diagnostics, and Plotly animation.
# Function comments use roxygen-style tags for package-like readability.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(plotly)
})

# ============================================================
# Small helpers
# ============================================================

#' Return a fallback value when the first value is NULL.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param x Numeric vector, coordinate, or value used by the helper.
#' @param y Fallback value or coordinate component.
#' @return Either x or y.
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

nrm <- function(x) sqrt(sum(x^2))

#' Normalize a numeric vector to unit length, with a guard against near-zero vectors.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param x Numeric vector, coordinate, or value used by the helper.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @param label Input argument used by this helper.
#' @return A normalized numeric vector.
#' @keywords internal

unit2 <- function(x, tol = 1e-12, label = "vector") {
  d <- nrm(x)
  if (!is.finite(d) || d < tol) {
    stop(sprintf(
      "Cannot normalize near-zero %s: [%g, %g, %g], norm=%g",
      label, x[1], x[2], x[3], d
    ))
  }
  x / d
}

#' Compute Euclidean distance between two 3D points.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param a First vector, point, or object used in the calculation.
#' @param b Second vector, point, or object used in the calculation.
#' @return A single numeric distance.
#' @keywords internal
dist3 <- function(a, b) nrm(a - b)

cross3 <- function(a, b) {
  c(
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1]
  )
}

#' Reflect a 3D point across the x-z plane by reversing its y coordinate.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param p Point to evaluate.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
reflect_xz <- function(p) c(p[1], -p[2], p[3])

angle_3pt <- function(A, B, C, deg = TRUE, tol = 1e-12) {
  u <- A - B
  v <- C - B
  nu <- nrm(u)
  nv <- nrm(v)
  if (nu < tol || nv < tol) return(NA_real_)
  cs <- sum(u * v) / (nu * nv)
  cs <- max(-1, min(1, cs))
  ang <- acos(cs)
  if (deg) ang * 180 / pi else ang
}

#' Test whether an object is a finite numeric 3D point.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param x Numeric vector, coordinate, or value used by the helper.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
is_xyz_point <- function(x) {
  is.numeric(x) && length(x) == 3 && all(is.finite(x))
}

#' Convert a named list of 3D points to a coordinate tibble.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param pts Named list or vector of point definitions.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
pts_to_coords <- function(pts) {
  nm_pts <- names(pts)[vapply(pts, is_xyz_point, logical(1))]
  tibble(
    name = nm_pts,
    x = map_dbl(nm_pts, ~ pts[[.x]][1]),
    y = map_dbl(nm_pts, ~ pts[[.x]][2]),
    z = map_dbl(nm_pts, ~ pts[[.x]][3])
  )
}

# ============================================================
# Rotation helpers
# ============================================================

#' Rotate a vector around an arbitrary 3D axis using Rodrigues rotation.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param v Input argument used by this helper.
#' @param axis Axis to align when recasting coordinates.
#' @param angle Rotation angle in radians.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
rotate_about_axis <- function(v, axis, angle, tol = 1e-12) {
  axis <- unit2(axis, tol, "rotation axis")
  v * cos(angle) +
    cross3(axis, v) * sin(angle) +
    axis * sum(axis * v) * (1 - cos(angle))
}

#' Rotate a point around a 3D line defined by two points.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param P 3D point.
#' @param P1 First point defining a line or rotation axis.
#' @param P2 Second point defining a line or rotation axis.
#' @param angle Rotation angle in radians.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
rotate_point_about_line <- function(P, P1, P2, angle, tol = 1e-12) {
  ax <- P2 - P1
  if (nrm(ax) < tol) stop("Rotation axis is degenerate.")
  P1 + rotate_about_axis(P - P1, ax, angle, tol = tol)
}

#' Rotate all named 3D points in a list around a line, optionally excluding points.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param pts Named list or vector of point definitions.
#' @param P1 First point defining a line or rotation axis.
#' @param P2 Second point defining a line or rotation axis.
#' @param angle Rotation angle in radians.
#' @param exclude Character vector of named points to leave unrotated.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
rotate_named_points_about_line <- function(pts,
                                           P1,
                                           P2,
                                           angle,
                                           exclude = character(),
                                           tol = 1e-12) {
  out <- pts
  nm_keep <- setdiff(names(out), exclude)
  
  for (nm in nm_keep) {
    if (is_xyz_point(out[[nm]])) {
      out[[nm]] <- rotate_point_about_line(out[[nm]], P1, P2, angle, tol = tol)
    }
  }
  
  out
}

#' Create a function that rotates one axis direction onto another.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param axis_from Original axis direction.
#' @param axis_to Target axis direction.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
make_axis_alignment_rotation <- function(axis_from, axis_to, tol = 1e-12) {
  axis_from <- unit2(axis_from, tol, "axis_from")
  axis_to   <- unit2(axis_to, tol, "axis_to")
  
  ax <- cross3(axis_from, axis_to)
  s  <- nrm(ax)
  c0 <- sum(axis_from * axis_to)
  
  if (s < tol) {
    if (c0 > 0) {
      return(function(v) v)
    } else {
      trial <- c(1, 0, 0)
      if (abs(sum(trial * axis_from)) > 0.9) trial <- c(0, 1, 0)
      ax2 <- unit2(cross3(axis_from, trial), tol, "180deg helper axis")
      return(function(v) rotate_about_axis(v, ax2, pi, tol))
    }
  }
  
  ax <- ax / s
  ang <- atan2(s, c0)
  function(v) rotate_about_axis(v, ax, ang, tol)
}

# ============================================================
# 2D circle-circle helper
# ============================================================

#' Find the intersection point(s) of two circles in 2D.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param c1 Center of the first circle.
#' @param r1 Radius of the first circle.
#' @param c2 Center of the second circle.
#' @param r2 Radius of the second circle.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A list with ok, pts, and reason fields.
#' @keywords internal
circle_circle_2d <- function(c1, r1, c2, r2, tol = 1e-10) {
  d <- sqrt(sum((c2 - c1)^2))
  
  if (d > r1 + r2 + tol) {
    return(list(ok = FALSE, pts = NULL, reason = "separate circles"))
  }
  if (d < abs(r1 - r2) - tol) {
    return(list(ok = FALSE, pts = NULL, reason = "one circle inside another"))
  }
  if (d < tol && abs(r1 - r2) < tol) {
    return(list(ok = FALSE, pts = NULL, reason = "coincident circles"))
  }
  
  a  <- (r1^2 - r2^2 + d^2) / (2 * d)
  h2 <- r1^2 - a^2
  
  if (h2 < -tol) {
    return(list(ok = FALSE, pts = NULL, reason = "no real intersection"))
  }
  
  h  <- sqrt(max(0, h2))
  ex <- (c2 - c1) / d
  ey <- c(-ex[2], ex[1])
  p0 <- c1 + a * ex
  
  pts <- if (h < tol) {
    matrix(p0, nrow = 1)
  } else {
    rbind(p0 + h * ey, p0 - h * ey)
  }
  
  list(ok = TRUE, pts = pts, reason = NULL)
}

# ============================================================
# Model builder
# Required names: A, A_ref, B, B_ref, C, C_ref, D, E
# ============================================================

#' Build the symmetric linkage model and fixed geometric constraints from named coordinates.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param coords Data frame or tibble containing named coordinates with columns name, x, y, and z.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A sym_linkage_model object containing starting points, fixed link lengths, and projection helpers.
#' @keywords internal
make_sym_model <- function(coords, tol = 1e-10) {
  dat <- as_tibble(coords)
  
  req_cols <- c("name", "x", "y", "z")
  if (!all(req_cols %in% names(dat))) {
    stop("coords must contain columns: ", paste(req_cols, collapse = ", "))
  }
  
  get_pt <- function(nm) {
    ii <- match(nm, dat$name)
    if (is.na(ii)) stop("Point not found: ", nm)
    as.numeric(dat[ii, c("x", "y", "z")])
  }
  
  pts0 <- setNames(
    lapply(seq_len(nrow(dat)), function(i) as.numeric(dat[i, c("x", "y", "z")])),
    dat$name
  )
  
  A      <- get_pt("A")
  A_ref  <- get_pt("A_ref")
  B0     <- get_pt("B")
  B0_ref <- get_pt("B_ref")
  C0     <- get_pt("C")
  D0     <- get_pt("D")
  E0     <- get_pt("E")
  
  e1 <- unit2(C0 - A, tol, "C0 - A")
  nC <- unit2(cross3(A_ref - A, C0 - A), tol, "hinge normal")
  e2 <- unit2(cross3(nC, e1), tol, "hinge basis e2")
  
  to_C_plane2d <- function(P) {
    v <- P - A
    c(sum(v * e1), sum(v * e2))
  }
  
  from_C_plane2d <- function(uv) {
    A + uv[1] * e1 + uv[2] * e2
  }
  
  project_to_C_plane <- function(P) {
    v <- P - A
    A + sum(v * e1) * e1 + sum(v * e2) * e2
  }
  
  axis_AC0 <- unit2(A - C0, tol, "A - C0")
  relB0    <- B0 - C0
  
  O0 <- c(B0[1], 0, B0[3])
  rE0 <- sqrt((E0[1] - O0[1])^2 + (E0[3] - O0[3])^2)
  phi0 <- atan2(E0[3] - O0[3], E0[1] - O0[1])
  
  structure(
    list(
      pts0 = pts0,
      A = A,
      A_ref = A_ref,
      B0 = B0,
      B0_ref = B0_ref,
      C0 = C0,
      D0 = D0,
      E0 = E0,
      AC = dist3(A, C0),
      DC = dist3(D0, C0),
      BC = dist3(B0, C0),
      BE = dist3(B0, E0),
      DE = dist3(D0, E0),
      axis_AC0 = axis_AC0,
      relB0 = relB0,
      E_phase0 = phi0,
      E_r0 = rE0,
      to_C_plane2d = to_C_plane2d,
      from_C_plane2d = from_C_plane2d,
      project_to_C_plane = project_to_C_plane,
      tol = tol
    ),
    class = "sym_linkage_model"
  )
}

# ============================================================
# Global pose
# ============================================================

#' Apply a global rotation to model points before solving a frame.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param model Linkage model object produced by the corresponding model-construction function.
#' @param global_rot_deg Global rotation angle in degrees.
#' @param exclude_from_global Character vector of model point names excluded from global rotation.
#' @param global_axis Two point names defining the global rotation axis.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
apply_global_pose <- function(model,
                              global_rot_deg = 0,
                              exclude_from_global = character(),
                              global_axis = c("J", "J_ref"),
                              tol = NULL) {
  tol <- tol %||% model$tol
  pts <- model$pts0
  
  if (abs(global_rot_deg) > 0) {
    if (!all(global_axis %in% names(pts))) {
      stop("global_axis points not found in model$pts0.")
    }
    
    pts <- rotate_named_points_about_line(
      pts = pts,
      P1 = pts[[global_axis[1]]],
      P2 = pts[[global_axis[2]]],
      angle = global_rot_deg * pi / 180,
      exclude = unique(exclude_from_global),
      tol = tol
    )
  }
  
  pts
}

#' Rebuild a symmetric linkage model from a named point list.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param pts Named list or vector of point definitions.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
rebuild_model_from_points <- function(pts, tol = 1e-10) {
  make_sym_model(pts_to_coords(pts), tol = tol)
}

# PATCH: build a posed solver model that keeps original linkage lengths
# while updating orientation / plane geometry from the posed coordinates.
#' Create a posed solver model while preserving the original linkage lengths.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param model Linkage model object produced by the corresponding model-construction function.
#' @param pts_pose Named 3D points after global pose adjustment.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
make_posed_solver_model <- function(model, pts_pose, tol = NULL) {
  tol <- tol %||% model$tol
  
  pose_model <- rebuild_model_from_points(pts_pose, tol = tol)
  
  # Preserve original linkage constraints
  pose_model$AC <- model$AC
  pose_model$DC <- model$DC
  pose_model$BC <- model$BC
  pose_model$BE <- model$BE
  pose_model$DE <- model$DE
  
  pose_model
}

# ============================================================
# D path
# ============================================================

#' Build a frame-by-frame path for the driver point D.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param model Linkage model object produced by the corresponding model-construction function.
#' @param x_shift Driver-point x translation; either increments or absolute values.
#' @param z_shift Driver-point z translation; either increments or absolute values.
#' @param driver_mode Whether driver translations are stepwise increments or absolute offsets.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
build_D_path <- function(model,
                         x_shift,
                         z_shift = 0,
                         driver_mode = c("step", "absolute")) {
  driver_mode <- match.arg(driver_mode)
  
  n <- length(x_shift)
  if (length(z_shift) == 1) z_shift <- rep(z_shift, n)
  if (length(z_shift) != n) stop("z_shift must be length 1 or length(x_shift)")
  
  D_path <- vector("list", n)
  D_prev <- model$D0
  
  for (i in seq_len(n)) {
    if (driver_mode == "step") {
      D_now <- c(D_prev[1] + x_shift[i], model$D0[2], D_prev[3] + z_shift[i])
    } else {
      D_now <- c(model$D0[1] + x_shift[i], model$D0[2], model$D0[3] + z_shift[i])
    }
    D_path[[i]] <- D_now
    D_prev <- D_now
  }
  
  D_path
}

# ============================================================
# Geometric C path
# ============================================================

#' Solve the path of point C from the prescribed D path and fixed link lengths.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param model Linkage model object produced by the corresponding model-construction function.
#' @param D_path List of prescribed D-point coordinates through time.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
solve_C_path_geometric <- function(model, D_path, tol = NULL) {
  tol <- tol %||% model$tol
  
  prev_C <- model$C0
  C_path <- vector("list", length(D_path))
  
  for (i in seq_along(D_path)) {
    D <- D_path[[i]]
    
    Dp <- model$project_to_C_plane(D)
    h  <- dist3(D, Dp)
    
    rho2 <- model$DC^2 - h^2
    if (rho2 < -tol) {
      stop(sprintf("No real C solution at frame %d: D too far from C-motion plane.", i))
    }
    rho <- sqrt(max(0, rho2))
    
    A2 <- c(0, 0)
    D2 <- model$to_C_plane2d(Dp)
    
    sol <- circle_circle_2d(A2, model$AC, D2, rho, tol)
    if (!sol$ok) {
      stop(sprintf("No C solution at frame %d: %s", i, sol$reason))
    }
    
    Ccands <- lapply(seq_len(nrow(sol$pts)), function(j) {
      model$from_C_plane2d(sol$pts[j, ])
    })
    
    keep <- which.min(sapply(Ccands, function(x) dist3(x, prev_C)))
    C <- Ccands[[keep]]
    
    C_path[[i]] <- C
    prev_C <- C
  }
  
  C_path
}

# ============================================================
# B and E builders
# B rides rigidly with C
# E rotates independently about B
# ============================================================

#' Place point B from point C using the model reference geometry.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param model Linkage model object produced by the corresponding model-construction function.
#' @param C 3D point C, typically a linkage joint.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
build_B_from_C <- function(model, C, tol = NULL) {
  tol <- tol %||% model$tol
  
  axis_now  <- unit2(model$A - C, tol, "A - C current")
  rot_align <- make_axis_alignment_rotation(model$axis_AC0, axis_now, tol)
  
  relB_now <- rot_align(model$relB0)
  C + relB_now
}

#' Place point E by applying a prescribed rotation about B.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param model Linkage model object produced by the corresponding model-construction function.
#' @param B 3D point B, typically a linkage joint or angle vertex.
#' @param E_rot_deg Vector of prescribed rotations for point E in degrees.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
build_E_from_rotation <- function(model, B, E_rot_deg, tol = NULL) {
  tol <- tol %||% model$tol
  
  if (abs(B[2]) > model$BE + tol) {
    stop("No E solution: |B_y| > |BE|, so E cannot remain on y = 0.")
  }
  
  r <- sqrt(max(0, model$BE^2 - B[2]^2))
  phi <- model$E_phase0 + E_rot_deg * pi / 180
  O <- c(B[1], 0, B[3])
  
  E <- O + c(r * cos(phi), 0, r * sin(phi))
  E[2] <- 0
  E
}

# ============================================================
# Main path solver
# PATCH:
# Global pose updates the plane/orientation geometry, but original
# linkage lengths remain the solver constraints.
# ============================================================

#' Solve the full symmetric linkage path from prescribed driver translations and E rotations.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param model Linkage model object produced by the corresponding model-construction function.
#' @param x_shift Driver-point x translation; either increments or absolute values.
#' @param z_shift Driver-point z translation; either increments or absolute values.
#' @param E_rot_deg Vector of prescribed rotations for point E in degrees.
#' @param global_rot_deg Global rotation angle in degrees.
#' @param exclude_from_global Character vector of model point names excluded from global rotation.
#' @param global_axis Two point names defining the global rotation axis.
#' @param driver_mode Whether driver translations are stepwise increments or absolute offsets.
#' @param global_driver_mode Input argument used by this helper.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A list of solved linkage states.
#' @keywords internal
#' 
#' 
solve_dual5_path <- function(model,
                             x_shift,
                             z_shift = 0,
                             E_rot_deg = 0,
                             global_rot_deg = 0,
                             exclude_from_global = character(),
                             global_axis = c("J", "J_ref"),
                             driver_mode = c("step", "absolute"),
                             global_driver_mode = c("absolute", "step"),
                             tol = NULL) {
  tol <- tol %||% model$tol
  driver_mode <- match.arg(driver_mode)
  global_driver_mode <- match.arg(global_driver_mode)
  
  n <- length(x_shift)
  if (length(z_shift) == 1) z_shift <- rep(z_shift, n)
  if (length(E_rot_deg) == 1) E_rot_deg <- rep(E_rot_deg, n)
  if (length(global_rot_deg) == 1) global_rot_deg <- rep(global_rot_deg, n)
  
  if (length(z_shift) != n) stop("z_shift must be length 1 or length(x_shift)")
  if (length(E_rot_deg) != n) stop("E_rot_deg must be length 1 or length(x_shift)")
  if (length(global_rot_deg) != n) stop("global_rot_deg must be length 1 or length(x_shift)")
  
  global_rot_use <- if (global_driver_mode == "step") cumsum(global_rot_deg) else global_rot_deg
  
  sols <- vector("list", n)
  
  for (i in seq_len(n)) {
    pts_pose <- apply_global_pose(
      model = model,
      global_rot_deg = global_rot_use[i],
      exclude_from_global = exclude_from_global,
      global_axis = global_axis,
      tol = tol
    )
    
    pose_model <- make_posed_solver_model(
      model = model,
      pts_pose = pts_pose,
      tol = tol
    )
    
    D_path_i <- build_D_path(
      model = pose_model,
      x_shift = x_shift[seq_len(i)],
      z_shift = z_shift[seq_len(i)],
      driver_mode = driver_mode
    )
    
    D <- D_path_i[[i]]
    
    C_path_i <- solve_C_path_geometric(
      model = pose_model,
      D_path = D_path_i,
      tol = tol
    )
    
    C <- C_path_i[[i]]
    B <- build_B_from_C(pose_model, C, tol = tol)
    E <- build_E_from_rotation(pose_model, B, E_rot_deg[i], tol = tol)
    
    pts <- pts_pose
    pts$A <- pose_model$A
    pts$A_ref <- pose_model$A_ref
    pts$B <- B
    pts$B_ref <- reflect_xz(B)
    pts$C <- C
    pts$C_ref <- reflect_xz(C)
    pts$D <- D
    pts$E <- E
    
    pts$E_rot_deg <- E_rot_deg[i]
    pts$global_rot_deg <- global_rot_use[i]
    pts$AC <- dist3(pts$A, pts$C)
    pts$DC <- dist3(pts$D, pts$C)
    pts$BC <- dist3(pts$B, pts$C)
    pts$BE <- dist3(pts$B, pts$E)
    pts$DE <- dist3(pts$D, pts$E)
    pts$BBref <- dist3(pts$B, pts$B_ref)
    pts$ACB <- angle_3pt(pts$A, pts$C, pts$B, deg = TRUE)
    
    sols[[i]] <- pts
  }
  
  structure(
    list(
      solutions = sols
    ),
    class = "sym_path_result"
  )
}

# ============================================================
# Pipe-friendly wrapper
# ============================================================

#' Pipe-friendly wrapper for solving a path from a tibble of frame inputs.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param .data Input tibble containing frame-level driver variables.
#' @param model Linkage model object produced by the corresponding model-construction function.
#' @param global_rot_deg Global rotation angle in degrees.
#' @param exclude_from_global Character vector of model point names excluded from global rotation.
#' @param global_axis Two point names defining the global rotation axis.
#' @param driver_mode Whether driver translations are stepwise increments or absolute offsets.
#' @param global_driver_mode Input argument used by this helper.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return The input tibble with solution, ok, and diagnostic columns added.
#' @keywords internal
solve_path_tbl <- function(.data,
                           model,
                           global_rot_deg = NULL,
                           exclude_from_global = character(),
                           global_axis = c("J", "J_ref"),
                           driver_mode = c("step", "absolute"),
                           global_driver_mode = c("absolute", "step"),
                           tol = NULL) {
  driver_mode <- match.arg(driver_mode)
  global_driver_mode <- match.arg(global_driver_mode)
  dat <- as_tibble(.data)
  
  if (!"x_shift" %in% names(dat)) {
    stop("Input tibble must contain column `x_shift`.")
  }
  if (!"z_shift" %in% names(dat)) dat$z_shift <- 0
  if (!"E_rot_deg" %in% names(dat)) dat$E_rot_deg <- 0
  
  if (is.null(global_rot_deg)) {
    if (!"global_rot_deg" %in% names(dat)) dat$global_rot_deg <- 0
  } else {
    dat$global_rot_deg <- global_rot_deg
  }
  
  solve_dual5_path(
    model = model,
    x_shift = dat$x_shift,
    z_shift = dat$z_shift,
    E_rot_deg = dat$E_rot_deg,
    global_rot_deg = dat$global_rot_deg,
    exclude_from_global = exclude_from_global,
    global_axis = global_axis,
    driver_mode = driver_mode,
    global_driver_mode = global_driver_mode,
    tol = tol
  )
}

# ============================================================
# Output helpers
# ============================================================

#' Classify a point as original, reflected, or merged based on its name and y position.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param nm Name of a point to extract.
#' @param y Fallback value or coordinate component.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
infer_side <- function(nm, y, tol = 1e-10) {
  if (grepl("_ref$", nm)) return("reflected")
  if (abs(y) < tol) return("merged")
  "original"
}

#' Convert a list of solved linkage states to a long data frame of point coordinates.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param solutions List of solved linkage states, one element per frame.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
solutions_to_df <- function(solutions) {
  bind_rows(lapply(seq_along(solutions), function(i) {
    s <- solutions[[i]]
    nm_pts <- names(s)[vapply(s, is_xyz_point, logical(1))]
    
    bind_rows(lapply(nm_pts, function(nm) {
      P <- s[[nm]]
      tibble(
        frame = i,
        point = sub("_ref$", "", nm),
        name  = nm,
        x = P[1],
        y = P[2],
        z = P[3],
        side = infer_side(nm, P[2])
      )
    }))
  }))
}

#' Summarize key link distances and angles across solved frames for diagnostics.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param solutions List of solved linkage states, one element per frame.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
constraint_series <- function(solutions) {
  tibble(
    frame = seq_along(solutions),
    AC = sapply(solutions, function(s) dist3(s$A, s$C)),
    DC = sapply(solutions, function(s) dist3(s$D, s$C)),
    BC = sapply(solutions, function(s) dist3(s$B, s$C)),
    BE = sapply(solutions, function(s) dist3(s$B, s$E)),
    DE = sapply(solutions, function(s) dist3(s$D, s$E)),
    BBref = sapply(solutions, function(s) dist3(s$B, s$B_ref)),
    ACB = sapply(solutions, function(s) angle_3pt(s$A, s$C, s$B, deg = TRUE)),
    Ey = sapply(solutions, function(s) s$E[2]),
    E_rot_deg = sapply(solutions, function(s) s$E_rot_deg %||% NA_real_),
    global_rot_deg = sapply(solutions, function(s) s$global_rot_deg %||% NA_real_)
  )
}

#' Convert one solved linkage state to line-segment coordinates for plotting.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @param link_map Named list specifying which point pairs form plotted links.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
segments_from_solution <- function(sol, link_map) {
  bind_rows(lapply(names(link_map), function(nm) {
    pts <- link_map[[nm]]
    
    if (length(pts) != 2) return(NULL)
    if (!all(pts %in% names(sol))) return(NULL)
    if (!is_xyz_point(sol[[pts[1]]]) || !is_xyz_point(sol[[pts[2]]])) return(NULL)
    
    p1 <- sol[[pts[1]]]
    p2 <- sol[[pts[2]]]
    
    tibble(
      link = nm,
      x = c(p1[1], p2[1], NA_real_),
      y = c(p1[2], p2[2], NA_real_),
      z = c(p1[3], p2[3], NA_real_)
    )
  }))
}

#' Convert multiple solved linkage states to line-segment coordinates for plotting or animation.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param solutions List of solved linkage states, one element per frame.
#' @param link_map Named list specifying which point pairs form plotted links.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
segments_from_solutions <- function(solutions, link_map) {
  bind_rows(lapply(seq_along(solutions), function(i) {
    segments_from_solution(solutions[[i]], link_map) %>%
      mutate(frame = i)
  }))
}

#' Reset the starting spacing between C and C_ref while preserving symmetry of paired points.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param coords Data frame or tibble containing named coordinates with columns name, x, y, and z.
#' @param target_CCref Target starting distance between C and C_ref.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
set_start_C <- function(coords, target_CCref, tol = 1e-10) {
  dat <- as_tibble(coords)
  
  req_cols <- c("name", "x", "y", "z")
  if (!all(req_cols %in% names(dat))) {
    stop("coords must contain columns: name, x, y, z")
  }
  
  get_idx <- function(nm) {
    ii <- match(nm, dat$name)
    if (is.na(ii)) stop("Point not found: ", nm)
    ii
  }
  
  get_pt <- function(nm) {
    ii <- get_idx(nm)
    as.numeric(dat[ii, c("x", "y", "z")])
  }
  
  set_pt <- function(nm, val) {
    ii <- get_idx(nm)
    dat[ii, c("x", "y", "z")] <<- as.list(val)
  }
  
  pivot <- get_pt("A")
  B0    <- get_pt("B")
  C0    <- get_pt("C")
  axis  <- c(1, 0, 0)
  
  target_Cy <- target_CCref / 2
  
  k <- unit2(axis, tol, "set_start_C axis")
  v <- C0 - pivot
  
  v_par  <- k * sum(k * v)
  v_perp <- v - v_par
  
  e_y <- c(0, 1, 0)
  
  a  <- sum(e_y * v_perp)
  b  <- sum(e_y * cross3(k, v_perp))
  c0 <- pivot[2] + sum(e_y * v_par)
  
  R   <- sqrt(a^2 + b^2)
  rhs <- target_Cy - c0
  
  if (R < tol) {
    stop("Rotation about x-axis cannot change the y-coordinate of C.")
  }
  if (abs(rhs) > R + tol) {
    stop("Requested target_CCref is unreachable by rigid rotation from the current pose.")
  }
  
  rhs <- max(min(rhs, R), -R)
  phi <- atan2(b, a)
  alpha <- acos(rhs / R)
  
  theta_candidates <- c(phi + alpha, phi - alpha)
  
  rotate_line <- function(P, theta) {
    rotate_point_about_line(
      P = P,
      P1 = pivot,
      P2 = pivot + axis,
      angle = theta,
      tol = tol
    )
  }
  
  C_cands <- lapply(theta_candidates, function(th) rotate_line(C0, th))
  keep <- which.min(vapply(C_cands, function(x) sum((x - C0)^2), numeric(1)))
  theta <- theta_candidates[keep]
  
  B_new <- rotate_line(B0, theta)
  C_new <- rotate_line(C0, theta)
  
  set_pt("B", B_new)
  set_pt("C", C_new)
  
  if ("B_ref" %in% dat$name) set_pt("B_ref", reflect_xz(B_new))
  if ("C_ref" %in% dat$name) set_pt("C_ref", reflect_xz(C_new))
  
  attr(dat, "init_theta_rad") <- theta
  attr(dat, "init_theta_deg") <- theta * 180 / pi
  attr(dat, "target_CCref") <- target_CCref
  
  dat
}

# ============================================================
# Static 3D plotter
# ============================================================

#' Plot original and reflected coordinates with optional linkage segments in 3D.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param all_coords Coordinate data frame containing original/reflected/merged points.
#' @param link_map Named list specifying which point pairs form plotted links.
#' @param link_col Input argument used by this helper.
#' @param point_size Input argument used by this helper.
#' @param line_width Input argument used by this helper.
#' @param show_labels Input argument used by this helper.
#' @param pitch Camera pitch angle in degrees.
#' @param yaw Camera yaw angle in degrees.
#' @param zoom Camera distance or scale factor.
#' @return A Plotly figure.
#' @keywords internal
plot_symmetry_3d <- function(all_coords,
                             link_map = NULL,
                             link_col = NULL,
                             point_size = 5,
                             line_width = 6,
                             show_labels = TRUE,
                             pitch = 15,
                             yaw = 110,
                             zoom = 2) {
  dat <- as_tibble(all_coords)
  
  req_cols <- c("name", "x", "y", "z")
  if (!all(req_cols %in% names(dat))) {
    stop("all_coords must contain columns: name, x, y, z")
  }
  
  if (!"side" %in% names(dat)) {
    dat <- dat %>%
      mutate(side = case_when(
        grepl("_ref$", name) ~ "reflected",
        abs(y) < 1e-10 ~ "merged",
        TRUE ~ "original"
      ))
  }
  
  link_df <- NULL
  if (!is.null(link_map)) {
    link_df <- bind_rows(lapply(names(link_map), function(lnk) {
      pts <- link_map[[lnk]]
      if (length(pts) != 2) return(NULL)
      if (!all(pts %in% dat$name)) return(NULL)
      
      p1 <- dat %>% filter(name == pts[1]) %>% slice(1)
      p2 <- dat %>% filter(name == pts[2]) %>% slice(1)
      
      tibble(
        link = lnk,
        x = c(p1$x, p2$x, NA_real_),
        y = c(p1$y, p2$y, NA_real_),
        z = c(p1$z, p2$z, NA_real_)
      )
    }))
  }
  
  pitch_r <- pitch * pi / 180
  yaw_r   <- yaw * pi / 180
  
  eye <- list(
    x = zoom * cos(pitch_r) * cos(yaw_r),
    y = zoom * cos(pitch_r) * sin(yaw_r),
    z = zoom * sin(pitch_r)
  )
  
  p <- plotly::plot_ly()
  
  for (s in unique(dat$side)) {
    dsub <- dat %>% filter(side == s)
    
    p <- plotly::add_markers(
      p,
      data = dsub,
      x = ~x, y = ~y, z = ~z,
      split = ~side,
      text = ~name,
      hovertemplate = paste0(
        "<b>%{text}</b><br>",
        "x: %{x:.3f}<br>",
        "y: %{y:.3f}<br>",
        "z: %{z:.3f}<br>",
        "side: ", s, "<extra></extra>"
      ),
      marker = list(size = point_size),
      showlegend = TRUE
    )
    
    if (show_labels) {
      p <- plotly::add_text(
        p,
        data = dsub,
        x = ~x, y = ~y, z = ~z,
        text = ~name,
        textposition = "top center",
        showlegend = FALSE
      )
    }
  }
  
  if (!is.null(link_df) && nrow(link_df) > 0) {
    if (is.null(link_col)) {
      p <- plotly::add_trace(
        p,
        data = link_df,
        x = ~x, y = ~y, z = ~z,
        type = "scatter3d",
        mode = "lines",
        split = ~link,
        line = list(width = line_width),
        hovertemplate = "<b>%{fullData.name}</b><extra></extra>",
        showlegend = TRUE
      )
    } else {
      for (lnk in unique(link_df$link)) {
        dsub <- link_df %>% filter(link == lnk)
        p <- plotly::add_trace(
          p,
          data = dsub,
          x = ~x, y = ~y, z = ~z,
          type = "scatter3d",
          mode = "lines",
          name = lnk,
          line = list(width = line_width, color = unname(link_col[lnk] %||% "black")),
          hovertemplate = paste0("<b>", lnk, "</b><extra></extra>"),
          showlegend = TRUE
        )
      }
    }
  }
  
  p %>%
    plotly::layout(
      scene = list(
        xaxis = list(title = "X"),
        yaxis = list(title = "Y"),
        zaxis = list(title = "Z"),
        aspectmode = "data",
        camera = list(eye = eye)
      ),
      legend = list(orientation = "h")
    )
}

# ============================================================
# Animation helpers
# ============================================================

#' Compute padded x, y, and z ranges for plotting a set of solutions.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param solutions List of solved linkage states, one element per frame.
#' @param omit_points Optional point names to omit from range calculations or plots.
#' @param pad Proportional padding added to plot ranges.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
get_xyz_ranges <- function(solutions, omit_points = NULL, pad = 0.05) {
  pts_df <- solutions_to_df(solutions)
  
  if (!is.null(omit_points)) {
    pts_df <- pts_df %>% dplyr::filter(!name %in% omit_points)
  }
  
  rx <- range(pts_df$x, na.rm = TRUE)
  ry <- range(pts_df$y, na.rm = TRUE)
  rz <- range(pts_df$z, na.rm = TRUE)
  
  padx <- diff(rx) * pad
  pady <- diff(ry) * pad
  padz <- diff(rz) * pad
  
  if (padx == 0) padx <- 1
  if (pady == 0) pady <- 1
  if (padz == 0) padz <- 1
  
  list(
    x = c(rx[1] - padx, rx[2] + padx),
    y = c(ry[1] - pady, ry[2] + pady),
    z = c(rz[1] - padz, rz[2] + padz)
  )
}

#' Convert a camera roll angle to a Plotly up-vector.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param roll Camera roll angle in degrees.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
roll_to_up <- function(roll = 0) {
  rr <- roll * pi / 180
  list(
    x = -sin(rr),
    y =  cos(rr),
    z =  0
  )
}

#' Animate a symmetric linkage solution path in Plotly.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param solutions List of solved linkage states, one element per frame.
#' @param link_map Named list specifying which point pairs form plotted links.
#' @param link_col Input argument used by this helper.
#' @param omit_points Optional point names to omit from range calculations or plots.
#' @param label_points Input argument used by this helper.
#' @param point_col Input argument used by this helper.
#' @param point_size Input argument used by this helper.
#' @param line_width Input argument used by this helper.
#' @param show_labels Input argument used by this helper.
#' @param pitch Camera pitch angle in degrees.
#' @param yaw Camera yaw angle in degrees.
#' @param roll Camera roll angle in degrees.
#' @param zoom Camera distance or scale factor.
#' @param frame_duration Input argument used by this helper.
#' @param axis_pad Input argument used by this helper.
#' @param remove_axes Logical; if TRUE, hide 3D axes.
#' @param bgcolor Scene background color.
#' @return A Plotly animation object.
#' @keywords internal
animate_sym_linkage <- function(solutions,
                                link_map = NULL,
                                link_col = NULL,
                                omit_points = NULL,
                                label_points = NULL,
                                point_col = "black",
                                point_size = 4,
                                line_width = 6,
                                show_labels = TRUE,
                                pitch = 15,
                                yaw = 110,
                                roll = 0,
                                zoom = 2,
                                frame_duration = 50,
                                axis_pad = 0.05,
                                remove_axes = FALSE,
                                bgcolor = "white") {
  pts_df <- solutions_to_df(solutions)
  
  if (!is.null(omit_points)) {
    pts_df <- pts_df %>%
      dplyr::filter(!name %in% omit_points)
  }
  
  pts_df <- pts_df %>%
    dplyr::mutate(
      label = if (isTRUE(show_labels)) {
        if (is.null(label_points)) name else ifelse(name %in% label_points, name, "")
      } else {
        ""
      }
    )
  
  link_map_use <- link_map
  if (!is.null(link_map_use) && !is.null(omit_points)) {
    keep_link <- vapply(
      link_map_use,
      function(x) !any(x %in% omit_points),
      logical(1)
    )
    link_map_use <- link_map_use[keep_link]
  }
  
  link_df <- NULL
  if (!is.null(link_map_use) && length(link_map_use) > 0) {
    link_df <- segments_from_solutions(solutions, link_map_use)
  }
  
  bounds <- get_xyz_ranges(solutions, omit_points = omit_points, pad = axis_pad)
  
  pitch_r <- pitch * pi / 180
  yaw_r   <- yaw * pi / 180
  
  eye <- list(
    x = zoom * cos(pitch_r) * cos(yaw_r),
    y = zoom * cos(pitch_r) * sin(yaw_r),
    z = zoom * sin(pitch_r)
  )
  
  axis_x <- if (remove_axes) {
    list(title = "", showbackground = FALSE, showgrid = FALSE, zeroline = FALSE,
         showticklabels = FALSE, ticks = "", showline = FALSE, range = bounds$x)
  } else {
    list(showbackground = FALSE, range = bounds$x)
  }
  
  axis_y <- if (remove_axes) {
    list(title = "", showbackground = FALSE, showgrid = FALSE, zeroline = FALSE,
         showticklabels = FALSE, ticks = "", showline = FALSE, range = bounds$y)
  } else {
    list(showbackground = FALSE, range = bounds$y)
  }
  
  axis_z <- if (remove_axes) {
    list(title = "", showbackground = FALSE, showgrid = FALSE, zeroline = FALSE,
         showticklabels = FALSE, ticks = "", showline = FALSE, range = bounds$z)
  } else {
    list(showbackground = FALSE, range = bounds$z)
  }
  
  p <- plotly::plot_ly()
  
  if (!is.null(link_df) && nrow(link_df) > 0) {
    for (lnk in unique(link_df$link)) {
      dsub <- link_df %>% dplyr::filter(link == lnk)
      
      p <- plotly::add_trace(
        p,
        data = dsub,
        x = ~x, y = ~y, z = ~z,
        frame = ~frame,
        type = "scatter3d",
        mode = "lines",
        name = lnk,
        line = if (is.null(link_col)) {
          list(width = line_width)
        } else {
          list(width = line_width, color = unname(link_col[lnk] %||% "black"))
        },
        hovertemplate = paste0("<b>", lnk, "</b><extra></extra>"),
        showlegend = FALSE
      )
    }
  }
  
  p <- plotly::add_trace(
    p,
    data = pts_df,
    x = ~x, y = ~y, z = ~z,
    frame = ~frame,
    type = "scatter3d",
    mode = if (show_labels) "markers+text" else "markers",
    text = ~label,
    customdata = ~name,
    textposition = "top center",
    marker = list(size = point_size, color = point_col),
    hovertemplate = paste0(
      "<b>%{customdata}</b><br>",
      "x: %{x:.3f}<br>",
      "y: %{y:.3f}<br>",
      "z: %{z:.3f}<extra></extra>"
    ),
    showlegend = FALSE
  )
  
  p %>%
    plotly::layout(
      scene = list(
        xaxis = axis_x,
        yaxis = axis_y,
        zaxis = axis_z,
        aspectmode = "data",
        camera = list(
          eye = eye,
          up = roll_to_up(roll)
        ),
        bgcolor = bgcolor
      ),
      showlegend = FALSE,
      uirevision = TRUE
    ) %>%
    plotly::animation_opts(
      frame = frame_duration,
      transition = 0,
      redraw = TRUE
    )
}