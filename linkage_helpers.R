# ============================================================
# General linkage helpers, plotting utilities, and geometry tools
#
# This file contains reusable geometry, plotting, plane-casting, volume,
# and trajectory helpers used by the linkage workflow. Function comments
# use roxygen-style tags so the helper file can be migrated into an R
# package with minimal restructuring.
# ============================================================

suppressPackageStartupMessages({
  library(plotly)
  library(tidyverse)
})

# ============================================================
# Basic helpers
# ============================================================

#' Compute the Euclidean norm of a numeric vector.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param x Numeric vector, coordinate, or value used by the helper.
#' @return A single numeric norm.
#' @keywords internal
nrm <- function(x) sqrt(sum(x^2))

unit <- function(x, tol = 1e-12) {
  d <- nrm(x)
  if (d < tol) stop("Cannot normalize near-zero vector.")
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
dist3 <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  sqrt(sum((a - b)^2))
}

#' Compute the 3D cross product of two vectors.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param a First vector, point, or object used in the calculation.
#' @param b Second vector, point, or object used in the calculation.
#' @return A length-3 numeric vector.
#' @keywords internal
cross3 <- function(a, b) {
  c(
    a[2]*b[3] - a[3]*b[2],
    a[3]*b[1] - a[1]*b[3],
    a[1]*b[2] - a[2]*b[1]
  )
}

#' Return a fallback value when the first value is NULL.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param x Numeric vector, coordinate, or value used by the helper.
#' @param y Fallback value or coordinate component.
#' @return Either x or y.
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

# ============================================================
# Circle-circle intersection in 2D
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
  
  if (h < tol) {
    pts <- matrix(p0, nrow = 1)
  } else {
    pts <- rbind(
      p0 + h * ey,
      p0 - h * ey
    )
  }
  
  list(ok = TRUE, pts = pts, reason = NULL)
}

# ============================================================
# Build model
# ============================================================

#' Build a geometric dual five-bar linkage model from starting landmark coordinates.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param A 3D point A, typically a fixed or anatomical landmark.
#' @param D 3D point D, typically a driver or fixed landmark.
#' @param B0 Starting/reference coordinate for point B.
#' @param C0 Starting/reference coordinate for point C.
#' @param E0 Starting/reference coordinate for point E.
#' @param M0 Starting/reference coordinate for point M.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A list with class dual5bar_model containing starting geometry, link lengths, and projection helpers.
#' @keywords internal
make_dual5bar_model <- function(A, D, B0, C0, E0, M0, tol = 1e-10) {
  
  # fixed plate plane from A, D, B0
  e1 <- unit(D - A, tol)
  n_plane <- cross3(D - A, B0 - A)
  if (nrm(n_plane) < tol) stop("A, D, B0 are collinear.")
  e3 <- unit(n_plane, tol)
  e2 <- unit(cross3(e3, e1), tol)
  
  to_plane2d <- function(P) {
    v <- P - A
    c(sum(v * e1), sum(v * e2))
  }
  
  from_plane2d <- function(xy) {
    A + xy[1] * e1 + xy[2] * e2
  }
  
  project_to_plane <- function(P) {
    v <- P - A
    A + sum(v * e1) * e1 + sum(v * e2) * e2
  }
  
  model <- list(
    A = A, D = D,
    B0 = B0, C0 = C0, E0 = E0, M0 = M0,
    
    plate_e1 = e1,
    plate_e2 = e2,
    plate_e3 = e3,
    
    to_plane2d = to_plane2d,
    from_plane2d = from_plane2d,
    project_to_plane = project_to_plane,
    
    # fixed bars
    AB = dist3(A, B0),
    DC = dist3(D, C0),
    
    # M links
    BM = dist3(B0, M0),
    CM = dist3(C0, M0),
    
    # E links
    BE = dist3(B0, E0),
    CE = dist3(C0, E0),
    
    tol = tol
  )
  
  class(model) <- "dual5bar_model"
  model
}

# ============================================================
# Solve B from M
# ============================================================

#' Solve the position of point B from a prescribed driver point M.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param model Linkage model object produced by the corresponding model-construction function.
#' @param M Current 3D coordinate of driver point M.
#' @param prev_B Previous B coordinate used to choose the continuous solution branch.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
solve_B_from_M <- function(model, M, prev_B = NULL, tol = NULL) {
  tol <- tol %||% model$tol
  
  Mp <- model$project_to_plane(M)
  h  <- dist3(M, Mp)
  
  rho2 <- model$BM^2 - h^2
  if (rho2 < -tol) stop("No real B solution: M too far from plate plane.")
  rho <- sqrt(max(0, rho2))
  
  A2  <- model$to_plane2d(model$A)
  Mp2 <- model$to_plane2d(Mp)
  
  sol <- circle_circle_2d(A2, model$AB, Mp2, rho, tol)
  if (!sol$ok) stop(paste("No B solution:", sol$reason))
  
  Bcands <- lapply(seq_len(nrow(sol$pts)), function(i) model$from_plane2d(sol$pts[i, ]))
  
  if (length(Bcands) == 1) return(Bcands[[1]])
  
  target <- prev_B %||% model$B0
  d <- sapply(Bcands, function(x) dist3(x, target))
  Bcands[[which.min(d)]]
}

# ============================================================
# Solve C from M
# ============================================================

#' Solve points B, C, and E from a prescribed driver point M.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param model Linkage model object produced by the corresponding model-construction function.
#' @param M Current 3D coordinate of driver point M.
#' @param prev Previous solved frame used to choose the continuous solution branch.
#' @param y_sym Y coordinate of the imposed symmetry plane; default is 0.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
solve_BCE_from_M <- function(model, M, prev = NULL, y_sym = 0, tol = NULL) {
  tol <- tol %||% model$tol
  
  B <- solve_B_from_M(
    model, M,
    prev_B = if (!is.null(prev)) prev$B else NULL,
    tol = tol
  )
  
  C <- solve_C_from_M(
    model, M,
    prev_C = if (!is.null(prev)) prev$C else NULL,
    tol = tol
  )
  
  E <- solve_E_from_BC_sym(
    model, B, C,
    prev_E = if (!is.null(prev)) prev$E else NULL,
    tol = tol
  )
  
  list(
    A = model$A,
    D = model$D,
    B = B,
    C = C,
    E = E,
    M = M
  )
}
# ============================================================
# Solve E from B and C, constrained to y = y_sym
# Default is symmetry plane y = 0
# ============================================================

#' Solve the midline point E from the current B and C positions under symmetry constraints.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param model Linkage model object produced by the corresponding model-construction function.
#' @param B 3D point B, typically a linkage joint or angle vertex.
#' @param C 3D point C, typically a linkage joint.
#' @param y_sym Y coordinate of the imposed symmetry plane; default is 0.
#' @param prev_E Previous E coordinate used to choose the continuous solution branch.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
solve_E_from_BC_sym <- function(model, B, C, y_sym = 0, prev_E = NULL, tol = NULL) {
  tol <- tol %||% model$tol
  
  # Intersect:
  # (x-Bx)^2 + (y_sym-By)^2 + (z-Bz)^2 = BE^2
  # (x-Cx)^2 + (y_sym-Cy)^2 + (z-Cz)^2 = CE^2
  #
  # This is a circle-circle problem in the x-z plane at fixed y = y_sym.
  
  B2 <- c(B[1], B[3])
  C2 <- c(C[1], C[3])
  
  rB2 <- model$BE^2 - (y_sym - B[2])^2
  rC2 <- model$CE^2 - (y_sym - C[2])^2
  
  if (rB2 < -tol) stop("No real E solution from B: symmetry-plane slice misses sphere.")
  if (rC2 < -tol) stop("No real E solution from C: symmetry-plane slice misses sphere.")
  
  rB <- sqrt(max(0, rB2))
  rC <- sqrt(max(0, rC2))
  
  sol <- circle_circle_2d(B2, rB, C2, rC, tol)
  if (!sol$ok) stop(paste("No E solution:", sol$reason))
  
  Ecands <- lapply(seq_len(nrow(sol$pts)), function(i) {
    c(sol$pts[i, 1], y_sym, sol$pts[i, 2])
  })
  
  if (length(Ecands) == 1) return(Ecands[[1]])
  
  target <- prev_E %||% model$E0
  d <- sapply(Ecands, function(x) dist3(x, target))
  Ecands[[which.min(d)]]
}

# ============================================================
# Solve full state from M
# ============================================================

#' Solve the midline point E from the current B and C positions under symmetry constraints.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param model Linkage model object produced by the corresponding model-construction function.
#' @param B 3D point B, typically a linkage joint or angle vertex.
#' @param C 3D point C, typically a linkage joint.
#' @param prev_E Previous E coordinate used to choose the continuous solution branch.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
solve_E_from_BC_sym <- function(model, B, C, prev_E = NULL, tol = NULL) {
  tol <- tol %||% model$tol
  
  # midpoint of BC
  G  <- 0.5 * (B + C)
  G0 <- 0.5 * (model$B0 + model$C0)
  
  # symmetry-plane direction from midpoint to E in reference pose
  v0 <- model$E0 - G0
  
  # project that reference direction into y = 0 plane
  dir <- c(v0[1], 0, v0[3])
  if (nrm(dir) < tol) stop("Reference midpoint-to-E direction is degenerate.")
  dir <- unit(dir, tol)
  
  # Because E lies in y = 0 and |BE| is fixed:
  #   (Ex-Bx)^2 + (0-By)^2 + (Ez-Bz)^2 = BE^2
  # Let r be the in-plane distance from G to E in the y=0 plane.
  # Then:
  #   r^2 + By^2 = BE^2
  #
  # Using B only is enough by symmetry.
  r2 <- model$BE^2 - B[2]^2
  if (r2 < -tol) stop("No real E solution: |By| exceeds |BE|.")
  r <- sqrt(max(0, r2))
  
  # two symmetric possibilities along the reference direction
  E1 <- G + r * dir
  E2 <- G - r * dir
  
  # choose branch by continuity, otherwise closest to reference side
  if (!is.null(prev_E)) {
    E <- if (dist3(E1, prev_E) <= dist3(E2, prev_E)) E1 else E2
  } else {
    E <- if (dist3(E1, model$E0) <= dist3(E2, model$E0)) E1 else E2
  }
  
  E[2] <- 0
  E
}
# ============================================================
# Solve path
# ============================================================

#' Solve a sequence of linkage states from a path of driver-point M positions.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param model Linkage model object produced by the corresponding model-construction function.
#' @param M_path Matrix or list of driver-point M coordinates through time.
#' @param y_sym Y coordinate of the imposed symmetry plane; default is 0.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
solve_path_from_M <- function(model, M_path, y_sym = 0, tol = NULL) {
  tol <- tol %||% model$tol
  
  if (is.vector(M_path)) M_path <- matrix(M_path, ncol = 3)
  
  sols <- vector("list", nrow(M_path))
  prev <- NULL
  
  for (i in seq_len(nrow(M_path))) {
    sols[[i]] <- solve_BCE_from_M(
      model,
      M = M_path[i, ],
      prev = prev,
      y_sym = y_sym,
      tol = tol
    )
    prev <- sols[[i]]
  }
  
  sols
}

# ============================================================
# Data helpers
# ============================================================

#' Convert a list of solved linkage states to a long data frame of point coordinates.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param solutions List of solved linkage states, one element per frame.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
solutions_to_df <- function(solutions) {
  do.call(rbind, lapply(seq_along(solutions), function(i) {
    s <- solutions[[i]]
    rbind(
      data.frame(frame = i, point = "A", x = s$A[1], y = s$A[2], z = s$A[3]),
      data.frame(frame = i, point = "D", x = s$D[1], y = s$D[2], z = s$D[3]),
      data.frame(frame = i, point = "B", x = s$B[1], y = s$B[2], z = s$B[3]),
      data.frame(frame = i, point = "C", x = s$C[1], y = s$C[2], z = s$C[3]),
      data.frame(frame = i, point = "E", x = s$E[1], y = s$E[2], z = s$E[3]),
      data.frame(frame = i, point = "M", x = s$M[1], y = s$M[2], z = s$M[3])
    )
  }))
}

#' Convert one solved linkage state to line-segment coordinates for plotting.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @param link_map Named list specifying which point pairs form plotted links.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
segments_from_solution <- function(sol, link_map = NULL) {
  
  pts <- list(
    A = sol$A, D = sol$D, B = sol$B,
    C = sol$C, E = sol$E, M = sol$M
  )
  
  # default links
  if (is.null(link_map)) {
    link_map <- list(
      AB = c("A","B"),
      DC = c("D","C"),
      BM = c("B","M"),
      CM = c("C","M"),
      BE = c("B","E"),
      CE = c("C","E"),
      AD = c("A","D")
    )
  }
  
  do.call(rbind, lapply(names(link_map), function(nm) {
    p1 <- pts[[link_map[[nm]][1]]]
    p2 <- pts[[link_map[[nm]][2]]]
    
    data.frame(
      link = nm,
      x = c(p1[1], p2[1], NA),
      y = c(p1[2], p2[2], NA),
      z = c(p1[3], p2[3], NA)
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
segments_from_solutions <- function(solutions, link_map = NULL) {
  do.call(rbind, lapply(seq_along(solutions), function(i) {
    d <- segments_from_solution(solutions[[i]], link_map = link_map)
    d$frame <- i
    d
  }))
}


# ============================================================
# Plot helpers
# ============================================================
#' Create a Plotly 3D camera specification from pitch, yaw, roll, and zoom values.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param pitch Camera pitch angle in degrees.
#' @param yaw Camera yaw angle in degrees.
#' @param roll Camera roll angle in degrees.
#' @param zoom Camera distance or scale factor.
#' @param center Camera center as a length-3 numeric vector.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
camera_from_euler <- function(pitch = 20,
                              yaw = 35,
                              roll = 0,
                              zoom = 1.8,
                              center = c(0, 0, 0)) {
  # angles in degrees
  
  deg2rad <- pi / 180
  pitch <- pitch * deg2rad
  yaw   <- yaw   * deg2rad
  roll  <- roll  * deg2rad
  
  # eye from spherical coords
  ex <- zoom * cos(pitch) * cos(yaw)
  ey <- zoom * cos(pitch) * sin(yaw)
  ez <- zoom * sin(pitch)
  
  # rolled "up" vector
  up0 <- c(0, 0, 1)
  view <- c(ex, ey, ez)
  view <- view / sqrt(sum(view^2))
  
  # rotate up-vector around view axis by roll
  cross3_local <- function(a, b) {
    c(
      a[2]*b[3] - a[3]*b[2],
      a[3]*b[1] - a[1]*b[3],
      a[1]*b[2] - a[2]*b[1]
    )
  }
  
  rotate_about_axis <- function(v, axis, theta) {
    axis <- axis / sqrt(sum(axis^2))
    v*cos(theta) +
      cross3_local(axis, v)*sin(theta) +
      axis*sum(axis*v)*(1 - cos(theta))
  }
  
  up <- rotate_about_axis(up0, view, roll)
  
  list(
    eye = list(x = ex, y = ey, z = ez),
    up = list(x = up[1], y = up[2], z = up[3]),
    center = list(x = center[1], y = center[2], z = center[3])
  )
}

#' Compute padded plotting bounds from a set of linkage solutions.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param solutions List of solved linkage states, one element per frame.
#' @param pad Proportional padding added to plot ranges.
#' @param omit_points Optional point names to omit from range calculations or plots.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
get_bounds <- function(solutions, pad = 0.05, omit_points = NULL) {
  pts <- solutions_to_df(solutions)
  
  if (!is.null(omit_points)) {
    pts <- pts %>% dplyr::filter(!name %in% omit_points)
  }
  
  expand <- function(r) {
    d <- diff(r)
    if (!is.finite(d) || d == 0) d <- 1
    c(r[1] - pad * d, r[2] + pad * d)
  }
  
  list(
    x = expand(range(pts$x, na.rm = TRUE)),
    y = expand(range(pts$y, na.rm = TRUE)),
    z = expand(range(pts$z, na.rm = TRUE))
  )
}
#' Return a simple default Plotly 3D camera position.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param x Numeric vector, coordinate, or value used by the helper.
#' @param y Fallback value or coordinate component.
#' @param z Input argument used by this helper.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
default_camera <- function(x = 1.5, y = 1.5, z = 1.2) {
  list(eye = list(x = x, y = y, z = z))
}

#' Build a reusable Plotly scene object with optional fixed camera and axis settings.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param bounds Optional x/y/z plot bounds.
#' @param remove_axes Logical; if TRUE, hide 3D axes.
#' @param fix_camera Logical; if TRUE, use a fixed camera.
#' @param camera Optional Plotly camera object.
#' @param pitch Camera pitch angle in degrees.
#' @param yaw Camera yaw angle in degrees.
#' @param roll Camera roll angle in degrees.
#' @param zoom Camera distance or scale factor.
#' @param bgcolor Scene background color.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
build_scene <- function(bounds = NULL,
                        remove_axes = TRUE,
                        fix_camera = TRUE,
                        camera = NULL,
                        pitch = NULL,
                        yaw = NULL,
                        roll = NULL,
                        zoom = 1.8,
                        bgcolor = "white") {
  
  if (is.null(camera) && (!is.null(pitch) || !is.null(yaw) || !is.null(roll))) {
    camera <- camera_from_euler(
      pitch = pitch %||% 20,
      yaw   = yaw   %||% 35,
      roll  = roll  %||% 0,
      zoom  = zoom
    )
  }
  
  if (is.null(camera)) {
    camera <- default_camera()
  }
  
  ax <- list(
    visible = !remove_axes,
    showbackground = FALSE,
    zeroline = FALSE,
    showgrid = FALSE,
    showticklabels = FALSE
  )
  
  scene <- list(
    xaxis = ax,
    yaxis = ax,
    zaxis = ax,
    aspectmode = "data",
    bgcolor = bgcolor
  )
  
  if (!is.null(bounds)) {
    scene$xaxis$range <- bounds$x
    scene$yaxis$range <- bounds$y
    scene$zaxis$range <- bounds$z
  }
  
  if (fix_camera) {
    scene$camera <- camera
  }
  
  scene
}


#' Plot a single linkage configuration as a static Plotly 3D object.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @param link_colors Named vector of colors for links.
#' @param link_widths Named vector of line widths for links.
#' @param title Plot title.
#' @param show_axes Logical; if TRUE, draw reference axes.
#' @param axes_origin Optional origin for plotted reference axes.
#' @param axes_length Length of plotted reference axes.
#' @param axes_colors Named vector of axis colors.
#' @param axes_width Line width for plotted axes.
#' @param axes_alpha Transparency for plotted axes.
#' @param axes_labels Logical; if TRUE, label plotted axes.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
plot_linkage_static <- function(sol,
                                link_colors = NULL,
                                link_widths = NULL,
                                title = "Dual 5-bar linkage",
                                show_axes = FALSE,
                                axes_origin = NULL,
                                axes_length = NULL,
                                axes_colors = c(
                                  x = "red",
                                  y = "forestgreen",
                                  z = "dodgerblue"
                                ),
                                axes_width = 4,
                                axes_alpha = 0.5,
                                axes_labels = TRUE) {
  
  segs <- segments_from_solution(sol)
  
  # default colors
  if (is.null(link_colors)) {
    link_colors <- c(
      AB = "black",
      DC = "black",
      BM = "red",
      CM = "red",
      BE = "blue",
      CE = "blue",
      AD = "gray"
    )
  }
  
  if (is.null(link_widths)) {
    link_widths <- setNames(rep(6, length(link_colors)), names(link_colors))
  }
  
  pts <- data.frame(
    point = c("A","D","B","C","E","M"),
    x = c(sol$A[1], sol$D[1], sol$B[1], sol$C[1], sol$E[1], sol$M[1]),
    y = c(sol$A[2], sol$D[2], sol$B[2], sol$C[2], sol$E[2], sol$M[2]),
    z = c(sol$A[3], sol$D[3], sol$B[3], sol$C[3], sol$E[3], sol$M[3])
  )
  
  p <- plot_ly()
  
  # ------------------------------------------------------------
  # Optional xyz axes in background
  # ------------------------------------------------------------
  if (show_axes) {
    
    xr <- range(pts$x, na.rm = TRUE)
    yr <- range(pts$y, na.rm = TRUE)
    zr <- range(pts$z, na.rm = TRUE)
    
    span <- max(diff(xr), diff(yr), diff(zr))
    if (!is.finite(span) || span <= 0) span <- 1
    
    if (is.null(axes_length)) {
      axes_length <- 0.18 * span
    }
    
    if (is.null(axes_origin)) {
      axes_origin <- c(
        xr[1] - 0.12 * span,
        yr[1] - 0.12 * span,
        zr[1] - 0.12 * span
      )
    }
    
    axes_df <- rbind(
      data.frame(
        axis = "x",
        x = c(axes_origin[1], axes_origin[1] + axes_length, NA),
        y = c(axes_origin[2], axes_origin[2], NA),
        z = c(axes_origin[3], axes_origin[3], NA)
      ),
      data.frame(
        axis = "y",
        x = c(axes_origin[1], axes_origin[1], NA),
        y = c(axes_origin[2], axes_origin[2] + axes_length, NA),
        z = c(axes_origin[3], axes_origin[3], NA)
      ),
      data.frame(
        axis = "z",
        x = c(axes_origin[1], axes_origin[1], NA),
        y = c(axes_origin[2], axes_origin[2], NA),
        z = c(axes_origin[3], axes_origin[3] + axes_length, NA)
      )
    )
    
    for (ax in c("x", "y", "z")) {
      d <- axes_df[axes_df$axis == ax, ]
      
      p <- add_trace(
        p,
        data = d,
        x = ~x, y = ~y, z = ~z,
        type = "scatter3d",
        mode = "lines",
        line = list(
          color = scales::alpha(axes_colors[[ax]], axes_alpha),
          width = axes_width
        ),
        showlegend = FALSE,
        hoverinfo = "none"
      )
    }
    
    if (axes_labels) {
      axes_lab <- data.frame(
        axis = c("x", "y", "z"),
        x = c(axes_origin[1] + axes_length,
              axes_origin[1],
              axes_origin[1]),
        y = c(axes_origin[2],
              axes_origin[2] + axes_length,
              axes_origin[2]),
        z = c(axes_origin[3],
              axes_origin[3],
              axes_origin[3] + axes_length)
      )
      
      p <- add_trace(
        p,
        data = axes_lab,
        x = ~x, y = ~y, z = ~z,
        type = "scatter3d",
        mode = "text",
        text = ~axis,
        textposition = "top center",
        textfont = list(color = "black", size = 12),
        showlegend = FALSE,
        hoverinfo = "none"
      )
    }
  }
  
  # ------------------------------------------------------------
  # Add each linkage separately
  # ------------------------------------------------------------
  for (lk in unique(segs$link)) {
    d <- segs[segs$link == lk, ]
    
    p <- add_trace(
      p,
      data = d,
      x = ~x, y = ~y, z = ~z,
      type = "scatter3d",
      mode = "lines",
      line = list(color = link_colors[lk], width = link_widths[lk]),
      showlegend = FALSE,
      hoverinfo = "none"
    )
  }
  
  p <- add_trace(
    p,
    data = pts,
    x = ~x, y = ~y, z = ~z,
    type = "scatter3d",
    mode = "markers+text",
    text = ~point,
    textposition = "top center",
    marker = list(size = 5, color = "black"),
    showlegend = FALSE
  )
  
  p %>% layout(
    title = title,
    scene = list(
      aspectmode = "data",
      
      xaxis = list(
        visible = FALSE,
        showgrid = FALSE,
        zeroline = FALSE,
        showbackground = FALSE,
        showticklabels = FALSE
      ),
      yaxis = list(
        visible = FALSE,
        showgrid = FALSE,
        zeroline = FALSE,
        showbackground = FALSE,
        showticklabels = FALSE
      ),
      zaxis = list(
        visible = FALSE,
        showgrid = FALSE,
        zeroline = FALSE,
        showbackground = FALSE,
        showticklabels = FALSE
      )
    )
  )
}


#' Plot a single linkage configuration using an alternate static Plotly layout.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @param link_map Named list specifying which point pairs form plotted links.
#' @param link_colors Named vector of colors for links.
#' @param link_widths Named vector of line widths for links.
#' @param remove_axes Logical; if TRUE, hide 3D axes.
#' @param fix_camera Logical; if TRUE, use a fixed camera.
#' @param camera Optional Plotly camera object.
#' @param pitch Camera pitch angle in degrees.
#' @param yaw Camera yaw angle in degrees.
#' @param roll Camera roll angle in degrees.
#' @param zoom Camera distance or scale factor.
#' @param show_link_labels Input argument used by this helper.
#' @param link_label_size Input argument used by this helper.
#' @param link_label_color Input argument used by this helper.
#' @param title Plot title.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
plot_linkage_static2 <- function(sol,
                                 link_map = NULL,
                                 link_colors = NULL,
                                 link_widths = NULL,
                                 remove_axes = TRUE,
                                 fix_camera = TRUE,
                                 camera = NULL,
                                 pitch = NULL,
                                 yaw = NULL,
                                 roll = NULL,
                                 zoom = 1.8,
                                 show_link_labels = TRUE,
                                 link_label_size = 10,
                                 link_label_color = "black",
                                 title = "Linkage") {
  
  segs <- segments_from_solution(sol, link_map = link_map)
  
  if (is.null(link_colors)) {
    link_colors <- setNames(rep("black", length(unique(segs$link))), unique(segs$link))
    defaults <- c(
      AB = "black", DC = "black",
      BM = "red",   CM = "red",
      BE = "blue",  CE = "blue",
      AD = "gray"
    )
    link_colors[names(defaults)[names(defaults) %in% names(link_colors)]] <-
      defaults[names(defaults) %in% names(link_colors)]
  }
  
  if (is.null(link_widths)) {
    link_widths <- setNames(rep(6, length(unique(segs$link))), unique(segs$link))
  }
  
  pts <- data.frame(
    point = c("A","D","B","C","E","M"),
    x = c(sol$A[1], sol$D[1], sol$B[1], sol$C[1], sol$E[1], sol$M[1]),
    y = c(sol$A[2], sol$D[2], sol$B[2], sol$C[2], sol$E[2], sol$M[2]),
    z = c(sol$A[3], sol$D[3], sol$B[3], sol$C[3], sol$E[3], sol$M[3])
  )
  
  p <- plot_ly()
  
  for (lk in unique(segs$link)) {
    d <- segs[segs$link == lk, ]
    
    p <- add_trace(
      p,
      data = d,
      x = ~x, y = ~y, z = ~z,
      type = "scatter3d",
      mode = "lines",
      line = list(
        color = unname(link_colors[[lk]] %||% "black"),
        width = unname(link_widths[[lk]] %||% 6)
      ),
      showlegend = FALSE,
      hoverinfo = "none"
    )
  }
  
  if (show_link_labels) {
    label_df <- do.call(rbind, lapply(unique(segs$link), function(lk) {
      d <- segs[segs$link == lk, ]
      p1 <- d[1, ]
      p2 <- d[2, ]
      
      data.frame(
        link = lk,
        x = mean(c(p1$x, p2$x)),
        y = mean(c(p1$y, p2$y)),
        z = mean(c(p1$z, p2$z))
      )
    }))
    
    p <- add_trace(
      p,
      data = label_df,
      x = ~x, y = ~y, z = ~z,
      type = "scatter3d",
      mode = "text",
      text = ~link,
      textposition = "middle center",
      textfont = list(size = link_label_size, color = link_label_color),
      showlegend = FALSE,
      hoverinfo = "none"
    )
  }
  
  p <- add_trace(
    p,
    data = pts,
    x = ~x, y = ~y, z = ~z,
    type = "scatter3d",
    mode = "markers+text",
    text = ~point,
    textposition = "top center",
    marker = list(size = 5, color = "black"),
    showlegend = FALSE
  )
  
  p %>%
    layout(
      title = title,
      scene = build_scene(
        bounds = NULL,
        remove_axes = remove_axes,
        fix_camera = fix_camera,
        camera = camera,
        pitch = pitch,
        yaw = yaw,
        roll = roll,
        zoom = zoom
      )
    )
}

#' Animate a sequence of linkage configurations in Plotly.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param solutions List of solved linkage states, one element per frame.
#' @param link_map Named list specifying which point pairs form plotted links.
#' @param link_colors Named vector of colors for links.
#' @param link_widths Named vector of line widths for links.
#' @param remove_axes Logical; if TRUE, hide 3D axes.
#' @param fix_camera Logical; if TRUE, use a fixed camera.
#' @param camera Optional Plotly camera object.
#' @param pitch Camera pitch angle in degrees.
#' @param yaw Camera yaw angle in degrees.
#' @param roll Camera roll angle in degrees.
#' @param zoom Camera distance or scale factor.
#' @param show_link_labels Input argument used by this helper.
#' @param link_label_size Input argument used by this helper.
#' @param link_label_color Input argument used by this helper.
#' @param title Plot title.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
animate_linkage_plotly2 <- function(solutions,
                                    link_map = NULL,
                                    link_colors = NULL,
                                    link_widths = NULL,
                                    remove_axes = TRUE,
                                    fix_camera = TRUE,
                                    camera = NULL,
                                    pitch = NULL,
                                    yaw = NULL,
                                    roll = NULL,
                                    zoom = 1.8,
                                    show_link_labels = TRUE,
                                    link_label_size = 10,
                                    link_label_color = "black",
                                    title = "Linkage animation") {
  
  segs <- segments_from_solutions(solutions, link_map = link_map)
  pts  <- solutions_to_df(solutions)
  bounds <- get_bounds(solutions)
  
  if (is.null(link_colors)) {
    link_colors <- setNames(rep("black", length(unique(segs$link))), unique(segs$link))
    defaults <- c(
      AB = "black", DC = "black",
      BM = "red",   CM = "red",
      BE = "blue",  CE = "blue",
      AD = "gray"
    )
    link_colors[names(defaults)[names(defaults) %in% names(link_colors)]] <-
      defaults[names(defaults) %in% names(link_colors)]
  }
  
  if (is.null(link_widths)) {
    link_widths <- setNames(rep(6, length(unique(segs$link))), unique(segs$link))
  }
  
  p <- plot_ly()
  
  for (lk in unique(segs$link)) {
    d <- segs[segs$link == lk, ]
    
    p <- add_trace(
      p,
      data = d,
      x = ~x, y = ~y, z = ~z,
      frame = ~frame,
      type = "scatter3d",
      mode = "lines",
      line = list(
        color = unname(link_colors[[lk]] %||% "black"),
        width = unname(link_widths[[lk]] %||% 6)
      ),
      showlegend = FALSE,
      hoverinfo = "none"
    )
  }
  
  if (show_link_labels) {
    label_df <- do.call(rbind, lapply(split(segs, segs$frame), function(dd) {
      do.call(rbind, lapply(unique(dd$link), function(lk) {
        d <- dd[dd$link == lk, ]
        p1 <- d[1, ]
        p2 <- d[2, ]
        
        data.frame(
          frame = p1$frame,
          link = lk,
          x = mean(c(p1$x, p2$x)),
          y = mean(c(p1$y, p2$y)),
          z = mean(c(p1$z, p2$z))
        )
      }))
    }))
    
    p <- add_trace(
      p,
      data = label_df,
      x = ~x, y = ~y, z = ~z,
      frame = ~frame,
      type = "scatter3d",
      mode = "text",
      text = ~link,
      textposition = "middle center",
      textfont = list(size = link_label_size, color = link_label_color),
      showlegend = FALSE,
      hoverinfo = "none"
    )
  }
  
  p <- add_trace(
    p,
    data = pts,
    x = ~x, y = ~y, z = ~z,
    frame = ~frame,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 4, color = "black"),
    showlegend = FALSE
  )
  
  p %>%
    layout(
      title = title,
      scene = build_scene(
        bounds = bounds,
        remove_axes = remove_axes,
        fix_camera = fix_camera,
        camera = camera,
        pitch = pitch,
        yaw = yaw,
        roll = roll,
        zoom = zoom
      ),
      updatemenus = list(
        list(
          type = "buttons",
          buttons = list(
            list(
              method = "animate",
              args = list(
                NULL,
                list(
                  frame = list(duration = 50, redraw = TRUE),
                  transition = list(duration = 0)
                )
              ),
              label = "Play"
            ),
            list(
              method = "animate",
              args = list(
                NULL,
                list(
                  frame = list(duration = 0, redraw = FALSE),
                  transition = list(duration = 0)
                )
              ),
              label = "Pause"
            )
          )
        )
      )
    )
}



# ============================================================
# Basic 3d to 2d  helpers
# ============================================================

#' Compute the Euclidean norm of a numeric vector.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param x Numeric vector, coordinate, or value used by the helper.
#' @return A single numeric norm.
#' @keywords internal
nrm <- function(x) sqrt(sum(x^2))

rot_z <- function(theta) {
  matrix(c(
    cos(theta), -sin(theta), 0,
    sin(theta),  cos(theta), 0,
    0,           0,          1
  ), 3, 3, byrow = TRUE)
}

#' Create a rotation matrix for rotation about the y axis.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param theta Rotation angle in radians.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
rot_y <- function(theta) {
  matrix(c(
    cos(theta), 0, sin(theta),
    0,          1, 0,
    -sin(theta), 0, cos(theta)
  ), 3, 3, byrow = TRUE)
}

#' Create a rotation matrix for rotation about the x axis.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param theta Rotation angle in radians.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
rot_x <- function(theta) {
  matrix(c(
    1, 0,           0,
    0, cos(theta), -sin(theta),
    0, sin(theta),  cos(theta)
  ), 3, 3, byrow = TRUE)
}

# ============================================================
# Build point and segment tables from one solution
# ============================================================


#' Convert one solver output object to a data frame of named 3D points.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
solution_points_df <- function(sol) {
  lapply(1:length(sol),function(x){
  data.frame(
    frame=x,
    point = c("A","D","B","C","E"),
    x = c(sol[[x]]$A[1], sol[[x]]$D[1], sol[[x]]$B[1], sol[[x]]$C[1], sol[[x]]$E[1]),
    y = c(sol[[x]]$A[2], sol[[x]]$D[2], sol[[x]]$B[2], sol[[x]]$C[2], sol[[x]]$E[2]),
    z = c(sol[[x]]$A[3], sol[[x]]$D[3], sol[[x]]$B[3], sol[[x]]$C[3], sol[[x]]$E[3])
  )}
)%>% bind_rows()
}

#' Convert one solver output object to a data frame of plotted segments.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
solution_segments_df <- function(sol) {
  pts <- list(
    A = sol$A, D = sol$D, B = sol$B,
    C = sol$C, E = sol$E, M = sol$M
  )
  
  links <- list(
    AB = c("A","B"),
    DC = c("D","C"),
    BM = c("B","M"),
    CM = c("C","M"),
    BE = c("B","E"),
    CE = c("C","E"),
    AD = c("A","D")
  )
  
  do.call(rbind, lapply(names(links), function(nm) {
    p1 <- pts[[links[[nm]][1]]]
    p2 <- pts[[links[[nm]][2]]]
    data.frame(
      link = nm,
      x1 = p1[1], y1 = p1[2], z1 = p1[3],
      x2 = p2[1], y2 = p2[2], z2 = p2[3]
    )
  }))
}

# ============================================================
# Camera projection
# yaw   = turn around vertical axis
# pitch = tilt up/down
# roll  = optional in-plane rotation
# ============================================================

#' Project 3D points into a rotated 2D camera view for ggplot-style figures.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param df Data frame of coordinates or segments.
#' @param yaw Camera yaw angle in degrees.
#' @param pitch Camera pitch angle in degrees.
#' @param roll Camera roll angle in degrees.
#' @param center Camera center as a length-3 numeric vector.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
project_camera <- function(df, yaw = 0, pitch = 0, roll = 0,
                           center = NULL) {
  X <- as.matrix(df[, c("x","y","z")])
  
  if (is.null(center)) {
    center <- colMeans(X)
  }
  
  Xc <- sweep(X, 2, center, "-")
  
  R <- rot_z(roll) %*% rot_x(pitch) %*% rot_y(yaw)
  Xr <- Xc %*% t(R)
  
  out <- df
  out$x2 <- Xr[, 1]
  out$y2 <- Xr[, 2]
  out$depth <- Xr[, 3]
  out
}

#' Project 3D line segments into a rotated 2D camera view.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param seg_df Data frame of line segments.
#' @param yaw Camera yaw angle in degrees.
#' @param pitch Camera pitch angle in degrees.
#' @param roll Camera roll angle in degrees.
#' @param center Camera center as a length-3 numeric vector.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
project_segments_camera <- function(seg_df, yaw = 0, pitch = 0, roll = 0,
                                    center = NULL) {
  p1 <- data.frame(x = seg_df$x1, y = seg_df$y1, z = seg_df$z1)
  p2 <- data.frame(x = seg_df$x2, y = seg_df$y2, z = seg_df$z2)
  
  if (is.null(center)) {
    all_xyz <- rbind(
      cbind(seg_df$x1, seg_df$y1, seg_df$z1),
      cbind(seg_df$x2, seg_df$y2, seg_df$z2)
    )
    center <- colMeans(all_xyz)
  }
  
  p1p <- project_camera(p1, yaw = yaw, pitch = pitch, roll = roll, center = center)
  p2p <- project_camera(p2, yaw = yaw, pitch = pitch, roll = roll, center = center)
  
  out <- seg_df
  out$x1p <- p1p$x2
  out$y1p <- p1p$y2
  out$z1p <- p1p$depth
  
  out$x2p <- p2p$x2
  out$y2p <- p2p$y2
  out$z2p <- p2p$depth
  
  out$depth <- 0.5 * (out$z1p + out$z2p)
  out
}

# ============================================================
# Default lateral camera for your linkage
# Tries to make E left, A/D right
# ============================================================

#' Estimate a lateral camera orientation from a solved linkage configuration.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
default_lateral_camera <- function(sol) {
  pts <- solution_points_df(sol)
  
  # center on whole linkage
  ctr <- colMeans(pts[, c("x","y","z")])
  
  # default side-ish view:
  # look mostly along +/- y so the y spread becomes depth
  # with a slight pitch so z is visible
  cam <- list(
    yaw = 0,
    pitch = -0.25,
    roll = 0,
    center = ctr
  )
  
  # project and flip horizontally if needed so E is leftmost
  pp <- project_camera(pts, yaw = cam$yaw, pitch = cam$pitch, roll = cam$roll, center = cam$center)
  
  xE  <- pp$x2[pp$point == "E"]
  xAD <- mean(pp$x2[pp$point %in% c("A","D")])
  
  if (xE > xAD) {
    cam$yaw <- cam$yaw + pi
  }
  
  cam
}

# ============================================================
# 2D plot
# Correct layer order: far to near
# ============================================================
#' Draw a linkage configuration as a 2D ggplot projection.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @param yaw Camera yaw angle in degrees.
#' @param pitch Camera pitch angle in degrees.
#' @param roll Camera roll angle in degrees.
#' @param center Camera center as a length-3 numeric vector.
#' @param use_default_camera Input argument used by this helper.
#' @param link_colors Named vector of colors for links.
#' @param link_sizes Input argument used by this helper.
#' @param link_alpha Input argument used by this helper.
#' @param point_size Input argument used by this helper.
#' @param label_size Input argument used by this helper.
#' @param show_labels Input argument used by this helper.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
gg_plot_linkage <- function(sol,
                            yaw = NULL,
                            pitch = NULL,
                            roll = 0,
                            center = NULL,
                            use_default_camera = TRUE,
                            link_colors = NULL,
                            link_sizes = NULL,
                            link_alpha = NULL,
                            point_size = 2.5,
                            label_size = 3.5,
                            show_labels = TRUE) {
  
  pts0 <- solution_points_df(sol)
  seg0 <- solution_segments_df(sol)
  
  if (use_default_camera) {
    cam <- default_lateral_camera(sol)
    yaw    <- yaw    %||% cam$yaw
    pitch  <- pitch  %||% cam$pitch
    center <- center %||% cam$center
  } else {
    center <- center %||% colMeans(pts0[, c("x", "y", "z")])
    yaw <- yaw %||% 0
    pitch <- pitch %||% 0
  }
  
  default_link_colors <- c(
    AB = "black",
    DC = "black",
    BM = "red3",
    CM = "red3",
    BE = "blue3",
    CE = "blue3",
    AD = "gray50"
  )
  
  default_link_sizes <- c(
    AB = 1,
    DC = 1,
    BM = 1.2,
    CM = 1.2,
    BE = 1,
    CE = 1,
    AD = 0.8
  )
  
  default_link_alpha <- c(
    AB = 1,
    DC = 1,
    BM = 1,
    CM = 1,
    BE = 1,
    CE = 1,
    AD = 1
  )
  
  pts <- project_camera(
    pts0,
    yaw = yaw, pitch = pitch, roll = roll,
    center = center
  )
  
  seg <- project_segments_camera(
    seg0,
    yaw = yaw, pitch = pitch, roll = roll,
    center = center
  )
  
  seg <- seg[order(seg$depth), ]
  pts <- pts[order(pts$depth), ]
  
  g <- ggplot()
  
  for (i in seq_len(nrow(seg))) {
    lk <- as.character(seg$link[i])
    
    col <- link_colors[[lk]] %||% default_link_colors[[lk]] %||% "black"
    sz  <- link_sizes[[lk]]  %||% default_link_sizes[[lk]]  %||% 1
    alp <- link_alpha[[lk]]  %||% default_link_alpha[[lk]]  %||% 1
    
    g <- g + geom_segment(
      data = seg[i, ],
      aes(x = x1p, y = y1p, xend = x2p, yend = y2p),
      linewidth = sz,
      color = scales::alpha(col, alp),
      lineend = "round"
    )
  }
  
  g <- g + geom_point(
    data = pts,
    aes(x = x2, y = y2),
    size = point_size
  )
  
  if (show_labels) {
    y_rng <- range(pts$y2, na.rm = TRUE)
    y_nudge <- 0.02 * diff(y_rng)
    if (!is.finite(y_nudge) || y_nudge == 0) y_nudge <- 0.02
    
    g <- g + geom_text(
      data = pts,
      aes(x = x2, y = y2, label = point),
      size = label_size,
      nudge_y = y_nudge
    )
  }
  
  g + coord_equal() + theme_void()
}
# ============================================================
# Easy viewpoint modifier
# ============================================================

#' Create left and right projected views of a solved linkage configuration.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @param side Input argument used by this helper.
#' @param tilt Input argument used by this helper.
#' @param spin Input argument used by this helper.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
view_left_right <- function(sol,
                            side = c("lateral", "front", "top"),
                            tilt = -0.25,
                            spin = 0) {
  side <- match.arg(side)
  ctr <- colMeans(solution_points_df(sol)[, c("x","y","z")])
  
  if (side == "lateral") {
    yaw <- 0 + spin
    pitch <- tilt
  } else if (side == "front") {
    yaw <- pi/2 + spin
    pitch <- tilt
  } else {
    yaw <- 0 + spin
    pitch <- -pi/2
  }
  
  list(yaw = yaw, pitch = pitch, roll = 0, center = ctr)
}

# ============================================================
# Save directly to PDF
# ============================================================

#' Save a projected linkage plot to a PDF file.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @param file Output file path.
#' @param yaw Camera yaw angle in degrees.
#' @param pitch Camera pitch angle in degrees.
#' @param roll Camera roll angle in degrees.
#' @param center Camera center as a length-3 numeric vector.
#' @param use_default_lateral Input argument used by this helper.
#' @param width Input argument used by this helper.
#' @param height Input argument used by this helper.
#' @param ... Additional arguments passed to downstream helpers.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
gg_save_linkage_pdf <- function(sol, file,
                                 yaw = NULL, pitch = NULL, roll = 0,
                                 center = NULL,
                                 use_default_lateral = TRUE,
                                 width = 6, height = 5,
                                 ...) {
  p <- gg_plot_linkage(
    sol,
    yaw = yaw, pitch = pitch, roll = roll,
    center = center,
    use_default_lateral = use_default_lateral,
    ...
  )
  
  ggsave(file, p, width = width, height = height, device = cairo_pdf)
}
# ============================================================
# Basic geometry helpers
# ============================================================

#' Compute a dot product for two vectors.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param a First vector, point, or object used in the calculation.
#' @param b Second vector, point, or object used in the calculation.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
dot3 <- function(a, b) sum(a * b)

dist_point_to_line_2d <- function(P, A, B) {
  AB <- B - A
  AP <- P - A
  abs(AB[1] * AP[2] - AB[2] * AP[1]) / sqrt(sum(AB^2))
}

#' Construct a local 2D basis from four 3D points defining a plane or plate.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param A 3D point A, typically a fixed or anatomical landmark.
#' @param B 3D point B, typically a linkage joint or angle vertex.
#' @param C 3D point C, typically a linkage joint.
#' @param D 3D point D, typically a driver or fixed landmark.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
plane_basis_from_ABCD <- function(A, B, C, D, tol = 1e-10) {
  e1 <- unit(D - A, tol)
  n  <- cross3(D - A, B - A)
  if (nrm(n) < tol) stop("ABCD does not define a plane.")
  e3 <- unit(n, tol)
  e2 <- unit(cross3(e3, e1), tol)
  
  list(
    origin = A,
    e1 = e1,
    e2 = e2,
    e3 = e3
  )
}

#' Project a 3D point into a supplied local 2D plane basis.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param P 3D point.
#' @param basis Plane basis object returned by a basis-construction helper.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
to_plane2d_basis <- function(P, basis) {
  v <- P - basis$origin
  c(dot3(v, basis$e1), dot3(v, basis$e2))
}

#' Convert local 2D plane coordinates back to 3D coordinates.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param xy Length-2 local plane coordinate.
#' @param basis Plane basis object returned by a basis-construction helper.
#' @param z_off Offset normal to the plane basis.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
from_plane2d_basis <- function(xy, basis, z_off = 0) {
  basis$origin + xy[1] * basis$e1 + xy[2] * basis$e2 + z_off * basis$e3
}

# ============================================================
# Big base inside ABCD
# Uses centroid of the quadrilateral in plane coordinates
# and radius = minimum distance from centroid to edges
# ============================================================

#' Estimate the larger base geometry used for frustum volume calculations from a solution.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
frustum_big_base_from_sol <- function(sol, tol = 1e-10) {
  A <- sol$A; B <- sol$B; C <- sol$C; D <- sol$D
  
  basis <- plane_basis_from_ABCD(A, B, C, D, tol = tol)
  
  A2 <- to_plane2d_basis(A, basis)
  B2 <- to_plane2d_basis(B, basis)
  C2 <- to_plane2d_basis(C, basis)
  D2 <- to_plane2d_basis(D, basis)
  
  quad2 <- rbind(A2, B2, C2, D2)
  
  # simple centroid of vertices; for your nearly symmetric plate this is fine
  G2 <- colMeans(quad2)
  
  # inscribed radius at that center
  dists <- c(
    dist_point_to_line_2d(G2, A2, B2),
    dist_point_to_line_2d(G2, B2, C2),
    dist_point_to_line_2d(G2, C2, D2),
    dist_point_to_line_2d(G2, D2, A2)
  )
  
  R <- min(dists)
  
  list(
    basis = basis,
    center2 = G2,
    center3 = from_plane2d_basis(G2, basis, z_off = 0),
    radius = R,
    quad2 = quad2
  )
}

# ============================================================
# Build frustum mesh
# r_small must be supplied (or you can set a default)
# ============================================================

#' Construct a mesh approximation of a frustum from a solved linkage configuration.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @param r_small Small-base radius used in frustum construction.
#' @param n Number of samples, frames, or mesh divisions, depending on context.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
make_frustum_mesh <- function(sol, r_small = 0.08, n = 60, tol = 1e-10) {
  fb <- frustum_big_base_from_sol(sol, tol = tol)
  
  G <- fb$center3
  E <- sol$E
  R <- fb$radius
  r <- r_small
  
  basis <- fb$basis
  e1 <- basis$e1
  e2 <- basis$e2
  e3 <- basis$e3
  
  # perpendicular height from E to plane ABCD
  h <- abs(dot3(E - G, e3))
  
  theta <- seq(0, 2 * pi, length.out = n + 1)[-(n + 1)]
  
  big <- t(sapply(theta, function(tt) {
    G + R * cos(tt) * e1 + R * sin(tt) * e2
  }))
  
  small <- t(sapply(theta, function(tt) {
    E + r * cos(tt) * e1 + r * sin(tt) * e2
  }))
  
  # mesh triangles
  x <- c(big[,1], small[,1])
  y <- c(big[,2], small[,2])
  z <- c(big[,3], small[,3])
  
  ii <- c(); jj <- c(); kk <- c()
  for (k in seq_len(n)) {
    k2 <- if (k < n) k + 1 else 1
    
    b1 <- k
    b2 <- k2
    s1 <- n + k
    s2 <- n + k2
    
    # two triangles per quad strip
    ii <- c(ii, b1, b2)
    jj <- c(jj, b2, s2)
    kk <- c(kk, s2, s1)
  }
  
  list(
    big = big,
    small = small,
    mesh = list(x = x, y = y, z = z, i = ii - 1, j = jj - 1, k = kk - 1),
    center_big = G,
    center_small = E,
    R = R,
    r = r,
    h = h,
    volume = pi * h * (R^2 + R * r + r^2) / 3
  )
}

# ============================================================
# Volume through time
# ============================================================

#' Estimate frustum volume across a sequence of linkage solutions.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param solutions List of solved linkage states, one element per frame.
#' @param r_small Small-base radius used in frustum construction.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
frustum_volume_path <- function(solutions, r_small = 0.08, tol = 1e-10) {
  do.call(rbind, lapply(seq_along(solutions), function(i) {
    fr <- make_frustum_mesh(solutions[[i]], r_small = r_small, tol = tol)
    data.frame(
      frame = i,
      R = fr$R,
      r = fr$r,
      h = fr$h,
      volume = fr$volume,
      Ex = solutions[[i]]$E[1],
      Ey = solutions[[i]]$E[2],
      Ez = solutions[[i]]$E[3],
      Mx = solutions[[i]]$M[1],
      My = solutions[[i]]$M[2],
      Mz = solutions[[i]]$M[3]
    )
  }))
}

#' Plot a static linkage configuration with the frustum mesh overlaid.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @param r_small Small-base radius used in frustum construction.
#' @param link_colors Named vector of colors for links.
#' @param link_widths Named vector of line widths for links.
#' @param remove_axes Logical; if TRUE, hide 3D axes.
#' @param fix_camera Logical; if TRUE, use a fixed camera.
#' @param camera Optional Plotly camera object.
#' @param frustum_color Input argument used by this helper.
#' @param rim_color Input argument used by this helper.
#' @param title Plot title.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
plot_linkage_static_frustum <- function(sol,
                                        r_small = 0.08,
                                        link_colors = NULL,
                                        link_widths = NULL,
                                        remove_axes = TRUE,
                                        fix_camera = TRUE,
                                        camera = NULL,
                                        frustum_color = "rgba(100,150,255,0.35)",
                                        rim_color = "rgba(40,80,180,0.9)",
                                        title = "Linkage + frustum") {
  
  segs <- segments_from_solution(sol)
  
  if (is.null(link_colors)) {
    link_colors <- c(
      AB="black", DC="black",
      BM="red",   CM="red",
      BE="blue",  CE="blue",
      AD="gray"
    )
  }
  
  if (is.null(link_widths)) {
    link_widths <- setNames(rep(6, length(link_colors)), names(link_colors))
  }
  
  pts <- data.frame(
    point = c("A","D","B","C","E","M"),
    x = c(sol$A[1], sol$D[1], sol$B[1], sol$C[1], sol$E[1], sol$M[1]),
    y = c(sol$A[2], sol$D[2], sol$B[2], sol$C[2], sol$E[2], sol$M[2]),
    z = c(sol$A[3], sol$D[3], sol$B[3], sol$C[3], sol$E[3], sol$M[3])
  )
  
  fr <- make_frustum_mesh(sol, r_small = r_small)
  
  p <- plot_ly()
  
  # frustum surface
  p <- add_trace(
    p,
    type = "mesh3d",
    x = fr$mesh$x,
    y = fr$mesh$y,
    z = fr$mesh$z,
    i = fr$mesh$i,
    j = fr$mesh$j,
    k = fr$mesh$k,
    color = frustum_color,
    opacity = 1,
    hoverinfo = "skip",
    showlegend = FALSE
  )
  
  # rims
  rim_df <- rbind(
    data.frame(group = "big",   x = c(fr$big[,1],   fr$big[1,1]),
               y = c(fr$big[,2],   fr$big[1,2]),
               z = c(fr$big[,3],   fr$big[1,3])),
    data.frame(group = "small", x = c(fr$small[,1], fr$small[1,1]),
               y = c(fr$small[,2], fr$small[1,2]),
               z = c(fr$small[,3], fr$small[1,3]))
  )
  
  for (gg in unique(rim_df$group)) {
    d <- rim_df[rim_df$group == gg, ]
    p <- add_trace(
      p,
      data = d,
      x = ~x, y = ~y, z = ~z,
      type = "scatter3d",
      mode = "lines",
      line = list(color = rim_color, width = 5),
      showlegend = FALSE,
      hoverinfo = "none"
    )
  }
  
  # linkage
  for (lk in unique(segs$link)) {
    d <- segs[segs$link == lk, ]
    p <- add_trace(
      p,
      data = d,
      x = ~x, y = ~y, z = ~z,
      type = "scatter3d",
      mode = "lines",
      line = list(color = link_colors[lk], width = link_widths[lk]),
      showlegend = FALSE,
      hoverinfo = "none"
    )
  }
  
  p <- add_trace(
    p,
    data = pts,
    x = ~x, y = ~y, z = ~z,
    type = "scatter3d",
    mode = "markers+text",
    text = ~point,
    textposition = "top center",
    marker = list(size = 5, color = "black"),
    showlegend = FALSE
  )
  
  p %>%
    layout(
      title = paste0(title, "  |  V = ", signif(fr$volume, 4)),
      scene = build_scene(
        bounds = NULL,
        remove_axes = remove_axes,
        fix_camera = fix_camera,
        camera = camera
      )
    )
}

#find angles



#' Compute the angle formed by three named points in one solver output.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @param pts Named list or vector of point definitions.
#' @param degrees Logical; if TRUE, return angles in degrees.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
angle_points <- function(sol, pts =c("B","E","B_ref"),degrees = TRUE, tol = 1e-12) {
  A <- sol[[pts[1]]]
  B <- sol[[pts[2]]]
  C <- sol[[pts[3]]]
  
  v1 <- A - B
  v2 <- C - B
  
  n1 <- sqrt(sum(v1^2))
  n2 <- sqrt(sum(v2^2))
  
  if (n1 < tol || n2 < tol) {
    stop("One of the vectors BA or BC has near-zero length.")
  }
  
  cosang <- sum(v1 * v2) / (n1 * n2)
  cosang <- max(-1, min(1, cosang))
  
  ang <- acos(cosang)
  
  if (degrees) ang <- ang * 180 / pi
  
  ang
}


#' Compute a three-point angle over a list of solver outputs.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param solutions List of solved linkage states, one element per frame.
#' @param degrees Logical; if TRUE, return angles in degrees.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @param ... Additional arguments passed to downstream helpers.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
angles_points <- function(solutions, degrees = TRUE, tol = 1e-12,...) {
  data.frame(
    frame = seq_along(solutions),
    angle = sapply(solutions, function(x) angle_points(x, degrees = degrees, tol = tol,...))
  )
}



#' Recast landmark coordinates into a reference plane and optionally generate reflected points.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param dat Input data frame containing coordinate data.
#' @param pts Named list or vector of point definitions.
#' @param plane Plane definition, usually named points or coordinate axes.
#' @param origin Origin definition for recasting coordinates.
#' @param axis_pair Pair of landmarks used to align one plane axis.
#' @param axis Axis to align when recasting coordinates.
#' @param plane_axes Names of output axes used to span the recast plane.
#' @param reflect Input argument used by this helper.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
recast_to_plane_named <- function(dat,
                                  pts,
                                  plane,
                                  origin = c("first", "centroid"),
                                  axis_pair = NULL,
                                  axis = c("y", "x", "z"),
                                  plane_axes = c("x", "y"),
                                  reflect = TRUE,
                                  tol = 1e-10) {
  
  origin <- match.arg(origin)
  axis   <- match.arg(axis)
  
  # validate plane_axes manually
  if (!is.character(plane_axes) || length(plane_axes) != 2) {
    stop("plane_axes must be a character vector of length 2, e.g. c('y','z')")
  }
  if (!all(plane_axes %in% c("x", "y", "z"))) {
    stop("plane_axes must contain only 'x', 'y', and/or 'z'")
  }
  if (length(unique(plane_axes)) != 2) {
    stop("plane_axes must contain two different axes")
  }
  
  all_axes <- c("x", "y", "z")
  normal_axis <- setdiff(all_axes, plane_axes)
  
  dat <- as.data.frame(dat)
  if (!all(c("x","y","z") %in% names(dat))) {
    stop("dat must have columns x, y, z")
  }
  
  X <- as.matrix(dat[, c("x","y","z")])
  storage.mode(X) <- "double"
  
  nrm <- function(v) sqrt(sum(v^2))
  
  unit <- function(v) {
    d <- nrm(v)
    if (d < tol) stop("Zero-length vector")
    v / d
  }
  
  cross3 <- function(a, b) {
    c(
      a[2]*b[3] - a[3]*b[2],
      a[3]*b[1] - a[1]*b[3],
      a[1]*b[2] - a[2]*b[1]
    )
  }
  
  get_pt <- function(name) {
    if (!name %in% names(pts)) stop("Point name not found in pts: ", name)
    X[pts[[name]], ]
  }
  
  plane_names <- names(plane)
  if (length(plane_names) != 3) {
    stop("plane must name exactly 3 points")
  }
  
  pt_idx <- unlist(pts, use.names = FALSE)
  
  P1 <- get_pt(plane_names[1])
  P2 <- get_pt(plane_names[2])
  P3 <- get_pt(plane_names[3])
  
  O <- switch(
    origin,
    first    = P1,
    centroid = colMeans(rbind(P1, P2, P3))
  )
  
  u <- unit(P2 - P1)
  n <- unit(cross3(P2 - P1, P3 - P1))
  v <- unit(cross3(n, u))
  
  if (!is.null(axis_pair)) {
    nm <- names(axis_pair)
    if (length(nm) != 2) {
      stop("axis_pair must have exactly 2 named points")
    }
    
    pair_vec <- get_pt(nm[2]) - get_pt(nm[1])
    
    if (axis %in% plane_axes) {
      pair_proj <- pair_vec - sum(pair_vec * n) * n
      if (nrm(pair_proj) < tol) {
        stop("Bad axis_pair: projected vector in plane is too small")
      }
      
      pair_proj <- unit(pair_proj)
      
      if (axis == plane_axes[1]) {
        u <- pair_proj
        v <- unit(cross3(n, u))
      } else if (axis == plane_axes[2]) {
        v <- pair_proj
        u <- unit(cross3(v, n))
      }
      
    } else if (axis == normal_axis) {
      s <- sum(pair_vec * n)
      if (abs(s) < tol) {
        stop("Bad axis_pair for normal axis: point-pair is nearly perpendicular to plane normal")
      }
      if (s < 0) {
        n <- -n
        v <- unit(cross3(n, u))
      }
    }
  }
  
  basis <- list(x = NULL, y = NULL, z = NULL)
  basis[[plane_axes[1]]] <- u
  basis[[plane_axes[2]]] <- v
  basis[[normal_axis]]   <- n
  
  R <- rbind(
    basis[["x"]],
    basis[["y"]],
    basis[["z"]]
  )
  
  Xc <- sweep(X, 2, O, "-")
  Xnew <- Xc %*% t(R)
  
  coords <- data.frame(
    point = names(pts),
    name  = names(pts),
    x     = Xnew[pt_idx, 1],
    y     = Xnew[pt_idx, 2],
    z     = Xnew[pt_idx, 3],
    side  = "original",
    stringsAsFactors = FALSE
  )
  
  coords[[normal_axis]][coords$point %in% plane_names] <- 0
  
  reflected <- NULL
  if (reflect) {
    reflected <- coords[!(coords$point %in% plane_names), , drop = FALSE]
    reflected[[normal_axis]] <- -reflected[[normal_axis]]
    reflected$side <- "reflected"
    reflected$name <- paste0(reflected$point, "_ref")
  }
  
  all_coords <- rbind(coords, reflected)
  rownames(all_coords) <- NULL
  
  list(
    all_coords = all_coords,
    coords = coords,
    reflected_coords = reflected,
    rotation = R,
    origin = O,
    axes = basis,
    plane = plane_names,
    plane_axes = plane_axes,
    normal_axis = normal_axis
  )
}
#' Plot original and reflected coordinates with optional linkage segments in 3D.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param all_coords Coordinate data frame containing original/reflected/merged points.
#' @param link_map Named list specifying which point pairs form plotted links.
#' @param pitch Camera pitch angle in degrees.
#' @param yaw Camera yaw angle in degrees.
#' @param zoom Camera distance or scale factor.
#' @param show_labels Input argument used by this helper.
#' @param point_size Input argument used by this helper.
#' @param line_width Input argument used by this helper.
#' @return A Plotly figure.
#' @keywords internal
plot_symmetry_3d <- function(all_coords,
                             link_map = NULL,
                             pitch = 15,
                             yaw = 110,
                             zoom = 2,
                             show_labels = TRUE,
                             point_size = 6,
                             line_width = 5) {
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required.")
  }
  
  dat <- as.data.frame(all_coords)
  
  req_cols <- c("name", "x", "y", "z", "side")
  if (!all(req_cols %in% names(dat))) {
    stop("all_coords must contain columns: ", paste(req_cols, collapse = ", "))
  }
  
  dat$name <- as.character(dat$name)
  dat$side <- as.character(dat$side)
  
  # -----------------------------
  # helper: make link dataframe
  # -----------------------------
  make_link_df <- function(link_map, dat) {
    if (is.null(link_map) || length(link_map) == 0) return(NULL)
    
    out <- lapply(names(link_map), function(link_nm) {
      pair <- link_map[[link_nm]]
      
      if (length(pair) != 2) {
        stop("Each element of link_map must be a character vector of length 2.")
      }
      
      i1 <- match(pair[1], dat$name)
      i2 <- match(pair[2], dat$name)
      
      if (any(is.na(c(i1, i2)))) {
        missing_pts <- pair[is.na(c(i1, i2))]
        stop("These points in link_map were not found in all_coords$name: ",
             paste(missing_pts, collapse = ", "))
      }
      
      data.frame(
        link = link_nm,
        x = c(dat$x[i1], dat$x[i2], NA_real_),
        y = c(dat$y[i1], dat$y[i2], NA_real_),
        z = c(dat$z[i1], dat$z[i2], NA_real_),
        stringsAsFactors = FALSE
      )
    })
    
    do.call(rbind, out)
  }
  
  link_df <- make_link_df(link_map, dat)
  
  # -----------------------------
  # camera from pitch/yaw/zoom
  # -----------------------------
  deg2rad <- function(x) x * pi / 180
  
  pitch_r <- deg2rad(pitch)
  yaw_r   <- deg2rad(yaw)
  
  eye <- list(
    x = zoom * cos(pitch_r) * cos(yaw_r),
    y = zoom * cos(pitch_r) * sin(yaw_r),
    z = zoom * sin(pitch_r)
  )
  
  # -----------------------------
  # build plot
  # -----------------------------
  p <- plotly::plot_ly()
  
  # points
  sides <- unique(dat$side)
  for (s in sides) {
    dsub <- dat[dat$side == s, , drop = FALSE]
    
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
  
  # links
  if (!is.null(link_df) && nrow(link_df) > 0) {
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
  }
  
  # layout
  p <- plotly::layout(
    p,
    scene = list(
      xaxis = list(title = "X"),
      yaxis = list(title = "Y"),
      zaxis = list(title = "Z"),
      aspectmode = "data",
      camera = list(eye = eye)
    ),
    legend = list(orientation = "h")
  )
  
  p
}


#' Merge selected reflected/original coordinate pairs into midline coordinates.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param coords Data frame or tibble containing named coordinates with columns name, x, y, and z.
#' @param merge_pts Character vector of paired points to merge to the midline.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
merge_coords <- function(coords=NULL,merge_pts=c("E","D")){
  
  d_unmerge <- coords %>% filter(!point%in%merge_pts)
  d_merge <- coords %>% filter(point%in%merge_pts) %>% 
    group_by(point) %>% 
    mutate(across(c(x, y, z), mean)) %>% 
    mutate(side="merged",name=point) %>% 
    distinct()
  
  rbind(d_unmerge,d_merge) %>% arrange(point,name)
}


# distance from point p to plane through a, b, c
#' Compute signed perpendicular distance from a point to a 3D plane.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param a First vector, point, or object used in the calculation.
#' @param b Second vector, point, or object used in the calculation.
#' @param c Third point defining a plane or angle.
#' @param p Point to evaluate.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
point_plane_dist_3d <- function(a, b, c, p) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  c <- as.numeric(c)
  p <- as.numeric(p)
  
  ab <- b - a
  ac <- c - a
  
  n <- c(
    ab[2] * ac[3] - ab[3] * ac[2],
    ab[3] * ac[1] - ab[1] * ac[3],
    ab[1] * ac[2] - ab[2] * ac[1]
  )
  
  denom <- sqrt(sum(n^2))
  if (denom == 0) return(NA_real_)
  
  abs(sum(n * (p - a))) / denom
}

#' Compute the distance between two named points over a sequence of solutions.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param solutions List of solved linkage states, one element per frame.
#' @param p1 First point name or coordinate.
#' @param p2 Second point name or coordinate.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
dist_series <- function(solutions, p1, p2) {
  
  nrm <- function(x) sqrt(sum(x^2))
  
  sapply(solutions, function(sol) {
    nrm(sol[[p1]] - sol[[p2]])
  })
}

#' Compute point-to-plane distance over a sequence of solutions.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param solutions List of solved linkage states, one element per frame.
#' @param point Point name to evaluate or transform.
#' @param plane_pts Three points defining a plane.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
point_plane_series <- function(solutions, point, plane_pts) {
  
  resolve_pt <- function(sol, p) {
    if (is.character(p) && length(p) == 1) {
      if (!p %in% names(sol)) {
        stop("Point name not found in solution: ", p)
      }
      return(sol[[p]])
    }
    
    if (is.numeric(p) && length(p) == 3) {
      return(as.numeric(p))
    }
    
    stop("Each point must be either a single character name or a numeric vector of length 3.")
  }
  
  if (length(plane_pts) != 3) {
    stop("plane_pts must contain exactly 3 points.")
  }
  
  sapply(solutions, function(sol) {
    P <- resolve_pt(sol, point)
    A <- resolve_pt(sol, plane_pts[[1]])
    B <- resolve_pt(sol, plane_pts[[2]])
    C <- resolve_pt(sol, plane_pts[[3]])
    
    point_plane_dist_3d(A, B, C, P)
  })
}

#' Generate a smooth logistic-style trajectory profile.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param n Number of samples, frames, or mesh divisions, depending on context.
#' @param min_val Minimum value of the generated profile.
#' @param max_val Maximum value of the generated profile.
#' @param k Steepness of the logistic profile.
#' @param mid Relative position of the profile midpoint.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
track_profile <- function(n,min_val=0,max_val, k = 10, mid = 0.5) {
  t <- seq(min_val, 1, length.out = n)
  
  raw <- 1 / (1 + exp(-k * (t - mid)))
  
  # rescale so it starts at 0 and ends at x_total
  prof <- (raw - min(raw)) / (max(raw) - min(raw))
  max_val * prof
}


#' Generate the trajectory of one model point rotated around another point in a chosen plane.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param model Linkage model object produced by the corresponding model-construction function.
#' @param point Point name to evaluate or transform.
#' @param center Camera center as a length-3 numeric vector.
#' @param angles_deg Vector of rotation angles in degrees.
#' @param plane Plane definition, usually named points or coordinate axes.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
point_trajectory <- function(model,
                                           point,
                                           center,
                                           angles_deg,
                                           plane = c("xz", "xy", "yz")) {
  plane <- match.arg(plane)
  
  if (is.null(model$pts0)) {
    stop("model must contain a pts0 element with named 3D points.")
  }
  
  if (!point %in% names(model$pts0)) {
    stop("Point not found in model$pts0: ", point)
  }
  
  if (!center %in% names(model$pts0)) {
    stop("Center point not found in model$pts0: ", center)
  }
  
  P <- model$pts0[[point]]
  C <- model$pts0[[center]]
  
  idx <- switch(
    plane,
    xy = c(1, 2),
    xz = c(1, 3),
    yz = c(2, 3)
  )
  
  rotate_one <- function(angle_deg) {
    th <- angle_deg * pi / 180
    
    p2 <- P[idx]
    c2 <- C[idx]
    
    d1 <- p2[1] - c2[1]
    d2 <- p2[2] - c2[2]
    
    c(
      x = c2[1] + d1 * cos(th) - d2 * sin(th),
      y = c2[2] + d1 * sin(th) + d2 * cos(th)
    )
  }
  
  out <- t(sapply(angles_deg, rotate_one))
  
  data.frame(
    step = seq_along(angles_deg),
    angle_deg = angles_deg,
    x = out[, 1],
    y = out[, 2]
  )
}

## ------------------------------------------------------------
## Centroid size helpers for selected 3D points in a solution list
## ------------------------------------------------------------

#' Compute centroid size for a set of coordinates.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param coords Data frame or tibble containing named coordinates with columns name, x, y, and z.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
centroid_size <- function(coords) {
  coords <- as.matrix(coords)
  
  if (!is.numeric(coords)) {
    stop("coords must be numeric.")
  }
  
  if (ncol(coords) != 3) {
    stop("coords must have 3 columns: x, y, z.")
  }
  
  ctr <- colMeans(coords)
  
  sqrt(sum(rowSums((coords - matrix(ctr, nrow(coords), 3, byrow = TRUE))^2)))
}


#' Extract selected named points from one solver output as a coordinate matrix.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @param points Character vector of point names to extract or use.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
get_points_matrix <- function(sol, points) {
  missing_pts <- setdiff(points, names(sol))
  
  if (length(missing_pts) > 0) {
    stop(
      "These points are missing from the solution: ",
      paste(missing_pts, collapse = ", ")
    )
  }
  
  coords <- do.call(rbind, sol[points])
  
  coords <- as.data.frame(coords)
  names(coords) <- c("x", "y", "z")
  coords$point <- rownames(coords)
  
  coords[, c("point", "x", "y", "z")]
}


#' Compute centroid size for selected points in one solver output.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @param points Character vector of point names to extract or use.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
centroid_size_from_solution <- function(sol, points) {
  coords <- get_points_matrix(sol, points)
  
  centroid_size(coords[, c("x", "y", "z")])
}

#' Compute centroid size for selected points across a sequence of solver outputs.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param solutions List of solved linkage states, one element per frame.
#' @param points Character vector of point names to extract or use.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
centroid_size_over_solutions <- function(solutions, points) {
  out <- lapply(seq_along(solutions), function(i) {
    cs <- centroid_size_from_solution(solutions[[i]], points)
    
    data.frame(
      frame = i,
      centroid_size = cs
    )
  })
  
  do.call(rbind, out)
}

