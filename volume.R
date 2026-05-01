# ============================================================
# Volume-estimation helpers for linkage solutions
#
# Functions in this file resolve named and derived hull vertices,
# estimate convex-hull volumes for buccal compartments, and summarize
# volume changes along the anteroposterior axis. Function comments use
# roxygen-style tags so the file can be converted into package
# documentation later if desired.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(plotly)
  library(geometry)
  library(purrr)
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
#' @return A normalized numeric vector.
#' @keywords internal
unit <- function(x, tol = 1e-12) {
  d <- nrm(x)
  if (d < tol) stop("Cannot normalize near-zero vector.")
  x / d
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
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1]
  )
}

#' Remove duplicate rows from a coordinate matrix using a numerical tolerance.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param mat Numeric matrix.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
unique_rows_tol <- function(mat, tol = 1e-10) {
  mat <- as.matrix(mat)
  key <- apply(round(mat / tol) * tol, 1, paste, collapse = "_")
  mat[!duplicated(key), , drop = FALSE]
}

#' Extract a named 3D point from one solver output object.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @param nm Name of a point to extract.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
get_sol_pt <- function(sol, nm) {
  if (!nm %in% names(sol)) stop("Point not found in solution: ", nm)
  p <- sol[[nm]]
  if (!is.numeric(p) || length(p) != 3) stop("Element is not a 3D point: ", nm)
  as.numeric(p)
}

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

# ============================================================
# Plane helpers
# ============================================================

#' Define a 3D plane from three non-collinear points.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param p1 First point name or coordinate.
#' @param p2 Second point name or coordinate.
#' @param p3 Input argument used by this helper.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
plane_from_3pts <- function(p1, p2, p3, tol = 1e-12) {
  n <- cross3(p2 - p1, p3 - p1)
  if (nrm(n) < tol) stop("Degenerate plane definition.")
  n <- unit(n, tol = tol)
  d <- -sum(n * p1)
  list(normal = n, d = d)
}

#' Evaluate the signed plane equation value for a point.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param p Point to evaluate.
#' @param plane Plane definition, usually named points or coordinate axes.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
plane_value <- function(p, plane) {
  sum(plane$normal * p) + plane$d
}

#' Find the intersection of a line segment with a plane.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param p1 First point name or coordinate.
#' @param p2 Second point name or coordinate.
#' @param plane Plane definition, usually named points or coordinate axes.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
line_plane_intersection_segment <- function(p1, p2, plane, tol = 1e-10) {
  v1 <- plane_value(p1, plane)
  v2 <- plane_value(p2, plane)
  
  if (abs(v1) < tol) return(list(point = p1, t = 0, on_segment = TRUE))
  if (abs(v2) < tol) return(list(point = p2, t = 1, on_segment = TRUE))
  if (abs(v1 - v2) < tol) stop("Segment is parallel to plane or lies in plane.")
  
  t <- v1 / (v1 - v2)
  p <- p1 + t * (p2 - p1)
  
  list(
    point = p,
    t = t,
    on_segment = (t >= -tol && t <= 1 + tol)
  )
}

# ============================================================
# Solver-output helpers
# ============================================================

#' Convert one solver output object to a data frame of named 3D points.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
solution_points_df <- function(sol) {
  nm_pts <- names(sol)[vapply(sol, function(x) is.numeric(x) && length(x) == 3, logical(1))]
  
  bind_rows(lapply(nm_pts, function(nm) {
    P <- sol[[nm]]
    tibble(
      point = sub("_ref$", "", nm),
      name = nm,
      x = P[1],
      y = P[2],
      z = P[3],
      side = infer_side(nm, P[2])
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
segments_from_solution <- function(sol, link_map) {
  bind_rows(lapply(names(link_map), function(nm) {
    pts <- link_map[[nm]]
    if (length(pts) != 2) return(NULL)
    if (!all(pts %in% names(sol))) return(NULL)
    if (!is.numeric(sol[[pts[1]]]) || length(sol[[pts[1]]]) != 3) return(NULL)
    if (!is.numeric(sol[[pts[2]]]) || length(sol[[pts[2]]]) != 3) return(NULL)
    
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

# ============================================================
# General volume-spec system
# ============================================================

# A volume spec is a list whose elements are either:
# 1. a character string naming an existing point in sol, e.g. "B"
# 2. a list describing a derived point from a line-plane intersection:
#    list(
#      type  = "line_plane_intersection",
#      name  = "M",
#      line  = c("A", "C"),
#      plane = c("L", "G", "L_ref")
#    )
#
# Example default specs matching your original functions:
#' Return the default point specification for volume compartment 1.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
default_vol1_spec <- function() {
  list(
    list(
      type  = "line_plane_intersection",
      name  = "M",
      line  = c("A", "C"),
      plane = c("L", "G", "L_ref")
    ),
    list(
      type  = "line_plane_intersection",
      name  = "M_ref",
      line  = c("A_ref", "C_ref"),
      plane = c("L", "G", "L_ref")
    ),
    "C",
    "C_ref",
    "B",
    "B_ref",
    "D",
    "L",
    "L_ref"
  )
}

#' Return the default point specification for volume compartment 2.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
default_vol2_spec <- function() {
  list(
    "B",
    "B_ref",
    "D",
    "K",
    "E",
    "L",
    "L_ref"
  )
}

#' Resolve direct and derived point specifications into hull vertices.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @param spec Volume specification listing direct and derived hull points.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
resolve_volume_spec <- function(sol, spec, tol = 1e-10) {
  if (is.null(spec) || length(spec) == 0) {
    stop("Volume spec is empty.")
  }
  
  pts <- vector("list", length(spec))
  nm_out <- character(length(spec))
  
  for (i in seq_along(spec)) {
    item <- spec[[i]]
    
    # direct point name
    if (is.character(item) && length(item) == 1) {
      nm_out[i] <- item
      pts[[i]] <- get_sol_pt(sol, item)
      next
    }
    
    # derived point
    if (is.list(item) && !is.null(item$type)) {
      if (identical(item$type, "line_plane_intersection")) {
        if (is.null(item$name)) stop("Derived spec is missing `name`.")
        if (is.null(item$line) || length(item$line) != 2) {
          stop("line_plane_intersection spec must have `line = c(pt1, pt2)`.")
        }
        if (is.null(item$plane) || length(item$plane) != 3) {
          stop("line_plane_intersection spec must have `plane = c(p1, p2, p3)`.")
        }
        
        p1 <- get_sol_pt(sol, item$line[1])
        p2 <- get_sol_pt(sol, item$line[2])
        q1 <- get_sol_pt(sol, item$plane[1])
        q2 <- get_sol_pt(sol, item$plane[2])
        q3 <- get_sol_pt(sol, item$plane[3])
        
        pl <- plane_from_3pts(q1, q2, q3, tol = tol)
        int <- line_plane_intersection_segment(p1, p2, pl, tol = tol)
        
        if (!isTRUE(int$on_segment)) {
          stop(
            "Plane ", paste(item$plane, collapse = "-"),
            " does not intersect segment ",
            paste(item$line, collapse = "-"),
            " within the segment."
          )
        }
        
        nm_out[i] <- item$name
        pts[[i]] <- int$point
        next
      }
      
      stop("Unknown derived point type: ", item$type)
    }
    
    stop("Each volume spec element must be a point name or a valid derived-point spec.")
  }
  
  mat <- do.call(rbind, pts)
  colnames(mat) <- c("x", "y", "z")
  rownames(mat) <- nm_out
  unique_rows_tol(mat, tol = tol)
}

# Backward-compatible wrappers
#' Return hull vertices for volume compartment 1.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @param spec Volume specification listing direct and derived hull points.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
get_vol1_vertices <- function(sol, tol = 1e-10, spec = default_vol1_spec()) {
  resolve_volume_spec(sol, spec = spec, tol = tol)
}

#' Return hull vertices for volume compartment 2.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @param spec Volume specification listing direct and derived hull points.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
get_vol2_vertices <- function(sol, tol = 1e-10, spec = default_vol2_spec()) {
  resolve_volume_spec(sol, spec = spec, tol = tol)
}

# General getter by volume name
#' Return hull vertices for a named volume compartment.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @param volume Volume name to retrieve, such as vol1 or vol2.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @param vol1_spec Point specification for volume compartment 1.
#' @param vol2_spec Point specification for volume compartment 2.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
get_volume_vertices <- function(sol,
                                volume = c("vol1", "vol2"),
                                tol = 1e-10,
                                vol1_spec = default_vol1_spec(),
                                vol2_spec = default_vol2_spec()) {
  volume <- match.arg(volume)
  if (volume == "vol1") return(get_vol1_vertices(sol, tol = tol, spec = vol1_spec))
  if (volume == "vol2") return(get_vol2_vertices(sol, tol = tol, spec = vol2_spec))
}

# ============================================================
# Volume estimation
# ============================================================

#' Estimate the volume enclosed by the convex hull of a point cloud.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param pts Named list or vector of point definitions.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
convex_hull_volume <- function(pts, tol = 1e-10) {
  pts <- unique_rows_tol(pts, tol = tol)
  if (nrow(pts) < 4) return(NA_real_)
  geometry::convhulln(pts, options = "FA")$vol
}

#' Estimate volume compartment 1 from one solver output.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @param spec Volume specification listing direct and derived hull points.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
estimate_vol1_from_solution <- function(sol,
                                        tol = 1e-10,
                                        spec = default_vol1_spec()) {
  pts <- get_vol1_vertices(sol, tol = tol, spec = spec)
  list(volume = convex_hull_volume(pts, tol = tol), hull_points = pts)
}

#' Estimate volume compartment 2 from one solver output.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @param spec Volume specification listing direct and derived hull points.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
estimate_vol2_from_solution <- function(sol,
                                        tol = 1e-10,
                                        spec = default_vol2_spec()) {
  pts <- get_vol2_vertices(sol, tol = tol, spec = spec)
  list(volume = convex_hull_volume(pts, tol = tol), hull_points = pts)
}

#' Estimate volume compartments across an entire solution path.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param solutions List of solved linkage states, one element per frame.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @param quiet Logical; if TRUE, suppress per-frame volume failure messages.
#' @param vol1_spec Point specification for volume compartment 1.
#' @param vol2_spec Point specification for volume compartment 2.
#' @return A tibble of per-frame volume estimates.
#' @keywords internal
estimate_meshvol_over_path <- function(solutions,
                                       tol = 1e-10,
                                       quiet = FALSE,
                                       vol1_spec = default_vol1_spec(),
                                       vol2_spec = default_vol2_spec()) {
  out <- vector("list", length(solutions))
  
  for (i in seq_along(solutions)) {
    ans <- tryCatch({
      v1 <- estimate_vol1_from_solution(solutions[[i]], tol = tol, spec = vol1_spec)
      v2 <- estimate_vol2_from_solution(solutions[[i]], tol = tol, spec = vol2_spec)
      list(vol1 = v1$volume, vol2 = v2$volume, total_volume = v1$volume + v2$volume)
    }, error = function(e) e)
    
    if (inherits(ans, "error")) {
      if (!quiet) message("Frame ", i, " failed: ", conditionMessage(ans))
      out[[i]] <- list(
        vol1 = NA_real_,
        vol2 = NA_real_,
        total_volume = NA_real_,
        error = conditionMessage(ans)
      )
    } else {
      out[[i]] <- ans
    }
  }
  
  tibble(
    frame = seq_along(out),
    vol1 = purrr::map_dbl(out, ~ .x$vol1 %||% NA_real_),
    vol2 = purrr::map_dbl(out, ~ .x$vol2 %||% NA_real_),
    total_volume = purrr::map_dbl(out, ~ .x$total_volume %||% NA_real_),
    error = purrr::map_chr(out, ~ .x$error %||% NA_character_)
  )
}

# ============================================================
# Plotly static plot
# ============================================================

#' Plot a linkage configuration with selected volume hulls overlaid.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param sol One solved linkage state, usually a named list of 3D points.
#' @param link_map Named list specifying which point pairs form plotted links.
#' @param volumes Input argument used by this helper.
#' @param vol1_spec Point specification for volume compartment 1.
#' @param vol2_spec Point specification for volume compartment 2.
#' @param link_col Input argument used by this helper.
#' @param omit_points Optional point names to omit from range calculations or plots.
#' @param label_points Input argument used by this helper.
#' @param point_col Input argument used by this helper.
#' @param point_size Input argument used by this helper.
#' @param line_width Input argument used by this helper.
#' @param show_labels Input argument used by this helper.
#' @param show_vertices Input argument used by this helper.
#' @param opacity_vol1 Input argument used by this helper.
#' @param opacity_vol2 Input argument used by this helper.
#' @param pitch Camera pitch angle in degrees.
#' @param yaw Camera yaw angle in degrees.
#' @param roll Camera roll angle in degrees.
#' @param zoom Camera distance or scale factor.
#' @param axis_pad Input argument used by this helper.
#' @param remove_axes Logical; if TRUE, hide 3D axes.
#' @param bgcolor Scene background color.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A Plotly figure.
#' @keywords internal
plot_linkage_volumes_static <- function(sol,
                                        link_map,
                                        volumes = c("vol1", "vol2"),
                                        vol1_spec = default_vol1_spec(),
                                        vol2_spec = default_vol2_spec(),
                                        link_col = NULL,
                                        omit_points = NULL,
                                        label_points = NULL,
                                        point_col = "black",
                                        point_size = 4,
                                        line_width = 6,
                                        show_labels = TRUE,
                                        show_vertices = FALSE,
                                        opacity_vol1 = 0.20,
                                        opacity_vol2 = 0.20,
                                        pitch = 15,
                                        yaw = 110,
                                        roll = 0,
                                        zoom = 2,
                                        axis_pad = 0.05,
                                        remove_axes = FALSE,
                                        bgcolor = "white",
                                        tol = 1e-10) {
  volumes <- match.arg(volumes, choices = c("vol1", "vol2"), several.ok = TRUE)
  
  if (is.null(link_col)) {
    link_col <- setNames(rep("black", length(link_map)), names(link_map))
  } else {
    if (is.null(names(link_col))) {
      stop("link_col must be a named vector with names matching link_map.")
    }
    missing_cols <- setdiff(names(link_map), names(link_col))
    if (length(missing_cols) > 0) {
      fill_cols <- setNames(rep("black", length(missing_cols)), missing_cols)
      link_col <- c(link_col, fill_cols)
    }
    link_col <- link_col[names(link_map)]
  }
  
  filter_pts_df <- function(df, omit_points = NULL) {
    if (is.null(omit_points)) return(df)
    df %>% dplyr::filter(!name %in% omit_points)
  }
  
  filter_link_map <- function(link_map, omit_points = NULL) {
    if (is.null(omit_points)) return(link_map)
    keep <- vapply(link_map, function(x) !any(x %in% omit_points), logical(1))
    link_map[keep]
  }
  
  filter_vertices_mat <- function(verts, omit_points = NULL) {
    if (is.null(verts)) return(NULL)
    if (is.null(omit_points)) return(verts)
    rn <- rownames(verts)
    if (is.null(rn)) return(verts)
    verts[!rn %in% omit_points, , drop = FALSE]
  }
  
  get_volume_vertices_local <- function(sol, vv, tol = 1e-10, omit_points = NULL) {
    verts <- get_volume_vertices(
      sol,
      volume = vv,
      tol = tol,
      vol1_spec = vol1_spec,
      vol2_spec = vol2_spec
    )
    filter_vertices_mat(verts, omit_points = omit_points)
  }
  
  link_map_use <- filter_link_map(link_map, omit_points = omit_points)
  
  pts_df <- solution_points_df(sol) %>%
    filter_pts_df(omit_points = omit_points) %>%
    mutate(label = if (isTRUE(show_labels)) {
      if (is.null(label_points)) name else ifelse(name %in% label_points, name, "")
    } else {
      ""
    })
  
  seg_df <- segments_from_solution(sol, link_map_use)
  
  expand <- function(r, pad = 0.05) {
    d <- diff(r)
    if (!is.finite(d) || d == 0) d <- 1
    c(r[1] - pad * d, r[2] + pad * d)
  }
  
  bounds <- list(
    x = expand(range(pts_df$x, na.rm = TRUE), axis_pad),
    y = expand(range(pts_df$y, na.rm = TRUE), axis_pad),
    z = expand(range(pts_df$z, na.rm = TRUE), axis_pad)
  )
  
  scene <- build_scene(
    bounds = bounds,
    remove_axes = remove_axes,
    fix_camera = TRUE,
    pitch = pitch,
    yaw = yaw,
    roll = roll,
    zoom = zoom,
    bgcolor = bgcolor
  )
  
  p <- plotly::plot_ly()
  
  if (length(link_map_use) > 0) {
    for (lnk in names(link_map_use)) {
      dsub <- seg_df %>% dplyr::filter(link == lnk)
      p <- plotly::add_trace(
        p,
        type = "scatter3d",
        mode = "lines",
        x = dsub$x,
        y = dsub$y,
        z = dsub$z,
        name = lnk,
        line = list(width = line_width, color = unname(link_col[[lnk]])),
        hovertemplate = paste0("<b>", lnk, "</b><extra></extra>"),
        showlegend = FALSE
      )
    }
  }
  
  p <- plotly::add_trace(
    p,
    data = pts_df,
    x = ~x, y = ~y, z = ~z,
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
  
  if (show_vertices) {
    for (vv in volumes) {
      verts <- tryCatch(
        get_volume_vertices_local(sol, vv, tol = tol, omit_points = omit_points),
        error = function(e) NULL
      )
      
      if (!is.null(verts) && nrow(verts) > 0) {
        p <- plotly::add_trace(
          p,
          type = "scatter3d",
          mode = "markers",
          x = verts[, 1],
          y = verts[, 2],
          z = verts[, 3],
          name = paste0(vv, "_verts"),
          marker = list(size = 3, color = point_col),
          hoverinfo = "skip",
          showlegend = FALSE
        )
      }
    }
  }
  
  for (vv in volumes) {
    verts <- tryCatch(
      get_volume_vertices_local(sol, vv, tol = tol, omit_points = omit_points),
      error = function(e) NULL
    )
    
    if (!is.null(verts) && nrow(verts) >= 4) {
      tri <- geometry::convhulln(verts, options = "Qt") - 1L
      p <- plotly::add_trace(
        p,
        type = "mesh3d",
        x = verts[, 1],
        y = verts[, 2],
        z = verts[, 3],
        i = tri[, 1],
        j = tri[, 2],
        k = tri[, 3],
        name = vv,
        opacity = if (vv == "vol1") opacity_vol1 else opacity_vol2,
        hoverinfo = "skip",
        showscale = FALSE,
        showlegend = FALSE
      )
    }
  }
  
  p %>%
    plotly::layout(
      scene = scene,
      showlegend = FALSE,
      uirevision = TRUE
    )
}

# ============================================================
# x-slice area and slab volume profiles
# ============================================================

#' Compute the scalar 2D cross product.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param a First vector, point, or object used in the calculation.
#' @param b Second vector, point, or object used in the calculation.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
cross2 <- function(a, b) a[1] * b[2] - a[2] * b[1]

poly_area_2d <- function(mat2d) {
  if (is.null(mat2d) || nrow(mat2d) < 3) return(0)
  x <- mat2d[, 1]
  y <- mat2d[, 2]
  xp <- c(x[-1], x[1])
  yp <- c(y[-1], y[1])
  abs(sum(x * yp - y * xp)) / 2
}

#' Order unordered 2D points around their centroid.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param mat2d Two-column matrix of 2D points.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
order_polygon_points_2d <- function(mat2d) {
  if (is.null(mat2d) || nrow(mat2d) < 3) return(mat2d)
  ctr <- colMeans(mat2d)
  ang <- atan2(mat2d[, 2] - ctr[2], mat2d[, 1] - ctr[1])
  mat2d[order(ang), , drop = FALSE]
}

#' Intersect a 3D hull with a plane perpendicular to the x axis.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param verts Matrix of 3D hull vertices.
#' @param x0 X coordinate of the slicing plane.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
slice_hull_at_x <- function(verts, x0, tol = 1e-10) {
  verts <- as.matrix(verts)
  if (nrow(verts) < 4) {
    return(list(points3d = NULL, points2d = NULL, area = 0))
  }
  
  tri <- tryCatch(
    geometry::convhulln(verts, options = "Qt"),
    error = function(e) NULL
  )
  if (is.null(tri)) {
    return(list(points3d = NULL, points2d = NULL, area = 0))
  }
  
  pts <- list()
  k <- 1L
  
  for (i in seq_len(nrow(tri))) {
    face <- tri[i, ]
    P <- verts[face, , drop = FALSE]
    
    edges <- list(
      c(1, 2),
      c(2, 3),
      c(3, 1)
    )
    
    face_pts <- matrix(NA_real_, nrow = 0, ncol = 3)
    
    for (e in edges) {
      p1 <- P[e[1], ]
      p2 <- P[e[2], ]
      
      dx1 <- p1[1] - x0
      dx2 <- p2[1] - x0
      
      if (abs(dx1) < tol && abs(dx2) < tol) {
        face_pts <- rbind(face_pts, p1, p2)
      } else if (abs(dx1) < tol) {
        face_pts <- rbind(face_pts, p1)
      } else if (abs(dx2) < tol) {
        face_pts <- rbind(face_pts, p2)
      } else if (dx1 * dx2 < 0) {
        t <- (x0 - p1[1]) / (p2[1] - p1[1])
        pint <- p1 + t * (p2 - p1)
        face_pts <- rbind(face_pts, pint)
      }
    }
    
    if (nrow(face_pts) > 0) {
      face_pts <- unique_rows_tol(face_pts, tol = tol)
      pts[[k]] <- face_pts
      k <- k + 1L
    }
  }
  
  if (length(pts) == 0) {
    return(list(points3d = NULL, points2d = NULL, area = 0))
  }
  
  pts3d <- do.call(rbind, pts)
  pts3d <- unique_rows_tol(pts3d, tol = tol)
  
  if (nrow(pts3d) < 3) {
    return(list(points3d = pts3d, points2d = NULL, area = 0))
  }
  
  pts2d <- pts3d[, 2:3, drop = FALSE]
  pts2d <- order_polygon_points_2d(pts2d)
  area <- poly_area_2d(pts2d)
  
  list(points3d = pts3d, points2d = pts2d, area = area)
}

#' Create a common x-axis grid spanning volume hulls across all frames.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param solutions List of solved linkage states, one element per frame.
#' @param volumes Input argument used by this helper.
#' @param n_slices Input argument used by this helper.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @param vol1_spec Point specification for volume compartment 1.
#' @param vol2_spec Point specification for volume compartment 2.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
get_global_x_grid <- function(solutions,
                              volumes = c("vol1", "vol2"),
                              n_slices = 60,
                              tol = 1e-10,
                              vol1_spec = default_vol1_spec(),
                              vol2_spec = default_vol2_spec()) {
  all_x <- c()
  
  for (i in seq_along(solutions)) {
    for (v in volumes) {
      verts <- tryCatch(
        get_volume_vertices(
          solutions[[i]],
          volume = v,
          tol = tol,
          vol1_spec = vol1_spec,
          vol2_spec = vol2_spec
        ),
        error = function(e) NULL
      )
      if (!is.null(verts) && nrow(verts) > 0) {
        all_x <- c(all_x, verts[, 1])
      }
    }
  }
  
  if (length(all_x) < 2) stop("Could not determine x range from solutions.")
  xr <- range(all_x, na.rm = TRUE)
  seq(xr[1], xr[2], length.out = n_slices)
}

#' Estimate cross-sectional area or segment-volume profiles along the x axis.
#'
#' This helper is part of the linkage/volume workflow. It assumes coordinates are
#' numeric and that named points use the same conventions as the rest of the scripts.
#' @param solutions List of solved linkage states, one element per frame.
#' @param volumes Input argument used by this helper.
#' @param n_slices Input argument used by this helper.
#' @param tol Numerical tolerance used for degeneracy checks and floating-point comparisons.
#' @param vol1_spec Point specification for volume compartment 1.
#' @param vol2_spec Point specification for volume compartment 2.
#' @return A value, data frame, list, or plot object used by the linkage workflow.
#' @keywords internal
estimate_xslice_profiles <- function(solutions,
                                     volumes = c("vol1", "vol2"),
                                     n_slices = 60,
                                     tol = 1e-10,
                                     vol1_spec = default_vol1_spec(),
                                     vol2_spec = default_vol2_spec()) {
  x_grid <- get_global_x_grid(
    solutions = solutions,
    volumes = volumes,
    n_slices = n_slices,
    tol = tol,
    vol1_spec = vol1_spec,
    vol2_spec = vol2_spec
  )
  
  out <- vector("list", length(solutions) * length(volumes))
  kk <- 1L
  
  for (i in seq_along(solutions)) {
    sol <- solutions[[i]]
    
    for (v in volumes) {
      verts <- tryCatch(
        get_volume_vertices(
          sol,
          volume = v,
          tol = tol,
          vol1_spec = vol1_spec,
          vol2_spec = vol2_spec
        ),
        error = function(e) NULL
      )
      
      if (is.null(verts) || nrow(verts) < 4) {
        dat <- tibble(
          frame = i,
          volume = v,
          x = x_grid,
          area_yz = 0
        )
      } else {
        dat <- tibble(
          frame = i,
          volume = v,
          x = x_grid,
          area_yz = map_dbl(x_grid, ~ slice_hull_at_x(verts, .x, tol = tol)$area)
        )
      }
      
      out[[kk]] <- dat
      kk <- kk + 1L
    }
  }
  
  area_dat <- bind_rows(out) %>%
    group_by(volume, frame) %>%
    arrange(x, .by_group = TRUE) %>%
    mutate(
      x_next = lead(x),
      area_next = lead(area_yz),
      dx = x_next - x,
      slab_volume = if_else(
        !is.na(dx),
        (area_yz + area_next) / 2 * dx,
        NA_real_
      ),
      x_mid = if_else(!is.na(x_next), (x + x_next) / 2, NA_real_)
    ) %>%
    ungroup()
  
  total_area <- area_dat %>%
    group_by(frame, x) %>%
    summarise(area_yz = sum(area_yz, na.rm = TRUE), .groups = "drop") %>%
    group_by(frame) %>%
    arrange(x, .by_group = TRUE) %>%
    mutate(
      x_next = lead(x),
      area_next = lead(area_yz),
      dx = x_next - x,
      slab_volume = if_else(
        !is.na(dx),
        (area_yz + area_next) / 2 * dx,
        NA_real_
      ),
      x_mid = if_else(!is.na(x_next), (x + x_next) / 2, NA_real_)
    ) %>%
    ungroup() %>%
    mutate(volume = "total") %>%
    relocate(frame, volume, x)
  
  area_all <- bind_rows(area_dat, total_area) %>%
    arrange(volume, frame, x)
  
  slab_all <- area_all %>%
    filter(!is.na(x_mid), !is.na(slab_volume)) %>%
    group_by(volume, x_mid) %>%
    arrange(frame, .by_group = TRUE) %>%
    mutate(
      dV_frame = slab_volume - lag(slab_volume),
      dV_frame = replace_na(dV_frame, 0)
    ) %>%
    ungroup()
  
  list(
    x_grid = x_grid,
    area = area_all,
    slab = slab_all
  )
}