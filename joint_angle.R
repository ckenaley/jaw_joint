# ============================================================
# Live-feeding 3D kinematics and joint-angle summaries
#
# This script computes gape, symphyseal angle, lateral expansion
# from digitized 3D landmarks. The processed curves
# are smoothed and exported for downstream linkage-model comparisons.
# ============================================================

library(ggbeeswarm)
library(ggthemes)
library(signal)
library(cowplot)
library(tidyverse)
library(lme4)
library(nationalparkcolors)

# ============================================================
# helpers
# ============================================================
 
source("linkage_helpers.R")
angle_3D <- function(A, B, C, degrees = TRUE, tol = 1e-12) {
  v1 <- A - B
  v2 <- C - B
  
  n1 <- sqrt(sum(v1^2))
  n2 <- sqrt(sum(v2^2))
  
  if (n1 < tol || n2 < tol) {
    return(NA_real_)
  }
  
  cosang <- sum(v1 * v2) / (n1 * n2)
  cosang <- max(-1, min(1, cosang))
  
  ang <- acos(cosang)
  
  if (degrees) {
    ang <- ang * 180 / pi
  }
  
  ang
}

kinematics_vars <- c(
  "sym_angle",
  "lat_exp",
  "ph_strain",
  "jaw_angle"
)

panel_levels <- c(
  "Gape\nangle (°)",
  "Symphyseal\nangle (°)",
  "Lateral\nexpansion (%)"
)

panel_labs <- c(
  jaw_angle = "Gape\nangle (°)",
  sym_angle = "Symphyseal\nangle (°)",
  lat_exp = "Lateral\nexpansion (%)"
)

# ============================================================
# read data
# ============================================================

kin_dat <- read_csv("data/bass_kinematic_data.csv", show_col_types = FALSE)

# ============================================================
# compute kinematic variables
# ============================================================


ang_dat <- kin_dat %>%
  drop_na() %>%
  group_by(fish,trial,frame) %>%
  dplyr::filter(length(landmark)==8) %>%
  group_by(fish,trial) %>%
  mutate(frame=frame-first(frame)) %>% 
  dplyr::select(-pt_num) %>%
  pivot_wider(
    values_from = x:z,
    names_from = landmark
  ) %>%
  rowwise() %>%
  mutate(
    sym_angle = angle_3D(
      c(x_quad, y_quad, z_quad),
      c(x_lj,   y_lj,   z_lj),
      c(x_basi, y_basi, z_basi)
    ) * 2,
    
    jaw_angle = angle_3D(
      c(x_pm,   y_pm,   z_pm),
      c(x_quad, y_quad, z_quad),
      c(x_lj,   y_lj,   z_lj)
    ),
    
    lat_exp = point_plane_dist_3d(
      c(x_lj,   y_lj,   z_lj),
      c(x_basi, y_basi, z_basi),
      c(x_pm,   y_pm,   z_pm),
      c(x_so,   y_so,   z_so)
    ),
    
    ph_strain = dist3(
      c(x_lj,   y_lj,   z_lj),
      c(x_basi, y_basi, z_basi)
    )
  ) %>%
  ungroup() %>%
  group_by(fish, trial) %>%
  mutate(
    sym_angle = sym_angle - first(sym_angle),
    jaw_angle = jaw_angle - first(jaw_angle),
    ph_strain = (ph_strain - first(ph_strain)) / first(ph_strain) * 100,
    lat_exp   = (lat_exp   - first(lat_exp))   / first(lat_exp)   * 100
  ) %>%
  ungroup() %>%
  mutate(
    fish = recode(
      fish,
      bass_1 = "Bass 1",
      bass_2 = "Bass 2",
      bass_4 = "Bass 3"
    )
  ) %>%
  group_by(fish, trial) %>%
  mutate(
    per_open = plyr::round_any(
      time_ms / time_ms[which.max(jaw_angle)],
      0.025
    ),
    per_abduct = plyr::round_any(
      time_ms / time_ms[which.max(sym_angle)],
      0.025
    )
  ) %>%
  ungroup()

# ============================================================
# jaw lengths
# ============================================================

jls <- tibble(
  fish = unique(ang_dat$fish),
  jl = c(3.01, 2.44, 2.87)
)

jl <- mean(jls$jl)

kin_dat %>%
  drop_na() %>%
  mutate(
    fish = recode(
      fish,
      bass_1 = "Bass 1",
      bass_2 = "Bass 2",
      bass_4 = "Bass 3"
    )
  ) %>% 
  group_by(fish,trial,frame) %>%
  dplyr::filter(length(landmark)==8) %>%
  group_by(fish,trial) %>%
  mutate(frame=frame-first(frame)) %>% 
  dplyr::select(-pt_num) %>%
  pivot_wider(
    values_from = x:z,
    names_from = landmark
  ) %>%
  dplyr::filter(frame==min(frame)) %>% 
  rowwise() %>% 
  group_by(fish, trial) %>% 
  mutate(
    lat_exp = point_plane_dist_3d(
      c(x_lj,   y_lj,   z_lj),
      c(x_basi, y_basi, z_basi),
      c(x_pm,   y_pm,   z_pm),
      c(x_so,   y_so,   z_so)
    )*2
  ) %>% 
  dplyr::select(fish,trial,lat_exp) %>% 
  left_join(jls) %>% 
  mutate(lat_exp=lat_exp/jl) %>% 
  ungroup() %>% 
  summarize(m_lat_exp=mean(lat_exp),sd=sd(lat_exp))


# ============================================================
# summaries
# ============================================================

ang_sum <- ang_dat %>%
  dplyr::select(fish, trial, per_open, per_abduct, all_of(kinematics_vars)) %>%
  pivot_longer(
    cols = all_of(kinematics_vars),
    values_to = "value"
  ) %>%
  group_by(fish, per_open, per_abduct, name) %>%
  summarise(
    m = mean(value, na.rm = TRUE),
    se = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
    .groups = "drop"
  ) %>%
  dplyr::filter(per_open <= 1.4)


# ============================================================
# smoothed species-average trajectories
# ============================================================

ang_sum2 <- ang_sum %>%
  dplyr::filter(per_open <= 1.4) %>%
  mutate(per_open = plyr::round_any(per_open, 0.05)) %>%
  group_by(per_open, name) %>%
  summarise(
    angle = mean(m, na.rm = TRUE),
    .groups = "drop"
  )

spec_data_for_model <- ang_sum2 %>%
  group_by(name) %>%
  group_modify(~{
    fit <- smooth.spline(
      x = .x$per_open,
      y = .x$angle,
      spar = 0.5
    )
    pred <- predict(fit, x = .x$per_open)
    
    tibble(
      per_open = pred$x,
      angle = pred$y
    )
  }) %>%
  ungroup()

# optional exploratory checks
ang_sum2 %>%
  ggplot(aes(per_open, angle, col = name)) +
  geom_point() +
  geom_point(shape = 3)


spec_data_for_model %>%
  ggplot(aes(per_open, angle, col = name)) +
  geom_point() +
  geom_point(data = ang_sum2, shape = 3)

saveRDS(spec_data_for_model, "data/spec_data_for_model.RDS")

# ============================================================
# trim to opening window for downstream summaries
# ============================================================

ang_dat <- ang_dat %>%
  dplyr::filter(per_open <= 1.3)

# summary at maximum gape
spec_sum <- ang_dat %>%
  dplyr::select(fish, trial, per_open, per_abduct, sym_angle,lat_exp, jaw_angle) %>%
  pivot_longer(sym_angle:jaw_angle) %>%
  dplyr::filter(per_open == 1) %>%
  group_by(fish, trial, name) %>%
  summarise(
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(fish, name) %>%
  summarise(
    mean = mean(max, na.rm = TRUE),
    sd = sd(max, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(name)

spec_sum %>%
  group_by(name) %>%
  summarise(
    m = mean(mean, na.rm = TRUE),
    sd = sd(mean, na.rm = TRUE),
    .groups = "drop"
  )

# long format version
ang_dat_long <- ang_dat %>%
  dplyr::select(fish, trial, per_open, sym_angle, lat_exp, jaw_angle) %>%
  pivot_longer(
    cols = sym_angle:jaw_angle,
    values_to = "angle"
  )


# ============================================================
# plot data for main figure
# ============================================================

ang_sum_open <- ang_dat %>%
  dplyr::select(fish, trial, per_open, sym_angle, lat_exp, jaw_angle) %>%
  pivot_longer(sym_angle:jaw_angle) %>%
  group_by(fish, per_open, name) %>%
  summarise(
    m = mean(value, na.rm = TRUE),
    se = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
    .groups = "drop"
  )

per_open_full <- ang_sum_open %>%
  mutate(per_open=per_open-1) %>% 
  dplyr::filter(per_open == 0)

plot_dat <- ang_sum_open %>%
  mutate(per_open=per_open-1) %>% 
  dplyr::filter(name %in% c("jaw_angle", "sym_angle", "lat_exp")) %>%
  mutate(
    panel = recode(name, !!!panel_labs),
    panel = factor(panel, levels = panel_levels),
    se_mult = 1
  )

vline_dat <- per_open_full %>%
  dplyr::filter(name %in% c("jaw_angle", "sym_angle", "lat_exp")) %>%
  mutate(
    panel = recode(name, !!!panel_labs),
    panel = factor(panel, levels = panel_levels),
    se_mult = 1
  )

# ============================================================
# plot
# ============================================================

p_onecol <-plot_dat %>%  ggplot() +
  geom_errorbar(
    aes(
      x = per_open,
      y = m,
      ymin = m - se * se_mult,
      ymax = m + se * se_mult
    ),
    alpha = 0.5
  ) +
  geom_point(
    aes(x = per_open, y = m),
    size = 2,
    alpha = 0.6
  ) +
  geom_vline(
    data = vline_dat,
    aes(xintercept = per_open),
    colour = "gray40",
    linetype = 4
  ) +
  facet_grid(
    panel ~ fish,
    scales = "free_y",
    switch = "y"
  ) +
  theme_classic(15) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y = element_text(size = 12),
    axis.text.y = element_text(size = 10)
  ) +
  labs(
    x = "Relative Time to Maximum Gape",
    y = NULL
  )

p_onecol

ggsave(
  filename = "manuscript/livekin.pdf",
  plot = p_onecol,
  height = 9,
  width = 7
)

# ============================================================
# lme model fitting sym angle ~ lat expansion
# ============================================================


fit <- lmer(
  sym_angle ~ lat_exp +
    (1 + lat_exp || fish), 
  data = ang_dat,
  REML = T,
  control = lmerControl(optimizer = "bobyqa")
)


anova(fit)
summary(fit)
car::Anova(fit)
coef(fit)

pal <- nationalparkcolors::park_palette("Badlands",5)

p_sym_lat <- ang_dat %>% 
  mutate(sym_ang_fit=predict(fit)) %>% 
  ggplot(aes(lat_exp,sym_angle,col=trial))+
  geom_point(alpha=0.6)+
  geom_smooth(alpha=0.6,method="lm",se=F)+
  geom_line(aes(lat_exp,sym_ang_fit),col="black",linewidth = 1.5,alpha=0.8)+
  facet_grid(
    .~ fish,
    scales = "free_y",
    switch = "y"
  )+
  scale_color_manual(values=pal)+
  theme_classic(15)+
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    legend.position = c(0.9,0.25)
  )+
  xlab("Lateral Expansion (%)")+
  ylab("Symphyseal nangle (°)")

ggsave(
  filename = "manuscript/symlat.pdf",
  plot = p_sym_lat,
  height = 5,
  width = 7
)


