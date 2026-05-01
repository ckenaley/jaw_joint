# ============================================================
# Largemouth bass linkage model workflow
#
# This script ties together the helper, solver, and volume functions. It
# reads CT landmarks, recasts them into a reference plane, scales the
# geometry to live specimens, drives the linkage with empirical feeding
# kinematics, and compares modeled outputs with live data.
# ============================================================

source("linkage_helpers.R")
source("linkage_solver.R")
source("volume.R")

suppressPackageStartupMessages({
  library(plotly)
  library(geomorph)
  library(rgl)
  #use devtools::install_github("ckenaley/geomorphcompanion")
  library(geomorphcompanion) 
  library(cowplot)
  library(ggthemes)
  library(ggrepel)
  library(tidyverse)
})

# ------------------------------------------------------------------
# Read specimen and digitize landmarks
# ------------------------------------------------------------------
if(F){
spec <- read.ply2(
  file = "data/MCZ48917_Micropterus_nigricans.ply",
  ShowSpecimen = FALSE
)

# Landmark definitions
# 1/A = Suspensorium at neurocranium
# 2/B = Quadrate at lower jaw
# 3/C = Suspensorium at interhyal
# 4/D = Ceratohyal symphysis
# 5/E = Lower jaw at lower jaw joint
# 6/F = Dorsoposterior tip of supraoccipital crest
# 7/G = Parasphenoid
# 8/H = Ethmoid
# 9/I = Anterior tip of ventral cleithrum
# 10/J = Neurocranium at cleithral joint
# 11/K = Premax (for volume)
# 12L = Anterior Pterygoid (for volume)
# 13/M = basioaccipital (for alignment)

digit.fixed2(spec = spec, out.dir = "data")
}


bass_points <- read_csv("data/MCZ48917_Micropterus_nigricans.ply.csv", show_col_types = FALSE)

pts <- list(
  A = 1,
  B = 2,
  C = 3,
  D = 4,
  E = 5,
  `F` = 6,
  G = 7,
  H = 8,
  I = 9,
  J = 10,
  K = 11,
  L = 12,
  M = 13
)

bass_points <- bass_points %>% mutate(name=names(pts))

plane_pts <- pts[c("M", "G", "H")]
axis_pair <- pts[c("M", "H")]

# ------------------------------------------------------------------
# Recast specimen to reference plane and merge selected coordinates
# ------------------------------------------------------------------

spec_cast <- recast_to_plane_named(
  bass_points %>% select(-name),
  pts = pts,
  plane = plane_pts,
  origin = "first",
  plane_axes = c("x", "z"),
  axis_pair = axis_pair,
  axis = "x"
)

coords <- coords_orig <- merge_coords(spec_cast$all_coords, c("E", "D", "I","K","M","G"))

# ------------------------------------------------------------------
# Quick visualization of reflected linkage geometry
# ------------------------------------------------------------------

link_map_plot <- list(
  N  = c("A", "A_ref"),
  Sl = c("A_ref", "C_ref"),
  Jl = c("B_ref", "E"),
  Hl = c("D", "C_ref"),
  Sr = c("A", "C"),
  Jr = c("B", "E"),
  H  = c("D", "C")
)

plot_symmetry_3d(
  coords,
  link_map = link_map_plot,
  pitch = 15,
  yaw = 110,
  zoom = 2
)

# ------------------------------------------------------------------
# Scale coordinates of CT specimen based on man jaw length of live specimens
# ------------------------------------------------------------------
bl_live <- mean(c(16.3, 20.4, 22.08))
jl_live <- 2.7

get_pt <- function(dat, nm) {
  ii <- match(nm, dat$name)
  if (is.na(ii)) stop("Point not found: ", nm)
  as.numeric(dat[ii, c("x", "y", "z")])
}

# define CT jaw length on the ORIGINAL geometry
# change these names if your intended jaw-length pair is different
jl_ct <- dist3(get_pt(coords, "B"), get_pt(coords, "E"))

# target starting C-C_ref spacing from CT jaw length
start_C <- 0.45 * jl_ct

# reposition C/B first on the unscaled CT geometry
coords <- set_start_C(coords, target_CCref = start_C)

jl_ct2 <- dist3(get_pt(coords, "B"), get_pt(coords, "E"))

# confirm new C-C_ref spacing before scaling
ct_start <- dist3(get_pt(coords, "C"), get_pt(coords, "C_ref"))
print(round(ct_start,4)==round(start_C,4))


plot_symmetry_3d(
  coords,
  link_map = link_map_plot,
  pitch = 15,
  yaw = 110,
  zoom = 2
)

# NOW recompute CT jaw length on the adjusted starting geometry
jl_ct_adj <- dist3(get_pt(coords, "B"), get_pt(coords, "E"))

# scale factor based on adjusted geometry
scale_factor <- jl_ct_adj / jl_live

# scale all coordinates
coords <- coords %>%
  mutate(
    x = x / scale_factor,
    y = y / scale_factor,
    z = z / scale_factor
  )

# check scaling
jl_ct2 <- dist3(get_pt(coords, "B"), get_pt(coords, "E"))
print(round(jl_ct2, 2) == round(jl_live, 2))
print(jl_ct2)



#Body lengths from Camp and Brainer (2014)
camp_bl <- c(24.2, 28.1, 27.7) %>% mean


camp_rat <- bl_live/camp_bl
# ------------------------------------------------------------------
# Pull feeding data
# ------------------------------------------------------------------

spec_dat <- readRDS("data/spec_data_for_model.RDS")

E_angles <- spec_dat %>%
  filter(name == "jaw_angle") %>%
  pull(angle) * -1

ph_strain <- spec_dat %>%
  filter(name == "ph_strain") %>%
  pull(angle)

sym_ang <- spec_dat %>%
  filter(name == "sym_angle") %>%
  pull(angle)


n_frames <- length(E_angles)

#Camp and Brainerd (2014) urohyal translations, x=0.85 and z=0.9 mm

retract <- camp_rat*c(0.85,0.9)


#smooth driver translations
trans_x <- track_profile(
  n = n_frames,
  max_val = retract[1],
  k = 8,
  mid = 0.75
)*0.33

trans_z <- track_profile(
  n = n_frames,
  max_val = retract[2],
  k = 8,
  mid = 0.75
)*0.33

plot(trans_z)

# Neurocranial rotation profile
neuro_rot <- track_profile(
  n = n_frames,
  max_val = 10,
  k = 8,
  mid = 0.75
)


# Optional diagnostics
plot(trans_x,type = "l")
points(trans_z,col="blue",type = "l")
plot(neuro_rot, type = "l", xlab = "Frame",ylab= expression(paste("Rotation (", degree, ")")))


# ------------------------------------------------------------------
# Build model
# ------------------------------------------------------------------


model <- make_sym_model(coords)

# Optional trajectory preview for point D rotated about C in xz plane
traj <- point_trajectory(
  model = model,
  point = "D",
  center = "C",
  angles_deg = seq(0, max(E_angles), length.out = n_frames),
  plane = "xz"
)



# ------------------------------------------------------------------
# Assemble solver inputs
# ------------------------------------------------------------------

inputs <- tibble(
  x_shift        = -c(0, diff(trans_x)),
  z_shift        = -c(0, diff(trans_z)),
  E_rot_deg      = E_angles,
  global_rot_deg = c(0,diff(neuro_rot))
)  

# ------------------------------------------------------------------
# Solve linkage path
# ------------------------------------------------------------------

res <- inputs %>%
  solve_path_tbl(
    model = model,
    exclude_from_global = c("I", "D"),
    global_axis = c("J", "J_ref"),
    driver_mode = "step"
      
  )



sols <- res$solutions
print(all(res$ok))

# ------------------------------------------------------------------
# Animate result
# ------------------------------------------------------------------

link_map_anim <- list(
  N   = c("A", "A_ref"),
  Sr  = c("A", "C"),
  Sl  = c("A_ref", "C_ref"),
  Hl  = c("D", "C_ref"),
  Hr  = c("D", "C"),
  Jl  = c("B_ref", "E"),
  Jr  = c("B", "E"),
  Qr  = c("C", "B"),
  Ql  = c("C_ref", "B_ref"),
  ShM = c("D", "I"),
  PhM = c("D", "E")
)

link_col <- setNames(
  c(
    "gray30",
    "#b5594bff",
    "#b5594bff",
    "#5766c0ff",
    "#5766c0ff",
    "#F0CF29",
    "#F0CF29",
    "#b5594bff",
    "#b5594bff",
    "black",
    "black"
  ),
  names(link_map_anim)
)

p_an <- animate_sym_linkage(
  sols,point_col = "black",axis_pad = F,
  link_map = link_map_anim,
  link_col = link_col,point_size = 5,
  show_labels = T,omit_points = c("L","L_ref","K","H","F","G"),
  label_points = LETTERS[1:12],
  zoom=-1.5,
  yaw=-90,
  pitch = 90
)

p_an

panel_levels <- c(
  "Symphyseal\nangle (°)",
  "Lateral\nexpansion (%)"
)

# ------------------------------------------------------------------
# Compare to live bass
# ------------------------------------------------------------------

comp_dat <- spec_dat %>% 
  filter(name%in%c("sym_angle","lat_exp")) %>% 
  mutate(data="live") %>% 
  rename(value=angle) %>% 
  bind_rows(tibble(
    per_open=spec_dat$per_open %>% unique,
    sym_angle=angles_points(sols,pts=c("B","E","B_ref"))$angle,
    lat_exp =dist_series(sols,"C","C_ref"),
    data="linkage model"
  ) %>% mutate(lat_exp=(lat_exp-lat_exp[1])/lat_exp[1]*100,
               sym_angle=(sym_angle-sym_angle[1])) %>% 
    pivot_longer(sym_angle:lat_exp)) %>%
  mutate(name = dplyr::recode(name, sym_ang=panel_levels[1],lat_exp=panel_levels[2]),
         panel = factor(name, levels = panel_levels),
  ) %>% 
  mutate(per_open=per_open-1)


lines_dat <- comp_dat %>% 
  filter(per_open==0|per_open==min(per_open)) %>% 
  group_by(data,panel) %>% 
  mutate(value=max(value)) %>% 
  bind_rows(
    comp_dat %>% 
      filter(per_open==0) %>% 
      group_by(data,panel) %>% 
      mutate(value=0)
  ) 

  
comp_dat %>% 
  ggplot(aes(per_open,value,linetype =data))+
  geom_line(data=lines_dat,aes(per_open,value),col="gray90")+
  geom_line()+
  facet_grid(
    name ~ .,
    scales = "free_y",
    switch = "y"
  )+theme_classic()+
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    legend.position = c(0.25,0.9)
  )+labs(
    linetype=NULL,
    x = "Relative Time to Maximum Gape",
    y = NULL
  )

ggsave("manuscript/modelcomp.pdf",h=6,w=4)
  

# ------------------------------------------------------------------
# Estimate volume along sequence of frames
# ------------------------------------------------------------------

#static plot
p_static_vol <- plot_linkage_volumes_static(
  sol = sols[[1]],
  link_map = link_map_anim,
  volumes = c("vol1", "vol2"),
  link_col = link_col,
  label_points = LETTERS[1:12],
  point_col = "black",
  show_labels = TRUE,
  show_vertices = FALSE,
  pitch = 14,
  yaw = 100,
  roll = -20,
  zoom = 2,
  remove_axes = FALSE
)
p_static_vol

#animation
p_vol <- animate_linkage_volumes(
  solutions = res$solutions,
  link_map = link_map_anim,
  volumes = c("vol1", "vol2"),
  link_col = link_col,
  omit_points = c("I"),
  label_points = LETTERS[1:13],
  point_size = 4,
  line_width = 6,
  show_labels = TRUE,
  show_vertices = FALSE,
  opacity_vol1 = 0.20,
  opacity_vol2 = 0.20,
  pitch = -20,
  yaw = 90,
  roll = 0,
  zoom = 2,
  frame_duration = 50,
  axis_pad = 0.05,
  remove_axes = FALSE,
  bgcolor = "white"
)

p_vol

#save animation
htmlwidgets::saveWidget(
  p_vol,
  "manuscript/sym_linkage_volume_animation.html",
  selfcontained = TRUE
)

#figures
vols_table <- estimate_meshvol_over_path(
  solutions = res$solutions,
  vol1_spec = my_vol1_spec,
  vol2_spec = my_vol2_spec
)
vol_tbl <- estimate_meshvol_over_path(sols)%>% 
  left_join(spec_dat %>% dplyr::select(per_open)%>% 
              mutate(frame=1:n(),
                     per_open=per_open-1)
            )

vol_tbl  %>% 
  pivot_longer(vol1:vol2) %>%
  mutate(name=ifelse(name=="vol1","Suspensorial","Mandibular")) %>% 
  group_by(name,per_open) %>% 
  mutate(per_volume=(value/total_volume)) %>%
  group_by(name) %>% 
  mutate(value=value-value[1],
         total_volume=total_volume-total_volume[1]
         ) %>%
  mutate(d_volume=(value-lag(value))) %>% 
  ggplot(aes(per_open,per_volume,col=name))+
  geom_line()+theme_classic(16)+labs(color="Volume")+theme(legend.position = c(.6,.2))+
  ylab("Prop. Increase")+xlab("Relative Time to Maximum Gape")


vol2 <- vol_tbl %>% 
  pivot_longer(vol1:vol2) %>%
  mutate(name = ifelse(name == "vol1", "Suspensorial", "Mandibular")) %>% 
  group_by(name, per_open) %>% 
  mutate(per_volume = value / total_volume) %>%
  ungroup() %>%
  group_by(name) %>% 
  mutate(
    value = value - value[1],
    total_volume = total_volume - total_volume[1]
  ) %>%
  ungroup()

sf <- max(vol2$value, na.rm = TRUE) / max(vol2$per_volume, na.rm = TRUE)

p_vol1 <- ggplot(vol2, aes(per_open, color = name)) +
  geom_line(aes(y = value), linewidth = 1) +
  geom_line(aes(y = per_volume * sf), linetype = 2, linewidth = 1) +
  scale_y_continuous(
    name = "Prop. Increase",
    sec.axis = sec_axis(~ . / sf, name = "Prop. Total. Vol.")
  ) +
  theme_classic(16) +
  labs(x = "Relative Time to Maximum Gape", color = "") +
  theme(legend.position = c(.7, .2),legend.text = element_text(size=10),legend.background = element_blank())+
  scale_color_manual(values = c("black","gray60"))


# ------------------------------------------------------------------
# Estimate dV/dt with respect to x
# ------------------------------------------------------------------

#define volume 1
my_vol1_spec <- list(
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

my_vol2_spec <- list(
  "B", "B_ref", "D", "E", "K","L","L_ref"
)




vol_x <- estimate_xslice_profiles(
  solutions = res$solutions,
  volumes = c("vol1", "vol2"),
  n_slices = 50,
  vol1_spec = my_vol1_spec,
  vol2_spec = my_vol2_spec
)

p_vol2 <- vol_x$slab %>%
  left_join(spec_dat %>% dplyr::select(per_open)%>% 
              mutate(frame=1:n(),
                     per_open=per_open-1)
  ) %>% 
  mutate(volume=recode(volume,
                       total="Total",
                       vol1="Suspensorial",
                       vol2="Mandibular")) %>% 
  filter(area_yz!=0) %>% 
  ggplot(aes(x = x, y = per_open, fill = dV_frame)) +
  geom_tile() +
  scale_y_reverse(position = "right") +
  scale_x_reverse() +
  facet_wrap(~ volume, ncol = 1) +
  labs(
    x = xlab("x position \n(anterior \u2190\u2192 posterior)"),
    y = "Relative Time to Maximum Gape",
    fill = expression("Segment volume ("*Delta~cm^2/Delta~t*")")
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "bottom"
    )
  )+
  theme_minimal()+
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick"
  )+
  theme(legend.position = "bottom",axis.title.y = element_text(size=16))


#combind figures and save

plot_grid(p_vol1,p_vol2,labels = c("A","B"),ncol = 2,rel_heights = c(.7,1))
ggsave("manuscript/volume.pdf",h=6,w=10)


# ------------------------------------------------------------------
# Model inputs
# ------------------------------------------------------------------

input_dat <- tibble(frame=1:length(trans_x),
                    x=trans_x,
                    z=trans_z,
                    Neurocranium=neuro_rot,
                    `Lower jaw`=E_angles) %>% 
  pivot_longer(x:`Lower jaw`)

p_in1 <- input_dat %>% 
  filter(name%in%c("x","z")) %>%
  ggplot(aes(frame,value,shape=name)) +
  geom_point()+
  labs(shape="")+
  geom_line()+theme_few(15)+theme(legend.position = c(0.8,.2))+
  ylab("Sagittal axis translation (cm)")+
  xlab("Frame")

p_in2 <- input_dat %>% 
  filter(grepl("aw|Neur",name)) %>%
  ggplot(aes(frame,abs(value),shape=name)) +
  geom_point()+
  labs(shape="")+
  geom_line()+theme_few(15)+theme(legend.position = c(0.2,.8))+
  ylab(expression(paste("Rotation (", degree, ")")))+
  xlab("Frame")

plot_grid(p_in2,p_in1,labels = c("A","B"),ncol = 1)
ggsave("manuscript/modelinput.pdf",h=8,w=5)


 # ============================================================
 # Sensitivity analysis: attenuation of Camp & Brainerd driver
 # ============================================================
 
 atten_factors <- c(1,0.75,0.5,0.33,0.2,0.1)
 
 # helper: extract C coordinates from a solution list
 get_C_path <- function(solutions) {
   bind_rows(lapply(seq_along(solutions), function(i) {
     C <- solutions[[i]]$C
     tibble(
       frame = i,
       Cx = C[1],
       Cy = C[2],
       Cz = C[3]
     )
   }))
 }
 
 # helper: extract D coordinates from a solution list
 get_D_path <- function(solutions) {
   bind_rows(lapply(seq_along(solutions), function(i) {
     D <- solutions[[i]]$D
     tibble(
       frame = i,
       Dx = D[1],
       Dy = D[2],
       Dz = D[3]
     )
   }))
 }
 
 # run one attenuation level
 run_attenuation <- function(att,
                             camp_rat,
                             n_frames,
                             E_angles,
                             neuro_rot,
                             model) {
   
   # attenuated driver endpoints
   retract_x_att <- (-0.85 * camp_rat) * att
   retract_z_att <- (-0.90 * camp_rat) * att
   
   # same temporal profile, different magnitude
   trans_x_att <- track_profile(
     n = n_frames,
     max_val = retract_x_att,
     k = 7,
     mid = 0.7
   )
   
   trans_z_att <- track_profile(
     n = n_frames,
     max_val = retract_z_att,
     k = 7,
     mid = 0.7
   )
   
   inputs_att <- tibble(
     x_shift        = c(0, diff(trans_x_att)),
     z_shift        = c(0, diff(trans_z_att)),
     E_rot_deg      = E_angles,
     global_rot_deg = c(0, diff(neuro_rot))
   )
   
   res_att <- inputs_att %>%
     solve_path_tbl(
       model = model,
       exclude_from_global = c("I", "D"),
       global_axis = c("J", "J_ref"),
       driver_mode = "step"
     )
   
   diag_att <- constraint_series(res_att$solutions) %>%
     mutate(
       attenuation = 1-att,
       atten_label = sprintf("%.2f", att)
     )
   
   C_att <- get_C_path(res_att$solutions) %>%
     mutate(
       attenuation = att,
       atten_label = sprintf("%.2f", att)
     ) %>%
     group_by(attenuation, atten_label) %>%
     mutate(
       Cy_change = Cy - first(Cy),
       C_disp = sqrt((Cx - first(Cx))^2 + (Cy - first(Cy))^2 + (Cz - first(Cz))^2)
     ) %>%
     ungroup()
   
   D_att <- get_D_path(res_att$solutions) %>%
     mutate(
       attenuation = att,
       atten_label = sprintf("%.2f", att)
     ) %>%
     group_by(attenuation, atten_label) %>%
     mutate(
       D_disp = sqrt((Dx - first(Dx))^2 + (Dy - first(Dy))^2 + (Dz - first(Dz))^2)
     ) %>%
     ungroup()
   
   list(
     attenuation = 1-att,
     retract_x = retract_x_att,
     retract_z = retract_z_att,
     trans_x = trans_x_att,
     trans_z = trans_z_att,
     inputs = inputs_att,
     res = res_att,
     diag = diag_att,
     C_path = C_att,
     D_path = D_att
   )
 }
 
 # run all attenuation levels
 sens_runs <- lapply(
   atten_factors,
   run_attenuation,
   camp_rat = camp_rat,
   n_frames = n_frames,
   E_angles = E_angles,
   neuro_rot = 0,
   model = model
 )
 
 names(sens_runs) <- paste0("att_", atten_factors)
 
 # bind diagnostics
 C_sens    <- bind_rows(lapply(sens_runs, `[[`, "C_path"))
 
 
 C_sens <- C_sens %>%ungroup %>% 
   group_by(attenuation, atten_label) %>%
   mutate(
     lat_exp = (Cy-Cy[1])/Cy[1]
   ) %>%
   left_join(spec_dat %>% mutate(frame=1:n()) %>% dplyr::select(frame,per_open)) %>% 
   ungroup()
 

#plot and figures
 
 lab_dat <- C_sens %>%
   group_by(atten_label) %>%
   dplyr::filter(per_open== max(per_open)) %>% 
   ungroup()
 
 spec_latex <- spec_dat %>% 
   dplyr::filter(name=="lat_exp") %>% 
   pivot_wider(names_from = "name",values_from = "angle")
 
 C_sens %>% ggplot( aes(per_open, lat_exp*100, group = atten_label, colour = atten_label)) +
   geom_line(data=spec_latex,col="red",aes(per_open, lat_exp),inherit.aes = F)+
   geom_line(linewidth = 1.1, show.legend = FALSE) +
   geom_text_repel(
     data = lab_dat,
     aes(label = atten_label),
     direction = "y",
     hjust = -1,
     nudge_x = 0.01,
     nudge_y = 0.085,
     segment.color = NA,
     size = 4,
     show.legend = FALSE,
     color="black"
   ) +
   coord_cartesian(clip = "off") +
   theme_classic( 16) +
   theme(
     plot.margin = margin(5.5, 80, 5.5, 5.5)
   )+
   labs(
     x = "Relative Time to Maximum Gape",
     y = "Lateral Epansion (%)",
     color = "Attenuation"
   )+
   scale_color_manual(values=gray.colors(n=length(atten_factors)))
 

ggsave("manuscript/latexpatten.pdf",h=4,w=5) 

 
 


 
