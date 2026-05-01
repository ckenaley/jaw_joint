# ============================================================
# Angular stiffness analysis
#
# This script reads force/displacement trials, converts angular
# displacement and moment-arm data into torque, estimates angular
# stiffness from trial-level regressions, generates manuscript plots,
# and runs species-level contrasts against the rat comparison group.
# ============================================================

library(ggthemes)
library(multcomp)
library(lme4)
library(tidyverse)
# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------

moments <- read_csv("data/moment_lengths.csv", show_col_types = FALSE)

name_key <- tibble(
  species = c("CP", "LB", "RODENT", "WP", "YP"),
  name = c("Chain pickerel", "LM bass", "Rat", "White perch", "Yellow perch")
)

jaw_files <- list.files("data/stiffness", pattern = "JAW", full.names = TRUE)

name_levels <- c("Rat", "LM bass", "Yellow perch", "White perch", "Chain pickerel")

# ------------------------------------------------------------
# Read and process raw data
# ------------------------------------------------------------

dat <- read_csv(jaw_files, id = "file", show_col_types = FALSE) %>%
  mutate(
    fish = gsub("JAW_(.*)_.*\\.csv", "\\1", basename(file)),
    species = gsub("(\\D+)\\d+$", "\\1", fish),
    trial = gsub("JAW_.*_(.*)\\.csv", "\\1", basename(file))
  ) %>%
  group_by(file) %>%
  filter(ms < ms[which.max(g)]) %>%
  mutate(
    g = g - first(g),
    n = row_number(),
    ang = n / max(n) * 30,
    ang = plyr::round_any(ang, 0.5)
  ) %>%
  ungroup() %>%
  group_by(fish, species, trial, ang) %>%
  summarise(g = mean(g), .groups = "drop") %>%
  arrange(ang, fish, species, trial) %>%
  left_join(name_key, by = "species") %>%
  left_join(moments, by = "fish") %>%
  mutate(
    specimen = stringr::str_extract(fish, "\\d+$"),
    f = g / 1000 * 9.81,
    moment_length_m = moment_length / 100,
    theta_rad = ang * pi / 180,
    torque = f * moment_length_m * cos(theta_rad),
    name = factor(name, levels = name_levels),
    fish = factor(fish),
    specimen = factor(specimen)
  ) %>% filter(ang <= 20)
# Optional check
dat %>% filter(fish == "LB4", trial == "01")

# ------------------------------------------------------------
# Trial-level angular stiffness
# ------------------------------------------------------------

stiff_dat <- dat %>%
  mutate(
    position = case_when(
      ang < 5 ~ "low",
      ang >= 5 & ang <= 20 ~ "mid"
    )
  ) %>%
  filter(!is.na(position)) %>%
  group_by(fish, name, trial, specimen, position) %>%
  group_modify(~{
    fit <- lm(torque ~ theta_rad, data = .x)
    tibble(
      stiff = unname(coef(fit)[["theta_rad"]]),
      r2 = summary(fit)$r.squared
    )
  }) %>%
  ungroup() %>%
  mutate(
    name = factor(name, levels = name_levels),
    fish = factor(fish),
    specimen = factor(specimen),
    position = factor(position, levels = c("low", "mid"))
  )

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

ex_dat <- dat %>% 
  filter(fish=="LB3", trial=="01") %>% 
  rowwise() %>%   
  mutate(low=ang<5,
         mid=ang>=5&ang<=20,
         position=ifelse(low,"low",NA),
         position=ifelse(mid,"mid",position))

ex_fit_l <- lm(torque~theta_rad,ex_dat %>% filter(position=="low"))
ex_fit_m <- lm(torque~theta_rad,ex_dat %>% filter(position=="mid"))

ex_pred_l <- tibble(theta_rad=seq(0,5*pi/180,0.001),torque=predict(ex_fit_l,newdata = data.frame(theta_rad=seq(0,5*pi/180,0.001))))
ex_pred_m <- tibble(theta_rad=seq(5*pi/180,20*pi/180,0.001),torque=predict(ex_fit_m,newdata = data.frame(theta_rad=seq(5*pi/180,20*pi/180,0.001))))

p_ex <- ex_dat%>% 
  ggplot(aes(theta_rad,torque)) +
  geom_line(data=ex_pred_l,linewidth = 5,alpha=0.6,col="red")+
  geom_line(data=ex_pred_m,linewidth = 5,alpha=0.6,col="blue")+
  geom_line(col="gray30",linewidth = 3)+
  theme_classic(24) + 
  labs(
    y = "Torque (N m)",
    x = "Angle (rad.)"
  )

p_ex

ggsave("manuscript/stiffexamp.pdf", plot = p_ex, height = 8, width = 6)


p_ex+theme(axis.text = element_blank())

ggsave("manuscript/stiffexamp2.pdf", height = 4, width = 5)


rat_median <- stiff_dat %>%
  filter(name == "Rat") %>%
  group_by(position) %>% 
  summarise(med = median(stiff, na.rm = TRUE)) %>% 
  mutate(
    position = recode(position,
                      low = "Initial",
                      mid = "Mid"
    ),
    position = factor(position, levels = c("Initial", "Mid"))
  ) 
p_stiff <- stiff_dat %>%
  mutate(
    position = recode(position,
                      low = "Initial",
                      mid = "Mid"
    ),
    position = factor(position, levels = c("Initial", "Mid"))
  ) %>%
  ggplot(aes(name, stiff)) +
  geom_hline(data = rat_median, aes(yintercept = med), linetype = 2) +
  geom_boxplot() +
  theme_few(20) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    axis.title.y = element_text(size = 20)
  ) +
  labs(
    x = "",
    y = expression("Angular stiffness (N m " ~ rad^{-1} * ")")
  ) +
  facet_wrap(. ~ position, ncol = 1, scales = "free_y")
p_stiff

ggsave("manuscript/angstiff.pdf", plot = p_stiff, height = 6, width = 8)

# ------------------------------------------------------------
# Dunnett contrasts
# ------------------------------------------------------------

stiff_dat <- stiff_dat %>%
  mutate(
    name = relevel(factor(name), ref = "Rat"),
    fish = factor(fish)
  )

stiff_dat_sum <- stiff_dat%>% 
  group_by(name,position) %>% 
  summarise(m=mean(stiff),se=sd(stiff)/sqrt(length(stiff)))
mod <- lmer(stiff~ position*name + (1 | fish), data = stiff_dat )
car::Anova(mod)

mod_l <- lmer(stiff ~ name + (1 | fish), data = stiff_dat %>% filter(position=="low"))
dun_l <- glht(mod_l, linfct = mcp(name = "Dunnett"))
summary(dun_l)

mod_m <- lmer(stiff ~ name + (1 | fish), data = stiff_dat %>% filter(position=="mid"))
dun_m <- glht(mod_m, linfct = mcp(name = "Dunnett"))
summary(dun_m)


library(emmeans)

emm <- emmeans(mod, ~ name | position)
pairs(emm, adjust = "dunnett", ref = "Rat")
