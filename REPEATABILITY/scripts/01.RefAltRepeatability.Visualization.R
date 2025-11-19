
# Install if needed:
# install.packages(c("ggplot2", "dplyr", "scales"))

library(ggplot2)
library(dplyr)
library(scales)

# Genome and parameters
MTDNA_LENGTH <- 16568
radius <- 1

# Load your repeat data (edit path as needed)
repeats <- read.csv("/home/popadin/Documents/REPEATABILITY/data/2_derived/pos8473_TtoC/my_mtDNA_repeat_all_repeats.csv")

# Filter: keep only 16 bp and 15 bp repeats
repeats <- repeats %>% filter(EffectiveLength %in% c(16, 15))

if (nrow(repeats) == 0) stop("No 15 or 16 bp repeats found in your file!")

# Helper: mtDNA position (1 at top, clockwise)
angle_for_pos <- function(pos) (0.5 - (pos-1)/MTDNA_LENGTH) * 2 * pi

repeats <- repeats %>%
  mutate(
    angle_pos = angle_for_pos(pos),
    angle_rep = angle_for_pos((repeat.start + repeat.end) / 2),
    x_pos = cos(angle_pos) * radius,
    y_pos = sin(angle_pos) * radius,
    x_rep = cos(angle_rep) * radius,
    y_rep = sin(angle_rep) * radius,
    arc_color = ifelse(RefAlt == "Ref", "#3182bd", "#de2d26"),
    arc_thick = ifelse(EffectiveLength == 16, 7, 4),
    arc_type  = case_when(
      RefAlt == "Ref" & EffectiveLength == 16 ~ "Ref, 16bp",
      RefAlt == "Alt" & EffectiveLength == 16 ~ "Alt, 16bp",
      RefAlt == "Ref" & EffectiveLength == 15 ~ "Ref, 15bp",
      RefAlt == "Alt" & EffectiveLength == 15 ~ "Alt, 15bp"
    ),
    arc_label = paste0(arc_type, "<br>repeat: ", repeat.start, "-", repeat.end)
  )

# Genome circle
n_points <- 400
circle_path <- data.frame(x = cos(seq(0, 2*pi, length.out=n_points)) * radius,
                          y = sin(seq(0, 2*pi, length.out=n_points)) * radius)

# Where is position 1?
angle_1 <- angle_for_pos(1)
x1 <- cos(angle_1) * radius
y1 <- sin(angle_1) * radius

# Where is "pos"?
variant_angle <- angle_for_pos(repeats$pos[1])
x_pos <- cos(variant_angle) * radius
y_pos <- sin(variant_angle) * radius

# Create arcs as quadratic Béziers (visually pleasing and simple)
make_arc_df <- function(x0, y0, x1, y1, arc_color, arc_thick, arc_type, arc_label, n=50, curve=0.37) {
  xm <- (x0 + x1)/2
  ym <- (y0 + y1)/2
  # Outward control point
  dx <- x1 - x0; dy <- y1 - y0
  ctrl_x <- xm - dy * curve
  ctrl_y <- ym + dx * curve
  t <- seq(0, 1, length.out = n)
  data.frame(
    x = (1 - t)^2 * x0 + 2 * (1 - t) * t * ctrl_x + t^2 * x1,
    y = (1 - t)^2 * y0 + 2 * (1 - t) * t * ctrl_y + t^2 * y1,
    arc_color = arc_color, arc_thick = arc_thick, arc_type = arc_type, arc_label = arc_label,
    group = paste0(x0, "_", y0, "_", x1, "_", y1, "_", arc_color)
  )
}
arcs_df <- bind_rows(lapply(1:nrow(repeats), function(i) {
  make_arc_df(
    repeats$x_pos[i], repeats$y_pos[i], repeats$x_rep[i], repeats$y_rep[i],
    arc_color = repeats$arc_color[i], arc_thick = repeats$arc_thick[i],
    arc_type = repeats$arc_type[i], arc_label = repeats$arc_label[i])
}))

# Final Plot
ggplot() +
  # mtDNA circle
  geom_path(data = circle_path, aes(x = x, y = y), color = "gray30", linewidth = 2) +
  # Position 1
  geom_point(aes(x=x1, y=y1), color="black", size=5) +
  geom_text(aes(x=x1, y=y1), label = "1", vjust = -1.3, fontface = "bold", size = 5) +
  # Variant pos (big golden dot and position number)
  geom_point(aes(x=x_pos, y=y_pos), color="goldenrod", size=7) +
  geom_text(aes(x=x_pos, y=y_pos), label = unique(repeats$pos), nudge_y=0.11, fontface="bold", size=5.1, color="black") +
  # Repeat arcs
  geom_path(
    data = arcs_df,
    aes(x = x, y = y, group = group, color = arc_type, linewidth = arc_thick),
    lineend = "round", alpha = 0.82, show.legend = TRUE
  ) +
  scale_color_manual(
    values = c("Ref, 16bp" = "#3182bd", "Alt, 16bp" = "#de2d26",
               "Ref, 15bp" = "#6baed6", "Alt, 15bp" = "#fc9272"),
    name = "Repeat type"
  ) +
  scale_linewidth_identity() +
  guides(color = guide_legend(override.aes = list(linewidth = c(7,7,4,4)), ncol=1)) +
  coord_equal(expand=FALSE) +
  theme_void(base_size = 16) +
  labs(
    title = "mtDNA: 15/16-bp Repeats Connected To Variant Position",
    subtitle = "Genome circle 1–16568 (top=1). Blue=Ref, Red=Alt. Thickness: Repeat Length. Gold dot marks variant.",
    x=NULL, y=NULL
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )
