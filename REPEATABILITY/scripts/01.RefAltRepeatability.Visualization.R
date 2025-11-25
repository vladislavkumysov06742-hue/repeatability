
## Install if needed:
## install.packages(c("ggplot2", "dplyr", "scales"))

## Подключаем библиотеки
library(ggplot2)  # графики
library(dplyr)    # манипуляции с таблицами
library(scales)   # вспомогательные функции для визуализации

## Параметры для круга генома
MTDNA_LENGTH <- 16568  # ожидаемая длина человеческой митохондриальной ДНК
radius <- 1            # радиус круга для отрисовки (произвольная единица)

## Загружаем результаты поиска повторов
## ЗАМЕНА ПУТИ: в оригинальном скрипте путь был абсолютный и специфичен для машины разработчика.
## Здесь использован относительный путь, подходящий для этого репозитория.
repeats <- read.csv("../data/2_derived/pos8473_TtoC/my_mtDNA_repeat_all_repeats.csv")

## Фильтруем по длине EffectiveLength — оставляем только 15 и 16 (пример визуализации)
repeats <- repeats %>% filter(EffectiveLength %in% c(16, 15))

if (nrow(repeats) == 0) stop("No 15 or 16 bp repeats found in your file!")

## angle_for_pos переводит позицию в угол (радианы) на круге, где позиция 1 располагается сверху
angle_for_pos <- function(pos) (0.5 - (pos-1)/MTDNA_LENGTH) * 2 * pi

## Добавляем вспомогательные координаты: углы и x/y для позиции и центра повтора
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

## Подготавливаем траекторию круга
n_points <- 400
circle_path <- data.frame(x = cos(seq(0, 2*pi, length.out=n_points)) * radius,
                          y = sin(seq(0, 2*pi, length.out=n_points)) * radius)

## Отмечаем позицию 1 на круге
angle_1 <- angle_for_pos(1)
x1 <- cos(angle_1) * radius
y1 <- sin(angle_1) * radius

## Координаты варианта (предполагаем, что во всех строках repeats одна и та же pos)
variant_angle <- angle_for_pos(repeats$pos[1])
x_pos <- cos(variant_angle) * radius
y_pos <- sin(variant_angle) * radius

## Функция для построения квадратичной кривой (Bézier) между двумя точками
make_arc_df <- function(x0, y0, x1, y1, arc_color, arc_thick, arc_type, arc_label, n=50, curve=0.37) {
  # xm, ym — середина отрезка; затем строим контрольную точку для кривой
  xm <- (x0 + x1)/2
  ym <- (y0 + y1)/2
  dx <- x1 - x0; dy <- y1 - y0
  # сдвиг контрольной точки в наружную сторону круга для более "пышной" дуги
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

## Собираем все дуги в один датафрейм
arcs_df <- bind_rows(lapply(1:nrow(repeats), function(i) {
  make_arc_df(
    repeats$x_pos[i], repeats$y_pos[i], repeats$x_rep[i], repeats$y_rep[i],
    arc_color = repeats$arc_color[i], arc_thick = repeats$arc_thick[i],
    arc_type = repeats$arc_type[i], arc_label = repeats$arc_label[i])
}))

## Рисуем итоговую фигуру
ggplot() +
  geom_path(data = circle_path, aes(x = x, y = y), color = "gray30", linewidth = 2) +
  geom_point(aes(x=x1, y=y1), color="black", size=5) +
  geom_text(aes(x=x1, y=y1), label = "1", vjust = -1.3, fontface = "bold", size = 5) +
  geom_point(aes(x=x_pos, y=y_pos), color="goldenrod", size=7) +
  geom_text(aes(x=x_pos, y=y_pos), label = unique(repeats$pos), nudge_y=0.11, fontface="bold", size=5.1, color="black") +
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
