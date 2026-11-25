# 13_figures.R — 正文图 1–6 与补充图 1–4（ggplot2）| Article Figs 1–6 and Supplementary Figs 1–4 (ggplot2)
#
# 需要 01–02、07（12 个配置）、08、09、11 与 12 的结果；本脚本只读取结果，不运行模型。
# Needs the results of steps 01–02, 07 (12 configurations), 08, 09, 11 and 12; reads results only and runs no model.
# 版式：每个面板有完整外框与灰色标题条，图例在面板下方；分类色为深蓝、深红、赭黄、紫、绿，连续色用 ColorBrewer 色带。
# Style: each panel has a full frame and a grey strip naming it, and one legend sits below the panels; categories
# use navy, deep red, ochre, purple and green, and continuous scales use ColorBrewer palettes.
# 图 1a 的 EEZ 外框为可选项：设定 EEZ_SHP 指向 Marine Regions EEZ v12（0–360 经度版）并安装 sf。
# EEZ outlines in Fig. 1a are optional: set EEZ_SHP to the Marine Regions EEZ v12 shapefile (0–360 longitudes)
# and install sf.
#   EEZ_SHP=/path/to/eez_v12_0_360.shp Rscript R/13_figures.R

source("R/00_settings.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(ggtext); library(grid) })

# ---- 0. 颜色、主题与图例 | colours, theme and legends ----
NAVY <- "#1F5673"; RED <- "#A31515"; OCHRE <- "#C98A00"; PURPLE <- "#6B4C8A"; GREEN <- "#2E7D4F"
STRIP <- "grey85"; BAND <- "grey95"; DARKGREY <- "#4D4D4D"
MP_LONG <- c("MP-A WCPFC TAC", "MP-B joint TAC", "MP-C TAC with deduction", "MP-D overlap effort cap")
MP <- c("MP-A", "MP-B", "MP-C", "MP-D")
MP_COL <- setNames(c(NAVY, RED, OCHRE, PURPLE), MP)
MP_SHAPE <- setNames(c(16, 15, 17, 18), MP)        # 圆、方、三角、菱形 | circle, square, triangle, diamond
SB <- "SB/SB<sub>F=0</sub>"                        # 标题条与轴标题按 Markdown 显示下标 | strips and axis titles render the subscript
FONT <- "Arial"

theme_paper <- theme_bw(base_size = 8, base_family = FONT) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        strip.background = element_rect(fill = STRIP, colour = "black", linewidth = 0.6),
        strip.text = element_markdown(face = "bold", size = 8.5, margin = margin(3.5, 2, 3.5, 2)),
        axis.text = element_text(colour = "black", size = 7.5),
        axis.title.x = element_markdown(size = 8.5), axis.title.y = element_markdown(size = 8.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.4),
        legend.position = "none",
        plot.tag = element_text(face = "bold", size = 10),
        plot.margin = margin(3, 5, 3, 3))
# 每个面板的灰色标题条用单格分面画出：geom_blank 提供分面变量 lab = "x"，labeller 把它换成面板名称。
# Each panel's grey strip is a one-panel facet: geom_blank supplies the facet variable lab = "x", and the labeller
# turns it into the panel name.

# 底部图例用 grid 画，线、点与色块可放在同一行；先列出每张图的图例项，再在一个循环里画好。
# Bottom legends are drawn with grid so that lines, points and patches share one row; the keys of every figure
# are listed first and drawn in one loop.
key_mp <- list(list(label = "MP-A", col = NAVY, fill = NAVY, lty = 1, pch = 16), list(label = "MP-B", col = RED, fill = RED, lty = 1, pch = 15),
               list(label = "MP-C", col = OCHRE, fill = OCHRE, lty = 1, pch = 17), list(label = "MP-D", col = PURPLE, fill = PURPLE, lty = 1, pch = 18))
legend_keys <- list(
  Fig1 = list(list(label = "WCPFC", col = NAVY, lty = 1), list(label = "IATTC", col = RED, lty = 2),
              list(label = "Hooks", fill = PURPLE, rect = TRUE), list(label = "Albacore weight", fill = OCHRE, rect = TRUE),
              list(label = "Albacore number", fill = GREEN, rect = TRUE)),
  Fig2 = list(list(label = "Catch scalar", col = NAVY, lty = 1), list(label = "Control points of CMM 2025-01", col = "grey20", pch = 16),
              list(label = "Target reference point", col = "black", lty = 2, lwd = 1.1), list(label = "Implementation factor", fill = "white", rect = TRUE)),
  Fig3 = list(list(label = "MP-A", col = NAVY, lty = 1, lwd = 2.2), list(label = "MP-B", col = RED, lty = 1, lwd = 2.2),
              list(label = "Target reference point", col = "black", lty = 3, lwd = 1.5), list(label = "Quota", col = "black", lty = 2, lwd = 1.3)),
  FigS1 = list(list(label = "MP-C", col = OCHRE, lty = 1, lwd = 2.2), list(label = "MP-D", col = PURPLE, lty = 1, lwd = 2.2),
               list(label = "Target reference point", col = "black", lty = 3, lwd = 1.5), list(label = "Quota", col = "black", lty = 2, lwd = 1.3)),
  Fig5 = list(list(label = "Shares updated from reported catch", col = RED, fill = RED, lty = 1, pch = 15),
              list(label = "Fixed shares", col = "black", fill = "white", pch = 22),
              list(label = "Fleet W", fill = DARKGREY, rect = TRUE), list(label = "Fleet OW", fill = STRIP, rect = TRUE)),
  MP = key_mp)                                       # 图 4、图 6、补充图 3–4 | Figs 4 and 6, Supplementary Figs 3–4
sym <- 0.30; gap <- 0.07; space <- 0.24              # 符号宽、符号与文字的间距、两项的间距（英寸）| symbol width and gaps (inches)
gp <- gpar(fontsize = 7.5, fontfamily = FONT)
measure <- tempfile(fileext = ".png")
ragg::agg_png(measure, width = 7.2, height = 1, units = "in", res = 72)   # 只用来量文字宽度 | only for measuring text widths
legends <- list()
for (fig in names(legend_keys)) {
  keys <- legend_keys[[fig]]
  text_w <- numeric(length(keys))
  for (i in seq_along(keys)) text_w[i] <- convertWidth(grobWidth(textGrob(keys[[i]]$label, gp = gp)), "in", valueOnly = TRUE)
  x <- -(sum(sym + gap + text_w) + space * (length(keys) - 1)) / 2   # 整行置中 | centre the row
  y <- unit(0.5, "npc")
  grobs <- list()
  for (i in seq_along(keys)) {
    k <- keys[[i]]
    x0 <- unit(0.5, "npc") + unit(x, "in"); x1 <- unit(0.5, "npc") + unit(x + sym, "in"); xm <- unit(0.5, "npc") + unit(x + sym / 2, "in")
    if (!is.null(k$lty)) grobs[[length(grobs) + 1]] <- segmentsGrob(x0, y, x1, y, gp = gpar(col = k$col, lty = k$lty, lwd = if (is.null(k$lwd)) 1.8 else k$lwd, lineend = "butt"))
    if (!is.null(k$pch)) grobs[[length(grobs) + 1]] <- pointsGrob(xm, y, pch = k$pch, size = unit(5, "pt"), gp = gpar(col = k$col, fill = k$fill))
    if (!is.null(k$rect)) grobs[[length(grobs) + 1]] <- rectGrob(xm, y, width = unit(0.17, "in"), height = unit(0.115, "in"), gp = gpar(fill = k$fill, col = "black", lwd = 0.8))
    grobs[[length(grobs) + 1]] <- textGrob(k$label, x = unit(0.5, "npc") + unit(x + sym + gap, "in"), y = y, just = "left", gp = gp)
    x <- x + sym + gap + text_w[i] + space
  }
  legends[[fig]] <- wrap_elements(full = gTree(children = do.call(gList, grobs)))
}
invisible(dev.off()); unlink(measure)

# ================================================================ 图 1 | Figure 1
# a) 2015–2024 每个 5° 格的钩数与两个公约区边界；b、c) 重疊区份额；d) 五个船旗两组织数据的比值
# a) hooks per 5° cell in 2015–2024 with the two convention-area boundaries; b, c) overlap shares; d) ratios of the
#    two datasets for five flags
public <- read_csv(file.path(OUT, "02_overlap_share_by_year.csv"), show_col_types = FALSE)
flags <- read_csv(file.path(OUT, "02_flag_ratio_2015_2024.csv"), show_col_types = FALSE)
cells <- read_csv(file.path(OUT, "01_ll_long.csv"), col_select = c(rfmo, year, lat_sw, lon_sw, hooks), show_col_types = FALSE) |>
  filter(year >= 2015, year <= 2024) |> group_by(rfmo, lat_sw, lon_sw) |> summarise(hooks = sum(hooks), .groups = "drop") |>
  mutate(lon = lon_sw %% 360) |> select(rfmo, lat = lat_sw, lon, hooks) |>
  pivot_wider(names_from = rfmo, values_from = hooks, values_fill = 0) |> filter(WCPFC > 0 | IATTC > 0)
# 两个数据集都有的格分成两个三角形：左上 WCPFC，右下 IATTC | a cell in both datasets is split: WCPFC upper left, IATTC lower right
upper <- cbind(c(0, 0, 5), c(0, 5, 5)); lower <- cbind(c(0, 5, 5), c(0, 0, 5)); square <- cbind(c(0, 5, 5, 0), c(0, 0, 5, 5))
polygons <- list()
for (i in seq_len(nrow(cells))) {
  lat <- cells$lat[i]; lon <- cells$lon[i]; hw <- cells$WCPFC[i]; hi <- cells$IATTC[i]
  if (hw > 0 && hi > 0) {
    polygons[[length(polygons) + 1]] <- data.frame(id = paste0(i, "w"), x = lon + upper[, 1], y = lat + upper[, 2], hooks = hw / 1e6)
    polygons[[length(polygons) + 1]] <- data.frame(id = paste0(i, "i"), x = lon + lower[, 1], y = lat + lower[, 2], hooks = hi / 1e6)
  } else {
    polygons[[length(polygons) + 1]] <- data.frame(id = paste0(i, "s"), x = lon + square[, 1], y = lat + square[, 2], hooks = max(hw, hi) / 1e6)
  }
}
polygons <- bind_rows(polygons)
# EEZ 外框（可选）| EEZ outlines (optional)
outline <- NULL
eez_file <- Sys.getenv("EEZ_SHP", "")
if (nzchar(eez_file) && file.exists(eez_file) && requireNamespace("sf", quietly = TRUE)) {
  sf::sf_use_s2(FALSE)                               # 0–360 经度按平面坐标处理 | treat 0–360 longitudes as planar
  eez <- sf::st_read(eez_file, quiet = TRUE)
  box <- sapply(sf::st_geometry(eez), sf::st_bbox)
  eez <- eez[!(box["xmax", ] < 140 | box["xmin", ] > 290 | box["ymax", ] < -55 | box["ymin", ] > 5), ]   # 地图范围内 | within the map
  outline <- as.data.frame(sf::st_coordinates(sf::st_cast(sf::st_geometry(eez), "MULTILINESTRING")))
  outline$piece <- interaction(outline[, setdiff(names(outline), c("X", "Y"))], drop = TRUE)              # 每段线一组 | one group per line piece
}
# 公约区边界：WCPFC 实线、IATTC 虚线，交集是重疊区 | convention-area boundaries: WCPFC solid, IATTC dashed; the intersection is the overlap area
RFMO_COL <- c(WCPFC = NAVY, IATTC = RED); RFMO_LTY <- c(WCPFC = "solid", IATTC = "dashed")
boundary <- rbind(data.frame(rfmo = "WCPFC", x = c(210, 210, 230, 230), y = c(0, -4, -4, -55)),
                  data.frame(rfmo = "IATTC", x = c(210, 210, 290), y = c(5, -50, -50)))
p1a <- ggplot() + geom_polygon(data = polygons, aes(x, y, group = id, fill = hooks))
if (!is.null(outline)) p1a <- p1a + geom_path(data = outline, aes(X, Y, group = piece), colour = "grey55", linewidth = 0.12)
p1a <- p1a +
  geom_path(data = boundary, aes(x, y, colour = rfmo, linetype = rfmo), linewidth = 0.85) +
  annotate("text", x = 220, y = -27, label = "Overlap\narea", fontface = "bold", size = 2.7, family = FONT, lineheight = 0.95) +
  scale_fill_distiller(palette = "Blues", direction = 1, trans = "sqrt", limits = c(0, NA), name = "Hooks (millions)",   # 平方根色标 | square-root scale
                       guide = guide_colourbar(barheight = unit(1.75, "in"), barwidth = unit(0.11, "in"), frame.colour = "black", ticks.colour = "black",
                                               title.position = "right", title.theme = element_text(angle = 90, size = 8.5, hjust = 0.5, family = FONT))) +
  scale_colour_manual(values = RFMO_COL, guide = "none") + scale_linetype_manual(values = RFMO_LTY, guide = "none") +
  scale_x_continuous(breaks = c(150, 180, 210, 240, 270), labels = c("150°E", "180°", "150°W", "120°W", "90°W"), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(-50, -25, 0), labels = c("50°S", "25°S", "0°"), expand = c(0, 0)) +
  coord_fixed(xlim = c(140, 290), ylim = c(-55, 5)) + labs(x = NULL, y = NULL, tag = "a)") +
  geom_blank(data = data.frame(lab = "x"), inherit.aes = FALSE) + facet_wrap(~lab, labeller = as_labeller(c(x = "Longline hooks, 2015–2024"))) +
  theme_paper + theme(legend.position = "right", legend.margin = margin(0, 0, 0, 2))
# b、c) 面积加权的重疊区份额 | b, c) area-weighted overlap shares
p1bc <- list()
for (j in 1:2) {
  d <- public |> filter(scheme == "area") |> mutate(value = 100 * .data[[c("share_hooks", "share_albmt")[j]]])
  p1bc[[j]] <- ggplot(d, aes(year, value, colour = rfmo, linetype = rfmo)) + geom_line(linewidth = 0.8) +
    scale_colour_manual(values = RFMO_COL) + scale_linetype_manual(values = RFMO_LTY) +
    scale_x_continuous(breaks = c(1990, 2000, 2010, 2020), limits = c(1990, 2024), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) + labs(x = "Year", y = "Overlap share (%)", tag = c("b)", "c)")[j]) +
    geom_blank(data = data.frame(lab = "x"), inherit.aes = FALSE) + facet_wrap(~lab, labeller = as_labeller(c(x = c("Hooks", "Albacore weight")[j]))) +
    theme_paper
}
# d) 十年在两个数据集中都有钩数的五个船旗：IATTC / WCPFC | d) five flags with hooks in both datasets in all ten years: IATTC / WCPFC
five <- c("Chinese Taipei", "French Polynesia", "China", "Korea", "Vanuatu")
QUANTITY <- c(hooks = "Hooks", alb_mt = "Albacore weight", alb_n = "Albacore number")
ratio <- flags |> filter(flag_label %in% five, !is.na(ratio)) |>
  mutate(row = match(flag_label, five), offset = c(hooks = -0.26, alb_mt = 0, alb_n = 0.26)[metric], quantity = factor(QUANTITY[metric], QUANTITY))
p1d <- ggplot(ratio) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = seq(1, 5, 2) - 0.46, ymax = seq(1, 5, 2) + 0.46, fill = BAND) +   # 交替底色只为对齐各行 | bands only align rows
  geom_rect(aes(xmin = pmin(1, ratio), xmax = pmax(1, ratio), ymin = row + offset - 0.13, ymax = row + offset + 0.13, fill = quantity), colour = "black", linewidth = 0.25) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4) +
  scale_fill_manual(values = setNames(c(PURPLE, OCHRE, GREEN), QUANTITY)) +
  scale_x_continuous(breaks = seq(0.8, 1.8, 0.2), limits = c(0.8, 1.9), expand = c(0, 0)) +
  scale_y_reverse(breaks = 1:5, labels = five, limits = c(5.5, 0.5), expand = c(0, 0)) +
  labs(x = "Ratio of IATTC to WCPFC totals in the overlap area", y = NULL, tag = "d)") +
  geom_blank(data = data.frame(lab = "x"), inherit.aes = FALSE) + facet_wrap(~lab, labeller = as_labeller(c(x = "Five flags, 2015–2024"))) +
  theme_paper + theme(axis.ticks.y = element_blank())
fig1 <- p1a / (p1bc[[1]] | p1bc[[2]]) / p1d / legends$Fig1 + plot_layout(heights = c(2.45, 1.5, 1.75, 0.16))
ggsave(file.path(FIG, "Fig1_public_data.png"), fig1, width = 7.2, height = 7.7, dpi = 450, device = ragg::agg_png, bg = "white")
ggsave(file.path(FIG, "Fig1_public_data.pdf"), fig1, width = 7.2, height = 7.7, device = cairo_pdf, bg = "white")

# ================================================================ 图 2 | Figure 2
# a) 三个区域与四支船队；b) 每年的闭环；c) 捕捞控制规则 | a) regions and fleets; b) annual closed loop; c) harvest control rule
cards <- rbind(data.frame(x = c(0.8, 3.9, 7.0), y = 5.5, w = 2.25, head = c("R1  WCPFC only", "R2  Overlap area", "R3  EPO only"), body = c("Fleet W", "Fleets OW and OI", "Fleet E")),
               data.frame(x = c(0.6, 6.5, 6.5, 0.6), y = c(3.2, 3.2, 0.6, 0.6), w = 2.9,
                          head = c("Operating model", "Observations", "Estimation method", "Management and fishing"),
                          body = c("Survival, recruitment and catch", "Reported catch and abundance indices", "Age-structured population update", "HCR, TAC, allocation, catch")))
band <- 0.42                                         # 方框顶部灰条的高度 | height of the grey band on top of each box
arrows <- data.frame(x = c(3.5, 7.95, 6.5, 2.05), y = c(3.7, 3.2, 1.1, 1.6), xend = c(6.5, 7.95, 3.5, 2.05), yend = c(3.7, 1.6, 1.1, 3.2))
links <- data.frame(x = c(3.05, 6.15), xend = c(3.9, 7.0), y = 5.85)          # 区域之间的鱼群移动 | movement between regions
tags <- data.frame(x = c(3.475, 6.575, 4.3, 7.4, 5.0), y = c(6.2, 6.2, 6.78, 6.78, 4.02), label = c("S3", "S3", "S1", "S4", "S2"))   # 实施因子 | implementation factors
notes <- data.frame(x = c(5, 8.1, 5, 2.2), y = c(3.47, 2.4, 0.85, 2.4), label = c("Catch and indices", "Reported catch\nand indices", "Estimated stock status", "Removals"), hjust = c(0.5, 0, 0.5, 0))
p2ab <- ggplot() +
  geom_rect(data = cards, aes(xmin = x, xmax = x + w, ymin = y, ymax = y + 1), fill = "white", colour = "black", linewidth = 0.45) +
  geom_rect(data = cards, aes(xmin = x, xmax = x + w, ymin = y + 1 - band, ymax = y + 1), fill = STRIP, colour = "black", linewidth = 0.45) +
  geom_segment(data = cards, aes(x = x, xend = x + w, y = y + 1 - band, yend = y + 1 - band), linewidth = 0.8) +
  geom_text(data = cards, aes(x + w / 2, y + 1 - band / 2, label = head), fontface = "bold", size = 3.0, family = FONT) +
  geom_text(data = cards, aes(x + w / 2, y + (1 - band) / 2, label = body), size = 2.65, family = FONT) +
  geom_segment(data = arrows, aes(x, y, xend = xend, yend = yend), linewidth = 0.45, arrow = arrow(length = unit(1.8, "mm"), type = "closed")) +
  geom_segment(data = links, aes(x, y, xend = xend, yend = y), linewidth = 0.45, arrow = arrow(length = unit(1.6, "mm"), ends = "both", type = "closed")) +
  geom_label(data = tags, aes(x, y, label = label), fill = "white", colour = "black", fontface = "bold", size = 2.65, family = FONT,
             label.padding = unit(0.16, "lines"), label.r = unit(0, "lines"), linewidth = 0.3) +
  geom_text(data = notes, aes(x, y, label = label, hjust = hjust), size = 2.65, family = FONT, lineheight = 0.95) +
  annotate("text", x = 0, y = c(6.85, 4.7), label = c("a)", "b)"), fontface = "bold", size = 3.55, family = FONT, hjust = 0) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0.4, 7), expand = FALSE) + theme_void()
# HCR（CMM 2025-01）：0.37 以下 0.20；0.37–0.94 曲率 10 的指数段；0.94–1.29 为 1.09；1.59 时 1.29
# HCR (CMM 2025-01): 0.20 below 0.37; exponential section with curvature 10 from 0.37 to 0.94; 1.09 from 0.94 to 1.29; 1.29 at 1.59
curve <- data.frame(status = seq(0, 2, length.out = 801))
curve$scalar <- ifelse(curve$status < 0.37, 0.20,
                ifelse(curve$status < 0.94, 0.20 + 0.89 * (1 - exp(-10 * (curve$status - 0.37))) / (1 - exp(-5.7)),
                ifelse(curve$status <= 1.29, 1.09,
                ifelse(curve$status < 1.59, 1.09 + 0.20 * (curve$status - 1.29) / 0.30, 1.29))))
control <- data.frame(status = c(0.37, 0.94, 1.29, 1.59), scalar = c(0.20, 1.09, 1.09, 1.29))
p2c <- ggplot(curve, aes(status, scalar)) + geom_vline(xintercept = 0.96, linetype = "dashed", linewidth = 0.4) +
  geom_line(colour = NAVY, linewidth = 0.85) + geom_point(data = control, colour = "grey20", size = 1.9) +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) + scale_y_continuous(limits = c(0, 1.45), expand = c(0, 0)) +
  labs(x = "Relative stock status", y = "Catch scalar", tag = "c)") +
  geom_blank(data = data.frame(lab = "x"), inherit.aes = FALSE) + facet_wrap(~lab, labeller = as_labeller(c(x = "Harvest control rule"))) +
  theme_paper
fig2 <- p2ab / (plot_spacer() | p2c | plot_spacer()) / legends$Fig2 + plot_layout(heights = c(3.75, 1.5, 0.16)) &
  theme(plot.background = element_rect(fill = "white", colour = NA))
fig2[[2]] <- fig2[[2]] + plot_layout(widths = c(0.55, 1, 0.55))   # HCR 面板置中 | centre the HCR panel
ggsave(file.path(FIG, "Fig2_design.png"), fig2, width = 7.2, height = 6.0, dpi = 450, device = ragg::agg_png, bg = "white")
ggsave(file.path(FIG, "Fig2_design.pdf"), fig2, width = 7.2, height = 6.0, device = cairo_pdf, bg = "white")

# ================================================================ 图 3 与补充图 1 | Figure 3 and Supplementary Figure 1
# 参考情景 C05 第 1 个重复：SB/SB_F=0 与目标参考点，W、OW 船队的渔获与配额；该情景是确定性的
# Replicate 1 of the reference scenario C05: SB/SB_F=0 with the target reference point, and catch and quota of
# fleets W and OW; the scenario is deterministic
traj <- read_csv(file.path(OUT, "07_reference", "first_replicate_trajectories.csv"), show_col_types = FALSE) |> filter(case_id == "C05")
for (fig in c("Fig3", "FigS1")) {
  idx <- if (fig == "Fig3") 1:2 else 3:4             # 图 3：MP-A、MP-B；补充图 1：MP-C、MP-D | Fig. 3: MP-A, MP-B; Supp. Fig. 1: MP-C, MP-D
  panels <- list()
  for (i in idx) {
    d <- traj |> filter(mp == MP_LONG[i]) |> arrange(year)
    xlab <- if (i == idx[2]) "Projection year" else NULL                  # 只有下一行写 x 轴标题 | x title on the lower row only
    k <- (match(i, idx) - 1) * 3                                          # 该行第一个面板之前的面板数 | panels before this row
    panels[[k + 1]] <- ggplot(d, aes(year)) + geom_line(aes(y = 100 * target_D), linetype = "dotted", linewidth = 0.7) +
      geom_line(aes(y = 100 * D), colour = MP_COL[[i]], linewidth = 1) +
      scale_x_continuous(breaks = c(1, 5, 10, 15, 20), limits = c(1, 20), expand = c(0, 0)) +
      scale_y_continuous(limits = c(40.2, 46), expand = c(0, 0)) +
      labs(x = xlab, y = paste0(SB, " (%)"), tag = paste0(letters[k + 1], ")")) +
      geom_blank(data = data.frame(lab = "x"), inherit.aes = FALSE) + facet_wrap(~lab, labeller = as_labeller(c(x = paste0(MP[i], ", ", SB)))) +
      theme_paper
    for (j in 1:2) {
      fleet <- c("W", "OW")[j]
      # 实线 = 渔获，虚线阶梯 = 配额，都以历史末年总渔获的百分比表示 | solid = catch, dashed steps = quota, % of final historical catch
      panels[[k + 1 + j]] <- ggplot(d, aes(year)) +
        geom_line(aes(y = 100 * .data[[paste0("C_", fleet)]] / last_total_catch), colour = MP_COL[[i]], linewidth = 1) +
        geom_step(aes(y = 100 * .data[[paste0("attempted_TAC_", fleet)]] / last_total_catch), linetype = "dashed", linewidth = 0.6) +
        scale_x_continuous(breaks = c(1, 5, 10, 15, 20), limits = c(1, 20), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, c(90, 12)[j]), expand = c(0, 0)) +
        labs(x = xlab, y = if (j == 1) "% of final historical catch" else NULL, tag = paste0(letters[k + 1 + j], ")")) +
        geom_blank(data = data.frame(lab = "x"), inherit.aes = FALSE) + facet_wrap(~lab, labeller = as_labeller(c(x = paste0(MP[i], ", fleet ", fleet)))) +
        theme_paper
    }
  }
  figure <- wrap_plots(panels, ncol = 3) / legends[[fig]] + plot_layout(heights = c(1, 0.045))
  name <- if (fig == "Fig3") "Fig3_reference" else "FigS1_reference_CD"
  ggsave(file.path(FIG, paste0(name, ".png")), figure, width = 7.2, height = 5.0, dpi = 450, device = ragg::agg_png, bg = "white")
  ggsave(file.path(FIG, paste0(name, ".pdf")), figure, width = 7.2, height = 5.0, device = cairo_pdf, bg = "white")
}

# ================================================================ 图 4 | Figure 4
# 参考配置下实施因子对 SB/SB_F=0 与渔获的效应，含无偏误差（12）| effects of the implementation factors in the reference case,
# with the unbiased error of step 12
marg <- read_csv(file.path(OUT, "08_analysis", "marginal.csv"), show_col_types = FALSE) |> filter(config_id == "reference") |>
  select(mp, contrast, metric, estimate, mcse)
unbiased <- read_csv(file.path(OUT, "12_unbiased_error", "effects.csv"), show_col_types = FALSE) |> filter(level == "unbiased_error") |>
  transmute(mp, contrast = "unbiased", metric, estimate, mcse)
factor_rows <- c("S1+" = "S1 Choice of applicable measures", "S2-" = "S2 Under-reporting", "S2+" = "S2 Over-reporting", "unbiased" = "S2 Unbiased error",
                 "S3+" = "S3 Time-varying movement", "S4+" = "S4 Catch outside the MP\n(fleet E capacity +3% per year)")
effects4 <- bind_rows(marg, unbiased) |> filter(contrast %in% names(factor_rows), metric %in% c("mean_D", "relative_catch")) |>
  mutate(row = match(contrast, names(factor_rows)), mp = MP[match(mp, MP_LONG)], estimate = 100 * estimate, mcse = 100 * mcse)
p4 <- list()
for (j in 1:2) {
  df <- filter(effects4, metric == c("mean_D", "relative_catch")[j])
  df$y <- df$row + (as.integer(factor(df$mp, MP)) - 2.5) * 0.19       # 一行内四个 MP 上下错开 | the four MPs offset within a row
  # 由零到符号的线段、符号与 ±1.96 MCSE | a line from zero to the symbol, the symbol and ±1.96 MCSE
  p4[[j]] <- ggplot(df) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = seq(1, 6, 2) - 0.46, ymax = seq(1, 6, 2) + 0.46, fill = BAND) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = estimate, y = y, yend = y, colour = mp), linewidth = 0.75, lineend = "butt") +
    geom_errorbar(aes(xmin = estimate - 1.96 * mcse, xmax = estimate + 1.96 * mcse, y = y, colour = mp), orientation = "y", width = 0.08, linewidth = 0.4) +
    geom_point(aes(estimate, y, colour = mp, shape = mp), size = 1.7) +
    scale_colour_manual(values = MP_COL) + scale_shape_manual(values = MP_SHAPE) +
    scale_y_reverse(breaks = 1:6, labels = if (j == 1) factor_rows else NULL, limits = c(6.5, 0.5), expand = c(0, 0)) +
    labs(x = c("Change (percentage points)", "Change (% of final historical catch)")[j], y = NULL, tag = c("a)", "b)")[j]) +
    geom_blank(data = data.frame(lab = "x"), inherit.aes = FALSE) + facet_wrap(~lab, labeller = as_labeller(c(x = c(SB, "Catch")[j]))) +
    theme_paper + theme(axis.ticks.y = element_blank())
}
fig4 <- (p4[[1]] | p4[[2]]) / legends$MP + plot_layout(heights = c(1, 0.05))
ggsave(file.path(FIG, "Fig4_factor_effects.png"), fig4, width = 7.2, height = 4.3, dpi = 450, device = ragg::agg_png, bg = "white")
ggsave(file.path(FIG, "Fig4_factor_effects.pdf"), fig4, width = 7.2, height = 4.3, device = cairo_pdf, bg = "white")
write_csv(effects4 |> transmute(panel = ifelse(metric == "mean_D", "a", "b"), factor = contrast, mp, estimate, mcse), file.path(FIG, "Fig4_factor_effects_data.csv"))

# ================================================================ 图 5 | Figure 5
# MP-B，S1、S3、S4 都关：四个上报水平；实心 = 按上报渔获更新的份额（07、12），空心 = 固定份额（09）
# MP-B with S1, S3 and S4 off: four reporting levels; filled = shares updated from reported catch (steps 07, 12),
# open = fixed shares (step 09)
levels5 <- read_csv(file.path(OUT, "12_unbiased_error", "levels_other_factors_off.csv"), show_col_types = FALSE) |> filter(mp == MP_LONG[2])
fixed <- read_csv(file.path(OUT, "09_fixed_shares", "replicate_metrics.csv.gz"), show_col_types = FALSE) |> filter(S1 == 0, S3 == 0, S4 == 0) |>
  group_by(S2) |> summarise(mean_D = mean(mean_D), relative_catch = mean(relative_catch), .groups = "drop") |>
  mutate(level = c("under_reporting", "accurate", "over_reporting")[match(S2, c(-1, 0, 1))])
level_names <- c(under_reporting = "Under", accurate = "Accurate", unbiased_error = "Unbiased error", over_reporting = "Over")   # x 轴顺序 | x-axis order
level_axis <- scale_x_continuous(breaks = 1:4, labels = level_names, limits = c(0.5, 4.5), expand = c(0, 0), guide = guide_axis(n.dodge = 2))
p5 <- list()
for (j in 1:2) {
  metric5 <- c("mean_D", "relative_catch")[j]
  d <- levels5 |> filter(metric == metric5) |> transmute(x = match(level, names(level_names)), value = 100 * mean, se = 100 * mcse) |> arrange(x)
  f <- fixed |> transmute(x = match(level, names(level_names)), value = 100 * .data[[metric5]])
  span <- diff(range(d$value))
  p5[[j]] <- ggplot(d, aes(x, value)) + geom_line(colour = RED, linewidth = 0.75) +
    geom_errorbar(aes(ymin = value - 1.96 * se, ymax = value + 1.96 * se), width = 0.08, colour = RED, linewidth = 0.4) +
    geom_point(colour = RED, shape = 15, size = 2.3) + geom_point(data = f, shape = 22, fill = "white", colour = "black", size = 2.1, stroke = 0.55) +
    # 数值标在点上方，最低点标在右侧 | values above the points, the lowest to the right
    geom_text(aes(label = sprintf("%.1f", value), hjust = ifelse(value == min(value), -0.35, 0.5), vjust = ifelse(value == min(value), 0.4, -1.1)), size = 2.65, family = FONT) +
    level_axis +
    scale_y_continuous(limits = c(min(d$value, f$value) - 0.25 * span, max(d$value) + 0.38 * span), expand = c(0, 0)) +
    labs(x = "Reporting of overlap catch", y = c(paste0(SB, " (%)"), "% of final historical catch")[j], tag = c("a)", "b)")[j]) +
    geom_blank(data = data.frame(lab = "x"), inherit.aes = FALSE) + facet_wrap(~lab, labeller = as_labeller(c(x = c(SB, "Catch")[j]))) +
    theme_paper
}
# c) W、OW 船队捕满配额的年份比例；接近整数的值不写小数 | c) years in which fleets W and OW take their full quota; near-integers without decimals
quota <- levels5 |> filter(metric %in% c("quota_fully_caught_W", "quota_fully_caught_OW")) |>
  transmute(level, fleet = factor(ifelse(metric == "quota_fully_caught_W", "W", "OW"), c("W", "OW")), value = 100 * mean) |>
  mutate(x = match(level, names(level_names)) + ifelse(fleet == "W", -0.18, 0.18),
         label = ifelse(abs(value - round(value)) < 0.05, sprintf("%.0f", value), sprintf("%.1f", value))) |> arrange(x)
p5c <- ggplot(quota) + geom_rect(aes(xmin = x - 0.18, xmax = x + 0.18, ymin = 0, ymax = value, fill = fleet), colour = "black", linewidth = 0.3) +
  geom_text(aes(x, value, label = label, hjust = ifelse(fleet == "W", 0.7, 0.3)), vjust = -0.5, size = 2.5, family = FONT) +
  scale_fill_manual(values = c(W = DARKGREY, OW = STRIP)) + level_axis +
  scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 115), expand = c(0, 0)) +
  labs(x = "Reporting of overlap catch", y = "% of years 11–20", tag = "c)") +
  geom_blank(data = data.frame(lab = "x"), inherit.aes = FALSE) + facet_wrap(~lab, labeller = as_labeller(c(x = "Quota fully caught"))) +
  theme_paper
fig5 <- (p5[[1]] | p5[[2]] | p5c) / legends$Fig5 + plot_layout(heights = c(1, 0.07))
ggsave(file.path(FIG, "Fig5_reporting.png"), fig5, width = 7.2, height = 3.6, dpi = 450, device = ragg::agg_png, bg = "white")
ggsave(file.path(FIG, "Fig5_reporting.pdf"), fig5, width = 7.2, height = 3.6, device = cairo_pdf, bg = "white")
fig5_data <- levels5 |> filter(metric %in% c("mean_D", "relative_catch", "quota_fully_caught_W", "quota_fully_caught_OW")) |>
  transmute(level, shares = "updated", metric, percent = 100 * mean) |>
  bind_rows(fixed |> pivot_longer(c(mean_D, relative_catch), names_to = "metric") |> transmute(level, shares = "fixed", metric, percent = 100 * value))
write_csv(fig5_data, file.path(FIG, "Fig5_reporting_data.csv"))

# ================================================================ 图 6 | Figure 6
# MP 范围外的渔获：E 船队（重疊区以外的东太）与 T 船队（赤道至 10°S 的中西太延绳钓），容量每年 −3% 或 +3%
# Catch outside the scope of the MP: fleet E (EPO outside the overlap area) and fleet T (WCPFC longline fishery
# between the equator and 10°S), capacity −3% or +3% per year
ext_rows <- data.frame(E_level = c(-1, 1, 0, 0, -1, 1), T_level = c(0, 0, -1, 1, -1, 1),
                       label = c("Fleet E, −3% per year", "Fleet E, +3% per year", "Fleet T, −3% per year", "Fleet T, +3% per year",
                                 "Both fleets, −3% per year", "Both fleets, +3% per year"))
ext_rows$row <- seq_len(nrow(ext_rows))
ext <- read_csv(file.path(OUT, "11_external_catch", "marginal_effects.csv"), show_col_types = FALSE) |> inner_join(ext_rows, by = c("E_level", "T_level")) |>
  mutate(estimate = 100 * effect, mcse = 100 * mcse)
p6 <- list()
for (j in 1:3) {
  df <- filter(ext, metric == c("mean_D", "relative_catch_MP_fleets", "relative_catch")[j])
  df$y <- df$row + (as.integer(factor(df$mp, MP)) - 2.5) * 0.19
  p6[[j]] <- ggplot(df) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = seq(1, 6, 2) - 0.46, ymax = seq(1, 6, 2) + 0.46, fill = BAND) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = estimate, y = y, yend = y, colour = mp), linewidth = 0.75, lineend = "butt") +
    geom_errorbar(aes(xmin = estimate - 1.96 * mcse, xmax = estimate + 1.96 * mcse, y = y, colour = mp), orientation = "y", width = 0.08, linewidth = 0.4) +
    geom_point(aes(estimate, y, colour = mp, shape = mp), size = 1.7) +
    scale_colour_manual(values = MP_COL) + scale_shape_manual(values = MP_SHAPE) +
    scale_y_reverse(breaks = 1:6, labels = if (j == 1) ext_rows$label else NULL, limits = c(6.5, 0.5), expand = c(0, 0)) +
    labs(x = c("Change<br>(percentage points)", "Change<br>(% of final historical catch)", "Change<br>(% of final historical catch)")[j], y = NULL,
         tag = c("a)", "b)", "c)")[j]) +
    geom_blank(data = data.frame(lab = "x"), inherit.aes = FALSE) +
    facet_wrap(~lab, labeller = as_labeller(c(x = c(SB, "Catch of fleets W and OW", "Total catch")[j]))) +
    theme_paper + theme(axis.ticks.y = element_blank())
}
fig6 <- (p6[[1]] | p6[[2]] | p6[[3]]) / legends$MP + plot_layout(heights = c(1, 0.055))
ggsave(file.path(FIG, "Fig6_external_catch.png"), fig6, width = 7.2, height = 4.1, dpi = 450, device = ragg::agg_png, bg = "white")
ggsave(file.path(FIG, "Fig6_external_catch.pdf"), fig6, width = 7.2, height = 4.1, device = cairo_pdf, bg = "white")
write_csv(ext |> filter(metric %in% c("mean_D", "relative_catch_MP_fleets", "relative_catch")) |> transmute(metric, level = label, mp, estimate, mcse),
          file.path(FIG, "Fig6_external_catch_data.csv"))

# ================================================================ 补充图 2 | Supplementary Figure 2
# 参考配置的全部 24 个情景，相对同一 MP 的参考情景之差（配对，300 个重复的均值）
# All 24 scenarios of the reference case as differences from the reference scenario of the same MP (paired, mean of 300 replicates)
paired <- read_csv(file.path(OUT, "08_analysis", "paired.csv"), show_col_types = FALSE)
REPORTING <- c("-1" = "under-reporting", "0" = "accurate", "1" = "over-reporting")
scen <- paired |> filter(config_id == "reference", metric %in% c("mean_D", "relative_catch")) |>
  mutate(mp = factor(MP[match(mp, MP_LONG)], MP), estimate = 100 * estimate)
case_order <- scen |> distinct(case_id, S1, S2, S3, S4) |> arrange(S3, S4, S2, S1) |>
  mutate(label = paste0("S1 ", ifelse(S1 == 1, "on", "off"), ", S2 ", REPORTING[as.character(S2)], ", S3 ", ifelse(S3 == 1, "on", "off"),
                        ", S4 ", ifelse(S4 == 1, "on", "off"), ifelse(S1 == 0 & S2 == 0 & S3 == 0 & S4 == 0, " (reference scenario)", "")),
         row = row_number())
scen <- scen |> inner_join(select(case_order, case_id, row, label), by = "case_id")
pS2 <- list()
for (j in 1:2) {
  d <- filter(scen, metric == c("mean_D", "relative_catch")[j])
  lim <- max(abs(d$estimate))                        # 两个面板各自对称的色标 | a symmetric colour scale for each panel
  pS2[[j]] <- ggplot(d, aes(mp, row, fill = estimate)) + geom_tile() +
    scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(-lim, lim), name = c("Change (percentage points)", "Change (% of final historical catch)")[j],
                         guide = guide_colourbar(barwidth = unit(1.7, "in"), barheight = unit(0.09, "in"), title.position = "bottom", frame.colour = "black",
                                                 ticks.colour = "black", title.theme = element_text(size = 8.5, hjust = 0.5, family = FONT))) +
    scale_x_discrete(expand = c(0, 0)) + scale_y_reverse(breaks = case_order$row, labels = if (j == 1) case_order$label else NULL, expand = c(0, 0)) +
    labs(x = NULL, y = NULL, tag = c("a)", "b)")[j]) +
    geom_blank(data = data.frame(lab = "x"), inherit.aes = FALSE) + facet_wrap(~lab, labeller = as_labeller(c(x = c(SB, "Catch")[j]))) +
    theme_paper + theme(legend.position = "bottom", axis.ticks.y = element_blank(), legend.margin = margin(0, 0, 0, 0))
}
figS2 <- pS2[[1]] | pS2[[2]]
ggsave(file.path(FIG, "FigS2_all_scenarios.png"), figS2, width = 7.2, height = 6.4, dpi = 450, device = ragg::agg_png, bg = "white")
ggsave(file.path(FIG, "FigS2_all_scenarios.pdf"), figS2, width = 7.2, height = 6.4, device = cairo_pdf, bg = "white")

# ================================================================ 补充图 3 | Supplementary Figure 3
# 合并情景（S1 开、多报、S3 开、S4 关，C23）在 12 个配置中相对 C05 的效应；前 6 行为参考配置与生物学配置
# The combined scenario (S1 on, over-reporting, S3 on, S4 off; C23) against C05 in the 12 configurations; the first
# six rows are the reference case and the biological configurations
config_rows <- c(reference = "Reference case", adult_M_0.24 = "Adult M = 0.24 per year", adult_M_0.36 = "Adult M = 0.36 per year",
                 adult_M_0.54 = "Adult M = 0.54 per year", steepness_0.65 = "Steepness h = 0.65", steepness_0.95 = "Steepness h = 0.95",
                 connectivity_0.5 = "Connectivity × 0.5", connectivity_2 = "Connectivity × 2", selectivity_2.5 = "Selectivity a50 = 2.5",
                 selectivity_4.5 = "Selectivity a50 = 4.5", errors = "Errors in history\nand projection",
                 movement_information = "Movement information\nin estimation")
combined <- paired |> filter(S1 == 1, S2 == 1, S3 == 1, S4 == 0, metric %in% c("mean_D", "relative_catch")) |>
  mutate(row = match(config_id, names(config_rows)), mp = MP[match(mp, MP_LONG)], estimate = 100 * estimate, mcse = 100 * mcse)
pS3 <- list()
for (j in 1:2) {
  df <- filter(combined, metric == c("mean_D", "relative_catch")[j])
  df$y <- df$row + (as.integer(factor(df$mp, MP)) - 2.5) * 0.17
  pS3[[j]] <- ggplot(df) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = seq(1, 12, 2) - 0.46, ymax = seq(1, 12, 2) + 0.46, fill = BAND) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = 6.5, linewidth = 0.4) +    # 生物学配置与其他配置之间 | between the biological and the other configurations
    geom_segment(aes(x = 0, xend = estimate, y = y, yend = y, colour = mp), linewidth = 0.75, lineend = "butt") +
    geom_errorbar(aes(xmin = estimate - 1.96 * mcse, xmax = estimate + 1.96 * mcse, y = y, colour = mp), orientation = "y", width = 0.08, linewidth = 0.4) +
    geom_point(aes(estimate, y, colour = mp, shape = mp), size = 1.7) +
    scale_colour_manual(values = MP_COL) + scale_shape_manual(values = MP_SHAPE) +
    scale_y_reverse(breaks = 1:12, labels = if (j == 1) config_rows else NULL, limits = c(12.5, 0.5), expand = c(0, 0)) +
    labs(x = c("Change (percentage points)", "Change (% of final historical catch)")[j], y = NULL, tag = c("a)", "b)")[j]) +
    geom_blank(data = data.frame(lab = "x"), inherit.aes = FALSE) + facet_wrap(~lab, labeller = as_labeller(c(x = c(SB, "Catch")[j]))) +
    theme_paper + theme(axis.ticks.y = element_blank())
}
figS3 <- (pS3[[1]] | pS3[[2]]) / legends$MP + plot_layout(heights = c(1, 0.035))
ggsave(file.path(FIG, "FigS3_configurations.png"), figS3, width = 7.2, height = 6.4, dpi = 450, device = ragg::agg_png, bg = "white")
ggsave(file.path(FIG, "FigS3_configurations.pdf"), figS3, width = 7.2, height = 6.4, device = cairo_pdf, bg = "white")

# ================================================================ 补充图 4 | Supplementary Figure 4
# S3 开时，估计方法得到实际移动减去使用平均移动（对 S1、S2、S4 在重复内平均）
# With S3 on, the estimation method given the realized movement minus mean movement (averaged over S1, S2 and S4 within replicate)
info <- read_csv(file.path(OUT, "08_analysis", "information.csv"), show_col_types = FALSE) |> filter(S3 == 1) |>
  mutate(mp = MP[match(mp, MP_LONG)], row = match(mp, MP), estimate = 100 * estimate, mcse = 100 * mcse)
pS4 <- list()
for (j in 1:2) {
  df <- filter(info, metric == c("mean_D", "estimation_rmse_fixed")[j])
  df$y <- df$row                                     # 每行一个 MP | one MP per row
  pS4[[j]] <- ggplot(df) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = seq(1, 4, 2) - 0.46, ymax = seq(1, 4, 2) + 0.46, fill = BAND) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = estimate, y = y, yend = y, colour = mp), linewidth = 0.75, lineend = "butt") +
    geom_errorbar(aes(xmin = estimate - 1.96 * mcse, xmax = estimate + 1.96 * mcse, y = y, colour = mp), orientation = "y", width = 0.08, linewidth = 0.4) +
    geom_point(aes(estimate, y, colour = mp, shape = mp), size = 1.7) +
    scale_colour_manual(values = MP_COL) + scale_shape_manual(values = MP_SHAPE) +
    scale_y_reverse(breaks = 1:4, labels = MP, limits = c(4.5, 0.5), expand = c(0, 0)) +
    labs(x = "Change (percentage points)", y = NULL, tag = c("a)", "b)")[j]) +
    geom_blank(data = data.frame(lab = "x"), inherit.aes = FALSE) +
    facet_wrap(~lab, labeller = as_labeller(c(x = c(SB, paste0("Error of estimated ", SB))[j]))) +
    theme_paper + theme(axis.ticks.y = element_blank())
}
figS4 <- (pS4[[1]] | pS4[[2]]) / legends$MP + plot_layout(heights = c(1, 0.09))
ggsave(file.path(FIG, "FigS4_movement_information.png"), figS4, width = 7.2, height = 2.9, dpi = 450, device = ragg::agg_png, bg = "white")
ggsave(file.path(FIG, "FigS4_movement_information.pdf"), figS4, width = 7.2, height = 2.9, device = cairo_pdf, bg = "white")
cat("13 done:", FIG, "\n")
