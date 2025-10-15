library(ape)
library(phylolm)
library(dplyr)
library(caper)
library(gghalves)
library(phytools)
library(ggpubr)
calculate_pgls_residuals <- function(tree, df, label = NULL) {
  if (is.null(label)) {
    label <- deparse(substitute(df))
    label <- sub(".*_", "", label)  # 提取最后一个下划线后的部分
  }
  # 构建comparative data
  stats <- comparative.data(tree, df, Species, vcv = TRUE, vcv.dim = 3)
  # 进行PGLS拟合
  model <- pgls(log10(beak_width) ~ log10(mass), data = stats, lambda = "ML")
  # 提取残差
  residuals_mod <- residuals(model)
  residuals_mod <- as.data.frame(residuals_mod)
  residuals_mod$Species <- rownames(residuals_mod)
  residuals_mod <- residuals_mod[,c(2,1)]
  colnames(residuals_mod) <- c("Species", "residual")
  # 合并原数据和残差
  merged_data <- merge(df, residuals_mod, by = "Species", all = TRUE)
  # 添加分组标签
  # merged_data$Group <- label
  # return(merged_data)
  return(list(data = merged_data, model = model))
}

tree <- read.tree("./1034_newick.txt")
tree$node.label <- NULL
data <- read.csv("./1034_morph_group.csv") 
data <- data %>% dplyr::select(species_name,group,mass,beak_length,beak_width,beak_depth)
colnames(data) <- c("Species","pgls.group","mass","beak_length","beak_width","beak_depth")

shared_species <- intersect(data$Species, tree$tip.label)
filtered_data <- data[data$Species %in% shared_species, ]
pruned_tree <- ape::keep.tip(tree, shared_species)


df <- filtered_data %>% dplyr::select(Species,pgls.group,mass,beak_width)
df <- na.omit(df)

# 创建空列表存结果
result_list <- list()

# 开始循环
for (i in 1:13) {
  # 筛选数据
  cat(i)
  df_i <- df %>% filter(pgls.group == i)
  
  # 修剪树
  tree_i <- keep.tip(pruned_tree, df_i$Species)
  
  # 按树的tip顺序重新排列df
  df_i <- df_i[match(tree_i$tip.label, df_i$Species), ]
  
  # 计算pgls residuals
  result_i <- calculate_pgls_residuals(tree_i, df_i)
  
  # 保存到列表
  result_list[[i]] <- result_i
}
# 合并rbs
combined_df <- do.call(rbind, lapply(seq_along(result_list), function(i) {
  data_i <- result_list[[i]]$data
  data_i
}))

# 整合获取斜率和截距
regression_lines <- do.call(rbind, lapply(seq_along(result_list), function(i) {
  coef_i <- summary(result_list[[i]]$model)$coefficients
  data.frame(
    pgls.group = i,
    intercept = coef_i[1],
    slope = coef_i[2],
    intercept_SE = coef_i[3],
    slope_SE = coef_i[4]
  )
}))




# 先计算每组的 log10(Mass) 范围
range_df <- combined_df %>%
  group_by(pgls.group) %>%
  summarise(
    x_min = min(log10(mass), na.rm = TRUE),
    x_max = max(log10(mass), na.rm = TRUE)
  )

# 把范围和回归参数合并
regression_lines <- left_join(regression_lines, range_df, by = "pgls.group")

# 根据 y = intercept + slope * x 算出 y 的起止点
regression_lines <- regression_lines %>%
  mutate(
    y_min = intercept + slope * x_min,
    y_max = intercept + slope * x_max
  )

combined_df_copy <- combined_df

combined_df$plotgroup <- ifelse(combined_df$pgls.group %in% c(1,10,2,5), "Low",
                                ifelse(combined_df$pgls.group %in% c(8,11), "Intermediate_low",
                                       ifelse(combined_df$pgls.group %in% c(4,13,12), "Intermediate_high", 
                                              ifelse(combined_df$pgls.group %in% c(3,6,9,7), "High", NA))))
#备份一下

regression_lines_copy <- regression_lines

# 数据分组
group_points_1 <- combined_df %>% filter(pgls.group == 1)
group_points_2 <- combined_df %>% filter(pgls.group == 2)
group_points_3 <- combined_df %>% filter(pgls.group == 3)
group_points_4 <- combined_df %>% filter(pgls.group == 4)
group_points_5 <- combined_df %>% filter(pgls.group == 5)
group_points_6 <- combined_df %>% filter(pgls.group == 6)
group_points_7 <- combined_df %>% filter(pgls.group == 7)
group_points_8 <- combined_df %>% filter(pgls.group == 8)
group_points_9 <- combined_df %>% filter(pgls.group == 9)
group_points_10 <- combined_df %>% filter(pgls.group == 10)
group_points_11 <- combined_df %>% filter(pgls.group == 11)
group_points_12 <- combined_df %>% filter(pgls.group == 12)
group_points_13 <- combined_df %>% filter(pgls.group == 13)

group_line_1 <- regression_lines %>% filter(pgls.group == 1)
group_line_2 <- regression_lines %>% filter(pgls.group == 2)
group_line_3 <- regression_lines %>% filter(pgls.group == 3)
group_line_4 <- regression_lines %>% filter(pgls.group == 4)
group_line_5 <- regression_lines %>% filter(pgls.group == 5)
group_line_6 <- regression_lines %>% filter(pgls.group == 6)
group_line_7 <- regression_lines %>% filter(pgls.group == 7)
group_line_8 <- regression_lines %>% filter(pgls.group == 8)
group_line_9 <- regression_lines %>% filter(pgls.group == 9)
group_line_10 <- regression_lines %>% filter(pgls.group == 10)
group_line_11 <- regression_lines %>% filter(pgls.group == 11)
group_line_12 <- regression_lines %>% filter(pgls.group == 12)
group_line_13 <- regression_lines %>% filter(pgls.group == 13)

#### redo plot #####
custom_colors <- c("1" = "grey80",
                   "2" = "grey80",
                   "3" = "grey80",
                   "4" = "grey80",
                   "5" = "grey80",
                   "6" = "grey80",
                   "7" = "grey80",
                   "8" = "grey80",
                   "9" = "grey80",
                   "10" = "#3B4E76",
                   "11" = "#B3CFDD",
                   "12" = "#EBD22A",
                   "13" = "#A0D568")
# 
p1 <- ggplot() +
  theme_minimal(base_size = 7) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    legend.title = element_blank()
  ) +
  labs(
    x = "log10(mass)",
    y = "log10(beak_width)"
  )+
  # point 
  # 1
  geom_point(data = group_points_1,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["1"], size = 0.5) +
  # 2
  geom_point(data = group_points_2,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["2"], size = 0.5) +
  # 3
  geom_point(data = group_points_3,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["3"], size = 0.5) +
  # 4
  geom_point(data = group_points_4,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["4"], size = 0.5) +
  # 5
  geom_point(data = group_points_5,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["5"], size = 0.5) +
  # 6
  geom_point(data = group_points_6,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["6"], size = 0.5) +
  # 7
  geom_point(data = group_points_7,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["7"], size = 0.5) +
  # 8
  geom_point(data = group_points_8,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["8"], size = 0.5) +
  # 9
  geom_point(data = group_points_9,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["9"], size = 0.5) +
  # line
  # 1
  geom_segment(data = group_line_1,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["1"], size = 0.5)+
  # 2
  geom_segment(data = group_line_2,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["2"], size = 0.5) +
  # 3
  geom_segment(data = group_line_3,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["3"], size = 0.5)+
  # 4
  geom_segment(data = group_line_4,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["4"], size = 0.5)+
  # 5
  geom_segment(data = group_line_5,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["5"], size = 0.5)+
  # 6
  geom_segment(data = group_line_6,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["6"], size = 0.5) +
  # 7
  geom_segment(data = group_line_7,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["7"], size = 0.5)+
  # 8
  geom_segment(data = group_line_8,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["8"], size = 0.5)+
  # 9
  geom_segment(data = group_line_9,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["9"], size = 0.5)+
  # 12
  geom_point(data = group_points_12,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["12"], size = 0.5) +
  # 11
  geom_point(data = group_points_11,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["11"], size = 0.5) +
  # 13
  geom_point(data = group_points_13,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["13"], size = 0.5) +
  # 10
  geom_point(data = group_points_10,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["10"], size = 0.5) +
  # 12
  geom_segment(data = group_line_12,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["12"], size = 0.5)+
  # 11
  geom_segment(data = group_line_11,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["11"], size = 0.5)+
  # 13
  geom_segment(data = group_line_13,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["13"], size = 0.5)+
  # 10
  geom_segment(data = group_line_10,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["10"], size = 0.5)
  
p1


custom_colors <- c("1" = "grey80",
                   "2" = "#7A01E6",
                   "3" = "#DE9525",
                   "4" = "grey80",
                   "5" = "grey80",
                   "6" = "grey80",
                   "7" = "grey80",
                   "8" = "#7097E4",
                   "9" = "grey80",
                   "10" = "grey80",
                   "11" = "grey80",
                   "12" = "grey80",
                   "13" = "grey80")
# 
p2 <- ggplot() +
  theme_minimal(base_size = 7) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    legend.title = element_blank()
  ) +
  labs(
    x = "log10(mass)",
    y = "log10(beak_width)"
  )+
  # point 
  # 1
  geom_point(data = group_points_1,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["1"], size = 0.5) +
  
  # 4
  geom_point(data = group_points_4,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["4"], size = 0.5) +
  # 5
  geom_point(data = group_points_5,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["5"], size = 0.5) +
  # 6
  geom_point(data = group_points_6,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["6"], size = 0.5) +
  # 7
  geom_point(data = group_points_7,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["7"], size = 0.5) +
  # 9
  geom_point(data = group_points_9,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["9"], size = 0.5) +
  # 10
  geom_point(data = group_points_10,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["10"], size = 0.5) +
  # 11
  geom_point(data = group_points_11,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["11"], size = 0.5) +
  # 12
  geom_point(data = group_points_12,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["12"], size = 0.5) +
  # 13
  geom_point(data = group_points_13,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["13"], size = 0.5) +
  # line
  # 1
  geom_segment(data = group_line_1,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["1"], size = 0.5)+
  
  # 4
  geom_segment(data = group_line_4,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["4"], size = 0.5)+
  # 5
  geom_segment(data = group_line_5,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["5"], size = 0.5)+
  # 6
  geom_segment(data = group_line_6,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["6"], size = 0.5) +
  # 7
  geom_segment(data = group_line_7,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["7"], size = 0.5)+
  
  # 9
  geom_segment(data = group_line_9,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["9"], size = 0.5)+
  # 10
  geom_segment(data = group_line_10,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["10"], size = 0.5) +
  # 11
  geom_segment(data = group_line_11,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["11"], size = 0.5)+
  # 12
  geom_segment(data = group_line_12,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["12"], size = 0.5)+
  # 13
  geom_segment(data = group_line_13,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["13"], size = 0.5)+
  # 2
  geom_point(data = group_points_2,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["2"], size = 0.5) +
  # 3
  geom_point(data = group_points_3,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["3"], size = 0.5) +
  # 8
  geom_point(data = group_points_8,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["8"], size = 0.5) +
  # 2
  geom_segment(data = group_line_2,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["2"], size = 0.5) +
  # 3
  geom_segment(data = group_line_3,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["3"], size = 0.5)+
  # 8
  geom_segment(data = group_line_8,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["8"], size = 0.5)
p2


custom_colors <- c("1" = "grey80",
                   "2" = "grey80",
                   "3" = "grey80",
                   "4" = "#437920",
                   "5" = "#7B8CDE",
                   "6" = "grey80",
                   "7" = "grey80",
                   "8" = "grey80",
                   "9" = "#CF3D3E",
                   "10" = "grey80",
                   "11" = "grey80",
                   "12" = "grey80",
                   "13" = "grey80")
# 
p3 <- ggplot() +
  theme_minimal(base_size = 7) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    legend.title = element_blank()
  ) +
  labs(
    x = "log10(mass)",
    y = "log10(beak_width)"
  )+
  # point 
  # 1
  geom_point(data = group_points_1,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["1"], size = 0.5) +
  # 2
  geom_point(data = group_points_2,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["2"], size = 0.5) +
  # 3
  geom_point(data = group_points_3,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["3"], size = 0.5) +
  # 6
  geom_point(data = group_points_6,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["6"], size = 0.5) +
  # 7
  geom_point(data = group_points_7,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["7"], size = 0.5) +
  # 8
  geom_point(data = group_points_8,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["8"], size = 0.5) +
  
  # 10
  geom_point(data = group_points_10,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["10"], size = 0.5) +
  # 11
  geom_point(data = group_points_11,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["11"], size = 0.5) +
  # 12
  geom_point(data = group_points_12,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["12"], size = 0.5) +
  # 13
  geom_point(data = group_points_13,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["13"], size = 0.5) +
  # line
  # 1
  geom_segment(data = group_line_1,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["1"], size = 0.5)+
  # 2
  geom_segment(data = group_line_2,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["2"], size = 0.5) +
  # 3
  geom_segment(data = group_line_3,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["3"], size = 0.5)+
  
  
  # 6
  geom_segment(data = group_line_6,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["6"], size = 0.5) +
  # 7
  geom_segment(data = group_line_7,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["7"], size = 0.5)+
  # 8
  geom_segment(data = group_line_8,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["8"], size = 0.5)+
  
  # 10
  geom_segment(data = group_line_10,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["10"], size = 0.5) +
  # 11
  geom_segment(data = group_line_11,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["11"], size = 0.5)+
  # 12
  geom_segment(data = group_line_12,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["12"], size = 0.5)+
  # 13
  geom_segment(data = group_line_13,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["13"], size = 0.5)+
  # 5
  geom_point(data = group_points_5,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["5"], size = 0.5) +
  # 4
  geom_point(data = group_points_4,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["4"], size = 0.5) +
  # 9
  geom_point(data = group_points_9,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["9"], size = 0.5) +
  # 5
  geom_segment(data = group_line_5,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["5"], size = 0.5)+
  # 4
  geom_segment(data = group_line_4,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["4"], size = 0.5)+
  # 9
  geom_segment(data = group_line_9,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["9"], size = 0.5)
  
p3


custom_colors <- c("1" = "#5A5A5A",
                   "2" = "grey80",
                   "3" = "grey80",
                   "4" = "grey80",
                   "5" = "grey80",
                   "6" = "#C05454",
                   "7" = "#911417",
                   "8" = "grey80",
                   "9" = "grey80",
                   "10" = "grey80",
                   "11" = "grey80",
                   "12" = "grey80",
                   "13" = "grey80")
# 
p4 <- ggplot() +
  theme_minimal(base_size = 7) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    legend.title = element_blank()
  ) +
  labs(
    x = "log10(mass)",
    y = "log10(beak_width)"
  )+
  # point 
  
  # 2
  geom_point(data = group_points_2,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["2"], size = 0.5) +
  # 3
  geom_point(data = group_points_3,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["3"], size = 0.5) +
  # 4
  geom_point(data = group_points_4,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["4"], size = 0.5) +
  # 5
  geom_point(data = group_points_5,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["5"], size = 0.5) +
  # 8
  geom_point(data = group_points_8,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["8"], size = 0.5) +
  # 9
  geom_point(data = group_points_9,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["9"], size = 0.5) +
  # 10
  geom_point(data = group_points_10,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["10"], size = 0.5) +
  # 11
  geom_point(data = group_points_11,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["11"], size = 0.5) +
  # 12
  geom_point(data = group_points_12,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["12"], size = 0.5) +
  # 13
  geom_point(data = group_points_13,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["13"], size = 0.5) +
  # line
  
  # 2
  geom_segment(data = group_line_2,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["2"], size = 0.5) +
  # 3
  geom_segment(data = group_line_3,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["3"], size = 0.5)+
  # 4
  geom_segment(data = group_line_4,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["4"], size = 0.5)+
  # 5
  geom_segment(data = group_line_5,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["5"], size = 0.5)+
  
  # 8
  geom_segment(data = group_line_8,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["8"], size = 0.5)+
  # 9
  geom_segment(data = group_line_9,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["9"], size = 0.5)+
  # 10
  geom_segment(data = group_line_10,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["10"], size = 0.5) +
  # 11
  geom_segment(data = group_line_11,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["11"], size = 0.5)+
  # 12
  geom_segment(data = group_line_12,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["12"], size = 0.5)+
  # 13
  geom_segment(data = group_line_13,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["13"], size = 0.5)+
  # 6
  geom_point(data = group_points_6,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["6"], size = 0.5) +
  # 1
  geom_point(data = group_points_1,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["1"], size = 0.5) +
  # 7
  geom_point(data = group_points_7,
             aes(x = log10(mass), y = log10(beak_width)),
             color = custom_colors["7"], size = 0.5) +
  # 6
  geom_segment(data = group_line_6,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["6"], size = 0.5) +
  # 1
  geom_segment(data = group_line_1,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["1"], size = 0.5)+
  # 7
  geom_segment(data = group_line_7,
               aes(x = x_min - 0.1, xend = x_max + 0.1, 
                   y = y_min - 0.1, yend = y_max + 0.1),
               color = custom_colors["7"], size = 0.5)
p4
p5 <- p1+p2+p3+p4
