colors_pam50 <- gg_color_hue(6)
names(colors_pam50) <- c(
    "basal",
    "her2", 
    "lumb", 
    "luma", 
    "claudin-low", 
    "normal"
)
selected_colors <- intersect(names(colors_pam50), unique(selected_data$pam50))
selected_colors <- colors_pam50[selected_colors]