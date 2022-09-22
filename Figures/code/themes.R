theme_Publication <- function(base_size=14) {
library(grid)
library(ggthemes)
(theme_foundation(base_size=base_size)
+ theme(plot.title = element_text(face = "bold",
size = rel(1.2), hjust = 0.5),
text = element_text(),
panel.border = element_rect(colour = NA),
axis.title = element_text(size = rel(1)),
axis.title.y = element_text(angle=90,vjust =2),
axis.title.x = element_text(vjust = -0.2),
axis.text = element_text(),
axis.line = element_line(colour="black"),
axis.ticks = element_line(),
panel.grid.minor = element_blank(),
panel.grid.major = element_blank(),
panel.background = element_rect(fill = "transparent",colour = NA),
plot.background = element_rect(fill = "transparent",colour = NA),
legend.key = element_rect(colour = NA),
legend.position = "right",
legend.direction = "vertical",
legend.key.size= unit(0.2, "cm"),
legend.title = element_text(face="italic"),
strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
strip.text = element_text(face="bold")
))
}
scale_colour_Publication <- function(...){
load('data/color_abb_codes.Robj')
scale_color_manual(values=color_code$color)
}
scale_fill_Publication <- function(...){
load('data/color_abb_codes.Robj')
scale_fill_manual(values=color_code$color)
}
