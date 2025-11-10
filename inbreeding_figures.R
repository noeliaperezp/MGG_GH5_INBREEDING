# LIBRARIES
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("ggpmisc", quietly = TRUE)) install.packages("ggpmisc")
library(ggplot2)
library(ggpmisc)


################################################################################
# INBREEDING DEPRESSION (ID)
################################################################################

# Data
phenotypes <- read.table("pop_ID.phe")
inbreeding <- read.table("pop_ID_pruned.inbreeding", header=TRUE)
TABLE <- data.frame(inbreeding,"phenotype"=phenotypes$V1)

# Fhat3
Fhat3 <- ggplot() + 
  geom_point(data=TABLE, aes(x=Fhat3, y=phenotype), alpha=0.8, color="black") +
  
  stat_poly_eq(data=TABLE,aes(x=Fhat3, y=phenotype,
               label = paste(..eq.label.., sep = "~~~")), label.x = "center", label.y = 0.98,
               eq.with.lhs = "italic(y)~`=`~",  eq.x.rhs = "~italic(x)", parse = TRUE, size = 4.6, color="black") +

  geom_smooth(data=TABLE,
              aes(x=Fhat3, y=phenotype),method = "lm", se = F, color="black", size=0.5) +
  
  labs(x=expression(bolditalic(F["hat3"])), y=expression(bold(Log)~(bolditalic(P)))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold",vjust = 0.5, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        axis.text = element_text(size=12, face = "bold"),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.3, 0.1), "cm"))

ggsave("pop_ID_pruned.Fhat3.jpg", plot = Fhat3, width = 4.4, height = 3.5, dpi=600)

# FROH01
FROH01 <- ggplot() + 
  geom_point(data=TABLE, aes(x=FROH01, y=phenotype), alpha=0.8, color="black") +
  
  stat_poly_eq(data=TABLE,aes(x=FROH01, y=phenotype,
                              label = paste(..eq.label.., sep = "~~~")), label.x = "center", label.y = 0.98,
               eq.with.lhs = "italic(y)~`=`~",  eq.x.rhs = "~italic(x)", parse = TRUE, size = 4.6, color="black") +
  
  geom_smooth(data=TABLE,
              aes(x=FROH01, y=phenotype),method = "lm", se = F, color="black", size=0.5) +
  
  labs(x=expression(bolditalic(F["ROH01"])), y=expression(bold(Log)~(bolditalic(P)))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold",vjust = 0.5, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        axis.text = element_text(size=12, face = "bold"),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.3, 0.1), "cm"))

ggsave("pop_ID_pruned.FROH01.jpg", plot = FROH01, width = 4.4, height = 3.5, dpi=600)

# FROH1
FROH1 <- ggplot() + 
  geom_point(data=TABLE, aes(x=FROH1, y=phenotype), alpha=0.8, color="black") +
  
  stat_poly_eq(data=TABLE,aes(x=FROH1, y=phenotype,
                              label = paste(..eq.label.., sep = "~~~")), label.x = "center", label.y = 0.98,
               eq.with.lhs = "italic(y)~`=`~",  eq.x.rhs = "~italic(x)", parse = TRUE, size = 4.6, color="black") +
  
  geom_smooth(data=TABLE,
              aes(x=FROH1, y=phenotype),method = "lm", se = F, color="black", size=0.5) +
  
  labs(x=expression(bolditalic(F["ROH1"])), y=expression(bold(Log)~(bolditalic(P)))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold",vjust = 0.5, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        axis.text = element_text(size=12, face = "bold"),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.3, 0.1), "cm"))

ggsave("pop_ID_pruned.FROH1.jpg", plot = FROH1, width = 4.4, height = 3.5, dpi=600)

# FROH5
FROH5 <- ggplot() + 
  geom_point(data=TABLE, aes(x=FROH5, y=phenotype), alpha=0.8, color="black") +
  
  stat_poly_eq(data=TABLE,aes(x=FROH5, y=phenotype,
                              label = paste(..eq.label.., sep = "~~~")), label.x = "center", label.y = 0.98,
               eq.with.lhs = "italic(y)~`=`~",  eq.x.rhs = "~italic(x)", parse = TRUE, size = 4.6, color="black") +
  
  geom_smooth(data=TABLE,
              aes(x=FROH5, y=phenotype),method = "lm", se = F, color="black", size=0.5) +
  
  labs(x=expression(bolditalic(F["ROH5"])), y=expression(bold(Log)~(bolditalic(P)))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14,face="bold"),
        axis.title.y = element_text(size = 14,face="bold",vjust = 0.5, margin = margin(t = 0, r = 6, b = 0, l = 0)),
        axis.text = element_text(size=12, face = "bold"),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.3, 0.1, 0.3, 0.1), "cm"))

ggsave("pop_ID_pruned.FROH5.jpg", plot = FROH5, width = 4.4, height = 3.5, dpi=600)


