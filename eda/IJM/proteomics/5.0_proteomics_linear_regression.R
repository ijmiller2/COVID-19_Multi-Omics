
library(cowplot)
library(ggplot2)
library(ggpmisc)

# load proteomics data with clinical metadata
proteomics_data_path = "~/Desktop/UW_2020/COVID19_multi-omics/COVID-19_Multi-Omics/data/proteomics_measurements_w_clinical_metadata.csv"

proteomics_df <- read.table(proteomics_data_path, sep=",", header=T)
covid_df <- subset(proteomics_df, COVID == 1)
no_covid_df <- subset(proteomics_df, COVID == 0)

# combine all datasets
my.formula <- y ~ x
a <- ggplot(proteomics_df, aes(x=Hospital_free_days_45,y=X7974)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula= my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +  
  labs(title = "Hospital Free Days", 
       y = "Gelsolin log2(int.)",
       x = "") +
  theme(plot.title = element_text(color = "black", size = 18, angle = 0, hjust = .5, vjust = .5, face = "bold"),
        axis.text.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = .5, face = "bold"),
        axis.text.y = element_text(color = "grey40", size = 20, angle = 0, hjust = 1, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 22, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "grey20", size = 22, angle = 90, hjust = .5, vjust = 1.5, face = "bold"),
        axis.ticks.x = element_line(size=1),
        axis.ticks.y = element_line(size=1),
        axis.ticks.length=unit(.25, "cm"),
        axis.line = element_line(colour = "grey20", size = 1),
        panel.grid.major = element_blank(), 
        panel.grid.major.y = element_line(color = "grey80", size=1),
        panel.background = element_blank()) +
  facet_grid(rows = vars(COVID))
a

# combine all datasets
my.formula <- y ~ x
b <- ggplot(proteomics_df, aes(x=Vent_free_days,y=X7974)) +
  geom_point(size=3) +
  geom_smooth(method='lm', formula= my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +  
  labs(title = "Vent Free Days", 
       y = "Gelsolin log2(int.)",
       x = "") +
  theme(plot.title = element_text(color = "black", size = 18, angle = 0, hjust = .5, vjust = .5, face = "bold"),
        axis.text.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = .5, face = "bold"),
        axis.text.y = element_text(color = "grey40", size = 20, angle = 0, hjust = 1, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 22, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "grey20", size = 22, angle = 90, hjust = .5, vjust = 1.5, face = "bold"),
        axis.ticks.x = element_line(size=1),
        axis.ticks.y = element_line(size=1),
        axis.ticks.length=unit(.25, "cm"),
        axis.line = element_line(colour = "grey20", size = 1),
        panel.grid.major = element_blank(), 
        panel.grid.major.y = element_line(color = "grey80", size=1),
        panel.background = element_blank()) +
  facet_grid(rows = vars(COVID))
b

### Box plots with ANOVA ###
c <- ggplot(proteomics_df, aes(x=color_by, y=X7974, fill=color_by)) + 
  #geom_boxplot(aes(fill=color_by), show.legend = F) +
  geom_jitter(aes(color=color_by), position=position_dodge(0.5), show.legend = F) +
  geom_boxplot(outlier.shape=NA) + 
  #geom_point(position=position_dodge(width=2)) +
  labs(title = "", 
       y = "Gelsolin log2(int.)",
       x = "Group") +
  theme(plot.title = element_text(color = "black", size = 18, angle = 0, hjust = .5, vjust = .5, face = "bold"),
        axis.text.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = .5, face = "bold"),
        axis.text.y = element_text(color = "grey40", size = 20, angle = 0, hjust = 1, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 22, angle = 0, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "grey20", size = 22, angle = 90, hjust = .5, vjust = 1.5, face = "bold"),
        axis.ticks.x = element_line(size=1),
        axis.ticks.y = element_line(size=1),
        axis.ticks.length=unit(.25, "cm"),
        axis.line = element_line(colour = "grey20", size = 1),
        panel.grid.major = element_blank(), 
        panel.grid.major.y = element_line(color = "grey80", size=1),
        panel.background = element_blank()) 
c


linearMod <- lm(X7974 ~ Hospital_free_days_45, data=covid_df)  # build linear regression model on full data
print(linearMod)
summary(linearMod)

linearMod <- lm(X7974 ~ Hospital_free_days_45, data=no_covid_df)  # build linear regression model on full data
print(linearMod)
summary(linearMod)

### with confounders / covariates ###
linearMod <- lm(X7974 ~ Hospital_free_days_45 + Age_less_than_90 + Gender, data=covid_df)
print(linearMod)
summary(linearMod)