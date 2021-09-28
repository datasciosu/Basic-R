
## set working directory
setwd("~/Desktop")

## Packages 

#install.packages("package_name")
#library(ggplot2)

packages = c("ggplot2", "dplyr", "tidyr", "lemon", "ggpubr", "tibble", "scales")

for ( i in seq_along(packages) ){
  
  suppressPackageStartupMessages(library(packages[i],character.only = T))
  
}

## read in data
data = read.delim("example_data_928.txt", sep = "\t", header = TRUE)
head(data)

rownames(data) = data$stable_ID
data = subset(data, select = -c(stable_ID))

head(data)
dim(data)

## Compare gene expresssion between wild type and mutant

#1. barplot
head(data)

data.stacked = data %>% gather(key = "sample", value = "rpm")
head(data.stacked)

data.stacked.grouped = data.stacked %>% group_by(sample) %>% summarise(rpm = mean(rpm))
data.stacked.grouped

barplot(data.stacked.grouped$rpm ~ data.stacked.grouped$sample)

plot_bar = ggplot(data = data.stacked.grouped, aes(x = sample, y = rpm)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", fill = "grey60", lwd = 1) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        axis.ticks.length.y = unit(-.25, "cm"),
        axis.ticks.length.x = unit(-.25, "cm"),
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.y = element_text(margin=unit(c(0,0.5,0.5,0.5), "cm")),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        text = element_text(size=20),
        plot.title = element_text(size = 20, hjust = 0.5, vjust = -0.5),
        axis.text = element_text(size = 20, color = "black")) + 
  coord_cartesian() + 
  coord_capped_cart(bottom = "both", left = "both") + 
  labs(x = " ", y = "Reads / million")

plot_bar

#2. Boxplot 

# unlogged
plot_boxplot = ggplot(data = data.stacked, aes(x = sample, y = rpm)) + 
  geom_boxplot(width = 0.7, color = "black", fill = "grey60", lwd = 1) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        axis.ticks.length.y = unit(-.25, "cm"),
        axis.ticks.length.x = unit(-.25, "cm"),
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.y = element_text(margin=unit(c(0,0.5,0.5,0.5), "cm")),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        text = element_text(size=20),
        plot.title = element_text(size = 20, hjust = 0.5, vjust = -0.5),
        axis.text = element_text(size = 20, color = "black")) + 
  coord_cartesian() + 
  coord_capped_cart(bottom = "both", left = "both") + 
  labs(x = " ", y = "Reads / million")

plot_boxplot


# log2 transformation
plot_boxplot = ggplot(data = data.stacked, aes(x = sample, y = rpm)) + 
  geom_boxplot(width = 0.7, color = "black", fill = "grey70", lwd = 1) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        axis.ticks.length.y = unit(-.25, "cm"),
        axis.ticks.length.x = unit(-.25, "cm"),
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.y = element_text(margin=unit(c(0,0.5,0.5,0.5), "cm")),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        text = element_text(size=20),
        plot.title = element_text(size = 20, hjust = 0.5, vjust = -0.5),
        axis.text = element_text(size = 20, color = "black")) + 
  coord_cartesian() + 
  coord_capped_cart(bottom = "both", left = "both") + 
  labs(x = " ", y = "Reads / million") + 
  scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))

plot_boxplot


## Does gene expression go up significantly in the mutant
mutants = c("mut_1", "mut_2")
wildtype = c("wt_1", "wt_2")

head(data)
data['mut_avg'] = rowMeans(data[,mutants])
data['wt_avg'] = rowMeans(data[,wildtype])

# test for normality
# shapiro.test()

hist(data$wt_avg)
hist(log2(data$wt_avg))

hist(data$mut_avg)
hist(log2(data$mut_avg))

# compare median change in abundance using wilcoxon rank sum test
wilcox.test(data$wt_avg, data$mut_avg)

data.avg = subset(data, select = c(wt_avg, mut_avg))
data.avg.stacked = data.avg %>% gather(key = "sample", value = "rpm")

plot_boxplot = ggplot(data = data.avg.stacked, aes(x = sample, y = rpm)) + 
  geom_boxplot(width = 0.7, color = "black", fill = "grey70", lwd = 1) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        axis.ticks.length.y = unit(-.25, "cm"),
        axis.ticks.length.x = unit(-.25, "cm"),
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.y = element_text(margin=unit(c(0,0.5,0.5,0.5), "cm")),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        text = element_text(size=20),
        plot.title = element_text(size = 20, hjust = 0.5, vjust = -0.5),
        axis.text = element_text(size = 20, color = "black")) + 
  coord_cartesian() + 
  coord_capped_cart(bottom = "both", left = "both") + 
  labs(x = " ", y = "Reads / million") + 
  scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)))
plot_boxplot


## compare per gene the abundance in wild type compared to mutant 
head(data)

#1. we can look at the fold change (WT / MUT)
#2. look at if the change in abundance (WT vs MUT) is significant (p-value)

#1. Fold change calculation
data['fold_change'] = (data$mut_avg + 0.01) / (data$wt_avg + 0.01)
data['log2FoldChange'] = log2(data$fold_change)
head(data)

#2. P value calculation 
# t.test( [array1], [array2], alternative = "two.sided", var.equal = T) -> Students 2 sample t-test
# t.test( [array1], [array2], alternative = "two.sided", var.equal = F) -> Welch test

# apply()
data$gene = rownames(data)
data = data %>% 
        rowwise() %>% 
        mutate(pvalue = t.test( c(wt_1, wt_2), c(mut_1, mut_2), alternative = "two.sided", var.equal = T )$p.value) %>% 
        ungroup()
    
head(data)
view(data)

## scatter plot to look at the abundance of genes in wild type vs mutant (average of 2 biological replicates)

data = data %>%
  mutate(diff_express = ifelse( log2FoldChange >= .38 & pvalue < 0.05, "up",
                                ifelse( log2FoldChange <= -.38 & pvalue < 0.05, "down", "none")))


# look at genes go up 
up = subset(data, diff_express == "up")
down = subset(data, diff_express == "down")
none = subset(data, diff_express == "none")

write.table(data, "mut-vs-wt-analysis.csv", sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)

plot_scatter = ggplot() + 
  geom_point(data = data %>% filter(diff_express == "none"), aes(x = wt_avg, y = mut_avg), size = 0.6, color = "grey60") + 
  geom_point(data = data %>% filter(diff_express == "up"), aes(x = wt_avg, y = mut_avg), size = 0.8, color = "green3") +
  geom_point(data = data %>% filter(diff_express == "down"), aes(x = wt_avg, y = mut_avg), size = 0.6, color = "red") +
  theme_classic() + 
  theme(aspect.ratio = 1,
        axis.ticks.length.y = unit(-.25, "cm"),
        axis.ticks.length.x = unit(-.25, "cm"),
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.y = element_text(margin=unit(c(0,0.5,0.5,0.5), "cm")),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        text = element_text(size=20),
        plot.title = element_text(size = 20, hjust = 0.5, vjust = -0.5),
        axis.text = element_text(size = 20, color = "black")) + 
  coord_cartesian() + 
  coord_capped_cart(bottom = "both", left = "both") + 
  labs(x = "Wild Type Reads / million", y = "Mutant Reads / million") + 
  scale_y_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)),
                     limits = c(2^-10, 2^10)) + 
  scale_x_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x)),
                     limits = c(2^-10, 2^10)) + 
  geom_abline(slope = 1, intercept = 0, color = "grey80", lty = 2) + 
  geom_abline(slope = 1, intercept = 1, color = "grey80", lty = 2) + 
  geom_abline(slope = 1, intercept = -1, color = "grey80", lty = 2) 

plot_scatter  
ggsave(plot_scatter, filename = "pirna_diff_express_scatter.pdf", dpi = 300, height = 10, width = 10)  


plot_volcano = ggplot() + 
  geom_point(data = data %>% filter(diff_express == "none"), aes(x = log2FoldChange, y = -log10(pvalue)), size = 0.6, color = "grey60") + 
  geom_point(data = data %>% filter(diff_express == "up"),  aes(x = log2FoldChange, y = -log10(pvalue)), size = 0.8, color = "green3") +
  geom_point(data = data %>% filter(diff_express == "down"),  aes(x = log2FoldChange, y = -log10(pvalue)), size = 0.6, color = "red") +
  theme_classic() + 
  theme(aspect.ratio = 1,
        axis.ticks.length.y = unit(-.25, "cm"),
        axis.ticks.length.x = unit(-.25, "cm"),
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.y = element_text(margin=unit(c(0,0.5,0.5,0.5), "cm")),
        axis.line = element_line(size = 0.8),
        axis.ticks = element_line(size = 0.8),
        text = element_text(size=20),
        plot.title = element_text(size = 20, hjust = 0.5, vjust = -0.5),
        axis.text = element_text(size = 20, color = "black")) + 
  coord_cartesian() + 
  coord_capped_cart(bottom = "both", left = "both") + 
  labs(x = "log2(Fold Change)", y = "-log10(p-value)") + 
  ylim(0,6) + 
  geom_hline(yintercept = 1.3, color = "blue", lty = 2) + 
  geom_vline(xintercept = .38, color = "blue", lty = 2) + 
  geom_vline(xintercept = -.38, color = "blue", lty = 2) 

plot_volcano








