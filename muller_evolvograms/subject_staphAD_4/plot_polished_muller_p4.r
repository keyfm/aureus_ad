# Rscript to generate polished muller plot for subject 4
library(ggmuller)
library(ggplot2)

edges <- read.table('results_nested/tables/muller_extended_tbl.edges.tsv',sep="\t",stringsAsFactors=F,header=TRUE)
pop = read.table('results_nested/tables/muller_extended_tbl.populations.tsv',sep="\t",header=TRUE,stringsAsFactors=F)
muller_df <- get_Muller_df(edges, pop)
# Muller_plot(muller_df, add_legend = TRUE)

pdf('ggmuller_plot_pat4.pdf', width=14, height=4.5, onefile = TRUE) 
my_palette <- c('white','#fcae91','#fb6a4a','#de2d26','#a50f15') #fee5d9
mc=my_palette#[1:5]
names(mc)=c('genotype-0','genotype-2','genotype-3','genotype-1','genotype-4')
ggplot(muller_df, aes_string(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) + 
    theme_classic() +
    geom_area() +
    theme(legend.position = "none",text = element_text(size=30,color="black"),axis.line = element_line(color = "black"),axis.text.x = element_text(color = "black"),axis.text.y = element_text(color = "black"),axis.ticks=element_line(color = "black")) + 
    guides(linetype = FALSE, color = FALSE) + 
    scale_y_continuous(breaks = c(0,0.5,1),labels = 50 * (0:2), name = "Frequency") +
    scale_x_continuous(breaks = c(0,100,200,1200),labels=c("0\n3","1\n14","2\n60","12\n50"), name="Time (in months)",limits = c(0, 1400)) + 
    scale_fill_manual(name = "Identity", values = mc) +
    scale_color_manual(values = mc)   
dev.off()

