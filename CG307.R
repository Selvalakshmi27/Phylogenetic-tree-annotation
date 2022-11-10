pacman::p_load(
  tidyverse,
  ape,
  ggtree,
  rhierbaps,
  phytools,
  ggnewscale,
  ggtreeExtra,
  ggpubr,
  scales, 
  ggplot2,
  RColorBrewer
)

meta_data1 <- read.csv("CG307.csv")


rownames(df1) <- df1$X
df1 <- as.data.frame(meta_data1)
df1[1] <- NULL
attach(df1)





gg_tree
tree<- read.tree("core_gene_alignment.aln.treefile")
tree <- midpoint.root(tree)




gg_tree<-ggtree(tree) + geom_treescale(fontsize = 2, linesize = 0.01, family = "Arial Bold")

p <- gg_tree %<+% df1 + geom_tippoint(aes(color=Location), size=3)+ scale_color_manual(values=c("blue", "red")) + 
  geom_tiplab(size=3, aes(color =Location)) + geom_nodelab() 
p

table(df1$ESBL_gene)

p2 <-gheatmap(p, df1[, "ESBL_gene", drop=FALSE], offset=0.000003, width=0.02, hjust = 1, colnames_angle = -45,
              colnames_position = "top", color = "white", font.size = 2.5, family="Arial Bold", custom_column_labels = c("ESBL")) +
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5) +
  ggplot2::ylim(-1.0, NA) + 
  coord_cartesian(clip = "off") + 
  scale_fill_manual(values = c("#F781BF"),
                    name = "ESBL",
                    breaks =c("CTX-M-15"),na.value="white")
p2
p3<-p2+new_scale_fill()
p4 <-gheatmap(p3, df1[, "Plasmid", drop=FALSE], offset=0.000004, width=0.02, hjust = 1, colnames_angle = -45,
              colnames_position = "top", color = "white", font.size = 2.5, family="Arial Bold", custom_column_labels = c("plasmid")) +
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5) +
  ggplot2::ylim(-1.0, NA) + 
  coord_cartesian(clip = "off") + 
  scale_fill_manual(values = c("#66C2A5", "#FC8D62" ,"#8DA0CB","#F781BF"),
                    name = "Plasmid",
                    breaks =c("FIA", "FIB(K)", "FIB(pQil)", "FIB(K),FIB(pQil),FIA"), na.value="white")
p4
p5<-p4+new_scale_fill()
table(df1$Chromosomal.Integration)
df1$Chromosomal.Integration[df1$Chromosomal.Integration == "No"] <- NA
p6 <-gheatmap(p5, df1[, "Chromosomal.Integration", drop=FALSE], offset=0.000005, width=0.02, hjust = 1, colnames_angle = -45,
              colnames_position = "top", color = "white", font.size = 2.5, family="Arial Bold", custom_column_labels = c("Chr")) +
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5) +
  ggplot2::ylim(-1.0, NA) + 
  coord_cartesian(clip = "off") + 
  scale_fill_manual(values = c("red"),
                    name = "Chromosome",
                    breaks =c("Yes"), na.value="white")
p6
p7<- p6+new_scale_fill()
p8<- gheatmap(p7, df1[, "GyrA.83I", drop=FALSE], offset=0.000008, width=0.02, hjust = 1, colnames_angle = -45,
         colnames_position = "top", color = "white", font.size = 2.5, family="Arial Bold", custom_column_labels = c("GyrA-83I")) +
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5) +
  ggplot2::ylim(-1.0, NA) + 
  coord_cartesian(clip = "off") + 
  scale_fill_manual(values = c("#377EB8"),
                    name = "GyrA-83I",
                    breaks =c("GyrA-83I"), na.value="white")
p9<- p8+new_scale_fill()
p10<- gheatmap(p9, df1[, "GyrA.87N", drop=FALSE], offset=0.000009, width=0.02, hjust = 1, colnames_angle = -45,
               colnames_position = "top", color = "white", font.size = 2.5, family="Arial Bold", custom_column_labels = c("GyrA-87N")) +
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5) +
  ggplot2::ylim(-1.0, NA) + 
  coord_cartesian(clip = "off") + 
  scale_fill_manual(values = c("#4DAF4A"),
                    name = "GyrA-87N",
                    breaks =c("GyrA-87N"), na.value="white")

p11<- p10+new_scale_fill()
df1$qnr[df1$qnr== "qnrB1"] <- 'qnrB'
p12<-  gheatmap(p11, df1[, "qnr", drop=FALSE], offset=0.000010, width=0.017, hjust = 1, colnames_angle = -45,
                colnames_position = "top", color = "white", font.size = 3, family="Arial Bold", custom_column_labels = c("qnr")) +
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5) +
  ggplot2::ylim(-1.0, NA) + 
  coord_cartesian(clip = "off") + 
  scale_fill_manual(values = c("#66C2A5"),
                    name = "qnr",
                    breaks =c( "qnrB"), na.value="white")
p12
df1$State
p13<- p12+new_scale_fill()
p14<-  gheatmap(p13, df1[, "State", drop=FALSE], offset=0.000020, width=0.017, hjust = 1, colnames_angle = -45,
                colnames_position = "top", color = "white", font.size = 3, family="Arial Bold", custom_column_labels = c("State")) +
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5) +
  ggplot2::ylim(-1.0, NA) + 
  coord_cartesian(clip = "off") + 
  scale_fill_manual(values = c("#66C2A5", "#F781BF", "#377EB8", "#4DAF4A"),
                    name = "State",
                    breaks =c( "LA","TX", "Chile","UAE"), na.value="white")

p14
table(df1$qnr)

ggsave("CG307_qnr1.png",width = 15.36, height = 8.14, dpi = 900)


