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

meta_data1 <- read.csv("Kleborate_main.csv")


rownames(df1) <- df1$X
df1 <- as.data.frame(meta_data1)
df1[1] <- NULL
attach(df1)
df1
tree<- read.tree("core_gene_alignment_main.aln.treefile")
tree <- midpoint.root(tree)



gg_tree<-ggtree(tree)  + 
  geom_hilight(node=60, fill="gold") + geom_hilight(node=79, fill="#c8b6ff")  + geom_hilight(node=74, fill="hot pink")+
  geom_hilight(node=48, fill="#90e0ef") + geom_hilight(node=56, fill="#80ed99",show.legend =TRUE) + geom_treescale(fontsize = 2, linesize = 0.01, family = "Arial Bold") 
gg_tree
p1<- gg_tree + geom_tippoint(color="black",size = 1)


p1
table(df1$ESBL_gene)

df1$ESBL_gene[df1$ESBL_gene == ""] <- NA
df1$ESBL_gene
p3 <-gheatmap(p1, df1[, "ESBL_gene", drop=FALSE], offset=0.0001, width=0.017, hjust = 1, colnames_angle = -45,
              colnames_position = "top", color = "white", font.size = 2.5, family="Arial Bold", custom_column_labels = c("ESBL")) +
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5) +
  ggplot2::ylim(-1.0, NA) + 
  coord_cartesian(clip = "off") + 
  scale_fill_manual(values = c("#F781BF", "#377EB8", "#4DAF4A"),
                    name = "ESBL",
                    breaks =c("CTX-M-15", "CTX-M-27", "SHV-31"),na.value="white")
p3
p4<-p3+new_scale_fill()

df1$Chromosomal.Integration[df1$Chromosomal.Integration == "No"] <- NA
table(df1$Plasmid)


p5 <-gheatmap(p4, df1[, "Plasmid", drop=FALSE], offset=0.00018, width=0.017, hjust = 1, colnames_angle = -45,
              colnames_position = "top", color = "white", font.size = 2.5, family="Arial Bold", custom_column_labels = c("plasmid")) +
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5) +
  ggplot2::ylim(-1.0, NA) + 
  coord_cartesian(clip = "off") + 
  scale_fill_manual(values = c("#66C2A5", "#FC8D62" ,"#8DA0CB"),
                    name = "Plasmid",
                    breaks =c("FIA", "FIB(K)", "FIB(pQil)"), na.value="white")
p5
p6<-p5+new_scale_fill()

p7<-p6+new_scale_fill()
table(df1$GyrA.83I)
df1$GyrA.87N[df1$GyrA.87N== ""] <- NA

p8<- gheatmap(p7, df1[, "GyrA.83I", drop=FALSE], offset=0.0004, width=0.017, hjust = 1, colnames_angle = -45,
              colnames_position = "top", color = "white", font.size = 2.8, family="Arial Bold", custom_column_labels = c("GyrA-83")) +
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5) +
  ggplot2::ylim(-1.0, NA) + 
  coord_cartesian(clip = "off") + 
  scale_fill_manual(values = c("#f13c77","#377EB8","#4DAF4A","#984EA3","white"),
                    name = "GyrA-83",
                    breaks =c("GyrA-83F", "GyrA-83I", "GyrA-83Y", "GyrA-87G","-"), na.value="white")

p8

p9<-p8+new_scale_fill()
table(df1$GyrA.87N)

p10<- gheatmap(p9, df1[, "GyrA.87N", drop=FALSE], offset=0.00048, width=0.017, hjust = 1, colnames_angle = -45,
              colnames_position = "top", color = "white", font.size = 2.8, family="Arial Bold", custom_column_labels = c("GyrA-87")) +
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5) +
  ggplot2::ylim(-1.0, NA) + 
  coord_cartesian(clip = "off") + 
  scale_fill_manual(values = c("#f13c77","#377EB8","#4DAF4A"),
                    name = "GyrA-87",
                    breaks =c("GyrA-87A", "GyrA-87G", "GyrA-87N"), na.value="white")

p10
table(df1$qnr)

df1$qnr[df1$qnr== "qnrA1"] <- 'qnrA'
p11<-p10+new_scale_fill()
p12<- gheatmap(p11, df1[, "qnr", drop=FALSE], offset=0.00056, width=0.017, hjust = 1, colnames_angle = -45,
               colnames_position = "top", color = "white", font.size = 3, family="Arial Bold", custom_column_labels = c("qnr")) +
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5) +
  ggplot2::ylim(-1.0, NA) + 
  coord_cartesian(clip = "off") + 
  scale_fill_manual(values = c("#f13c77","#377EB8","#66C2A5","red"),
                    name = "qnr",
                    breaks =c("qepA2", "qnrA", "qnrB", "qnrS1"), na.value="white")


p12
p13<- p12+new_scale_fill()

df1$P.like.fimbriae[df1$P.like.fimbriae == "?-type"] <- "π-type"
df1$X.1.fimbriae[df1$X.1.fimbriae=="?1-fimbriae"] <- "γ1-fimbriae" 

df1$Bla_Carb_acquired
df1$Bla_Carb_acquired[df1$Bla_Carb_acquired== "OXA-181"] <- NA


p14<- gheatmap(p13, df1[, "X.1.fimbriae", drop=FALSE], offset=0.00074, width=0.017, hjust = 1, colnames_angle = -45,
               colnames_position = "top", color = "white", font.size = 2.8, family="Arial Bold", custom_column_labels = c("γ1")) +
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5) +
  ggplot2::ylim(-1.0, NA) + 
  coord_cartesian(clip = "off") + 
  scale_fill_manual(values = c("gold"),
                    name = "Type-1-fimbriae",
                    breaks =c("γ1-fimbriae"), na.value="white")

p14








p14<- gheatmap(p13, df1[, "P.like.fimbriae", drop=FALSE], offset=0.00064, width=0.017, hjust = 1, colnames_angle = -45,
               colnames_position = "top", color = "white", font.size = 2.8, family="Arial Bold", custom_column_labels = c("π-type")) +
  theme(legend.text=element_text(size=8, face="bold"),
        legend.title=element_text(size=8, face="bold"),
        legend.key.size = unit(0.5, 'cm'),
        legend.title.align=0.5) +
  ggplot2::ylim(-1.0, NA) + 
  coord_cartesian(clip = "off") + 
  scale_fill_manual(values = c("red"),
                    name = "π-type Fimbriae",
                    breaks =c("π-type"), na.value="white")


p14
p5
ggsave("Fimbriae-Pi.png",width = 15.36, height = 8.14, dpi = 900)
