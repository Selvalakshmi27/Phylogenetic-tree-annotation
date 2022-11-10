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
  ggplot2
)



meta_data1 <- read.csv("Kleborate_results.csv")
meta_data1

df1 <- as.data.frame(meta_data1)
rownames(df1) <- df1$strain
df1[1] <- NULL
df1$ST
df1
tree<- read.tree("core_gene_alignment.aln.treefile")
tree <- midpoint.root(tree)


attach(df1)


gg_tree<-ggtree(tree, color="black",layout="circular", branch.length = 20)  + 
  geom_hilight(node=60, fill="gold") + geom_hilight(node=79, fill="#c8b6ff")  + geom_hilight(node=74, fill="hot pink")+
  geom_hilight(node=48, fill="#90e0ef") + geom_hilight(node=56, fill="#80ed99",show.legend =TRUE)+ 
  geom_cladelab(60,"CG307", barsize=1, barcolor="white", 
                hjust=.5,angle=0,offset=0.00003,offset.text=0.0004,
                hjust=0.5,fontcolor='red', fontsize=5, fontface="bold")+
  geom_cladelab(79,"CG20", barsize=1, barcolor="white", 
                hjust=.5,angle=0,offset=0.0001,offset.text=0.0006,
                hjust=0.5,fontcolor='red', fontsize=5, fontface="bold") +
  geom_cladelab(74,"CG29", barsize=1, barcolor="white", 
                hjust=.5,angle=0,offset=0.0002,offset.text=0.0005,
                hjust=0.5,fontcolor='red', fontsize=5, fontface="bold") +
  geom_cladelab(48,"CG15", barsize=1, barcolor="white", 
                hjust=.5,angle=0,offset=0.0001,offset.text=0.00027,
                hjust=0.5,fontcolor='red', fontsize=5, fontface="bold") +
  geom_cladelab(56,"CG258", barsize=1, barcolor="white", 
                hjust=.5,angle=0,offset=0.0001,offset.text=0.0006,
                hjust=0.5,fontcolor='red', fontsize=5, fontface="bold") +
  geom_treescale(fontsize = 2, linesize = 0.01, family = "Arial Bold") 
 

gg_tree %<+% df1 + 
  geom_tippoint(aes(color=Clonal.group), show.legend =FALSE) + theme(legend.position = "none")

ggsave("p1.png", dpi=600)






