miR_target_filter<-miR_target[miR_target$genesymbol %in% colnames(target_heatmap),] 
connect <- miR_target_filter[c(1,3,4)]
connect$mirnaid<-gsub("hsa-","",connect$mirnaid)
names(connect)[3]<-"value"
c( as.character(connect$mirnaid), as.character(connect$genesymbol)) %>%
  as.tibble() %>%
  group_by(value) %>%
  summarize(n=n()) -> coauth
colnames(coauth) <- c("name", "n")
dim(coauth)
mygraph <- graph_from_data_frame(connect, vertices = coauth, directed = FALSE )
com <- walktrap.community(mygraph)
max(com$membership)
coauth<-rbind(coauth[7:27,],coauth[c(1:6),],coauth[c(28:38),])
mygraph <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )

mycolor <- colormap(colormap=colormaps$jet, nshades=3)
mycolor <- sample(mycolor, length(mycolor))
coauth$color<-c(rep(mycolor[2],21),rep(mycolor[3],17))

ggraph(mygraph, layout="linear") + 
  geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
  geom_node_point(aes(size=n,color = coauth$color, fill = coauth$color), alpha=0.5) +
  scale_size_continuous(range=c(0.5,8)) +
  scale_color_manual(values=mycolor) +
  geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0.4,0), "null"),
    panel.spacing=unit(c(0,0,3.4,0), "null")
  ) +
  expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2))

miR_heapmap<-miR_heapmap[,colnames(miR_heapmap) %in% coauth$name]
pheatmap(miR_heapmap,cluster_rows = F, cluster_cols = F,treeheight_col = 0,
         show_colnames =T,show_rownames = T, 
         color = colorRampPalette(color.key)(50),
         cutree_cols = 0,
         #scale = "row",
         na_col = "white"
         #annotation_colors = color.annotation,
         #annotation_col = annotation_col,
         #cellwidth = 20,
         #cellheight = 10
         #file = "/pheatmap.pdf",plot = "pdf"
)

intersect(miRsummary)