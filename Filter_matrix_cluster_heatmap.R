library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dplyr)
library(ggplot2)
library(reshape2)




#Expression data table 
data <- read.table("OC_site_all_expressions.csv",header = TRUE,row.names=1,sep=",")
#list of genes to filter data by
genelist <- read.csv('C:/Users/Mason/PycharmProjects/gene_stat5332/selected_features.csv')





#function to flip datatable 
tdt <- function(dt) {
  names<- dt[,1]
  dt.T <- as.data.frame(as.matrix(t(dt[,-1])))
  colnames(dt.T) <-names
  return(dt.T)
}
#function to count gene matrix clusters to determine if viable for vizualizing clear 
# it averages the genes inside the cluster. Using its Pearson correlation number. 
# .3 is a number i have selected. I can not tell you why it is important i do not remember
# but it helps understand how well divereged the cluster is
avg_count<- function(x,number){
  clust_by_row<- filter(x, cluster==number) 
  cluster_list<- rownames(clust_by_row)
  col.num <- which(colnames(clust_by_row) %in% cluster_list)
  cluster_matrix <- clust_by_row[,sort(c(col.num))]
  Cluster_gene_avg <- colMeans(cluster_matrix)
  names(Cluster_gene_avg)<- NULL
  cluster_count <- length(Cluster_gene_avg[Cluster_gene_avg>.3])
  return(cluster_count)
  
}

#filters data data by genelist
filtered_data <- merge(genelist,data,by.x=1,by.y=0)
filtered_data <- filtered_data%>% select_if(~ !any(is.na(.)))
rownames(filtered_data) <- filtered_data[,1]
filtered_data <-tdt(filtered_data)

#removes sample names 
df<-Filter(function(x) sd(x) != 0, filtered_data)

# creates matrix table and adds gene name to a column
cormat <-cor(df,method="pearson")
df.cormat <-as.data.frame(cormat)
df <- tibble::rownames_to_column(df.cormat, "Gene")
df<- na.omit(df)
# Dissimilarity matrix
d <- dist(df, method = "euclidean")

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "ward.D" )
hc1$labels
hc1$labels <-df$Gene
gene.order <- hc1$order
#labels cluster order as gene order
gene.order<- hc1$labels[hc1$order]
gene.order.column<- append(gene.order,'cluster')

hc1.order<- hc1$order

# Cut tree into # of clusters
sub_grp <- cutree(hc1, k = 6)

# Number of members in each cluster
table(sub_grp)

# attaches cluster numbers to the genes. 
df.clusters <- df %>% mutate(cluster = sub_grp) 
df.ordered <- df.clusters %>% slice(match(gene.order.column,Gene))
df.matrix <-df.ordered[c(gene.order)]
df.matrix2 <- df.ordered[c(gene.order.column)]
rownames(df.matrix2) <- c(gene.order)
df.matrix3 <- as.matrix(df.matrix)
write.csv(df.matrix2,"ordered_gene_matrix.csv")

# plots denogram plus coloring 
plot(hc1, cex = 0.6, hang = -1)
rect.hclust(hc1, k = 6, border = 2:5)

# makes the matrix usable for ggplot
melted_df.matrix <- melt(df.matrix3)

ggplot(data = melted_df.matrix, aes(x=Var1, y=Var2, fill=value))+ 
  geom_tile() +
  scale_fill_gradient2(low="blue", high="red") +
  scale_x_discrete("Genes") +
  scale_y_discrete("")

cluster_1_avg <-avg_count(df.matrix2,1)
print(cluster_1_avg)
cluster_2_avg <-avg_count(df.matrix2,2)
print(cluster_2_avg)
cluster_3_avg <-avg_count(df.matrix2,3)
print(cluster_3_avg)
cluster_4_avg <-avg_count(df.matrix2,4)
print(cluster_4_avg)
cluster_5_avg <-avg_count(df.matrix2,5)
print(cluster_5_avg)
cluster_6_avg <-avg_count(df.matrix2,6)
print(cluster_6_avg)

cluster_gene_list <- function(x,number){
  clust_by_row<- filter(x, cluster==number) 
  cluster_list<- rownames(clust_by_row)
  return(c(cluster_list))
}
cluster_6_list <- cluster_gene_list(df.matrix2,6)
cluster_5_list <- cluster_gene_list(df.matrix2,5)
cluster_4_list <- cluster_gene_list(df.matrix2,4)
cluster_3_list <- cluster_gene_list(df.matrix2,3)
cluster_2_list <- cluster_gene_list(df.matrix2,2)
cluster_1_list <- cluster_gene_list(df.matrix2,1)

# writes the gene lists from clusters 
write.csv(cluster_6_list,'Foldchange_cluster_6_gene_list.csv')
write.csv(cluster_5_list,'Foldchange_cluster_5_gene_list.csv')
write.csv(cluster_4_list,'Foldchange_cluster_4_gene_list.csv')
write.csv(cluster_3_list,'Foldchange_cluster_3_gene_list.csv')
write.csv(cluster_2_list,'Foldchange_cluster_2_gene_list.csv')
write.csv(cluster_1_list,'Foldchange_cluster_1_gene_list.csv')




















