


#Use the function below to make a graphic showing group-level correlations based on the top N most variable genes in the dataset

make_corr_heatmap = function(object, slot = 'counts', ngenes = 2000, display_nums = FALSE, assay = "RNA", method = "pearson", textsize = 10, group_by = 'seurat_clusters') {
  Idents(object) = object[[group_by]]
  vargenes = FindVariableFeatures(object, nfeatures = ngenes, assay = assay)@assays$RNA@var.features
  av.exp <- AverageExpression(object, slot = slot)$RNA
  av.exp <- av.exp[vargenes,]
  cor.exp <- as.data.frame(cor(av.exp, method = method))
  cor.exp$x <- rownames(cor.exp)
  cor.df <- tidyr::gather(data = cor.exp, y, correlation, names(table(Idents(object))))
  return(pheatmap(cor.exp[1:length(cor.exp) - 1], display_numbers = display_nums, fontsize = textsize))
}


#Use the below function to either get enrichment results from escape, or add them directly to the seurat object in a new assay.
#This one can take awhile if you don't have many cores to use. 

Do_escape = function(object, category, subcategory = FALSE, group_by = "seurat_clusters", sample_n_cells = Inf, cores = 2, groups = 1000, species = "Mus musculus", add_to_object = FALSE) {
  require(GSEABase)
  require(msigdb)
  require(ExperimentHub)
  require(Seurat)
  require(escape)
  require(SingleCellExperiment)
  require(msigdbr)
  gc()
  #https://bioconductor.org/packages/release/bioc/vignettes/escape/inst/doc/vignette.html
  Idents(object) = object[[group_by]]
  object = subset(object, downsample = sample_n_cells)
  print(object)
  
  if (subcategory != FALSE) {
  gene_set = msigdbr(species = species, category = category, subcategory = subcategory)
  } else {
    gene_set = msigdbr(species = species, category = category)
  }
  gene_set = gene_set %>% split(x = .$gene_symbol, f = .$gs_name)
  print(head(gene_set))
  ES.seurat <- enrichIt(obj = object, 
                      gene.sets = gene_set, 
                      groups = groups, cores = cores)
  if (add_to_object) {
    object = Seurat::AddMetaData(object, ES.seurat)
    object[["escape_enrichment"]] = CreateAssayObject(data = t(as.matrix(ES.seurat)))
    return(object)
  }else{
    return(ES.seurat)
  }
}


#The below function changes the level order of an object's identitiy class
reorder_idents = function(object, group_by, new_levels) {
  require(Seurat)
  object[[group_by]][[1]] <- factor(x = object[[group_by]][[1]], levels = new_levels)
  return(object)
}

#The below function allows one to quickly get the genes involved with a set of gene sets. Output is a list, so can function as a custom gene set for enrichment analysis. 
make_custom_gset_list = function(organism_type = 'mm', id_type = 'SYM', version = '7.4', brownse_gene_sets = FALSE, bahareh_gene_sets) {
  #http://127.0.0.1:10378/library/msigdb/doc/msigdb.html
  #https://www.bioconductor.org/packages/release/bioc/vignettes/GSEABase/inst/doc/GSEABase.pdf
  require(GSEABase)
  require(msigdb)
  require(ExperimentHub)
  require(GSEA)
  require(gage)
  require(GSVA)
  library(fgsea)
  eh = ExperimentHub()
  if (brownse_gene_sets) {
    print(query(eh, 'msigdb'))
    return()
    }
  gene_set = msigdb::getMsigdb(org = organism_type, id = id_type, version = version)
  gene_set = appendKEGG(gene_set)
    
  gset_list = list()
  for (b_gset in bahareh_gene_sets) {
    specific_gene_set = gene_set[[b_gset]]
    print(specific_gene_set)
    gset_list[[b_gset]] = specific_gene_set@geneIds
  }
  return(gset_list)
}


run_dorothea_and_viper_on_seurat = function(object, cores = 1) {
  require(dorothea)
  require(ggplot2)
  require(tidyverse)
  data(dorothea_mm, package = "dorothea")
  # run viper
  tf_activities <- run_viper(object, dorothea_mm,
  options = list(method = "scale", minsize = 4,
  eset.filter = FALSE, cores = cores,
  verbose = FALSE))
  
  return(tf_activities)
}



make_cell_type_bar_plot = function(object, xident, yident, make_prop = TRUE, convert_small_yidents_to_other = FALSE, small_ident_thresh = 0.05, fontsize = 5, xname = deparse(substitute(xident)), fillname = deparse(substitute(yident)), yname = ifelse(test = make_prop, yes = "proportion of total cells", no = "cell counts")) {
  require(tidyverse)
  print(xname)
  
  meta_data = object@meta.data
  meta_data_2 <- meta_data  %>% dplyr::select({{xident}}, {{yident}}) %>% group_by({{xident}}) %>% table()
  meta_data_2 = data.frame(t(data.frame(rbind(meta_data_2))))
  
  if (convert_small_yidents_to_other) {
    count_totals = data.frame(apply(X = meta_data_2, MARGIN = 1, FUN = sum))
    count_props = count_totals/sum(count_totals)
    
    names_of_lowbies = rownames(subset(count_props, apply.X...meta_data_2..MARGIN...1..FUN...sum. < small_ident_thresh))
    names_of_highbies = rownames(subset(count_props, apply.X...meta_data_2..MARGIN...1..FUN...sum. >= small_ident_thresh))
    meta_data_lowbie = meta_data_2[names_of_lowbies,]
    other_vals = apply(meta_data_lowbie, MARGIN = 2, FUN = sum)
    meta_data_2 = meta_data_2[names_of_highbies,]
    print("other vals")
    print(other_vals)
    meta_data_2 = rbind(meta_data_2, other = other_vals)
    print(meta_data_2)
    
  }
  
  if (make_prop) {
    find_proportion = function(input_set) {
      return(input_set/sum(input_set))
    }
    meta_data_2 = data.frame(apply(meta_data_2, MARGIN = 2, FUN = find_proportion))
  }
  print(meta_data_2)
  #meta_data_2$xident = rownames(meta_data_2)
  columns_to_elongate = colnames(meta_data_2)
  meta_data_2$yident = rownames(meta_data_2)
  meta_data_2 = meta_data_2 %>% pivot_longer(cols = columns_to_elongate, names_to = 'xident', values_to = 'count_or_prop')
  
  #bar_plot = ggplot(meta_data_2, aes(x = rownames(meta_data_2), y = ))
  
  plot_object = ggplot(meta_data_2, aes(x = xident, y = count_or_prop, fill = yident)) + geom_bar(stat = 'identity')+ theme(text=element_text(size=fontsize)) + labs(x = xname, y = yname, fill = fillname)
  print(plot_object)
  
}



make_volcano = function(object, de, log_fc_cutoff = 0.5, p_val_cutoff = NULL, groupby = 'seurat_clusters', graph_title = 'add a title idiot', overlap_metric = 17, point_font_size = ) {
    require(Seurat)
    require(ggplot2)
    require(ggrepel)
    Idents(object) = object[[groupby]]
    
    de$reg = ""
    de$reg[de$p_val < p_val_cutoff & abs(de$avg_log2FC) > log_fc_cutoff & de$avg_log2FC > 0] <- "UP"
    de$reg[de$p_val < p_val_cutoff & abs(de$avg_log2FC) > log_fc_cutoff & de$avg_log2FC < 0] <- "DOWN"
    de$name = rownames(de)
    de$name[de$reg == ""] <- ""
    
   
    plot = ggplot(data=de, aes(x=avg_log2FC, y=-log10(p_val_adj), col=reg, label=name)) + 
        geom_point(color = 'black', size = 2.5) + 
        theme_minimal() +
        ggrepel::geom_text_repel(color = 'Black', max.overlaps = overlap_metric, ) +
        scale_color_manual(breaks = c("DOWN", "", "UP"),values=c("Blue", "Gray", "Yellow")) + geom_point(size = 2) + 
        ggtitle(graph_title)
    print(plot)
    

}

generate_3d_plot = function(object, num_PCs = 20, downsample_size = Inf, group.by = 'seurat_clusters') {
  
  object = subset(object, downsample = downsample_size)

  combined.sub <- RunUMAP(object, dims = 1:num_PCs, n.components = 3L)
  PC_embeddings = Embeddings(object = combined.sub, reduction = 'umap')
  plot.data <- FetchData(object = combined.sub, vars = c("UMAP_1", "UMAP_2", "UMAP_3", group.by))
  plot.data$label <- paste(rownames(plot.data))
  print(head(plot.data))
  
  figfull = plot_ly(x=plot.data$UMAP_1, y=plot.data$UMAP_2, z=plot.data$UMAP_3, type="scatter3d", mode="markers", color = plot.data$tissue, size = 15)
  figfull = figfull %>% layout(scene = list(xaxis = list(title = "umap 1"), yaxis = list(title = "umap 2"), zaxis = list(title = "umap 3")))
  
  
  
  

  plot.data.av = plot.data %>% group_by_({{group.by}}) %>% summarise(average_UMAP1 = mean(UMAP_1), average_UMAP2 = mean(UMAP_2), average_UMAP3 = mean(UMAP_3), grouping = group.by)
  fig = plot_ly(x=plot.data.av$average_UMAP1, y=plot.data.av$average_UMAP2, z=plot.data.av$average_UMAP3, type="scatter3d", mode="markers", color = plot.data.av$grouping, size = 15)
  
  
  fig = fig %>% layout(scene = list(xaxis = list(title = "umap 1"), yaxis = list(title = "umap 2"), zaxis = list(title = "umap 3")))
  
  
  
  
}




















