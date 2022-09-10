


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









