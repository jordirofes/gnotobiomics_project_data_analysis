group_distribution_plots <- function(methrix_obj, group_var, scale_dt = TRUE, 
                                     plot_title = "group values distribution"){
    
    if(scale_dt){
        cpg_dt <- t(scale(t(assay(methrix_obj))))    
      } else{
        cpg_dt <- assay(methrix_obj)
      }
      plot_dt <- lapply(unique(methrix_obj[[group_var]]), function(ctrs_id){
        data.frame("CpG_Means" = rowMeans2(cpg_dt, cols = methrix_obj[[group_var]] == ctrs_id, na.rm = TRUE),
                    "Group" = ctrs_id)   
      })
      plot_dt <- do.call(rbind, plot_dt)
      p_1 <- ggplot(plot_dt, aes(x = CpG_Means, fill = Group)) + 
                        geom_histogram(position = "identity") + theme_minimal() +
                        ggtitle(paste("Mean", plot_title))
      
      
      plot_dt <- lapply(unique(methrix_sub[[params$group_var]]), function(ctrs_id){
        data.frame("CpG_Means" = matrixStats::rowMedians(cpg_dt, cols = methrix_sub[[group_var]] == ctrs_id), "Group" = ctrs_id)
      })
      plot_dt <- do.call(rbind, plot_dt)
      p_2 <- ggplot(plot_dt, aes(x = CpG_Means, fill = Group)) + geom_histogram(position = "identity") + theme_minimal()+
                        ggtitle(paste("Median", plot_title))
      
      return(list(p_1, p_2))
}

sample_distribution_plot <- function(methrix_obj, group_var, scale_dt = TRUE, 
                                     plot_title = "Sample values distribution (10000 random CpGs)"){
  if(scale_dt){
    cpg_dt <- t(scale(t(assay(methrix_obj))))    
  } else{
    cpg_dt <- assay(methrix_obj)
  }
  plot_dt2 <- lapply(colnames(methrix_sub), function(sample_id){
    data.frame("Sample_Values" = cpg_dt[,sample_id], 
               "Sample_ID" = sample_id, "Sample_Group" = methrix_sub[,sample_id][[group_var]])
  })
  plot_dt2 <- do.call(rbind, plot_dt2)
  ggplot(plot_dt2, aes(x = Sample_Values, fill = Sample_ID)) + 
             geom_histogram(position = "identity", bins = 100) +
             theme_minimal() + facet_grid(rows = vars(Sample_Group))+
             ggtitle(plot_title)
}


### Annotation Plot

dmr_annotation_plot <- function(omic_dt, 
                                annotation_labels = c("Promoter", "Exon", 
                                                      "5' UTR", "3' UTR", 
                                                      "Distal Intergenic", 
                                                      "Intron", "Downstream"), 
                                annotation_metadata_var = "annotation", 
                                x_lab = "Region Proportions"){
  if(length(annotation_labels) <= 8){
    color_palette <- scale_fill_brewer(palette = "Set2")
  } else{
    color_palette <- NULL
  }
  
  factor_tables <- lapply(omic_dt, function(exp_dt){
    if(is(exp_dt, "SummarizedExperiment")){
      if(annotation_metadata_var == "seqnames"){
        annotation_variable <- as.vector(seqnames(rowRanges(exp_dt)))
      } else{
        annotation_variable <- rowData(exp_dt)[[annotation_metadata_var]]    
      }
    } else if(is(exp_dt, "GRanges")){
      if(annotation_metadata_var == "seqnames"){
        annotation_variable <- as.vector(seqnames(exp_dt))
      } else{
        annotation_variable <- mcols(exp_dt)[[annotation_metadata_var]]
      }
    }
    exp_annotation <- sapply(as.character(annotation_labels), function(x){
      found_anno <- grepl(pattern = x, annotation_variable)
      sum(found_anno)
    })
    exp_annotation
    # table(unlist(exp_annotation))
  })
  factor_tables <- do.call(rbind, factor_tables)
  factor_annotation_dt <-factor_tables/rowSums(factor_tables)
  
  plot_dt <- reshape2::melt(factor_annotation_dt)
  
  plot_dt$Var1 <- factor(plot_dt$Var1, levels = rev(levels(plot_dt$Var1)))
  ggplot(plot_dt, aes(x = value, y = Var1, fill = Var2)) + geom_col() + color_palette + 
    theme_minimal() + xlab(x_lab) + ylab(NULL) + theme(legend.title = element_blank())
}