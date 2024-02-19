
plot_pc <- function(plot_dt, group_var, x_lab, y_lab, loading_dt = NA, 
                    allElipse = TRUE, ellipse_var = NA, interactive = TRUE, level = 0.95){
    pc_plot <- ggplot(plot_dt) + 
        geom_point(aes(x = .data[[colnames(plot_dt)[1]]], 
                    y = .data[[colnames(plot_dt)[2]]], colour = group_var)) + 
        theme_classic() + xlab(x_lab) + ylab(y_lab) +
        theme(axis.line = element_line(colour = "black", size = 0.5,
                                       linetype = "solid"),
              legend.title = element_blank(),
              plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(vjust = 0.5, hjust = 1)) + geom_hline(yintercept = 0) + 
        geom_vline(xintercept = 0)
    if(allElipse){
        pc_plot <- pc_plot + stat_ellipse(aes(x = .data[[colnames(plot_dt)[1]]], 
                                        y = .data[[colnames(plot_dt)[2]]]),
                                        colour = "grey7", type = "t", 
                                        linetype = 1,level = level)
    }
    
    if(!any(is.na(ellipse_var))){
        pc_plot <- pc_plot + stat_ellipse(aes(x = .data[[colnames(plot_dt)[1]]], 
                                        y = .data[[colnames(plot_dt)[2]]], 
                                        color = ellipse_var),
                                        type = "t", linetype = 1, level = level)
    }
    
    arrow_scale <- max(abs(plot_dt))
    
    if(!all(is.na(loading_dt))){
        pc_plot <- pc_plot + 
            annotate(geom = "segment", x = 0, y = 0, 
                            xend = loading_dt[,1]*(arrow_scale + arrow_scale/2), 
                            yend = loading_dt[,2]*(arrow_scale + arrow_scale/2), size = 0.2,  
                        arrow = arrow(length = unit(0.2, "cm"), ends = "last")) + 
            annotate(geom = "text",
                    x = loading_dt[,1]*(arrow_scale+ arrow_scale/2 + 1), 
                    y = loading_dt[,2]*(arrow_scale+ arrow_scale/2 + 1), 
                    label = rownames(loading_dt), 
                    size = 2.5, color = "blue")
    }
    
    if(interactive){
        pc_plot <- ggplotly(pc_plot)
    }
    return(pc_plot)
}

pca_plot <- function(pca_data, pc_x = 1, pc_y = 2, group_var, ellipse_var = NA, 
                     allElipse = TRUE, interactive = TRUE, level = 0.95){
    pc_var <- pca_data$sdev^2
    pc_var_abs <- pc_var/sum(pc_var)*100
    
    x_lab <- paste0("PC", pc_x, " (", round(pc_var_abs[pc_x]), "%)")
    y_lab <- paste0("PC", pc_y, " (",round(pc_var_abs[pc_y]), "%)")
    
    plot_dt <- as.data.frame(pca_data$x[,c(pc_x, pc_y)])
    
    loading_dt <- pca_data$rotation[,c(pc_x, pc_y)]
    
    loading_effect_sum <- apply(loading_dt, 1, function(x){sum(abs(x))})
    loading_names <- names(loading_effect_sum)[order(loading_effect_sum, 
                                                            decreasing = TRUE)]
    loading_dt <- loading_dt[loading_names[1:c(min(length(loading_names), 10))],]
    
    return(plot_pc(plot_dt, group_var, x_lab, y_lab, loading_dt, allElipse, 
                     ellipse_var, interactive, level))
}


var_plot <- function(data, nfilter = 500){
    
    vars_var <- rowVars(data)
    if(nfilter > vars_var){
        nfilter <- length(vars_var)
    }
    top_var_data <- data[order(vars_var, decreasing = TRUE)[1:nfilter],]
    pca_data <- prcomp(top_var_data)
    
    pc_var <- pca_data$sdev^2
    var_vec_1 <- cumsum(pc_var)/sum(pc_var)*100
    variance_plot <- qplot(y = var_vec_1) + theme_minimal() + 
        xlab("Number of components") + ylab("Cummulative variance (%)") + 
        ggtitle("Scree Plot") + theme(plot.title = element_text(hjust = 0.5))
    return(variance_plot)
}

quantile_breaks <- function(dt, length_out){
    dt_quantiles <- quantile(x = dt, seq(0, 1, length.out = length_out))
    return(dt_quantiles[!duplicated(dt_quantiles)])
}

heatMapFun <- function(dt_list, metadata, stratificationVar, cdQuantile = 0, 
                       plot_name = NULL, show_colnames, 
                       clustering_distance_rows, ...){
    group_var <- metadata[stratificationVar]
    dt_merged <- t(dt_list)
    if(any(is.na(dt_merged))){
        dt_merged <- t(impute::impute.knn(t(dt_merged))$data)
    }
    det_coef <- apply(dt_merged, 2, function(x){
        sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)
    })
    det_coef <- sort(det_coef, decreasing = TRUE)
    sel_vars <- names(det_coef)[det_coef >= quantile(det_coef, cdQuantile)]
    
    dt_quantiles <- quantile_breaks(dt = dt_merged[,sel_vars], 
                                    length_out = 101)
    
    if(cdQuantile != 0){
        plot_name <- paste0(plot_name, " Top ", 
                                     length(sel_vars),
                                     " variables (>=", cdQuantile*100, "%) CV")
    } else if(cdQuantile != 0 & is.null(plot_name)){
        plot_name <- paste("Top", length(sel_vars), "variables")
    }
    
    plot_colors <- viridis::viridis(n = length(dt_quantiles) - 1)
    cat("\n")
    pheatmap::pheatmap(dt_merged[,sel_vars], color = plot_colors,
                       annotation_row = as.data.frame(group_var), 
                       show_colnames = show_colnames, fontsize = 5, 
                       fontsize_col = 6, fontsize_row = 6, show_rownames = FALSE, 
                       clustering_method = "average", cluster_rows = TRUE, 
                       cluster_cols = TRUE, clustering_distance_rows = clustering_distance_rows, 
                       main = plot_name, 
                       breaks = dt_quantiles)
    cat("\n")
}


pcaFun <- function(dt_list, metadata, stratificationVar){
    patient_metadata_full_num <-  do.call(cbind, dt_list)
    if(any(is.na(patient_metadata_full_num))){
        pca_dt <- impute::impute.knn(t(patient_metadata_full_num))$data    
    } else{
        pca_dt <- t(patient_metadata_full_num)
    }
    
    pca_dt <- prcomp(t(pca_dt), scale. = TRUE)
    for(group_var_name in stratificationVar){
        group_var <-  metadata[[group_var_name]]
        cat("##### ", group_var_name, " \n")
        cat(knit_print(pca_plot(pca_data = pca_dt, pc_x = 1, pc_y = 2, 
                                group_var = group_var, 
                 ellipse_var = group_var)), " \n")
        cat("\n")
        pcvar <- round(pca_dt$sdev^2/sum(pca_dt$sdev^2)*100, 2)
    
        d3plot <- plot_ly(x = pca_dt$x[,1], y = pca_dt$x[,2], 
                          z = pca_dt$x[,3], color = group_var)
        d3plot <- d3plot %>% add_markers()
        d3plot <- d3plot %>% layout(scene = list(xaxis = list(title = paste0('PC1', "(", pcvar[1], "%)")),
                            yaxis = list(title = paste0('PC2', "(", pcvar[2], "%)")),
                            zaxis = list(title = paste0('PC3', "(", pcvar[3], "%)"))))
        cat(knit_print(d3plot), " \n")
        cat("\n")
    }
}

boxplotFun <- function(dt_list, metadata, stratificationVar){
    test_boxplots_list <- lapply(dt_list, function(dt){
        test_boxplots <- lapply(names(dt), function(vars_name){
            ggplot() + geom_boxplot(aes(x = metadata[[stratificationVar]], 
                                        y = dt[[vars_name]]), 
                                    fill = brewer.pal(n = 2, name = "Dark2")[1:2]) + 
                xlab("Recurrencia") + ylab(vars_name) + theme_minimal()
        })
        names(test_boxplots) <- names(dt)
        test_boxplots
    })

    for(dt_name in names(test_boxplots_list)){
        cat("##### ", dt_name, "{.tabset} \n")
        for(var_name in names(test_boxplots_list[[dt_name]])){
            cat("###### ", var_name, " \n")
            print(test_boxplots_list[[dt_name]][[var_name]])
            cat("\n")
            cat("\n")
        }
        cat("\n")
    }
}

significant_tables_chunks <- function(contrs_result_sig_list, contrs_names, template = c(
    "##### `r contrs_names[{{i}}]` \n",
    "```{r, echo = FALSE}\n",
    "datatable(contrs_result_sig_list[[{{i}}]], options = list(scrollX = TRUE)) %>% formatRound(columns = colnames(contrs_result_sig_list[[{{i}}]])[-7], digits = 3)",
    "```\n",
    "\n"
)){
    if(missing(contrs_names)){
        contrs_names <- names(contrs_result_sig_list)
    }
    tables_chunks <- lapply(seq_along(contrs_result_sig_list), function(i){
        knit_expand(text = template)
    })
    return(tables_chunks)
}


annotate_data <- function(deseqRes, annotationPackage){
    ensembl_id <- rownames(deseqRes)
    entrez_ids <- AnnotationDbi::mapIds(annotationPackage, keys = ensembl_id, column = c("ENTREZID"), 
                                             keytype = "ENSEMBL", multiVals = "first")
    symbols <- AnnotationDbi::mapIds(annotationPackage, keys = ensembl_id, column = c("SYMBOL"), 
                                             keytype = "ENSEMBL", multiVals = "first")
    gene_name <- AnnotationDbi::mapIds(annotationPackage, keys = ensembl_id, column = c("GENENAME"), 
                                            keytype = "ENSEMBL", multiVals = "first")
    annotated_data_frame <- cbind(deseqRes, entrez_ids, symbols, gene_name)
    return(annotated_data_frame)
}

enrich_filter <- function(dt_list){
    filtered_data <- lapply(dt_list, function(x){
        x$entrez_ids[!is.na(x$entrez_ids)]
    })
    return(filtered_data)
}

multi_enrich <- function(dt_to_pathway, universe, OrgDb = "org.Hs.eg.db", 
                        keyType = keyType, qvalueCutoff = 0.01, 
                        pAdjustMthd = "BH", pvalueCutoff = 0.01, 
                        simplify_res = FALSE, simplify_cutoff = 0.6, 
                        simplify_col = "p.adjust", simplify_fun = min,
                        enrichFun = c("ora", "gsea", "kegg", "wp", 
                                      "david", "msigdbrORA", "gseKegg", 
                                      "gseWp", "msigdbrGSEA"), 
                        keggOrg, keggKeyType, davidOrg, wpOrg, 
                        msigdbCategory, msigdbSpc, ...){
    
    if(!is(enrichFun, "function")){
        enrichFun <- match.arg(enrichFun)
        enrich_fun <- switch (enrichFun,
                              "ora" = enrichGO,
                              "gsea" = gseGO,
                              "kegg" = enrichKEGG,
                              "gseKegg" = gseKEGG,
                              "gseWp" = gseWP,
                              "wp" = enrichWP,
                              "david" = enrichDAVID,
                              "msigdbrORA" = enricher,
                              "msigdbrGSEA" = GSEA
        )
        switch(enrichFun, "ora" = {
            param_list <- list(keyType = keyType, OrgDb = OrgDb, 
                               pAdjustMethod = pAdjustMthd, 
                               pvalueCutoff = pvalueCutoff,
                               qvalueCutoff = qvalueCutoff)
            
            gene_input_name <- "gene"
            param_list$universe <- universe
            param_list$readable <- TRUE
        }, "gsea" = {
            param_list <- list(keyType = keyType, OrgDb = OrgDb, 
                               pAdjustMethod = pAdjustMthd, 
                               pvalueCutoff = pvalueCutoff)
            
            gene_input_name <- "geneList"
            dt_to_pathway <- lapply(dt_to_pathway, sort, decreasing = TRUE)
        }, "kegg" = {
            param_list <- list(keyType = keggKeyType,
                               pAdjustMethod = pAdjustMthd, 
                               pvalueCutoff = pvalueCutoff)
            
            gene_input_name <- "gene"
            dt_to_pathway <- lapply(dt_to_pathway, mapIds, x = OrgDb, 
                                    column = "ENTREZID", keytype = keyType, 
                                    multiVals = "first")
            param_list$organism <- keggOrg
            
            param_list$universe <- universe
            param_list$qvalueCutoff  <- qvalueCutoff
            simplify_res <- FALSE
        }, "wp" = {
            param_list <- list(pAdjustMethod = pAdjustMthd, 
                               pvalueCutoff = pvalueCutoff)
            
            wpOrg <- match.arg(wpOrg, get_wp_organisms())
            param_list$organism <- wpOrg
            dt_to_pathway <- lapply(dt_to_pathway, mapIds, x = OrgDb, 
                                    column = "ENTREZID", keytype = keyType, 
                                    multiVals = "first")
            gene_input_name <- "gene"
            param_list$qvalueCutoff  <- qvalueCutoff
            simplify_res <- FALSE
            
            
        }, "david" = {
            dt_to_pathway <- lapply(dt_to_pathway, mapIds, x = OrgDb, 
                                    column = "ENTREZID", keytype = keyType, 
                                    multiVals = FALSE)
            param_list$idType <- "ENTREZ_GENE_ID"
        }, "msigdbrORA"= {
            param_list <- list(pAdjustMethod = pAdjustMthd, 
                               pvalueCutoff = pvalueCutoff,
                               qvalueCutoff = qvalueCutoff)
            dt_to_pathway <- lapply(dt_to_pathway, mapIds, x = OrgDb, 
                                    column = "ENTREZID", keytype = keyType, 
                                    multiVals = "first")
            C3_t2g <- msigdbr(species = msigdbSpc, category = msigdbCategory) %>% 
                            dplyr::select(gs_name, entrez_gene)
            
            enrichFun <- paste0(enrichFun, "_", msigdbCategory)
            
            param_list$TERM2GENE <- C3_t2g
            gene_input_name <- "gene"
            simplify_res <- FALSE
            
        }, "msigdbrGSEA" = {
            param_list <- list(pAdjustMethod = pAdjustMthd, 
                               pvalueCutoff = pvalueCutoff)
            C3_t2g <- msigdbr(species = msigdbSpc, category = msigdbCategory) %>% 
                            dplyr::select(gs_name, entrez_gene)
            
            enrichFun <- paste0(enrichFun, "_", msigdbCategory)
            
            dt_to_pathway <- lapply(dt_to_pathway, function(x){
                names(x) <- mapIds(x = OrgDb, keys = names(x),
                        column = "ENTREZID", 
                        keytype = keyType, 
                        multiVals = "first")
                x
            })
            
            param_list$TERM2GENE <- C3_t2g
            gene_input_name <- "geneList"
            dt_to_pathway <- lapply(dt_to_pathway, sort, decreasing = TRUE)
            simplify_res <- FALSE

        }, 
        "gseWp" = {
            param_list <- list(pAdjustMethod = pAdjustMthd, 
                               pvalueCutoff = pvalueCutoff)
            wpOrg <- match.arg(wpOrg, get_wp_organisms())
            param_list$organism <- wpOrg
            dt_to_pathway <- lapply(dt_to_pathway, function(x){
                names(x) <- mapIds(x = OrgDb, column = "ENTREZID", 
                                    keys = names(x), keytype = keyType, 
                                    multiVals = "first")
                x
            })
            dt_to_pathway <- lapply(dt_to_pathway, sort, decreasing = TRUE)
            
            gene_input_name <- "geneList"
            simplify_res <- FALSE
            
        }, 
        "gseKegg" = {
            param_list <- list(keyType = keggKeyType,
                               pAdjustMethod = pAdjustMthd, 
                               pvalueCutoff = pvalueCutoff)
            
            gene_input_name <- "geneList"
            dt_to_pathway <- lapply(dt_to_pathway, function(x){
                names(x) <- mapIds(x = OrgDb, column = "ENTREZID", 
                                    keys = names(x), keytype = keyType, 
                                    multiVals = "first")
                x
            })
            dt_to_pathway <- lapply(dt_to_pathway, sort, decreasing = TRUE)
            
            param_list$organism <- keggOrg
            
            param_list$universe <- universe
            simplify_res <- FALSE
        })
        
    }
    enrich_data <- list()
    for(i in seq_along(dt_to_pathway)){
        
        param_list[[gene_input_name]] <- dt_to_pathway[[i]]
        
        comparison <- names(dt_to_pathway)[i]
        
        if(enrichFun %in% c("ora", "gsea")){
            param_list$ont <- "BP"
            gohyperBP <- do.call(enrich_fun, param_list)
            
            param_list$ont <- "CC"
            gohyperCC <- do.call(enrich_fun, param_list)
            
            param_list$ont <- "MF"
            gohyperMF <- do.call(enrich_fun, param_list)
            enrich_data[[i]] <- list("BP" = gohyperBP,"CC" = gohyperCC,"MF" = gohyperMF)    
        } else{
            enrichRes <- do.call(enrich_fun, param_list)
            enrich_data[[i]] <- list(enrichRes)
            names(enrich_data[[i]]) <- enrichFun
        }
        
        names(enrich_data)[i] <- comparison
        if(simplify_res){
            enrich_data[[i]][!sapply(enrich_data[[i]], is.null)] <- lapply(enrich_data[[i]][!sapply(enrich_data[[i]], is.null)], 
                                                                   simplify, cutoff = simplify_cutoff)
        }
    }
    
    return(enrich_data)
}



enrich_plots <- function(enrich_data, fold_change_list = NULL, 
                         node_label = "category", legend_name = "Fold Change (log2)"){
    enrich_plots_res <- lapply(seq_along(enrich_data), function(path_data_id){
        path_data <- enrich_data[[path_data_id]]
        fold_change_dt <- fold_change_list[[path_data_id]]
        enrich_plot <- lapply(path_data, function(component){
            if(is.null(component)){return(NA)}
            if(nrow(as.data.frame(component)) == 0){return(NA)}
            dot_plot <- enrichplot::dotplot(component, showCategory = 30, 
                                            font.size = 8)
            network_plot <- enrichplot::cnetplot(component, 
                categorySize = "geneNum", foldChange = fold_change_dt, 
                showCategory = 15, node_label= node_label, 
                cex.params = list(category_label = 0.4)) + 
                guides(colour=guide_colorbar(title = legend_name))
            return(list(dot_plot, network_plot))
        })
        return(enrich_plot)
    })
    names(enrich_plots_res) <- names(enrich_data)
    return(enrich_plots_res)
}


enrich_chunks <- function(enrichData_all, enrich_plot_data){
    for(i in 1:length(enrich_plot_data)){
        cat("#### ", names(enrich_plot_data)[i], "{.tabset} \n");
        for(z in 1:length(enrichData_all[[i]])){
            dt_to_table <- as.data.frame(enrichData_all[[i]][[z]])
            num_vars <- sapply(dt_to_table, is.numeric)
            
            cat("#####", names(enrichData_all[[i]])[z], " {.tabset} \n")
            if(!any(is.na(enrich_plot_data[[i]][[z]]))){
                # knit_print(enrich_plot_data[[i]][[z]][[1]])
                # knit_print(enrich_plot_data[[i]][[z]][[2]])
                lapply(enrich_plot_data[[i]][[z]], knit_print)
            }
            cat(knit_print( datatable(dt_to_table,
                                      options = list(scrollX = TRUE, pageLength = 5
                                      )) %>%
                                formatRound(columns = colnames(dt_to_table)[num_vars], digits = 3)), "\n");
            cat("\n")
        }
    }
}

calculate_enrich_fc <- function(enrich_table){
    g_ratio <- enrich_table@result$GeneRatio
    b_ratio <- enrich_table@result$BgRatio
    
    g_ratio <- sapply(g_ratio, function(x){eval(parse(text = x))})
    b_ratio <- sapply(b_ratio, function(x){eval(parse(text = x))})
    return(g_ratio/b_ratio)
}


### Annotation Plot

dmr_annotation_plot <- function(omic_dt, 
                                annotation_labels = c("Intron", "Promoter", 
                                                    "Distal Intergenic", 
                                                    "3' UTR", "5' UTR", 
                                                    "Exon", "Downstream"), 
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

### Utilities for multivariate normalization
count_vst_deseq2 <- function(omic_dt, round_dt = TRUE, multiplier = 1, ...){
    omic_dt_temp <- assay(omic_dt)
    if(round_dt){
        omic_dt_temp <- round(omic_dt_temp)
    }
    omic_dt_temp <- omic_dt_temp * multiplier
    deseq_dt <- DESeqDataSetFromMatrix(omic_dt_temp, colData = colData(omic_dt), design =  ~ group)
    assay(omic_dt) <- as.data.frame(assay(vst(deseq_dt, blind = FALSE)))
    omic_dt
}

clr_proc <- function(omic_dt, ...){
    omic_dt_temp <- omic_dt
    # if(!("relabundance" %in% names(assays(omic_dt)))){
    #     omic_dt_temp <- transformAssay(omic_dt_temp, assay.type = "counts", method = "relabundance", pseudocount = 1)
    # }
    # pseudo_c <- ifelse(any(assay(omic_dt_temp, "relabundance") == 0),
    #                     min(assay(omic_dt_temp, "relabundance")[assay(omic_dt_temp, "relabundance") != 0])*0.1,
    #                     0)
    # omic_dt_temp <- transformAssay(x = omic_dt_temp, 
    #                                 assay.type = "relabundance", 
    #                                 method = "clr", name = "clr", 
    #                                 pseudocount = pseudo_c,)
    if(!("relabundance" %in% names(assays(omic_dt)))){
        omic_dt_temp <- transformAssay(omic_dt_temp, assay.type = "counts", 
                                       method = "relabundance")
    }
    omic_dt_temp <- transformAssay(x = omic_dt_temp, 
                                    assay.type = "relabundance", 
                                    method = "rclr", name = "clr")

    assays(omic_dt) <- list("clr" = assay(omic_dt_temp, "clr"))
    omic_dt
}

transformToMvalues <- function(omic_dt, ...) {
    omic_assay <- assay(omic_dt)
    
    omic_assay[omic_assay == 1] <- 0.99
    omic_assay[omic_assay == 0] <- 0.01
    M <- log2(omic_assay / (1 - omic_assay))
    assay(omic_dt) <- M
    omic_dt
}

scale_assay <- function(omic_dt, ...){
    assay(omic_dt) <- as.data.frame(t(scale(t(assay(omic_dt)))))
    omic_dt
}
impute_knn_assay <- function(omic_dt, k = 10, ...){
    assay(omic_dt) <- as.data.frame(impute.knn(as.matrix(assay(omic_dt)), k = k)$data)
    omic_dt
}

ComBat_norm <- function(omic_dt, covariate, ...){
    assay(omic_dt) <- as.data.frame(ComBat(dat = as.matrix(assay(omic_dt)), 
                                            batch = omic_dt[[covariate]]))
    omic_dt
}


go_similarity_matrix <- function(go_list1, go_list2, semData, measure){
    sim_mat <- matrix(nrow = length(go_list1), ncol = length(go_list2))
    for(go1 in seq_along(go_list1)){
        for(go2 in seq_along(go_list2)){
            sim_mat[go1, go2] <- goSim(go_list1[go1], go_list2[go2], 
                                        semData = semData, measure = measure)
        }
    }
    rownames(sim_mat) <- go_list1
    colnames(sim_mat) <- go_list2
    return(sim_mat)
}
