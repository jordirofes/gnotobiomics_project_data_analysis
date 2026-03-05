data_pca <- function(dt_list, metadata, stratificationVar, pc_x = 1, pc_y = 2, 
                    scale. = FALSE, center = FALSE, plot_name, color_vect,
                    loading_arrows = FALSE, ...){
    group_var <- metadata[[stratificationVar]]
    if(any(is.na(dt_list))){
        dt_list <- t(impute::impute.knn(as.matrix(dt_list))$data)    
    } else{
        dt_list <- t(dt_list)
    }
    
    pca_dt <- prcomp(dt_list, center = center, scale. = scale.)
    
    pca_plot(pca_dt, pc_x = pc_x, pc_y = pc_y, group_var = group_var, 
            ellipse_var = group_var, interactive = FALSE, 
            color_vect = color_vect, loading_arrows = loading_arrows) + 
        ggtitle(plot_name)
}

plot_pc <- function(plot_dt, group_var, x_lab, y_lab, loading_dt = NA, color_vect, 
                    allElipse = TRUE, ellipse_var = NA, interactive = TRUE, 
                    level = 0.95, loading_arrows = FALSE){
    pc_plot <- ggplot(plot_dt) + 
        geom_point(aes(x = .data[[colnames(plot_dt)[1]]], 
                    y = .data[[colnames(plot_dt)[2]]], colour = group_var)) + 
        theme_minimal() + xlab(x_lab) + ylab(y_lab) +
        theme(legend.title = element_blank(),
              plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
              plot.title = element_text(hjust = 0.5), 
              legend.text = element_text(size = 14), 
              axis.title.x = element_text(hjust = 0.5, size = 16),
              axis.title.y = element_text(size = 16)) + 
        geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
        
    if(!missing(color_vect)){
        if(!is.numeric(group_var)){
            pc_plot <- pc_plot + scale_color_manual(values = color_vect)
        } else{
            pc_plot <- pc_plot + scale_color_continuous("viridis")
        }    
    }
    
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
    
    if(!all(is.na(loading_dt)) & loading_arrows){
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
                    allElipse = TRUE, interactive = TRUE, level = 0.95, 
                    color_vect, loading_arrows = FALSE){
    pc_var <- pca_data$sdev^2
    pc_var_abs <- pc_var/sum(pc_var)*100
    
    if(missing(color_vect) & length(unique(group_var)) < 9){
        color_vect <- RColorBrewer::brewer.pal(length(unique(group_var)), name = "Set1")    
    } else if(missing(color_vect) & length(unique(group_var)) >= 9){
        color_vect <- RCy3::paletteColorRandom(value.count = length(unique(group_var)))
    }
    
    x_lab <- paste0("PC", pc_x, " (", round(pc_var_abs[pc_x]), "%)")
    y_lab <- paste0("PC", pc_y, " (",round(pc_var_abs[pc_y]), "%)")
    
    plot_dt <- as.data.frame(pca_data$x[,c(pc_x, pc_y)])
    
    loading_dt <- pca_data$rotation[,c(pc_x, pc_y)]
    
    loading_effect_sum <- apply(loading_dt, 1, function(x){sum(abs(x))})
    loading_names <- names(loading_effect_sum)[order(loading_effect_sum, 
                                                     decreasing = TRUE)]
    loading_dt <- loading_dt[loading_names[1:c(min(length(loading_names), 10))],]
    
    return(plot_pc(plot_dt, group_var, x_lab, y_lab, loading_dt, color_vect, 
                   allElipse, ellipse_var, interactive, level, loading_arrows))
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
                       plot_name = "Data Heatmap", show_colnames, show_rownames = FALSE, 
                       clustering_method = "average", fontsize_col = 6, fontsize_row = 6,
                       clustering_distance_rows, color_list = NA, use_cat = TRUE, 
                       ...){
    group_var <- metadata[stratificationVar]
    if(any(is.na(dt_list))){
        dt_list <- impute::impute.knn(as.matrix(dt_list))$data
    }
    dt_list <- t(dt_list)
    
    det_coef <- apply(dt_list, 2, function(x){
        sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)
    })
    det_coef <- sort(det_coef, decreasing = TRUE)
    sel_vars <- names(det_coef)[det_coef >= quantile(det_coef, cdQuantile)]
    
    dt_quantiles <- quantile_breaks(dt = dt_list[,sel_vars], 
                                    length_out = 101)
    
    if(cdQuantile != 0){
        plot_name <- paste0(plot_name, " Top ", 
                            length(sel_vars),
                            " variables (>=", cdQuantile*100, "%) CV")
    } else if(cdQuantile != 0 & is.null(plot_name)){
        plot_name <- paste("Top", length(sel_vars), "variables")
    }
    
    plot_colors <- viridis::viridis(n = length(dt_quantiles) - 1)
    
    if(use_cat){
        
        cat("\n")
        pheatmap::pheatmap(dt_list[,sel_vars], color = plot_colors,
                           annotation_row = as.data.frame(group_var),
                           show_colnames = show_colnames, fontsize = 14,
                           fontsize_col = fontsize_col, fontsize_row = fontsize_row, show_rownames = show_rownames,
                           clustering_method = clustering_method, cluster_rows = TRUE,
                           cluster_cols = TRUE, angle_col = 45,
                           clustering_distance_rows = clustering_distance_rows,
                           main = plot_name,
                           breaks = dt_quantiles, annotation_colors = color_list, 
                           annotation_legend = FALSE, annotation_names_row = FALSE)
        cat("\n")
    } else{
        pheatmap::pheatmap(dt_list[,sel_vars], color = plot_colors,
                           annotation_row = as.data.frame(group_var),
                           show_colnames = show_colnames, fontsize = 5,
                           fontsize_col = 6, fontsize_row = 6, show_rownames = show_rownames,
                           clustering_method = clustering_method, cluster_rows = TRUE,
                           cluster_cols = TRUE,
                           clustering_distance_rows = clustering_distance_rows,
                           main = plot_name,
                           breaks = dt_quantiles, annotation_colors = color_list, 
                           annotation_names_row = FALSE)
    }
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


volcanoFun <- function(comp_dt, metadata, stratificationVar, 
                        adj.pvalue = TRUE, jitterseed = 1234){
    if(adj.pvalue){
        comp_dt <- comptable[,c("logFC", "adj.P.Val")]
        dt_thr <- comp_dt$adj.P.Val <= 0.05 & abs(comp_dt$logFC) > 1
        groups <- ifelse(dt_thr, "royalblue4", "grey40")
        comp_dt$adj.P.Val <- -log10(comp_dt$adj.P.Val)
    } else{
        comp_dt <- comptable[,c("logFC", "P.Value")]
        dt_thr <- comp_dt$P.Value <= 0.05 & abs(comp_dt$logFC) > 1
        groups <- ifelse(dt_thr, "royalblue4", "grey40")
        comp_dt$P.Value <- -log10(comp_dt$P.Value)
    }
    dt_labs <- rownames(comptable)
    comp_dt <- as.list.data.frame(comp_dt)
    ggplot(comp_dt) + geom_scatter(aes(x = logFC, y = adj.p.val)) +
        geom_vline(xintercept = -1, color = "lightgrey") + 
        geom_vline(xintercept = 1, color = "lightgrey") +
        geom_hline(yintercept = -log10(0.05),color = "lightgrey") +
        ggplot2::annotate(geom = "point",x = comp_dt[[1]], 
                            y = comp_dt[[2]], colour = groups) +
        geom_text(aes(x = comp_dt[[1]][dt_thr], y = comp_dt[[2]][dt_thr], 
                    label = dt_labs[dt_thr]), 
                position = position_jitter(width = 0.1, height = 0.1, 
                                           seed = jitterseed))
    
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
                                       "gseWp", "msigdbrGSEA", "custom_ora", 
                                       "custom_gsea"), 
                         keggOrg, keggKeyType, davidOrg, wpOrg, 
                         msigdbCategory, msigdbSpc, 
                         TERM2GENE, TERM2NAME, ...){
    
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
                              "msigdbrGSEA" = GSEA,
                              "custom_ora" = enricher,
                              "custom_gsea" = GSEA
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
            
            # wpOrg <- match.arg(wpOrg, get_wp_organisms())
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
            # wpOrg <- match.arg(wpOrg, get_wp_organisms())
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
        },
        "custom_ora" = {
            param_list <- list(pAdjustMethod = pAdjustMthd, 
                               pvalueCutoff = pvalueCutoff)
            
            gene_input_name <- "gene"
            dt_to_pathway <- lapply(dt_to_pathway, mapIds, x = OrgDb, 
                                    column = "ENTREZID", keytype = keyType, 
                                    multiVals = "first")
            param_list$TERM2GENE <- TERM2GENE
            param_list$TERM2NAME <- TERM2NAME
            
            param_list$universe <- universe
            param_list$qvalueCutoff  <- qvalueCutoff
            simplify_res <- FALSE
        },
        "custom_gsea" = {
            param_list <- list(pAdjustMethod = pAdjustMthd, 
                               pvalueCutoff = pvalueCutoff)
            
            gene_input_name <- "geneList"
            dt_to_pathway <- lapply(dt_to_pathway, function(x){
                names(x) <- mapIds(x = OrgDb, column = "ENTREZID", 
                                   keys = names(x), keytype = keyType, 
                                   multiVals = "first")
                x
            })
            dt_to_pathway <- lapply(dt_to_pathway, sort, decreasing = TRUE)
            
            param_list$TERM2GENE <- TERM2GENE
            param_list$TERM2NAME <- TERM2NAME
            
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
            enrich_data[[i]] <- list("BP" = gohyperBP, 
                                     "CC" = gohyperCC,
                                    "MF" = gohyperMF)    
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

enrich_plots_list <- function(enrich_list, fold_change_list = NULL, 
                         node_label = "category", legend_name = "Fold Change (log2)"){
    enrich_plots_res <- lapply(seq_along(enrich_list), function(path_data_id){
        path_data <- enrich_list[[path_data_id]]
        if(!is.null(fold_change_list)){
            fold_change_dt <- fold_change_list[[path_data_id]]
        } else{
            fold_change_dt <- NULL
        }
        
        enrich_plot <- lapply(path_data, enrich_plots, 
                              fold_change_dt = fold_change_dt, 
                              node_label = node_label, legend_name = legend_name)
        return(enrich_plot)
    })
    names(enrich_plots_res) <- names(enrich_list)
    return(enrich_plots_res)
}

enrich_plots <- function(enrich_data, fold_change_dt, node_label, legend_name){
    if(is.null(enrich_data)){return(NA)}
    if(nrow(as.data.frame(enrich_data)) == 0){return(NA)}
    
    plot_dt_list <- list()
    
    plot_dt_list$dot_plot <- enrichplot::dotplot(enrich_data, showCategory = 30, 
                                    font.size = 8)
    plot_dt_list$network_plot <- enrichplot::cnetplot(enrich_data, 
                                        categorySize = "geneNum", 
                                        foldChange = fold_change_dt, 
                                        showCategory = 15, node_label= node_label, 
                                        cex.params = list(category_label = 0.4)) + 
        guides(colour=guide_colorbar(title = legend_name))
    
    go_types <- c("BP", "CC", "MF")
    
    ont_slot <- grep("ontology|setType", slotNames(enrich_data), value = TRUE) 
    
    if(length(ont_slot) & slot(enrich_data, ont_slot) %in% go_types & nrow(enrich_data) >= 5){
        enrich_data <- pairwise_termsim(x = enrich_data)
        
        plot_dt_list$tree_plot <- enrichplot::treeplot(enrich_data, 
                                                       cex.params = list(category_label = 0.4))
        plot_dt_list$emaplot <- enrichplot::emapplot(enrich_data, cex.params = list(category_label = 0.4),
                                                    cluster.params = list(cluster = TRUE, 
                                                                        legend = TRUE))
    }
    
    return(plot_dt_list)
}

enrich_chunks <- function(enrichData_all, enrich_plot_data){
    
    if(is.null(enrichData_all) || length(enrichData_all) == 0){return(NULL)}
    for(i in 1:length(enrich_plot_data)){
        
        if(is.null(enrichData_all[[i]]) || length(enrichData_all[[i]]) == 0){next}
        
        cat("#### ", names(enrich_plot_data)[i], "{.tabset} \n");
        for(z in 1:length(enrichData_all[[i]])){
            dt_to_table <- as.data.frame(enrichData_all[[i]][[z]])
            if(is.null(dt_to_table) || nrow(dt_to_table) == 0){next}
            num_vars <- sapply(dt_to_table, is.numeric)
            
            cat("#####", names(enrichData_all[[i]])[z], " {.tabset} \n")
            if(!any(is.na(enrich_plot_data[[i]][[z]]))){
                lapply(enrich_plot_data[[i]][[z]], knit_print)
            }
            cat(knit_print(datatable(dt_to_table,
                                    options = list(scrollX = TRUE, 
                                                    pageLength = 5)) %>% 
                                formatRound(columns = colnames(dt_to_table)[num_vars], 
                                            digits = 3)), "\n");
            cat("\n")
        }
    }
}

calculate_enrich_fc <- function(enrich_table){
    if(is(enrich_table, "data.frame")){
        g_ratio <- enrich_table$GeneRatio
        b_ratio <- enrich_table$BgRatio    
    } else{
        g_ratio <- enrich_table@result$GeneRatio
        b_ratio <- enrich_table@result$BgRatio    
    }
    
    
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
count_vst_deseq2 <- function(omic_dt, round_dt = TRUE, multiplier = 1, design = ~ group + 0, ...){
    omic_dt_temp <- as.matrix(assay(omic_dt))
    if(round_dt){
        mode(omic_dt_temp) <- "integer"
    }
    if(any(is.na(omic_dt_temp))){
        omic_dt_temp <- impute.knn(omic_dt_temp, k = 10)$data
        mode(omic_dt_temp) <- "integer"
    }
    
    omic_dt_temp <- omic_dt_temp * multiplier
    deseq_dt <- DESeqDataSetFromMatrix(omic_dt_temp, colData = colData(omic_dt), design = design)
    assay(omic_dt) <- as.data.frame(assay(varianceStabilizingTransformation(deseq_dt, blind = FALSE)))
    omic_dt
}

count_rlog_deseq2 <- function(omic_dt, round_dt = TRUE, multiplier = 1, ...){
    omic_dt_temp <- assay(omic_dt)
    if(round_dt){
        omic_dt_temp <- round(omic_dt_temp)
    }
    omic_dt_temp <- omic_dt_temp * multiplier
    deseq_dt <- DESeqDataSetFromMatrix(omic_dt_temp, colData = colData(omic_dt), design =  ~ group + 0)
    assay(omic_dt) <- as.data.frame(assay(rlog(deseq_dt, blind = FALSE)))
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
                                       method = "relabundance", pseudocount = 1)
    }
    omic_dt_temp <- transformAssay(x = omic_dt_temp, 
                                   assay.type = "relabundance", 
                                   method = "clr", name = "clr")
    
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


network_heatmap_fun <- function(plot_dt, labels_dt, omic_color_dict, 
                                x_colors, plot_title){
    dt_quantiles <- quantile_breaks(dt = plot_dt[,"value"], 
                                    length_out = 101)
    
    p <- ggplot(plot_dt, aes(x = sample_name, fill = value, 
                             y = var_id)) + 
        geom_tile() + 
        theme_minimal() + 
        scale_fill_viridis_c(values = scales::rescale(dt_quantiles)) +
        theme(axis.text.x = element_text(angle = 90, 
                                         colour = x_colors),
              axis.text.y = element_text(size = 5),
              legend.position = "right", 
              legend.justification = "top", 
              legend.title = element_blank(), 
              plot.title = element_text(hjust = 0.5),
              legend.spacing.y = unit(0.1, "cm")) + 
        ggtitle(plot_title) + 
        xlab(NULL) + ylab(NULL) + 
        ggnewscale::new_scale_fill() + 
        geom_tile(inherit.aes = FALSE, 
                  aes(x = -1, y = var_id, 
                      fill = omic_name), data = labels_dt) +
        scale_fill_manual(values = omic_color_dict) 
    
}

network_lineplot_fun <- function(plot_dt, labels_dt, omic_color_dict, x_colors, 
                                 legend_colors, plot_title){
    p <- ggplot(plot_dt, aes(x = sample_name, y = value, color = omic_name,
                             group = var_id, shape = var_id)) +
        geom_point() + geom_smooth(aes(color = omic_name, fill = omic_name),
                                   se = FALSE, span = 0.5, alpha = 0.1) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, colour = x_colors),
              legend.position = "right",
              legend.justification = "top",
              legend.title = element_blank(),
              plot.title = element_text(hjust = 0.5),
              legend.spacing.y = unit(0.1, "cm")) +
        ggtitle(plot_title) +
        xlab(NULL) + ylab("Normalized Values") +
        scale_color_manual(values = omic_color_dict)  +
        guides(shape = guide_legend(label.theme = element_text(size = 7.5),
                                    order = 2,  byrow = TRUE,
                                    override.aes = list(
                                        color = legend_colors
                                    ),
                                    keyheight = 0.05),
               col = guide_legend(order = 1), fill = "none") +
        scale_shape_manual(values = rep(16, length(unique(plot_dt$var_id))))
    return(p)
}

plot_nw_term <- function(omic_dt, feature_ids, term_ids, group_var = "group", 
                         sample_var = "samplename", feature_ids_dict = NULL, 
                         omic_names_dict = NULL, omic_color_dict, 
                         group_color_dict = NULL,
                         var_order_col = "Factor_2", plot_title, 
                         plot_type = c("Heatmap", "Lineplot")){
    
    # Preparing some dictionary data
    if(!is.null(params$omic_names_dict)){
        names(omic_color_dict) <- omic_names_dict[as.character(names(omic_color_dict))]
    }
    
    group_var_dict <- colData(omic_dt)[[group_var]]
    names(group_var_dict) <- colData(omic_dt)[[sample_var]]
    
    # Preparing melted plot data for ggplot
    plot_dt <- lapply(names(experiments(omic_dt)), function(omic_id){
        omic_filt <- omic_dt[[omic_id]]
        found <- paste0(omic_id, ".", rownames(omic_filt)) %in% feature_ids
        if(any(found)){
            return(omic_filt[found,])
        } else{
            return(NA)
        }
    })
    names(plot_dt) <- names(experiments(omic_dt))
    plot_dt <- plot_dt[!is.na(plot_dt)]
    
    # Melt data to plot and add omic name and group variable
    plot_dt <- lapply(names(plot_dt), function(omic_exp_name){
        omic_exp <- plot_dt[[omic_exp_name]]
        melted_dt <- melt(as.matrix(assay(omic_exp)))
        
        colnames(melted_dt)[1:2] <- c("var_id", "sample_name")
        
        melted_dt$factor_weight <- rowData(omic_exp)[melted_dt$var_id, var_order_col]
        melted_dt$factor_weight_pos <- ifelse(melted_dt$factor_weight >= 0,
                                              "Positive", "Negative")
        
        # melted_dt$var_id <- factor(gsub(paste0(omic_exp_name, "."), "", melted_dt$var_id))
        melted_dt$var_id <- factor(paste(omic_exp_name, 
                                         melted_dt$var_id, sep = "."))
        melted_dt$omic_name <- factor(omic_exp_name)
        
        melted_dt$group <- factor(group_var_dict[as.character(melted_dt$sample_name)])
        
        if(!is.null(omic_names_dict)){
            melted_dt$omic_name <- factor(omic_names_dict[as.character(melted_dt$omic_name)])
        }
        if(!is.null(feature_ids_dict)){
            melted_dt$var_id <- factor(feature_ids_dict[[omic_exp_name]][as.character(melted_dt$var_id)])
        }
        
        melted_dt
    })
    
    plot_dt <- do.call(rbind, plot_dt)
    
    # plot_dt$var_id <- factor(abbreviate(as.character(plot_dt$var_id), 75, ))
    
    plot_dt$var_id <- factor(plot_dt$var_id, levels = levels(plot_dt$var_id), 
                             labels = base::make.unique(str_trunc(levels(plot_dt$var_id), 75)))
    
    plot_dt$var_id <- factor(str_wrap(plot_dt$var_id, width = 40))
    
    
    # Sample factor re-ordering and colors
    sample_levels <- levels(plot_dt$sample_name)
    reordered_sample_levels <- sample_levels[order(group_var_dict[as.character(sample_levels)])]
    
    plot_dt$sample_name <- reorder(plot_dt$sample_name, 
                                   new.order = reordered_sample_levels)
    
    x_colors <- group_color_dict[as.character(group_var_dict[levels(plot_dt$sample_name)])]
    
    # Variable re-ordering and colors
    labels_dt <- distinct(plot_dt[,c("var_id", "omic_name", "factor_weight")])
    # labels_dt <- labels_dt[order(labels_dt$omic_name),]
    
    var_omic_dict <- labels_dt$omic_name
    names(var_omic_dict) <- labels_dt$var_id
    
    var_levels <- levels(plot_dt$var_id)
    
    # reordered_var_levels <- var_levels[order(var_omic_dict[as.character(var_levels)])]
    reordered_var_levels <- labels_dt$var_id[order(labels_dt$factor_weight, 
                                                   decreasing = TRUE)]
    
    
    plot_dt$var_id <- reorder(plot_dt$var_id, 
                              new.order = reordered_var_levels)
    
    legend_colors <- omic_color_dict[as.character(sapply(levels(plot_dt$var_id), function(x){
        labels_dt$omic_name[as.character(labels_dt$var_id) == x]
    }))]
    
    names(legend_colors) <- levels(plot_dt$var_id)
    labels_dt$omic_colors <- legend_colors[labels_dt$var_id]
    
    legend_colors
    
    # Plotting
    plot_type <- match.arg(plot_type)
    
    switch(plot_type, 
           "Heatmap" = {
               p <- network_heatmap_fun(plot_dt = plot_dt, labels_dt = labels_dt, 
                                        omic_color_dict = omic_color_dict,
                                        x_colors = x_colors, plot_title = plot_title)
           }, "Lineplot" = {
               p <- network_lineplot_fun(plot_dt = plot_dt, labels_dt = labels_dt, 
                                         legend_colors = legend_colors, 
                                         omic_color_dict = omic_color_dict,
                                         x_colors = x_colors, plot_title = plot_title)
           })
    
    
    return(p)
}


setReadable_custom <- function(enrich_res, featureDict, dt_type = "ora"){
    
    dt_type <- switch(dt_type, "ora" = "geneID", "gsea" = "core_enrichment")
    if((!is(enrich_res, "enrichResult") & !is(enrich_res, "gseaResult")) || nrow(enrich_res@result) == 0){return(enrich_res)}
    geneID_list <- strsplit(enrich_res@result[[dt_type]], "/")
    enrich_res@result[[dt_type]] <- sapply(geneID_list, function(pathway_ids){
        
        found <- names(featureDict)[sapply(featureDict, function(ids_dict){any(ids_dict %in% pathway_ids)})]
        
        paste(unique(unlist(found)), collapse = "/")
    })
    enrich_res
}

enrich_multi_omic_clusters <- function(graph, feature_dict, 
                                       membership_attribute, universe, 
                                       min_n = 25){
    
    membership_vector <- vertex_attr(graph, membership_attribute)
    clu_vertex <- lapply(unique(membership_vector), function(cl_id){
        names(V(graph))[membership_vector == cl_id]
    })
    names(clu_vertex) <- unique(membership_vector)
    clu_vertex <- clu_vertex[lengths(clu_vertex) >= min_n]
    rna_seq_tag <- "RNA_Sequencing\\."
    dna_seq_tag <- "WGB_DNA_Methylation\\."
    
    general_tag <- paste0(rna_seq_tag, "|", dna_seq_tag)
    
    enrich_res <- lapply(clu_vertex, function(clu_ids){
        merged_ids <- grep(general_tag, clu_ids, value = TRUE)
        merged_ids <- feature_dict[merged_ids]
        enrich_go <- clusterProfiler::enrichGO(
            gene = merged_ids, OrgDb = org.Mm.eg.db,
            keyType = "ENTREZID", ont = "BP",
            pvalueCutoff = 0.05, pAdjustMethod = "fdr",
            universe = universe)
        enrich_go <- clusterProfiler::simplify(enrich_go)
        
        enrich_go <- setReadable_custom(enrich_go, featureDict = merged_ids)
        enrich_go
    })
    enrich_res
}



lollipop_compare_enrich <- function(enrich_1, enrich_2, 
                                    dt_colors = c("steelblue1", "indianred1"),
                                    plot_title = ""){
    plot_data <- rbind(enrich_1, enrich_2)
    
    if (dim(plot_data)[1] == 0) {
        return(NULL)
    }
    plot_data$Description <- gsub(" - Mus musculus \\(house mouse)", replacement = "", x = plot_data$Description)
    ctrs_dt <- lapply(unique(plot_data$type), function(x) {
        ctrs_dt <- plot_data[plot_data$type == x, ]
        ctrs_dt$Description[order(ctrs_dt$enrich_log2fc, 
                                  decreasing = TRUE)]
    })
    names(ctrs_dt) <- unique(plot_data$type)
    if (length(ctrs_dt) == 0) {
        return(NULL)
    }
    if (length(ctrs_dt) == 1) {
        ctrs_dt[[2]] <- character(0)
        dt_colors[2] <- "white"
    }
    set_1 <- setdiff(ctrs_dt[[2]], ctrs_dt[[1]])
    set_2 <- setdiff(ctrs_dt[[1]], ctrs_dt[[2]])
    set_intersect <- intersect(ctrs_dt[[1]], ctrs_dt[[2]])
    ordered_dt <- c(rev(set_1), rev(set_2), rev(set_intersect))
    
    plot_data$Description_factor <- factor(plot_data$Description, 
                                           levels = ordered_dt, ordered = TRUE)
    
    color_ctrs <- dt_colors
    names(color_ctrs) <- unique(plot_data$type)

    fc_ceiling <- ceiling(max(plot_data$enrich_log2fc)*2)/2
    
    p_d <- ggplot(plot_data, aes(x = enrich_log2fc, y = Description_factor)) + 
        geom_point(aes(color = type), size = 8) + 
        geom_segment(aes(x = 0, xend = enrich_log2fc, y = Description_factor,
                         yend = Description_factor, color = type), 
                     alpha = 0.9, size = 1.5) + 
        geom_vline(xintercept = 0) + 
        theme_minimal() + 
        theme(legend.position = "bottom",
              plot.title = element_text(size = 22, 
                                        hjust = 0.5), 
              axis.text.y.left = element_text(size = 22, 
                                              margin = margin(t = 100, b = 100, 
                                                              r = 10, l = 10)),
              legend.text = element_text(size = 22),
              legend.title = element_text(size = 22),
              axis.title.x  = element_text(size = 22),
              axis.title.y = element_text(size = 22), 
              axis.text.x = element_text(size = 22)) + 
        xlab("ORA Fold-Enrich") + 
        ylab("KEGG Term Description") + 
        scale_color_manual(name = "Contrasts", 
                           values = color_ctrs) + 
        ggtitle(plot_title) + 
        scale_x_continuous(breaks = seq(0, fc_ceiling, 0.5),
                           limits = c(0, fc_ceiling))
    
    venn_dt <- lapply(unique(plot_data$type), function(x) {
        plot_data$Description_factor[plot_data$type == 
                                         x]
    })
    names(venn_dt) <- unique(plot_data$type)
    rotation_degree <- ifelse(length(set_1) < length(set_2), 
                              180, 0)
    venn_plot <- draw.pairwise.venn(area1 = length(ctrs_dt[[1]]), 
                                    area2 = length(ctrs_dt[[2]]), cross.area = length(set_intersect), 
                                    fill = dt_colors, lty = "blank", cex = 0, label.col = c(dt_colors[1], 
                                                                                            "thistle", dt_colors[2]), 
                                    rotation.degree = rotation_degree, 
                                    ind = FALSE)
    list("plot_p" = p_d, "venn_diagram" = venn_plot)
    # plot_grid(p_d, plot_grid(venn_plot, nrow = 8, ncol = 1), 
    #           rel_widths = c(7, 1))
}