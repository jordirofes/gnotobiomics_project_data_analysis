library(readxl)

lc_ms_meta <- read_excel(path = "./worklist_untargeted.xlsx", sheet = 3)
gf_lima_meta <- read.csv("./MERGED_gfcv_LIMA_metadata.csv")

lc_ms_meta <- lc_ms_meta[-c(1, 4:12, 105, 106),]
lc_ms_meta[1:2,1] <- c("Blank_00", "Blank_01")
b <- do.call(rbind,lapply(strsplit(gf_lima_meta$samplename[1:36], split = ""), function(x){ paste(x[1:3], collapse = "")}))

dt_to_join <- lapply(lc_ms_meta$`injection order`, function(x){
    found <- gf_lima_meta[which(sapply(b, function(y){grepl(pattern = y, x = x)})),]
    if(nrow(found) == 0){
        found <- as.data.frame(lapply(found, function(x){x <- NA}))
    }
    return(found)
})
dt_to_join <- do.call(rbind, dt_to_join)
rownames(dt_to_join) <- 1:95
dt_to_join$injection_order <- lc_ms_meta$`injection order`
dt_to_join[dt_to_join$injection_order == "20_686_INT", 1:6] <- gf_lima_meta[32,]

for(x in seq_along(dt_to_join$injection_order)){
    dt <- dt_to_join$injection_order[x]
    found_qc <- grepl(pattern = "QC", x = dt) 
    found_blank <- grepl(pattern = "Blank", x = dt)
    if(found_qc){
        dt_to_join[x,1:6] <- c(dt_to_join[x,7], rep("RP_QC", 4), dt_to_join[x,7])     
    }
    if(found_blank){
        dt_to_join[x,1:6] <- c(dt_to_join[x,7], rep("RP_Blank", 4), dt_to_join[x,7])
    }

}
v1 <- substr(gf_lima_meta$samplename, start = 1, stop = 2)
v2 <- substr(gf_lima_meta$samplename, start = 3, stop = 4)
gf_lima_dt <- data.frame(v1, v2)

v1 <- dt_to_join$injection_order
v1 <- lapply(strsplit(v1, split = "_"), function(x){
    if(x[2] == "CV" | x[2] == "GF"){
        return(x)
    } else{
        return(rep(0, 4))
    }
})
v1 <- lapply(v1, function(x){
    if(any(x != 0)){
        group <- substr(x[2], start = 1, stop = 1)
        gender_num <- x[3]
        batch <- ifelse(x[4] == "a", 1, 2)    
        return(paste0(group, gender_num, batch, collapse = ""))
    } else{
        return(0)
    }
    
})


for(x in seq_along(v1)){
    dt <- v1[[x]]
    found <- dt == gf_lima_meta$samplename
    if(any(found)){
        dt_to_join[x, 1:6] <- gf_lima_meta[found, 1:6]
    }
}
dt_to_join[dt_to_join$injection_order == "3_CV_F1_b", 1:6] <- c("CF12", "CVF", 2, "F", "C", "CF12")
dt_to_join[dt_to_join$injection_order == "29_GF_F5_b", 1:6] <- c("GF52", "GFF", 2, "F", "G", "GF52")


write.csv(dt_to_join, file = "D:/Untargeted_LIMMA_GF/RP_Positive/LC_MS_metadata.csv", row.names = FALSE)
# for(x in seq_along(dt_to_join$injection_order)){
#     dt <- dt_to_join$injection_order[x]
#     found <- strsplit(x = dt, split = "_")[[1]]
#     if(substr(found[2], start = 1, stop = 1) == "C"){
#         gf_lima_dt$v1 == "CV" & gf_lima_dt$v2 == found[3]
#     }
#     if(substr(found[2], start = 1, stop = 1) == "G"){
#         
#     }
#     print(found)
# }

## Change HP lab for RP lab or another

gf_lima_meta_hn <- read.csv("D:/Untargeted_LIMMA_GF/Hilic_Negative/LC_MS_metadata.csv")
gf_lima_meta_hn <- gf_lima_meta
gf_lima_meta_hn_1 <- lapply(gf_lima_meta_hn, function(x){
    x <- sub(pattern = "HP", replacement = "HN", x = x)
    return(x)
})
gf_lima_meta_hn_1 <- as.data.frame(gf_lima_meta_hn_1)

write.csv(gf_lima_meta_hn_1, file = "D:/Untargeted_LIMMA_GF/Hilic_Negative/LC_MS_metadata.csv")
