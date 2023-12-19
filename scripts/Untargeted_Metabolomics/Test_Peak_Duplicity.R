feature_definitions <- as.data.frame(featureDefinitions(ms_data))
chrom_peaks <- chromPeaks(ms_data)

i <- sapply(feature_definitions$peakidx, 
            function(pk){table(table(chrom_peaks[pk, "sample"]))})

j <- unlist(i)
hist(j[names(j) == "1"])
