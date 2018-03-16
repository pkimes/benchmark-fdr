
download.file("https://www.dropbox.com/s/8hf0ccl27ircodz/dsdObjects.RData?dl=1",
              destfile="data/dsdObjects.RData")
load("data/dsdObjects.RData")
dsd1 <- dsd1[,colData( dsd1 )$tissue %in% 
  c("Brain - Nucleus accumbens (basal ganglia)", "Brain - Putamen (basal ganglia)")]
colData(dsd1)$tissue <- droplevels( colData( dsd1 )$tissue )
colData(dsd1)$individual <- droplevels( colData( dsd1 )$individual )
dsd1 <- dsd1[,colData(dsd1)$sex == "female"]
colData(dsd1)$sex <- droplevels( colData( dsd1 )$sex )

colData(dsd1)$tissue <- factor(gsub(" |Brain -|\\(basal ganglia\\)", "", 
                                    as.character( colData(dsd1)$tissue )))
design(dsd1) <- ~ tissue
saveRDS(dsd1, file="data/brain-counts.rds" )

