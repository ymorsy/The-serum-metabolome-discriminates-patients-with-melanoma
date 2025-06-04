pacman::p_load(MetaboAnalystR, rio, tidyverse, glue)

# create the results directories
dirs <- c("00_Data scaling", "01_PCA", "02_PLSDA", "03_Correlation", "04_Stats", "05_Heatmap")
lapply(dirs, dir.create)
df_name <- "met_anno_lab.txt"

mSet <- InitDataObjects("pktable", "stat", FALSE)
mSet <- Read.TextData(mSet, df_name, "rowu", "disc")
summary(mSet$dataSet$cls)

# mSet$dataSet$cls <- factor(mSet$dataSet$cls,  levels = c("Control","Melanoma_Baseline"))
# summary(mSet$dataSet$cls)
capture.output(summary(mSet$dataSet$cls), file = "summary.txt")
mSet <- SanityCheckData(mSet)
mSet <- ReplaceMin(mSet)
mSet <- ContainMissing(mSet)
mSet <- FilterVariable(mSet, "iqr", "F", 25)
mSet <- PreparePrenormData(mSet)
mSet <- Normalization(mSet, "MedianNorm", "NULL", "AutoNorm", ratio = FALSE, ratioNum = 20)

data_scaled <- mSet[["dataSet"]][["norm"]]
data_scaled <- rownames_to_column(data_scaled, "Samples")
data_scaled <- data.table::transpose(data_scaled, make.names = 1, keep.names = "Metabolites")

write.table(data_scaled, "00_Data scaling/data_scaled.txt", row.names = F, sep = "\t")

# plot the normalized features and samples

mSet <- PlotNormSummary(mSet, "00_Data scaling/Features_normalized", "pdf", 100, width = NA)
mSet <- PlotSampleNormSummary(mSet, "00_Data scaling/Samples_normalized", "png", 100, width = NA)


# perform PCA and PLSDA


mSet <- PCA.Anal(mSet)
# green , red , blue , purple, aqua , fushia , Yellow "#00FFFF","#FF00FF" , "#FFFF00"
# colVec<-c("#008000", "#FF0000","#0000FF","#800080","#00FFFF","#FF00FF" , "#FFFF00")
colVec <- c("blue3", "red3")
# shapeVec<-c(19,1)
tryCatch(
    {
        mSet <- UpdateGraphSettings(mSet, colVec)
    },
    error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
    }
)

mSet <- PlotPCA2DScore(mSet, "01_PCA/PCA_", "pdf", 100, width = NA, 1, 2, 0.95, 0, 0)
mSet <- PlotPCA2DScore(mSet, "01_PCA/PCA2_", "pdf", 100, width = NA, 1, 2)




mSet <- PLSR.Anal(mSet, reg = TRUE)
mSet <- PlotPLS2DScore(mSet, "02_PLSDA/PCA2_", "pdf", 100, width = NA, 1, 2)
mSet <- PlotPLS2DScore(mSet, "02_PLSDA/PCA_", "pdf", 100, width = NA, 1, 2, 0.95, 0, 0)
# Perform correlation
# setwd("03_Correlation")

mSet <- PlotCorrHeatMap(mSet, "03_Correlation/corr_sample_detail_", "pdf", 72, width = NA, "row", "pearson", "bwm", "detail", F, F, "0")
mSet <- PlotCorrHeatMap(mSet, "03_Correlation/corr_sample_overview_", "pdf", 72, width = NA, "row", "pearson", "bwm", "overview", F, F, "0")

mSet <- PlotCorrHeatMap(mSet, "03_Correlation/corr_feature_detail_", "pdf", 72, width = NA, "col", "pearson", "bwm", "detail", F, F, "0")
mSet <- PlotCorrHeatMap(mSet, "03_Correlation/corr_feature_overview_", "pdf", 72, width = NA, "col", "pearson", "bwm", "overview", F, F, "0")


# mSet<-FeatureCorrelation(mSet, "pearson", "NaBu")
# mSet<-PlotCorr(mSet, "NaBu_pattern_", "pdf", 100, width=NA)
# setwd("..")

# Ttest_test
setwd("04_Stats")
Ttests.Anal(mSet, )
mSet <- Ttests.Anal(mSet, nonpar = F, 0.05, paired = FALSE, equal.var = T, "fdr")

mSet <- mSet <- Ttests.Anal(mSet, nonpar = F, 0.05, paired = FALSE, equal.var = T, "fdr", all_results = TRUE)


if (file.exists("t_test.csv")) {
    res <- import("t_test.csv") %>%
        rename(ID = V1)
    num <- nrow(res)

    anno <- import("../anno_hmdb_ids.xlsx")

    res2 <- left_join(res, anno) %>%
        relocate(Name, .after = ID)
    export(res2, "Results.xlsx")
    setwd("..")

    # heatmap


    if (num < 50) {
        mSet <- PlotSubHeatMap(mSet, glue("05_Heatmap/Top_{num}_heatmap_overview_clustered_"), "pdf", 72, width = NA, "norm", "row", "euclidean", "ward.D", "bwm", "tanova", num, "overview", T, T, T, F)
        mSet <- PlotSubHeatMap(mSet, glue("05_Heatmap/Top_{num}_heatmap_overview_"), "pdf", 72, width = NA, "norm", "row", "euclidean", "ward.D", "bwm", "tanova", num, "overview", F, T, T, F)
    } else {
        mSet <- PlotSubHeatMap(mSet, "05_Heatmap/Top_50_heatmap_overview_clustered_", "pdf", 72, width = NA, "norm", "row", "euclidean", "ward.D", "bwm", "tanova", 50, "overview", T, T, T, F)
        mSet <- PlotSubHeatMap(mSet, "05_Heatmap/Top_50_heatmap_overview_", "pdf", 72, width = NA, "norm", "row", "euclidean", "ward.D", "bwm", "tanova", 50, "overview", F, T, T, F)


        mSet <- PlotSubHeatMap(mSet, glue("05_Heatmap/Top_{num}_heatmap_overview_clustered_"), "pdf", 72, width = NA, "norm", "row", "euclidean", "ward.D", "bwm", "tanova", num, "overview", T, T, T, F)
        mSet <- PlotSubHeatMap(mSet, glue("05_Heatmap/Top_{num}_heatmap_overview_"), "pdf", 72, width = NA, "norm", "row", "euclidean", "ward.D", "bwm", "tanova", num, "overview", F, T, T, F)
        mSet <- PlotSubHeatMap(mSet, glue("05_Heatmap/Top_{num}_heatmap_detail_clustered_"), "pdf", 72, width = NA, "norm", "row", "euclidean", "ward.D", "bwm", "tanova", num, "detail", T, T, T, F)
        mSet <- PlotSubHeatMap(mSet, glue("05_Heatmap/Top_{num}_heatmap_detail_"), "pdf", 72, width = NA, "norm", "row", "euclidean", "ward.D", "bwm", "tanova", num, "detail", F, T, T, F)
    }
}
