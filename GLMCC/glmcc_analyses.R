library(pzfx)
library(PairedData)
library(ggplot2)
library(ggpubr)
library(plot.matrix)
library(cluster)
library(RColorBrewer)
library(psych)
library(network)
library(igraph)
library(readxl)
library(circlize)
library(ComplexHeatmap)
library("ggvenn")

days <- c("171019", "171207", "171208", "171213", "180110", "180111", "180131", "180221", "180228", 
          "180302", "180419", "180420", "180423");
daynight <- c(1,2,2,2,2,2,1,2,2,1,1,1,1)


regions_match_s2 <- read_xlsx("Library/CloudStorage/OneDrive-CNR/tim_brown/regions_match_revision_May2023.xlsx", 
                                  sheet = "Sheet4")

colnames(regions_match_s2) <- c("Rec_id", "Regions", "Neurons")
regions_match_s2$Rec_id <- as.factor(regions_match_s2$Rec_id)
regions_match_s2$Regions <- as.factor(regions_match_s2$Regions)
png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/macroregions_pie_v2.png", width = 2560, height = 1440, units = "px")
ggplot(regions_match_s2, aes(x="", y=Neurons, fill=Regions)) + geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +  theme_void() + theme(text = element_text(size = 20))
dev.off()

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/macroregions_dist_v2.png", width = 2560, height = 1440, units = "px")
ggboxplot(regions_match_s2, x="Regions", y="Neurons", fill="Regions") + theme(text = element_text(size = 20)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + xlab("") + theme(legend.title = element_blank()) + 
  theme(legend.position = "none")
dev.off()

df_171019 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171019_deg_wei")
regions_171019 <- character(nrow(df_171019));
regions_171019[seq(1,114)] = "ventromedial thalamus"
regions_171019[c(115, 116)] = "paraventricular hypothalamus"
regions_171019[seq(117,126)] = "ventromedial thalamus"
nodeid <- rep(1:nrow(df_171019),2);
condition <- c(rep("Ongoing",nrow(df_171019)), rep("LightON",nrow(df_171019)));
dat_171019 <- data.frame(c(df_171019$Ongoing, df_171019$LightON), condition, nodeid, rep(regions_171019,2));
colnames(dat_171019) <- c("NodeDegree", "condition", "nodeid", "region")
dat_171019$condition <- factor(dat_171019$condition, levels = c("Ongoing", "LightON"))
dat_171019$region <- factor(dat_171019$region)
dat_171019$NodeDegree <- as.numeric(sub(",", ".", c(df_171019$Ongoing, df_171019$LightON), fixed = TRUE))
model_171019 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_171019)
res_aov_171019 <- aov(model_171019)
summary(res_aov_171019)
cm=compare_means(NodeDegree ~ condition, data = dat_171019, group.by = "region")

tmpdf_degree = aggregate(NodeDegree ~ region + condition, data = dat_171019, mean);
tmpdf_degree$expID = rep("171019",1,length(levels(as.factor(tmpdf_degree$region)))*2);
df_degree = tmpdf_degree;

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/degree_2way_171019_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171019, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171019") + ylab("Node Strength")
dev.off()

df_171019 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171019_clus_wei")
nodeid <- rep(1:nrow(df_171019), 2);
condition <- c(rep("Ongoing",nrow(df_171019)), rep("LightON", nrow(df_171019)));
dat_171019 <- data.frame(c(df_171019$Ongoing, df_171019$LightON), condition, nodeid, rep(regions_171019,2));
colnames(dat_171019) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_171019$condition <- factor(dat_171019$condition, levels = c("Ongoing", "LightON"))
dat_171019$ClusCoeff <- as.numeric(sub(",", ".", c(df_171019$Ongoing, df_171019$LightON), fixed = TRUE))
model_171019 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_171019)
res_aov_171019 <- aov(model_171019)
summary(res_aov_171019)
cm=compare_means(ClusCoeff ~ condition, data = dat_171019, group.by = "region")

tmpdf_cluscoeff = aggregate(ClusCoeff ~ region + condition, data = dat_171019, mean);
tmpdf_cluscoeff$expID = rep("171019",1,length(levels(as.factor(tmpdf_cluscoeff$region)))*2);
df_cluscoeff = tmpdf_cluscoeff;

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/cluscoeff_2way_171019_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171019, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171207") + ylab("Clustering Coefficient")
dev.off()


df_171019 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171019_localeff_wei")
nodeid <- rep(1:nrow(df_171019), 2);
condition <- c(rep("Ongoing",nrow(df_171019)), rep("LightON", nrow(df_171019)));
dat_171019 <- data.frame(c(df_171019$Ongoing, df_171019$LightON), condition, nodeid, rep(regions_171019,2));
colnames(dat_171019) <- c("LocEff", "condition", "nodeid", "region")
dat_171019$condition <- factor(dat_171019$condition, levels = c("Ongoing", "LightON"))
dat_171019$LocEff <- as.numeric(sub(",", ".", c(df_171019$Ongoing, df_171019$LightON), fixed = TRUE))
model_171019 <- lm(formula = "LocEff ~ condition * region + (1|nodeid)", data=dat_171019)
res_aov_171019 <- aov(model_171019)
summary(res_aov_171019)
cm=compare_means(LocEff ~ condition, data = dat_171019, group.by = "region")

tmpdf_loceff = aggregate(LocEff ~ region + condition, data = dat_171019, mean);
tmpdf_loceff$expID = rep("171019",1,length(levels(as.factor(tmpdf_loceff$region)))*2);
df_loceff = tmpdf_loceff;  
  
png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/localEff_2way_171019_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171019, x="region", y="LocEff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171207") + ylab("Local Efficiency")
dev.off()

###################################################################################################
df_171207 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171207_deg_wei")
regions_171207 <- character(nrow(df_171207));
regions_171207[c(1,2,3, 12,13,28, 18,23,24,25,26,27,31)] = "ventromedial thalamus"
regions_171207[c(4,5,6,7,11,14,15,16,29,30)] = "paraventricular hypothalamus"
regions_171207[c(8,9,10,17)] = "anterior hypothalamus"
regions_171207[c(19,20,21,22)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_171207),2);
condition <- c(rep("Ongoing",nrow(df_171207)), rep("LightON",nrow(df_171207)));
dat_171207 <- data.frame(c(df_171207$Ongoing, df_171207$LightON), condition, nodeid, rep(regions_171207,2));
colnames(dat_171207) <- c("NodeDegree", "condition", "nodeid", "region")
dat_171207$condition <- factor(dat_171207$condition, levels = c("Ongoing", "LightON"))
dat_171207$region <- factor(dat_171207$region)
dat_171207$NodeDegree <- as.numeric(sub(",", ".", c(df_171207$Ongoing, df_171207$LightON), fixed = TRUE))
model_171207 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_171207)
res_aov_171207 <- aov(model_171207)
summary(res_aov_171207)
cm=compare_means(NodeDegree ~ condition, data = dat_171207, group.by = "region")

tmpdf_degree = aggregate(NodeDegree ~ region + condition, data = dat_171207, mean);
tmpdf_degree$expID = rep("171207",1,length(levels(as.factor(tmpdf_degree$region)))*2);
df_degree = rbind(df_degree, tmpdf_degree);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/degree_2way_171207_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171207, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171207") + ylab("Node Strength")
dev.off()

df_171207 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171207_clus_wei")
nodeid <- rep(1:nrow(df_171207), 2);
condition <- c(rep("Ongoing",nrow(df_171207)), rep("LightON", nrow(df_171207)));
dat_171207 <- data.frame(c(df_171207$Ongoing, df_171207$LightON), condition, nodeid, rep(regions_171207,2));
colnames(dat_171207) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_171207$condition <- factor(dat_171207$condition, levels = c("Ongoing", "LightON"))
dat_171207$ClusCoeff <- as.numeric(sub(",", ".", c(df_171207$Ongoing, df_171207$LightON), fixed = TRUE))
model_171207 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_171207)
res_aov_171207 <- aov(model_171207)
summary(res_aov_171207)
cm=compare_means(ClusCoeff ~ condition, data = dat_171207, group.by = "region")

tmpdf_cluscoeff = aggregate(ClusCoeff ~ region + condition, data = dat_171207, mean);
tmpdf_cluscoeff$expID = rep("171207",1,length(levels(as.factor(tmpdf_cluscoeff$region)))*2);
df_cluscoeff = rbind(df_cluscoeff, tmpdf_cluscoeff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/cluscoeff_2way_171207_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171207, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171207") + ylab("Clustering Coefficient")
dev.off()


df_171207 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171207_localeff_wei")
nodeid <- rep(1:nrow(df_171207), 2);
condition <- c(rep("Ongoing",nrow(df_171207)), rep("LightON", nrow(df_171207)));
dat_171207 <- data.frame(c(df_171207$Ongoing, df_171207$LightON), condition, nodeid, rep(regions_171207,2));
colnames(dat_171207) <- c("LocEff", "condition", "nodeid", "region")
dat_171207$condition <- factor(dat_171207$condition, levels = c("Ongoing", "LightON"))
dat_171207$LocEff <- as.numeric(sub(",", ".", c(df_171207$Ongoing, df_171207$LightON), fixed = TRUE))
model_171207 <- lm(formula = "LocEff ~ condition * region + (1|nodeid)", data=dat_171207)
res_aov_171207 <- aov(model_171207)
summary(res_aov_171207)
cm=compare_means(LocEff ~ condition, data = dat_171207, group.by = "region")

tmpdf_loceff = aggregate(LocEff ~ region + condition, data = dat_171207, mean);
tmpdf_loceff$expID = rep("171207",1,length(levels(as.factor(tmpdf_loceff$region)))*2);
df_loceff = rbind(df_loceff, tmpdf_loceff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/localEff_2way_171207_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171207, x="region", y="LocEff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171207") + ylab("Local Efficiency")
dev.off()


###############################################################################
df_171208 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171208_deg_wei")
regions_171208 <- character(nrow(df_171208));
regions_171208[c(1,2,3,4,5,7,20,21,22,23,24,25,28)] = "ventromedial thalamus"
regions_171208[c(6, 8,9,11,13,16,17)] = "ventromedial hypothalamus"
regions_171208[c(10,12, 32)] = "paraventricular hypothalamus"
regions_171208[c(14,15,18,19)] = "dorsomedial hypothalamus"
regions_171208[c(26, 27, 29, 30, 31,33, 34)] = "anterior hypothalamus"
nodeid <- rep(1:nrow(df_171208),2);
condition <- c(rep("Ongoing",nrow(df_171208)), rep("LightON",nrow(df_171208)));
dat_171208 <- data.frame(c(df_171208$Ongoing, df_171208$LightON), condition, nodeid, rep(regions_171208,2));
colnames(dat_171208) <- c("NodeDegree", "condition", "nodeid", "region")
dat_171208$condition <- factor(dat_171208$condition, levels = c("Ongoing", "LightON"))
dat_171208$region <- factor(dat_171208$region)
dat_171208$NodeDegree <- as.numeric(sub(",", ".", c(df_171208$Ongoing, df_171208$LightON), fixed = TRUE))
model_171208 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_171208)
res_aov_171208 <- aov(model_171208)
summary(res_aov_171208)
cm=compare_means(NodeDegree ~ condition, data = dat_171208, group.by = "region")

tmpdf_degree = aggregate(NodeDegree ~ region + condition, data = dat_171208, mean);
tmpdf_degree$expID = rep("171208",1,length(levels(as.factor(tmpdf_degree$region)))*2);
df_degree = rbind(df_degree, tmpdf_degree);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/degree_2way_171208_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171208, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171208") + ylab("Node Strength")
dev.off()

df_171208 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171208_clus_wei")
nodeid <- rep(1:nrow(df_171208), 2);
condition <- c(rep("Ongoing",nrow(df_171208)), rep("LightON", nrow(df_171208)));
dat_171208 <- data.frame(c(df_171208$Ongoing, df_171208$LightON), condition, nodeid, rep(regions_171208,2));
colnames(dat_171208) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_171208$condition <- factor(dat_171208$condition, levels = c("Ongoing", "LightON"))
dat_171208$ClusCoeff <- as.numeric(sub(",", ".", c(df_171208$Ongoing, df_171208$LightON), fixed = TRUE))
model_171208 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_171208)
res_aov_171208 <- aov(model_171208)
summary(res_aov_171208)
cm=compare_means(ClusCoeff ~ condition, data = dat_171208, group.by = "region")

tmpdf_cluscoeff = aggregate(ClusCoeff ~ region + condition, data = dat_171208, mean);
tmpdf_cluscoeff$expID = rep("171208",1,length(levels(as.factor(tmpdf_cluscoeff$region)))*2);
df_cluscoeff = rbind(df_cluscoeff, tmpdf_cluscoeff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/cluscoeff_2way_171208_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171208, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171208") + ylab("Clustering Coefficient")
dev.off()


df_171208 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171208_localeff_wei")
nodeid <- rep(1:nrow(df_171208), 2);
condition <- c(rep("Ongoing",nrow(df_171208)), rep("LightON", nrow(df_171208)));
dat_171208 <- data.frame(c(df_171208$Ongoing, df_171208$LightON), condition, nodeid, rep(regions_171208,2));
colnames(dat_171208) <- c("LocEff", "condition", "nodeid", "region")
dat_171208$condition <- factor(dat_171208$condition, levels = c("Ongoing", "LightON"))
dat_171208$LocEff <- as.numeric(sub(",", ".", c(df_171208$Ongoing, df_171208$LightON), fixed = TRUE))
model_171208 <- lm(formula = "LocEff ~ condition * region + (1|nodeid)", data=dat_171208)
res_aov_171208 <- aov(model_171208)
summary(res_aov_171208)
cm=compare_means(LocEff ~ condition, data = dat_171208, group.by = "region")

tmpdf_loceff = aggregate(LocEff ~ region + condition, data = dat_171208, mean);
tmpdf_loceff$expID = rep("171208",1,length(levels(as.factor(tmpdf_loceff$region)))*2);
df_loceff = rbind(df_loceff, tmpdf_loceff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/localEff_2way_171208_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171208, x="region", y="LocEff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171213") + ylab("Local Efficiency")
dev.off()

##################################################################################################
df_171213 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171213_deg_wei")
regions_171213 <- character(nrow(df_171213));
regions_171213[c(11,12,13,20,21,22,23,24,25,26,27,29,30,31,32)] = "ventromedial thalamus"
regions_171213[c(1,2,3,4,5,6,7,8,9,10,14,15,16,17,18,19,28,33)] = "paraventricular hypothalamus"
nodeid <- rep(1:nrow(df_171213),2);
condition <- c(rep("Ongoing",nrow(df_171213)), rep("LightON",nrow(df_171213)));
dat_171213 <- data.frame(c(df_171213$Ongoing, df_171213$LightON), condition, nodeid, rep(regions_171213,2));
colnames(dat_171213) <- c("NodeDegree", "condition", "nodeid", "region")
dat_171213$condition <- factor(dat_171213$condition, levels = c("Ongoing", "LightON"))
dat_171213$region <- factor(dat_171213$region)
dat_171213$NodeDegree <- as.numeric(sub(",", ".", c(df_171213$Ongoing, df_171213$LightON), fixed = TRUE))
model_171213 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_171213)
res_aov_171213 <- aov(model_171213)
summary(res_aov_171213)
cm=compare_means(NodeDegree ~ condition, data = dat_171213, group.by = "region")

tmpdf_degree = aggregate(NodeDegree ~ region + condition, data = dat_171213, mean);
tmpdf_degree$expID = rep("171213",1,length(levels(as.factor(tmpdf_degree$region)))*2);
df_degree = rbind(df_degree, tmpdf_degree);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/degree_2way_171213_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171213, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171213") + ylab("Node Strength")
dev.off()

df_171213 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171213_clus_wei")
nodeid <- rep(1:nrow(df_171213), 2);
condition <- c(rep("Ongoing",nrow(df_171213)), rep("LightON", nrow(df_171213)));
dat_171213 <- data.frame(c(df_171213$Ongoing, df_171213$LightON), condition, nodeid, rep(regions_171213,2));
colnames(dat_171213) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_171213$condition <- factor(dat_171213$condition, levels = c("Ongoing", "LightON"))
dat_171213$ClusCoeff <- as.numeric(sub(",", ".", c(df_171213$Ongoing, df_171213$LightON), fixed = TRUE))
model_171213 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_171213)
res_aov_171213 <- aov(model_171213)
summary(res_aov_171213)
cm=compare_means(ClusCoeff ~ condition, data = dat_171213, group.by = "region")

tmpdf_cluscoeff = aggregate(ClusCoeff ~ region + condition, data = dat_171213, mean);
tmpdf_cluscoeff$expID = rep("171213",1,length(levels(as.factor(tmpdf_cluscoeff$region)))*2);
df_cluscoeff = rbind(df_cluscoeff, tmpdf_cluscoeff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/cluscoeff_2way_171213_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171213, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171213") + ylab("Clustering Coefficient")
dev.off()


df_171213 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171213_localeff_wei")
nodeid <- rep(1:nrow(df_171213), 2);
condition <- c(rep("Ongoing",nrow(df_171213)), rep("LightON", nrow(df_171213)));
dat_171213 <- data.frame(c(df_171213$Ongoing, df_171213$LightON), condition, nodeid, rep(regions_171213,2));
colnames(dat_171213) <- c("LocEff", "condition", "nodeid", "region")
dat_171213$condition <- factor(dat_171213$condition, levels = c("Ongoing", "LightON"))
dat_171213$LocEff <- as.numeric(sub(",", ".", c(df_171213$Ongoing, df_171213$LightON), fixed = TRUE))
model_171213 <- lm(formula = "LocEff ~ condition * region + (1|nodeid)", data=dat_171213)
res_aov_171213 <- aov(model_171213)
summary(res_aov_171213)
cm=compare_means(LocEff ~ condition, data = dat_171213, group.by = "region")

tmpdf_loceff = aggregate(LocEff ~ region + condition, data = dat_171213, mean);
tmpdf_loceff$expID = rep("171213",1,length(levels(as.factor(tmpdf_loceff$region)))*2);
df_loceff = rbind(df_loceff, tmpdf_loceff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/localEff_2way_171213_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171213, x="region", y="LocEff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171213") + ylab("Local Efficiency")
dev.off()


################################################################################################################################
df_180110 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180110_deg_wei")
regions_180110 <- character(nrow(df_180110));
regions_180110[c(1,2,11,18,20,21,23)] = "ventromedial thalamus"
regions_180110[c(3,4,5,6,7,9,10, 19,22,24,25,26,27,28,29,30,31,32,33,34,35)] = "paraventricular hypothalamus"
regions_180110[c(8)] = "anterior hypothalamus"
regions_180110[c(12,13,14,15,16,17)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_180110),2);
condition <- c(rep("Ongoing",nrow(df_180110)), rep("LightON",nrow(df_180110)));
dat_180110 <- data.frame(c(df_180110$Ongoing, df_180110$LightON), condition, nodeid, rep(regions_180110,2));
colnames(dat_180110) <- c("NodeDegree", "condition", "nodeid", "region")
dat_180110$condition <- factor(dat_180110$condition, levels = c("Ongoing", "LightON"))
dat_180110$region <- factor(dat_180110$region)
dat_180110$NodeDegree <- as.numeric(sub(",", ".", c(df_180110$Ongoing, df_180110$LightON), fixed = TRUE))
model_180110 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_180110)
res_aov_180110 <- aov(model_180110)
summary(res_aov_180110)
cm=compare_means(NodeDegree ~ condition, data = dat_180110, group.by = "region")

tmpdf_degree = aggregate(NodeDegree ~ region + condition, data = dat_180110, mean);
tmpdf_degree$expID = rep("180110",1,length(levels(as.factor(tmpdf_degree$region)))*2);
df_degree = rbind(df_degree, tmpdf_degree);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/degree_2way_180110_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180110, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180110") + ylab("Node Strength")
dev.off()

df_180110 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180110_clus_wei")
nodeid <- rep(1:nrow(df_180110), 2);
condition <- c(rep("Ongoing",nrow(df_180110)), rep("LightON", nrow(df_180110)));
dat_180110 <- data.frame(c(df_180110$Ongoing, df_180110$LightON), condition, nodeid, rep(regions_180110,2));
colnames(dat_180110) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_180110$condition <- factor(dat_180110$condition, levels = c("Ongoing", "LightON"))
dat_180110$ClusCoeff <- as.numeric(sub(",", ".", c(df_180110$Ongoing, df_180110$LightON), fixed = TRUE))
model_180110 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_180110)
res_aov_180110 <- aov(model_180110)
summary(res_aov_180110)
cm=compare_means(ClusCoeff ~ condition, data = dat_180110, group.by = "region")

tmpdf_cluscoeff = aggregate(ClusCoeff ~ region + condition, data = dat_180110, mean);
tmpdf_cluscoeff$expID = rep("180110",1,length(levels(as.factor(tmpdf_cluscoeff$region)))*2);
df_cluscoeff = rbind(df_cluscoeff, tmpdf_cluscoeff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/cluscoeff_2way_180110_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180110, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180110") + ylab("Clustering Coefficient")
dev.off()


df_180110 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180110_localeff_wei")
nodeid <- rep(1:nrow(df_180110), 2);
condition <- c(rep("Ongoing",nrow(df_180110)), rep("LightON", nrow(df_180110)));
dat_180110 <- data.frame(c(df_180110$Ongoing, df_180110$LightON), condition, nodeid, rep(regions_180110,2));
colnames(dat_180110) <- c("LocEff", "condition", "nodeid", "region")
dat_180110$condition <- factor(dat_180110$condition, levels = c("Ongoing", "LightON"))
dat_180110$LocEff <- as.numeric(sub(",", ".", c(df_180110$Ongoing, df_180110$LightON), fixed = TRUE))
model_180110 <- lm(formula = "LocEff ~ condition * region + (1|nodeid)", data=dat_180110)
res_aov_180110 <- aov(model_180110)
summary(res_aov_180110)
cm=compare_means(LocEff ~ condition, data = dat_180110, group.by = "region")

tmpdf_loceff = aggregate(LocEff ~ region + condition, data = dat_180110, mean);
tmpdf_loceff$expID = rep("180110",1,length(levels(as.factor(tmpdf_loceff$region)))*2);
df_loceff = rbind(df_loceff, tmpdf_loceff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/localEff_2way_180110_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180110, x="region", y="LocEff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180110") + ylab("Local Efficiency")
dev.off()

##############################################################################################
df_180111 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180111_deg_wei")
regions_180111 <- character(nrow(df_180111));
regions_180111[c(1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,22,23,24)] = "ventromedial thalamus"
regions_180111[c(14)] = "paraventricular hypothalamus"
regions_180111[c(19,20,21)] = "ventromedial hypothalamus"
nodeid <- rep(1:nrow(df_180111),2);
condition <- c(rep("Ongoing",nrow(df_180111)), rep("LightON",nrow(df_180111)));
dat_180111 <- data.frame(c(df_180111$Ongoing, df_180111$LightON), condition, nodeid, rep(regions_180111,2));
colnames(dat_180111) <- c("NodeDegree", "condition", "nodeid", "region")
dat_180111$condition <- factor(dat_180111$condition, levels = c("Ongoing", "LightON"))
dat_180111$region <- factor(dat_180111$region)
dat_180111$NodeDegree <- as.numeric(sub(",", ".", c(df_180111$Ongoing, df_180111$LightON), fixed = TRUE))
model_180111 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_180111)
res_aov_180111 <- aov(model_180111)
summary(res_aov_180111)
cm=compare_means(NodeDegree ~ condition, data = dat_180111, group.by = "region")

tmpdf_degree = aggregate(NodeDegree ~ region + condition, data = dat_180111, mean);
tmpdf_degree$expID = rep("180111",1,length(levels(as.factor(tmpdf_degree$region)))*2);
df_degree = rbind(df_degree, tmpdf_degree);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/degree_2way_180111_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180111, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180111") + ylab("Node Strength")
dev.off()

df_180111 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180111_clus_wei")
nodeid <- rep(1:nrow(df_180111), 2);
condition <- c(rep("Ongoing",nrow(df_180111)), rep("LightON", nrow(df_180111)));
dat_180111 <- data.frame(c(df_180111$Ongoing, df_180111$LightON), condition, nodeid, rep(regions_180111,2));
colnames(dat_180111) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_180111$condition <- factor(dat_180111$condition, levels = c("Ongoing", "LightON"))
dat_180111$ClusCoeff <- as.numeric(sub(",", ".", c(df_180111$Ongoing, df_180111$LightON), fixed = TRUE))
model_180111 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_180111)
res_aov_180111 <- aov(model_180111)
summary(res_aov_180111)
cm=compare_means(ClusCoeff ~ condition, data = dat_180111, group.by = "region")

tmpdf_cluscoeff = aggregate(ClusCoeff ~ region + condition, data = dat_180111, mean);
tmpdf_cluscoeff$expID = rep("180111",1,length(levels(as.factor(tmpdf_cluscoeff$region)))*2);
df_cluscoeff = rbind(df_cluscoeff, tmpdf_cluscoeff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/cluscoeff_2way_180111_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180111, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180111") + ylab("Clustering Coefficient")
dev.off()


df_180111 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180111_localeff_wei")
nodeid <- rep(1:nrow(df_180111), 2);
condition <- c(rep("Ongoing",nrow(df_180111)), rep("LightON", nrow(df_180111)));
dat_180111 <- data.frame(c(df_180111$Ongoing, df_180111$LightON), condition, nodeid, rep(regions_180111,2));
colnames(dat_180111) <- c("LocEff", "condition", "nodeid", "region")
dat_180111$condition <- factor(dat_180111$condition, levels = c("Ongoing", "LightON"))
dat_180111$LocEff <- as.numeric(sub(",", ".", c(df_180111$Ongoing, df_180111$LightON), fixed = TRUE))
model_180111 <- lm(formula = "LocEff ~ condition * region + (1|nodeid)", data=dat_180111)
res_aov_180111 <- aov(model_180111)
summary(res_aov_180111)
cm=compare_means(LocEff ~ condition, data = dat_180111, group.by = "region")

tmpdf_loceff = aggregate(LocEff ~ region + condition, data = dat_180111, mean);
tmpdf_loceff$expID = rep("180111",1,length(levels(as.factor(tmpdf_loceff$region)))*2);
df_loceff = rbind(df_loceff, tmpdf_loceff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/localEff_2way_180111_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180111, x="region", y="LocEff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180111") + ylab("Local Efficiency")
dev.off()

##############################################################################################################
df_180131 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180131_deg_wei")
regions_180131 <- character(nrow(df_180131));
regions_180131[c(1,2,4,5,6,14,15,16, 26)] = "ventromedial thalamus"
regions_180131[c(3,13,18,19,20,21,22,23,24, 27,29)] = "posterior hypothalamus"
regions_180131[c(7,8,9,10,11,12,25,28,30,31)] = "dorsomedial hypothalamus"
regions_180131[c(17)] = "paraventricular hypothalamus"
nodeid <- rep(1:nrow(df_180131),2);
condition <- c(rep("Ongoing",nrow(df_180131)), rep("LightON",nrow(df_180131)));
dat_180131 <- data.frame(c(df_180131$Ongoing, df_180131$LightON), condition, nodeid, rep(regions_180131,2));
colnames(dat_180131) <- c("NodeDegree", "condition", "nodeid", "region")
dat_180131$condition <- factor(dat_180131$condition, levels = c("Ongoing", "LightON"))
dat_180131$region <- factor(dat_180131$region)
dat_180131$NodeDegree <- as.numeric(sub(",", ".", c(df_180131$Ongoing, df_180131$LightON), fixed = TRUE))
model_180131 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_180131)
res_aov_180131 <- aov(model_180131)
summary(res_aov_180131)
cm=compare_means(NodeDegree ~ condition, data = dat_180131, group.by = "region")

tmpdf_degree = aggregate(NodeDegree ~ region + condition, data = dat_180131, mean);
tmpdf_degree$expID = rep("180131",1,length(levels(as.factor(tmpdf_degree$region)))*2);
df_degree = rbind(df_degree, tmpdf_degree);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/degree_2way_180131_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180131, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180131") + ylab("Node Strength")
dev.off()

df_180131 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180131_clus_wei")
nodeid <- rep(1:nrow(df_180131), 2);
condition <- c(rep("Ongoing",nrow(df_180131)), rep("LightON", nrow(df_180131)));
dat_180131 <- data.frame(c(df_180131$Ongoing, df_180131$LightON), condition, nodeid, rep(regions_180131,2));
colnames(dat_180131) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_180131$condition <- factor(dat_180131$condition, levels = c("Ongoing", "LightON"))
dat_180131$ClusCoeff <- as.numeric(sub(",", ".", c(df_180131$Ongoing, df_180131$LightON), fixed = TRUE))
model_180131 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_180131)
res_aov_180131 <- aov(model_180131)
summary(res_aov_180131)
cm=compare_means(ClusCoeff ~ condition, data = dat_180131, group.by = "region")

tmpdf_cluscoeff = aggregate(ClusCoeff ~ region + condition, data = dat_180131, mean);
tmpdf_cluscoeff$expID = rep("180131",1,length(levels(as.factor(tmpdf_cluscoeff$region)))*2);
df_cluscoeff = rbind(df_cluscoeff, tmpdf_cluscoeff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/cluscoeff_2way_180131_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180131, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180131") + ylab("Clustering Coefficient")
dev.off()


df_180131 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180131_localeff_wei")
nodeid <- rep(1:nrow(df_180131), 2);
condition <- c(rep("Ongoing",nrow(df_180131)), rep("LightON", nrow(df_180131)));
dat_180131 <- data.frame(c(df_180131$Ongoing, df_180131$LightON), condition, nodeid, rep(regions_180131,2));
colnames(dat_180131) <- c("LocEff", "condition", "nodeid", "region")
dat_180131$condition <- factor(dat_180131$condition, levels = c("Ongoing", "LightON"))
dat_180131$LocEff <- as.numeric(sub(",", ".", c(df_180131$Ongoing, df_180131$LightON), fixed = TRUE))
model_180131 <- lm(formula = "LocEff ~ condition * region + (1|nodeid)", data=dat_180131)
res_aov_180131 <- aov(model_180131)
summary(res_aov_180131)
cm=compare_means(LocEff ~ condition, data = dat_180131, group.by = "region")

tmpdf_loceff = aggregate(LocEff ~ region + condition, data = dat_180131, mean);
tmpdf_loceff$expID = rep("180131",1,length(levels(as.factor(tmpdf_loceff$region)))*2);
df_loceff = rbind(df_loceff, tmpdf_loceff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/localEff_2way_180131_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180131, x="region", y="LocEff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180131") + ylab("Local Efficiency")
dev.off()



########################################################################################################################
df_180221 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180221_deg_wei")
regions_180221 <- character(nrow(df_180221));
regions_180221[c(1, 35, 36, 48,49, 50, 2,3,4,12,17,18, 21, 22, 26, 31, 32, 33, 39, 40, 41, 42, 43, 44, 45, 46, 47)] = "ventromedial thalamus"
regions_180221[c(23, 24, 25)] = "ventromedial hypothalamus"
regions_180221[c(5,6,7,8,9,10,11,13,14,15,16, 27, 28, 29, 30, 37, 38)] = "dorsomedial hypothalamus"
regions_180221[c(19, 20, 34, 51, 52, 53, 54)] = "paraventricular hypothalamus"
nodeid <- rep(1:nrow(df_180221),2);
condition <- c(rep("Ongoing",nrow(df_180221)), rep("LightON",nrow(df_180221)));
dat_180221 <- data.frame(c(df_180221$Ongoing, df_180221$LightON), condition, nodeid, rep(regions_180221,2));
colnames(dat_180221) <- c("NodeDegree", "condition", "nodeid", "region")
dat_180221$condition <- factor(dat_180221$condition, levels = c("Ongoing", "LightON"))
dat_180221$region <- factor(dat_180221$region)
dat_180221$NodeDegree <- as.numeric(sub(",", ".", c(df_180221$Ongoing, df_180221$LightON), fixed = TRUE))
model_180221 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_180221)
res_aov_180221 <- aov(model_180221)
summary(res_aov_180221)
cm=compare_means(NodeDegree ~ condition, data = dat_180221, group.by = "region")

tmpdf_degree = aggregate(NodeDegree ~ region + condition, data = dat_180221, mean);
tmpdf_degree$expID = rep("180221",1,length(levels(as.factor(tmpdf_degree$region)))*2);
df_degree = rbind(df_degree, tmpdf_degree);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/degree_2way_180221_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180221, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180221") + ylab("Node Strength")
dev.off()

df_180221 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180221_clus_wei")
nodeid <- rep(1:nrow(df_180221), 2);
condition <- c(rep("Ongoing",nrow(df_180221)), rep("LightON", nrow(df_180221)));
dat_180221 <- data.frame(c(df_180221$Ongoing, df_180221$LightON), condition, nodeid, rep(regions_180221,2));
colnames(dat_180221) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_180221$condition <- factor(dat_180221$condition, levels = c("Ongoing", "LightON"))
dat_180221$ClusCoeff <- as.numeric(sub(",", ".", c(df_180221$Ongoing, df_180221$LightON), fixed = TRUE))
model_180221 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_180221)
res_aov_180221 <- aov(model_180221)
summary(res_aov_180221)
cm=compare_means(ClusCoeff ~ condition, data = dat_180221, group.by = "region")

tmpdf_cluscoeff = aggregate(ClusCoeff ~ region + condition, data = dat_180221, mean);
tmpdf_cluscoeff$expID = rep("180221",1,length(levels(as.factor(tmpdf_cluscoeff$region)))*2);
df_cluscoeff = rbind(df_cluscoeff, tmpdf_cluscoeff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/cluscoeff_2way_180221_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180221, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180221") + ylab("Clustering Coefficient")
dev.off()


df_180221 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180221_localeff_wei")
nodeid <- rep(1:nrow(df_180221), 2);
condition <- c(rep("Ongoing",nrow(df_180221)), rep("LightON", nrow(df_180221)));
dat_180221 <- data.frame(c(df_180221$Ongoing, df_180221$LightON), condition, nodeid, rep(regions_180221,2));
colnames(dat_180221) <- c("LocEff", "condition", "nodeid", "region")
dat_180221$condition <- factor(dat_180221$condition, levels = c("Ongoing", "LightON"))
dat_180221$LocEff <- as.numeric(sub(",", ".", c(df_180221$Ongoing, df_180221$LightON), fixed = TRUE))
model_180221 <- lm(formula = "LocEff ~ condition * region + (1|nodeid)", data=dat_180221)
res_aov_180221 <- aov(model_180221)
summary(res_aov_180221)
cm=compare_means(LocEff ~ condition, data = dat_180221, group.by = "region")

tmpdf_loceff = aggregate(LocEff ~ region + condition, data = dat_180221, mean);
tmpdf_loceff$expID = rep("180221",1,length(levels(as.factor(tmpdf_loceff$region)))*2);
df_loceff = rbind(df_loceff, tmpdf_loceff);


png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/localEff_2way_180221_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180221, x="region", y="LocEff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180221") + ylab("Local Efficiency")
dev.off()

######################################################################################################

df_180228 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180228_deg_wei")
regions_180228 <- character(nrow(df_180228));
regions_180228[c(1,4,5,6, 24, 25, 74)] = "ventromedial thalamus"
regions_180228[c(2, 11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23,26, 40, 42, 62, 63, 64)] = "posterior hypothalamus"
regions_180228[c(3,7,8,10, 19, 20, 27, 28, 29, 30, 31, 54, 56, 57, 58, 59, 60, 61, 65, 66,67,68,69, 73, 75, 76, 77,78,79)] = "dorsomedial hypothalamus"
regions_180228[c(9,46, 47, 48, 49, 50)] = "ventromedial hypothalamus"
regions_180228[c(32, 33, 34, 35, 36, 37, 38,39, 41)] = "arcuate hypothalamus"
regions_180228[c(43,44,45,51,52,53,55,70,71, 72)] = "zona incerta"
nodeid <- rep(1:nrow(df_180228),2);
condition <- c(rep("Ongoing",nrow(df_180228)), rep("LightON",nrow(df_180228)));
dat_180228 <- data.frame(c(df_180228$Ongoing, df_180228$LightON), condition, nodeid, rep(regions_180228,2));
colnames(dat_180228) <- c("NodeDegree", "condition", "nodeid", "region")
dat_180228$condition <- factor(dat_180228$condition, levels = c("Ongoing", "LightON"))
dat_180228$region <- factor(dat_180228$region)
dat_180228$NodeDegree <- as.numeric(sub(",", ".", c(df_180228$Ongoing, df_180228$LightON), fixed = TRUE))
model_180228 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_180228)
res_aov_180228 <- aov(model_180228)
summary(res_aov_180228)
cm=compare_means(NodeDegree ~ condition, data = dat_180228, group.by = "region")

tmpdf_degree = aggregate(NodeDegree ~ region + condition, data = dat_180228, mean);
tmpdf_degree$expID = rep("180228",1,length(levels(as.factor(tmpdf_degree$region)))*2);
df_degree = rbind(df_degree, tmpdf_degree);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/degree_2way_180228_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180228, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180228") + ylab("Node Strength")
dev.off()

df_180228 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180228_clus_wei")
nodeid <- rep(1:nrow(df_180228), 2);
condition <- c(rep("Ongoing",nrow(df_180228)), rep("LightON", nrow(df_180228)));
dat_180228 <- data.frame(c(df_180228$Ongoing, df_180228$LightON), condition, nodeid, rep(regions_180228,2));
colnames(dat_180228) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_180228$condition <- factor(dat_180228$condition, levels = c("Ongoing", "LightON"))
dat_180228$ClusCoeff <- as.numeric(sub(",", ".", c(df_180228$Ongoing, df_180228$LightON), fixed = TRUE))
model_180228 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_180228)
res_aov_180228 <- aov(model_180228)
summary(res_aov_180228)
cm=compare_means(ClusCoeff ~ condition, data = dat_180228, group.by = "region")

tmpdf_cluscoeff = aggregate(ClusCoeff ~ region + condition, data = dat_180228, mean);
tmpdf_cluscoeff$expID = rep("180228",1,length(levels(as.factor(tmpdf_cluscoeff$region)))*2);
df_cluscoeff = rbind(df_cluscoeff, tmpdf_cluscoeff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/cluscoeff_2way_180228_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180228, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180228") + ylab("Clustering Coefficient")
dev.off()


df_180228 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180228_localeff_wei")
nodeid <- rep(1:nrow(df_180228), 2);
condition <- c(rep("Ongoing",nrow(df_180228)), rep("LightON", nrow(df_180228)));
dat_180228 <- data.frame(c(df_180228$Ongoing, df_180228$LightON), condition, nodeid, rep(regions_180228,2));
colnames(dat_180228) <- c("LocEff", "condition", "nodeid", "region")
dat_180228$condition <- factor(dat_180228$condition, levels = c("Ongoing", "LightON"))
dat_180228$LocEff <- as.numeric(sub(",", ".", c(df_180228$Ongoing, df_180228$LightON), fixed = TRUE))
model_180228 <- lm(formula = "LocEff ~ condition * region + (1|nodeid)", data=dat_180228)
res_aov_180228 <- aov(model_180228)
summary(res_aov_180228)
cm=compare_means(LocEff ~ condition, data = dat_180228, group.by = "region")

tmpdf_loceff = aggregate(LocEff ~ region + condition, data = dat_180228, mean);
tmpdf_loceff$expID = rep("180228",1,length(levels(as.factor(tmpdf_loceff$region)))*2);
df_loceff = rbind(df_loceff, tmpdf_loceff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/localEff_2way_180228_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180228, x="region", y="LocEff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180228") + ylab("Local Efficiency")
dev.off()

###################################################################################################

df_180302 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180302_deg_wei")
regions_180302 <- character(nrow(df_180302));
#regions_180302[c(2,3,4,5,6, 13, 14, 15,30,36)] = "tuberomammillary nucleus"
regions_180302[c(2,3,4,5,6, 13, 14, 15,30,36, 9, 10, 41,42, 18, 19, 20, 21)] = "mammillary complex"
regions_180302[c(1,11,12,7, 8,16, 17, 23, 31, 37, 38,39,40, 47)] = "posterior hypothalamus"
#regions_180302[c(9, 10, 41,42)] = "premammillary nucleus"
#regions_180302[c(18, 19, 20, 21)] = "supramammillary nucleus"
regions_180302[c(26, 27, 28, 43,44,45,46)] = "arcuate hypothalamus"
regions_180302[c(22, 24, 25, 29)] = "paraventricular hypothalamus"
regions_180302[c(32,33,34,35)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_180302),2);
condition <- c(rep("Ongoing",nrow(df_180302)), rep("LightON",nrow(df_180302)));
dat_180302 <- data.frame(c(df_180302$Ongoing, df_180302$LightON), condition, nodeid, rep(regions_180302,2));
colnames(dat_180302) <- c("NodeDegree", "condition", "nodeid", "region")
dat_180302$condition <- factor(dat_180302$condition, levels = c("Ongoing", "LightON"))
dat_180302$region <- factor(dat_180302$region)
dat_180302$NodeDegree <- as.numeric(sub(",", ".", c(df_180302$Ongoing, df_180302$LightON), fixed = TRUE))
model_180302 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_180302)
res_aov_180302 <- aov(model_180302)
summary(res_aov_180302)
cm=compare_means(NodeDegree ~ condition, data = dat_180302, group.by = "region")

tmpdf_degree = aggregate(NodeDegree ~ region + condition, data = dat_180302, mean);
tmpdf_degree$expID = rep("180302",1,length(levels(as.factor(tmpdf_degree$region)))*2);
df_degree = rbind(df_degree, tmpdf_degree);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/degree_2way_180302_wei_v2.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180302, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180302") + ylab("Node Strength")
dev.off()

df_180302 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180302_clus_wei")
nodeid <- rep(1:nrow(df_180302), 2);
condition <- c(rep("Ongoing",nrow(df_180302)), rep("LightON", nrow(df_180302)));
dat_180302 <- data.frame(c(df_180302$Ongoing, df_180302$LightON), condition, nodeid, rep(regions_180302,2));
colnames(dat_180302) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_180302$condition <- factor(dat_180302$condition, levels = c("Ongoing", "LightON"))
dat_180302$ClusCoeff <- as.numeric(sub(",", ".", c(df_180302$Ongoing, df_180302$LightON), fixed = TRUE))
model_180302 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_180302)
res_aov_180302 <- aov(model_180302)
summary(res_aov_180302)
cm=compare_means(ClusCoeff ~ condition, data = dat_180302, group.by = "region")

tmpdf_cluscoeff = aggregate(ClusCoeff ~ region + condition, data = dat_180302, mean);
tmpdf_cluscoeff$expID = rep("180302",1,length(levels(as.factor(tmpdf_cluscoeff$region)))*2);
df_cluscoeff = rbind(df_cluscoeff, tmpdf_cluscoeff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/cluscoeff_2way_180302_wei_v2.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180302, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180302") + ylab("Clustering Coefficient")
dev.off()


df_180302 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180302_localeff_wei")
nodeid <- rep(1:nrow(df_180302), 2);
condition <- c(rep("Ongoing",nrow(df_180302)), rep("LightON", nrow(df_180302)));
dat_180302 <- data.frame(c(df_180302$Ongoing, df_180302$LightON), condition, nodeid, rep(regions_180302,2));
colnames(dat_180302) <- c("LocEff", "condition", "nodeid", "region")
dat_180302$condition <- factor(dat_180302$condition, levels = c("Ongoing", "LightON"))
dat_180302$LocEff <- as.numeric(sub(",", ".", c(df_180302$Ongoing, df_180302$LightON), fixed = TRUE))
model_180302 <- lm(formula = "LocEff ~ condition * region + (1|nodeid)", data=dat_180302)
res_aov_180302 <- aov(model_180302)
summary(res_aov_180302)
cm=compare_means(LocEff ~ condition, data = dat_180302, group.by = "region")

tmpdf_loceff = aggregate(LocEff ~ region + condition, data = dat_180302, mean);
tmpdf_loceff$expID = rep("180302",1,length(levels(as.factor(tmpdf_loceff$region)))*2);
df_loceff = rbind(df_loceff, tmpdf_loceff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/localEff_2way_180302_wei_v2.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180302, x="region", y="LocEff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180302") + ylab("Local Efficiency")
dev.off()

################################################################################################################################
df_180419 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180419_deg_wei")
regions_180419 <- character(nrow(df_180419));
regions_180419[c(1,2,11, 12,9, 10)] = "ventromedial hypothalamus"
regions_180419[c(4,5,6,7,8)] = "ventromedial thalamus"
regions_180419[c(36,37,38,39,51,52,53,54,55,56,57,58,59,60,61)] = "posterior hypothalamus"
regions_180419[c(3, 40, 41, 42, 43, 44,45,46,47,48,49,50, 13,14,15,16,17,18,23, 19,20,21,22,24,25,26,27,28,29,30,31,32,33,34,35)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_180419),2);
condition <- c(rep("Ongoing",nrow(df_180419)), rep("LightON",nrow(df_180419)));
dat_180419 <- data.frame(c(df_180419$Ongoing, df_180419$LightON), condition, nodeid, rep(regions_180419,2));
colnames(dat_180419) <- c("NodeDegree", "condition", "nodeid", "region")
dat_180419$condition <- factor(dat_180419$condition, levels = c("Ongoing", "LightON"))
dat_180419$region <- factor(dat_180419$region)
dat_180419$NodeDegree <- as.numeric(sub(",", ".", c(df_180419$Ongoing, df_180419$LightON), fixed = TRUE))
model_180419 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_180419)
res_aov_180419 <- aov(model_180419)
summary(res_aov_180419)
cm=compare_means(NodeDegree ~ condition, data = dat_180419, group.by = "region")

tmpdf_degree = aggregate(NodeDegree ~ region + condition, data = dat_180419, mean);
tmpdf_degree$expID = rep("180419",1,length(levels(as.factor(tmpdf_degree$region)))*2);
df_degree = rbind(df_degree, tmpdf_degree);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/degree_2way_180419_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180419, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180419") + ylab("Node Strength")
dev.off()

df_180419 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180419_clus_wei")
nodeid <- rep(1:nrow(df_180419), 2);
condition <- c(rep("Ongoing",nrow(df_180419)), rep("LightON", nrow(df_180419)));
dat_180419 <- data.frame(c(df_180419$Ongoing, df_180419$LightON), condition, nodeid, rep(regions_180419,2));
colnames(dat_180419) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_180419$condition <- factor(dat_180419$condition, levels = c("Ongoing", "LightON"))
dat_180419$ClusCoeff <- as.numeric(sub(",", ".", c(df_180419$Ongoing, df_180419$LightON), fixed = TRUE))
model_180419 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_180419)
res_aov_180419 <- aov(model_180419)
summary(res_aov_180419)
cm=compare_means(ClusCoeff ~ condition, data = dat_180419, group.by = "region")

tmpdf_cluscoeff = aggregate(ClusCoeff ~ region + condition, data = dat_180419, mean);
tmpdf_cluscoeff$expID = rep("180419",1,length(levels(as.factor(tmpdf_cluscoeff$region)))*2);
df_cluscoeff = rbind(df_cluscoeff, tmpdf_cluscoeff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/cluscoeff_2way_180419_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180419, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180419") + ylab("Clustering Coefficient")
dev.off()


df_180419 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180419_localeff_wei")
nodeid <- rep(1:nrow(df_180419), 2);
condition <- c(rep("Ongoing",nrow(df_180419)), rep("LightON", nrow(df_180419)));
dat_180419 <- data.frame(c(df_180419$Ongoing, df_180419$LightON), condition, nodeid, rep(regions_180419,2));
colnames(dat_180419) <- c("LocEff", "condition", "nodeid", "region")
dat_180419$condition <- factor(dat_180419$condition, levels = c("Ongoing", "LightON"))
dat_180419$LocEff <- as.numeric(sub(",", ".", c(df_180419$Ongoing, df_180419$LightON), fixed = TRUE))
model_180419 <- lm(formula = "LocEff ~ condition * region + (1|nodeid)", data=dat_180419)
res_aov_180419 <- aov(model_180419)
summary(res_aov_180419)
cm=compare_means(LocEff ~ condition, data = dat_180419, group.by = "region")

tmpdf_loceff = aggregate(LocEff ~ region + condition, data = dat_180419, mean);
tmpdf_loceff$expID = rep("180419",1,length(levels(as.factor(tmpdf_loceff$region)))*2);
df_loceff = rbind(df_loceff, tmpdf_loceff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/localEff_2way_180419_wei.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180419, x="region", y="LocEff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180419") + ylab("Local Efficiency")
dev.off()

########################################################################################################################
df_180420 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180420_deg_wei")
regions_180420 <- character(nrow(df_180420));
#regions_180420[c(51,52,53,54,55,56,57,58,59,60,61,62,63)] = "supramammillary nucleus"
regions_180420[c(51,52,53,54,55,56,57,58,59,60,61,62,63,19, 20, 21, 22, 23, 24, 25,26,27,28,29,30,31)] = "mammillary complex"
regions_180420[c(32, 33, 34, 35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50)] = "paraventricular hypothalamus"
#regions_180420[c(19, 20, 21, 22, 23, 24, 25,26,27,28,29,30,31)] = "premammillary nucleus"
regions_180420[c(2, 3, 4, 5, 6,7,8,9,10,11,12,13,14,15,16,17,18,64,65,66,67,68,69,70,71)] = "posterior hypothalamus"
regions_180420[c(1)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_180420),2);
condition <- c(rep("Ongoing",nrow(df_180420)), rep("LightON",nrow(df_180420)));
dat_180420 <- data.frame(c(df_180420$Ongoing, df_180420$LightON), condition, nodeid, rep(regions_180420,2));
colnames(dat_180420) <- c("NodeDegree", "condition", "nodeid", "region")
dat_180420$condition <- factor(dat_180420$condition, levels = c("Ongoing", "LightON"))
dat_180420$region <- factor(dat_180420$region)
dat_180420$NodeDegree <- as.numeric(sub(",", ".", c(df_180420$Ongoing, df_180420$LightON), fixed = TRUE))
model_180420 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_180420)
res_aov_180420 <- aov(model_180420)
summary(res_aov_180420)
cm=compare_means(NodeDegree ~ condition, data = dat_180420, group.by = "region")

tmpdf_degree = aggregate(NodeDegree ~ region + condition, data = dat_180420, mean);
tmpdf_degree$expID = rep("180420",1,length(levels(as.factor(tmpdf_degree$region)))*2);
df_degree = rbind(df_degree, tmpdf_degree);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/degree_2way_180420_wei_v2.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180420, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180420") + ylab("Node Strength")
dev.off()

df_180420 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180420_clus_wei")
nodeid <- rep(1:nrow(df_180420), 2);
condition <- c(rep("Ongoing",nrow(df_180420)), rep("LightON", nrow(df_180420)));
dat_180420 <- data.frame(c(df_180420$Ongoing, df_180420$LightON), condition, nodeid, rep(regions_180420,2));
colnames(dat_180420) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_180420$condition <- factor(dat_180420$condition, levels = c("Ongoing", "LightON"))
dat_180420$ClusCoeff <- as.numeric(sub(",", ".", c(df_180420$Ongoing, df_180420$LightON), fixed = TRUE))
model_180420 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_180420)
res_aov_180420 <- aov(model_180420)
summary(res_aov_180420)
cm=compare_means(ClusCoeff ~ condition, data = dat_180420, group.by = "region")

tmpdf_cluscoeff = aggregate(ClusCoeff ~ region + condition, data = dat_180420, mean);
tmpdf_cluscoeff$expID = rep("180420",1,length(levels(as.factor(tmpdf_cluscoeff$region)))*2);
df_cluscoeff = rbind(df_cluscoeff, tmpdf_cluscoeff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/cluscoeff_2way_180420_wei_v2.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180420, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180420") + ylab("Clustering Coefficient")
dev.off()


df_180420 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180420_localeff_wei")
nodeid <- rep(1:nrow(df_180420), 2);
condition <- c(rep("Ongoing",nrow(df_180420)), rep("LightON", nrow(df_180420)));
dat_180420 <- data.frame(c(df_180420$Ongoing, df_180420$LightON), condition, nodeid, rep(regions_180420,2));
colnames(dat_180420) <- c("LocEff", "condition", "nodeid", "region")
dat_180420$condition <- factor(dat_180420$condition, levels = c("Ongoing", "LightON"))
dat_180420$LocEff <- as.numeric(sub(",", ".", c(df_180420$Ongoing, df_180420$LightON), fixed = TRUE))
model_180420 <- lm(formula = "LocEff ~ condition * region + (1|nodeid)", data=dat_180420)
res_aov_180420 <- aov(model_180420)
summary(res_aov_180420)
cm=compare_means(LocEff ~ condition, data = dat_180420, group.by = "region")

tmpdf_loceff = aggregate(LocEff ~ region + condition, data = dat_180420, mean);
tmpdf_loceff$expID = rep("180420",1,length(levels(as.factor(tmpdf_loceff$region)))*2);
df_loceff = rbind(df_loceff, tmpdf_loceff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/localEff_2way_180420_wei_v2.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180420, x="region", y="LocEff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180420") + ylab("Local Efficiency")
dev.off()

#################################################################################################################

df_180423 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180423_deg_wei")
regions_180423 <- character(nrow(df_180423));
regions_180423[c(1,2,3,4)] = "ventromedial hypothalamus"
#regions_180423[c(5)] = "submammillothalamic nucleus"
regions_180423[c(6,7,8,9,10,18,19,20,26,27,28,29,30,31,32,33,34,35,36)] = "dorsomedial hypothalamus"
regions_180423[c(11, 12,13,14,15,16,17,21,22,23,24,25,51,52,53,54,55,56)] = "posterior hypothalamus"
#regions_180423[c(37, 38,39,40,41,42,43,44,45,46,47,48,49,50)] = "premammillary nucleus"
regions_180423[c(37, 38,39,40,41,42,43,44,45,46,47,48,49,50,5)] = "mammillary complex"
nodeid <- rep(1:nrow(df_180423),2);
condition <- c(rep("Ongoing",nrow(df_180423)), rep("LightON",nrow(df_180423)));
dat_180423 <- data.frame(c(df_180423$Ongoing, df_180423$LightON), condition, nodeid, rep(regions_180423,2));
colnames(dat_180423) <- c("NodeDegree", "condition", "nodeid", "region")
dat_180423$condition <- factor(dat_180423$condition, levels = c("Ongoing", "LightON"))
dat_180423$region <- factor(dat_180423$region)
dat_180423$NodeDegree <- as.numeric(sub(",", ".", c(df_180423$Ongoing, df_180423$LightON), fixed = TRUE))
model_180423 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_180423)
res_aov_180423 <- aov(model_180423)
summary(res_aov_180423)
cm=compare_means(NodeDegree ~ condition, data = dat_180423, group.by = "region")

tmpdf_degree = aggregate(NodeDegree ~ region + condition, data = dat_180423, mean);
tmpdf_degree$expID = rep("180423",1,length(levels(as.factor(tmpdf_degree$region)))*2);
df_degree = rbind(df_degree, tmpdf_degree);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/degree_2way_180423_wei_v2.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180423, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180423") + ylab("Node Strength")
dev.off()

df_180423 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180423_clus_wei")
nodeid <- rep(1:nrow(df_180423), 2);
condition <- c(rep("Ongoing",nrow(df_180423)), rep("LightON", nrow(df_180423)));
dat_180423 <- data.frame(c(df_180423$Ongoing, df_180423$LightON), condition, nodeid, rep(regions_180423,2));
colnames(dat_180423) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_180423$condition <- factor(dat_180423$condition, levels = c("Ongoing", "LightON"))
dat_180423$ClusCoeff <- as.numeric(sub(",", ".", c(df_180423$Ongoing, df_180423$LightON), fixed = TRUE))
model_180423 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_180423)
res_aov_180423 <- aov(model_180423)
summary(res_aov_180423)
cm=compare_means(ClusCoeff ~ condition, data = dat_180423, group.by = "region")

tmpdf_cluscoeff = aggregate(ClusCoeff ~ region + condition, data = dat_180423, mean);
tmpdf_cluscoeff$expID = rep("180423",1,length(levels(as.factor(tmpdf_cluscoeff$region)))*2);
df_cluscoeff = rbind(df_cluscoeff, tmpdf_cluscoeff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/cluscoeff_2way_180423_wei_v2.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180423, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180423") + ylab("Clustering Coefficient")
dev.off()


df_180423 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180423_localeff_wei")
nodeid <- rep(1:nrow(df_180423), 2);
condition <- c(rep("Ongoing",nrow(df_180423)), rep("LightON", nrow(df_180423)));
dat_180423 <- data.frame(c(df_180423$Ongoing, df_180423$LightON), condition, nodeid, rep(regions_180423,2));
colnames(dat_180423) <- c("LocEff", "condition", "nodeid", "region")
dat_180423$condition <- factor(dat_180423$condition, levels = c("Ongoing", "LightON"))
dat_180423$LocEff <- as.numeric(sub(",", ".", c(df_180423$Ongoing, df_180423$LightON), fixed = TRUE))
model_180423 <- lm(formula = "LocEff ~ condition * region + (1|nodeid)", data=dat_180423)
res_aov_180423 <- aov(model_180423)
summary(res_aov_180423)
cm=compare_means(LocEff ~ condition, data = dat_180423, group.by = "region")

tmpdf_loceff = aggregate(LocEff ~ region + condition, data = dat_180423, mean);
tmpdf_loceff$expID = rep("180423",1,length(levels(as.factor(tmpdf_loceff$region)))*2);
df_loceff = rbind(df_loceff, tmpdf_loceff);

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/localEff_2way_180423_wei_v2.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180423, x="region", y="LocEff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180423") + ylab("Local Efficiency")
dev.off()

#################################################################################################################
df_degree$expID <- as.factor(df_degree$expID)
df_degree$Day = rep(0,nrow(df_degree))
df_degree$Day[which(df_degree$expID == "171019")] = "Day";
df_degree$Day[which(df_degree$expID == "171207")] = "Night";
df_degree$Day[which(df_degree$expID == "171208")] = "Night";
df_degree$Day[which(df_degree$expID == "171213")] = "Night";
df_degree$Day[which(df_degree$expID == "180110")] = "Night";
df_degree$Day[which(df_degree$expID == "180111")] = "Night";
df_degree$Day[which(df_degree$expID == "180131")] = "Day";
df_degree$Day[which(df_degree$expID == "180221")] = "Night";
df_degree$Day[which(df_degree$expID == "180228")] = "Night";
df_degree$Day[which(df_degree$expID == "180302")] = "Day";
df_degree$Day[which(df_degree$expID == "180419")] = "Day";
df_degree$Day[which(df_degree$expID == "180420")] = "Day";
df_degree$Day[which(df_degree$expID == "180423")] = "Day";
model_degree = lm(data = df_degree, formula = "NodeDegree ~ condition + region + Day" )
res_aov_degree = aov(model_degree)
summary(res_aov_degree)
emm1 = emmeans(model_degree, ~ condition + Day | region, adjust = 'sidak')

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/strength_3way.png", width = 1240, height = 2840, units = "px")
plot(emm1, xlab = "Node Strength", ylab="") + theme(text = element_text(size = 34))
dev.off()

write.csv2(df_degree[which(df_degree$region == "paraventricular hypothalamus"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/pv_ht_3way.csv")
write.csv2(df_degree[which(df_degree$region == "ventromedial hypothalamus"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/vm_ht_3way.csv")
write.csv2(df_degree[which(df_degree$region == "posterior hypothalamus"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/p_ht_3way.csv")
write.csv2(df_degree[which(df_degree$region == "arcuate hypothalamus"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/arc_ht_3way.csv")
write.csv2(df_degree[which(df_degree$region == "zona incerta"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/zi_3way.csv")
write.csv2(df_degree[which(df_degree$region == "mammillary complex"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/mm_complex_3way.csv")
write.csv2(df_degree[which(df_degree$region == "anterior hypothalamus"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/a_ht_3way.csv")
write.csv2(df_degree[which(df_degree$region == "dorsomedial hypothalamus"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/dm_ht_3way.csv")
write.csv2(df_degree[which(df_degree$region == "ventromedial thalamus"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/vm_t_3way.csv")

con1 <- pairs(emm1)
cm = compare_means(NodeDegree ~ condition, data = df_degree, group.by = "region")
cmd = compare_means(NodeDegree ~ condition, data = df_degree, group.by = "Day")
cmd2 = compare_means(NodeDegree ~ condition, data = df_degree, group.by = c("region", "Day"), method = 'anova')
ggboxplot(df_degree, y="NodeDegree", x = "region", fill = "region", facet.by = "condition") + 
  stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
  theme(text = element_text(size = 28)) 

df_cluscoeff$expID <- as.factor(df_cluscoeff$expID)
df_cluscoeff$Day = df_degree$Day
model_cluscoeff = lm(data = df_cluscoeff, formula = "ClusCoeff ~ condition + region + Day" )
res_aov_cluscoeff = aov(model_cluscoeff)
summary(res_aov_cluscoeff)
cm = compare_means(ClusCoeff ~ condition, data = df_cluscoeff, group.by = "region")
cmd = compare_means(ClusCoeff ~ condition, data = df_cluscoeff, group.by = "Day")

write.csv2(df_cluscoeff[which(df_cluscoeff$region == "paraventricular hypothalamus"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/pv_ht_3way_cc.csv")
write.csv2(df_cluscoeff[which(df_cluscoeff$region == "ventromedial hypothalamus"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/vm_ht_3way_cc.csv")
write.csv2(df_cluscoeff[which(df_cluscoeff$region == "posterior hypothalamus"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/p_ht_3way_cc.csv")
write.csv2(df_cluscoeff[which(df_cluscoeff$region == "arcuate hypothalamus"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/arc_ht_3way_cc.csv")
write.csv2(df_cluscoeff[which(df_cluscoeff$region == "zona incerta"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/zi_3way_cc.csv")
write.csv2(df_cluscoeff[which(df_cluscoeff$region == "mammillary complex"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/mm_complex_3way_cc.csv")
write.csv2(df_cluscoeff[which(df_cluscoeff$region == "anterior hypothalamus"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/a_ht_3way_cc.csv")
write.csv2(df_cluscoeff[which(df_cluscoeff$region == "dorsomedial hypothalamus"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/dm_ht_3way_cc.csv")
write.csv2(df_cluscoeff[which(df_cluscoeff$region == "ventromedial thalamus"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/vm_t_3way_cc.csv")


df_loceff$expID <- as.factor(df_loceff$expID)
df_loceff$Day = df_cluscoeff$Day
model_loceff = lm(data = df_loceff, formula = "LocEff ~ condition + region + Day" )
res_aov_loceff = aov(model_loceff)
summary(res_aov_loceff)
cm = compare_means(LocEff ~ condition, data = df_loceff, group.by = "region")
cmd = compare_means(LocEff ~ condition, data = df_loceff, group.by = "Day")

write.csv2(df_loceff[which(df_loceff$region == "paraventricular hypothalamus"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/pv_ht_3way_le.csv")
write.csv2(df_loceff[which(df_loceff$region == "ventromedial hypothalamus"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/vm_ht_3way_le.csv")
write.csv2(df_loceff[which(df_loceff$region == "posterior hypothalamus"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/p_ht_3way_le.csv")
write.csv2(df_loceff[which(df_loceff$region == "arcuate hypothalamus"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/arc_ht_3way_le.csv")
write.csv2(df_loceff[which(df_loceff$region == "zona incerta"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/zi_3way_le.csv")
write.csv2(df_loceff[which(df_loceff$region == "mammillary complex"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/mm_complex_3way_le.csv")
write.csv2(df_loceff[which(df_loceff$region == "anterior hypothalamus"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/a_ht_3way_le.csv")
write.csv2(df_loceff[which(df_loceff$region == "dorsomedial hypothalamus"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/dm_ht_3way_le.csv")
write.csv2(df_loceff[which(df_loceff$region == "ventromedial thalamus"),], "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/vm_t_3way_le.csv")

##################################################################################################################


clus_coefs_o = c()
clus_coefs_e = c()
strengths_o = c()
strengths_e = c()
localeffs_o = c()
localeffs_e = c()
i = 1;
for(day in days){
  adj_glmcc <- read.csv(paste("~/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_", day, "_pycharm.csv", sep = "") , header=FALSE)
  tmp_o <- read.csv(paste("~/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/", day, "_clus_ongoing_wei.csv", sep = ""), header = FALSE)
  tmp_e <- read.csv(paste("~/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/", day, "_clus_evoked_wei.csv", sep = ""), header = FALSE)
  clus_coef <- data.frame("Clustering" = c(t(tmp_o), t(tmp_e)), "Condition" = c(rep("Ongoing",nrow(tmp_o)), 
                                                                                  rep("LightON",nrow(tmp_e))))
  clus_coefs_o[i] = mean(clus_coef$Clustering[clus_coef$Condition == "Ongoing" ])  
  clus_coefs_e[i] = mean(clus_coef$Clustering[clus_coef$Condition == "LightON" ])
  #write_pzfx(data.frame("Ongoing" = tmp_o, "lightON" = tmp_e),  "~/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/tmp.pzfx")
  #lf <- lm(formula = "Clustering ~ Condition", data=clus_coef)
  #res_aov <- aov(lf)
  #summary(res_aov)
  #sss = wilcox.test(formula = Clustering ~ Condition, data = clus_coef, paired = TRUE)
  #pvalues_cc = c(sss$p.value, pvalues_cc)
  #print(paste("~/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/", day, "_clus_wei_comparison.png", sep = ""))
  # png(filename = paste("~/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/", day, "_clus_wei_comparison.png", sep = ""),
  #     width = 680, height = 1020, units = "px")
  # ggboxplot(data = clus_coef, x = "Condition", y = "Clustering", color = "Condition", fill = "Condition") + 
  #   stat_compare_means() + xlab("") + ylab("Clustering Coefficient")  + xlab("") + theme(legend.title = element_blank()) + 
  #   theme(legend.position = "none") + theme(text = element_text(size = 20))
  # dev.off()
  
  tmp_o <- read.csv(paste("~/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/", day, "_deg_ongoing_wei.csv", sep = ""), header = FALSE)
  tmp_e <- read.csv(paste("~/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/", day, "_deg_evoked_wei.csv", sep = ""), header = FALSE)
  strength <- data.frame("Strength" = c(t(tmp_o), t(tmp_e)), "Condition" = c(rep("Ongoing",nrow(tmp_o)), 
                                                                             rep("LightON",nrow(tmp_e))))
  strengths_o[i] = mean(strength$Strength[strength$Condition == "Ongoing" ])  
  strengths_e[i] = mean(strength$Strength[strength$Condition == "LightON" ])
  # png(filename = paste("~/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/", day, "_strength_wei_comparison.png", sep = ""),
  #     width = 680, height = 1020, units = "px")
  # ggboxplot(data = strength, x = "Condition", y = "Strength", color = "Condition", fill = "Condition") + 
  #   stat_compare_means() + xlab("") + ylab("Strength")  + xlab("") + theme(legend.title = element_blank()) + 
  #   theme(legend.position = "none") + theme(text = element_text(size = 20))
  # dev.off()
  
  tmp_o <- read.csv(paste("~/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/", day, "_localeff_ongoing_wei.csv", sep = ""), header = FALSE)
  tmp_e <- read.csv(paste("~/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/", day, "_localeff_evoked_wei.csv", sep = ""), header = FALSE)
  localeff <- data.frame("LocalEff" = c(t(tmp_o), t(tmp_e)), "Condition" = c(rep("Ongoing",nrow(tmp_o)), 
                                                                             rep("LightON",nrow(tmp_e))))
  localeffs_o[i] = mean(localeff$LocalEff[localeff$Condition == "Ongoing" ])  
  localeffs_e[i] = mean(localeff$LocalEff[localeff$Condition == "LightON" ])
  # png(filename = paste("~/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/", day, "_localeff_wei_comparison.png", sep = ""),
  #     width = 680, height = 1020, units = "px")
  # ggboxplot(data = localeff, x = "Condition", y = "LocalEff", color = "Condition", fill = "Condition") + 
  #   stat_compare_means() + xlab("") + ylab("Local Efficiency")  + xlab("") + theme(legend.title = element_blank()) + 
  #   theme(legend.position = "none") + theme(text = element_text(size = 20))
  # dev.off()
  i = i +1;
}

tt = c()
for (i in 1:dim(adj_glmcc)[1]){
  tt = c(tt, toString(i))
}


day = days[1]
adj_glmcc <- read.csv(paste("~/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_", day, "_pycharm.csv", sep = "") , header=FALSE)
adj_glmcc <- as.matrix(adj_glmcc)
net <- network(adj_glmcc>quantile(adj_glmcc,0.97))
plot(net, label = as.character(tt), vertex.cex=1.95, label.pos=5)


#####Second levels
cc_df = data.frame("Clustering")




####################################################################################
######################## 171019 #########################
adj_pearson_171019 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_171019_pycharm.csv", header=FALSE)
adj_pearson_171019 <- as.matrix(adj_pearson_171019);
colnames(adj_pearson_171019) <- regions_171019;
rownames(adj_pearson_171019) <- regions_171019;
plot(adj_pearson_171019)
net <- network(adj_pearson_171019, loop=TRUE)
#plot(net)

region_names_171019 <- unique(regions_171019);
nregions_171019 = length(region_names_171019);
conn_171019 = matrix(0, nregions_171019, nregions_171019);
colnames(conn_171019) <- rownames(conn_171019) <- region_names_171019;
cx = 1;
for (i in region_names_171019){
  idx = which(regions_171019 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_171019){
    idy = which(regions_171019 == j, arr.ind = FALSE);
    conn_171019[cx,cy] = sum(adj_pearson_171019[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_171019, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_171019.csv");
sel_mat <- which(conn_171019 > quantile(conn_171019, .75), arr.ind = TRUE);
net_171019 <- network(conn_171019, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_171019.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_171019, label = unique(rownames(conn_171019)), edge.label = conn_171019, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171019_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171019_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171019_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171019_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_171019)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_171019, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_ongoing_171019.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1:2], track.index = c(2,2), text = "PV HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[3:126], track.index = c(2,2), text = "VM T", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))


for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_171019);
idd = which(adj_pearson_171019 > quantile(adj_pearson_171019, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_171019)[1]){
  for(jj in 1:dim(adj_pearson_171019)[2]){
    if (adj_pearson_171019[ii,jj]> quantile(adj_pearson_171019, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_171019[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()



################################### EVOKED ############################################################
adj_pearson_171019 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_171019_pycharm.csv", header=FALSE)
adj_pearson_171019 <- as.matrix(adj_pearson_171019);
colnames(adj_pearson_171019) <- regions_171019;
rownames(adj_pearson_171019) <- regions_171019;
plot(adj_pearson_171019)
net <- network(adj_pearson_171019, loop=TRUE)
#plot(net)

region_names_171019 <- unique(regions_171019);
nregions_171019 = length(region_names_171019);
conn_171019 = matrix(0, nregions_171019, nregions_171019);
colnames(conn_171019) <- rownames(conn_171019) <- region_names_171019;
cx = 1;
for (i in region_names_171019){
  idx = which(regions_171019 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_171019){
    idy = which(regions_171019 == j, arr.ind = FALSE);
    conn_171019[cx,cy] = sum(adj_pearson_171019[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_171019, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_171019.csv");
sel_mat <- which(conn_171019 > quantile(conn_171019, .7), arr.ind = TRUE);
net_171019 <- network(conn_171019, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_171019.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_171019, label = unique(rownames(conn_171019)), edge.label = conn_171019, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171019_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171019_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171019_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171019_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_171019)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_171019, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_evoked_171019.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1:2], track.index = c(2,2), text = "PV HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[3:126], track.index = c(2,2), text = "VM T", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))



for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_171019);
idd = which(adj_pearson_171019 > quantile(adj_pearson_171019, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_171019)[1]){
  for(jj in 1:dim(adj_pearson_171019)[2]){
    if (adj_pearson_171019[ii,jj]> quantile(adj_pearson_171019, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_171019[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()


####################################################################################
######################## 171207 #########################
adj_pearson_171207 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_171207_pycharm.csv", header=FALSE)
adj_pearson_171207 <- as.matrix(adj_pearson_171207);
colnames(adj_pearson_171207) <- regions_171207;
rownames(adj_pearson_171207) <- regions_171207;
plot(adj_pearson_171207)

idx = which(abs(adj_pearson_171207) < 0.7, arr.ind = TRUE);
for (idd in  idx){ 
  if(adj_pearson_171207[idd] > 0) 
    { 
    adj_pearson_171207[idd] = 0
  }
}
#adj_pearson_171207[idx] = 0;
colnames(adj_pearson_171207) <- seq(1,31)
net <- network(adj_pearson_171207, loop = TRUE)
col_mat = matrix(0,31,31)
col_mat_idx_neg = which(adj_pearson_171207 < 0, arr.ind = TRUE);
col_mat_idx_pos = which(adj_pearson_171207 > 0, arr.ind = TRUE);
col_mat[col_mat_idx_neg] = "blue";
col_mat[col_mat_idx_pos] = "red";
adj_pearson_171207[col_mat_idx_neg] = adj_pearson_171207[col_mat_idx_neg] * 15;
png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/circle_plot_glmcc_toy.png", 
    width=1200, height=1200, units = "px")
plot.network(net, mode = "circle", displaylabels = TRUE, edge.lwd = abs(adj_pearson_171207), vertex.cex = 3, label.cex = 2, 
             label.pos = 5, edge.col = array(col_mat), vertex.col = "white")
dev.off()

circos.clear()
tt = c()
for (i in 1:dim(adj_pearson_171207)[1]){
  tt = c(tt, toString(i))
}
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
  niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

for(ii in 1:dim(adj_pearson_171207)[1]){
  for(jj in 1:dim(adj_pearson_171207)[2]){
    if (adj_pearson_171207[ii,jj]> quantile(adj_pearson_171207, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 1, lwd=0.1 )
    }
  }
}







region_names_171207 <- unique(regions_171207);
nregions_171207 = length(region_names_171207);
conn_171207 = matrix(0, nregions_171207, nregions_171207);
colnames(conn_171207) <- rownames(conn_171207) <- region_names_171207;
cx = 1;
for (i in region_names_171207){
  idx = which(regions_171207 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_171207){
    idy = which(regions_171207 == j, arr.ind = FALSE);
    conn_171207[cx,cy] = sum(adj_pearson_171207[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_171207, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_171207.csv");
sel_mat <- which(conn_171207 > quantile(conn_171207, .75), arr.ind = TRUE);
net_171207 <- network(conn_171207, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_171207.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_171207, label = unique(rownames(conn_171207)), edge.label = conn_171207, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171207_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171207_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171207_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171207_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_171207)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_171207, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_ongoing_171207.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1:4], track.index = c(2,2), text = "A HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[5:8], track.index = c(2,2), text = "DM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[9:18], track.index = c(2,2), text = "PV HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[19:31], track.index = c(2,2), text = "VM T", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_171207);
idd = which(adj_pearson_171207 > quantile(adj_pearson_171207, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_171207)[1]){
  for(jj in 1:dim(adj_pearson_171207)[2]){
    if (adj_pearson_171207[ii,jj]> quantile(adj_pearson_171207, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_171207[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()



################################### EVOKED ############################################################
adj_pearson_171207 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_171207_pycharm.csv", header=FALSE)
adj_pearson_171207 <- as.matrix(adj_pearson_171207);
colnames(adj_pearson_171207) <- regions_171207;
rownames(adj_pearson_171207) <- regions_171207;
plot(adj_pearson_171207)
net <- network(adj_pearson_171207, loop=TRUE)
#plot(net)

region_names_171207 <- unique(regions_171207);
nregions_171207 = length(region_names_171207);
conn_171207 = matrix(0, nregions_171207, nregions_171207);
colnames(conn_171207) <- rownames(conn_171207) <- region_names_171207;
cx = 1;
for (i in region_names_171207){
  idx = which(regions_171207 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_171207){
    idy = which(regions_171207 == j, arr.ind = FALSE);
    conn_171207[cx,cy] = sum(adj_pearson_171207[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_171207, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_171207.csv");
sel_mat <- which(conn_171207 > quantile(conn_171207, .75), arr.ind = TRUE);
net_171207 <- network(conn_171207, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_171207.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_171207, label = unique(rownames(conn_171207)), edge.label = conn_171207, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171207_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171207_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171207_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171207_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_171207)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_171207, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_evoked_171207.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1:4], track.index = c(2,2), text = "A HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[5:8], track.index = c(2,2), text = "DM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[9:18], track.index = c(2,2), text = "PV HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[19:31], track.index = c(2,2), text = "VM T", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_171207);
idd = which(adj_pearson_171207 > quantile(adj_pearson_171207, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_171207)[1]){
  for(jj in 1:dim(adj_pearson_171207)[2]){
    if (adj_pearson_171207[ii,jj]> quantile(adj_pearson_171207, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_171207[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()

####################################################################################
######################## 171208 #########################
adj_pearson_171208 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_171208_pycharm.csv", header=FALSE)
adj_pearson_171208 <- as.matrix(adj_pearson_171208);
colnames(adj_pearson_171208) <- regions_171208;
rownames(adj_pearson_171208) <- regions_171208;
plot(adj_pearson_171208)
net <- network(adj_pearson_171208, loop=TRUE)
#plot(net)

region_names_171208 <- unique(regions_171208);
nregions_171208 = length(region_names_171208);
conn_171208 = matrix(0, nregions_171208, nregions_171208);
colnames(conn_171208) <- rownames(conn_171208) <- region_names_171208;
cx = 1;
for (i in region_names_171208){
  idx = which(regions_171208 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_171208){
    idy = which(regions_171208 == j, arr.ind = FALSE);
    conn_171208[cx,cy] = sum(adj_pearson_171208[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_171208, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_171208.csv");
sel_mat <- which(conn_171208 > quantile(conn_171208, .75), arr.ind = TRUE);
net_171208 <- network(conn_171208, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_171208.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_171208, label = unique(rownames(conn_171208)), edge.label = conn_171208, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171208_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171208_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171208_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171208_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_171208)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_171208, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_ongoing_171208.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1:7], track.index = c(2,2), text = "A HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[8:11], track.index = c(2,2), text = "DM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[12:14], track.index = c(2,2), text = "PV HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[15:21], track.index = c(2,2), text = "VM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[22:34], track.index = c(2,2), text = "VM T", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_171208);
idd = which(adj_pearson_171208 > quantile(adj_pearson_171208, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_171208)[1]){
  for(jj in 1:dim(adj_pearson_171208)[2]){
    if (adj_pearson_171208[ii,jj]> quantile(adj_pearson_171208, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_171208[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()



################################### EVOKED ############################################################
adj_pearson_171208 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_171208_pycharm.csv", header=FALSE)
adj_pearson_171208 <- as.matrix(adj_pearson_171208);
colnames(adj_pearson_171208) <- regions_171208;
rownames(adj_pearson_171208) <- regions_171208;
plot(adj_pearson_171208)
net <- network(adj_pearson_171208, loop=TRUE)
#plot(net)

region_names_171208 <- unique(regions_171208);
nregions_171208 = length(region_names_171208);
conn_171208 = matrix(0, nregions_171208, nregions_171208);
colnames(conn_171208) <- rownames(conn_171208) <- region_names_171208;
cx = 1;
for (i in region_names_171208){
  idx = which(regions_171208 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_171208){
    idy = which(regions_171208 == j, arr.ind = FALSE);
    conn_171208[cx,cy] = sum(adj_pearson_171208[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_171208, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_171208.csv");
sel_mat <- which(conn_171208 > quantile(conn_171208, .75), arr.ind = TRUE);
net_171208 <- network(conn_171208, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_171208.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_171208, label = unique(rownames(conn_171208)), edge.label = conn_171208, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171208_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171208_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171208_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171208_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_171208)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_171208, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_evoked_171208.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1:7], track.index = c(2,2), text = "A HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[8:11], track.index = c(2,2), text = "DM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[12:14], track.index = c(2,2), text = "PV HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[15:21], track.index = c(2,2), text = "VM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[22:34], track.index = c(2,2), text = "VM T", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_171208);
idd = which(adj_pearson_171208 > quantile(adj_pearson_171208, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_171208)[1]){
  for(jj in 1:dim(adj_pearson_171208)[2]){
    if (adj_pearson_171208[ii,jj]> quantile(adj_pearson_171208, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_171208[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()





####################################################################################
######################## 171213 #########################
adj_pearson_171213 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_171213_pycharm.csv", header=FALSE)
adj_pearson_171213 <- as.matrix(adj_pearson_171213);
colnames(adj_pearson_171213) <- regions_171213;
rownames(adj_pearson_171213) <- regions_171213;
plot(adj_pearson_171213)
net <- network(adj_pearson_171213, loop=TRUE)
#plot(net)

region_names_171213 <- unique(regions_171213);
nregions_171213 = length(region_names_171213);
conn_171213 = matrix(0, nregions_171213, nregions_171213);
colnames(conn_171213) <- rownames(conn_171213) <- region_names_171213;
cx = 1;
for (i in region_names_171213){
  idx = which(regions_171213 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_171213){
    idy = which(regions_171213 == j, arr.ind = FALSE);
    conn_171213[cx,cy] = sum(adj_pearson_171213[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_171213, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_171213.csv");
sel_mat <- which(conn_171213 > quantile(conn_171213, .75), arr.ind = TRUE);
net_171213 <- network(conn_171213, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_171213.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_171213, label = unique(rownames(conn_171213)), edge.label = conn_171213, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171213_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171213_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171213_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171213_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_171213)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_171213, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_ongoing_171213.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1:18], track.index = c(2,2), text = "PV HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[19:33], track.index = c(2,2), text = "VM T", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_171213);
idd = which(adj_pearson_171213 > quantile(adj_pearson_171213, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_171213)[1]){
  for(jj in 1:dim(adj_pearson_171213)[2]){
    if (adj_pearson_171213[ii,jj]> quantile(adj_pearson_171213, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_171213[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()



################################### EVOKED ############################################################
adj_pearson_171213 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_171213_pycharm.csv", header=FALSE)
adj_pearson_171213 <- as.matrix(adj_pearson_171213);
colnames(adj_pearson_171213) <- regions_171213;
rownames(adj_pearson_171213) <- regions_171213;
plot(adj_pearson_171213)
net <- network(adj_pearson_171213, loop=TRUE)
#plot(net)

region_names_171213 <- unique(regions_171213);
nregions_171213 = length(region_names_171213);
conn_171213 = matrix(0, nregions_171213, nregions_171213);
colnames(conn_171213) <- rownames(conn_171213) <- region_names_171213;
cx = 1;
for (i in region_names_171213){
  idx = which(regions_171213 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_171213){
    idy = which(regions_171213 == j, arr.ind = FALSE);
    conn_171213[cx,cy] = sum(adj_pearson_171213[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_171213, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_171213.csv");
sel_mat <- which(conn_171213 > quantile(conn_171213, .75), arr.ind = TRUE);
net_171213 <- network(conn_171213, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_171213.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_171213, label = unique(rownames(conn_171213)), edge.label = conn_171213, label.cex=2, edge.label.cex=2)
dev.off()


col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171213_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171213_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171213_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171213_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_171213)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_171213, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_evoked_171213.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1:18], track.index = c(2,2), text = "PV HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[19:33], track.index = c(2,2), text = "VM T", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_171213);
idd = which(adj_pearson_171213 > quantile(adj_pearson_171213, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_171213)[1]){
  for(jj in 1:dim(adj_pearson_171213)[2]){
    if (adj_pearson_171213[ii,jj]> quantile(adj_pearson_171213, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_171213[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()







####################################################################################
######################## 180110 #########################
adj_pearson_180110 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_180110_pycharm.csv", header=FALSE)
adj_pearson_180110 <- as.matrix(adj_pearson_180110);
colnames(adj_pearson_180110) <- regions_180110;
rownames(adj_pearson_180110) <- regions_180110;
plot(adj_pearson_180110)
net <- network(adj_pearson_180110, loop=TRUE)
#plot(net)

region_names_180110 <- unique(regions_180110);
nregions_180110 = length(region_names_180110);
conn_180110 = matrix(0, nregions_180110, nregions_180110);
colnames(conn_180110) <- rownames(conn_180110) <- region_names_180110;
cx = 1;
for (i in region_names_180110){
  idx = which(regions_180110 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180110){
    idy = which(regions_180110 == j, arr.ind = FALSE);
    conn_180110[cx,cy] = sum(adj_pearson_180110[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_180110, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_180110.csv");
sel_mat <- which(conn_180110 > quantile(conn_180110, .75), arr.ind = TRUE);
net_180110 <- network(conn_180110, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_180110.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180110, label = unique(rownames(conn_180110)), edge.label = conn_180110, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180110_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180110_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180110_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180110_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180110)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180110, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_ongoing_180110.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1], track.index = c(2,2), text = "A HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[2:7], track.index = c(2,2), text = "DM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[8:28], track.index = c(2,2), text = "PV HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[29:35], track.index = c(2,2), text = "VM T", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180110);
idd = which(adj_pearson_180110 > quantile(adj_pearson_180110, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180110)[1]){
  for(jj in 1:dim(adj_pearson_180110)[2]){
    if (adj_pearson_180110[ii,jj]> quantile(adj_pearson_180110, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_180110[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()



################################### EVOKED ############################################################
adj_pearson_180110 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_180110_pycharm.csv", header=FALSE, sep = "\t")
adj_pearson_180110 <- as.matrix(adj_pearson_180110);
colnames(adj_pearson_180110) <- regions_180110;
rownames(adj_pearson_180110) <- regions_180110;
plot(adj_pearson_180110)
net <- network(adj_pearson_180110, loop=TRUE)
#plot(net)

region_names_180110 <- unique(regions_180110);
nregions_180110 = length(region_names_180110);
conn_180110 = matrix(0, nregions_180110, nregions_180110);
colnames(conn_180110) <- rownames(conn_180110) <- region_names_180110;
cx = 1;
for (i in region_names_180110){
  idx = which(regions_180110 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180110){
    idy = which(regions_180110 == j, arr.ind = FALSE);
    conn_180110[cx,cy] = sum(adj_pearson_180110[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_180110, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_180110.csv");
sel_mat <- which(conn_180110 > quantile(conn_180110, .75), arr.ind = TRUE);
net_180110 <- network(conn_180110, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_180110.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180110, label = unique(rownames(conn_180110)), edge.label = conn_180110, label.cex=2, edge.label.cex=2 )
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180110_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180110_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180110_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180110_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180110)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180110, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_evoked_180110.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1], track.index = c(2,2), text = "A HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[2:7], track.index = c(2,2), text = "DM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[8:28], track.index = c(2,2), text = "PV HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[29:35], track.index = c(2,2), text = "VM T", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180110);
idd = which(adj_pearson_180110 > quantile(adj_pearson_180110, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180110)[1]){
  for(jj in 1:dim(adj_pearson_180110)[2]){
    if (adj_pearson_180110[ii,jj]> quantile(adj_pearson_180110, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_180110[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()










####################################################################################
######################## 180111 #########################
adj_pearson_180111 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_180111_pycharm.csv", header=FALSE)
adj_pearson_180111 <- as.matrix(adj_pearson_180111);
colnames(adj_pearson_180111) <- regions_180111;
rownames(adj_pearson_180111) <- regions_180111;
plot(adj_pearson_180111)
net <- network(adj_pearson_180111, loop=TRUE)
#plot(net)

region_names_180111 <- unique(regions_180111);
nregions_180111 = length(region_names_180111);
conn_180111 = matrix(0, nregions_180111, nregions_180111);
colnames(conn_180111) <- rownames(conn_180111) <- region_names_180111;
cx = 1;
for (i in region_names_180111){
  idx = which(regions_180111 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180111){
    idy = which(regions_180111 == j, arr.ind = FALSE);
    conn_180111[cx,cy] = sum(adj_pearson_180111[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_180111, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_180111.csv");
sel_mat <- which(conn_180111 > quantile(conn_180111, .75), arr.ind = TRUE);
net_180111 <- network(conn_180111, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_180111.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180111, label = unique(rownames(conn_180111)), edge.label = conn_180111, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180111_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180111_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180111_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180111_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180111)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180111, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_ongoing_180111.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1], track.index = c(2,2), text = "PV HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[2:4], track.index = c(2,2), text = "VM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[5:24], track.index = c(2,2), text = "VM T", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180111);
idd = which(adj_pearson_180111 > quantile(adj_pearson_180111, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180111)[1]){
  for(jj in 1:dim(adj_pearson_180111)[2]){
    if (adj_pearson_180111[ii,jj]> quantile(adj_pearson_180111, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_180111[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()



################################### EVOKED ############################################################
adj_pearson_180111 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_180111_pycharm.csv", header=FALSE, sep = ",")
adj_pearson_180111 <- as.matrix(adj_pearson_180111);
colnames(adj_pearson_180111) <- regions_180111;
rownames(adj_pearson_180111) <- regions_180111;
plot(adj_pearson_180111)
net <- network(adj_pearson_180111, loop=TRUE)
#plot(net)

region_names_180111 <- unique(regions_180111);
nregions_180111 = length(region_names_180111);
conn_180111 = matrix(0, nregions_180111, nregions_180111);
colnames(conn_180111) <- rownames(conn_180111) <- region_names_180111;
cx = 1;
for (i in region_names_180111){
  idx = which(regions_180111 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180111){
    idy = which(regions_180111 == j, arr.ind = FALSE);
    conn_180111[cx,cy] = sum(adj_pearson_180111[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_180111, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_180111.csv");
sel_mat <- which(conn_180111 > quantile(conn_180111, .75), arr.ind = TRUE);
net_180111 <- network(conn_180111, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_180111.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180111, label = unique(rownames(conn_180111)), edge.label = conn_180111, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180111_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180111_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180111_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180111_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180111)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180111, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_evoked_180111.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1], track.index = c(2,2), text = "PV HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[2:4], track.index = c(2,2), text = "VM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[5:24], track.index = c(2,2), text = "VM T", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))


for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180111);
idd = which(adj_pearson_180111 > quantile(adj_pearson_180111, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180111)[1]){
  for(jj in 1:dim(adj_pearson_180111)[2]){
    if (adj_pearson_180111[ii,jj]> quantile(adj_pearson_180111, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_180111[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()









####################################################################################
######################## 180131 #########################
adj_pearson_180131 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_180131_pycharm.csv", header=FALSE)
adj_pearson_180131 <- as.matrix(adj_pearson_180131);
colnames(adj_pearson_180131) <- regions_180131;
rownames(adj_pearson_180131) <- regions_180131;
plot(adj_pearson_180131)
net <- network(adj_pearson_180131, loop=TRUE)
#plot(net)

region_names_180131 <- unique(regions_180131);
nregions_180131 = length(region_names_180131);
conn_180131 = matrix(0, nregions_180131, nregions_180131);
colnames(conn_180131) <- rownames(conn_180131) <- region_names_180131;
cx = 1;
for (i in region_names_180131){
  idx = which(regions_180131 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180131){
    idy = which(regions_180131 == j, arr.ind = FALSE);
    conn_180131[cx,cy] = sum(adj_pearson_180131[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_180131, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_180131.csv");
sel_mat <- which(conn_180131 > quantile(conn_180131, .75), arr.ind = TRUE);
net_180131 <- network(conn_180131, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_180131.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180131, label = unique(rownames(conn_180131)), edge.label = conn_180131, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180131_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180131_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180131_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180131_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180131)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180131, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_ongoing_180131.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1:10], track.index = c(2,2), text = "DM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[11], track.index = c(2,2), text = "PV HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[12:22], track.index = c(2,2), text = "P HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[23:31], track.index = c(2,2), text = "VM T", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180131);
idd = which(adj_pearson_180131 > quantile(adj_pearson_180131, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180131)[1]){
  for(jj in 1:dim(adj_pearson_180131)[2]){
    if (adj_pearson_180131[ii,jj]> quantile(adj_pearson_180131, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_180131[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()



################################### EVOKED ############################################################
adj_pearson_180131 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_180131_pycharm.csv", header=FALSE, sep = ",")
adj_pearson_180131 <- as.matrix(adj_pearson_180131);
colnames(adj_pearson_180131) <- regions_180131;
rownames(adj_pearson_180131) <- regions_180131;
plot(adj_pearson_180131)
net <- network(adj_pearson_180131, loop=TRUE)
#plot(net)

region_names_180131 <- unique(regions_180131);
nregions_180131 = length(region_names_180131);
conn_180131 = matrix(0, nregions_180131, nregions_180131);
colnames(conn_180131) <- rownames(conn_180131) <- region_names_180131;
cx = 1;
for (i in region_names_180131){
  idx = which(regions_180131 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180131){
    idy = which(regions_180131 == j, arr.ind = FALSE);
    conn_180131[cx,cy] = sum(adj_pearson_180131[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_180131, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_180131.csv");
sel_mat <- which(conn_180131 > quantile(conn_180131, .75), arr.ind = TRUE);
net_180131 <- network(conn_180131, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_180131.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180131, label = unique(rownames(conn_180131)), edge.label = conn_180131, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180131_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180131_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180131_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180131_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180131)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180131, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_evoked_180131.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1:10], track.index = c(2,2), text = "DM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[11], track.index = c(2,2), text = "PV HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[12:22], track.index = c(2,2), text = "P HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[23:31], track.index = c(2,2), text = "VM T", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180131);
idd = which(adj_pearson_180131 > quantile(adj_pearson_180131, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180131)[1]){
  for(jj in 1:dim(adj_pearson_180131)[2]){
    if (adj_pearson_180131[ii,jj]> quantile(adj_pearson_180131, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_180131[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()













####################################################################################
######################## 180221 #########################
adj_pearson_180221 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_180221_pycharm.csv", header=FALSE)
adj_pearson_180221 <- as.matrix(adj_pearson_180221);
colnames(adj_pearson_180221) <- regions_180221;
rownames(adj_pearson_180221) <- regions_180221;
plot(adj_pearson_180221)
net <- network(adj_pearson_180221, loop=TRUE)
#plot(net)

region_names_180221 <- unique(regions_180221);
nregions_180221 = length(region_names_180221);
conn_180221 = matrix(0, nregions_180221, nregions_180221);
colnames(conn_180221) <- rownames(conn_180221) <- region_names_180221;
cx = 1;
for (i in region_names_180221){
  idx = which(regions_180221 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180221){
    idy = which(regions_180221 == j, arr.ind = FALSE);
    conn_180221[cx,cy] = sum(adj_pearson_180221[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_180221, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_180221.csv");
sel_mat <- which(conn_180221 > quantile(conn_180221, .75), arr.ind = TRUE);
net_180221 <- network(conn_180221, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_180221.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180221, label = unique(rownames(conn_180221)), edge.label = conn_180221, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180221_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180221_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180221_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180221_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180221)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180221, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_ongoing_180221.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1:17], track.index = c(2,2), text = "DM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[18:24], track.index = c(2,2), text = "PV HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[25:27], track.index = c(2,2), text = "VM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[28:54], track.index = c(2,2), text = "VM T", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180221);
idd = which(adj_pearson_180221 > quantile(adj_pearson_180221, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180221)[1]){
  for(jj in 1:dim(adj_pearson_180221)[2]){
    if (adj_pearson_180221[ii,jj]> quantile(adj_pearson_180221, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_180221[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()



################################### EVOKED ############################################################
adj_pearson_180221 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_180221_pycharm.csv", header=FALSE, sep = ",")
adj_pearson_180221 <- as.matrix(adj_pearson_180221);
colnames(adj_pearson_180221) <- regions_180221;
rownames(adj_pearson_180221) <- regions_180221;
plot(adj_pearson_180221)
net <- network(adj_pearson_180221, loop=TRUE)
#plot(net)

region_names_180221 <- unique(regions_180221);
nregions_180221 = length(region_names_180221);
conn_180221 = matrix(0, nregions_180221, nregions_180221);
colnames(conn_180221) <- rownames(conn_180221) <- region_names_180221;
cx = 1;
for (i in region_names_180221){
  idx = which(regions_180221 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180221){
    idy = which(regions_180221 == j, arr.ind = FALSE);
    conn_180221[cx,cy] = sum(adj_pearson_180221[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_180221, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_180221.csv");
sel_mat <- which(conn_180221 > quantile(conn_180221, .75), arr.ind = TRUE);
net_180221 <- network(conn_180221, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_180221.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180221, label = unique(rownames(conn_180221)), edge.label = conn_180221, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180221_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180221_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180221_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180221_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180221)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180221, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_evoked_180221.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1:17], track.index = c(2,2), text = "DM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[18:24], track.index = c(2,2), text = "PV HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[25:27], track.index = c(2,2), text = "VM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[28:54], track.index = c(2,2), text = "VM T", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180221);
idd = which(adj_pearson_180221 > quantile(adj_pearson_180221, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180221)[1]){
  for(jj in 1:dim(adj_pearson_180221)[2]){
    if (adj_pearson_180221[ii,jj]> quantile(adj_pearson_180221, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_180221[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()







####################################################################################
######################## 180228 #########################
adj_pearson_180228 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_180228_pycharm.csv", header=FALSE)
adj_pearson_180228 <- as.matrix(adj_pearson_180228);
colnames(adj_pearson_180228) <- regions_180228;
rownames(adj_pearson_180228) <- regions_180228;
plot(adj_pearson_180228)
net <- network(adj_pearson_180228, loop=TRUE)
#plot(net)

region_names_180228 <- unique(regions_180228);
nregions_180228 = length(region_names_180228);
conn_180228 = matrix(0, nregions_180228, nregions_180228);
colnames(conn_180228) <- rownames(conn_180228) <- region_names_180228;
cx = 1;
for (i in region_names_180228){
  idx = which(regions_180228 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180228){
    idy = which(regions_180228 == j, arr.ind = FALSE);
    conn_180228[cx,cy] = sum(adj_pearson_180228[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_180228, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_180228.csv");
sel_mat <- which(conn_180228 > quantile(conn_180228, .75), arr.ind = TRUE);
net_180228 <- network(conn_180228, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_180228.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180228, label = unique(rownames(conn_180228)), edge.label = conn_180228, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180228_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180228_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180228_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180228_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180228)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180228, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_ongoing_180228.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1:9], track.index = c(2,2), text = "Arc HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[10:38], track.index = c(2,2), text = "DM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[39:56], track.index = c(2,2), text = "P HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[57:62], track.index = c(2,2), text = "VM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[63:69], track.index = c(2,2), text = "VM T", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[70:79], track.index = c(2,2), text = "ZI", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180228);
idd = which(adj_pearson_180228 > quantile(adj_pearson_180228, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180228)[1]){
  for(jj in 1:dim(adj_pearson_180228)[2]){
    if (adj_pearson_180228[ii,jj]> quantile(adj_pearson_180228, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_180228[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()



################################### EVOKED ############################################################
adj_pearson_180228 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_180228_pycharm.csv", header=FALSE, sep = ",")
adj_pearson_180228 <- as.matrix(adj_pearson_180228);
colnames(adj_pearson_180228) <- regions_180228;
rownames(adj_pearson_180228) <- regions_180228;
plot(adj_pearson_180228)
net <- network(adj_pearson_180228, loop=TRUE)
#plot(net)

region_names_180228 <- unique(regions_180228);
nregions_180228 = length(region_names_180228);
conn_180228 = matrix(0, nregions_180228, nregions_180228);
colnames(conn_180228) <- rownames(conn_180228) <- region_names_180228;
cx = 1;
for (i in region_names_180228){
  idx = which(regions_180228 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180228){
    idy = which(regions_180228 == j, arr.ind = FALSE);
    conn_180228[cx,cy] = sum(adj_pearson_180228[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_180228, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_180228.csv");
sel_mat <- which(conn_180228 > quantile(conn_180228, .75), arr.ind = TRUE);
net_180228 <- network(conn_180228, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_180228.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180228, label = unique(rownames(conn_180228)), edge.label = conn_180228, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180228_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180228_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180228_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180228_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180228)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180228, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_evoked_180228.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1:9], track.index = c(2,2), text = "Arc HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[10:38], track.index = c(2,2), text = "DM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[39:56], track.index = c(2,2), text = "P HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[57:62], track.index = c(2,2), text = "VM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[63:69], track.index = c(2,2), text = "VM T", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[70:79], track.index = c(2,2), text = "ZI", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180228);
idd = which(adj_pearson_180228 > quantile(adj_pearson_180228, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180228)[1]){
  for(jj in 1:dim(adj_pearson_180228)[2]){
    if (adj_pearson_180228[ii,jj]> quantile(adj_pearson_180228, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_180228[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()






####################################################################################
######################## 180302 #########################
adj_pearson_180302 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_180302_pycharm.csv", header=FALSE)
adj_pearson_180302 <- as.matrix(adj_pearson_180302);
colnames(adj_pearson_180302) <- regions_180302;
rownames(adj_pearson_180302) <- regions_180302;
plot(adj_pearson_180302)
net <- network(adj_pearson_180302, loop=TRUE)
#plot(net)

region_names_180302 <- unique(regions_180302);
nregions_180302 = length(region_names_180302);
conn_180302 = matrix(0, nregions_180302, nregions_180302);
colnames(conn_180302) <- rownames(conn_180302) <- region_names_180302;
cx = 1;
for (i in region_names_180302){
  idx = which(regions_180302 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180302){
    idy = which(regions_180302 == j, arr.ind = FALSE);
    conn_180302[cx,cy] = sum(adj_pearson_180302[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_180302, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_180302.csv");
sel_mat <- which(conn_180302 > quantile(conn_180302, .75), arr.ind = TRUE);
net_180302 <- network(conn_180302, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_180302.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180302, label = unique(rownames(conn_180302)), edge.label = conn_180302, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180302_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180302_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180302_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180302_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180302)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180302, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_ongoing_180302.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1:7], track.index = c(2,2), text = "Arc HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[8:11], track.index = c(2,2), text = "DM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[12:29], track.index = c(2,2), text = "MM Complex", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[30:33], track.index = c(2,2), text = "PV HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[34:47], track.index = c(2,2), text = "P HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))


for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180302);
idd = which(adj_pearson_180302 > quantile(adj_pearson_180302, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180302)[1]){
  for(jj in 1:dim(adj_pearson_180302)[2]){
    if (adj_pearson_180302[ii,jj]> quantile(adj_pearson_180302, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_180302[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()



################################### EVOKED ############################################################
adj_pearson_180302 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_180302_pycharm.csv", header=FALSE, sep = ",")
adj_pearson_180302 <- as.matrix(adj_pearson_180302);
colnames(adj_pearson_180302) <- regions_180302;
rownames(adj_pearson_180302) <- regions_180302;
plot(adj_pearson_180302)
net <- network(adj_pearson_180302, loop=TRUE)
#plot(net)

region_names_180302 <- unique(regions_180302);
nregions_180302 = length(region_names_180302);
conn_180302 = matrix(0, nregions_180302, nregions_180302);
colnames(conn_180302) <- rownames(conn_180302) <- region_names_180302;
cx = 1;
for (i in region_names_180302){
  idx = which(regions_180302 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180302){
    idy = which(regions_180302 == j, arr.ind = FALSE);
    conn_180302[cx,cy] = sum(adj_pearson_180302[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_180302, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_180302.csv");
sel_mat <- which(conn_180302 > quantile(conn_180302, .75), arr.ind = TRUE);
net_180302 <- network(conn_180302, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_180302.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180302, label = unique(rownames(conn_180302)), edge.label = conn_180302, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180302_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180302_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180302_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180302_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180302)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180302, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_evoked_180302.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1:7], track.index = c(2,2), text = "Arc HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[8:11], track.index = c(2,2), text = "DM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[12:29], track.index = c(2,2), text = "MM Complex", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[30:33], track.index = c(2,2), text = "PV HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[34:47], track.index = c(2,2), text = "P HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180302);
idd = which(adj_pearson_180302 > quantile(adj_pearson_180302, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180302)[1]){
  for(jj in 1:dim(adj_pearson_180302)[2]){
    if (adj_pearson_180302[ii,jj]> quantile(adj_pearson_180302, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_180302[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()








####################################################################################
######################## 180419 #########################
adj_pearson_180419 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_180419_pycharm.csv", header=FALSE)
adj_pearson_180419 <- as.matrix(adj_pearson_180419);
colnames(adj_pearson_180419) <- regions_180419;
rownames(adj_pearson_180419) <- regions_180419;
plot(adj_pearson_180419)
net <- network(adj_pearson_180419, loop=TRUE)
#plot(net)

region_names_180419 <- unique(regions_180419);
nregions_180419 = length(region_names_180419);
conn_180419 = matrix(0, nregions_180419, nregions_180419);
colnames(conn_180419) <- rownames(conn_180419) <- region_names_180419;
cx = 1;
for (i in region_names_180419){
  idx = which(regions_180419 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180419){
    idy = which(regions_180419 == j, arr.ind = FALSE);
    conn_180419[cx,cy] = sum(adj_pearson_180419[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_180419, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_180419.csv");
sel_mat <- which(conn_180419 > quantile(conn_180419, .75), arr.ind = TRUE);
net_180419 <- network(conn_180419, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_180419.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180419, label = unique(rownames(conn_180419)), edge.label = conn_180419, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180419_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180419_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180419_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180419_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180419)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180419, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_ongoing_180419.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1:35], track.index = c(2,2), text = "DM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[36:50], track.index = c(2,2), text = "P HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[51:56], track.index = c(2,2), text = "VM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[57:61], track.index = c(2,2), text = "VM T", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180419);
idd = which(adj_pearson_180419 > quantile(adj_pearson_180419, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180419)[1]){
  for(jj in 1:dim(adj_pearson_180419)[2]){
    if (adj_pearson_180419[ii,jj]> quantile(adj_pearson_180419, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_180419[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()



################################### EVOKED ############################################################
adj_pearson_180419 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_180419_pycharm.csv", header=FALSE, sep = ",")
adj_pearson_180419 <- as.matrix(adj_pearson_180419);
colnames(adj_pearson_180419) <- regions_180419;
rownames(adj_pearson_180419) <- regions_180419;
plot(adj_pearson_180419)
net <- network(adj_pearson_180419, loop=TRUE)
#plot(net)

region_names_180419 <- unique(regions_180419);
nregions_180419 = length(region_names_180419);
conn_180419 = matrix(0, nregions_180419, nregions_180419);
colnames(conn_180419) <- rownames(conn_180419) <- region_names_180419;
cx = 1;
for (i in region_names_180419){
  idx = which(regions_180419 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180419){
    idy = which(regions_180419 == j, arr.ind = FALSE);
    conn_180419[cx,cy] = sum(adj_pearson_180419[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_180419, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_180419.csv");
sel_mat <- which(conn_180419 > quantile(conn_180419, .75), arr.ind = TRUE);
net_180419 <- network(conn_180419, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_180419.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180419, label = unique(rownames(conn_180419)), edge.label = conn_180419, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180419_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180419_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180419_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180419_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180419)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180419, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_evoked_180419.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1:35], track.index = c(2,2), text = "DM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[36:50], track.index = c(2,2), text = "P HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[51:56], track.index = c(2,2), text = "VM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[57:61], track.index = c(2,2), text = "VM T", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180419);
idd = which(adj_pearson_180419 > quantile(adj_pearson_180419, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180419)[1]){
  for(jj in 1:dim(adj_pearson_180419)[2]){
    if (adj_pearson_180419[ii,jj]> quantile(adj_pearson_180419, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_180419[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()








####################################################################################
######################## 180420 #########################
adj_pearson_180420 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_180420_pycharm.csv", header=FALSE)
adj_pearson_180420 <- as.matrix(adj_pearson_180420);
colnames(adj_pearson_180420) <- regions_180420;
rownames(adj_pearson_180420) <- regions_180420;
plot(adj_pearson_180420)
net <- network(adj_pearson_180420, loop=TRUE)
#plot(net)

region_names_180420 <- unique(regions_180420);
nregions_180420 = length(region_names_180420);
conn_180420 = matrix(0, nregions_180420, nregions_180420);
colnames(conn_180420) <- rownames(conn_180420) <- region_names_180420;
cx = 1;
for (i in region_names_180420){
  idx = which(regions_180420 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180420){
    idy = which(regions_180420 == j, arr.ind = FALSE);
    conn_180420[cx,cy] = sum(adj_pearson_180420[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_180420, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_180420.csv");
sel_mat <- which(conn_180420 > quantile(conn_180420, .75), arr.ind = TRUE);
net_180420 <- network(conn_180420, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_180420.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180420, label = unique(rownames(conn_180420)), edge.label = conn_180420, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180420_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180420_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180420_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180420_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180420)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180420, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_ongoing_180420.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1], track.index = c(2,2), text = "DM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[2:27], track.index = c(2,2), text = "MM Complex", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[28:46], track.index = c(2,2), text = "PV HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[47:71], track.index = c(2,2), text = "P HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180420);
idd = which(adj_pearson_180420 > quantile(adj_pearson_180420, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180420)[1]){
  for(jj in 1:dim(adj_pearson_180420)[2]){
    if (adj_pearson_180420[ii,jj]> quantile(adj_pearson_180420, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_180420[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()



################################### EVOKED ############################################################
adj_pearson_180420 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_180420_pycharm.csv", header=FALSE, sep = ",")
adj_pearson_180420 <- as.matrix(adj_pearson_180420);
colnames(adj_pearson_180420) <- regions_180420;
rownames(adj_pearson_180420) <- regions_180420;
plot(adj_pearson_180420)
net <- network(adj_pearson_180420, loop=TRUE)
#plot(net)

region_names_180420 <- unique(regions_180420);
nregions_180420 = length(region_names_180420);
conn_180420 = matrix(0, nregions_180420, nregions_180420);
colnames(conn_180420) <- rownames(conn_180420) <- region_names_180420;
cx = 1;
for (i in region_names_180420){
  idx = which(regions_180420 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180420){
    idy = which(regions_180420 == j, arr.ind = FALSE);
    conn_180420[cx,cy] = sum(adj_pearson_180420[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_180420, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_180420.csv");
sel_mat <- which(conn_180420 > quantile(conn_180420, .75), arr.ind = TRUE);
net_180420 <- network(conn_180420, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_180420.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180420, label = unique(rownames(conn_180420)) , edge.label = conn_180420, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180420_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180420_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180420_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180420_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180420)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180420, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_evoked_180420.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1], track.index = c(2,2), text = "DM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[2:27], track.index = c(2,2), text = "MM Complex", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[28:46], track.index = c(2,2), text = "PV HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[47:71], track.index = c(2,2), text = "P HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180420);
idd = which(adj_pearson_180420 > quantile(adj_pearson_180420, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180420)[1]){
  for(jj in 1:dim(adj_pearson_180420)[2]){
    if (adj_pearson_180420[ii,jj]> quantile(adj_pearson_180420, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_180420[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()




####################################################################################
######################## 180423 #########################
adj_pearson_180423 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_180423_pycharm.csv", header=FALSE)
adj_pearson_180423 <- as.matrix(adj_pearson_180423);
colnames(adj_pearson_180423) <- regions_180423;
rownames(adj_pearson_180423) <- regions_180423;
plot(adj_pearson_180423)
net <- network(adj_pearson_180423, loop=TRUE)
#plot(net)

region_names_180423 <- unique(regions_180423);
nregions_180423 = length(region_names_180423);
conn_180423 = matrix(0, nregions_180423, nregions_180423);
colnames(conn_180423) <- rownames(conn_180423) <- region_names_180423;
cx = 1;
for (i in region_names_180423){
  idx = which(regions_180423 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180423){
    idy = which(regions_180423 == j, arr.ind = FALSE);
    conn_180423[cx,cy] = sum(adj_pearson_180423[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_180423, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_180423.csv");
sel_mat <- which(conn_180423 > quantile(conn_180423, .75), arr.ind = TRUE);
net_180423 <- network(conn_180423, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_180423.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180423, label = unique(rownames(conn_180423)), edge.label = conn_180423, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180423_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180423_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180423_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180423_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180423)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180423, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_ongoing_180423.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1:19], track.index = c(2,2), text = "DM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[20:34], track.index = c(2,2), text = "MM Complex", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[35:52], track.index = c(2,2), text = "P HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[53:56], track.index = c(2,2), text = "VM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180423);
idd = which(adj_pearson_180423 > quantile(adj_pearson_180423, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180423)[1]){
  for(jj in 1:dim(adj_pearson_180423)[2]){
    if (adj_pearson_180423[ii,jj]> quantile(adj_pearson_180423, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_180423[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()



################################### EVOKED ############################################################
adj_pearson_180423 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_180423_pycharm.csv", header=FALSE, sep = ",")
adj_pearson_180423 <- as.matrix(adj_pearson_180423);
colnames(adj_pearson_180423) <- regions_180423;
rownames(adj_pearson_180423) <- regions_180423;
plot(adj_pearson_180423)
net <- network(adj_pearson_180423, loop=TRUE)
#plot(net)

region_names_180423 <- unique(regions_180423);
nregions_180423 = length(region_names_180423);
conn_180423 = matrix(0, nregions_180423, nregions_180423);
colnames(conn_180423) <- rownames(conn_180423) <- region_names_180423;
cx = 1;
for (i in region_names_180423){
  idx = which(regions_180423 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180423){
    idy = which(regions_180423 == j, arr.ind = FALSE);
    conn_180423[cx,cy] = sum(adj_pearson_180423[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
write.csv(conn_180423, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_180423.csv");
sel_mat <- which(conn_180423 > quantile(conn_180423, .75), arr.ind = TRUE);
net_180423 <- network(conn_180423, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_180423.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180423, label = unique(rownames(conn_180423)) , edge.label = conn_180423, label.cex=2, edge.label.cex=2)
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180423_fr")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180423_clus_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180423_deg_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180423_localeff_wei")
df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180423)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180423, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/regionwise_evoked_180423.png"), units = "in",  width = 16, height = 16, res = 600)
circos.clear()
sectors = factor(tt, levels = tt)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = fr_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = clus_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = deg_col, cell.padding = c(0.01, 0.01, 0.01, 0.01))
set_track_gap(cm_h(0.001))
circos.track(sectors, ylim = c(0, 2), bg.border = "black", track.height = cm_h(0.4), bg.lwd = 0.1, bg.col = loceff_col, cell.padding = c(0.01, 0.1, 0.01, 0.1))

highlight.sector(tt[1:19], track.index = c(2,2), text = "DM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[20:34], track.index = c(2,2), text = "MM Complex", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[35:52], track.index = c(2,2), text = "P HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[53:56], track.index = c(2,2), text = "VM HT", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180423);
idd = which(adj_pearson_180423 > quantile(adj_pearson_180423, 0.99), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180423)[1]){
  for(jj in 1:dim(adj_pearson_180423)[2]){
    if (adj_pearson_180423[ii,jj]> quantile(adj_pearson_180423, 0.99)){
      circos.link(tt[ii], c(0, 1), tt[jj] , c(0, 1), directional = 0, col = col_fun_edges(adj_pearson_180423[ii,jj]) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "RdBu")), transparency = 0)
leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
               col_fun = col_fun_edges, title = "Edges", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
col_fun_metric = colorRamp2( seq(0,1, 0.1), rev(brewer.pal(n = 11, name = "PRGn")), transparency = 0)
leg2 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"),
               col_fun = col_fun_metric, title = "Metrics", grid_width = unit(0.4, "in"),
               labels_gp = gpar(col = "black", fontsize = 14),
               title_gp = gpar(col = "black", fontsize = 18),
               title_gap = unit(4, "mm"))
ComplexHeatmap::draw(packLegend(leg1, leg2, direction = "horizontal"), x = unit(2, "in"), y = unit(2, "in")) #just = c("right", "top")

dev.off()

########################## FIX STRENGTH ALL ############################################################################################
strength_all = data.frame("Ongoing" = rep(0,13), "LightON" = rep(0,13))
i = 1;
for (day in days){
  df <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", 
                  table=paste0(day,"_deg_wei", sep = ""))
  df$Ongoing <- as.numeric(sub(",", ".",df$Ongoing, fixed=TRUE))
  df$LightON <- as.numeric(sub(",", ".",df$LightON, fixed=TRUE))
  strength_all$Ongoing[i] = mean(df$Ongoing)
  strength_all$LightON[i] = mean(df$LightON)
  i = i + 1;
}
####OK 
########################################################################################################################################
############### DIFFERENCE NETWORKS #####################################################################################################
adj_pearson_171019 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_171019_pycharm.csv", header=FALSE)
adj_pearson_171019 <- as.matrix(adj_pearson_171019);
colnames(adj_pearson_171019) <- regions_171019;
rownames(adj_pearson_171019) <- regions_171019;
net <- network(adj_pearson_171019, loop=TRUE)


region_names_171019 <- unique(regions_171019);
nregions_171019 = length(region_names_171019);
conn_171019 = matrix(0, nregions_171019, nregions_171019);
colnames(conn_171019) <- rownames(conn_171019) <- region_names_171019;
cx = 1;
for (i in region_names_171019){
  idx = which(regions_171019 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_171019){
    idy = which(regions_171019 == j, arr.ind = FALSE);
    conn_171019[cx,cy] = sum(adj_pearson_171019[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_o_171019 = conn_171019;
adj_pearson_171019 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_171019_pycharm.csv", header=FALSE)
adj_pearson_171019 <- as.matrix(adj_pearson_171019);
colnames(adj_pearson_171019) <- regions_171019;
rownames(adj_pearson_171019) <- regions_171019;
plot(adj_pearson_171019)
net <- network(adj_pearson_171019, loop=TRUE)

region_names_171019 <- unique(regions_171019);
nregions_171019 = length(region_names_171019);
conn_171019 = matrix(0, nregions_171019, nregions_171019);
colnames(conn_171019) <- rownames(conn_171019) <- region_names_171019;
cx = 1;
for (i in region_names_171019){
  idx = which(regions_171019 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_171019){
    idy = which(regions_171019 == j, arr.ind = FALSE);
    conn_171019[cx,cy] = sum(adj_pearson_171019[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_e_171019 = conn_171019;

conn_diff_171019 = conn_e_171019 - conn_o_171019
net_171019 <- network(conn_diff_171019, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_171019.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_171019, label = unique(rownames(conn_diff_171019)), edge.label = conn_diff_171019, label.cex=2, edge.label.cex=2,
             edge.lwd = conn_diff_171019/50)
dev.off()



adj_pearson_171207 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_171207_pycharm.csv", header=FALSE)
adj_pearson_171207 <- as.matrix(adj_pearson_171207);
colnames(adj_pearson_171207) <- regions_171207;
rownames(adj_pearson_171207) <- regions_171207;
#plot(adj_pearson_171207)
net <- network(adj_pearson_171207, loop=TRUE)
#plot(net)

region_names_171207 <- unique(regions_171207);
nregions_171207 = length(region_names_171207);
conn_171207 = matrix(0, nregions_171207, nregions_171207);
colnames(conn_171207) <- rownames(conn_171207) <- region_names_171207;
cx = 1;
for (i in region_names_171207){
  idx = which(regions_171207 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_171207){
    idy = which(regions_171207 == j, arr.ind = FALSE);
    conn_171207[cx,cy] = sum(adj_pearson_171207[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_o_171207 = conn_171207;
adj_pearson_171207 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_171207_pycharm.csv", header=FALSE)
adj_pearson_171207 <- as.matrix(adj_pearson_171207);
colnames(adj_pearson_171207) <- regions_171207;
rownames(adj_pearson_171207) <- regions_171207;
plot(adj_pearson_171207)
net <- network(adj_pearson_171207, loop=TRUE)
#plot(net)

region_names_171207 <- unique(regions_171207);
nregions_171207 = length(region_names_171207);
conn_171207 = matrix(0, nregions_171207, nregions_171207);
colnames(conn_171207) <- rownames(conn_171207) <- region_names_171207;
cx = 1;
for (i in region_names_171207){
  idx = which(regions_171207 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_171207){
    idy = which(regions_171207 == j, arr.ind = FALSE);
    conn_171207[cx,cy] = sum(adj_pearson_171207[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_e_171207 = conn_171207;

conn_diff_171207 = conn_e_171207 - conn_o_171207
net_171207 <- network(conn_diff_171207, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_171207.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_171207, label = unique(rownames(conn_diff_171207)), edge.label = conn_diff_171207, label.cex = 2, edge.label.cex = 2
             ,edge.lwd = conn_diff_171207/20)
dev.off()




adj_pearson_171208 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_171208_pycharm.csv", header=FALSE)
adj_pearson_171208 <- as.matrix(adj_pearson_171208);
colnames(adj_pearson_171208) <- regions_171208;
rownames(adj_pearson_171208) <- regions_171208;
plot(adj_pearson_171208)
net <- network(adj_pearson_171208, loop=TRUE)
#plot(net)

region_names_171208 <- unique(regions_171208);
nregions_171208 = length(region_names_171208);
conn_171208 = matrix(0, nregions_171208, nregions_171208);
colnames(conn_171208) <- rownames(conn_171208) <- region_names_171208;
cx = 1;
for (i in region_names_171208){
  idx = which(regions_171208 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_171208){
    idy = which(regions_171208 == j, arr.ind = FALSE);
    conn_171208[cx,cy] = sum(adj_pearson_171208[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_o_171208 = conn_171208;
adj_pearson_171208 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_171208_pycharm.csv", header=FALSE)
adj_pearson_171208 <- as.matrix(adj_pearson_171208);
colnames(adj_pearson_171208) <- regions_171208;
rownames(adj_pearson_171208) <- regions_171208;
plot(adj_pearson_171208)
net <- network(adj_pearson_171208, loop=TRUE)
#plot(net)

region_names_171208 <- unique(regions_171208);
nregions_171208 = length(region_names_171208);
conn_171208 = matrix(0, nregions_171208, nregions_171208);
colnames(conn_171208) <- rownames(conn_171208) <- region_names_171208;
cx = 1;
for (i in region_names_171208){
  idx = which(regions_171208 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_171208){
    idy = which(regions_171208 == j, arr.ind = FALSE);
    conn_171208[cx,cy] = sum(adj_pearson_171208[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_e_171208 = conn_171208;

conn_diff_171208 = conn_e_171208 - conn_o_171208
net_171208 <- network(conn_diff_171208, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_171208.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_171208, label = unique(rownames(conn_diff_171208)), edge.label = conn_diff_171208, label.cex = 2, edge.label.cex = 2,
             edge.lwd = conn_diff_171208/10)
dev.off()


adj_pearson_171213 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_171213_pycharm.csv", header=FALSE)
adj_pearson_171213 <- as.matrix(adj_pearson_171213);
colnames(adj_pearson_171213) <- regions_171213;
rownames(adj_pearson_171213) <- regions_171213;
plot(adj_pearson_171213)
net <- network(adj_pearson_171213, loop=TRUE)
#plot(net)

region_names_171213 <- unique(regions_171213);
nregions_171213 = length(region_names_171213);
conn_171213 = matrix(0, nregions_171213, nregions_171213);
colnames(conn_171213) <- rownames(conn_171213) <- region_names_171213;
cx = 1;
for (i in region_names_171213){
  idx = which(regions_171213 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_171213){
    idy = which(regions_171213 == j, arr.ind = FALSE);
    conn_171213[cx,cy] = sum(adj_pearson_171213[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_o_171213 = conn_171213;
adj_pearson_171213 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_171213_pycharm.csv", header=FALSE)
adj_pearson_171213 <- as.matrix(adj_pearson_171213);
colnames(adj_pearson_171213) <- regions_171213;
rownames(adj_pearson_171213) <- regions_171213;
plot(adj_pearson_171213)
net <- network(adj_pearson_171213, loop=TRUE)
#plot(net)

region_names_171213 <- unique(regions_171213);
nregions_171213 = length(region_names_171213);
conn_171213 = matrix(0, nregions_171213, nregions_171213);
colnames(conn_171213) <- rownames(conn_171213) <- region_names_171213;
cx = 1;
for (i in region_names_171213){
  idx = which(regions_171213 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_171213){
    idy = which(regions_171213 == j, arr.ind = FALSE);
    conn_171213[cx,cy] = sum(adj_pearson_171213[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_e_171213 = conn_171213;

conn_diff_171213 = conn_e_171213 - conn_o_171213
net_171213 <- network(conn_diff_171213, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_171213.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_171213, label = unique(rownames(conn_diff_171213)), edge.label = conn_diff_171213, label.cex = 2, edge.label.cex = 2
             ,edge.lwd = conn_diff_171213/50)
dev.off()



adj_pearson_180110 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_180110_pycharm.csv", header=FALSE)
adj_pearson_180110 <- as.matrix(adj_pearson_180110);
colnames(adj_pearson_180110) <- regions_180110;
rownames(adj_pearson_180110) <- regions_180110;
plot(adj_pearson_180110)
net <- network(adj_pearson_180110, loop=TRUE)
#plot(net)

region_names_180110 <- unique(regions_180110);
nregions_180110 = length(region_names_180110);
conn_180110 = matrix(0, nregions_180110, nregions_180110);
colnames(conn_180110) <- rownames(conn_180110) <- region_names_180110;
cx = 1;
for (i in region_names_180110){
  idx = which(regions_180110 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180110){
    idy = which(regions_180110 == j, arr.ind = FALSE);
    conn_180110[cx,cy] = sum(adj_pearson_180110[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_o_180110 = conn_180110;
adj_pearson_180110 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_180110_pycharm.csv", header=FALSE, sep = "\t")
adj_pearson_180110 <- as.matrix(adj_pearson_180110);
colnames(adj_pearson_180110) <- regions_180110;
rownames(adj_pearson_180110) <- regions_180110;
plot(adj_pearson_180110)
net <- network(adj_pearson_180110, loop=TRUE)
#plot(net)

region_names_180110 <- unique(regions_180110);
nregions_180110 = length(region_names_180110);
conn_180110 = matrix(0, nregions_180110, nregions_180110);
colnames(conn_180110) <- rownames(conn_180110) <- region_names_180110;
cx = 1;
for (i in region_names_180110){
  idx = which(regions_180110 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180110){
    idy = which(regions_180110 == j, arr.ind = FALSE);
    conn_180110[cx,cy] = sum(adj_pearson_180110[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_e_180110 = conn_180110;

conn_diff_180110 = conn_e_180110 - conn_o_180110
net_180110 <- network(conn_diff_180110, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_180110.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180110, label = unique(rownames(conn_diff_180110)), edge.label = conn_diff_180110, label.cex = 2, edge.label.cex = 2,
             edge.lwd = conn_diff_180110/40)
dev.off()



adj_pearson_180111 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_180111_pycharm.csv", header=FALSE)
adj_pearson_180111 <- as.matrix(adj_pearson_180111);
colnames(adj_pearson_180111) <- regions_180111;
rownames(adj_pearson_180111) <- regions_180111;
plot(adj_pearson_180111)
net <- network(adj_pearson_180111, loop=TRUE)
#plot(net)

region_names_180111 <- unique(regions_180111);
nregions_180111 = length(region_names_180111);
conn_180111 = matrix(0, nregions_180111, nregions_180111);
colnames(conn_180111) <- rownames(conn_180111) <- region_names_180111;
cx = 1;
for (i in region_names_180111){
  idx = which(regions_180111 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180111){
    idy = which(regions_180111 == j, arr.ind = FALSE);
    conn_180111[cx,cy] = sum(adj_pearson_180111[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_o_180111 = conn_180111;
adj_pearson_180111 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_180111_pycharm.csv", header=FALSE)
adj_pearson_180111 <- as.matrix(adj_pearson_180111);
colnames(adj_pearson_180111) <- regions_180111;
rownames(adj_pearson_180111) <- regions_180111;
plot(adj_pearson_180111)
net <- network(adj_pearson_180111, loop=TRUE)
#plot(net)

region_names_180111 <- unique(regions_180111);
nregions_180111 = length(region_names_180111);
conn_180111 = matrix(0, nregions_180111, nregions_180111);
colnames(conn_180111) <- rownames(conn_180111) <- region_names_180111;
cx = 1;
for (i in region_names_180111){
  idx = which(regions_180111 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180111){
    idy = which(regions_180111 == j, arr.ind = FALSE);
    conn_180111[cx,cy] = sum(adj_pearson_180111[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_e_180111 = conn_180111;

conn_diff_180111 = conn_e_180111 - conn_o_180111
net_180111 <- network(conn_diff_180111, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_180111.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180111, label = unique(rownames(conn_diff_180111)), edge.label = conn_diff_180111, label.cex = 2, edge.label.cex = 2,
             edge.lwd = conn_diff_180111/20)
dev.off()



adj_pearson_180131 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_180131_pycharm.csv", header=FALSE)
adj_pearson_180131 <- as.matrix(adj_pearson_180131);
colnames(adj_pearson_180131) <- regions_180131;
rownames(adj_pearson_180131) <- regions_180131;
plot(adj_pearson_180131)
net <- network(adj_pearson_180131, loop=TRUE)
#plot(net)

region_names_180131 <- unique(regions_180131);
nregions_180131 = length(region_names_180131);
conn_180131 = matrix(0, nregions_180131, nregions_180131);
colnames(conn_180131) <- rownames(conn_180131) <- region_names_180131;
cx = 1;
for (i in region_names_180131){
  idx = which(regions_180131 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180131){
    idy = which(regions_180131 == j, arr.ind = FALSE);
    conn_180131[cx,cy] = sum(adj_pearson_180131[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_o_180131 = conn_180131;
adj_pearson_180131 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_180131_pycharm.csv", header=FALSE)
adj_pearson_180131 <- as.matrix(adj_pearson_180131);
colnames(adj_pearson_180131) <- regions_180131;
rownames(adj_pearson_180131) <- regions_180131;
plot(adj_pearson_180131)
net <- network(adj_pearson_180131, loop=TRUE)
#plot(net)

region_names_180131 <- unique(regions_180131);
nregions_180131 = length(region_names_180131);
conn_180131 = matrix(0, nregions_180131, nregions_180131);
colnames(conn_180131) <- rownames(conn_180131) <- region_names_180131;
cx = 1;
for (i in region_names_180131){
  idx = which(regions_180131 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180131){
    idy = which(regions_180131 == j, arr.ind = FALSE);
    conn_180131[cx,cy] = sum(adj_pearson_180131[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_e_180131 = conn_180131;

conn_diff_180131 = conn_e_180131 - conn_o_180131
net_180131 <- network(conn_diff_180131, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_180131.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180131, label = unique(rownames(conn_diff_180131)), edge.label = conn_diff_180131, label.cex = 2, edge.label.cex = 2,
             edge.lwd = conn_diff_180131/30)
dev.off()



adj_pearson_180221 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_180221_pycharm.csv", header=FALSE)
adj_pearson_180221 <- as.matrix(adj_pearson_180221);
colnames(adj_pearson_180221) <- regions_180221;
rownames(adj_pearson_180221) <- regions_180221;
plot(adj_pearson_180221)
net <- network(adj_pearson_180221, loop=TRUE)
#plot(net)

region_names_180221 <- unique(regions_180221);
nregions_180221 = length(region_names_180221);
conn_180221 = matrix(0, nregions_180221, nregions_180221);
colnames(conn_180221) <- rownames(conn_180221) <- region_names_180221;
cx = 1;
for (i in region_names_180221){
  idx = which(regions_180221 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180221){
    idy = which(regions_180221 == j, arr.ind = FALSE);
    conn_180221[cx,cy] = sum(adj_pearson_180221[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_o_180221 = conn_180221;
adj_pearson_180221 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_180221_pycharm.csv", header=FALSE)
adj_pearson_180221 <- as.matrix(adj_pearson_180221);
colnames(adj_pearson_180221) <- regions_180221;
rownames(adj_pearson_180221) <- regions_180221;
plot(adj_pearson_180221)
net <- network(adj_pearson_180221, loop=TRUE)
#plot(net)

region_names_180221 <- unique(regions_180221);
nregions_180221 = length(region_names_180221);
conn_180221 = matrix(0, nregions_180221, nregions_180221);
colnames(conn_180221) <- rownames(conn_180221) <- region_names_180221;
cx = 1;
for (i in region_names_180221){
  idx = which(regions_180221 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180221){
    idy = which(regions_180221 == j, arr.ind = FALSE);
    conn_180221[cx,cy] = sum(adj_pearson_180221[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_e_180221 = conn_180221;

conn_diff_180221 = conn_e_180221 - conn_o_180221
net_180221 <- network(conn_diff_180221, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_180221.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180221, label = unique(rownames(conn_diff_180221)), edge.label = conn_diff_180221, label.cex = 2, edge.label.cex = 2, 
             edge.lwd = conn_diff_180221/80)
dev.off()




adj_pearson_180228 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_180228_pycharm.csv", header=FALSE)
adj_pearson_180228 <- as.matrix(adj_pearson_180228);
colnames(adj_pearson_180228) <- regions_180228;
rownames(adj_pearson_180228) <- regions_180228;
plot(adj_pearson_180228)
net <- network(adj_pearson_180228, loop=TRUE)
#plot(net)

region_names_180228 <- unique(regions_180228);
nregions_180228 = length(region_names_180228);
conn_180228 = matrix(0, nregions_180228, nregions_180228);
colnames(conn_180228) <- rownames(conn_180228) <- region_names_180228;
cx = 1;
for (i in region_names_180228){
  idx = which(regions_180228 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180228){
    idy = which(regions_180228 == j, arr.ind = FALSE);
    conn_180228[cx,cy] = sum(adj_pearson_180228[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_o_180228 = conn_180228;
adj_pearson_180228 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_180228_pycharm.csv", header=FALSE)
adj_pearson_180228 <- as.matrix(adj_pearson_180228);
colnames(adj_pearson_180228) <- regions_180228;
rownames(adj_pearson_180228) <- regions_180228;
plot(adj_pearson_180228)
net <- network(adj_pearson_180228, loop=TRUE)
#plot(net)

region_names_180228 <- unique(regions_180228);
nregions_180228 = length(region_names_180228);
conn_180228 = matrix(0, nregions_180228, nregions_180228);
colnames(conn_180228) <- rownames(conn_180228) <- region_names_180228;
cx = 1;
for (i in region_names_180228){
  idx = which(regions_180228 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180228){
    idy = which(regions_180228 == j, arr.ind = FALSE);
    conn_180228[cx,cy] = sum(adj_pearson_180228[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_e_180228 = conn_180228;

conn_diff_180228 = conn_e_180228 - conn_o_180228
net_180228 <- network(conn_diff_180228, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_180228.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180228, label = unique(rownames(conn_diff_180228)), edge.label = conn_diff_180228, label.cex = 2, edge.label.cex = 2
             , edge.lwd = conn_diff_180228/50)
dev.off()



adj_pearson_180228 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_180228_pycharm.csv", header=FALSE)
adj_pearson_180228 <- as.matrix(adj_pearson_180228);
colnames(adj_pearson_180228) <- regions_180228;
rownames(adj_pearson_180228) <- regions_180228;
plot(adj_pearson_180228)
net <- network(adj_pearson_180228, loop=TRUE)
#plot(net)

region_names_180228 <- unique(regions_180228);
nregions_180228 = length(region_names_180228);
conn_180228 = matrix(0, nregions_180228, nregions_180228);
colnames(conn_180228) <- rownames(conn_180228) <- region_names_180228;
cx = 1;
for (i in region_names_180228){
  idx = which(regions_180228 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180228){
    idy = which(regions_180228 == j, arr.ind = FALSE);
    conn_180228[cx,cy] = sum(adj_pearson_180228[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_o_180228 = conn_180228;
adj_pearson_180228 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_180228_pycharm.csv", header=FALSE)
adj_pearson_180228 <- as.matrix(adj_pearson_180228);
colnames(adj_pearson_180228) <- regions_180228;
rownames(adj_pearson_180228) <- regions_180228;
plot(adj_pearson_180228)
net <- network(adj_pearson_180228, loop=TRUE)
#plot(net)

region_names_180228 <- unique(regions_180228);
nregions_180228 = length(region_names_180228);
conn_180228 = matrix(0, nregions_180228, nregions_180228);
colnames(conn_180228) <- rownames(conn_180228) <- region_names_180228;
cx = 1;
for (i in region_names_180228){
  idx = which(regions_180228 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180228){
    idy = which(regions_180228 == j, arr.ind = FALSE);
    conn_180228[cx,cy] = sum(adj_pearson_180228[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_e_180228 = conn_180228;

conn_diff_180228 = conn_e_180228 - conn_o_180228
net_180228 <- network(conn_diff_180228, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_180228.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180228, label = unique(rownames(conn_diff_180228)), edge.label = conn_diff_180228, label.cex = 2, edge.label.cex = 2
             , edge.lwd = conn_diff_180228/50)
dev.off()


adj_pearson_180302 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_180302_pycharm.csv", header=FALSE)
adj_pearson_180302 <- as.matrix(adj_pearson_180302);
colnames(adj_pearson_180302) <- regions_180302;
rownames(adj_pearson_180302) <- regions_180302;
plot(adj_pearson_180302)
net <- network(adj_pearson_180302, loop=TRUE)
#plot(net)

region_names_180302 <- unique(regions_180302);
nregions_180302 = length(region_names_180302);
conn_180302 = matrix(0, nregions_180302, nregions_180302);
colnames(conn_180302) <- rownames(conn_180302) <- region_names_180302;
cx = 1;
for (i in region_names_180302){
  idx = which(regions_180302 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180302){
    idy = which(regions_180302 == j, arr.ind = FALSE);
    conn_180302[cx,cy] = sum(adj_pearson_180302[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_o_180302 = conn_180302;
adj_pearson_180302 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_180302_pycharm.csv", header=FALSE)
adj_pearson_180302 <- as.matrix(adj_pearson_180302);
colnames(adj_pearson_180302) <- regions_180302;
rownames(adj_pearson_180302) <- regions_180302;
plot(adj_pearson_180302)
net <- network(adj_pearson_180302, loop=TRUE)
#plot(net)

region_names_180302 <- unique(regions_180302);
nregions_180302 = length(region_names_180302);
conn_180302 = matrix(0, nregions_180302, nregions_180302);
colnames(conn_180302) <- rownames(conn_180302) <- region_names_180302;
cx = 1;
for (i in region_names_180302){
  idx = which(regions_180302 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180302){
    idy = which(regions_180302 == j, arr.ind = FALSE);
    conn_180302[cx,cy] = sum(adj_pearson_180302[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_e_180302 = conn_180302;

conn_diff_180302 = conn_e_180302 - conn_o_180302
net_180302 <- network(conn_diff_180302, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_180302.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180302, label = unique(rownames(conn_diff_180302)), edge.label = conn_diff_180302, label.cex = 2, edge.label.cex = 2
             , edge.lwd = conn_diff_180302/30)
dev.off()



adj_pearson_180419 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_180419_pycharm.csv", header=FALSE)
adj_pearson_180419 <- as.matrix(adj_pearson_180419);
colnames(adj_pearson_180419) <- regions_180419;
rownames(adj_pearson_180419) <- regions_180419;
plot(adj_pearson_180419)
net <- network(adj_pearson_180419, loop=TRUE)
#plot(net)

region_names_180419 <- unique(regions_180419);
nregions_180419 = length(region_names_180419);
conn_180419 = matrix(0, nregions_180419, nregions_180419);
colnames(conn_180419) <- rownames(conn_180419) <- region_names_180419;
cx = 1;
for (i in region_names_180419){
  idx = which(regions_180419 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180419){
    idy = which(regions_180419 == j, arr.ind = FALSE);
    conn_180419[cx,cy] = sum(adj_pearson_180419[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_o_180419 = conn_180419;
adj_pearson_180419 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_180419_pycharm.csv", header=FALSE)
adj_pearson_180419 <- as.matrix(adj_pearson_180419);
colnames(adj_pearson_180419) <- regions_180419;
rownames(adj_pearson_180419) <- regions_180419;
plot(adj_pearson_180419)
net <- network(adj_pearson_180419, loop=TRUE)
#plot(net)

region_names_180419 <- unique(regions_180419);
nregions_180419 = length(region_names_180419);
conn_180419 = matrix(0, nregions_180419, nregions_180419);
colnames(conn_180419) <- rownames(conn_180419) <- region_names_180419;
cx = 1;
for (i in region_names_180419){
  idx = which(regions_180419 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180419){
    idy = which(regions_180419 == j, arr.ind = FALSE);
    conn_180419[cx,cy] = sum(adj_pearson_180419[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_e_180419 = conn_180419;

conn_diff_180419 = conn_e_180419 - conn_o_180419
net_180419 <- network(conn_diff_180419, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_180419.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180419, label = unique(rownames(conn_diff_180419)), edge.label = conn_diff_180419, label.cex = 2, edge.label.cex = 2, 
             edge.lwd = conn_diff_180419/50)
dev.off()


adj_pearson_180420 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_180420_pycharm.csv", header=FALSE)
adj_pearson_180420 <- as.matrix(adj_pearson_180420);
colnames(adj_pearson_180420) <- regions_180420;
rownames(adj_pearson_180420) <- regions_180420;
plot(adj_pearson_180420)
net <- network(adj_pearson_180420, loop=TRUE)
#plot(net)

region_names_180420 <- unique(regions_180420);
nregions_180420 = length(region_names_180420);
conn_180420 = matrix(0, nregions_180420, nregions_180420);
colnames(conn_180420) <- rownames(conn_180420) <- region_names_180420;
cx = 1;
for (i in region_names_180420){
  idx = which(regions_180420 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180420){
    idy = which(regions_180420 == j, arr.ind = FALSE);
    conn_180420[cx,cy] = sum(adj_pearson_180420[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_o_180420 = conn_180420;
adj_pearson_180420 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_180420_pycharm.csv", header=FALSE)
adj_pearson_180420 <- as.matrix(adj_pearson_180420);
colnames(adj_pearson_180420) <- regions_180420;
rownames(adj_pearson_180420) <- regions_180420;
plot(adj_pearson_180420)
net <- network(adj_pearson_180420, loop=TRUE)
#plot(net)

region_names_180420 <- unique(regions_180420);
nregions_180420 = length(region_names_180420);
conn_180420 = matrix(0, nregions_180420, nregions_180420);
colnames(conn_180420) <- rownames(conn_180420) <- region_names_180420;
cx = 1;
for (i in region_names_180420){
  idx = which(regions_180420 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180420){
    idy = which(regions_180420 == j, arr.ind = FALSE);
    conn_180420[cx,cy] = sum(adj_pearson_180420[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_e_180420 = conn_180420;

conn_diff_180420 = conn_e_180420 - conn_o_180420
net_180420 <- network(conn_diff_180420, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_180420.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180420, label = unique(rownames(conn_diff_180420)), edge.label = conn_diff_180420, 
             label.cex = 2, edge.label.cex = 2, edge.lwd = conn_diff_180420/100)
dev.off()



adj_pearson_180423 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/prestimulus_pycharm/result_GLMCC_prestimulus_180423_pycharm.csv", header=FALSE)
adj_pearson_180423 <- as.matrix(adj_pearson_180423);
colnames(adj_pearson_180423) <- regions_180423;
rownames(adj_pearson_180423) <- regions_180423;
plot(adj_pearson_180423)
net <- network(adj_pearson_180423, loop=TRUE)
#plot(net)

region_names_180423 <- unique(regions_180423);
nregions_180423 = length(region_names_180423);
conn_180423 = matrix(0, nregions_180423, nregions_180423);
colnames(conn_180423) <- rownames(conn_180423) <- region_names_180423;
cx = 1;
for (i in region_names_180423){
  idx = which(regions_180423 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180423){
    idy = which(regions_180423 == j, arr.ind = FALSE);
    conn_180423[cx,cy] = sum(adj_pearson_180423[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_o_180423 = conn_180423;
adj_pearson_180423 <- read.csv("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/lightON_pycharm/result_GLMCC_lightON_180423_pycharm.csv", header=FALSE)
adj_pearson_180423 <- as.matrix(adj_pearson_180423);
colnames(adj_pearson_180423) <- regions_180423;
rownames(adj_pearson_180423) <- regions_180423;
plot(adj_pearson_180423)
net <- network(adj_pearson_180423, loop=TRUE)
#plot(net)

region_names_180423 <- unique(regions_180423);
nregions_180423 = length(region_names_180423);
conn_180423 = matrix(0, nregions_180423, nregions_180423);
colnames(conn_180423) <- rownames(conn_180423) <- region_names_180423;
cx = 1;
for (i in region_names_180423){
  idx = which(regions_180423 == i, arr.ind = FALSE);
  cy = 1;
  for (j in region_names_180423){
    idy = which(regions_180423 == j, arr.ind = FALSE);
    conn_180423[cx,cy] = sum(adj_pearson_180423[idx,idy]);
    cy = cy + 1;
  }
  cx = cx + 1;
}
conn_e_180423 = conn_180423;

conn_diff_180423 = conn_e_180423 - conn_o_180423
net_180423 <- network(conn_diff_180423, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_180423.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180423, label = unique(rownames(conn_diff_180423)), #edge.label = conn_diff_180423[conn_diff_180423>100], 
             label.cex = 2, edge.label.cex = 2, edge.lwd = conn_diff_180423/50, usecurve=TRUE, arrowhead.cex = 1.2, vertex.cex = 1.5, pad =1)
dev.off()


write.csv(conn_diff_171019, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_171019.csv");
write.csv(conn_diff_171207, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_171207.csv");
write.csv(conn_diff_171208, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_171208.csv");
write.csv(conn_diff_171213, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_171213.csv");
write.csv(conn_diff_180110, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_180110.csv");
write.csv(conn_diff_180111, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_180111.csv");
write.csv(conn_diff_180131, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_180131.csv");
write.csv(conn_diff_180221, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_180221.csv");
write.csv(conn_diff_180228, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_180228.csv");
write.csv(conn_diff_180302, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_180302.csv");
write.csv(conn_diff_180419, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_180419.csv");
write.csv(conn_diff_180420, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_180420.csv");
write.csv(conn_diff_180423, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_180423.csv");


################################# All difference matrices summed up together ##########################
conn_diff_all = matrix(0, nrow = 9, ncol = 9)
lab = levels(regions_match_s2$Regions);
rownames(conn_diff_all) <- levels(regions_match_s2$Regions);
colnames(conn_diff_all) <- levels(regions_match_s2$Regions);
for (i in seq(1,nrow(conn_diff_171019)))
  {
  for (j in seq(1,nrow(conn_diff_171019))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_171019)[i])
      posy = which(lab == colnames(conn_diff_171019)[j])
      values = conn_diff_171019[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_171207)))
{
  for (j in seq(1,nrow(conn_diff_171207))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_171207)[i])
      posy = which(lab == colnames(conn_diff_171207)[j])
      values = conn_diff_171207[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_171208)))
{
  for (j in seq(1,nrow(conn_diff_171208))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_171208)[i])
      posy = which(lab == colnames(conn_diff_171208)[j])
      values = conn_diff_171208[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_171213)))
{
  for (j in seq(1,nrow(conn_diff_171213))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_171213)[i])
      posy = which(lab == colnames(conn_diff_171213)[j])
      values = conn_diff_171213[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_180110)))
{
  for (j in seq(1,nrow(conn_diff_180110))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_180110)[i])
      posy = which(lab == colnames(conn_diff_180110)[j])
      values = conn_diff_180110[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_180111)))
{
  for (j in seq(1,nrow(conn_diff_180111))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_180111)[i])
      posy = which(lab == colnames(conn_diff_180111)[j])
      values = conn_diff_180111[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_180131)))
{
  for (j in seq(1,nrow(conn_diff_180131))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_180131)[i])
      posy = which(lab == colnames(conn_diff_180131)[j])
      values = conn_diff_180131[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_180221)))
{
  for (j in seq(1,nrow(conn_diff_180221))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_180221)[i])
      posy = which(lab == colnames(conn_diff_180221)[j])
      values = conn_diff_180221[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_180228)))
{
  for (j in seq(1,nrow(conn_diff_180228))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_180228)[i])
      posy = which(lab == colnames(conn_diff_180228)[j])
      values = conn_diff_180228[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_180302)))
{
  for (j in seq(1,nrow(conn_diff_180302))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_180302)[i])
      posy = which(lab == colnames(conn_diff_180302)[j])
      values = conn_diff_180302[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_180419)))
{
  for (j in seq(1,nrow(conn_diff_180419))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_180419)[i])
      posy = which(lab == colnames(conn_diff_180419)[j])
      values = conn_diff_180419[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_180420)))
{
  for (j in seq(1,nrow(conn_diff_180420))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_180420)[i])
      posy = which(lab == colnames(conn_diff_180420)[j])
      values = conn_diff_180420[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_180423)))
{
  for (j in seq(1,nrow(conn_diff_180423))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_180423)[i])
      posy = which(lab == colnames(conn_diff_180423)[j])
      values = conn_diff_180423[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
net_all <- network(conn_diff_all, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_all.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_all, label = unique(rownames(conn_diff_all)), #edge.label = conn_diff_180423[conn_diff_180423>100], 
             label.cex = 2, edge.label.cex = 2, edge.lwd = conn_diff_all/100, usecurve=TRUE, arrowhead.cex = 1.2, vertex.cex = 1.5, pad =1)
dev.off()
write.csv(conn_diff_all, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_all.csv");



###### JUST for the DAY ###############################################################################################
conn_diff_all = matrix(0, nrow = 9, ncol = 9)
lab = levels(regions_match_s2$Regions);
rownames(conn_diff_all) <- levels(regions_match_s2$Regions);
colnames(conn_diff_all) <- levels(regions_match_s2$Regions);
for (i in seq(1,nrow(conn_diff_171019)))
{
  for (j in seq(1,nrow(conn_diff_171019))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_171019)[i])
      posy = which(lab == colnames(conn_diff_171019)[j])
      values = conn_diff_171019[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_180131)))
{
  for (j in seq(1,nrow(conn_diff_180131))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_180131)[i])
      posy = which(lab == colnames(conn_diff_180131)[j])
      values = conn_diff_180131[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_180302)))
{
  for (j in seq(1,nrow(conn_diff_180302))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_180302)[i])
      posy = which(lab == colnames(conn_diff_180302)[j])
      values = conn_diff_180302[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_180419)))
{
  for (j in seq(1,nrow(conn_diff_180419))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_180419)[i])
      posy = which(lab == colnames(conn_diff_180419)[j])
      values = conn_diff_180419[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_180420)))
{
  for (j in seq(1,nrow(conn_diff_180420))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_180420)[i])
      posy = which(lab == colnames(conn_diff_180420)[j])
      values = conn_diff_180420[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_180423)))
{
  for (j in seq(1,nrow(conn_diff_180423))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_180423)[i])
      posy = which(lab == colnames(conn_diff_180423)[j])
      values = conn_diff_180423[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
net_all_day_reduced <- network(conn_diff_all, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_all_day_reduced_circ.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_all_day_reduced, label = unique(rownames(conn_diff_all)), #edge.label = conn_diff_180423[conn_diff_180423>100], 
             label.cex = 2, edge.label.cex = 2, edge.lwd = conn_diff_all/100, usecurve=TRUE, arrowhead.cex = 1.2, vertex.cex = 1.5, pad =1, 
             mode = "circle", displayisolates = FALSE)
dev.off()
write.csv(conn_diff_all, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_all_day_reduced.csv");



###### JUST for the NIGHT ###############################################################################################
conn_diff_all = matrix(0, nrow = 9, ncol = 9)
lab = levels(regions_match_s2$Regions);
rownames(conn_diff_all) <- levels(regions_match_s2$Regions);
colnames(conn_diff_all) <- levels(regions_match_s2$Regions);
for (i in seq(1,nrow(conn_diff_171207)))
{
  for (j in seq(1,nrow(conn_diff_171207))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_171207)[i])
      posy = which(lab == colnames(conn_diff_171207)[j])
      values = conn_diff_171207[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_171208)))
{
  for (j in seq(1,nrow(conn_diff_171208))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_171208)[i])
      posy = which(lab == colnames(conn_diff_171208)[j])
      values = conn_diff_171208[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_171213)))
{
  for (j in seq(1,nrow(conn_diff_171213))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_171213)[i])
      posy = which(lab == colnames(conn_diff_171213)[j])
      values = conn_diff_171213[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_180110)))
{
  for (j in seq(1,nrow(conn_diff_180110))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_180110)[i])
      posy = which(lab == colnames(conn_diff_180110)[j])
      values = conn_diff_180110[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_180111)))
{
  for (j in seq(1,nrow(conn_diff_180111))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_180111)[i])
      posy = which(lab == colnames(conn_diff_180111)[j])
      values = conn_diff_180111[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_180221)))
{
  for (j in seq(1,nrow(conn_diff_180221))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_180221)[i])
      posy = which(lab == colnames(conn_diff_180221)[j])
      values = conn_diff_180221[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_diff_180228)))
{
  for (j in seq(1,nrow(conn_diff_180228))){
    if (i != j){
      posx = which(lab == rownames(conn_diff_180228)[i])
      posy = which(lab == colnames(conn_diff_180228)[j])
      values = conn_diff_180228[i,j]
      conn_diff_all[posx,posy] = conn_diff_all[posx,posy] + values
    }
  }
}
conn_diff_all[which(rownames(conn_diff_all) == "zona incerta"),] <- rep(0,9,1)
conn_diff_all[,which(rownames(conn_diff_all) == "zona incerta")] <- rep(0,9,1)
conn_diff_all[,which(rownames(conn_diff_all) == "mammillary complex")] <- rep(0,9,1)
conn_diff_all[which(rownames(conn_diff_all) == "mammillary complex"),] <- rep(0,9,1)
conn_diff_all[,which(rownames(conn_diff_all) == "anterior hypothalamus")] <- rep(0,9,1)
conn_diff_all[which(rownames(conn_diff_all) == "anterior hypothalamus"),] <- rep(0,9,1)
net_all_night <- network(conn_diff_all, directed = TRUE, matrix.type = "adjacency")
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_all_night_reduced_circ.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_all_night, label = unique(rownames(conn_diff_all)), #edge.label = conn_diff_180423[conn_diff_180423>100], 
             label.cex = 2, edge.label.cex = 2, edge.lwd = conn_diff_all/100, usecurve=TRUE, arrowhead.cex = 1.2, vertex.cex = 1.5, pad = 1,
             mode = "circle", displayisolates = FALSE)
dev.off()
write.csv(conn_diff_all, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_diff_all_night_reduced.csv");

###############################################################################################################################################
########################### Compute all ongoing in DAY experiments ####################################################################
conn_ongoing_day = matrix(0, nrow = 9, ncol = 9)
lab = levels(regions_match_s2$Regions);
rownames(conn_ongoing_day) <- levels(regions_match_s2$Regions);
colnames(conn_ongoing_day) <- levels(regions_match_s2$Regions);
for (i in seq(1,nrow(conn_o_171019)))
{
  for (j in seq(1,nrow(conn_o_171019))){
    if (i != j){
      posx = which(lab == rownames(conn_o_171019)[i])
      posy = which(lab == colnames(conn_o_171019)[j])
      values = conn_o_171019[i,j]
      conn_ongoing_day[posx,posy] = conn_ongoing_day[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_o_180131)))
{
  for (j in seq(1,nrow(conn_o_180131))){
    if (i != j){
      posx = which(lab == rownames(conn_o_180131)[i])
      posy = which(lab == colnames(conn_o_180131)[j])
      values = conn_o_180131[i,j]
      conn_ongoing_day[posx,posy] = conn_ongoing_day[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_o_180302)))
{
  for (j in seq(1,nrow(conn_o_180302))){
    if (i != j){
      posx = which(lab == rownames(conn_o_180302)[i])
      posy = which(lab == colnames(conn_o_180302)[j])
      values = conn_o_180302[i,j]
      conn_ongoing_day[posx,posy] = conn_ongoing_day[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_o_180419)))
{
  for (j in seq(1,nrow(conn_o_180419))){
    if (i != j){
      posx = which(lab == rownames(conn_o_180419)[i])
      posy = which(lab == colnames(conn_o_180419)[j])
      values = conn_o_180419[i,j]
      conn_ongoing_day[posx,posy] = conn_ongoing_day[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_o_180420)))
{
  for (j in seq(1,nrow(conn_o_180420))){
    if (i != j){
      posx = which(lab == rownames(conn_o_180420)[i])
      posy = which(lab == colnames(conn_o_180420)[j])
      values = conn_o_180420[i,j]
      conn_ongoing_day[posx,posy] = conn_ongoing_day[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_o_180423)))
{
  for (j in seq(1,nrow(conn_o_180423))){
    if (i != j){
      posx = which(lab == rownames(conn_o_180423)[i])
      posy = which(lab == colnames(conn_o_180423)[j])
      values = conn_o_180423[i,j]
      conn_ongoing_day[posx,posy] = conn_ongoing_day[posx,posy] + values
    }
  }
}
net_all_ongoing_day <- network(conn_ongoing_day, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_all_day_reduced_circ.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_all_ongoing_day, label = unique(rownames(conn_ongoing_day)), #edge.label = conn_diff_180423[conn_diff_180423>100], 
             label.cex = 2, edge.label.cex = 2, edge.lwd = conn_ongoing_day/100, usecurve=TRUE, arrowhead.cex = 1.2, vertex.cex = 1.5, pad =1, 
             mode = "circle", displayisolates = FALSE)
dev.off()
write.csv(conn_ongoing_day, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_all_day.csv");

###############################################################################################################################################
########################### Compute all ongoing in NIGHT experiments ####################################################################
conn_ongoing_night = matrix(0, nrow = 9, ncol = 9)
lab = levels(regions_match_s2$Regions);
rownames(conn_ongoing_night) <- levels(regions_match_s2$Regions);
colnames(conn_ongoing_night) <- levels(regions_match_s2$Regions);
for (i in seq(1,nrow(conn_o_171207)))
{
  for (j in seq(1,nrow(conn_o_171207))){
    if (i != j){
      posx = which(lab == rownames(conn_o_171207)[i])
      posy = which(lab == colnames(conn_o_171207)[j])
      values = conn_o_171207[i,j]
      conn_ongoing_night[posx,posy] = conn_ongoing_night[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_o_171208)))
{
  for (j in seq(1,nrow(conn_o_171208))){
    if (i != j){
      posx = which(lab == rownames(conn_o_171208)[i])
      posy = which(lab == colnames(conn_o_171208)[j])
      values = conn_o_171208[i,j]
      conn_ongoing_night[posx,posy] = conn_ongoing_night[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_o_171213)))
{
  for (j in seq(1,nrow(conn_o_171213))){
    if (i != j){
      posx = which(lab == rownames(conn_o_171213)[i])
      posy = which(lab == colnames(conn_o_171213)[j])
      values = conn_o_171213[i,j]
      conn_ongoing_night[posx,posy] = conn_ongoing_night[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_o_180110)))
{
  for (j in seq(1,nrow(conn_o_180110))){
    if (i != j){
      posx = which(lab == rownames(conn_o_180110)[i])
      posy = which(lab == colnames(conn_o_180110)[j])
      values = conn_o_180110[i,j]
      conn_ongoing_night[posx,posy] = conn_ongoing_night[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_o_180111)))
{
  for (j in seq(1,nrow(conn_o_180111))){
    if (i != j){
      posx = which(lab == rownames(conn_o_180111)[i])
      posy = which(lab == colnames(conn_o_180111)[j])
      values = conn_o_180111[i,j]
      conn_ongoing_night[posx,posy] = conn_ongoing_night[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_o_180221)))
{
  for (j in seq(1,nrow(conn_o_180221))){
    if (i != j){
      posx = which(lab == rownames(conn_o_180221)[i])
      posy = which(lab == colnames(conn_o_180221)[j])
      values = conn_o_180221[i,j]
      conn_ongoing_night[posx,posy] = conn_ongoing_night[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_o_180228)))
{
  for (j in seq(1,nrow(conn_o_180228))){
    if (i != j){
      posx = which(lab == rownames(conn_o_180228)[i])
      posy = which(lab == colnames(conn_o_180228)[j])
      values = conn_o_180228[i,j]
      conn_ongoing_night[posx,posy] = conn_ongoing_night[posx,posy] + values
    }
  }
}
net_ongoing_night <- network(conn_ongoing_night, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_all_night_reduced_circ.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_ongoing_night, label = unique(rownames(conn_ongoing_night)), #edge.label = conn_diff_180423[conn_diff_180423>100], 
             label.cex = 2, edge.label.cex = 2, edge.lwd = conn_ongoing_night/100, usecurve=TRUE, arrowhead.cex = 1.2, vertex.cex = 1.5, pad =1, 
             mode = "circle", displayisolates = FALSE)
dev.off()
write.csv(conn_ongoing_night, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_all_night.csv");


###############################################################################################################################################
########################### Compute all evoked in DAY experiments ####################################################################
conn_evoked_day = matrix(0, nrow = 9, ncol = 9)
lab = levels(regions_match_s2$Regions);
rownames(conn_evoked_day) <- levels(regions_match_s2$Regions);
colnames(conn_evoked_day) <- levels(regions_match_s2$Regions);
for (i in seq(1,nrow(conn_e_171019)))
{
  for (j in seq(1,nrow(conn_e_171019))){
    if (i != j){
      posx = which(lab == rownames(conn_e_171019)[i])
      posy = which(lab == colnames(conn_e_171019)[j])
      values = conn_e_171019[i,j]
      conn_evoked_day[posx,posy] = conn_evoked_day[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_e_180131)))
{
  for (j in seq(1,nrow(conn_e_180131))){
    if (i != j){
      posx = which(lab == rownames(conn_e_180131)[i])
      posy = which(lab == colnames(conn_e_180131)[j])
      values = conn_e_180131[i,j]
      conn_evoked_day[posx,posy] = conn_evoked_day[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_e_180302)))
{
  for (j in seq(1,nrow(conn_e_180302))){
    if (i != j){
      posx = which(lab == rownames(conn_e_180302)[i])
      posy = which(lab == colnames(conn_e_180302)[j])
      values = conn_e_180302[i,j]
      conn_evoked_day[posx,posy] = conn_evoked_day[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_e_180419)))
{
  for (j in seq(1,nrow(conn_e_180419))){
    if (i != j){
      posx = which(lab == rownames(conn_e_180419)[i])
      posy = which(lab == colnames(conn_e_180419)[j])
      values = conn_e_180419[i,j]
      conn_evoked_day[posx,posy] = conn_evoked_day[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_e_180420)))
{
  for (j in seq(1,nrow(conn_e_180420))){
    if (i != j){
      posx = which(lab == rownames(conn_e_180420)[i])
      posy = which(lab == colnames(conn_e_180420)[j])
      values = conn_e_180420[i,j]
      conn_evoked_day[posx,posy] = conn_evoked_day[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_e_180423)))
{
  for (j in seq(1,nrow(conn_e_180423))){
    if (i != j){
      posx = which(lab == rownames(conn_e_180423)[i])
      posy = which(lab == colnames(conn_e_180423)[j])
      values = conn_e_180423[i,j]
      conn_evoked_day[posx,posy] = conn_evoked_day[posx,posy] + values
    }
  }
}
net_all_evoked_day <- network(conn_evoked_day, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_all_night_reduced_circ.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_all_evoked_day, label = unique(rownames(conn_evoked_day)), #edge.label = conn_diff_180423[conn_diff_180423>100], 
             label.cex = 2, edge.label.cex = 2, edge.lwd = conn_evoked_day/100, usecurve=TRUE, arrowhead.cex = 1.2, vertex.cex = 1.5, pad =1, 
             mode = "circle", displayisolates = FALSE)
dev.off()
write.csv(conn_evoked_day, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_all_day.csv");


###############################################################################################################################################
########################### Compute all evoked in NIGHT experiments ####################################################################
conn_evoked_night = matrix(0, nrow = 9, ncol = 9)
lab = levels(regions_match_s2$Regions);
rownames(conn_evoked_night) <- levels(regions_match_s2$Regions);
colnames(conn_evoked_night) <- levels(regions_match_s2$Regions);
for (i in seq(1,nrow(conn_e_171207)))
{
  for (j in seq(1,nrow(conn_e_171207))){
    if (i != j){
      posx = which(lab == rownames(conn_e_171207)[i])
      posy = which(lab == colnames(conn_e_171207)[j])
      values = conn_e_171207[i,j]
      conn_evoked_night[posx,posy] = conn_evoked_night[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_e_171208)))
{
  for (j in seq(1,nrow(conn_e_171208))){
    if (i != j){
      posx = which(lab == rownames(conn_e_171208)[i])
      posy = which(lab == colnames(conn_e_171208)[j])
      values = conn_e_171208[i,j]
      conn_evoked_night[posx,posy] = conn_evoked_night[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_e_171213)))
{
  for (j in seq(1,nrow(conn_e_171213))){
    if (i != j){
      posx = which(lab == rownames(conn_e_171213)[i])
      posy = which(lab == colnames(conn_e_171213)[j])
      values = conn_e_171213[i,j]
      conn_evoked_night[posx,posy] = conn_evoked_night[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_e_180110)))
{
  for (j in seq(1,nrow(conn_e_180110))){
    if (i != j){
      posx = which(lab == rownames(conn_e_180110)[i])
      posy = which(lab == colnames(conn_e_180110)[j])
      values = conn_e_180110[i,j]
      conn_evoked_night[posx,posy] = conn_evoked_night[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_e_180111)))
{
  for (j in seq(1,nrow(conn_e_180111))){
    if (i != j){
      posx = which(lab == rownames(conn_e_180111)[i])
      posy = which(lab == colnames(conn_e_180111)[j])
      values = conn_e_180111[i,j]
      conn_evoked_night[posx,posy] = conn_evoked_night[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_e_180221)))
{
  for (j in seq(1,nrow(conn_e_180221))){
    if (i != j){
      posx = which(lab == rownames(conn_e_180221)[i])
      posy = which(lab == colnames(conn_e_180221)[j])
      values = conn_e_180221[i,j]
      conn_evoked_night[posx,posy] = conn_evoked_night[posx,posy] + values
    }
  }
}
for (i in seq(1,nrow(conn_e_180228)))
{
  for (j in seq(1,nrow(conn_e_180228))){
    if (i != j){
      posx = which(lab == rownames(conn_e_180228)[i])
      posy = which(lab == colnames(conn_e_180228)[j])
      values = conn_e_180228[i,j]
      conn_evoked_night[posx,posy] = conn_evoked_night[posx,posy] + values
    }
  }
}
net_evoked_night <- network(conn_evoked_night, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_all_night_reduced_circ.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_evoked_night, label = unique(rownames(conn_evoked_night)), #edge.label = conn_diff_180423[conn_diff_180423>100], 
             label.cex = 2, edge.label.cex = 2, edge.lwd = conn_evoked_night/100, usecurve=TRUE, arrowhead.cex = 1.2, vertex.cex = 1.5, pad =1, 
             mode = "circle", displayisolates = FALSE)
dev.off()
write.csv(conn_evoked_night, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_all_night.csv");

########################CONTRASTS MATRICES/GRAPHS####################
conn_ongoing_day_vs_night = conn_ongoing_day - conn_ongoing_night
net_ongoing_night_vs_day <- network(conn_ongoing_day_vs_night, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoign_night_vs_day_circ.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_ongoing_night_vs_day, label = unique(rownames(conn_ongoing_day_vs_night)), #edge.label = conn_diff_180423[conn_diff_180423>100], 
             label.cex = 2, edge.label.cex = 2, edge.lwd = conn_ongoing_day_vs_night/100, usecurve=TRUE, arrowhead.cex = 1.2, vertex.cex = 1.5, pad =1, 
             mode = "circle", displayisolates = FALSE)
dev.off()
write.csv(conn_ongoing_day_vs_night, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_ongoing_night_vs_day.csv");

conn_evoked_day_vs_night = conn_evoked_day - conn_evoked_night
net_evoked_night_vs_day <- network(conn_evoked_day_vs_night, directed = TRUE)
png(file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_night_vs_day_circ.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_evoked_night_vs_day, label = unique(rownames(conn_evoked_day_vs_night)), #edge.label = conn_diff_180423[conn_diff_180423>100], 
             label.cex = 2, edge.label.cex = 2, edge.lwd = conn_tmp/100, usecurve=TRUE, arrowhead.cex = 1.2, vertex.cex = 1.5, pad =1, 
             mode = "circle", displayisolates = FALSE)
dev.off()
write.csv(conn_evoked_day_vs_night, file = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/conn_evoked_night_vs_day.csv");

regions_reduced = colnames(conn_evoked_day_vs_night);
### Remove i) Anterior Hypothalamus, ii) Mammillary Complex and iii) Zona Incerta
regions_reduced = regions_reduced[c(-1, -4, -9),]



#########################################################################
df_171207 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171207_fr")
regions_171207 <- character(nrow(df_171207));
regions_171207[c(1,2,3, 12,13,28, 18,23,24,25,26,27,31)] = "ventromedial thalamus"
regions_171207[c(4,5,6,7,11,14,15,16,29,30)] = "paraventricular hypothalamus"
regions_171207[c(8,9,10,17)] = "anterior hypothalamus"
regions_171207[c(19,20,21,22)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_171207),2);
condition <- c(rep("Ongoing",nrow(df_171207)), rep("LightON",nrow(df_171207)));
ei_171207 <- array("Excitatory", length(regions_171207));
ei_171207[9] <- "Inhibitory";
dat_171207 <- data.frame(c(df_171207$Ongoing, df_171207$LightON), condition, nodeid, rep(regions_171207,2), ei_171207);
colnames(dat_171207) <- c("FiringRate", "condition", "nodeid", "region", "EI");
dat_171207$condition <- factor(dat_171207$condition, levels = c("Ongoing", "LightON"))
dat_171207$region <- factor(dat_171207$region)
dat_171207$FiringRate <- as.numeric(sub(",", ".", c(df_171207$Ongoing, df_171207$LightON), fixed = TRUE))
dat_171207$EI <- factor(ei_171207);
model_171207 <- lm(formula = "FiringRate ~ condition * region * EI + (1|nodeid)", data=dat_171207)
res_aov_171207 <- aov(model_171207)
summary(res_aov_171207)
cm=compare_means(FiringRate ~ condition, data = dat_171207, group.by = "region")

tmpdf_fr = aggregate(FiringRate ~ region + condition + EI, data = dat_171207, mean);
#tmpdf_fr$expID = rep("171207",1,length(levels(as.factor(dat_171207)))*2);
tmpdf_fr$expID = rep("171207",1,length(dat_171207)*2);
df_fr = tmpdf_fr;
#df_degree = rbind(df_degree, tmpdf_degree);

###
df_171019 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171019_fr")
regions_171019 <- character(nrow(df_171019));
regions_171019[seq(1,114)] = "ventromedial thalamus"
regions_171019[c(115, 116)] = "paraventricular hypothalamus"
regions_171019[seq(117,126)] = "ventromedial thalamus"
nodeid <- rep(1:nrow(df_171019),2);
condition <- c(rep("Ongoing",nrow(df_171019)), rep("LightON",nrow(df_171019)));
ei_171019 <- array("Excitatory", length(regions_171019));
ei_171019[c(60,125)] <- "Inhibitory";
dat_171019 <- data.frame(c(df_171019$Ongoing, df_171019$LightON), condition, nodeid, rep(regions_171019,2), ei_171019);
colnames(dat_171019) <- c("FiringRate", "condition", "nodeid", "region", "EI")
dat_171019$condition <- factor(dat_171019$condition, levels = c("Ongoing", "LightON"))
dat_171019$region <- factor(dat_171019$region)
dat_171019$FiringRate <- as.numeric(sub(",", ".", c(df_171019$Ongoing, df_171019$LightON), fixed = TRUE))
dat_171019$EI <- factor(ei_171019);
model_171019 <- lm(formula = "FiringRate ~ condition * region * EI + (1|nodeid)", data=dat_171019)
res_aov_171019 <- aov(model_171019)
summary(res_aov_171019)
cm=compare_means(FiringRate ~ condition, data = dat_171019, group.by = "region")

tmpdf_fr = aggregate(FiringRate ~ region + condition + EI, data = dat_171019, mean);
#tmpdf_fr$expID = rep("171019",1,length(levels(as.factor(tmpdf_fr$region)))*2);
tmpdf_fr$expID = rep("171019",1,nrow(tmpdf_fr));
df_fr = rbind(df_fr, tmpdf_fr);


###
df_171208 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171208_fr")
regions_171208 <- character(nrow(df_171208));
regions_171208[c(1,2,3,4,5,7,20,21,22,23,24,25,28)] = "ventromedial thalamus"
regions_171208[c(6, 8,9,11,13,16,17)] = "ventromedial hypothalamus"
regions_171208[c(10,12, 32)] = "paraventricular hypothalamus"
regions_171208[c(14,15,18,19)] = "dorsomedial hypothalamus"
regions_171208[c(26, 27, 29, 30, 31,33, 34)] = "anterior hypothalamus"
nodeid <- rep(1:nrow(df_171208),2);
condition <- c(rep("Ongoing",nrow(df_171208)), rep("LightON",nrow(df_171208)));
ei_171208 <- array("Excitatory", length(regions_171208));
dat_171208 <- data.frame(c(df_171208$Ongoing, df_171208$LightON), condition, nodeid, rep(regions_171208,2), ei_171208);
colnames(dat_171208) <- c("FiringRate", "condition", "nodeid", "region", "EI")
dat_171208$condition <- factor(dat_171208$condition, levels = c("Ongoing", "LightON"))
dat_171208$region <- factor(dat_171208$region)
dat_171208$FiringRate <- as.numeric(sub(",", ".", c(df_171208$Ongoing, df_171208$LightON), fixed = TRUE))
dat_171208$EI <- factor(ei_171208);
model_171208 <- lm(formula = "FiringRate ~ condition * region + (1|nodeid)", data=dat_171208)
res_aov_171208 <- aov(model_171208)
summary(res_aov_171208)
cm=compare_means(FiringRate ~ condition, data = dat_171208, group.by = "region")

tmpdf_fr = aggregate(FiringRate ~ region + condition + EI, data = dat_171208, mean);
#tmpdf_fr$expID = rep("171208",1,length(levels(as.factor(tmpdf_fr$region)))*2);
tmpdf_fr$expID = rep("171208",1,nrow(tmpdf_fr));
df_fr = rbind(df_fr, tmpdf_fr);



###
df_171213 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171213_fr")
regions_171213 <- character(nrow(df_171213));
regions_171213[c(11,12,13,20,21,22,23,24,25,26,27,29,30,31,32)] = "ventromedial thalamus"
regions_171213[c(1,2,3,4,5,6,7,8,9,10,14,15,16,17,18,19,28,33)] = "paraventricular hypothalamus"
nodeid <- rep(1:nrow(df_171213),2);
condition <- c(rep("Ongoing",nrow(df_171213)), rep("LightON",nrow(df_171213)));
ei_171213 <- array("Excitatory", length(regions_171213));
dat_171213 <- data.frame(c(df_171213$Ongoing, df_171213$LightON), condition, nodeid, rep(regions_171213,2), ei_171213);
colnames(dat_171213) <- c("FiringRate", "condition", "nodeid", "region", "EI")
dat_171213$condition <- factor(dat_171213$condition, levels = c("Ongoing", "LightON"))
dat_171213$region <- factor(dat_171213$region)
dat_171213$FiringRate <- as.numeric(sub(",", ".", c(df_171213$Ongoing, df_171213$LightON), fixed = TRUE))
dat_171213$EI <- factor(ei_171213);
model_171213 <- lm(formula = "FiringRate ~ condition * region + (1|nodeid)", data=dat_171213)
res_aov_171213 <- aov(model_171213)
summary(res_aov_171213)
cm=compare_means(FiringRate ~ condition, data = dat_171213, group.by = "region")

tmpdf_fr = aggregate(FiringRate ~ region + condition + EI, data = dat_171213, mean);
#tmpdf_fr$expID = rep("171213",1,length(levels(as.factor(tmpdf_fr$region)))*2);
tmpdf_fr$expID = rep("171213",1,nrow(tmpdf_fr));
df_fr = rbind(df_fr, tmpdf_fr);


###
df_180110 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180110_fr")
regions_180110 <- character(nrow(df_180110));
regions_180110[c(1,2,11,18,20,21,23)] = "ventromedial thalamus"
regions_180110[c(3,4,5,6,7,9,10, 19,22,24,25,26,27,28,29,30,31,32,33,34,35)] = "paraventricular hypothalamus"
regions_180110[c(8)] = "anterior hypothalamus"
regions_180110[c(12,13,14,15,16,17)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_180110),2);
condition <- c(rep("Ongoing",nrow(df_180110)), rep("LightON",nrow(df_180110)));
ei_180110 <- array("Excitatory", length(regions_180110));
ei_180110[16] <- "Inhibitory"
dat_180110 <- data.frame(c(df_180110$Ongoing, df_180110$LightON), condition, nodeid, rep(regions_180110,2), ei_180110);
colnames(dat_180110) <- c("FiringRate", "condition", "nodeid", "region", "EI")
dat_180110$condition <- factor(dat_180110$condition, levels = c("Ongoing", "LightON"))
dat_180110$region <- factor(dat_180110$region)
dat_180110$FiringRate <- as.numeric(sub(",", ".", c(df_180110$Ongoing, df_180110$LightON), fixed = TRUE))
dat_180110$EI <- factor(ei_180110)
model_180110 <- lm(formula = "FiringRate ~ condition * region * EI + (1|nodeid)", data=dat_180110)
res_aov_180110 <- aov(model_180110)
summary(res_aov_180110)
cm=compare_means(FiringRate ~ condition, data = dat_180110, group.by = "region")

tmpdf_fr = aggregate(FiringRate ~ region + condition + EI, data = dat_180110, mean);
#tmpdf_fr$expID = rep("180110",1,length(levels(as.factor(tmpdf_fr$region)))*2);
tmpdf_fr$expID = rep("180110",1,nrow(tmpdf_fr));
df_fr = rbind(df_fr, tmpdf_fr);

###
df_180111 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180111_fr")
regions_180111 <- character(nrow(df_180111));
regions_180111[c(1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,22,23,24)] = "ventromedial thalamus"
regions_180111[c(14)] = "paraventricular hypothalamus"
regions_180111[c(19,20,21)] = "ventromedial hypothalamus"
nodeid <- rep(1:nrow(df_180111),2);
condition <- c(rep("Ongoing",nrow(df_180111)), rep("LightON",nrow(df_180111)));
ei_180111 <- array("Excitatory", length(regions_180111));
dat_180111 <- data.frame(c(df_180111$Ongoing, df_180111$LightON), condition, nodeid, rep(regions_180111,2), ei_180111);
colnames(dat_180111) <- c("FiringRate", "condition", "nodeid", "region", "EI")
dat_180111$condition <- factor(dat_180111$condition, levels = c("Ongoing", "LightON"))
dat_180111$region <- factor(dat_180111$region)
dat_180111$FiringRate <- as.numeric(sub(",", ".", c(df_180111$Ongoing, df_180111$LightON), fixed = TRUE))
dat_180111$EI <- factor(ei_180111);
model_180111 <- lm(formula = "FiringRate ~ condition * region + (1|nodeid)", data=dat_180111)
res_aov_180111 <- aov(model_180111)
summary(res_aov_180111)
cm=compare_means(FiringRate ~ condition, data = dat_180111, group.by = "region")

tmpdf_fr = aggregate(FiringRate ~ region + condition + EI, data = dat_180111, mean);
#tmpdf_fr$expID = rep("180111",1,length(levels(as.factor(tmpdf_fr$region)))*2);
tmpdf_fr$expID = rep("180111",1,nrow(tmpdf_fr));
df_fr = rbind(df_fr, tmpdf_fr);


###
df_180131 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180131_fr")
regions_180131 <- character(nrow(df_180131));
regions_180131[c(1,2,4,5,6,14,15,16, 26)] = "ventromedial thalamus"
regions_180131[c(3,13,18,19,20,21,22,23,24, 27,29)] = "posterior hypothalamus"
regions_180131[c(7,8,9,10,11,12,25,28,30,31)] = "dorsomedial hypothalamus"
regions_180131[c(17)] = "paraventricular hypothalamus"
nodeid <- rep(1:nrow(df_180131),2);
condition <- c(rep("Ongoing",nrow(df_180131)), rep("LightON",nrow(df_180131)));
ei_180131 <- array("Excitatory", length(regions_180131));
ei_180131[c(5,6,8,9,12,17,19,23,24,25,26,28)] <- "Inhibitory";
dat_180131 <- data.frame(c(df_180131$Ongoing, df_180131$LightON), condition, nodeid, rep(regions_180131,2), ei_180131);
colnames(dat_180131) <- c("FiringRate", "condition", "nodeid", "region", "EI")
dat_180131$condition <- factor(dat_180131$condition, levels = c("Ongoing", "LightON"))
dat_180131$region <- factor(dat_180131$region)
dat_180131$FiringRate <- as.numeric(sub(",", ".", c(df_180131$Ongoing, df_180131$LightON), fixed = TRUE))
dat_180131$EI <- factor(ei_180131);
model_180131 <- lm(formula = "FiringRate ~ condition * region * EI + (1|nodeid)", data=dat_180131)
res_aov_180131 <- aov(model_180131)
summary(res_aov_180131)
cm=compare_means(FiringRate ~ condition, data = dat_180131, group.by = "region")

tmpdf_fr = aggregate(FiringRate ~ region + condition + EI, data = dat_180131, mean);
#tmpdf_fr$expID = rep("180131",1,length(levels(as.factor(tmpdf_fr$region)))*2);
tmpdf_fr$expID = rep("180131",1,nrow(tmpdf_fr));
df_fr = rbind(df_fr, tmpdf_fr);


###
df_180221 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180221_fr")
regions_180221 <- character(nrow(df_180221));
regions_180221[c(1, 35, 36, 48,49, 50, 2,3,4,12,17,18, 21, 22, 26, 31, 32, 33, 39, 40, 41, 42, 43, 44, 45, 46, 47)] = "ventromedial thalamus"
regions_180221[c(23, 24, 25)] = "ventromedial hypothalamus"
regions_180221[c(5,6,7,8,9,10,11,13,14,15,16, 27, 28, 29, 30, 37, 38)] = "dorsomedial hypothalamus"
regions_180221[c(19, 20, 34, 51, 52, 53, 54)] = "paraventricular hypothalamus"
nodeid <- rep(1:nrow(df_180221),2);
condition <- c(rep("Ongoing",nrow(df_180221)), rep("LightON",nrow(df_180221)));
ei_180221 <- array("Excitatory", length(regions_180221));
ei_180221[c(1,2,5,7,8,9,10,14,15,19,20,23,25,28,30,32,33,34,37,42,43,44,51,52,53,54)] <- "Inhibitory";
dat_180221 <- data.frame(c(df_180221$Ongoing, df_180221$LightON), condition, nodeid, rep(regions_180221,2), ei_180221);
colnames(dat_180221) <- c("FiringRate", "condition", "nodeid", "region", "EI")
dat_180221$condition <- factor(dat_180221$condition, levels = c("Ongoing", "LightON"))
dat_180221$region <- factor(dat_180221$region)
dat_180221$FiringRate <- as.numeric(sub(",", ".", c(df_180221$Ongoing, df_180221$LightON), fixed = TRUE))
dat_180221$EI <- factor(ei_180221);
model_180221 <- lm(formula = "FiringRate ~ condition * region *EI + (1|nodeid)", data=dat_180221)
res_aov_180221 <- aov(model_180221)
summary(res_aov_180221)
cm=compare_means(FiringRate ~ condition, data = dat_180221, group.by = "region")

tmpdf_fr = aggregate(FiringRate ~ region + condition + EI, data = dat_180221, mean);
#tmpdf_fr$expID = rep("180131",1,length(levels(as.factor(tmpdf_fr$region)))*2);
tmpdf_fr$expID = rep("180221",1,nrow(tmpdf_fr));
df_fr = rbind(df_fr, tmpdf_fr);


###
df_180228 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180228_fr")
regions_180228 <- character(nrow(df_180228));
regions_180228[c(1,4,5,6, 24, 25, 74)] = "ventromedial thalamus"
regions_180228[c(2, 11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23,26, 40, 42, 62, 63, 64)] = "posterior hypothalamus"
regions_180228[c(3,7,8,10, 19, 20, 27, 28, 29, 30, 31, 54, 56, 57, 58, 59, 60, 61, 65, 66,67,68,69, 73, 75, 76, 77,78,79)] = "dorsomedial hypothalamus"
regions_180228[c(9,46, 47, 48, 49, 50)] = "ventromedial hypothalamus"
regions_180228[c(32, 33, 34, 35, 36, 37, 38,39, 41)] = "arcuate hypothalamus"
regions_180228[c(43,44,45,51,52,53,55,70,71, 72)] = "zona incerta"
nodeid <- rep(1:nrow(df_180228),2);
condition <- c(rep("Ongoing",nrow(df_180228)), rep("LightON",nrow(df_180228)));
ei_180228 <- array("Excitatory", length(regions_180228));
ei_180228[c(1,4,14,16,17,18,19,20,21,23,24,25,28,29,30,32,33,34,35,36,39,40,46,47,50,51,52,54,56,59,60,61,62)] <- "Inhibitory";
dat_180228 <- data.frame(c(df_180228$Ongoing, df_180228$LightON), condition, nodeid, rep(regions_180228,2), ei_180228);
colnames(dat_180228) <- c("FiringRate", "condition", "nodeid", "region", "EI")
dat_180228$condition <- factor(dat_180228$condition, levels = c("Ongoing", "LightON"))
dat_180228$region <- factor(dat_180228$region)
dat_180228$FiringRate <- as.numeric(sub(",", ".", c(df_180228$Ongoing, df_180228$LightON), fixed = TRUE))
dat_180228$EI <- factor(ei_180228);
model_180228 <- lm(formula = "FiringRate ~ condition * region * EI + (1|nodeid)", data=dat_180228)
res_aov_180228 <- aov(model_180228)
summary(res_aov_180228)
cm=compare_means(FiringRate ~ condition, data = dat_180228, group.by = "region")

tmpdf_fr = aggregate(FiringRate ~ region + condition + EI, data = dat_180228, mean);
#tmpdf_fr$expID = rep("180131",1,length(levels(as.factor(tmpdf_fr$region)))*2);
tmpdf_fr$expID = rep("180228",1,nrow(tmpdf_fr));
df_fr = rbind(df_fr, tmpdf_fr);


####
df_180302 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180302_fr")
regions_180302 <- character(nrow(df_180302));
#regions_180302[c(2,3,4,5,6, 13, 14, 15,30,36)] = "tuberomammillary nucleus"
regions_180302[c(2,3,4,5,6, 13, 14, 15,30,36, 9, 10, 41,42, 18, 19, 20, 21)] = "mammillary complex"
regions_180302[c(1,11,12,7, 8,16, 17, 23, 31, 37, 38,39,40, 47)] = "posterior hypothalamus"
#regions_180302[c(9, 10, 41,42)] = "premammillary nucleus"
#regions_180302[c(18, 19, 20, 21)] = "supramammillary nucleus"
regions_180302[c(26, 27, 28, 43,44,45,46)] = "arcuate hypothalamus"
regions_180302[c(22, 24, 25, 29)] = "paraventricular hypothalamus"
regions_180302[c(32,33,34,35)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_180302),2);
condition <- c(rep("Ongoing",nrow(df_180302)), rep("LightON",nrow(df_180302)));
ei_180302 <- array("Excitatory", length(regions_180302));
ei_180302[c(1,3,8,9,11,15,16,17,19,20,25,26,28,29,30,33,34,35,36,37,40,41,42,47)] <- "Inhibitory";
dat_180302 <- data.frame(c(df_180302$Ongoing, df_180302$LightON), condition, nodeid, rep(regions_180302,2), ei_180302);
colnames(dat_180302) <- c("FiringRate", "condition", "nodeid", "region", "EI")
dat_180302$condition <- factor(dat_180302$condition, levels = c("Ongoing", "LightON"))
dat_180302$region <- factor(dat_180302$region)
dat_180302$FiringRate <- as.numeric(sub(",", ".", c(df_180302$Ongoing, df_180302$LightON), fixed = TRUE))
dat_180302$EI <- factor(ei_180302)
model_180302 <- lm(formula = "FiringRate ~ condition * region *EI + (1|nodeid)", data=dat_180302)
res_aov_180302 <- aov(model_180302)
summary(res_aov_180302)
cm=compare_means(FiringRate ~ condition, data = dat_180302, group.by = "region")

tmpdf_fr = aggregate(FiringRate ~ region + condition + EI, data = dat_180302, mean);
#tmpdf_fr$expID = rep("180131",4,length(levels(as.factor(tmpdf_fr$region)))*2);
tmpdf_fr$expID = rep("180302",1,nrow(tmpdf_fr));
df_fr = rbind(df_fr, tmpdf_fr);

###
df_180419 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180419_fr")
regions_180419 <- character(nrow(df_180419));
regions_180419[c(1,2,11, 12,9, 10)] = "ventromedial hypothalamus"
regions_180419[c(4,5,6,7,8)] = "ventromedial thalamus"
regions_180419[c(36,37,38,39,51,52,53,54,55,56,57,58,59,60,61)] = "posterior hypothalamus"
regions_180419[c(3, 40, 41, 42, 43, 44,45,46,47,48,49,50, 13,14,15,16,17,18,23, 19,20,21,22,24,25,26,27,28,29,30,31,32,33,34,35)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_180419),2);
condition <- c(rep("Ongoing",nrow(df_180419)), rep("LightON",nrow(df_180419)));
ei_180419 <- array("Excitatory", length(regions_180419));
ei_180419[c(1,3,4,5,6,7,9,10,11,15,16,17,19,20,22,33,34,41,42,43,44,47,48,55,58,59)] <- "Inhibitory";
dat_180419 <- data.frame(c(df_180419$Ongoing, df_180419$LightON), condition, nodeid, rep(regions_180419,2), ei_180419);
colnames(dat_180419) <- c("FiringRate", "condition", "nodeid", "region", "EI")
dat_180419$condition <- factor(dat_180419$condition, levels = c("Ongoing", "LightON"))
dat_180419$region <- factor(dat_180419$region)
dat_180419$FiringRate <- as.numeric(sub(",", ".", c(df_180419$Ongoing, df_180419$LightON), fixed = TRUE))
dat_180419$EI <- factor(ei_180419)
model_180419 <- lm(formula = "FiringRate ~ condition * region *EI + (1|nodeid)", data=dat_180419)
res_aov_180419 <- aov(model_180419)
summary(res_aov_180419)
cm=compare_means(FiringRate ~ condition, data = dat_180419, group.by = "region")

tmpdf_fr = aggregate(FiringRate ~ region + condition + EI, data = dat_180419, mean);
#tmpdf_fr$expID = rep("180131",4,length(levels(as.factor(tmpdf_fr$region)))*2);
tmpdf_fr$expID = rep("180419",1,nrow(tmpdf_fr));
df_fr = rbind(df_fr, tmpdf_fr);

####
df_180420 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180420_fr")
regions_180420 <- character(nrow(df_180420));
#regions_180420[c(51,52,53,54,55,56,57,58,59,60,61,62,63)] = "supramammillary nucleus"
regions_180420[c(51,52,53,54,55,56,57,58,59,60,61,62,63,19, 20, 21, 22, 23, 24, 25,26,27,28,29,30,31)] = "mammillary complex"
regions_180420[c(32, 33, 34, 35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50)] = "paraventricular hypothalamus"
#regions_180420[c(19, 20, 21, 22, 23, 24, 25,26,27,28,29,30,31)] = "premammillary nucleus"
regions_180420[c(2, 3, 4, 5, 6,7,8,9,10,11,12,13,14,15,16,17,18,64,65,66,67,68,69,70,71)] = "posterior hypothalamus"
regions_180420[c(1)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_180420),2);
condition <- c(rep("Ongoing",nrow(df_180420)), rep("LightON",nrow(df_180420)));
ei_180420 <- array("Excitatory", length(regions_180420));
ei_180420[c(4,5,6,7,9,11,12,13,15,16,17,18,19,20,28,29,30,31,41,42,49,51,52,53,59,62,64,68,69,70,71)] <- "Inhibitory";
dat_180420 <- data.frame(c(df_180420$Ongoing, df_180420$LightON), condition, nodeid, rep(regions_180420,2), ei_180420);
colnames(dat_180420) <- c("FiringRate", "condition", "nodeid", "region", "EI")
dat_180420$condition <- factor(dat_180420$condition, levels = c("Ongoing", "LightON"))
dat_180420$region <- factor(dat_180420$region)
dat_180420$FiringRate <- as.numeric(sub(",", ".", c(df_180420$Ongoing, df_180420$LightON), fixed = TRUE))
dat_180420$EI <- factor(ei_180420)
model_180420 <- lm(formula = "FiringRate ~ condition * region * EI + (1|nodeid)", data=dat_180420)
res_aov_180420 <- aov(model_180420)
summary(res_aov_180420)
cm=compare_means(FiringRate ~ condition, data = dat_180420, group.by = "region")

tmpdf_fr = aggregate(FiringRate ~ region + condition + EI, data = dat_180420, mean);
#tmpdf_fr$expID = rep("180131",4,length(levels(as.factor(tmpdf_fr$region)))*2);
tmpdf_fr$expID = rep("180420",1,nrow(tmpdf_fr));
df_fr = rbind(df_fr, tmpdf_fr);

###
df_180423 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180423_fr")
regions_180423 <- character(nrow(df_180423));
regions_180423[c(1,2,3,4)] = "ventromedial hypothalamus"
#regions_180423[c(5)] = "submammillothalamic nucleus"
regions_180423[c(6,7,8,9,10,18,19,20,26,27,28,29,30,31,32,33,34,35,36)] = "dorsomedial hypothalamus"
regions_180423[c(11, 12,13,14,15,16,17,21,22,23,24,25,51,52,53,54,55,56)] = "posterior hypothalamus"
#regions_180423[c(37, 38,39,40,41,42,43,44,45,46,47,48,49,50)] = "premammillary nucleus"
regions_180423[c(37, 38,39,40,41,42,43,44,45,46,47,48,49,50,5)] = "mammillary complex"
nodeid <- rep(1:nrow(df_180423),2);
condition <- c(rep("Ongoing",nrow(df_180423)), rep("LightON",nrow(df_180423)));
ei_180423 <- array("Excitatory", length(regions_180423));
ei_180423[c(1,2,3,4,8,12,14,16,19,26,27,29,31,32,35,38,40,41,45,47,48,49,50,51,52,54)] <- "Inhibitory";
dat_180423 <- data.frame(c(df_180423$Ongoing, df_180423$LightON), condition, nodeid, rep(regions_180423,2), ei_180423);
colnames(dat_180423) <- c("FiringRate", "condition", "nodeid", "region", "EI")
dat_180423$condition <- factor(dat_180423$condition, levels = c("Ongoing", "LightON"))
dat_180423$region <- factor(dat_180423$region)
dat_180423$FiringRate <- as.numeric(sub(",", ".", c(df_180423$Ongoing, df_180423$LightON), fixed = TRUE))
dat_180423$EI <- factor(ei_180423);
model_180423 <- lm(formula = "FiringRate ~ condition * region * EI + (1|nodeid)", data=dat_180423)
res_aov_180423 <- aov(model_180423)
summary(res_aov_180423)
cm=compare_means(FiringRate ~ condition, data = dat_180423, group.by = "region")

tmpdf_fr = aggregate(FiringRate ~ region + condition + EI, data = dat_180423, mean);
#tmpdf_fr$expID = rep("180131",4,length(levels(as.factor(tmpdf_fr$region)))*2);
tmpdf_fr$expID = rep("180423",1,nrow(tmpdf_fr));
df_fr = rbind(df_fr, tmpdf_fr);

df_fr$expID <- factor(df_fr$expID)
model_fr <- lm(formula = "FiringRate ~ condition * region * EI", data = df_fr)
res_aov_fr <- aov(model_fr)
summary(res_aov_fr)
posthoc <- TukeyHSD(res_aov_fr)
cm=compare_means(FiringRate ~ condition, data = df_fr, group.by = "region")
ttt <- posthoc$`region:EI`
which(ttt[,4] < 0.05)


df_fr$phase[df_fr$expID == "171019"] = rep("Day", length(df_fr$expID == "171019"))
df_fr$phase[df_fr$expID == "180131"] = rep("Day", length(df_fr$expID == "180131"))
df_fr$phase[df_fr$expID == "180302"] = rep("Day", length(df_fr$expID == "180302"))
df_fr$phase[df_fr$expID == "180131"] = rep("Day", length(df_fr$expID == "180131"))
df_fr$phase[df_fr$expID == "180419"] = rep("Day", length(df_fr$expID == "180419"))
df_fr$phase[df_fr$expID == "180420"] = rep("Day", length(df_fr$expID == "180420"))
df_fr$phase[df_fr$expID == "180423"] = rep("Day", length(df_fr$expID == "180423"))
df_fr$phase[is.na(df_fr$phase)] = "Night";


model_fr <- lm(formula = "FiringRate ~ condition * region * EI * phase", data = df_fr)
res_aov_fr <- aov(model_fr)
summary(res_aov_fr)
cm=compare_means(FiringRate ~ condition + region, data = df_fr, group.by = "region")
posthoc <- TukeyHSD(res_aov_fr)
cm=compare_means(FiringRate ~ condition, data = df_fr, group.by = "region")
ttt <- posthoc$`region:phase`
which(ttt[,4] < 0.05)

model_fr <- lm(formula = "FiringRate ~ region + condition + phase", data = df_fr)
res_aov_fr <- aov(model_fr)
summary(res_aov_fr)

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/region_fr_phase.png", width = 1860, height = 1440, units = "px")
ggboxplot(df_fr, x="region", y="FiringRate", color="EI", add = "jitter") + #stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("")  + ylab("Firing Rate (Hz)")
dev.off()







#### NODE STRENGTH
##################################################################################
#################################################################################
df_171207 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171207_deg_wei")
regions_171207 <- character(nrow(df_171207));
regions_171207[c(1,2,3, 12,13,28, 18,23,24,25,26,27,31)] = "ventromedial thalamus"
regions_171207[c(4,5,6,7,11,14,15,16,29,30)] = "paraventricular hypothalamus"
regions_171207[c(8,9,10,17)] = "anterior hypothalamus"
regions_171207[c(19,20,21,22)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_171207),2);
condition <- c(rep("Ongoing",nrow(df_171207)), rep("LightON",nrow(df_171207)));
ei_171207 <- array("Excitatory", length(regions_171207));
ei_171207[9] <- "Inhibitory";
dat_171207 <- data.frame(c(df_171207$Ongoing, df_171207$LightON), condition, nodeid, rep(regions_171207,2), ei_171207);
colnames(dat_171207) <- c("Strength", "condition", "nodeid", "region", "EI");
dat_171207$condition <- factor(dat_171207$condition, levels = c("Ongoing", "LightON"))
dat_171207$region <- factor(dat_171207$region)
dat_171207$Strength <- as.numeric(sub(",", ".", c(df_171207$Ongoing, df_171207$LightON), fixed = TRUE))
dat_171207$EI <- factor(ei_171207);
model_171207 <- lm(formula = "Strength ~ condition * region * EI + (1|nodeid)", data=dat_171207)
res_aov_171207 <- aov(model_171207)
summary(res_aov_171207)
cm=compare_means(Strength ~ condition, data = dat_171207, group.by = "region")

###searching for hubs###
potout = boxplot.stats(dat_171207$Strength[dat_171207$condition == "Ongoing"], 2.3)$out
potout_id = dat_171207$nodeid[which(dat_171207$Strength[dat_171207$condition == "Ongoing"] %in%  potout)]
potout_evoked = boxplot.stats(dat_171207$Strength[dat_171207$condition == "LightON"], 2.3)$out
potout_id_evoked = dat_171207$nodeid[which(dat_171207$Strength[dat_171207$condition == "LightON"] %in%  potout_evoked)]
x = list(ongoing = potout_id, lightON = potout_id_evoked);
png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/hubness_171207.png")
ggvenn(x)
dev.off()



tmpdf_deg_wei = aggregate(Strength ~ region + condition + EI, data = dat_171207, mean);
#tmpdf_deg_wei$expID = rep("171207",1,length(levels(as.factor(dat_171207)))*2);
tmpdf_deg_wei$expID = rep("171207",1,length(dat_171207)*2);
df_deg_wei = tmpdf_deg_wei;
#df_degree = rbind(df_degree, tmpdf_degree);

###
df_171019 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171019_deg_wei")
regions_171019 <- character(nrow(df_171019));
regions_171019[seq(1,114)] = "ventromedial thalamus"
regions_171019[c(115, 116)] = "paraventricular hypothalamus"
regions_171019[seq(117,126)] = "ventromedial thalamus"
nodeid <- rep(1:nrow(df_171019),2);
condition <- c(rep("Ongoing",nrow(df_171019)), rep("LightON",nrow(df_171019)));
ei_171019 <- array("Excitatory", length(regions_171019));
ei_171019[c(60,125)] <- "Inhibitory";
dat_171019 <- data.frame(c(df_171019$Ongoing, df_171019$LightON), condition, nodeid, rep(regions_171019,2), ei_171019);
colnames(dat_171019) <- c("Strength", "condition", "nodeid", "region", "EI")
dat_171019$condition <- factor(dat_171019$condition, levels = c("Ongoing", "LightON"))
dat_171019$region <- factor(dat_171019$region)
dat_171019$Strength <- as.numeric(sub(",", ".", c(df_171019$Ongoing, df_171019$LightON), fixed = TRUE))
dat_171019$EI <- factor(ei_171019);
model_171019 <- lm(formula = "Strength ~ condition * region * EI + (1|nodeid)", data=dat_171019)
res_aov_171019 <- aov(model_171019)
summary(res_aov_171019)
cm=compare_means(Strength ~ condition, data = dat_171019, group.by = "region")

potout = boxplot.stats(dat_171019$Strength[dat_171019$condition == "Ongoing"], 2.3)$out
potout_id = dat_171019$nodeid[which(dat_171019$Strength[dat_171019$condition == "Ongoing"] %in%  potout)]
potout_evoked = boxplot.stats(dat_171019$Strength[dat_171019$condition == "LightON"], 2.3)$out
potout_id_evoked = dat_171019$nodeid[which(dat_171019$Strength[dat_171019$condition == "LightON"] %in%  potout_evoked)]
x = list(ongoing = potout_id, lightON = potout_id_evoked);
png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/hubness_171019.png")
ggvenn(x)
dev.off()


tmpdf_deg_wei = aggregate(Strength ~ region + condition + EI, data = dat_171019, mean);
#tmpdf_deg_wei$expID = rep("171019",1,length(levels(as.factor(tmpdf_deg_wei$region)))*2);
tmpdf_deg_wei$expID = rep("171019",1,nrow(tmpdf_deg_wei));
df_deg_wei = rbind(df_deg_wei, tmpdf_deg_wei);


###
df_171208 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171208_deg_wei")
regions_171208 <- character(nrow(df_171208));
regions_171208[c(1,2,3,4,5,7,20,21,22,23,24,25,28)] = "ventromedial thalamus"
regions_171208[c(6, 8,9,11,13,16,17)] = "ventromedial hypothalamus"
regions_171208[c(10,12, 32)] = "paraventricular hypothalamus"
regions_171208[c(14,15,18,19)] = "dorsomedial hypothalamus"
regions_171208[c(26, 27, 29, 30, 31,33, 34)] = "anterior hypothalamus"
nodeid <- rep(1:nrow(df_171208),2);
condition <- c(rep("Ongoing",nrow(df_171208)), rep("LightON",nrow(df_171208)));
ei_171208 <- array("Excitatory", length(regions_171208));
dat_171208 <- data.frame(c(df_171208$Ongoing, df_171208$LightON), condition, nodeid, rep(regions_171208,2), ei_171208);
colnames(dat_171208) <- c("Strength", "condition", "nodeid", "region", "EI")
dat_171208$condition <- factor(dat_171208$condition, levels = c("Ongoing", "LightON"))
dat_171208$region <- factor(dat_171208$region)
dat_171208$Strength <- as.numeric(sub(",", ".", c(df_171208$Ongoing, df_171208$LightON), fixed = TRUE))
dat_171208$EI <- factor(ei_171208);
model_171208 <- lm(formula = "Strength ~ condition * region + (1|nodeid)", data=dat_171208)
res_aov_171208 <- aov(model_171208)
summary(res_aov_171208)
cm=compare_means(Strength ~ condition, data = dat_171208, group.by = "region")

potout = boxplot.stats(dat_171208$Strength[dat_171208$condition == "Ongoing"], 2.3)$out
potout_id = dat_171208$nodeid[which(dat_171208$Strength[dat_171208$condition == "Ongoing"] %in%  potout)]
potout_evoked = boxplot.stats(dat_171208$Strength[dat_171208$condition == "LightON"], 2.3)$out
potout_id_evoked = dat_171208$nodeid[which(dat_171208$Strength[dat_171208$condition == "LightON"] %in%  potout_evoked)]
x = list(ongoing = potout_id, lightON = potout_id_evoked);
png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/hubness_171208.png")
ggvenn(x)
dev.off()

tmpdf_deg_wei = aggregate(Strength ~ region + condition + EI, data = dat_171208, mean);
#tmpdf_deg_wei$expID = rep("171208",1,length(levels(as.factor(tmpdf_deg_wei$region)))*2);
tmpdf_deg_wei$expID = rep("171208",1,nrow(tmpdf_deg_wei));
df_deg_wei = rbind(df_deg_wei, tmpdf_deg_wei);



###
df_171213 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171213_deg_wei")
regions_171213 <- character(nrow(df_171213));
regions_171213[c(11,12,13,20,21,22,23,24,25,26,27,29,30,31,32)] = "ventromedial thalamus"
regions_171213[c(1,2,3,4,5,6,7,8,9,10,14,15,16,17,18,19,28,33)] = "paraventricular hypothalamus"
nodeid <- rep(1:nrow(df_171213),2);
condition <- c(rep("Ongoing",nrow(df_171213)), rep("LightON",nrow(df_171213)));
ei_171213 <- array("Excitatory", length(regions_171213));
dat_171213 <- data.frame(c(df_171213$Ongoing, df_171213$LightON), condition, nodeid, rep(regions_171213,2), ei_171213);
colnames(dat_171213) <- c("Strength", "condition", "nodeid", "region", "EI")
dat_171213$condition <- factor(dat_171213$condition, levels = c("Ongoing", "LightON"))
dat_171213$region <- factor(dat_171213$region)
dat_171213$Strength <- as.numeric(sub(",", ".", c(df_171213$Ongoing, df_171213$LightON), fixed = TRUE))
dat_171213$EI <- factor(ei_171213);
model_171213 <- lm(formula = "Strength ~ condition * region + (1|nodeid)", data=dat_171213)
res_aov_171213 <- aov(model_171213)
summary(res_aov_171213)
cm=compare_means(Strength ~ condition, data = dat_171213, group.by = "region")

potout = boxplot.stats(dat_171213$Strength[dat_171213$condition == "Ongoing"], 2.3)$out
potout_id = dat_171213$nodeid[which(dat_171213$Strength[dat_171213$condition == "Ongoing"] %in%  potout)]
potout_evoked = boxplot.stats(dat_171213$Strength[dat_171213$condition == "LightON"], 2.3)$out
potout_id_evoked = dat_171213$nodeid[which(dat_171213$Strength[dat_171213$condition == "LightON"] %in%  potout_evoked)]
x = list(ongoing = potout_id, lightON = potout_id_evoked);
png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/hubness_171213.png")
ggvenn(x)
dev.off()


tmpdf_deg_wei = aggregate(Strength ~ region + condition + EI, data = dat_171213, mean);
#tmpdf_deg_wei$expID = rep("171213",1,length(levels(as.factor(tmpdf_deg_wei$region)))*2);
tmpdf_deg_wei$expID = rep("171213",1,nrow(tmpdf_deg_wei));
df_deg_wei = rbind(df_deg_wei, tmpdf_deg_wei);


###
df_180110 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180110_deg_wei")
regions_180110 <- character(nrow(df_180110));
regions_180110[c(1,2,11,18,20,21,23)] = "ventromedial thalamus"
regions_180110[c(3,4,5,6,7,9,10, 19,22,24,25,26,27,28,29,30,31,32,33,34,35)] = "paraventricular hypothalamus"
regions_180110[c(8)] = "anterior hypothalamus"
regions_180110[c(12,13,14,15,16,17)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_180110),2);
condition <- c(rep("Ongoing",nrow(df_180110)), rep("LightON",nrow(df_180110)));
ei_180110 <- array("Excitatory", length(regions_180110));
ei_180110[16] <- "Inhibitory"
dat_180110 <- data.frame(c(df_180110$Ongoing, df_180110$LightON), condition, nodeid, rep(regions_180110,2), ei_180110);
colnames(dat_180110) <- c("Strength", "condition", "nodeid", "region", "EI")
dat_180110$condition <- factor(dat_180110$condition, levels = c("Ongoing", "LightON"))
dat_180110$region <- factor(dat_180110$region)
dat_180110$Strength <- as.numeric(sub(",", ".", c(df_180110$Ongoing, df_180110$LightON), fixed = TRUE))
dat_180110$EI <- factor(ei_180110)
model_180110 <- lm(formula = "Strength ~ condition * region * EI + (1|nodeid)", data=dat_180110)
res_aov_180110 <- aov(model_180110)
summary(res_aov_180110)
cm=compare_means(Strength ~ condition, data = dat_180110, group.by = "region")

potout = boxplot.stats(dat_180110$Strength[dat_180110$condition == "Ongoing"])$out
potout_id = dat_180110$nodeid[which(dat_180110$Strength[dat_180110$condition == "Ongoing"] %in%  potout)]
potout_evoked = boxplot.stats(dat_180110$Strength[dat_180110$condition == "LightON"])$out
potout_id_evoked = dat_180110$nodeid[which(dat_180110$Strength[dat_180110$condition == "LightON"] %in%  potout_evoked)]
x = list(ongoing = potout_id, lightON = potout_id_evoked);
png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/hubness_180110.png")
ggvenn(x)
dev.off()


tmpdf_deg_wei = aggregate(Strength ~ region + condition + EI, data = dat_180110, mean);
#tmpdf_deg_wei$expID = rep("180110",1,length(levels(as.factor(tmpdf_deg_wei$region)))*2);
tmpdf_deg_wei$expID = rep("180110",1,nrow(tmpdf_deg_wei));
df_deg_wei = rbind(df_deg_wei, tmpdf_deg_wei);

###
df_180111 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180111_deg_wei")
regions_180111 <- character(nrow(df_180111));
regions_180111[c(1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,22,23,24)] = "ventromedial thalamus"
regions_180111[c(14)] = "paraventricular hypothalamus"
regions_180111[c(19,20,21)] = "ventromedial hypothalamus"
nodeid <- rep(1:nrow(df_180111),2);
condition <- c(rep("Ongoing",nrow(df_180111)), rep("LightON",nrow(df_180111)));
ei_180111 <- array("Excitatory", length(regions_180111));
dat_180111 <- data.frame(c(df_180111$Ongoing, df_180111$LightON), condition, nodeid, rep(regions_180111,2), ei_180111);
colnames(dat_180111) <- c("Strength", "condition", "nodeid", "region", "EI")
dat_180111$condition <- factor(dat_180111$condition, levels = c("Ongoing", "LightON"))
dat_180111$region <- factor(dat_180111$region)
dat_180111$Strength <- as.numeric(sub(",", ".", c(df_180111$Ongoing, df_180111$LightON), fixed = TRUE))
dat_180111$EI <- factor(ei_180111);
model_180111 <- lm(formula = "Strength ~ condition * region + (1|nodeid)", data=dat_180111)
res_aov_180111 <- aov(model_180111)
summary(res_aov_180111)
cm=compare_means(Strength ~ condition, data = dat_180111, group.by = "region")

potout = boxplot.stats(dat_180111$Strength[dat_180111$condition == "Ongoing"])$out
potout_id = dat_180111$nodeid[which(dat_180111$Strength[dat_180111$condition == "Ongoing"] %in%  potout)]
potout_evoked = boxplot.stats(dat_180111$Strength[dat_180111$condition == "LightON"])$out
potout_id_evoked = dat_180111$nodeid[which(dat_180111$Strength[dat_180111$condition == "LightON"] %in%  potout_evoked)]
x = list(ongoing = potout_id, lightON = potout_id_evoked);
png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/hubness_180111.png")
ggvenn(x)
dev.off()

tmpdf_deg_wei = aggregate(Strength ~ region + condition + EI, data = dat_180111, mean);
#tmpdf_deg_wei$expID = rep("180111",1,length(levels(as.factor(tmpdf_deg_wei$region)))*2);
tmpdf_deg_wei$expID = rep("180111",1,nrow(tmpdf_deg_wei));
df_deg_wei = rbind(df_deg_wei, tmpdf_deg_wei);


###
df_180131 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180131_deg_wei")
regions_180131 <- character(nrow(df_180131));
regions_180131[c(1,2,4,5,6,14,15,16, 26)] = "ventromedial thalamus"
regions_180131[c(3,13,18,19,20,21,22,23,24, 27,29)] = "posterior hypothalamus"
regions_180131[c(7,8,9,10,11,12,25,28,30,31)] = "dorsomedial hypothalamus"
regions_180131[c(17)] = "paraventricular hypothalamus"
nodeid <- rep(1:nrow(df_180131),2);
condition <- c(rep("Ongoing",nrow(df_180131)), rep("LightON",nrow(df_180131)));
ei_180131 <- array("Excitatory", length(regions_180131));
ei_180131[c(5,6,8,9,12,17,19,23,24,25,26,28)] <- "Inhibitory";
dat_180131 <- data.frame(c(df_180131$Ongoing, df_180131$LightON), condition, nodeid, rep(regions_180131,2), ei_180131);
colnames(dat_180131) <- c("Strength", "condition", "nodeid", "region", "EI")
dat_180131$condition <- factor(dat_180131$condition, levels = c("Ongoing", "LightON"))
dat_180131$region <- factor(dat_180131$region)
dat_180131$Strength <- as.numeric(sub(",", ".", c(df_180131$Ongoing, df_180131$LightON), fixed = TRUE))
dat_180131$EI <- factor(ei_180131);
model_180131 <- lm(formula = "Strength ~ condition * region * EI + (1|nodeid)", data=dat_180131)
res_aov_180131 <- aov(model_180131)
summary(res_aov_180131)
cm=compare_means(Strength ~ condition, data = dat_180131, group.by = "region")

potout = boxplot.stats(dat_180131$Strength[dat_180131$condition == "Ongoing"])$out
potout_id = dat_180131$nodeid[which(dat_180131$Strength[dat_180131$condition == "Ongoing"] %in%  potout)]
potout_evoked = boxplot.stats(dat_180131$Strength[dat_180131$condition == "LightON"])$out
potout_id_evoked = dat_180131$nodeid[which(dat_180131$Strength[dat_180131$condition == "LightON"] %in%  potout_evoked)]
x = list(ongoing = potout_id, lightON = potout_id_evoked);
png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/hubness_180131.png")
ggvenn(x)
dev.off()


tmpdf_deg_wei = aggregate(Strength ~ region + condition + EI, data = dat_180131, mean);
#tmpdf_deg_wei$expID = rep("180131",1,length(levels(as.factor(tmpdf_deg_wei$region)))*2);
tmpdf_deg_wei$expID = rep("180131",1,nrow(tmpdf_deg_wei));
df_deg_wei = rbind(df_deg_wei, tmpdf_deg_wei);


###
df_180221 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180221_deg_wei")
regions_180221 <- character(nrow(df_180221));
regions_180221[c(1, 35, 36, 48,49, 50, 2,3,4,12,17,18, 21, 22, 26, 31, 32, 33, 39, 40, 41, 42, 43, 44, 45, 46, 47)] = "ventromedial thalamus"
regions_180221[c(23, 24, 25)] = "ventromedial hypothalamus"
regions_180221[c(5,6,7,8,9,10,11,13,14,15,16, 27, 28, 29, 30, 37, 38)] = "dorsomedial hypothalamus"
regions_180221[c(19, 20, 34, 51, 52, 53, 54)] = "paraventricular hypothalamus"
nodeid <- rep(1:nrow(df_180221),2);
condition <- c(rep("Ongoing",nrow(df_180221)), rep("LightON",nrow(df_180221)));
ei_180221 <- array("Excitatory", length(regions_180221));
ei_180221[c(1,2,5,7,8,9,10,14,15,19,20,23,25,28,30,32,33,34,37,42,43,44,51,52,53,54)] <- "Inhibitory";
dat_180221 <- data.frame(c(df_180221$Ongoing, df_180221$LightON), condition, nodeid, rep(regions_180221,2), ei_180221);
colnames(dat_180221) <- c("Strength", "condition", "nodeid", "region", "EI")
dat_180221$condition <- factor(dat_180221$condition, levels = c("Ongoing", "LightON"))
dat_180221$region <- factor(dat_180221$region)
dat_180221$Strength <- as.numeric(sub(",", ".", c(df_180221$Ongoing, df_180221$LightON), fixed = TRUE))
dat_180221$EI <- factor(ei_180221);
model_180221 <- lm(formula = "Strength ~ condition * region *EI + (1|nodeid)", data=dat_180221)
res_aov_180221 <- aov(model_180221)
summary(res_aov_180221)
cm=compare_means(Strength ~ condition, data = dat_180221, group.by = "region")

potout = boxplot.stats(dat_180221$Strength[dat_180221$condition == "Ongoing"])$out
potout_id = dat_180221$nodeid[which(dat_180221$Strength[dat_180221$condition == "Ongoing"] %in%  potout)]
potout_evoked = boxplot.stats(dat_180221$Strength[dat_180221$condition == "LightON"])$out
potout_id_evoked = dat_180221$nodeid[which(dat_180221$Strength[dat_180221$condition == "LightON"] %in%  potout_evoked)]
x = list(ongoing = potout_id, lightON = potout_id_evoked);
png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/hubness_180221.png")
ggvenn(x)
dev.off()


tmpdf_deg_wei = aggregate(Strength ~ region + condition + EI, data = dat_180221, mean);
#tmpdf_deg_wei$expID = rep("180131",1,length(levels(as.factor(tmpdf_deg_wei$region)))*2);
tmpdf_deg_wei$expID = rep("180221",1,nrow(tmpdf_deg_wei));
df_deg_wei = rbind(df_deg_wei, tmpdf_deg_wei);


###
df_180228 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180228_deg_wei")
regions_180228 <- character(nrow(df_180228));
regions_180228[c(1,4,5,6, 24, 25, 74)] = "ventromedial thalamus"
regions_180228[c(2, 11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23,26, 40, 42, 62, 63, 64)] = "posterior hypothalamus"
regions_180228[c(3,7,8,10, 19, 20, 27, 28, 29, 30, 31, 54, 56, 57, 58, 59, 60, 61, 65, 66,67,68,69, 73, 75, 76, 77,78,79)] = "dorsomedial hypothalamus"
regions_180228[c(9,46, 47, 48, 49, 50)] = "ventromedial hypothalamus"
regions_180228[c(32, 33, 34, 35, 36, 37, 38,39, 41)] = "arcuate hypothalamus"
regions_180228[c(43,44,45,51,52,53,55,70,71, 72)] = "zona incerta"
nodeid <- rep(1:nrow(df_180228),2);
condition <- c(rep("Ongoing",nrow(df_180228)), rep("LightON",nrow(df_180228)));
ei_180228 <- array("Excitatory", length(regions_180228));
ei_180228[c(1,4,14,16,17,18,19,20,21,23,24,25,28,29,30,32,33,34,35,36,39,40,46,47,50,51,52,54,56,59,60,61,62)] <- "Inhibitory";
dat_180228 <- data.frame(c(df_180228$Ongoing, df_180228$LightON), condition, nodeid, rep(regions_180228,2), ei_180228);
colnames(dat_180228) <- c("Strength", "condition", "nodeid", "region", "EI")
dat_180228$condition <- factor(dat_180228$condition, levels = c("Ongoing", "LightON"))
dat_180228$region <- factor(dat_180228$region)
dat_180228$Strength <- as.numeric(sub(",", ".", c(df_180228$Ongoing, df_180228$LightON), fixed = TRUE))
dat_180228$EI <- factor(ei_180228);
model_180228 <- lm(formula = "Strength ~ condition * region * EI + (1|nodeid)", data=dat_180228)
res_aov_180228 <- aov(model_180228)
summary(res_aov_180228)
cm=compare_means(Strength ~ condition, data = dat_180228, group.by = "region")

potout = boxplot.stats(dat_180228$Strength[dat_180228$condition == "Ongoing"])$out
potout_id = dat_180228$nodeid[which(dat_180228$Strength[dat_180228$condition == "Ongoing"] %in%  potout)]
potout_evoked = boxplot.stats(dat_180228$Strength[dat_180228$condition == "LightON"])$out
potout_id_evoked = dat_180228$nodeid[which(dat_180228$Strength[dat_180228$condition == "LightON"] %in%  potout_evoked)]
x = list(ongoing = potout_id, lightON = potout_id_evoked);
png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/hubness_180228.png")
ggvenn(x)
dev.off()


tmpdf_deg_wei = aggregate(Strength ~ region + condition + EI, data = dat_180228, mean);
#tmpdf_deg_wei$expID = rep("180131",1,length(levels(as.factor(tmpdf_deg_wei$region)))*2);
tmpdf_deg_wei$expID = rep("180228",1,nrow(tmpdf_deg_wei));
df_deg_wei = rbind(df_deg_wei, tmpdf_deg_wei);


####
df_180302 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180302_deg_wei")
regions_180302 <- character(nrow(df_180302));
#regions_180302[c(2,3,4,5,6, 13, 14, 15,30,36)] = "tuberomammillary nucleus"
regions_180302[c(2,3,4,5,6, 13, 14, 15,30,36, 9, 10, 41,42, 18, 19, 20, 21)] = "mammillary complex"
regions_180302[c(1,11,12,7, 8,16, 17, 23, 31, 37, 38,39,40, 47)] = "posterior hypothalamus"
#regions_180302[c(9, 10, 41,42)] = "premammillary nucleus"
#regions_180302[c(18, 19, 20, 21)] = "supramammillary nucleus"
regions_180302[c(26, 27, 28, 43,44,45,46)] = "arcuate hypothalamus"
regions_180302[c(22, 24, 25, 29)] = "paraventricular hypothalamus"
regions_180302[c(32,33,34,35)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_180302),2);
condition <- c(rep("Ongoing",nrow(df_180302)), rep("LightON",nrow(df_180302)));
ei_180302 <- array("Excitatory", length(regions_180302));
ei_180302[c(1,3,8,9,11,15,16,17,19,20,25,26,28,29,30,33,34,35,36,37,40,41,42,47)] <- "Inhibitory";
dat_180302 <- data.frame(c(df_180302$Ongoing, df_180302$LightON), condition, nodeid, rep(regions_180302,2), ei_180302);
colnames(dat_180302) <- c("Strength", "condition", "nodeid", "region", "EI")
dat_180302$condition <- factor(dat_180302$condition, levels = c("Ongoing", "LightON"))
dat_180302$region <- factor(dat_180302$region)
dat_180302$Strength <- as.numeric(sub(",", ".", c(df_180302$Ongoing, df_180302$LightON), fixed = TRUE))
dat_180302$EI <- factor(ei_180302)
model_180302 <- lm(formula = "Strength ~ condition * region *EI + (1|nodeid)", data=dat_180302)
res_aov_180302 <- aov(model_180302)
summary(res_aov_180302)
cm=compare_means(Strength ~ condition, data = dat_180302, group.by = "region")

potout = boxplot.stats(dat_180302$Strength[dat_180302$condition == "Ongoing"])$out
potout_id = dat_180302$nodeid[which(dat_180302$Strength[dat_180302$condition == "Ongoing"] %in%  potout)]
potout_evoked = boxplot.stats(dat_180302$Strength[dat_180302$condition == "LightON"])$out
potout_id_evoked = dat_180302$nodeid[which(dat_180302$Strength[dat_180302$condition == "LightON"] %in%  potout_evoked)]
x = list(ongoing = potout_id, lightON = potout_id_evoked);
png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/hubness_180302.png")
ggvenn(x)
dev.off()


tmpdf_deg_wei = aggregate(Strength ~ region + condition + EI, data = dat_180302, mean);
#tmpdf_deg_wei$expID = rep("180131",4,length(levels(as.factor(tmpdf_deg_wei$region)))*2);
tmpdf_deg_wei$expID = rep("180302",1,nrow(tmpdf_deg_wei));
df_deg_wei = rbind(df_deg_wei, tmpdf_deg_wei);

###
df_180419 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180419_deg_wei")
regions_180419 <- character(nrow(df_180419));
regions_180419[c(1,2,11, 12,9, 10)] = "ventromedial hypothalamus"
regions_180419[c(4,5,6,7,8)] = "ventromedial thalamus"
regions_180419[c(36,37,38,39,51,52,53,54,55,56,57,58,59,60,61)] = "posterior hypothalamus"
regions_180419[c(3, 40, 41, 42, 43, 44,45,46,47,48,49,50, 13,14,15,16,17,18,23, 19,20,21,22,24,25,26,27,28,29,30,31,32,33,34,35)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_180419),2);
condition <- c(rep("Ongoing",nrow(df_180419)), rep("LightON",nrow(df_180419)));
ei_180419 <- array("Excitatory", length(regions_180419));
ei_180419[c(1,3,4,5,6,7,9,10,11,15,16,17,19,20,22,33,34,41,42,43,44,47,48,55,58,59)] <- "Inhibitory";
dat_180419 <- data.frame(c(df_180419$Ongoing, df_180419$LightON), condition, nodeid, rep(regions_180419,2), ei_180419);
colnames(dat_180419) <- c("Strength", "condition", "nodeid", "region", "EI")
dat_180419$condition <- factor(dat_180419$condition, levels = c("Ongoing", "LightON"))
dat_180419$region <- factor(dat_180419$region)
dat_180419$Strength <- as.numeric(sub(",", ".", c(df_180419$Ongoing, df_180419$LightON), fixed = TRUE))
dat_180419$EI <- factor(ei_180419)
model_180419 <- lm(formula = "Strength ~ condition * region *EI + (1|nodeid)", data=dat_180419)
res_aov_180419 <- aov(model_180419)
summary(res_aov_180419)
cm=compare_means(Strength ~ condition, data = dat_180419, group.by = "region")

potout = boxplot.stats(dat_180419$Strength[dat_180419$condition == "Ongoing"])$out
potout_id = dat_180419$nodeid[which(dat_180419$Strength[dat_180419$condition == "Ongoing"] %in%  potout)]
potout_evoked = boxplot.stats(dat_180419$Strength[dat_180419$condition == "LightON"])$out
potout_id_evoked = dat_180419$nodeid[which(dat_180419$Strength[dat_180419$condition == "LightON"] %in%  potout_evoked)]
x = list(ongoing = potout_id, lightON = potout_id_evoked);
png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/hubness_180419.png")
ggvenn(x)
dev.off()


tmpdf_deg_wei = aggregate(Strength ~ region + condition + EI, data = dat_180419, mean);
#tmpdf_deg_wei$expID = rep("180131",4,length(levels(as.factor(tmpdf_deg_wei$region)))*2);
tmpdf_deg_wei$expID = rep("180419",1,nrow(tmpdf_deg_wei));
df_deg_wei = rbind(df_deg_wei, tmpdf_deg_wei);

####
df_180420 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180420_deg_wei")
regions_180420 <- character(nrow(df_180420));
#regions_180420[c(51,52,53,54,55,56,57,58,59,60,61,62,63)] = "supramammillary nucleus"
regions_180420[c(51,52,53,54,55,56,57,58,59,60,61,62,63,19, 20, 21, 22, 23, 24, 25,26,27,28,29,30,31)] = "mammillary complex"
regions_180420[c(32, 33, 34, 35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50)] = "paraventricular hypothalamus"
#regions_180420[c(19, 20, 21, 22, 23, 24, 25,26,27,28,29,30,31)] = "premammillary nucleus"
regions_180420[c(2, 3, 4, 5, 6,7,8,9,10,11,12,13,14,15,16,17,18,64,65,66,67,68,69,70,71)] = "posterior hypothalamus"
regions_180420[c(1)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_180420),2);
condition <- c(rep("Ongoing",nrow(df_180420)), rep("LightON",nrow(df_180420)));
ei_180420 <- array("Excitatory", length(regions_180420));
ei_180420[c(4,5,6,7,9,11,12,13,15,16,17,18,19,20,28,29,30,31,41,42,49,51,52,53,59,62,64,68,69,70,71)] <- "Inhibitory";
dat_180420 <- data.frame(c(df_180420$Ongoing, df_180420$LightON), condition, nodeid, rep(regions_180420,2), ei_180420);
colnames(dat_180420) <- c("Strength", "condition", "nodeid", "region", "EI")
dat_180420$condition <- factor(dat_180420$condition, levels = c("Ongoing", "LightON"))
dat_180420$region <- factor(dat_180420$region)
dat_180420$Strength <- as.numeric(sub(",", ".", c(df_180420$Ongoing, df_180420$LightON), fixed = TRUE))
dat_180420$EI <- factor(ei_180420)
model_180420 <- lm(formula = "Strength ~ condition * region * EI + (1|nodeid)", data=dat_180420)
res_aov_180420 <- aov(model_180420)
summary(res_aov_180420)
cm=compare_means(Strength ~ condition, data = dat_180420, group.by = "region")

potout = boxplot.stats(dat_180420$Strength[dat_180420$condition == "Ongoing"])$out
potout_id = dat_180420$nodeid[which(dat_180420$Strength[dat_180420$condition == "Ongoing"] %in%  potout)]
potout_evoked = boxplot.stats(dat_180420$Strength[dat_180420$condition == "LightON"])$out
potout_id_evoked = dat_180420$nodeid[which(dat_180420$Strength[dat_180420$condition == "LightON"] %in%  potout_evoked)]
x = list(ongoing = potout_id, lightON = potout_id_evoked);
png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/hubness_180420.png")
ggvenn(x)
dev.off()



tmpdf_deg_wei = aggregate(Strength ~ region + condition + EI, data = dat_180420, mean);
#tmpdf_deg_wei$expID = rep("180131",4,length(levels(as.factor(tmpdf_deg_wei$region)))*2);
tmpdf_deg_wei$expID = rep("180420",1,nrow(tmpdf_deg_wei));
df_deg_wei = rbind(df_deg_wei, tmpdf_deg_wei);

###
df_180423 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180423_deg_wei")
regions_180423 <- character(nrow(df_180423));
regions_180423[c(1,2,3,4)] = "ventromedial hypothalamus"
#regions_180423[c(5)] = "submammillothalamic nucleus"
regions_180423[c(6,7,8,9,10,18,19,20,26,27,28,29,30,31,32,33,34,35,36)] = "dorsomedial hypothalamus"
regions_180423[c(11, 12,13,14,15,16,17,21,22,23,24,25,51,52,53,54,55,56)] = "posterior hypothalamus"
#regions_180423[c(37, 38,39,40,41,42,43,44,45,46,47,48,49,50)] = "premammillary nucleus"
regions_180423[c(37, 38,39,40,41,42,43,44,45,46,47,48,49,50,5)] = "mammillary complex"
nodeid <- rep(1:nrow(df_180423),2);
condition <- c(rep("Ongoing",nrow(df_180423)), rep("LightON",nrow(df_180423)));
ei_180423 <- array("Excitatory", length(regions_180423));
ei_180423[c(1,2,3,4,8,12,14,16,19,26,27,29,31,32,35,38,40,41,45,47,48,49,50,51,52,54)] <- "Inhibitory";
dat_180423 <- data.frame(c(df_180423$Ongoing, df_180423$LightON), condition, nodeid, rep(regions_180423,2), ei_180423);
colnames(dat_180423) <- c("Strength", "condition", "nodeid", "region", "EI")
dat_180423$condition <- factor(dat_180423$condition, levels = c("Ongoing", "LightON"))
dat_180423$region <- factor(dat_180423$region)
dat_180423$Strength <- as.numeric(sub(",", ".", c(df_180423$Ongoing, df_180423$LightON), fixed = TRUE))
dat_180423$EI <- factor(ei_180423);
model_180423 <- lm(formula = "Strength ~ condition * region * EI + (1|nodeid)", data=dat_180423)
res_aov_180423 <- aov(model_180423)
summary(res_aov_180423)
cm=compare_means(Strength ~ condition, data = dat_180423, group.by = "region")

potout = boxplot.stats(dat_180423$Strength[dat_180423$condition == "Ongoing"])$out
potout_id = dat_180423$nodeid[which(dat_180423$Strength[dat_180423$condition == "Ongoing"] %in%  potout)]
potout_evoked = boxplot.stats(dat_180423$Strength[dat_180423$condition == "LightON"])$out
potout_id_evoked = dat_180423$nodeid[which(dat_180423$Strength[dat_180423$condition == "LightON"] %in%  potout_evoked)]
x = list(ongoing = potout_id, lightON = potout_id_evoked);
png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/hubness_180423.png")
ggvenn(x)
dev.off()



tmpdf_deg_wei = aggregate(Strength ~ region + condition + EI, data = dat_180423, mean);
#tmpdf_deg_wei$expID = rep("180131",4,length(levels(as.factor(tmpdf_deg_wei$region)))*2);
tmpdf_deg_wei$expID = rep("180423",1,nrow(tmpdf_deg_wei));
df_deg_wei = rbind(df_deg_wei, tmpdf_deg_wei);

df_deg_wei$expID <- factor(df_deg_wei$expID)
model_deg_wei <- lm(formula = "Strength ~ condition * region * EI", data = df_deg_wei)
res_aov_deg_wei <- aov(model_deg_wei)
summary(res_aov_deg_wei)
posthoc <- TukeyHSD(res_aov_deg_wei)
cm=compare_means(Strength ~ condition, data = df_deg_wei, group.by = "region")
ttt <- posthoc$`region:EI`
which(ttt[,4] < 0.05)


df_deg_wei$phase[df_deg_wei$expID == "171019"] = rep("Day", length(df_deg_wei$expID == "171019"))
df_deg_wei$phase[df_deg_wei$expID == "180131"] = rep("Day", length(df_deg_wei$expID == "180131"))
df_deg_wei$phase[df_deg_wei$expID == "180302"] = rep("Day", length(df_deg_wei$expID == "180302"))
df_deg_wei$phase[df_deg_wei$expID == "180131"] = rep("Day", length(df_deg_wei$expID == "180131"))
df_deg_wei$phase[df_deg_wei$expID == "180419"] = rep("Day", length(df_deg_wei$expID == "180419"))
df_deg_wei$phase[df_deg_wei$expID == "180420"] = rep("Day", length(df_deg_wei$expID == "180420"))
df_deg_wei$phase[df_deg_wei$expID == "180423"] = rep("Day", length(df_deg_wei$expID == "180423"))
df_deg_wei$phase[is.na(df_deg_wei$phase)] = "Night";
df_deg_wei$phase = as.factor(df_deg_wei$phase);

model_deg_wei <- lm(formula = "Strength ~ condition * region * EI * phase", data = df_deg_wei)
res_aov_deg_wei <- aov(model_deg_wei)
summary(res_aov_deg_wei)
cm=compare_means(Strength ~ condition + region, data = df_deg_wei, group.by = "region")
posthoc <- TukeyHSD(res_aov_deg_wei)
cm=compare_means(Strength ~ condition, data = df_deg_wei, group.by = "region")
ttt <- posthoc$`EI:phase`
which(ttt[,4] < 0.05)

model_deg_wei <- lm(formula = "Strength ~ phase * EI", data = df_deg_wei)
res_aov_deg_wei <- aov(model_deg_wei)
summary(res_aov_deg_wei)


png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/region_deg_wei_phase.png", width = 1860, height = 1440, units = "px")
ggboxplot(df_deg_wei, x="region", y="Strength", color="EI", add = "jitter") + #stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("")  + ylab("Node Strength")
dev.off()











#### Clustering Coefficient
##################################################################################
#################################################################################
df_171207 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171207_clus_wei")
regions_171207 <- character(nrow(df_171207));
regions_171207[c(1,2,3, 12,13,28, 18,23,24,25,26,27,31)] = "ventromedial thalamus"
regions_171207[c(4,5,6,7,11,14,15,16,29,30)] = "paraventricular hypothalamus"
regions_171207[c(8,9,10,17)] = "anterior hypothalamus"
regions_171207[c(19,20,21,22)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_171207),2);
condition <- c(rep("Ongoing",nrow(df_171207)), rep("LightON",nrow(df_171207)));
ei_171207 <- array("Excitatory", length(regions_171207));
ei_171207[9] <- "Inhibitory";
dat_171207 <- data.frame(c(df_171207$Ongoing, df_171207$LightON), condition, nodeid, rep(regions_171207,2), ei_171207);
colnames(dat_171207) <- c("LocalEff", "condition", "nodeid", "region", "EI");
dat_171207$condition <- factor(dat_171207$condition, levels = c("Ongoing", "LightON"))
dat_171207$region <- factor(dat_171207$region)
dat_171207$LocalEff <- as.numeric(sub(",", ".", c(df_171207$Ongoing, df_171207$LightON), fixed = TRUE))
dat_171207$EI <- factor(ei_171207);
model_171207 <- lm(formula = "LocalEff ~ condition * region * EI + (1|nodeid)", data=dat_171207)
res_aov_171207 <- aov(model_171207)
summary(res_aov_171207)
cm=compare_means(LocalEff ~ condition, data = dat_171207, group.by = "region")

tmpdf_clus_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_171207, mean);
#tmpdf_clus_wei$expID = rep("171207",1,length(levels(as.factor(dat_171207)))*2);
tmpdf_clus_wei$expID = rep("171207",1,length(dat_171207)*2);
df_clus_wei = tmpdf_clus_wei;
#df_degree = rbind(df_degree, tmpdf_degree);

###
df_171019 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171019_clus_wei")
regions_171019 <- character(nrow(df_171019));
regions_171019[seq(1,114)] = "ventromedial thalamus"
regions_171019[c(115, 116)] = "paraventricular hypothalamus"
regions_171019[seq(117,126)] = "ventromedial thalamus"
nodeid <- rep(1:nrow(df_171019),2);
condition <- c(rep("Ongoing",nrow(df_171019)), rep("LightON",nrow(df_171019)));
ei_171019 <- array("Excitatory", length(regions_171019));
ei_171019[c(60,125)] <- "Inhibitory";
dat_171019 <- data.frame(c(df_171019$Ongoing, df_171019$LightON), condition, nodeid, rep(regions_171019,2), ei_171019);
colnames(dat_171019) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_171019$condition <- factor(dat_171019$condition, levels = c("Ongoing", "LightON"))
dat_171019$region <- factor(dat_171019$region)
dat_171019$LocalEff <- as.numeric(sub(",", ".", c(df_171019$Ongoing, df_171019$LightON), fixed = TRUE))
dat_171019$EI <- factor(ei_171019);
model_171019 <- lm(formula = "LocalEff ~ condition * region * EI + (1|nodeid)", data=dat_171019)
res_aov_171019 <- aov(model_171019)
summary(res_aov_171019)
cm=compare_means(LocalEff ~ condition, data = dat_171019, group.by = "region")

tmpdf_clus_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_171019, mean);
#tmpdf_clus_wei$expID = rep("171019",1,length(levels(as.factor(tmpdf_clus_wei$region)))*2);
tmpdf_clus_wei$expID = rep("171019",1,nrow(tmpdf_clus_wei));
df_clus_wei = rbind(df_clus_wei, tmpdf_clus_wei);


###
df_171208 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171208_clus_wei")
regions_171208 <- character(nrow(df_171208));
regions_171208[c(1,2,3,4,5,7,20,21,22,23,24,25,28)] = "ventromedial thalamus"
regions_171208[c(6, 8,9,11,13,16,17)] = "ventromedial hypothalamus"
regions_171208[c(10,12, 32)] = "paraventricular hypothalamus"
regions_171208[c(14,15,18,19)] = "dorsomedial hypothalamus"
regions_171208[c(26, 27, 29, 30, 31,33, 34)] = "anterior hypothalamus"
nodeid <- rep(1:nrow(df_171208),2);
condition <- c(rep("Ongoing",nrow(df_171208)), rep("LightON",nrow(df_171208)));
ei_171208 <- array("Excitatory", length(regions_171208));
dat_171208 <- data.frame(c(df_171208$Ongoing, df_171208$LightON), condition, nodeid, rep(regions_171208,2), ei_171208);
colnames(dat_171208) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_171208$condition <- factor(dat_171208$condition, levels = c("Ongoing", "LightON"))
dat_171208$region <- factor(dat_171208$region)
dat_171208$LocalEff <- as.numeric(sub(",", ".", c(df_171208$Ongoing, df_171208$LightON), fixed = TRUE))
dat_171208$EI <- factor(ei_171208);
model_171208 <- lm(formula = "LocalEff ~ condition * region + (1|nodeid)", data=dat_171208)
res_aov_171208 <- aov(model_171208)
summary(res_aov_171208)
cm=compare_means(LocalEff ~ condition, data = dat_171208, group.by = "region")

tmpdf_clus_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_171208, mean);
#tmpdf_clus_wei$expID = rep("171208",1,length(levels(as.factor(tmpdf_clus_wei$region)))*2);
tmpdf_clus_wei$expID = rep("171208",1,nrow(tmpdf_clus_wei));
df_clus_wei = rbind(df_clus_wei, tmpdf_clus_wei);



###
df_171213 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171213_clus_wei")
regions_171213 <- character(nrow(df_171213));
regions_171213[c(11,12,13,20,21,22,23,24,25,26,27,29,30,31,32)] = "ventromedial thalamus"
regions_171213[c(1,2,3,4,5,6,7,8,9,10,14,15,16,17,18,19,28,33)] = "paraventricular hypothalamus"
nodeid <- rep(1:nrow(df_171213),2);
condition <- c(rep("Ongoing",nrow(df_171213)), rep("LightON",nrow(df_171213)));
ei_171213 <- array("Excitatory", length(regions_171213));
dat_171213 <- data.frame(c(df_171213$Ongoing, df_171213$LightON), condition, nodeid, rep(regions_171213,2), ei_171213);
colnames(dat_171213) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_171213$condition <- factor(dat_171213$condition, levels = c("Ongoing", "LightON"))
dat_171213$region <- factor(dat_171213$region)
dat_171213$LocalEff <- as.numeric(sub(",", ".", c(df_171213$Ongoing, df_171213$LightON), fixed = TRUE))
dat_171213$EI <- factor(ei_171213);
model_171213 <- lm(formula = "LocalEff ~ condition * region + (1|nodeid)", data=dat_171213)
res_aov_171213 <- aov(model_171213)
summary(res_aov_171213)
cm=compare_means(LocalEff ~ condition, data = dat_171213, group.by = "region")

tmpdf_clus_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_171213, mean);
#tmpdf_clus_wei$expID = rep("171213",1,length(levels(as.factor(tmpdf_clus_wei$region)))*2);
tmpdf_clus_wei$expID = rep("171213",1,nrow(tmpdf_clus_wei));
df_clus_wei = rbind(df_clus_wei, tmpdf_clus_wei);


###
df_180110 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180110_clus_wei")
regions_180110 <- character(nrow(df_180110));
regions_180110[c(1,2,11,18,20,21,23)] = "ventromedial thalamus"
regions_180110[c(3,4,5,6,7,9,10, 19,22,24,25,26,27,28,29,30,31,32,33,34,35)] = "paraventricular hypothalamus"
regions_180110[c(8)] = "anterior hypothalamus"
regions_180110[c(12,13,14,15,16,17)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_180110),2);
condition <- c(rep("Ongoing",nrow(df_180110)), rep("LightON",nrow(df_180110)));
ei_180110 <- array("Excitatory", length(regions_180110));
ei_180110[16] <- "Inhibitory"
dat_180110 <- data.frame(c(df_180110$Ongoing, df_180110$LightON), condition, nodeid, rep(regions_180110,2), ei_180110);
colnames(dat_180110) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_180110$condition <- factor(dat_180110$condition, levels = c("Ongoing", "LightON"))
dat_180110$region <- factor(dat_180110$region)
dat_180110$LocalEff <- as.numeric(sub(",", ".", c(df_180110$Ongoing, df_180110$LightON), fixed = TRUE))
dat_180110$EI <- factor(ei_180110)
model_180110 <- lm(formula = "LocalEff ~ condition * region * EI + (1|nodeid)", data=dat_180110)
res_aov_180110 <- aov(model_180110)
summary(res_aov_180110)
cm=compare_means(LocalEff ~ condition, data = dat_180110, group.by = "region")

tmpdf_clus_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_180110, mean);
#tmpdf_clus_wei$expID = rep("180110",1,length(levels(as.factor(tmpdf_clus_wei$region)))*2);
tmpdf_clus_wei$expID = rep("180110",1,nrow(tmpdf_clus_wei));
df_clus_wei = rbind(df_clus_wei, tmpdf_clus_wei);

###
df_180111 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180111_clus_wei")
regions_180111 <- character(nrow(df_180111));
regions_180111[c(1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,22,23,24)] = "ventromedial thalamus"
regions_180111[c(14)] = "paraventricular hypothalamus"
regions_180111[c(19,20,21)] = "ventromedial hypothalamus"
nodeid <- rep(1:nrow(df_180111),2);
condition <- c(rep("Ongoing",nrow(df_180111)), rep("LightON",nrow(df_180111)));
ei_180111 <- array("Excitatory", length(regions_180111));
dat_180111 <- data.frame(c(df_180111$Ongoing, df_180111$LightON), condition, nodeid, rep(regions_180111,2), ei_180111);
colnames(dat_180111) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_180111$condition <- factor(dat_180111$condition, levels = c("Ongoing", "LightON"))
dat_180111$region <- factor(dat_180111$region)
dat_180111$LocalEff <- as.numeric(sub(",", ".", c(df_180111$Ongoing, df_180111$LightON), fixed = TRUE))
dat_180111$EI <- factor(ei_180111);
model_180111 <- lm(formula = "LocalEff ~ condition * region + (1|nodeid)", data=dat_180111)
res_aov_180111 <- aov(model_180111)
summary(res_aov_180111)
cm=compare_means(LocalEff ~ condition, data = dat_180111, group.by = "region")

tmpdf_clus_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_180111, mean);
#tmpdf_clus_wei$expID = rep("180111",1,length(levels(as.factor(tmpdf_clus_wei$region)))*2);
tmpdf_clus_wei$expID = rep("180111",1,nrow(tmpdf_clus_wei));
df_clus_wei = rbind(df_clus_wei, tmpdf_clus_wei);


###
df_180131 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180131_clus_wei")
regions_180131 <- character(nrow(df_180131));
regions_180131[c(1,2,4,5,6,14,15,16, 26)] = "ventromedial thalamus"
regions_180131[c(3,13,18,19,20,21,22,23,24, 27,29)] = "posterior hypothalamus"
regions_180131[c(7,8,9,10,11,12,25,28,30,31)] = "dorsomedial hypothalamus"
regions_180131[c(17)] = "paraventricular hypothalamus"
nodeid <- rep(1:nrow(df_180131),2);
condition <- c(rep("Ongoing",nrow(df_180131)), rep("LightON",nrow(df_180131)));
ei_180131 <- array("Excitatory", length(regions_180131));
ei_180131[c(5,6,8,9,12,17,19,23,24,25,26,28)] <- "Inhibitory";
dat_180131 <- data.frame(c(df_180131$Ongoing, df_180131$LightON), condition, nodeid, rep(regions_180131,2), ei_180131);
colnames(dat_180131) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_180131$condition <- factor(dat_180131$condition, levels = c("Ongoing", "LightON"))
dat_180131$region <- factor(dat_180131$region)
dat_180131$LocalEff <- as.numeric(sub(",", ".", c(df_180131$Ongoing, df_180131$LightON), fixed = TRUE))
dat_180131$EI <- factor(ei_180131);
model_180131 <- lm(formula = "LocalEff ~ condition * region * EI + (1|nodeid)", data=dat_180131)
res_aov_180131 <- aov(model_180131)
summary(res_aov_180131)
cm=compare_means(LocalEff ~ condition, data = dat_180131, group.by = "region")

tmpdf_clus_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_180131, mean);
#tmpdf_clus_wei$expID = rep("180131",1,length(levels(as.factor(tmpdf_clus_wei$region)))*2);
tmpdf_clus_wei$expID = rep("180131",1,nrow(tmpdf_clus_wei));
df_clus_wei = rbind(df_clus_wei, tmpdf_clus_wei);


###
df_180221 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180221_clus_wei")
regions_180221 <- character(nrow(df_180221));
regions_180221[c(1, 35, 36, 48,49, 50, 2,3,4,12,17,18, 21, 22, 26, 31, 32, 33, 39, 40, 41, 42, 43, 44, 45, 46, 47)] = "ventromedial thalamus"
regions_180221[c(23, 24, 25)] = "ventromedial hypothalamus"
regions_180221[c(5,6,7,8,9,10,11,13,14,15,16, 27, 28, 29, 30, 37, 38)] = "dorsomedial hypothalamus"
regions_180221[c(19, 20, 34, 51, 52, 53, 54)] = "paraventricular hypothalamus"
nodeid <- rep(1:nrow(df_180221),2);
condition <- c(rep("Ongoing",nrow(df_180221)), rep("LightON",nrow(df_180221)));
ei_180221 <- array("Excitatory", length(regions_180221));
ei_180221[c(1,2,5,7,8,9,10,14,15,19,20,23,25,28,30,32,33,34,37,42,43,44,51,52,53,54)] <- "Inhibitory";
dat_180221 <- data.frame(c(df_180221$Ongoing, df_180221$LightON), condition, nodeid, rep(regions_180221,2), ei_180221);
colnames(dat_180221) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_180221$condition <- factor(dat_180221$condition, levels = c("Ongoing", "LightON"))
dat_180221$region <- factor(dat_180221$region)
dat_180221$LocalEff <- as.numeric(sub(",", ".", c(df_180221$Ongoing, df_180221$LightON), fixed = TRUE))
dat_180221$EI <- factor(ei_180221);
model_180221 <- lm(formula = "LocalEff ~ condition * region *EI + (1|nodeid)", data=dat_180221)
res_aov_180221 <- aov(model_180221)
summary(res_aov_180221)
cm=compare_means(LocalEff ~ condition, data = dat_180221, group.by = "region")

tmpdf_clus_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_180221, mean);
#tmpdf_clus_wei$expID = rep("180131",1,length(levels(as.factor(tmpdf_clus_wei$region)))*2);
tmpdf_clus_wei$expID = rep("180221",1,nrow(tmpdf_clus_wei));
df_clus_wei = rbind(df_clus_wei, tmpdf_clus_wei);


###
df_180228 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180228_clus_wei")
regions_180228 <- character(nrow(df_180228));
regions_180228[c(1,4,5,6, 24, 25, 74)] = "ventromedial thalamus"
regions_180228[c(2, 11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23,26, 40, 42, 62, 63, 64)] = "posterior hypothalamus"
regions_180228[c(3,7,8,10, 19, 20, 27, 28, 29, 30, 31, 54, 56, 57, 58, 59, 60, 61, 65, 66,67,68,69, 73, 75, 76, 77,78,79)] = "dorsomedial hypothalamus"
regions_180228[c(9,46, 47, 48, 49, 50)] = "ventromedial hypothalamus"
regions_180228[c(32, 33, 34, 35, 36, 37, 38,39, 41)] = "arcuate hypothalamus"
regions_180228[c(43,44,45,51,52,53,55,70,71, 72)] = "zona incerta"
nodeid <- rep(1:nrow(df_180228),2);
condition <- c(rep("Ongoing",nrow(df_180228)), rep("LightON",nrow(df_180228)));
ei_180228 <- array("Excitatory", length(regions_180228));
ei_180228[c(1,4,14,16,17,18,19,20,21,23,24,25,28,29,30,32,33,34,35,36,39,40,46,47,50,51,52,54,56,59,60,61,62)] <- "Inhibitory";
dat_180228 <- data.frame(c(df_180228$Ongoing, df_180228$LightON), condition, nodeid, rep(regions_180228,2), ei_180228);
colnames(dat_180228) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_180228$condition <- factor(dat_180228$condition, levels = c("Ongoing", "LightON"))
dat_180228$region <- factor(dat_180228$region)
dat_180228$LocalEff <- as.numeric(sub(",", ".", c(df_180228$Ongoing, df_180228$LightON), fixed = TRUE))
dat_180228$EI <- factor(ei_180228);
model_180228 <- lm(formula = "LocalEff ~ condition * region * EI + (1|nodeid)", data=dat_180228)
res_aov_180228 <- aov(model_180228)
summary(res_aov_180228)
cm=compare_means(LocalEff ~ condition, data = dat_180228, group.by = "region")

tmpdf_clus_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_180228, mean);
#tmpdf_clus_wei$expID = rep("180131",1,length(levels(as.factor(tmpdf_clus_wei$region)))*2);
tmpdf_clus_wei$expID = rep("180228",1,nrow(tmpdf_clus_wei));
df_clus_wei = rbind(df_clus_wei, tmpdf_clus_wei);


####
df_180302 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180302_clus_wei")
regions_180302 <- character(nrow(df_180302));
#regions_180302[c(2,3,4,5,6, 13, 14, 15,30,36)] = "tuberomammillary nucleus"
regions_180302[c(2,3,4,5,6, 13, 14, 15,30,36, 9, 10, 41,42, 18, 19, 20, 21)] = "mammillary complex"
regions_180302[c(1,11,12,7, 8,16, 17, 23, 31, 37, 38,39,40, 47)] = "posterior hypothalamus"
#regions_180302[c(9, 10, 41,42)] = "premammillary nucleus"
#regions_180302[c(18, 19, 20, 21)] = "supramammillary nucleus"
regions_180302[c(26, 27, 28, 43,44,45,46)] = "arcuate hypothalamus"
regions_180302[c(22, 24, 25, 29)] = "paraventricular hypothalamus"
regions_180302[c(32,33,34,35)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_180302),2);
condition <- c(rep("Ongoing",nrow(df_180302)), rep("LightON",nrow(df_180302)));
ei_180302 <- array("Excitatory", length(regions_180302));
ei_180302[c(1,3,8,9,11,15,16,17,19,20,25,26,28,29,30,33,34,35,36,37,40,41,42,47)] <- "Inhibitory";
dat_180302 <- data.frame(c(df_180302$Ongoing, df_180302$LightON), condition, nodeid, rep(regions_180302,2), ei_180302);
colnames(dat_180302) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_180302$condition <- factor(dat_180302$condition, levels = c("Ongoing", "LightON"))
dat_180302$region <- factor(dat_180302$region)
dat_180302$LocalEff <- as.numeric(sub(",", ".", c(df_180302$Ongoing, df_180302$LightON), fixed = TRUE))
dat_180302$EI <- factor(ei_180302)
model_180302 <- lm(formula = "LocalEff ~ condition * region *EI + (1|nodeid)", data=dat_180302)
res_aov_180302 <- aov(model_180302)
summary(res_aov_180302)
cm=compare_means(LocalEff ~ condition, data = dat_180302, group.by = "region")

tmpdf_clus_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_180302, mean);
#tmpdf_clus_wei$expID = rep("180131",4,length(levels(as.factor(tmpdf_clus_wei$region)))*2);
tmpdf_clus_wei$expID = rep("180302",1,nrow(tmpdf_clus_wei));
df_clus_wei = rbind(df_clus_wei, tmpdf_clus_wei);

###
df_180419 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180419_clus_wei")
regions_180419 <- character(nrow(df_180419));
regions_180419[c(1,2,11, 12,9, 10)] = "ventromedial hypothalamus"
regions_180419[c(4,5,6,7,8)] = "ventromedial thalamus"
regions_180419[c(36,37,38,39,51,52,53,54,55,56,57,58,59,60,61)] = "posterior hypothalamus"
regions_180419[c(3, 40, 41, 42, 43, 44,45,46,47,48,49,50, 13,14,15,16,17,18,23, 19,20,21,22,24,25,26,27,28,29,30,31,32,33,34,35)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_180419),2);
condition <- c(rep("Ongoing",nrow(df_180419)), rep("LightON",nrow(df_180419)));
ei_180419 <- array("Excitatory", length(regions_180419));
ei_180419[c(1,3,4,5,6,7,9,10,11,15,16,17,19,20,22,33,34,41,42,43,44,47,48,55,58,59)] <- "Inhibitory";
dat_180419 <- data.frame(c(df_180419$Ongoing, df_180419$LightON), condition, nodeid, rep(regions_180419,2), ei_180419);
colnames(dat_180419) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_180419$condition <- factor(dat_180419$condition, levels = c("Ongoing", "LightON"))
dat_180419$region <- factor(dat_180419$region)
dat_180419$LocalEff <- as.numeric(sub(",", ".", c(df_180419$Ongoing, df_180419$LightON), fixed = TRUE))
dat_180419$EI <- factor(ei_180419)
model_180419 <- lm(formula = "LocalEff ~ condition * region *EI + (1|nodeid)", data=dat_180419)
res_aov_180419 <- aov(model_180419)
summary(res_aov_180419)
cm=compare_means(LocalEff ~ condition, data = dat_180419, group.by = "region")

tmpdf_clus_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_180419, mean);
#tmpdf_clus_wei$expID = rep("180131",4,length(levels(as.factor(tmpdf_clus_wei$region)))*2);
tmpdf_clus_wei$expID = rep("180419",1,nrow(tmpdf_clus_wei));
df_clus_wei = rbind(df_clus_wei, tmpdf_clus_wei);

####
df_180420 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180420_clus_wei")
regions_180420 <- character(nrow(df_180420));
#regions_180420[c(51,52,53,54,55,56,57,58,59,60,61,62,63)] = "supramammillary nucleus"
regions_180420[c(51,52,53,54,55,56,57,58,59,60,61,62,63,19, 20, 21, 22, 23, 24, 25,26,27,28,29,30,31)] = "mammillary complex"
regions_180420[c(32, 33, 34, 35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50)] = "paraventricular hypothalamus"
#regions_180420[c(19, 20, 21, 22, 23, 24, 25,26,27,28,29,30,31)] = "premammillary nucleus"
regions_180420[c(2, 3, 4, 5, 6,7,8,9,10,11,12,13,14,15,16,17,18,64,65,66,67,68,69,70,71)] = "posterior hypothalamus"
regions_180420[c(1)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_180420),2);
condition <- c(rep("Ongoing",nrow(df_180420)), rep("LightON",nrow(df_180420)));
ei_180420 <- array("Excitatory", length(regions_180420));
ei_180420[c(4,5,6,7,9,11,12,13,15,16,17,18,19,20,28,29,30,31,41,42,49,51,52,53,59,62,64,68,69,70,71)] <- "Inhibitory";
dat_180420 <- data.frame(c(df_180420$Ongoing, df_180420$LightON), condition, nodeid, rep(regions_180420,2), ei_180420);
colnames(dat_180420) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_180420$condition <- factor(dat_180420$condition, levels = c("Ongoing", "LightON"))
dat_180420$region <- factor(dat_180420$region)
dat_180420$LocalEff <- as.numeric(sub(",", ".", c(df_180420$Ongoing, df_180420$LightON), fixed = TRUE))
dat_180420$EI <- factor(ei_180420)
model_180420 <- lm(formula = "LocalEff ~ condition * region * EI + (1|nodeid)", data=dat_180420)
res_aov_180420 <- aov(model_180420)
summary(res_aov_180420)
cm=compare_means(LocalEff ~ condition, data = dat_180420, group.by = "region")

tmpdf_clus_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_180420, mean);
#tmpdf_clus_wei$expID = rep("180131",4,length(levels(as.factor(tmpdf_clus_wei$region)))*2);
tmpdf_clus_wei$expID = rep("180420",1,nrow(tmpdf_clus_wei));
df_clus_wei = rbind(df_clus_wei, tmpdf_clus_wei);

###
df_180423 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180423_clus_wei")
regions_180423 <- character(nrow(df_180423));
regions_180423[c(1,2,3,4)] = "ventromedial hypothalamus"
#regions_180423[c(5)] = "submammillothalamic nucleus"
regions_180423[c(6,7,8,9,10,18,19,20,26,27,28,29,30,31,32,33,34,35,36)] = "dorsomedial hypothalamus"
regions_180423[c(11, 12,13,14,15,16,17,21,22,23,24,25,51,52,53,54,55,56)] = "posterior hypothalamus"
#regions_180423[c(37, 38,39,40,41,42,43,44,45,46,47,48,49,50)] = "premammillary nucleus"
regions_180423[c(37, 38,39,40,41,42,43,44,45,46,47,48,49,50,5)] = "mammillary complex"
nodeid <- rep(1:nrow(df_180423),2);
condition <- c(rep("Ongoing",nrow(df_180423)), rep("LightON",nrow(df_180423)));
ei_180423 <- array("Excitatory", length(regions_180423));
ei_180423[c(1,2,3,4,8,12,14,16,19,26,27,29,31,32,35,38,40,41,45,47,48,49,50,51,52,54)] <- "Inhibitory";
dat_180423 <- data.frame(c(df_180423$Ongoing, df_180423$LightON), condition, nodeid, rep(regions_180423,2), ei_180423);
colnames(dat_180423) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_180423$condition <- factor(dat_180423$condition, levels = c("Ongoing", "LightON"))
dat_180423$region <- factor(dat_180423$region)
dat_180423$LocalEff <- as.numeric(sub(",", ".", c(df_180423$Ongoing, df_180423$LightON), fixed = TRUE))
dat_180423$EI <- factor(ei_180423);
model_180423 <- lm(formula = "LocalEff ~ condition * region * EI + (1|nodeid)", data=dat_180423)
res_aov_180423 <- aov(model_180423)
summary(res_aov_180423)
cm=compare_means(LocalEff ~ condition, data = dat_180423, group.by = "region")

tmpdf_clus_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_180423, mean);
#tmpdf_clus_wei$expID = rep("180131",4,length(levels(as.factor(tmpdf_clus_wei$region)))*2);
tmpdf_clus_wei$expID = rep("180423",1,nrow(tmpdf_clus_wei));
df_clus_wei = rbind(df_clus_wei, tmpdf_clus_wei);

df_clus_wei$expID <- factor(df_clus_wei$expID)
model_clus_wei <- lm(formula = "LocalEff ~ condition * region * EI", data = df_clus_wei)
res_aov_clus_wei <- aov(model_clus_wei)
summary(res_aov_clus_wei)
posthoc <- TukeyHSD(res_aov_clus_wei)
cm=compare_means(LocalEff ~ condition, data = df_clus_wei, group.by = "region")
ttt <- posthoc$`region:EI`
which(ttt[,4] < 0.05)


df_clus_wei$phase[df_clus_wei$expID == "171019"] = rep("Day", length(df_clus_wei$expID == "171019"))
df_clus_wei$phase[df_clus_wei$expID == "180131"] = rep("Day", length(df_clus_wei$expID == "180131"))
df_clus_wei$phase[df_clus_wei$expID == "180302"] = rep("Day", length(df_clus_wei$expID == "180302"))
df_clus_wei$phase[df_clus_wei$expID == "180131"] = rep("Day", length(df_clus_wei$expID == "180131"))
df_clus_wei$phase[df_clus_wei$expID == "180419"] = rep("Day", length(df_clus_wei$expID == "180419"))
df_clus_wei$phase[df_clus_wei$expID == "180420"] = rep("Day", length(df_clus_wei$expID == "180420"))
df_clus_wei$phase[df_clus_wei$expID == "180423"] = rep("Day", length(df_clus_wei$expID == "180423"))
df_clus_wei$phase[is.na(df_clus_wei$phase)] = "Night";
df_clus_wei$phase = as.factor(df_clus_wei$phase);

model_clus_wei <- lm(formula = "LocalEff ~ condition * region * EI * phase", data = df_clus_wei)
res_aov_clus_wei <- aov(model_clus_wei)
summary(res_aov_clus_wei)
cm=compare_means(LocalEff ~ condition + region, data = df_clus_wei, group.by = "region")
posthoc <- TukeyHSD(res_aov_clus_wei)
cm=compare_means(LocalEff ~ condition, data = df_clus_wei, group.by = "region")
ttt <- posthoc$`phase`
which(ttt[,4] < 0.05)

model_clus_wei <- lm(formula = "LocalEff ~ phase * EI", data = df_clus_wei)
res_aov_clus_wei <- aov(model_clus_wei)
summary(res_aov_clus_wei)


png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/region_clus_wei_phase.png", width = 1860, height = 1440, units = "px")
ggboxplot(df_clus_wei, x="region", y="ClusCoef", color="EI", add = "jitter") + #stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("")  + ylab("Node ClusCoeff")
dev.off()








#### Local Efficiency
##################################################################################
#################################################################################
df_171207 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171207_localeff_wei")
regions_171207 <- character(nrow(df_171207));
regions_171207[c(1,2,3, 12,13,28, 18,23,24,25,26,27,31)] = "ventromedial thalamus"
regions_171207[c(4,5,6,7,11,14,15,16,29,30)] = "paraventricular hypothalamus"
regions_171207[c(8,9,10,17)] = "anterior hypothalamus"
regions_171207[c(19,20,21,22)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_171207),2);
condition <- c(rep("Ongoing",nrow(df_171207)), rep("LightON",nrow(df_171207)));
ei_171207 <- array("Excitatory", length(regions_171207));
ei_171207[9] <- "Inhibitory";
dat_171207 <- data.frame(c(df_171207$Ongoing, df_171207$LightON), condition, nodeid, rep(regions_171207,2), ei_171207);
colnames(dat_171207) <- c("LocalEff", "condition", "nodeid", "region", "EI");
dat_171207$condition <- factor(dat_171207$condition, levels = c("Ongoing", "LightON"))
dat_171207$region <- factor(dat_171207$region)
dat_171207$LocalEff <- as.numeric(sub(",", ".", c(df_171207$Ongoing, df_171207$LightON), fixed = TRUE))
dat_171207$EI <- factor(ei_171207);
model_171207 <- lm(formula = "LocalEff ~ condition * region * EI + (1|nodeid)", data=dat_171207)
res_aov_171207 <- aov(model_171207)
summary(res_aov_171207)
cm=compare_means(LocalEff ~ condition, data = dat_171207, group.by = "region")

tmpdf_localeff_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_171207, mean);
#tmpdf_localeff_wei$expID = rep("171207",1,length(levels(as.factor(dat_171207)))*2);
tmpdf_localeff_wei$expID = rep("171207",1,length(dat_171207)*2);
df_localeff_wei = tmpdf_localeff_wei;
#df_degree = rbind(df_degree, tmpdf_degree);

###
df_171019 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171019_localeff_wei")
regions_171019 <- character(nrow(df_171019));
regions_171019[seq(1,114)] = "ventromedial thalamus"
regions_171019[c(115, 116)] = "paraventricular hypothalamus"
regions_171019[seq(117,126)] = "ventromedial thalamus"
nodeid <- rep(1:nrow(df_171019),2);
condition <- c(rep("Ongoing",nrow(df_171019)), rep("LightON",nrow(df_171019)));
ei_171019 <- array("Excitatory", length(regions_171019));
ei_171019[c(60,125)] <- "Inhibitory";
dat_171019 <- data.frame(c(df_171019$Ongoing, df_171019$LightON), condition, nodeid, rep(regions_171019,2), ei_171019);
colnames(dat_171019) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_171019$condition <- factor(dat_171019$condition, levels = c("Ongoing", "LightON"))
dat_171019$region <- factor(dat_171019$region)
dat_171019$LocalEff <- as.numeric(sub(",", ".", c(df_171019$Ongoing, df_171019$LightON), fixed = TRUE))
dat_171019$EI <- factor(ei_171019);
model_171019 <- lm(formula = "LocalEff ~ condition * region * EI + (1|nodeid)", data=dat_171019)
res_aov_171019 <- aov(model_171019)
summary(res_aov_171019)
cm=compare_means(LocalEff ~ condition, data = dat_171019, group.by = "region")

tmpdf_localeff_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_171019, mean);
#tmpdf_localeff_wei$expID = rep("171019",1,length(levels(as.factor(tmpdf_localeff_wei$region)))*2);
tmpdf_localeff_wei$expID = rep("171019",1,nrow(tmpdf_localeff_wei));
df_localeff_wei = rbind(df_localeff_wei, tmpdf_localeff_wei);


###
df_171208 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171208_localeff_wei")
regions_171208 <- character(nrow(df_171208));
regions_171208[c(1,2,3,4,5,7,20,21,22,23,24,25,28)] = "ventromedial thalamus"
regions_171208[c(6, 8,9,11,13,16,17)] = "ventromedial hypothalamus"
regions_171208[c(10,12, 32)] = "paraventricular hypothalamus"
regions_171208[c(14,15,18,19)] = "dorsomedial hypothalamus"
regions_171208[c(26, 27, 29, 30, 31,33, 34)] = "anterior hypothalamus"
nodeid <- rep(1:nrow(df_171208),2);
condition <- c(rep("Ongoing",nrow(df_171208)), rep("LightON",nrow(df_171208)));
ei_171208 <- array("Excitatory", length(regions_171208));
dat_171208 <- data.frame(c(df_171208$Ongoing, df_171208$LightON), condition, nodeid, rep(regions_171208,2), ei_171208);
colnames(dat_171208) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_171208$condition <- factor(dat_171208$condition, levels = c("Ongoing", "LightON"))
dat_171208$region <- factor(dat_171208$region)
dat_171208$LocalEff <- as.numeric(sub(",", ".", c(df_171208$Ongoing, df_171208$LightON), fixed = TRUE))
dat_171208$EI <- factor(ei_171208);
model_171208 <- lm(formula = "LocalEff ~ condition * region + (1|nodeid)", data=dat_171208)
res_aov_171208 <- aov(model_171208)
summary(res_aov_171208)
cm=compare_means(LocalEff ~ condition, data = dat_171208, group.by = "region")

tmpdf_localeff_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_171208, mean);
#tmpdf_localeff_wei$expID = rep("171208",1,length(levels(as.factor(tmpdf_localeff_wei$region)))*2);
tmpdf_localeff_wei$expID = rep("171208",1,nrow(tmpdf_localeff_wei));
df_localeff_wei = rbind(df_localeff_wei, tmpdf_localeff_wei);



###
df_171213 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="171213_localeff_wei")
regions_171213 <- character(nrow(df_171213));
regions_171213[c(11,12,13,20,21,22,23,24,25,26,27,29,30,31,32)] = "ventromedial thalamus"
regions_171213[c(1,2,3,4,5,6,7,8,9,10,14,15,16,17,18,19,28,33)] = "paraventricular hypothalamus"
nodeid <- rep(1:nrow(df_171213),2);
condition <- c(rep("Ongoing",nrow(df_171213)), rep("LightON",nrow(df_171213)));
ei_171213 <- array("Excitatory", length(regions_171213));
dat_171213 <- data.frame(c(df_171213$Ongoing, df_171213$LightON), condition, nodeid, rep(regions_171213,2), ei_171213);
colnames(dat_171213) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_171213$condition <- factor(dat_171213$condition, levels = c("Ongoing", "LightON"))
dat_171213$region <- factor(dat_171213$region)
dat_171213$LocalEff <- as.numeric(sub(",", ".", c(df_171213$Ongoing, df_171213$LightON), fixed = TRUE))
dat_171213$EI <- factor(ei_171213);
model_171213 <- lm(formula = "LocalEff ~ condition * region + (1|nodeid)", data=dat_171213)
res_aov_171213 <- aov(model_171213)
summary(res_aov_171213)
cm=compare_means(LocalEff ~ condition, data = dat_171213, group.by = "region")

tmpdf_localeff_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_171213, mean);
#tmpdf_localeff_wei$expID = rep("171213",1,length(levels(as.factor(tmpdf_localeff_wei$region)))*2);
tmpdf_localeff_wei$expID = rep("171213",1,nrow(tmpdf_localeff_wei));
df_localeff_wei = rbind(df_localeff_wei, tmpdf_localeff_wei);


###
df_180110 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180110_localeff_wei")
regions_180110 <- character(nrow(df_180110));
regions_180110[c(1,2,11,18,20,21,23)] = "ventromedial thalamus"
regions_180110[c(3,4,5,6,7,9,10, 19,22,24,25,26,27,28,29,30,31,32,33,34,35)] = "paraventricular hypothalamus"
regions_180110[c(8)] = "anterior hypothalamus"
regions_180110[c(12,13,14,15,16,17)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_180110),2);
condition <- c(rep("Ongoing",nrow(df_180110)), rep("LightON",nrow(df_180110)));
ei_180110 <- array("Excitatory", length(regions_180110));
ei_180110[16] <- "Inhibitory"
dat_180110 <- data.frame(c(df_180110$Ongoing, df_180110$LightON), condition, nodeid, rep(regions_180110,2), ei_180110);
colnames(dat_180110) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_180110$condition <- factor(dat_180110$condition, levels = c("Ongoing", "LightON"))
dat_180110$region <- factor(dat_180110$region)
dat_180110$LocalEff <- as.numeric(sub(",", ".", c(df_180110$Ongoing, df_180110$LightON), fixed = TRUE))
dat_180110$EI <- factor(ei_180110)
model_180110 <- lm(formula = "LocalEff ~ condition * region * EI + (1|nodeid)", data=dat_180110)
res_aov_180110 <- aov(model_180110)
summary(res_aov_180110)
cm=compare_means(LocalEff ~ condition, data = dat_180110, group.by = "region")

tmpdf_localeff_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_180110, mean);
#tmpdf_localeff_wei$expID = rep("180110",1,length(levels(as.factor(tmpdf_localeff_wei$region)))*2);
tmpdf_localeff_wei$expID = rep("180110",1,nrow(tmpdf_localeff_wei));
df_localeff_wei = rbind(df_localeff_wei, tmpdf_localeff_wei);

###
df_180111 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180111_localeff_wei")
regions_180111 <- character(nrow(df_180111));
regions_180111[c(1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,22,23,24)] = "ventromedial thalamus"
regions_180111[c(14)] = "paraventricular hypothalamus"
regions_180111[c(19,20,21)] = "ventromedial hypothalamus"
nodeid <- rep(1:nrow(df_180111),2);
condition <- c(rep("Ongoing",nrow(df_180111)), rep("LightON",nrow(df_180111)));
ei_180111 <- array("Excitatory", length(regions_180111));
dat_180111 <- data.frame(c(df_180111$Ongoing, df_180111$LightON), condition, nodeid, rep(regions_180111,2), ei_180111);
colnames(dat_180111) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_180111$condition <- factor(dat_180111$condition, levels = c("Ongoing", "LightON"))
dat_180111$region <- factor(dat_180111$region)
dat_180111$LocalEff <- as.numeric(sub(",", ".", c(df_180111$Ongoing, df_180111$LightON), fixed = TRUE))
dat_180111$EI <- factor(ei_180111);
model_180111 <- lm(formula = "LocalEff ~ condition * region + (1|nodeid)", data=dat_180111)
res_aov_180111 <- aov(model_180111)
summary(res_aov_180111)
cm=compare_means(LocalEff ~ condition, data = dat_180111, group.by = "region")

tmpdf_localeff_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_180111, mean);
#tmpdf_localeff_wei$expID = rep("180111",1,length(levels(as.factor(tmpdf_localeff_wei$region)))*2);
tmpdf_localeff_wei$expID = rep("180111",1,nrow(tmpdf_localeff_wei));
df_localeff_wei = rbind(df_localeff_wei, tmpdf_localeff_wei);


###
df_180131 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180131_localeff_wei")
regions_180131 <- character(nrow(df_180131));
regions_180131[c(1,2,4,5,6,14,15,16, 26)] = "ventromedial thalamus"
regions_180131[c(3,13,18,19,20,21,22,23,24, 27,29)] = "posterior hypothalamus"
regions_180131[c(7,8,9,10,11,12,25,28,30,31)] = "dorsomedial hypothalamus"
regions_180131[c(17)] = "paraventricular hypothalamus"
nodeid <- rep(1:nrow(df_180131),2);
condition <- c(rep("Ongoing",nrow(df_180131)), rep("LightON",nrow(df_180131)));
ei_180131 <- array("Excitatory", length(regions_180131));
ei_180131[c(5,6,8,9,12,17,19,23,24,25,26,28)] <- "Inhibitory";
dat_180131 <- data.frame(c(df_180131$Ongoing, df_180131$LightON), condition, nodeid, rep(regions_180131,2), ei_180131);
colnames(dat_180131) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_180131$condition <- factor(dat_180131$condition, levels = c("Ongoing", "LightON"))
dat_180131$region <- factor(dat_180131$region)
dat_180131$LocalEff <- as.numeric(sub(",", ".", c(df_180131$Ongoing, df_180131$LightON), fixed = TRUE))
dat_180131$EI <- factor(ei_180131);
model_180131 <- lm(formula = "LocalEff ~ condition * region * EI + (1|nodeid)", data=dat_180131)
res_aov_180131 <- aov(model_180131)
summary(res_aov_180131)
cm=compare_means(LocalEff ~ condition, data = dat_180131, group.by = "region")

tmpdf_localeff_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_180131, mean);
#tmpdf_localeff_wei$expID = rep("180131",1,length(levels(as.factor(tmpdf_localeff_wei$region)))*2);
tmpdf_localeff_wei$expID = rep("180131",1,nrow(tmpdf_localeff_wei));
df_localeff_wei = rbind(df_localeff_wei, tmpdf_localeff_wei);


###
df_180221 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180221_localeff_wei")
regions_180221 <- character(nrow(df_180221));
regions_180221[c(1, 35, 36, 48,49, 50, 2,3,4,12,17,18, 21, 22, 26, 31, 32, 33, 39, 40, 41, 42, 43, 44, 45, 46, 47)] = "ventromedial thalamus"
regions_180221[c(23, 24, 25)] = "ventromedial hypothalamus"
regions_180221[c(5,6,7,8,9,10,11,13,14,15,16, 27, 28, 29, 30, 37, 38)] = "dorsomedial hypothalamus"
regions_180221[c(19, 20, 34, 51, 52, 53, 54)] = "paraventricular hypothalamus"
nodeid <- rep(1:nrow(df_180221),2);
condition <- c(rep("Ongoing",nrow(df_180221)), rep("LightON",nrow(df_180221)));
ei_180221 <- array("Excitatory", length(regions_180221));
ei_180221[c(1,2,5,7,8,9,10,14,15,19,20,23,25,28,30,32,33,34,37,42,43,44,51,52,53,54)] <- "Inhibitory";
dat_180221 <- data.frame(c(df_180221$Ongoing, df_180221$LightON), condition, nodeid, rep(regions_180221,2), ei_180221);
colnames(dat_180221) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_180221$condition <- factor(dat_180221$condition, levels = c("Ongoing", "LightON"))
dat_180221$region <- factor(dat_180221$region)
dat_180221$LocalEff <- as.numeric(sub(",", ".", c(df_180221$Ongoing, df_180221$LightON), fixed = TRUE))
dat_180221$EI <- factor(ei_180221);
model_180221 <- lm(formula = "LocalEff ~ condition * region *EI + (1|nodeid)", data=dat_180221)
res_aov_180221 <- aov(model_180221)
summary(res_aov_180221)
cm=compare_means(LocalEff ~ condition, data = dat_180221, group.by = "region")

tmpdf_localeff_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_180221, mean);
#tmpdf_localeff_wei$expID = rep("180131",1,length(levels(as.factor(tmpdf_localeff_wei$region)))*2);
tmpdf_localeff_wei$expID = rep("180221",1,nrow(tmpdf_localeff_wei));
df_localeff_wei = rbind(df_localeff_wei, tmpdf_localeff_wei);


###
df_180228 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180228_localeff_wei")
regions_180228 <- character(nrow(df_180228));
regions_180228[c(1,4,5,6, 24, 25, 74)] = "ventromedial thalamus"
regions_180228[c(2, 11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23,26, 40, 42, 62, 63, 64)] = "posterior hypothalamus"
regions_180228[c(3,7,8,10, 19, 20, 27, 28, 29, 30, 31, 54, 56, 57, 58, 59, 60, 61, 65, 66,67,68,69, 73, 75, 76, 77,78,79)] = "dorsomedial hypothalamus"
regions_180228[c(9,46, 47, 48, 49, 50)] = "ventromedial hypothalamus"
regions_180228[c(32, 33, 34, 35, 36, 37, 38,39, 41)] = "arcuate hypothalamus"
regions_180228[c(43,44,45,51,52,53,55,70,71, 72)] = "zona incerta"
nodeid <- rep(1:nrow(df_180228),2);
condition <- c(rep("Ongoing",nrow(df_180228)), rep("LightON",nrow(df_180228)));
ei_180228 <- array("Excitatory", length(regions_180228));
ei_180228[c(1,4,14,16,17,18,19,20,21,23,24,25,28,29,30,32,33,34,35,36,39,40,46,47,50,51,52,54,56,59,60,61,62)] <- "Inhibitory";
dat_180228 <- data.frame(c(df_180228$Ongoing, df_180228$LightON), condition, nodeid, rep(regions_180228,2), ei_180228);
colnames(dat_180228) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_180228$condition <- factor(dat_180228$condition, levels = c("Ongoing", "LightON"))
dat_180228$region <- factor(dat_180228$region)
dat_180228$LocalEff <- as.numeric(sub(",", ".", c(df_180228$Ongoing, df_180228$LightON), fixed = TRUE))
dat_180228$EI <- factor(ei_180228);
model_180228 <- lm(formula = "LocalEff ~ condition * region * EI + (1|nodeid)", data=dat_180228)
res_aov_180228 <- aov(model_180228)
summary(res_aov_180228)
cm=compare_means(LocalEff ~ condition, data = dat_180228, group.by = "region")

tmpdf_localeff_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_180228, mean);
#tmpdf_localeff_wei$expID = rep("180131",1,length(levels(as.factor(tmpdf_localeff_wei$region)))*2);
tmpdf_localeff_wei$expID = rep("180228",1,nrow(tmpdf_localeff_wei));
df_localeff_wei = rbind(df_localeff_wei, tmpdf_localeff_wei);


####
df_180302 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180302_localeff_wei")
regions_180302 <- character(nrow(df_180302));
#regions_180302[c(2,3,4,5,6, 13, 14, 15,30,36)] = "tuberomammillary nucleus"
regions_180302[c(2,3,4,5,6, 13, 14, 15,30,36, 9, 10, 41,42, 18, 19, 20, 21)] = "mammillary complex"
regions_180302[c(1,11,12,7, 8,16, 17, 23, 31, 37, 38,39,40, 47)] = "posterior hypothalamus"
#regions_180302[c(9, 10, 41,42)] = "premammillary nucleus"
#regions_180302[c(18, 19, 20, 21)] = "supramammillary nucleus"
regions_180302[c(26, 27, 28, 43,44,45,46)] = "arcuate hypothalamus"
regions_180302[c(22, 24, 25, 29)] = "paraventricular hypothalamus"
regions_180302[c(32,33,34,35)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_180302),2);
condition <- c(rep("Ongoing",nrow(df_180302)), rep("LightON",nrow(df_180302)));
ei_180302 <- array("Excitatory", length(regions_180302));
ei_180302[c(1,3,8,9,11,15,16,17,19,20,25,26,28,29,30,33,34,35,36,37,40,41,42,47)] <- "Inhibitory";
dat_180302 <- data.frame(c(df_180302$Ongoing, df_180302$LightON), condition, nodeid, rep(regions_180302,2), ei_180302);
colnames(dat_180302) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_180302$condition <- factor(dat_180302$condition, levels = c("Ongoing", "LightON"))
dat_180302$region <- factor(dat_180302$region)
dat_180302$LocalEff <- as.numeric(sub(",", ".", c(df_180302$Ongoing, df_180302$LightON), fixed = TRUE))
dat_180302$EI <- factor(ei_180302)
model_180302 <- lm(formula = "LocalEff ~ condition * region *EI + (1|nodeid)", data=dat_180302)
res_aov_180302 <- aov(model_180302)
summary(res_aov_180302)
cm=compare_means(LocalEff ~ condition, data = dat_180302, group.by = "region")

tmpdf_localeff_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_180302, mean);
#tmpdf_localeff_wei$expID = rep("180131",4,length(levels(as.factor(tmpdf_localeff_wei$region)))*2);
tmpdf_localeff_wei$expID = rep("180302",1,nrow(tmpdf_localeff_wei));
df_localeff_wei = rbind(df_localeff_wei, tmpdf_localeff_wei);

###
df_180419 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180419_localeff_wei")
regions_180419 <- character(nrow(df_180419));
regions_180419[c(1,2,11, 12,9, 10)] = "ventromedial hypothalamus"
regions_180419[c(4,5,6,7,8)] = "ventromedial thalamus"
regions_180419[c(36,37,38,39,51,52,53,54,55,56,57,58,59,60,61)] = "posterior hypothalamus"
regions_180419[c(3, 40, 41, 42, 43, 44,45,46,47,48,49,50, 13,14,15,16,17,18,23, 19,20,21,22,24,25,26,27,28,29,30,31,32,33,34,35)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_180419),2);
condition <- c(rep("Ongoing",nrow(df_180419)), rep("LightON",nrow(df_180419)));
ei_180419 <- array("Excitatory", length(regions_180419));
ei_180419[c(1,3,4,5,6,7,9,10,11,15,16,17,19,20,22,33,34,41,42,43,44,47,48,55,58,59)] <- "Inhibitory";
dat_180419 <- data.frame(c(df_180419$Ongoing, df_180419$LightON), condition, nodeid, rep(regions_180419,2), ei_180419);
colnames(dat_180419) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_180419$condition <- factor(dat_180419$condition, levels = c("Ongoing", "LightON"))
dat_180419$region <- factor(dat_180419$region)
dat_180419$LocalEff <- as.numeric(sub(",", ".", c(df_180419$Ongoing, df_180419$LightON), fixed = TRUE))
dat_180419$EI <- factor(ei_180419)
model_180419 <- lm(formula = "LocalEff ~ condition * region *EI + (1|nodeid)", data=dat_180419)
res_aov_180419 <- aov(model_180419)
summary(res_aov_180419)
cm=compare_means(LocalEff ~ condition, data = dat_180419, group.by = "region")

tmpdf_localeff_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_180419, mean);
#tmpdf_localeff_wei$expID = rep("180131",4,length(levels(as.factor(tmpdf_localeff_wei$region)))*2);
tmpdf_localeff_wei$expID = rep("180419",1,nrow(tmpdf_localeff_wei));
df_localeff_wei = rbind(df_localeff_wei, tmpdf_localeff_wei);

####
df_180420 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180420_localeff_wei")
regions_180420 <- character(nrow(df_180420));
#regions_180420[c(51,52,53,54,55,56,57,58,59,60,61,62,63)] = "supramammillary nucleus"
regions_180420[c(51,52,53,54,55,56,57,58,59,60,61,62,63,19, 20, 21, 22, 23, 24, 25,26,27,28,29,30,31)] = "mammillary complex"
regions_180420[c(32, 33, 34, 35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50)] = "paraventricular hypothalamus"
#regions_180420[c(19, 20, 21, 22, 23, 24, 25,26,27,28,29,30,31)] = "premammillary nucleus"
regions_180420[c(2, 3, 4, 5, 6,7,8,9,10,11,12,13,14,15,16,17,18,64,65,66,67,68,69,70,71)] = "posterior hypothalamus"
regions_180420[c(1)] = "dorsomedial hypothalamus"
nodeid <- rep(1:nrow(df_180420),2);
condition <- c(rep("Ongoing",nrow(df_180420)), rep("LightON",nrow(df_180420)));
ei_180420 <- array("Excitatory", length(regions_180420));
ei_180420[c(4,5,6,7,9,11,12,13,15,16,17,18,19,20,28,29,30,31,41,42,49,51,52,53,59,62,64,68,69,70,71)] <- "Inhibitory";
dat_180420 <- data.frame(c(df_180420$Ongoing, df_180420$LightON), condition, nodeid, rep(regions_180420,2), ei_180420);
colnames(dat_180420) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_180420$condition <- factor(dat_180420$condition, levels = c("Ongoing", "LightON"))
dat_180420$region <- factor(dat_180420$region)
dat_180420$LocalEff <- as.numeric(sub(",", ".", c(df_180420$Ongoing, df_180420$LightON), fixed = TRUE))
dat_180420$EI <- factor(ei_180420)
model_180420 <- lm(formula = "LocalEff ~ condition * region * EI + (1|nodeid)", data=dat_180420)
res_aov_180420 <- aov(model_180420)
summary(res_aov_180420)
cm=compare_means(LocalEff ~ condition, data = dat_180420, group.by = "region")

tmpdf_localeff_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_180420, mean);
#tmpdf_localeff_wei$expID = rep("180131",4,length(levels(as.factor(tmpdf_localeff_wei$region)))*2);
tmpdf_localeff_wei$expID = rep("180420",1,nrow(tmpdf_localeff_wei));
df_localeff_wei = rbind(df_localeff_wei, tmpdf_localeff_wei);

###
df_180423 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network_GLMCC_v2.pzfx", table="180423_localeff_wei")
regions_180423 <- character(nrow(df_180423));
regions_180423[c(1,2,3,4)] = "ventromedial hypothalamus"
#regions_180423[c(5)] = "submammillothalamic nucleus"
regions_180423[c(6,7,8,9,10,18,19,20,26,27,28,29,30,31,32,33,34,35,36)] = "dorsomedial hypothalamus"
regions_180423[c(11, 12,13,14,15,16,17,21,22,23,24,25,51,52,53,54,55,56)] = "posterior hypothalamus"
#regions_180423[c(37, 38,39,40,41,42,43,44,45,46,47,48,49,50)] = "premammillary nucleus"
regions_180423[c(37, 38,39,40,41,42,43,44,45,46,47,48,49,50,5)] = "mammillary complex"
nodeid <- rep(1:nrow(df_180423),2);
condition <- c(rep("Ongoing",nrow(df_180423)), rep("LightON",nrow(df_180423)));
ei_180423 <- array("Excitatory", length(regions_180423));
ei_180423[c(1,2,3,4,8,12,14,16,19,26,27,29,31,32,35,38,40,41,45,47,48,49,50,51,52,54)] <- "Inhibitory";
dat_180423 <- data.frame(c(df_180423$Ongoing, df_180423$LightON), condition, nodeid, rep(regions_180423,2), ei_180423);
colnames(dat_180423) <- c("LocalEff", "condition", "nodeid", "region", "EI")
dat_180423$condition <- factor(dat_180423$condition, levels = c("Ongoing", "LightON"))
dat_180423$region <- factor(dat_180423$region)
dat_180423$LocalEff <- as.numeric(sub(",", ".", c(df_180423$Ongoing, df_180423$LightON), fixed = TRUE))
dat_180423$EI <- factor(ei_180423);
model_180423 <- lm(formula = "LocalEff ~ condition * region * EI + (1|nodeid)", data=dat_180423)
res_aov_180423 <- aov(model_180423)
summary(res_aov_180423)
cm=compare_means(LocalEff ~ condition, data = dat_180423, group.by = "region")

tmpdf_localeff_wei = aggregate(LocalEff ~ region + condition + EI, data = dat_180423, mean);
#tmpdf_localeff_wei$expID = rep("180131",4,length(levels(as.factor(tmpdf_localeff_wei$region)))*2);
tmpdf_localeff_wei$expID = rep("180423",1,nrow(tmpdf_localeff_wei));
df_localeff_wei = rbind(df_localeff_wei, tmpdf_localeff_wei);

df_localeff_wei$expID <- factor(df_localeff_wei$expID)
model_localeff_wei <- lm(formula = "LocalEff ~ condition * region * EI", data = df_localeff_wei)
res_aov_localeff_wei <- aov(model_localeff_wei)
summary(res_aov_localeff_wei)
posthoc <- TukeyHSD(res_aov_localeff_wei)
cm=compare_means(LocalEff ~ condition, data = df_localeff_wei, group.by = "region")
ttt <- posthoc$`EI`
which(ttt[,4] < 0.05)


df_localeff_wei$phase[df_localeff_wei$expID == "171019"] = rep("Day", length(df_localeff_wei$expID == "171019"))
df_localeff_wei$phase[df_localeff_wei$expID == "180131"] = rep("Day", length(df_localeff_wei$expID == "180131"))
df_localeff_wei$phase[df_localeff_wei$expID == "180302"] = rep("Day", length(df_localeff_wei$expID == "180302"))
df_localeff_wei$phase[df_localeff_wei$expID == "180131"] = rep("Day", length(df_localeff_wei$expID == "180131"))
df_localeff_wei$phase[df_localeff_wei$expID == "180419"] = rep("Day", length(df_localeff_wei$expID == "180419"))
df_localeff_wei$phase[df_localeff_wei$expID == "180420"] = rep("Day", length(df_localeff_wei$expID == "180420"))
df_localeff_wei$phase[df_localeff_wei$expID == "180423"] = rep("Day", length(df_localeff_wei$expID == "180423"))
df_localeff_wei$phase[is.na(df_localeff_wei$phase)] = "Night";
df_localeff_wei$phase = as.factor(df_localeff_wei$phase);

model_localeff_wei <- lm(formula = "LocalEff ~ condition * region * EI * phase", data = df_localeff_wei)
res_aov_localeff_wei <- aov(model_localeff_wei)
summary(res_aov_localeff_wei)
cm=compare_means(LocalEff ~ condition + region, data = df_localeff_wei, group.by = "region")
posthoc <- TukeyHSD(res_aov_localeff_wei)
cm=compare_means(LocalEff ~ condition, data = df_localeff_wei, group.by = "region")
ttt <- posthoc$`EI`
which(ttt[,4] < 0.05)

model_localeff_wei <- lm(formula = "LocalEff ~ phase * EI", data = df_localeff_wei)
res_aov_localeff_wei <- aov(model_localeff_wei)
summary(res_aov_localeff_wei)


png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/GLMCC/region_localeff_wei_phase.png", width = 1860, height = 1440, units = "px")
ggboxplot(df_localeff_wei, x="region", y="LocalEff", color="EI", add = "jitter") + #stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("")  + ylab("Node LocalEff")
dev.off()



df_hubness = data.frame(Hubs = c(2,2,2,0,0,1,1,7,4,3,4,2,3,0,2,1,0,0,0,0,0,1,0,0,0,0,0,2,1,0,0,0,0,0,0,3,0,0,0 ), 
                        Membership = c(rep("Ongoing",1,13), rep("both",1,13), rep("LightON",1,13)), 
                        Phase = rep(c(1,2,2,2,2,2,1,2,2,1,1,1,1),3))
df_hubness$Membership <- as.factor(df_hubness$Membership)
df_hubness$Phase <- as.factor(df_hubness$Phase)
model_hubness <- lm(formula = "Hubs ~ Membership * Phase", data = df_hubness)
res_aov_hubness <- aov(model_hubness)
summary(res_aov_hubness)
posthoc <- TukeyHSD(res_aov_hubness)
