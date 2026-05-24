library(lme4)
library(rstatix)
library(pzfx)
library(ggplot2)
library(ggpubr)
library(memisc)
library(scales)
library(readxl)
library(circlize)
library(R.matlab)
library(plot.matrix)
library(cluster)
library(RColorBrewer)
library(psych)
library(network)
library(igraph)

regions_match_s2 <- read_excel("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/regions_match.xlsx", 
                               sheet = "Sheet2")
colnames(regions_match_s2) <- c("Rec_id", "Regions", "Neurons")
regions_match_s2$Rec_id <- as.factor(regions_match_s2$Rec_id)
regions_match_s2$Regions <- as.factor(regions_match_s2$Regions)
png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/regions_pie.png", width = 2560, height = 1440, units = "px")
ggplot(regions_match_s2, aes(x="", y=Neurons, fill=Regions)) + geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +  theme_void() + theme(text = element_text(size = 20))
dev.off()

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/regions_dist.png", width = 2560, height = 1440, units = "px")
ggboxplot(regions_match_s2, x="Regions", y="Neurons", fill="Regions", facet.by = "Rec_id", add = "jitter") + theme(text = element_text(size = 20)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + xlab("") + theme(legend.title = element_blank()) + theme(legend.position = "none")
dev.off()

manova <- read_excel("/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/manova.xlsx")
model = lm(cbind(manova$Degree, manova$FiringRate, manova$ClusteringCoeff, manova$LocalEfficiency) ~ manova$Condition, data=manova)
res.man <- manova(cbind(manova$Degree, manova$FiringRate, manova$ClusteringCoeff, manova$LocalEfficiency) ~ manova$Condition, data=manova)
summary(res.man)
summary.aov(res.man)

######################################################## 2-way anova Y \in {Degree, CLusCoef, LocalEff} ~ Regions + condition

################ DEGREE #######################
df_171019 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="171019_degree")
regions_171019 <- character(nrow(df_171019));
regions_171019[c(1,2,3,4,5,6,7,8,9,10, 12, 16, 19, 21, 23, 24, 26, 31, 33,43, 44, 46, 48, 49, 51, 52, 53, 54, 80, 81, 82, 83, 84, 85, 86, 87, 88)] = "interanteromedial thalamic nucleus"
regions_171019[c(11, 13, 14, 15, 17, 18, 20, 22, 25, 27, 28, 29, 30, 32, 34,35, 36, 37,38, 39, 40, 41, 42,45, 47, 50, 55, 56, 57, 58)] = "reuniens thalamic nucleus"
regions_171019[c(59, 60, 61, 62, 66, 67, 68, 69, 70, 71, 72, 73, 74)] = "mediodorsal thalamic nucleus"
regions_171019[c(63, 64, 65, 75, 76, 77, 78, 79, 95, 96, 97, 98, 99, 109, 110, 111, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126)] = "submedius thalamic nucleus"
regions_171019[c(89, 90, 91, 92, 93, 94, 100, 101, 102, 103, 104, 105, 106, 107, 108, 112, 113, 114)] = "central medial nucleus"
regions_171019[c(115,116)] = "paraventricular hypothalamic area, medial, parvicellular part"

nodeid <- rep(1:nrow(df_171019),2);
condition <- c(rep("Ongoing",nrow(df_171019)), rep("LightON",nrow(df_171019)));
dat_171019 <- data.frame(c(df_171019$Ongoing, df_171019$LightON), condition, nodeid, rep(regions_171019,2));
colnames(dat_171019) <- c("NodeDegree", "condition", "nodeid", "region")
dat_171019$condition <- factor(dat_171019$condition, levels = c("Ongoing", "LightON"))
model_171019 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_171019)
res_aov_171019 <- aov(model_171019)
summary(res_aov_171019)
cm=compare_means(NodeDegree ~ condition, data = dat_171019, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/degree_2way_171019.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171019, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171019") + ylab("Node Degree")
dev.off()

df_171207 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="171207_degree")
regions_171207 <- character(nrow(df_171207));
regions_171207[c(1,2,3, 12,13,28)] <- "reuniens thalamic nucleus"
regions_171207[c(4,5,11,14,15,16)] <- "paraventricular hypothalamic area, medial, magnocellular part"
regions_171207[c(6,7)] <- "paraventricular hypothalamic area, lateral magnocellular part"
regions_171207[c(8,9,10,17)] <- "anterior nucleus of hypothalamus"
regions_171207[c(18,23,24,25,26,27)] <- "submedius thalamic nucleus"
regions_171207[c(19,20,21,22)] <- "dorsomedial nucleus of hypothalamus"
regions_171207[c(29,30)] <- "paraventricular hypothalamic nucleus, posterior"
regions_171207[31] <- "tuber cinereum area"

nodeid <- rep(1:nrow(df_171207),2);
condition <- c(rep("Ongoing",nrow(df_171207)), rep("LightON",nrow(df_171207)));
dat_171207 <- data.frame(c(df_171207$Ongoing, df_171207$LightON), condition, nodeid, rep(regions_171207,2));
colnames(dat_171207) <- c("NodeDegree", "condition", "nodeid", "region")
dat_171207$condition <- factor(dat_171207$condition, levels = c("Ongoing", "LightON"))
model_171207 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_171207)
res_aov_171207 <- aov(model_171207)
summary(res_aov_171207)
cm=compare_means(NodeDegree ~ condition, data = dat_171207, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/degree_2way_171207.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171207, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171207") + ylab("Node Degree")
dev.off()

df_171208 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="171208_degree")
regions_171208 <- character(nrow(df_171208));
regions_171208[c(1,2,3,4,5,7,20,21,22,23,24,25,28)] <- "submedius thalamic nucleus"
regions_171208[6] <- "ventromedial hypothalamic nucleus"
regions_171208[c(8,9,11,13,16,17)] <- "ventromedial hypothalamic nucleus, dorsomedial part"
regions_171208[c(10,12,32)] <- "paraventricular hypothalamic nucleus, posterior"
regions_171208[c(14,15,18,19)] <- "dorsomedial nucleus of hypothalamus"
regions_171208[c(26, 27, 29, 30, 31,32,33, 34)] <- "anterior nucleus of hypothalamus"

nodeid <- rep(1:nrow(df_171208),2);
condition <- c(rep("Ongoing",nrow(df_171208)), rep("LightON",nrow(df_171208)));
dat_171208 <- data.frame(c(df_171208$Ongoing, df_171208$LightON), condition, nodeid, rep(regions_171208,2));
colnames(dat_171208) <- c("NodeDegree", "condition", "nodeid", "region")
dat_171208$condition <- factor(dat_171208$condition, levels = c("Ongoing", "LightON"))
model_171208 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_171208)
res_aov_171208 <- aov(model_171208)
summary(res_aov_171208)
cm=compare_means(NodeDegree ~ condition, data = dat_171208, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/degree_2way_171208.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171208, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171208") + ylab("Node Degree")
dev.off()


df_171213 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="171213_degree")
regions_171213 <- character(nrow(df_171213));
regions_171213[c(1,2,3,4,5,6,7,8,9,10,14,15,16,17,18,19,28,33)] <- "paraventricular hypothalamic area, medial, parvicellular part"
regions_171213[c(11,12,13,20,21,22,23,24,25,26,27,29,30,31,32)] <- "reuniens thalamic nucleus"

nodeid <- rep(1:nrow(df_171213),2);
condition <- c(rep("Ongoing", nrow(df_171213)), rep("LightON",nrow(df_171213)));
dat_171213 <- data.frame(c(df_171213$Ongoing, df_171213$LightON), condition, nodeid, rep(regions_171213,2));
colnames(dat_171213) <- c("NodeDegree", "condition", "nodeid", "region")
dat_171213$condition <- factor(dat_171213$condition, levels = c("Ongoing", "LightON"))
model_171213 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_171213)
res_aov_171213 <- aov(model_171213)
summary(res_aov_171213)
cm=compare_means(NodeDegree ~ condition, data = dat_171213, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/degree_2way_171213.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171213, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171213") + ylab("Node Degree")
dev.off()

df_180110 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180110_degree")
regions_180110 <- character(nrow(df_180110));
regions_180110[c(1,2,11,18,23)] <- "reuniens thalamic nucleus"
regions_180110[c(3,4,5,6,9,10)] <- "paraventricular hypothalamic area, medial, magnocellular part"
regions_180110[7] <- "paraventricular hypothalamic area, lateral magnocellular part"
regions_180110[8] <- "subparaventricular zone"
regions_180110[c(12,13,14,15,16,17)] <- "dorsomedial nucleus of hypothalamus"
regions_180110[c(19,22,24,25,26,27,28,29,30,31,32,33,34,35)] <- "paraventricular hypothalamic nucleus, posterior"
regions_180110[c(20,21)] <-"submedius thalamic nucleus"

nodeid <- rep(1:nrow(df_180110),2);
condition <- c(rep("Ongoing", nrow(df_180110)), rep("LightON",nrow(df_180110)));
dat_180110 <- data.frame(c(df_180110$Ongoing, df_180110$LightON), condition, nodeid, rep(regions_180110,2));
colnames(dat_180110) <- c("NodeDegree", "condition", "nodeid", "region")
dat_180110$condition <- factor(dat_180110$condition, levels = c("Ongoing", "LightON"))
model_180110 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_180110)
res_aov_180110 <- aov(model_180110)
summary(res_aov_180110)
cm=compare_means(NodeDegree ~ condition, data = dat_180110, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/degree_2way_180110.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180110, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180110") + ylab("Node Degree")
dev.off()


df_180111 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180111_degree")
regions_180111 <- character(nrow(df_180111));
regions_180111[c(1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,22,23,24)] <- "reuniens thalamic nucleus"
regions_180111[14] <- "paraventricular hypothalamic nucleus, posterior"
regions_180111[c(19,20,21)] <- "ventromedial hypothalamic nucleus, dorsomedial part"

nodeid <- rep(1:nrow(df_180111),2);
condition <- c(rep("Ongoing", nrow(df_180111)), rep("LightON",nrow(df_180111)));
dat_180111 <- data.frame(c(df_180111$Ongoing, df_180111$LightON), condition, nodeid, rep(regions_180111,2));
colnames(dat_180111) <- c("NodeDegree", "condition", "nodeid", "region")
dat_180111$condition <- factor(dat_180111$condition, levels = c("Ongoing", "LightON"))
model_180111 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_180111)
res_aov_180111 <- aov(model_180111)
summary(res_aov_180111)
cm=compare_means(NodeDegree ~ condition, data = dat_180111, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/degree_2way_180111.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180111, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180111") + ylab("Node Degree")
dev.off()


df_180131 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180131_degree")
regions_180131 <- character(nrow(df_180131));
regions_180131[c(1,2,4,5,6,14,15,16)] <- "subparafascicular thalamic nucleus"
regions_180131[c(3,13,18,19,20,21,22,23,24, 27,29)] <- "posterior hypothalamic area"
regions_180131[c(7,8,9,10,11,12,25)] <- "dorsomedial hypothalamic nucleus, compact part"
regions_180131[17] <- "parvicellular part of ventral posteromedial nucleus"
regions_180131[26] <- "reuniens thalamic nucleus"
regions_180131[c(28,30,31)] <- "dorsomedial hypothalamic nucleus, dorsal part"

nodeid <- rep(1:nrow(df_180131),2);
condition <- c(rep("Ongoing", nrow(df_180131)), rep("LightON",nrow(df_180131)));
dat_180131 <- data.frame(c(df_180131$Ongoing, df_180131$LightON), condition, nodeid, rep(regions_180131,2));
colnames(dat_180131) <- c("NodeDegree", "condition", "nodeid", "region")
dat_180131$condition <- factor(dat_180131$condition, levels = c("Ongoing", "LightON"))
model_180131 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_180131)
res_aov_180131 <- aov(model_180131)
summary(res_aov_180131)
cm=compare_means(NodeDegree ~ condition, data = dat_180131, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/degree_2way_180131.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180131, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180131") + ylab("Node Degree")
dev.off()


df_180221 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180221_degree")
regions_180221 <- character(nrow(df_180221));
regions_180221[c(1, 35, 36, 48,49, 50)] <- "submedius thalamic nucleus"
regions_180221[c(2,3,4,12,17,18, 21, 22, 26, 31, 32, 33, 39, 40, 41, 42, 43, 44, 45, 46, 47)] <- "reuniens thalamic nucleus"
regions_180221[c(5,6,7,8,9,10,11,13,14,15,16, 27, 28, 29, 30, 37, 38)] <- "dorsomedial nucleus of hypothalamus"
regions_180221[c(19, 20, 34, 51, 52, 53, 54)] <- "parvicellular part of ventral posteromedial nucleus"
regions_180221[c(23, 24, 25)] <- "ventromedial hypothalamic nucleus, dorsomedial part"

nodeid <- rep(1:nrow(df_180221),2);
condition <- c(rep("Ongoing", nrow(df_180221)), rep("LightON",nrow(df_180221)));
dat_180221 <- data.frame(c(df_180221$Ongoing, df_180221$LightON), condition, nodeid, rep(regions_180221,2));
colnames(dat_180221) <- c("NodeDegree", "condition", "nodeid", "region")
dat_180221$condition <- factor(dat_180221$condition, levels = c("Ongoing", "LightON"))
model_180221 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_180221)
res_aov_180221 <- aov(model_180221)
summary(res_aov_180221)
cm=compare_means(NodeDegree ~ condition, data = dat_180221, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/degree_2way_180221.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180221, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180221") + ylab("Node Degree")
dev.off()


df_180228 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180228_degree")
regions_180228 <- character(nrow(df_180228));
regions_180228[c(1,4,5,6, 24, 25, 74)] <- "ventromedial thalamic nucleus"
regions_180228[c(2, 11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23,26, 40, 42, 62, 63, 64)] <- "posterior hypothalamic area"
regions_180228[c(3,7,8,10, 19, 20, 27, 28, 29, 30, 31, 54, 56, 57, 58, 59, 60, 61, 65, 66,67,68,69, 73, 75, 76, 77,78,79)] <- "dorsomedial nucleus of hypothalamus"
regions_180228[9] <- "ventromedial hypothalamic nucleus, dorsomedial part"
regions_180228[c(32, 33, 34, 35, 36, 37, 38,39, 41)] <- "arcuate nucleus of hypothalamus"
regions_180228[c(43,44,45)] <- "subincertal nucleus"
regions_180228[c(46, 47, 48, 49, 50)] <- "ventromedial hypothalamic nucleus"
regions_180228[c(51,52,53,55,70,71, 72)] <- "zona incerta"

nodeid <- rep(1:nrow(df_180228),2);
condition <- c(rep("Ongoing", nrow(df_180228)), rep("LightON",nrow(df_180228)));
dat_180228 <- data.frame(c(df_180228$Ongoing, df_180228$LightON), condition, nodeid, rep(regions_180228,2));
colnames(dat_180228) <- c("NodeDegree", "condition", "nodeid", "region")
dat_180228$condition <- factor(dat_180228$condition, levels = c("Ongoing", "LightON"))
model_180228 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_180228)
res_aov_180228 <- aov(model_180228)
summary(res_aov_180228)
cm=compare_means(NodeDegree ~ condition, data = dat_180228, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/degree_2way_180228.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180228, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180228") + ylab("Node Degree")
dev.off()

df_180302 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180302_degree")
regions_180302 <- character(nrow(df_180302));
regions_180302[c(1,11,12)] <- "posterior hypothalamic area"
regions_180302[c(2,3,4,5,6, 13, 14, 15,30,36)] <- "dorsal tuberomammillary nucleus"
regions_180302[c(7, 8,16, 17, 23, 31, 37, 38,39,40, 47)] <- "posterior hypothalamic area"
regions_180302[c(9, 10, 41,42)] <- "premammillary nucleus dorsal"
regions_180302[c(18, 19, 20, 21)] <- "supramammillary nucleus"
regions_180302[c(22, 24, 25, 29)] <- "medial magnocellular hypothalamic area"
regions_180302[c(26, 27, 28, 43,44,45,46)] <- "arcuate nucleus of hypothalamus"
regions_180302[c(32,33,34,35)] <- "dorsomedial hypothalamic nucleus, compact part"

nodeid <- rep(1:nrow(df_180302),2);
condition <- c(rep("Ongoing", nrow(df_180302)), rep("LightON",nrow(df_180302)));
dat_180302 <- data.frame(c(df_180302$Ongoing, df_180302$LightON), condition, nodeid, rep(regions_180302,2));
colnames(dat_180302) <- c("NodeDegree", "condition", "nodeid", "region")
dat_180302$condition <- factor(dat_180302$condition, levels = c("Ongoing", "LightON"))
model_180302 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_180302)
res_aov_180302 <- aov(model_180302)
summary(res_aov_180302)
cm=compare_means(NodeDegree ~ condition, data = dat_180302, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/degree_2way_180302.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180302, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180302") + ylab("Node Degree")
dev.off()

df_180419 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180419_degree")
regions_180419 <- character(nrow(df_180419));
regions_180419[c(1,2,11,12)] <- "ventromedial hypothalamic nucleus, dorsomedial part"
regions_180419[c(3, 40, 41, 42, 43, 44,45,46,47,48,49,50)] <- "dorsomedial nucleus of hypothalamus"
regions_180419[c(4,5,6,7,8)] <- "submedius thalamic nucleus"
regions_180419[c(9, 10)] <- "ventromedial thalamic nucleus"
regions_180419[c(13,14,15,16,17,18,23)] <- "dorsomedial hypothalamic nucleus, ventral part"
regions_180419[c(19,20,21,22,24,25,26,27,28,29,30,31,32,33,34,35)] <- "dorsomedial hypothalamic nucleus, dorsal part"
regions_180419[c(36,37,38,39,51,52,53,54,55,56,57,58,59,60,61)] <- "posterior hypothalamic area"

nodeid <- rep(1:nrow(df_180419),2);
condition <- c(rep("Ongoing", nrow(df_180419)), rep("LightON",nrow(df_180419)));
dat_180419 <- data.frame(c(df_180419$Ongoing, df_180419$LightON), condition, nodeid, rep(regions_180419,2));
colnames(dat_180419) <- c("NodeDegree", "condition", "nodeid", "region")
dat_180419$condition <- factor(dat_180419$condition, levels = c("Ongoing", "LightON"))
model_180419 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_180419)
res_aov_180419 <- aov(model_180419)
summary(res_aov_180419)
cm=compare_means(NodeDegree ~ condition, data = dat_180419, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/degree_2way_180419.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180419, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180419") + ylab("Node Degree")
dev.off()

df_180420 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180420_degree")
regions_180420 <- character(nrow(df_180420));
regions_180420[1] <- "dorsomedial hypothalamic nucleus, compact part"
regions_180420[c(2, 3, 4, 5, 6,7,8,9,10,11,12,13,14,15,16,17,18,64,65,66,67,68,69,70,71)] <- "posterior hypothalamic area"
regions_180420[c(19, 20, 21, 22, 23, 24, 25,26,27,28,29,30,31)] <- "premammillary nucleus dorsal"
regions_180420[c(32, 33, 34, 35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50)] <- "medial magnocellular hypothalamic area"
regions_180420[c(51,52,53,54,55,56,57,58,59,60,61,62,63)] <- "supramammillary nucleus, medial part"

nodeid <- rep(1:nrow(df_180420),2);
condition <- c(rep("Ongoing", nrow(df_180420)), rep("LightON",nrow(df_180420)));
dat_180420 <- data.frame(c(df_180420$Ongoing, df_180420$LightON), condition, nodeid, rep(regions_180420,2));
colnames(dat_180420) <- c("NodeDegree", "condition", "nodeid", "region")
dat_180420$condition <- factor(dat_180420$condition, levels = c("Ongoing", "LightON"))
model_180420 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_180420)
res_aov_180420 <- aov(model_180420)
summary(res_aov_180420)
cm=compare_means(NodeDegree ~ condition, data = dat_180420, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/degree_2way_180420.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180420, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180420") + ylab("Node Degree")
dev.off()

df_180423 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180423_degree")
regions_180423 <- character(nrow(df_180423));
regions_180423[c(1,2,3,4)] <- "ventromedial hypothalamic nucleus"
regions_180423[5] <- "lateral hypothalamic area"
regions_180423[c(6,7,8,9,10,18,19,20,26,27,28,29,30,31,32,33,34,35,36)] <- "dorsomedial nucleus of hypothalamus"
regions_180423[c(11, 12,13,14,15,16,17,21,22,23,24,25,51,52,53,54,55,56)] <- "posterior hypothalamic area"
regions_180423[37] <- "premammillary nucleus, ventral part"
regions_180423[c(38,39,40,41,42,43,44,45,46,47,48,49,50)] <- "premammillary nucleus, dorsal part"

nodeid <- rep(1:nrow(df_180423),2);
condition <- c(rep("Ongoing", nrow(df_180423)), rep("LightON",nrow(df_180423)));
dat_180423 <- data.frame(c(df_180423$Ongoing, df_180423$LightON), condition, nodeid, rep(regions_180423,2));
colnames(dat_180423) <- c("NodeDegree", "condition", "nodeid", "region")
dat_180423$condition <- factor(dat_180423$condition, levels = c("Ongoing", "LightON"))
model_180423 <- lm(formula = "NodeDegree ~ condition * region + (1|nodeid)", data=dat_180423)
res_aov_180423 <- aov(model_180423)
summary(res_aov_180423)
cm=compare_means(NodeDegree ~ condition, data = dat_180423, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/degree_2way_180423.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180423, x="region", y="NodeDegree", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180423") + ylab("Node Degree")
dev.off()



############### CLUSTERING COEFFICIENT ####################
df_171019 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="171019_cluscoeff")
regions_171019 <- character(nrow(df_171019));
regions_171019[c(1,2,3,4,5,6,7,8,9,10, 12, 16, 19, 21, 23, 24, 26, 31, 33,43, 44, 46, 48, 49, 51, 52, 53, 54, 80, 81, 82, 83, 84, 85, 86, 87, 88)] = "interanteromedial thalamic nucleus"
regions_171019[c(11, 13, 14, 15, 17, 18, 20, 22, 25, 27, 28, 29, 30, 32, 34,35, 36, 37,38, 39, 40, 41, 42,45, 47, 50, 55, 56, 57, 58)] = "reuniens thalamic nucleus"
regions_171019[c(59, 60, 61, 62, 66, 67, 68, 69, 70, 71, 72, 73, 74)] = "mediodorsal thalamic nucleus"
regions_171019[c(63, 64, 65, 75, 76, 77, 78, 79, 95, 96, 97, 98, 99, 109, 110, 111, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126)] = "submedius thalamic nucleus"
regions_171019[c(89, 90, 91, 92, 93, 94, 100, 101, 102, 103, 104, 105, 106, 107, 108, 112, 113, 114)] = "central medial nucleus"
regions_171019[c(115,116)] = "paraventricular hypothalamic area, medial, parvicellular part"

nodeid <- rep(1:nrow(df_171019),2);
condition <- c(rep("Ongoing",nrow(df_171019)), rep("LightON",nrow(df_171019)));
dat_171019 <- data.frame(c(df_171019$Ongoing, df_171019$LightON), condition, nodeid, rep(regions_171019,2));
colnames(dat_171019) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_171019$condition <- factor(dat_171019$condition, levels = c("Ongoing", "LightON"))
model_171019 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_171019)
res_aov_171019 <- aov(model_171019)
summary(res_aov_171019)
cm=compare_means(ClusCoeff ~ condition, data = dat_171019, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/cluscoeff_2way_171019.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171019, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171019") + ylab("ClusCoeff")
dev.off()

df_171207 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="171207_cluscoeff")
regions_171207 <- character(nrow(df_171207));
regions_171207[c(1,2,3, 12,13,28)] <- "reuniens thalamic nucleus"
regions_171207[c(4,5,11,14,15,16)] <- "paraventricular hypothalamic area, medial, magnocellular part"
regions_171207[c(6,7)] <- "paraventricular hypothalamic area, lateral magnocellular part"
regions_171207[c(8,9,10,17)] <- "anterior nucleus of hypothalamus"
regions_171207[c(18,23,24,25,26,27)] <- "submedius thalamic nucleus"
regions_171207[c(19,20,21,22)] <- "dorsomedial nucleus of hypothalamus"
regions_171207[c(29,30)] <- "paraventricular hypothalamic nucleus, posterior"
regions_171207[31] <- "tuber cinereum area"

nodeid <- rep(1:nrow(df_171207),2);
condition <- c(rep("Ongoing",nrow(df_171207)), rep("LightON",nrow(df_171207)));
dat_171207 <- data.frame(c(df_171207$Ongoing, df_171207$LightON), condition, nodeid, rep(regions_171207,2));
colnames(dat_171207) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_171207$condition <- factor(dat_171207$condition, levels = c("Ongoing", "LightON"))
model_171207 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_171207)
res_aov_171207 <- aov(model_171207)
summary(res_aov_171207)
cm=compare_means(ClusCoeff ~ condition, data = dat_171207, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/cluscoeff_2way_171207.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171207, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171207") + ylab("ClusCoeff")
dev.off()

df_171208 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="171208_cluscoeff")
regions_171208 <- character(nrow(df_171208));
regions_171208[c(1,2,3,4,5,7,20,21,22,23,24,25,28)] <- "submedius thalamic nucleus"
regions_171208[6] <- "ventromedial hypothalamic nucleus"
regions_171208[c(8,9,11,13,16,17)] <- "ventromedial hypothalamic nucleus, dorsomedial part"
regions_171208[c(10,12,32)] <- "paraventricular hypothalamic nucleus, posterior"
regions_171208[c(14,15,18,19)] <- "dorsomedial nucleus of hypothalamus"
regions_171208[c(26, 27, 29, 30, 31,32,33, 34)] <- "anterior nucleus of hypothalamus"

nodeid <- rep(1:nrow(df_171208),2);
condition <- c(rep("Ongoing",nrow(df_171208)), rep("LightON",nrow(df_171208)));
dat_171208 <- data.frame(c(df_171208$Ongoing, df_171208$LightON), condition, nodeid, rep(regions_171208,2));
colnames(dat_171208) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_171208$condition <- factor(dat_171208$condition, levels = c("Ongoing", "LightON"))
model_171208 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_171208)
res_aov_171208 <- aov(model_171208)
summary(res_aov_171208)
cm=compare_means(ClusCoeff ~ condition, data = dat_171208, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/cluscoeff_2way_171208.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171208, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171208") + ylab("ClusCoeff")
dev.off()


df_171213 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="171213_cluscoeff")
regions_171213 <- character(nrow(df_171213));
regions_171213[c(1,2,3,4,5,6,7,8,9,10,14,15,16,17,18,19,28,33)] <- "paraventricular hypothalamic area, medial, parvicellular part"
regions_171213[c(11,12,13,20,21,22,23,24,25,26,27,29,30,31,32)] <- "reuniens thalamic nucleus"

nodeid <- rep(1:nrow(df_171213),2);
condition <- c(rep("Ongoing", nrow(df_171213)), rep("LightON",nrow(df_171213)));
dat_171213 <- data.frame(c(df_171213$Ongoing, df_171213$LightON), condition, nodeid, rep(regions_171213,2));
colnames(dat_171213) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_171213$condition <- factor(dat_171213$condition, levels = c("Ongoing", "LightON"))
model_171213 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_171213)
res_aov_171213 <- aov(model_171213)
summary(res_aov_171213)
cm=compare_means(ClusCoeff ~ condition, data = dat_171213, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/cluscoeff_2way_171213.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171213, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171213") + ylab("ClusCoeff")
dev.off()

df_180110 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180110_cluscoeff")
regions_180110 <- character(nrow(df_180110));
regions_180110[c(1,2,11,18,23)] <- "reuniens thalamic nucleus"
regions_180110[c(3,4,5,6,9,10)] <- "paraventricular hypothalamic area, medial, magnocellular part"
regions_180110[7] <- "paraventricular hypothalamic area, lateral magnocellular part"
regions_180110[8] <- "subparaventricular zone"
regions_180110[c(12,13,14,15,16,17)] <- "dorsomedial nucleus of hypothalamus"
regions_180110[c(19,22,24,25,26,27,28,29,30,31,32,33,34,35)] <- "paraventricular hypothalamic nucleus, posterior"
regions_180110[c(20,21)] <-"submedius thalamic nucleus"

nodeid <- rep(1:nrow(df_180110),2);
condition <- c(rep("Ongoing", nrow(df_180110)), rep("LightON",nrow(df_180110)));
dat_180110 <- data.frame(c(df_180110$Ongoing, df_180110$LightON), condition, nodeid, rep(regions_180110,2));
colnames(dat_180110) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_180110$condition <- factor(dat_180110$condition, levels = c("Ongoing", "LightON"))
model_180110 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_180110)
res_aov_180110 <- aov(model_180110)
summary(res_aov_180110)
cm=compare_means(ClusCoeff ~ condition, data = dat_180110, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/cluscoeff_2way_180110.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180110, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180110") + ylab("ClusCoeff")
dev.off()


df_180111 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180111_cluscoeff")
regions_180111 <- character(nrow(df_180111));
regions_180111[c(1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,22,23,24)] <- "reuniens thalamic nucleus"
regions_180111[14] <- "paraventricular hypothalamic nucleus, posterior"
regions_180111[c(19,20,21)] <- "ventromedial hypothalamic nucleus, dorsomedial part"

nodeid <- rep(1:nrow(df_180111),2);
condition <- c(rep("Ongoing", nrow(df_180111)), rep("LightON",nrow(df_180111)));
dat_180111 <- data.frame(c(df_180111$Ongoing, df_180111$LightON), condition, nodeid, rep(regions_180111,2));
colnames(dat_180111) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_180111$condition <- factor(dat_180111$condition, levels = c("Ongoing", "LightON"))
model_180111 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_180111)
res_aov_180111 <- aov(model_180111)
summary(res_aov_180111)
cm=compare_means(ClusCoeff ~ condition, data = dat_180111, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/cluscoeff_2way_180111.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180111, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180111") + ylab("ClusCoeff")
dev.off()


df_180131 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180131_cluscoeff")
regions_180131 <- character(nrow(df_180131));
regions_180131[c(1,2,4,5,6,14,15,16)] <- "subparafascicular thalamic nucleus"
regions_180131[c(3,13,18,19,20,21,22,23,24, 27,29)] <- "posterior hypothalamic area"
regions_180131[c(7,8,9,10,11,12,25)] <- "dorsomedial hypothalamic nucleus, compact part"
regions_180131[17] <- "parvicellular part of ventral posteromedial nucleus"
regions_180131[26] <- "reuniens thalamic nucleus"
regions_180131[c(28,30,31)] <- "dorsomedial hypothalamic nucleus, dorsal part"

nodeid <- rep(1:nrow(df_180131),2);
condition <- c(rep("Ongoing", nrow(df_180131)), rep("LightON",nrow(df_180131)));
dat_180131 <- data.frame(c(df_180131$Ongoing, df_180131$LightON), condition, nodeid, rep(regions_180131,2));
colnames(dat_180131) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_180131$condition <- factor(dat_180131$condition, levels = c("Ongoing", "LightON"))
model_180131 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_180131)
res_aov_180131 <- aov(model_180131)
summary(res_aov_180131)
cm=compare_means(ClusCoeff ~ condition, data = dat_180131, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/cluscoeff_2way_180131.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180131, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180131") + ylab("ClusCoeff")
dev.off()


df_180221 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180221_cluscoeff")
regions_180221 <- character(nrow(df_180221));
regions_180221[c(1, 35, 36, 48,49, 50)] <- "submedius thalamic nucleus"
regions_180221[c(2,3,4,12,17,18, 21, 22, 26, 31, 32, 33, 39, 40, 41, 42, 43, 44, 45, 46, 47)] <- "reuniens thalamic nucleus"
regions_180221[c(5,6,7,8,9,10,11,13,14,15,16, 27, 28, 29, 30, 37, 38)] <- "dorsomedial nucleus of hypothalamus"
regions_180221[c(19, 20, 34, 51, 52, 53, 54)] <- "parvicellular part of ventral posteromedial nucleus"
regions_180221[c(23, 24, 25)] <- "ventromedial hypothalamic nucleus, dorsomedial part"

nodeid <- rep(1:nrow(df_180221),2);
condition <- c(rep("Ongoing", nrow(df_180221)), rep("LightON",nrow(df_180221)));
dat_180221 <- data.frame(c(df_180221$Ongoing, df_180221$LightON), condition, nodeid, rep(regions_180221,2));
colnames(dat_180221) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_180221$condition <- factor(dat_180221$condition, levels = c("Ongoing", "LightON"))
model_180221 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_180221)
res_aov_180221 <- aov(model_180221)
summary(res_aov_180221)
cm=compare_means(ClusCoeff ~ condition, data = dat_180221, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/cluscoeff_2way_180221.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180221, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180221") + ylab("ClusCoeff")
dev.off()


df_180228 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180228_cluscoeff")
regions_180228 <- character(nrow(df_180228));
regions_180228[c(1,4,5,6, 24, 25, 74)] <- "ventromedial thalamic nucleus"
regions_180228[c(2, 11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23,26, 40, 42, 62, 63, 64)] <- "posterior hypothalamic area"
regions_180228[c(3,7,8,10, 19, 20, 27, 28, 29, 30, 31, 54, 56, 57, 58, 59, 60, 61, 65, 66,67,68,69, 73, 75, 76, 77,78,79)] <- "dorsomedial nucleus of hypothalamus"
regions_180228[9] <- "ventromedial hypothalamic nucleus, dorsomedial part"
regions_180228[c(32, 33, 34, 35, 36, 37, 38,39, 41)] <- "arcuate nucleus of hypothalamus"
regions_180228[c(43,44,45)] <- "subincertal nucleus"
regions_180228[c(46, 47, 48, 49, 50)] <- "ventromedial hypothalamic nucleus"
regions_180228[c(51,52,53,55,70,71, 72)] <- "zona incerta"

nodeid <- rep(1:nrow(df_180228),2);
condition <- c(rep("Ongoing", nrow(df_180228)), rep("LightON",nrow(df_180228)));
dat_180228 <- data.frame(c(df_180228$Ongoing, df_180228$LightON), condition, nodeid, rep(regions_180228,2));
colnames(dat_180228) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_180228$condition <- factor(dat_180228$condition, levels = c("Ongoing", "LightON"))
model_180228 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_180228)
res_aov_180228 <- aov(model_180228)
summary(res_aov_180228)
cm=compare_means(ClusCoeff ~ condition, data = dat_180228, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/cluscoeff_2way_180228.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180228, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180228") + ylab("ClusCoeff")
dev.off()

df_180302 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180302_cluscoeff")
regions_180302 <- character(nrow(df_180302));
regions_180302[c(1,11,12)] <- "posterior hypothalamic area"
regions_180302[c(2,3,4,5,6, 13, 14, 15,30,36)] <- "dorsal tuberomammillary nucleus"
regions_180302[c(7, 8,16, 17, 23, 31, 37, 38,39,40, 47)] <- "posterior hypothalamic area"
regions_180302[c(9, 10, 41,42)] <- "premammillary nucleus dorsal"
regions_180302[c(18, 19, 20, 21)] <- "supramammillary nucleus"
regions_180302[c(22, 24, 25, 29)] <- "medial magnocellular hypothalamic area"
regions_180302[c(26, 27, 28, 43,44,45,46)] <- "arcuate nucleus of hypothalamus"
regions_180302[c(32,33,34,35)] <- "dorsomedial hypothalamic nucleus, compact part"

nodeid <- rep(1:nrow(df_180302),2);
condition <- c(rep("Ongoing", nrow(df_180302)), rep("LightON",nrow(df_180302)));
dat_180302 <- data.frame(c(df_180302$Ongoing, df_180302$LightON), condition, nodeid, rep(regions_180302,2));
colnames(dat_180302) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_180302$condition <- factor(dat_180302$condition, levels = c("Ongoing", "LightON"))
model_180302 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_180302)
res_aov_180302 <- aov(model_180302)
summary(res_aov_180302)
cm=compare_means(ClusCoeff ~ condition, data = dat_180302, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/cluscoeff_2way_180302.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180302, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180302") + ylab("ClusCoeff")
dev.off()

df_180419 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180419_cluscoeff")
regions_180419 <- character(nrow(df_180419));
regions_180419[c(1,2,11,12)] <- "ventromedial hypothalamic nucleus, dorsomedial part"
regions_180419[c(3, 40, 41, 42, 43, 44,45,46,47,48,49,50)] <- "dorsomedial nucleus of hypothalamus"
regions_180419[c(4,5,6,7,8)] <- "submedius thalamic nucleus"
regions_180419[c(9, 10)] <- "ventromedial thalamic nucleus"
regions_180419[c(13,14,15,16,17,18,23)] <- "dorsomedial hypothalamic nucleus, ventral part"
regions_180419[c(19,20,21,22,24,25,26,27,28,29,30,31,32,33,34,35)] <- "dorsomedial hypothalamic nucleus, dorsal part"
regions_180419[c(36,37,38,39,51,52,53,54,55,56,57,58,59,60,61)] <- "posterior hypothalamic area"

nodeid <- rep(1:nrow(df_180419),2);
condition <- c(rep("Ongoing", nrow(df_180419)), rep("LightON",nrow(df_180419)));
dat_180419 <- data.frame(c(df_180419$Ongoing, df_180419$LightON), condition, nodeid, rep(regions_180419,2));
colnames(dat_180419) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_180419$condition <- factor(dat_180419$condition, levels = c("Ongoing", "LightON"))
model_180419 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_180419)
res_aov_180419 <- aov(model_180419)
summary(res_aov_180419)
cm=compare_means(ClusCoeff ~ condition, data = dat_180419, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/cluscoeff_2way_180419.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180419, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180419") + ylab("ClusCoeff")
dev.off()

df_180420 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180420_cluscoeff")
regions_180420 <- character(nrow(df_180420));
regions_180420[1] <- "dorsomedial hypothalamic nucleus, compact part"
regions_180420[c(2, 3, 4, 5, 6,7,8,9,10,11,12,13,14,15,16,17,18,64,65,66,67,68,69,70,71)] <- "posterior hypothalamic area"
regions_180420[c(19, 20, 21, 22, 23, 24, 25,26,27,28,29,30,31)] <- "premammillary nucleus dorsal"
regions_180420[c(32, 33, 34, 35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50)] <- "medial magnocellular hypothalamic area"
regions_180420[c(51,52,53,54,55,56,57,58,59,60,61,62,63)] <- "supramammillary nucleus, medial part"

nodeid <- rep(1:nrow(df_180420),2);
condition <- c(rep("Ongoing", nrow(df_180420)), rep("LightON",nrow(df_180420)));
dat_180420 <- data.frame(c(df_180420$Ongoing, df_180420$LightON), condition, nodeid, rep(regions_180420,2));
colnames(dat_180420) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_180420$condition <- factor(dat_180420$condition, levels = c("Ongoing", "LightON"))
model_180420 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_180420)
res_aov_180420 <- aov(model_180420)
summary(res_aov_180420)
cm=compare_means(ClusCoeff ~ condition, data = dat_180420, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/cluscoeff_2way_180420.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180420, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180420") + ylab("ClusCoeff")
dev.off()

df_180423 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180423_cluscoeff")
regions_180423 <- character(nrow(df_180423));
regions_180423[c(1,2,3,4)] <- "ventromedial hypothalamic nucleus"
regions_180423[5] <- "lateral hypothalamic area"
regions_180423[c(6,7,8,9,10,18,19,20,26,27,28,29,30,31,32,33,34,35,36)] <- "dorsomedial nucleus of hypothalamus"
regions_180423[c(11, 12,13,14,15,16,17,21,22,23,24,25,51,52,53,54,55,56)] <- "posterior hypothalamic area"
regions_180423[37] <- "premammillary nucleus, ventral part"
regions_180423[c(38,39,40,41,42,43,44,45,46,47,48,49,50)] <- "premammillary nucleus, dorsal part"

nodeid <- rep(1:nrow(df_180423),2);
condition <- c(rep("Ongoing", nrow(df_180423)), rep("LightON",nrow(df_180423)));
dat_180423 <- data.frame(c(df_180423$Ongoing, df_180423$LightON), condition, nodeid, rep(regions_180423,2));
colnames(dat_180423) <- c("ClusCoeff", "condition", "nodeid", "region")
dat_180423$condition <- factor(dat_180423$condition, levels = c("Ongoing", "LightON"))
model_180423 <- lm(formula = "ClusCoeff ~ condition * region + (1|nodeid)", data=dat_180423)
res_aov_180423 <- aov(model_180423)
summary(res_aov_180423)
cm=compare_means(ClusCoeff ~ condition, data = dat_180423, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/cluscoeff_2way_180423.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180423, x="region", y="ClusCoeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180423") + ylab("ClusCoeff")
dev.off()


######################### LCOAL EFFICIENCY ###################
df_171019 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="171019_localeff")
regions_171019 <- character(nrow(df_171019));
regions_171019[c(1,2,3,4,5,6,7,8,9,10, 12, 16, 19, 21, 23, 24, 26, 31, 33,43, 44, 46, 48, 49, 51, 52, 53, 54, 80, 81, 82, 83, 84, 85, 86, 87, 88)] = "interanteromedial thalamic nucleus"
regions_171019[c(11, 13, 14, 15, 17, 18, 20, 22, 25, 27, 28, 29, 30, 32, 34,35, 36, 37,38, 39, 40, 41, 42,45, 47, 50, 55, 56, 57, 58)] = "reuniens thalamic nucleus"
regions_171019[c(59, 60, 61, 62, 66, 67, 68, 69, 70, 71, 72, 73, 74)] = "mediodorsal thalamic nucleus"
regions_171019[c(63, 64, 65, 75, 76, 77, 78, 79, 95, 96, 97, 98, 99, 109, 110, 111, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126)] = "submedius thalamic nucleus"
regions_171019[c(89, 90, 91, 92, 93, 94, 100, 101, 102, 103, 104, 105, 106, 107, 108, 112, 113, 114)] = "central medial nucleus"
regions_171019[c(115,116)] = "paraventricular hypothalamic area, medial, parvicellular part"

nodeid <- rep(1:nrow(df_171019),2);
condition <- c(rep("Ongoing",nrow(df_171019)), rep("LightON",nrow(df_171019)));
dat_171019 <- data.frame(c(df_171019$Ongoing, df_171019$LightON), condition, nodeid, rep(regions_171019,2));
colnames(dat_171019) <- c("localeff", "condition", "nodeid", "region")
dat_171019$condition <- factor(dat_171019$condition, levels = c("Ongoing", "LightON"))
model_171019 <- lm(formula = "localeff ~ condition * region + (1|nodeid)", data=dat_171019)
res_aov_171019 <- aov(model_171019)
summary(res_aov_171019)
cm=compare_means(localeff ~ condition, data = dat_171019, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/localeff_2way_171019.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171019, x="region", y="localeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171019") + ylab("localeff")
dev.off()

df_171207 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="171207_localeff")
regions_171207 <- character(nrow(df_171207));
regions_171207[c(1,2,3, 12,13,28)] <- "reuniens thalamic nucleus"
regions_171207[c(4,5,11,14,15,16)] <- "paraventricular hypothalamic area, medial, magnocellular part"
regions_171207[c(6,7)] <- "paraventricular hypothalamic area, lateral magnocellular part"
regions_171207[c(8,9,10,17)] <- "anterior nucleus of hypothalamus"
regions_171207[c(18,23,24,25,26,27)] <- "submedius thalamic nucleus"
regions_171207[c(19,20,21,22)] <- "dorsomedial nucleus of hypothalamus"
regions_171207[c(29,30)] <- "paraventricular hypothalamic nucleus, posterior"
regions_171207[31] <- "tuber cinereum area"

nodeid <- rep(1:nrow(df_171207),2);
condition <- c(rep("Ongoing",nrow(df_171207)), rep("LightON",nrow(df_171207)));
dat_171207 <- data.frame(c(df_171207$Ongoing, df_171207$LightON), condition, nodeid, rep(regions_171207,2));
colnames(dat_171207) <- c("localeff", "condition", "nodeid", "region")
dat_171207$condition <- factor(dat_171207$condition, levels = c("Ongoing", "LightON"))
model_171207 <- lm(formula = "localeff ~ condition * region + (1|nodeid)", data=dat_171207)
res_aov_171207 <- aov(model_171207)
summary(res_aov_171207)
cm=compare_means(localeff ~ condition, data = dat_171207, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/localeff_2way_171207.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171207, x="region", y="localeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171207") + ylab("localeff")
dev.off()

df_171208 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="171208_localeff")
regions_171208 <- character(nrow(df_171208));
regions_171208[c(1,2,3,4,5,7,20,21,22,23,24,25,28)] <- "submedius thalamic nucleus"
regions_171208[6] <- "ventromedial hypothalamic nucleus"
regions_171208[c(8,9,11,13,16,17)] <- "ventromedial hypothalamic nucleus, dorsomedial part"
regions_171208[c(10,12,32)] <- "paraventricular hypothalamic nucleus, posterior"
regions_171208[c(14,15,18,19)] <- "dorsomedial nucleus of hypothalamus"
regions_171208[c(26, 27, 29, 30, 31,32,33, 34)] <- "anterior nucleus of hypothalamus"

nodeid <- rep(1:nrow(df_171208),2);
condition <- c(rep("Ongoing",nrow(df_171208)), rep("LightON",nrow(df_171208)));
dat_171208 <- data.frame(c(df_171208$Ongoing, df_171208$LightON), condition, nodeid, rep(regions_171208,2));
colnames(dat_171208) <- c("localeff", "condition", "nodeid", "region")
dat_171208$condition <- factor(dat_171208$condition, levels = c("Ongoing", "LightON"))
model_171208 <- lm(formula = "localeff ~ condition * region + (1|nodeid)", data=dat_171208)
res_aov_171208 <- aov(model_171208)
summary(res_aov_171208)
cm=compare_means(localeff ~ condition, data = dat_171208, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/localeff_2way_171208.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171208, x="region", y="localeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171208") + ylab("localeff")
dev.off()


df_171213 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="171213_localeff")
regions_171213 <- character(nrow(df_171213));
regions_171213[c(1,2,3,4,5,6,7,8,9,10,14,15,16,17,18,19,28,33)] <- "paraventricular hypothalamic area, medial, parvicellular part"
regions_171213[c(11,12,13,20,21,22,23,24,25,26,27,29,30,31,32)] <- "reuniens thalamic nucleus"

nodeid <- rep(1:nrow(df_171213),2);
condition <- c(rep("Ongoing", nrow(df_171213)), rep("LightON",nrow(df_171213)));
dat_171213 <- data.frame(c(df_171213$Ongoing, df_171213$LightON), condition, nodeid, rep(regions_171213,2));
colnames(dat_171213) <- c("localeff", "condition", "nodeid", "region")
dat_171213$condition <- factor(dat_171213$condition, levels = c("Ongoing", "LightON"))
model_171213 <- lm(formula = "localeff ~ condition * region + (1|nodeid)", data=dat_171213)
res_aov_171213 <- aov(model_171213)
summary(res_aov_171213)
cm=compare_means(localeff ~ condition, data = dat_171213, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/localeff_2way_171213.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_171213, x="region", y="localeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("171213") + ylab("localeff")
dev.off()

df_180110 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180110_localeff")
regions_180110 <- character(nrow(df_180110));
regions_180110[c(1,2,11,18,23)] <- "reuniens thalamic nucleus"
regions_180110[c(3,4,5,6,9,10)] <- "paraventricular hypothalamic area, medial, magnocellular part"
regions_180110[7] <- "paraventricular hypothalamic area, lateral magnocellular part"
regions_180110[8] <- "subparaventricular zone"
regions_180110[c(12,13,14,15,16,17)] <- "dorsomedial nucleus of hypothalamus"
regions_180110[c(19,22,24,25,26,27,28,29,30,31,32,33,34,35)] <- "paraventricular hypothalamic nucleus, posterior"
regions_180110[c(20,21)] <-"submedius thalamic nucleus"

nodeid <- rep(1:nrow(df_180110),2);
condition <- c(rep("Ongoing", nrow(df_180110)), rep("LightON",nrow(df_180110)));
dat_180110 <- data.frame(c(df_180110$Ongoing, df_180110$LightON), condition, nodeid, rep(regions_180110,2));
colnames(dat_180110) <- c("localeff", "condition", "nodeid", "region")
dat_180110$condition <- factor(dat_180110$condition, levels = c("Ongoing", "LightON"))
model_180110 <- lm(formula = "localeff ~ condition * region + (1|nodeid)", data=dat_180110)
res_aov_180110 <- aov(model_180110)
summary(res_aov_180110)
cm=compare_means(localeff ~ condition, data = dat_180110, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/localeff_2way_180110.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180110, x="region", y="localeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180110") + ylab("localeff")
dev.off()


df_180111 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180111_localeff")
regions_180111 <- character(nrow(df_180111));
regions_180111[c(1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,22,23,24)] <- "reuniens thalamic nucleus"
regions_180111[14] <- "paraventricular hypothalamic nucleus, posterior"
regions_180111[c(19,20,21)] <- "ventromedial hypothalamic nucleus, dorsomedial part"

nodeid <- rep(1:nrow(df_180111),2);
condition <- c(rep("Ongoing", nrow(df_180111)), rep("LightON",nrow(df_180111)));
dat_180111 <- data.frame(c(df_180111$Ongoing, df_180111$LightON), condition, nodeid, rep(regions_180111,2));
colnames(dat_180111) <- c("localeff", "condition", "nodeid", "region")
dat_180111$condition <- factor(dat_180111$condition, levels = c("Ongoing", "LightON"))
model_180111 <- lm(formula = "localeff ~ condition * region + (1|nodeid)", data=dat_180111)
res_aov_180111 <- aov(model_180111)
summary(res_aov_180111)
cm=compare_means(localeff ~ condition, data = dat_180111, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/localeff_2way_180111.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180111, x="region", y="localeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180111") + ylab("localeff")
dev.off()


df_180131 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180131_localeff")
regions_180131 <- character(nrow(df_180131));
regions_180131[c(1,2,4,5,6,14,15,16)] <- "subparafascicular thalamic nucleus"
regions_180131[c(3,13,18,19,20,21,22,23,24, 27,29)] <- "posterior hypothalamic area"
regions_180131[c(7,8,9,10,11,12,25)] <- "dorsomedial hypothalamic nucleus, compact part"
regions_180131[17] <- "parvicellular part of ventral posteromedial nucleus"
regions_180131[26] <- "reuniens thalamic nucleus"
regions_180131[c(28,30,31)] <- "dorsomedial hypothalamic nucleus, dorsal part"

nodeid <- rep(1:nrow(df_180131),2);
condition <- c(rep("Ongoing", nrow(df_180131)), rep("LightON",nrow(df_180131)));
dat_180131 <- data.frame(c(df_180131$Ongoing, df_180131$LightON), condition, nodeid, rep(regions_180131,2));
colnames(dat_180131) <- c("localeff", "condition", "nodeid", "region")
dat_180131$condition <- factor(dat_180131$condition, levels = c("Ongoing", "LightON"))
model_180131 <- lm(formula = "localeff ~ condition * region + (1|nodeid)", data=dat_180131)
res_aov_180131 <- aov(model_180131)
summary(res_aov_180131)
cm=compare_means(localeff ~ condition, data = dat_180131, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/localeff_2way_180131.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180131, x="region", y="localeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180131") + ylab("localeff")
dev.off()


df_180221 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180221_localeff")
regions_180221 <- character(nrow(df_180221));
regions_180221[c(1, 35, 36, 48,49, 50)] <- "submedius thalamic nucleus"
regions_180221[c(2,3,4,12,17,18, 21, 22, 26, 31, 32, 33, 39, 40, 41, 42, 43, 44, 45, 46, 47)] <- "reuniens thalamic nucleus"
regions_180221[c(5,6,7,8,9,10,11,13,14,15,16, 27, 28, 29, 30, 37, 38)] <- "dorsomedial nucleus of hypothalamus"
regions_180221[c(19, 20, 34, 51, 52, 53, 54)] <- "parvicellular part of ventral posteromedial nucleus"
regions_180221[c(23, 24, 25)] <- "ventromedial hypothalamic nucleus, dorsomedial part"

nodeid <- rep(1:nrow(df_180221),2);
condition <- c(rep("Ongoing", nrow(df_180221)), rep("LightON",nrow(df_180221)));
dat_180221 <- data.frame(c(df_180221$Ongoing, df_180221$LightON), condition, nodeid, rep(regions_180221,2));
colnames(dat_180221) <- c("localeff", "condition", "nodeid", "region")
dat_180221$condition <- factor(dat_180221$condition, levels = c("Ongoing", "LightON"))
model_180221 <- lm(formula = "localeff ~ condition * region + (1|nodeid)", data=dat_180221)
res_aov_180221 <- aov(model_180221)
summary(res_aov_180221)
cm=compare_means(localeff ~ condition, data = dat_180221, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/localeff_2way_180221.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180221, x="region", y="localeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180221") + ylab("localeff")
dev.off()


df_180228 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180228_localeff")
regions_180228 <- character(nrow(df_180228));
regions_180228[c(1,4,5,6, 24, 25, 74)] <- "ventromedial thalamic nucleus"
regions_180228[c(2, 11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23,26, 40, 42, 62, 63, 64)] <- "posterior hypothalamic area"
regions_180228[c(3,7,8,10, 19, 20, 27, 28, 29, 30, 31, 54, 56, 57, 58, 59, 60, 61, 65, 66,67,68,69, 73, 75, 76, 77,78,79)] <- "dorsomedial nucleus of hypothalamus"
regions_180228[9] <- "ventromedial hypothalamic nucleus, dorsomedial part"
regions_180228[c(32, 33, 34, 35, 36, 37, 38,39, 41)] <- "arcuate nucleus of hypothalamus"
regions_180228[c(43,44,45)] <- "subincertal nucleus"
regions_180228[c(46, 47, 48, 49, 50)] <- "ventromedial hypothalamic nucleus"
regions_180228[c(51,52,53,55,70,71, 72)] <- "zona incerta"

nodeid <- rep(1:nrow(df_180228),2);
condition <- c(rep("Ongoing", nrow(df_180228)), rep("LightON",nrow(df_180228)));
dat_180228 <- data.frame(c(df_180228$Ongoing, df_180228$LightON), condition, nodeid, rep(regions_180228,2));
colnames(dat_180228) <- c("localeff", "condition", "nodeid", "region")
dat_180228$condition <- factor(dat_180228$condition, levels = c("Ongoing", "LightON"))
model_180228 <- lm(formula = "localeff ~ condition * region + (1|nodeid)", data=dat_180228)
res_aov_180228 <- aov(model_180228)
summary(res_aov_180228)
cm=compare_means(localeff ~ condition, data = dat_180228, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/localeff_2way_180228.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180228, x="region", y="localeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180228") + ylab("localeff")
dev.off()

df_180302 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180302_localeff")
regions_180302 <- character(nrow(df_180302));
regions_180302[c(1,11,12)] <- "posterior hypothalamic area"
regions_180302[c(2,3,4,5,6, 13, 14, 15,30,36)] <- "dorsal tuberomammillary nucleus"
regions_180302[c(7, 8,16, 17, 23, 31, 37, 38,39,40, 47)] <- "posterior hypothalamic area"
regions_180302[c(9, 10, 41,42)] <- "premammillary nucleus dorsal"
regions_180302[c(18, 19, 20, 21)] <- "supramammillary nucleus"
regions_180302[c(22, 24, 25, 29)] <- "medial magnocellular hypothalamic area"
regions_180302[c(26, 27, 28, 43,44,45,46)] <- "arcuate nucleus of hypothalamus"
regions_180302[c(32,33,34,35)] <- "dorsomedial hypothalamic nucleus, compact part"

nodeid <- rep(1:nrow(df_180302),2);
condition <- c(rep("Ongoing", nrow(df_180302)), rep("LightON",nrow(df_180302)));
dat_180302 <- data.frame(c(df_180302$Ongoing, df_180302$LightON), condition, nodeid, rep(regions_180302,2));
colnames(dat_180302) <- c("localeff", "condition", "nodeid", "region")
dat_180302$condition <- factor(dat_180302$condition, levels = c("Ongoing", "LightON"))
model_180302 <- lm(formula = "localeff ~ condition * region + (1|nodeid)", data=dat_180302)
res_aov_180302 <- aov(model_180302)
summary(res_aov_180302)
cm=compare_means(localeff ~ condition, data = dat_180302, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/localeff_2way_180302.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180302, x="region", y="localeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180302") + ylab("localeff")
dev.off()

df_180419 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180419_localeff")
regions_180419 <- character(nrow(df_180419));
regions_180419[c(1,2,11,12)] <- "ventromedial hypothalamic nucleus, dorsomedial part"
regions_180419[c(3, 40, 41, 42, 43, 44,45,46,47,48,49,50)] <- "dorsomedial nucleus of hypothalamus"
regions_180419[c(4,5,6,7,8)] <- "submedius thalamic nucleus"
regions_180419[c(9, 10)] <- "ventromedial thalamic nucleus"
regions_180419[c(13,14,15,16,17,18,23)] <- "dorsomedial hypothalamic nucleus, ventral part"
regions_180419[c(19,20,21,22,24,25,26,27,28,29,30,31,32,33,34,35)] <- "dorsomedial hypothalamic nucleus, dorsal part"
regions_180419[c(36,37,38,39,51,52,53,54,55,56,57,58,59,60,61)] <- "posterior hypothalamic area"

nodeid <- rep(1:nrow(df_180419),2);
condition <- c(rep("Ongoing", nrow(df_180419)), rep("LightON",nrow(df_180419)));
dat_180419 <- data.frame(c(df_180419$Ongoing, df_180419$LightON), condition, nodeid, rep(regions_180419,2));
colnames(dat_180419) <- c("localeff", "condition", "nodeid", "region")
dat_180419$condition <- factor(dat_180419$condition, levels = c("Ongoing", "LightON"))
model_180419 <- lm(formula = "localeff ~ condition * region + (1|nodeid)", data=dat_180419)
res_aov_180419 <- aov(model_180419)
summary(res_aov_180419)
cm=compare_means(localeff ~ condition, data = dat_180419, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/localeff_2way_180419.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180419, x="region", y="localeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180419") + ylab("localeff")
dev.off()

df_180420 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180420_localeff")
regions_180420 <- character(nrow(df_180420));
regions_180420[1] <- "dorsomedial hypothalamic nucleus, compact part"
regions_180420[c(2, 3, 4, 5, 6,7,8,9,10,11,12,13,14,15,16,17,18,64,65,66,67,68,69,70,71)] <- "posterior hypothalamic area"
regions_180420[c(19, 20, 21, 22, 23, 24, 25,26,27,28,29,30,31)] <- "premammillary nucleus dorsal"
regions_180420[c(32, 33, 34, 35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50)] <- "medial magnocellular hypothalamic area"
regions_180420[c(51,52,53,54,55,56,57,58,59,60,61,62,63)] <- "supramammillary nucleus, medial part"

nodeid <- rep(1:nrow(df_180420),2);
condition <- c(rep("Ongoing", nrow(df_180420)), rep("LightON",nrow(df_180420)));
dat_180420 <- data.frame(c(df_180420$Ongoing, df_180420$LightON), condition, nodeid, rep(regions_180420,2));
colnames(dat_180420) <- c("localeff", "condition", "nodeid", "region")
dat_180420$condition <- factor(dat_180420$condition, levels = c("Ongoing", "LightON"))
model_180420 <- lm(formula = "localeff ~ condition * region + (1|nodeid)", data=dat_180420)
res_aov_180420 <- aov(model_180420)
summary(res_aov_180420)
cm=compare_means(localeff ~ condition, data = dat_180420, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/localeff_2way_180420.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180420, x="region", y="localeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180420") + ylab("localeff")
dev.off()

df_180423 <- read_pzfx(path = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/analysis_network.pzfx", table="180423_localeff")
regions_180423 <- character(nrow(df_180423));
regions_180423[c(1,2,3,4)] <- "ventromedial hypothalamic nucleus"
regions_180423[5] <- "lateral hypothalamic area"
regions_180423[c(6,7,8,9,10,18,19,20,26,27,28,29,30,31,32,33,34,35,36)] <- "dorsomedial nucleus of hypothalamus"
regions_180423[c(11, 12,13,14,15,16,17,21,22,23,24,25,51,52,53,54,55,56)] <- "posterior hypothalamic area"
regions_180423[37] <- "premammillary nucleus, ventral part"
regions_180423[c(38,39,40,41,42,43,44,45,46,47,48,49,50)] <- "premammillary nucleus, dorsal part"

nodeid <- rep(1:nrow(df_180423),2);
condition <- c(rep("Ongoing", nrow(df_180423)), rep("LightON",nrow(df_180423)));
dat_180423 <- data.frame(c(df_180423$Ongoing, df_180423$LightON), condition, nodeid, rep(regions_180423,2));
colnames(dat_180423) <- c("localeff", "condition", "nodeid", "region")
dat_180423$condition <- factor(dat_180423$condition, levels = c("Ongoing", "LightON"))
model_180423 <- lm(formula = "localeff ~ condition * region + (1|nodeid)", data=dat_180423)
res_aov_180423 <- aov(model_180423)
summary(res_aov_180423)
cm=compare_means(localeff ~ condition, data = dat_180423, group.by = "region")

png(filename = "/Users/antoniogiulianozippo/Library/CloudStorage/OneDrive-CNR/tim_brown/localeff_2way_180423.png", width = 1860, height = 1440, units = "px")
ggboxplot(dat_180423, x="region", y="localeff", fill="condition") + stat_compare_means(aes(group=condition), label = "p.signif", size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + theme(text = element_text(size = 28)) + 
  theme(legend.title = element_blank()) + xlab("") + ggtitle("180423") + ylab("localeff")
dev.off()


######################################################################################################################
adj_pearson_171019 <- read.csv("~/Library/CloudStorage/OneDrive-CNR/tim_brown/adj_pearson_171019.csv", header=FALSE)
adj_pearson_171019 <- as.matrix(adj_pearson_171019)
plot(adj_pearson_171019)
net <- network(adj_pearson_171019, loop=TRUE)
plot(net)

neuron_rec_id <- c(126,31,34,33,35,24,31,54,79,47,61,71,56)



