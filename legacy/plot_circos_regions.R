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
library(ComplexHeatmap)

df_171019 <- read_pzfx(path = "/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171019_degree")
regions_171019 <- character(nrow(df_171019));
regions_171019[c(1,2,3,4,5,6,7,8,9,10, 12, 16, 19, 21, 23, 24, 26, 31, 33,43, 44, 46, 48, 49, 51, 52, 53, 54, 80, 81, 82, 83, 84, 85, 86, 87, 88)] = "interanteromedial thalamic nucleus"
regions_171019[c(11, 13, 14, 15, 17, 18, 20, 22, 25, 27, 28, 29, 30, 32, 34,35, 36, 37,38, 39, 40, 41, 42,45, 47, 50, 55, 56, 57, 58)] = "reuniens thalamic nucleus"
regions_171019[c(59, 60, 61, 62, 66, 67, 68, 69, 70, 71, 72, 73, 74)] = "mediodorsal thalamic nucleus"
regions_171019[c(63, 64, 65, 75, 76, 77, 78, 79, 95, 96, 97, 98, 99, 109, 110, 111, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126)] = "submedius thalamic nucleus"
regions_171019[c(89, 90, 91, 92, 93, 94, 100, 101, 102, 103, 104, 105, 106, 107, 108, 112, 113, 114)] = "central medial nucleus"
regions_171019[c(115,116)] = "paraventricular hypothalamic area, medial, parvicellular part"

df_171207 <- read_pzfx(path = "/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171207_degree")
regions_171207 <- character(nrow(df_171207));
regions_171207[c(1,2,3, 12,13,28)] <- "reuniens thalamic nucleus"
regions_171207[c(4,5,11,14,15,16)] <- "paraventricular hypothalamic area, medial, magnocellular part"
regions_171207[c(6,7)] <- "paraventricular hypothalamic area, lateral magnocellular part"
regions_171207[c(8,9,10,17)] <- "anterior nucleus of hypothalamus"
regions_171207[c(18,23,24,25,26,27)] <- "submedius thalamic nucleus"
regions_171207[c(19,20,21,22)] <- "dorsomedial nucleus of hypothalamus"
regions_171207[c(29,30)] <- "paraventricular hypothalamic nucleus, posterior"
regions_171207[31] <- "tuber cinereum area"

df_171208 <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171208_degree")
regions_171208 <- character(nrow(df_171208));
regions_171208[c(1,2,3,4,5,7,20,21,22,23,24,25,28)] <- "submedius thalamic nucleus"
regions_171208[6] <- "ventromedial hypothalamic nucleus"
regions_171208[c(8,9,11,13,16,17)] <- "ventromedial hypothalamic nucleus, dorsomedial part"
regions_171208[c(10,12,32)] <- "paraventricular hypothalamic nucleus, posterior"
regions_171208[c(14,15,18,19)] <- "dorsomedial nucleus of hypothalamus"
regions_171208[c(26, 27, 29, 30, 31,32,33, 34)] <- "anterior nucleus of hypothalamus"


df_171213 <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171213_degree")
regions_171213 <- character(nrow(df_171213));
regions_171213[c(1,2,3,4,5,6,7,8,9,10,14,15,16,17,18,19,28,33)] <- "paraventricular hypothalamic area, medial, parvicellular part"
regions_171213[c(11,12,13,20,21,22,23,24,25,26,27,29,30,31,32)] <- "reuniens thalamic nucleus"


df_180110 <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180110_degree")
regions_180110 <- character(nrow(df_180110));
regions_180110[c(1,2,11,18,23)] <- "reuniens thalamic nucleus"
regions_180110[c(3,4,5,6,9,10)] <- "paraventricular hypothalamic area, medial, magnocellular part"
regions_180110[7] <- "paraventricular hypothalamic area, lateral magnocellular part"
regions_180110[8] <- "subparaventricular zone"
regions_180110[c(12,13,14,15,16,17)] <- "dorsomedial nucleus of hypothalamus"
regions_180110[c(19,22,24,25,26,27,28,29,30,31,32,33,34,35)] <- "paraventricular hypothalamic nucleus, posterior"
regions_180110[c(20,21)] <-"submedius thalamic nucleus"

df_180111 <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180111_degree")
regions_180111 <- character(nrow(df_180111));
regions_180111[c(1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,22,23,24)] <- "reuniens thalamic nucleus"
regions_180111[14] <- "paraventricular hypothalamic nucleus, posterior"
regions_180111[c(19,20,21)] <- "ventromedial hypothalamic nucleus, dorsomedial part"

df_180131 <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180131_degree")
regions_180131 <- character(nrow(df_180131));
regions_180131[c(1,2,4,5,6,14,15,16)] <- "subparafascicular thalamic nucleus"
regions_180131[c(3,13,18,19,20,21,22,23,24, 27,29)] <- "posterior hypothalamic area"
regions_180131[c(7,8,9,10,11,12,25)] <- "dorsomedial hypothalamic nucleus, compact part"
regions_180131[17] <- "parvicellular part of ventral posteromedial nucleus"
regions_180131[26] <- "reuniens thalamic nucleus"
regions_180131[c(28,30,31)] <- "dorsomedial hypothalamic nucleus, dorsal part"

df_180221 <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180221_degree")
regions_180221 <- character(nrow(df_180221));
regions_180221[c(1, 35, 36, 48,49, 50)] <- "submedius thalamic nucleus"
regions_180221[c(2,3,4,12,17,18, 21, 22, 26, 31, 32, 33, 39, 40, 41, 42, 43, 44, 45, 46, 47)] <- "reuniens thalamic nucleus"
regions_180221[c(5,6,7,8,9,10,11,13,14,15,16, 27, 28, 29, 30, 37, 38)] <- "dorsomedial nucleus of hypothalamus"
regions_180221[c(19, 20, 34, 51, 52, 53, 54)] <- "parvicellular part of ventral posteromedial nucleus"
regions_180221[c(23, 24, 25)] <- "ventromedial hypothalamic nucleus, dorsomedial part"

df_180228 <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180228_degree")
regions_180228 <- character(nrow(df_180228));
regions_180228[c(1,4,5,6, 24, 25, 74)] <- "ventromedial thalamic nucleus"
regions_180228[c(2, 11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23,26, 40, 42, 62, 63, 64)] <- "posterior hypothalamic area"
regions_180228[c(3,7,8,10, 19, 20, 27, 28, 29, 30, 31, 54, 56, 57, 58, 59, 60, 61, 65, 66,67,68,69, 73, 75, 76, 77,78,79)] <- "dorsomedial nucleus of hypothalamus"
regions_180228[9] <- "ventromedial hypothalamic nucleus, dorsomedial part"
regions_180228[c(32, 33, 34, 35, 36, 37, 38,39, 41)] <- "arcuate nucleus of hypothalamus"
regions_180228[c(43,44,45)] <- "subincertal nucleus"
regions_180228[c(46, 47, 48, 49, 50)] <- "ventromedial hypothalamic nucleus"
regions_180228[c(51,52,53,55,70,71, 72)] <- "zona incerta"

df_180302 <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180302_degree")
regions_180302 <- character(nrow(df_180302));
regions_180302[c(1,11,12)] <- "posterior hypothalamic area"
regions_180302[c(2,3,4,5,6, 13, 14, 15,30,36)] <- "dorsal tuberomammillary nucleus"
regions_180302[c(7, 8,16, 17, 23, 31, 37, 38,39,40, 47)] <- "posterior hypothalamic area"
regions_180302[c(9, 10, 41,42)] <- "premammillary nucleus dorsal"
regions_180302[c(18, 19, 20, 21)] <- "supramammillary nucleus"
regions_180302[c(22, 24, 25, 29)] <- "medial magnocellular hypothalamic area"
regions_180302[c(26, 27, 28, 43,44,45,46)] <- "arcuate nucleus of hypothalamus"
regions_180302[c(32,33,34,35)] <- "dorsomedial hypothalamic nucleus, compact part"

df_180419 <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180419_degree")
regions_180419 <- character(nrow(df_180419));
regions_180419[c(1,2,11,12)] <- "ventromedial hypothalamic nucleus, dorsomedial part"
regions_180419[c(3, 40, 41, 42, 43, 44,45,46,47,48,49,50)] <- "dorsomedial nucleus of hypothalamus"
regions_180419[c(4,5,6,7,8)] <- "submedius thalamic nucleus"
regions_180419[c(9, 10)] <- "ventromedial thalamic nucleus"
regions_180419[c(13,14,15,16,17,18,23)] <- "dorsomedial hypothalamic nucleus, ventral part"
regions_180419[c(19,20,21,22,24,25,26,27,28,29,30,31,32,33,34,35)] <- "dorsomedial hypothalamic nucleus, dorsal part"
regions_180419[c(36,37,38,39,51,52,53,54,55,56,57,58,59,60,61)] <- "posterior hypothalamic area"

df_180420 <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180420_degree")
regions_180420 <- character(nrow(df_180420));
regions_180420[1] <- "dorsomedial hypothalamic nucleus, compact part"
regions_180420[c(2, 3, 4, 5, 6,7,8,9,10,11,12,13,14,15,16,17,18,64,65,66,67,68,69,70,71)] <- "posterior hypothalamic area"
regions_180420[c(19, 20, 21, 22, 23, 24, 25,26,27,28,29,30,31)] <- "premammillary nucleus dorsal"
regions_180420[c(32, 33, 34, 35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50)] <- "medial magnocellular hypothalamic area"
regions_180420[c(51,52,53,54,55,56,57,58,59,60,61,62,63)] <- "supramammillary nucleus, medial part"

df_180423 <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180423_degree")
regions_180423 <- character(nrow(df_180423));
regions_180423[c(1,2,3,4)] <- "ventromedial hypothalamic nucleus"
regions_180423[5] <- "lateral hypothalamic area"
regions_180423[c(6,7,8,9,10,18,19,20,26,27,28,29,30,31,32,33,34,35,36)] <- "dorsomedial nucleus of hypothalamus"
regions_180423[c(11, 12,13,14,15,16,17,21,22,23,24,25,51,52,53,54,55,56)] <- "posterior hypothalamic area"
regions_180423[37] <- "premammillary nucleus, ventral part"
regions_180423[c(38,39,40,41,42,43,44,45,46,47,48,49,50)] <- "premammillary nucleus, dorsal part"


####################################################################################
######################## 171019 #########################
adj_pearson_171019 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_ongoing_171019.csv", header=FALSE)
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
write.csv(conn_171019, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_171019.csv");
sel_mat <- which(conn_171019 > quantile(conn_171019, .75), arr.ind = TRUE);
net_171019 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_171019.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_171019, label = unique(rownames(sel_mat)))
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171019_fr")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171019_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171019_degree")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171019_localeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_171019)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_171019, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_ongoing_171019.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[1:18], track.index = c(2,2), text = "CM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[19:55], track.index = c(2,2), text = "IAM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[56:68], track.index = c(2,2), text = "MD", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[69:70], track.index = c(2,2), text = "PaMP", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[71:100], track.index = c(2,2), text = "Re", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[101:126], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))


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
adj_pearson_171019 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_evoked_171019.csv", header=FALSE)
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
write.csv(conn_171019, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_171019.csv");
sel_mat <- which(conn_171019 > quantile(conn_171019, .7), arr.ind = TRUE);
net_171019 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_171019.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_171019, label = unique(rownames(sel_mat)))
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171019_fr")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171019_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171019_degree")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171019_localeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$LightON)

tt = c()
for (i in 1:dim(adj_pearson_171019)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_171019, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_evoked_171019.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[1:18], track.index = c(2,2), text = "CM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[19:55], track.index = c(2,2), text = "IAM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[56:68], track.index = c(2,2), text = "MD", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[69:70], track.index = c(2,2), text = "PaMP", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[71:100], track.index = c(2,2), text = "Re", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[101:126], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))


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
adj_pearson_171207 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_ongoing_171207.csv", header=FALSE)
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
write.csv(conn_171207, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_171207.csv");
sel_mat <- which(conn_171207 > quantile(conn_171207, .7), arr.ind = TRUE);
net_171207 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_171207.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_171207, label = unique(rownames(sel_mat)))
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171207_fr")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171207_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171207_degree")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171207_localeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_171207)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_171207, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_ongoing_171207.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[1:4], track.index = c(2,2), text = "AHC", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[5:8], track.index = c(2,2), text = "DM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[9:10], track.index = c(2,2), text = "PaLM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[11:16], track.index = c(2,2), text = "PaMM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[17:18], track.index = c(2,2), text = "PaPo", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[19:24], track.index = c(2,2), text = "Re", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[25:30], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))




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
adj_pearson_171207 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_evoked_171207.csv", header=FALSE)
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
write.csv(conn_171207, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_171207.csv");
sel_mat <- which(conn_171207 > quantile(conn_171207, .7), arr.ind = TRUE);
net_171207 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_171207.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_171207, label = unique(rownames(sel_mat)))
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171207_fr")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171207_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171207_degree")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171207_localeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$LightON)

tt = c()
for (i in 1:dim(adj_pearson_171207)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_171207, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_evoked_171207.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[1:4], track.index = c(2,2), text = "AHC", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[5:8], track.index = c(2,2), text = "DM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[9:10], track.index = c(2,2), text = "PaLM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[11:16], track.index = c(2,2), text = "PaMM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[17:18], track.index = c(2,2), text = "PaPo", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[19:24], track.index = c(2,2), text = "Re", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[25:30], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))



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
adj_pearson_171208 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_ongoing_171208.csv", header=FALSE)
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
write.csv(conn_171208, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_171208.csv");
sel_mat <- which(conn_171208 > quantile(conn_171208, .75), arr.ind = TRUE);
net_171208 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_171208.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_171208, label = unique(rownames(sel_mat)))
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171208_fr")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171208_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171208_degree")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171208_localeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_171208)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_171208, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_ongoing_171208.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[1:8], track.index = c(2,2), text = "AHC", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[9:12], track.index = c(2,2), text = "DM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[13:14], track.index = c(2,2), text = "PaPo", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[15:27], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[28:34], track.index = c(2,2), text = "VMHDM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[19:24], track.index = c(2,2), text = "Re", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[25:30], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))




for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_171208);
idd = which(adj_pearson_171208 > quantile(adj_pearson_171208, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_171208)[1]){
  for(jj in 1:dim(adj_pearson_171208)[2]){
    if (adj_pearson_171208[ii,jj]> quantile(adj_pearson_171208, 0.98)){
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
adj_pearson_171208 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_evoked_171208.csv", header=FALSE)
adj_pearson_171208 <- as.matrix(adj_pearson_171208);
colnames(adj_pearson_171208) <- regions_171208;
rownames(adj_pearson_171208) <- regions_171208;
plot(adj_pearson_171208)
net <- network(adj_pearson_171208, loop=TRUE)
plot(net)

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
write.csv(conn_171208, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_171208.csv");
sel_mat <- which(conn_171208 > quantile(conn_171208, .7), arr.ind = TRUE);
net_171208 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_171208.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_171208, label = unique(rownames(sel_mat)))
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171208_fr")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171208_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171208_degree")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171208_localeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$LightON)

tt = c()
for (i in 1:dim(adj_pearson_171208)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_171208, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_evoked_171208.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[1:8], track.index = c(2,2), text = "AHC", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[9:12], track.index = c(2,2), text = "DM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[13:14], track.index = c(2,2), text = "PaPo", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[15:27], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[28:34], track.index = c(2,2), text = "VMHDM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[19:24], track.index = c(2,2), text = "Re", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[25:30], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))



for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_171208);
idd = which(adj_pearson_171208 > quantile(adj_pearson_171208, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_171208)[1]){
  for(jj in 1:dim(adj_pearson_171208)[2]){
    if (adj_pearson_171208[ii,jj]> quantile(adj_pearson_171208, 0.98)){
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
adj_pearson_171213 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_ongoing_171213.csv", header=FALSE)
adj_pearson_171213 <- as.matrix(adj_pearson_171213);
colnames(adj_pearson_171213) <- regions_171213;
rownames(adj_pearson_171213) <- regions_171213;
plot(adj_pearson_171213)
net <- network(adj_pearson_171213, loop=TRUE)
plot(net)

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
write.csv(conn_171213, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_171213.csv");
sel_mat <- which(conn_171213 > quantile(conn_171213, .7), arr.ind = TRUE);
net_171213 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_171213.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_171213, label = unique(rownames(sel_mat)))
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171213_fr")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171213_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171213_degree")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171213_localeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_171213)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_171213, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_ongoing_171213.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[1:18], track.index = c(2,2), text = "PaMP", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[19:33], track.index = c(2,2), text = "Re", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[13:14], track.index = c(2,2), text = "PaPo", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[15:27], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[28:34], track.index = c(2,2), text = "VMHDM", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[19:24], track.index = c(2,2), text = "Re", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[25:30], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))




for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_171213);
idd = which(adj_pearson_171213 > quantile(adj_pearson_171213, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_171213)[1]){
  for(jj in 1:dim(adj_pearson_171213)[2]){
    if (adj_pearson_171213[ii,jj]> quantile(adj_pearson_171213, 0.98)){
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
adj_pearson_171213 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_evoked_171213.csv", header=FALSE)
adj_pearson_171213 <- as.matrix(adj_pearson_171213);
colnames(adj_pearson_171213) <- regions_171213;
rownames(adj_pearson_171213) <- regions_171213;
plot(adj_pearson_171213)
net <- network(adj_pearson_171213, loop=TRUE)
plot(net)

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
write.csv(conn_171213, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_171213.csv");
sel_mat <- which(conn_171213 > quantile(conn_171213, .7), arr.ind = TRUE);
net_171213 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_171213.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_171213, label = unique(rownames(sel_mat)))
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171213_fr")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171213_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171213_degree")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="171213_localeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$LightON)

tt = c()
for (i in 1:dim(adj_pearson_171213)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_171213, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_evoked_171213.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[1:18], track.index = c(2,2), text = "PaMP", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[19:33], track.index = c(2,2), text = "Re", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[13:14], track.index = c(2,2), text = "PaPo", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[15:27], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[28:34], track.index = c(2,2), text = "VMHDM", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[19:24], track.index = c(2,2), text = "Re", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[25:30], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))



for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_171213);
idd = which(adj_pearson_171213 > quantile(adj_pearson_171213, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_171213)[1]){
  for(jj in 1:dim(adj_pearson_171213)[2]){
    if (adj_pearson_171213[ii,jj]> quantile(adj_pearson_171213, 0.98)){
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
adj_pearson_180110 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_ongoing_180110.csv", header=FALSE)
adj_pearson_180110 <- as.matrix(adj_pearson_180110);
colnames(adj_pearson_180110) <- regions_180110;
rownames(adj_pearson_180110) <- regions_180110;
plot(adj_pearson_180110)
net <- network(adj_pearson_180110, loop=TRUE)
plot(net)

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
write.csv(conn_180110, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_180110.csv");
sel_mat <- which(conn_180110 > quantile(conn_180110, .7), arr.ind = TRUE);
net_180110 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_180110.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180110, label = unique(rownames(sel_mat)))
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180110_fr")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180110_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180110_degree")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180110_localeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180110)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180110, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_ongoing_180110.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[1:6], track.index = c(2,2), text = "DM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[7], track.index = c(2,2), text = "PaLM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[8:13], track.index = c(2,2), text = "PaMM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[14:27], track.index = c(2,2), text = "PaPo", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[28:32], track.index = c(2,2), text = "Re", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[33:35], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[25:30], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))



for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180110);
idd = which(adj_pearson_180110 > quantile(adj_pearson_180110, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180110)[1]){
  for(jj in 1:dim(adj_pearson_180110)[2]){
    if (adj_pearson_180110[ii,jj]> quantile(adj_pearson_180110, 0.98)){
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
adj_pearson_180110 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_evoked_180110.csv", header=FALSE)
adj_pearson_180110 <- as.matrix(adj_pearson_180110);
colnames(adj_pearson_180110) <- regions_180110;
rownames(adj_pearson_180110) <- regions_180110;
plot(adj_pearson_180110)
net <- network(adj_pearson_180110, loop=TRUE)
plot(net)

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
write.csv(conn_180110, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_180110.csv");
sel_mat <- which(conn_180110 > quantile(conn_180110, .75), arr.ind = TRUE);
net_180110 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_180110.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180110, label = unique(rownames(sel_mat)))
dev.off()


col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180110_fr")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180110_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180110_degree")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180110_localeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$LightON)

tt = c()
for (i in 1:dim(adj_pearson_180110)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180110, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_evoked_180110.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[1:6], track.index = c(2,2), text = "DM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[7], track.index = c(2,2), text = "PaLM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[8:13], track.index = c(2,2), text = "PaMM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[14:27], track.index = c(2,2), text = "PaPo", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[28:32], track.index = c(2,2), text = "Re", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[33:35], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[25:30], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))



for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180110);
idd = which(adj_pearson_180110 > quantile(adj_pearson_180110, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180110)[1]){
  for(jj in 1:dim(adj_pearson_180110)[2]){
    if (adj_pearson_180110[ii,jj]> quantile(adj_pearson_180110, 0.98)){
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
adj_pearson_180111 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_ongoing_180111.csv", header=FALSE)
adj_pearson_180111 <- as.matrix(adj_pearson_180111);
colnames(adj_pearson_180111) <- regions_180111;
rownames(adj_pearson_180111) <- regions_180111;
#plot(adj_pearson_180111)
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
write.csv(conn_180111, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_180111.csv");
sel_mat <- which(conn_180111 > quantile(conn_180111, .7), arr.ind = TRUE);
net_180111 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_180111.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180111, label = unique(rownames(sel_mat)))
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180111_fr")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180111_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180111_degree")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180111_localeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180111)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180111, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_ongoing_180111.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[1], track.index = c(2,2), text = "PaPo", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[2:21], track.index = c(2,2), text = "Re", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[22:24], track.index = c(2,2), text = "VMHDM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[15:27], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[28:34], track.index = c(2,2), text = "VMHDM", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[33:35], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[25:30], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))



for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180111);
idd = which(adj_pearson_180111 > quantile(adj_pearson_180111, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180111)[1]){
  for(jj in 1:dim(adj_pearson_180111)[2]){
    if (adj_pearson_180111[ii,jj]> quantile(adj_pearson_180111, 0.98)){
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
adj_pearson_180111 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_evoked_180111.csv", header=FALSE)
adj_pearson_180111 <- as.matrix(adj_pearson_180111);
colnames(adj_pearson_180111) <- regions_180111;
rownames(adj_pearson_180111) <- regions_180111;
#plot(adj_pearson_180111)
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
write.csv(conn_180111, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_180111.csv");
sel_mat <- which(conn_180111 > quantile(conn_180111, .75), arr.ind = TRUE);
net_180111 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_180111.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180111, label = unique(rownames(sel_mat)))
dev.off()


col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180111_fr")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180111_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180111_degree")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180111_localeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$LightON)

tt = c()
for (i in 1:dim(adj_pearson_180111)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180111, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_evoked_180111.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[1], track.index = c(2,2), text = "PaPo", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[2:21], track.index = c(2,2), text = "Re", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[22:24], track.index = c(2,2), text = "VMHDM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[15:27], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[28:34], track.index = c(2,2), text = "VMHDM", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[33:35], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[25:30], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))



for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180111);
idd = which(adj_pearson_180111 > quantile(adj_pearson_180111, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180111)[1]){
  for(jj in 1:dim(adj_pearson_180111)[2]){
    if (adj_pearson_180111[ii,jj]> quantile(adj_pearson_180111, 0.98)){
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
adj_pearson_180131 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_ongoing_180131.csv", header=FALSE)
adj_pearson_180131 <- as.matrix(adj_pearson_180131);
colnames(adj_pearson_180131) <- regions_180131;
rownames(adj_pearson_180131) <- regions_180131;
#plot(adj_pearson_180131)
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
write.csv(conn_180131, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_180131.csv");
sel_mat <- which(conn_180131 > quantile(conn_180131, .7), arr.ind = TRUE);
net_180131 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_180131.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180131, label = unique(rownames(sel_mat)))
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180131_fr")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180131_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180131_degree")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180131_localeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180131)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180131, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_ongoing_180131.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[1:7], track.index = c(2,2), text = "DMC", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[8:10], track.index = c(2,2), text = "DMD", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[11], track.index = c(2,2), text = "Gus", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[12:22], track.index = c(2,2), text = "PH", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[23], track.index = c(2,2), text = "Re", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[24:31], track.index = c(2,2), text = "SPF", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[25:30], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))



for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180131);
idd = which(adj_pearson_180131 > quantile(adj_pearson_180131, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180131)[1]){
  for(jj in 1:dim(adj_pearson_180131)[2]){
    if (adj_pearson_180131[ii,jj]> quantile(adj_pearson_180131, 0.98)){
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
adj_pearson_180131 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_evoked_180131.csv", header=FALSE)
adj_pearson_180131 <- as.matrix(adj_pearson_180131);
colnames(adj_pearson_180131) <- regions_180131;
rownames(adj_pearson_180131) <- regions_180131;
#plot(adj_pearson_180131)
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
write.csv(conn_180131, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_180131.csv");
sel_mat <- which(conn_180131 > quantile(conn_180131, .75), arr.ind = TRUE);
net_180131 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_180131.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180131, label = unique(rownames(sel_mat)))
dev.off()


col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180131_fr")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180131_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180131_degree")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180131_localeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$LightON)

tt = c()
for (i in 1:dim(adj_pearson_180131)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180131, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_evoked_180131.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[1:7], track.index = c(2,2), text = "DMC", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[8:10], track.index = c(2,2), text = "DMD", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[11], track.index = c(2,2), text = "Gus", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[12:22], track.index = c(2,2), text = "PH", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[23], track.index = c(2,2), text = "Re", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[24:31], track.index = c(2,2), text = "SPF", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180131);
idd = which(adj_pearson_180131 > quantile(adj_pearson_180131, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180131)[1]){
  for(jj in 1:dim(adj_pearson_180131)[2]){
    if (adj_pearson_180131[ii,jj]> quantile(adj_pearson_180131, 0.98)){
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
adj_pearson_180221 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_ongoing_180221.csv", header=FALSE)
adj_pearson_180221 <- as.matrix(adj_pearson_180221);
colnames(adj_pearson_180221) <- regions_180221;
rownames(adj_pearson_180221) <- regions_180221;
#plot(adj_pearson_180221)
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
write.csv(conn_180221, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_180221.csv");
sel_mat <- which(conn_180221 > quantile(conn_180221, .7), arr.ind = TRUE);
net_180221 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_180221.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180221, label = unique(rownames(sel_mat)))
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180221_fr")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180221_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180221_degree")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180221_localeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180221)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180221, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_ongoing_180221.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[1:17], track.index = c(2,2), text = "DM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[18:24], track.index = c(2,2), text = "Gus", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[25:45], track.index = c(2,2), text = "Re", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[46:51], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[52:54], track.index = c(2,2), text = "VMHDM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[24:31], track.index = c(2,2), text = "SPF", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[25:30], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))



for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180221);
idd = which(adj_pearson_180221 > quantile(adj_pearson_180221, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180221)[1]){
  for(jj in 1:dim(adj_pearson_180221)[2]){
    if (adj_pearson_180221[ii,jj]> quantile(adj_pearson_180221, 0.98)){
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
adj_pearson_180221 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_evoked_180221.csv", header=FALSE)
adj_pearson_180221 <- as.matrix(adj_pearson_180221);
colnames(adj_pearson_180221) <- regions_180221;
rownames(adj_pearson_180221) <- regions_180221;
#plot(adj_pearson_180221)
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
write.csv(conn_180221, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_180221.csv");
sel_mat <- which(conn_180221 > quantile(conn_180221, .75), arr.ind = TRUE);
net_180221 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_180221.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180221, label = unique(rownames(sel_mat)))
dev.off()


col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180221_fr")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180221_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180221_degree")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180221_localeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$LightON)

tt = c()
for (i in 1:dim(adj_pearson_180221)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180221, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_evoked_180221.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[1:17], track.index = c(2,2), text = "DM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[18:24], track.index = c(2,2), text = "Gus", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[25:45], track.index = c(2,2), text = "Re", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[46:51], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[52:54], track.index = c(2,2), text = "VMHDM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[24:31], track.index = c(2,2), text = "SPF", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[25:30], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))



for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180221);
idd = which(adj_pearson_180221 > quantile(adj_pearson_180221, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180221)[1]){
  for(jj in 1:dim(adj_pearson_180221)[2]){
    if (adj_pearson_180221[ii,jj]> quantile(adj_pearson_180221, 0.98)){
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
adj_pearson_180228 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_ongoing_180228.csv", header=FALSE)
adj_pearson_180228 <- as.matrix(adj_pearson_180228);
colnames(adj_pearson_180228) <- regions_180228;
rownames(adj_pearson_180228) <- regions_180228;
#plot(adj_pearson_180228)
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
write.csv(conn_180228, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_180228.csv");
sel_mat <- which(conn_180228 > quantile(conn_180228, .7), arr.ind = TRUE);
net_180228 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_180228.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180228, label = unique(rownames(sel_mat)))
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180228_fr")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180228_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180228_degree")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180228_localeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180228)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180228, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_ongoing_180228.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[1:9], track.index = c(2,2), text = "Arc", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[10:38], track.index = c(2,2), text = "DM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[39:56], track.index = c(2,2), text = "PH", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[57:59], track.index = c(2,2), text = "SubI", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[60:65], track.index = c(2,2), text = "VMH", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[66:72], track.index = c(2,2), text = "VM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[73:79], track.index = c(2,2), text = "ZI", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))



for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180228);
idd = which(adj_pearson_180228 > quantile(adj_pearson_180228, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180228)[1]){
  for(jj in 1:dim(adj_pearson_180228)[2]){
    if (adj_pearson_180228[ii,jj]> quantile(adj_pearson_180228, 0.98)){
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
adj_pearson_180228 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_evoked_180228.csv", header=FALSE)
adj_pearson_180228 <- as.matrix(adj_pearson_180228);
colnames(adj_pearson_180228) <- regions_180228;
rownames(adj_pearson_180228) <- regions_180228;
#plot(adj_pearson_180228)
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
write.csv(conn_180228, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_180228.csv");
sel_mat <- which(conn_180228 > quantile(conn_180228, .75), arr.ind = TRUE);
net_180228 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_180228.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180228, label = unique(rownames(sel_mat)))
dev.off()


col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180228_fr")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180228_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180228_degree")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180228_localeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$LightON)

tt = c()
for (i in 1:dim(adj_pearson_180228)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180228, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_evoked_180228.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[1:9], track.index = c(2,2), text = "Arc", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[10:38], track.index = c(2,2), text = "DM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[39:56], track.index = c(2,2), text = "PH", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[57:59], track.index = c(2,2), text = "SubI", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[60:65], track.index = c(2,2), text = "VMH", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[66:72], track.index = c(2,2), text = "VM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[73:79], track.index = c(2,2), text = "ZI", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))



for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180228);
idd = which(adj_pearson_180228 > quantile(adj_pearson_180228, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180228)[1]){
  for(jj in 1:dim(adj_pearson_180228)[2]){
    if (adj_pearson_180228[ii,jj]> quantile(adj_pearson_180228, 0.98)){
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
adj_pearson_180302 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_ongoing_180302.csv", header=FALSE)
adj_pearson_180302 <- as.matrix(adj_pearson_180302);
colnames(adj_pearson_180302) <- regions_180302;
rownames(adj_pearson_180302) <- regions_180302;
#plot(adj_pearson_180302)
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
write.csv(conn_180302, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_180302.csv");
sel_mat <- which(conn_180302 > quantile(conn_180302, .7), arr.ind = TRUE);
net_180302 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_180302.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180302, label = unique(rownames(sel_mat)))
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180302_fr")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180302_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180302_degree")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180302_localeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180302)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180302, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_ongoing_180302.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[1:7], track.index = c(2,2), text = "Arc", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[8:17], track.index = c(2,2), text = "DTM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[18:21], track.index = c(2,2), text = "DMC", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[22:25], track.index = c(2,2), text = "MM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[26:39], track.index = c(2,2), text = "PH", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[40:43], track.index = c(2,2), text = "PMD", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[44:47], track.index = c(2,2), text = "SuM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))



for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180302);
idd = which(adj_pearson_180302 > quantile(adj_pearson_180302, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180302)[1]){
  for(jj in 1:dim(adj_pearson_180302)[2]){
    if (adj_pearson_180302[ii,jj]> quantile(adj_pearson_180302, 0.98)){
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
adj_pearson_180302 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_evoked_180302.csv", header=FALSE)
adj_pearson_180302 <- as.matrix(adj_pearson_180302);
colnames(adj_pearson_180302) <- regions_180302;
rownames(adj_pearson_180302) <- regions_180302;
#plot(adj_pearson_180302)
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
write.csv(conn_180302, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_180302.csv");
sel_mat <- which(conn_180302 > quantile(conn_180302, .75), arr.ind = TRUE);
net_180302 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_180302.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180302, label = unique(rownames(sel_mat)))
dev.off()


col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180302_fr")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180302_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180302_degree")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180302_localeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$LightON)

tt = c()
for (i in 1:dim(adj_pearson_180302)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180302, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_evoked_180302.png"), units = "in",  width = 16, height = 16, res = 600)
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


highlight.sector(tt[1:7], track.index = c(2,2), text = "Arc", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[8:17], track.index = c(2,2), text = "DTM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[18:21], track.index = c(2,2), text = "DMC", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[22:25], track.index = c(2,2), text = "MM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[26:39], track.index = c(2,2), text = "PH", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[40:43], track.index = c(2,2), text = "PMD", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[44:47], track.index = c(2,2), text = "SuM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))



for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180302);
idd = which(adj_pearson_180302 > quantile(adj_pearson_180302, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180302)[1]){
  for(jj in 1:dim(adj_pearson_180302)[2]){
    if (adj_pearson_180302[ii,jj]> quantile(adj_pearson_180302, 0.98)){
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
adj_pearson_180419 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_ongoing_180419.csv", header=FALSE)
adj_pearson_180419 <- as.matrix(adj_pearson_180419);
colnames(adj_pearson_180419) <- regions_180419;
rownames(adj_pearson_180419) <- regions_180419;
#plot(adj_pearson_180419)
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
write.csv(conn_180419, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_180419.csv");
sel_mat <- which(conn_180419 > quantile(conn_180419, .7), arr.ind = TRUE);
net_180419 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_180419.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180419, label = unique(rownames(sel_mat)))
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180419_fr")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180419_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180419_degree")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180419_localeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180419)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180419, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_ongoing_180419.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[1:16], track.index = c(2,2), text = "DMD", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[17:23], track.index = c(2,2), text = "DMV", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[24:35], track.index = c(2,2), text = "DM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[36:50], track.index = c(2,2), text = "PH", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[51:55], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[56:59], track.index = c(2,2), text = "VMHDM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[60:61], track.index = c(2,2), text = "VM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))



for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180419);
idd = which(adj_pearson_180419 > quantile(adj_pearson_180419, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180419)[1]){
  for(jj in 1:dim(adj_pearson_180419)[2]){
    if (adj_pearson_180419[ii,jj]> quantile(adj_pearson_180419, 0.98)){
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
adj_pearson_180419 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_evoked_180419.csv", header=FALSE)
adj_pearson_180419 <- as.matrix(adj_pearson_180419);
colnames(adj_pearson_180419) <- regions_180419;
rownames(adj_pearson_180419) <- regions_180419;
#plot(adj_pearson_180419)
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
write.csv(conn_180419, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_180419.csv");
sel_mat <- which(conn_180419 > quantile(conn_180419, .75), arr.ind = TRUE);
net_180419 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_180419.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180419, label = unique(rownames(sel_mat)))
dev.off()


col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180419_fr")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180419_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180419_degree")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180419_localeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$LightON)

tt = c()
for (i in 1:dim(adj_pearson_180419)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180419, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_evoked_180419.png"), units = "in",  width = 16, height = 16, res = 600)
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


highlight.sector(tt[1:16], track.index = c(2,2), text = "DMD", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[17:23], track.index = c(2,2), text = "DMV", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[24:35], track.index = c(2,2), text = "DM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[36:50], track.index = c(2,2), text = "PH", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[51:55], track.index = c(2,2), text = "Sub", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[56:59], track.index = c(2,2), text = "VMHDM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[60:61], track.index = c(2,2), text = "VM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))


for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180419);
idd = which(adj_pearson_180419 > quantile(adj_pearson_180419, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180419)[1]){
  for(jj in 1:dim(adj_pearson_180419)[2]){
    if (adj_pearson_180419[ii,jj]> quantile(adj_pearson_180419, 0.98)){
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
adj_pearson_180420 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_ongoing_180420.csv", header=FALSE)
adj_pearson_180420 <- as.matrix(adj_pearson_180420);
colnames(adj_pearson_180420) <- regions_180420;
rownames(adj_pearson_180420) <- regions_180420;
#plot(adj_pearson_180420)
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
write.csv(conn_180420, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_180420.csv");
sel_mat <- which(conn_180420 > quantile(conn_180420, .7), arr.ind = TRUE);
net_180420 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_180420.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180420, label = unique(rownames(sel_mat)))
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180420_fr")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180420_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180420_degree")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180420_localeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180420)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180420, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_ongoing_180420.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[2:20], track.index = c(2,2), text = "MM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[21:45], track.index = c(2,2), text = "PH", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[46:58], track.index = c(2,2), text = "PMD", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[59:71], track.index = c(2,2), text = "SuMM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[51:55], track.index = c(2,2), text = "SuMM", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[56:59], track.index = c(2,2), text = "VMHDM", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[60:61], track.index = c(2,2), text = "VM", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))



for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180420);
idd = which(adj_pearson_180420 > quantile(adj_pearson_180420, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180420)[1]){
  for(jj in 1:dim(adj_pearson_180420)[2]){
    if (adj_pearson_180420[ii,jj]> quantile(adj_pearson_180420, 0.98)){
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
adj_pearson_180420 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_evoked_180420.csv", header=FALSE)
adj_pearson_180420 <- as.matrix(adj_pearson_180420);
colnames(adj_pearson_180420) <- regions_180420;
rownames(adj_pearson_180420) <- regions_180420;
#plot(adj_pearson_180420)
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
write.csv(conn_180420, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_180420.csv");
sel_mat <- which(conn_180420 > quantile(conn_180420, .75), arr.ind = TRUE);
net_180420 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_180420.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180420, label = unique(rownames(sel_mat)))
dev.off()


col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180420_fr")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180420_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180420_degree")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180420_localeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$LightON)

tt = c()
for (i in 1:dim(adj_pearson_180420)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180420, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_evoked_180420.png"), units = "in",  width = 16, height = 16, res = 600)
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


highlight.sector(tt[2:20], track.index = c(2,2), text = "MM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[21:45], track.index = c(2,2), text = "PH", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[46:58], track.index = c(2,2), text = "PMD", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[59:71], track.index = c(2,2), text = "SuMM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[51:55], track.index = c(2,2), text = "SuMM", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[56:59], track.index = c(2,2), text = "VMHDM", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[60:61], track.index = c(2,2), text = "VM", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))


for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180420);
idd = which(adj_pearson_180420 > quantile(adj_pearson_180420, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180420)[1]){
  for(jj in 1:dim(adj_pearson_180420)[2]){
    if (adj_pearson_180420[ii,jj]> quantile(adj_pearson_180420, 0.98)){
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
adj_pearson_180423 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_ongoing_180423.csv", header=FALSE)
adj_pearson_180423 <- as.matrix(adj_pearson_180423);
colnames(adj_pearson_180423) <- regions_180423;
rownames(adj_pearson_180423) <- regions_180423;
#plot(adj_pearson_180423)
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
write.csv(conn_180423, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_180423.csv");
sel_mat <- which(conn_180423 > quantile(conn_180423, .7), arr.ind = TRUE);
net_180423 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_ongoing_180423.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180423, label = unique(rownames(sel_mat)))
dev.off()

col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180423_fr")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180423_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180423_degree")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$Ongoing)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180423_localeff")
col_fun_metric = colorRamp2( seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10), 
                             rev(brewer.pal(n = length(seq(min(df$Ongoing),max(df$Ongoing), (max(df$Ongoing) - min(df$Ongoing))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$Ongoing)

tt = c()
for (i in 1:dim(adj_pearson_180423)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180423, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_ongoing_180423.png"), units = "in",  width = 16, height = 16, res = 600)
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

highlight.sector(tt[1:19], track.index = c(2,2), text = "DM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[20], track.index = c(2,2), text = "LH", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[21:38], track.index = c(2,2), text = "PH", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[39:51], track.index = c(2,2), text = "PMD", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[52], track.index = c(2,2), text = "PMV", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[53:56], track.index = c(2,2), text = "VMH", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[60:61], track.index = c(2,2), text = "VM", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))



for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180423);
idd = which(adj_pearson_180423 > quantile(adj_pearson_180423, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180423)[1]){
  for(jj in 1:dim(adj_pearson_180423)[2]){
    if (adj_pearson_180423[ii,jj]> quantile(adj_pearson_180423, 0.98)){
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
adj_pearson_180423 <- read.csv("C:/Users/antonio/OneDrive - CNR/tim_brown/adj_pearson_evoked_180423.csv", header=FALSE)
adj_pearson_180423 <- as.matrix(adj_pearson_180423);
colnames(adj_pearson_180423) <- regions_180423;
rownames(adj_pearson_180423) <- regions_180423;
#plot(adj_pearson_180423)
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
write.csv(conn_180423, file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_180423.csv");
sel_mat <- which(conn_180423 > quantile(conn_180423, .75), arr.ind = TRUE);
net_180423 <- network(sel_mat, directed = FALSE)
png(file = "C:/Users/antonio/OneDrive - CNR/tim_brown/conn_evoked_180423.png", res=600,  units = "in",  width = 26, height = 16)
plot.network(net_180423, label = unique(rownames(sel_mat)))
dev.off()


col_fun = colorRamp2( seq(-0.6,0.6,0.12), rev(brewer.pal(n = length(seq(-0.6,0.6,0.12)), name = "RdBu")), transparency = 0)

df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180423_fr")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
fr_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180423_cluscoeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
clus_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180423_degree")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
deg_col = col_fun_metric(df$LightON)
df <- read_pzfx(path = "C:/Users/antonio/OneDrive - CNR/tim_brown/analysis_network.pzfx", table="180423_localeff")
col_fun_metric = colorRamp2( seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10), 
                             rev(brewer.pal(n = length(seq(min(df$LightON),max(df$LightON), (max(df$LightON) - min(df$LightON))/10)), name = "PRGn")), transparency = 0)
loceff_col = col_fun_metric(df$LightON)

tt = c()
for (i in 1:dim(adj_pearson_180423)[1]){
  tt = c(tt, toString(i))
}
sss=sort(regions_180423, index.return=TRUE)
tt = tt[sss$ix]

png(file = paste0("C:/Users/antonio/OneDrive - CNR/tim_brown/regionwise_evoked_180423.png"), units = "in",  width = 16, height = 16, res = 600)
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


highlight.sector(tt[1:19], track.index = c(2,2), text = "DM", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "cadetblue1", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[20], track.index = c(2,2), text = "LH", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "violet", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[21:38], track.index = c(2,2), text = "PH", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "green", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[39:51], track.index = c(2,2), text = "PMD", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "gold", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[52], track.index = c(2,2), text = "PMV", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "burlywood", padding = c(-0.6, 0, 0, 0))
highlight.sector(tt[53:56], track.index = c(2,2), text = "VMH", facing = "bending.inside", 
                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "blueviolet", padding = c(-0.6, 0, 0, 0))
#highlight.sector(tt[60:61], track.index = c(2,2), text = "VM", facing = "bending.inside", 
#                 niceFacing = TRUE, text.vjust = 0, cex = 1.4, col = "red", padding = c(-0.6, 0, 0, 0))



for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 0.8, col = "black", sector.index = si, track.index = 2)
}

values = array(adj_pearson_180423);
idd = which(adj_pearson_180423 > quantile(adj_pearson_180423, 0.98), arr.ind = FALSE);
min_values = min(values[idd])
svalues = sort(values, index.return = TRUE);
col_fun_edges = colorRamp2( seq(mean(values),max(values), (max(values) - mean(values))/10 ), 
                            rev(brewer.pal(n = length(seq(mean(values),max(values), (max(values) - mean(values))/10 )), name = "RdBu")), transparency = 0.5)

for(ii in 1:dim(adj_pearson_180423)[1]){
  for(jj in 1:dim(adj_pearson_180423)[2]){
    if (adj_pearson_180423[ii,jj]> quantile(adj_pearson_180423, 0.98)){
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
