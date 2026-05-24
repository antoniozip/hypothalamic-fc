library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(cartography)
library(ggpubr)
library(pals)

all_values_pos = c()
values_raw = array(conn_ongoing_day)
all_values_pos = c(all_values_pos, values_raw)
values_raw = array(conn_ongoing_night)
all_values_pos = c(all_values_pos, values_raw)

all_values_vs = c()
values_raw = array(conn_ongoing_day_vs_night)
all_values_vs = c(all_values_vs, values_raw)

values_raw = array(conn_evoked_day)
all_values_pos = c(all_values_pos, values_raw)
values_raw = array(conn_evoked_night)
all_values_pos = c(all_values_pos, values_raw)

values_raw = array(conn_evoked_day_vs_night)
all_values_vs = c(all_values_vs, values_raw)




depth_values = 19;
min_values_pos = min(all_values_pos);
max_values_pos = max(all_values_pos);
seq_values_pos = seq(min(all_values_pos), max_values_pos, (max_values_pos - min(all_values_pos))/depth_values);
col_fun_edges_pos = colorRamp2(seq_values_pos, carto.pal(pal1 = "blue.pal" , n1 = depth_values+1), transparency = 0.5)

min_values_vs = min(all_values_vs);
max_values_vs = max(all_values_vs);
seq_values_vs = seq(min(all_values_vs), max_values_vs, (max_values_vs - min(all_values_vs))/depth_values);
col_fun_edges_vs = colorRamp2(seq_values_vs, rev(coolwarm(depth_values+1)), transparency = 0.5)


ref_conn = conn_ongoing_day
all_region_labels = rownames(ref_conn)
all_region_labels[all_region_labels == "anterior hypothalamus"] = "Ant HT"
all_region_labels[all_region_labels == "arcuate hypothalamus"] = "Arc HT"
all_region_labels[all_region_labels == "dorsomedial hypothalamus"] = "DM HT"
all_region_labels[all_region_labels == "mammillary complex"] = "MaMa"
all_region_labels[all_region_labels == "paraventricular hypothalamus"] = "PV HT"
all_region_labels[all_region_labels == "posterior hypothalamus"] = "Post HT"
all_region_labels[all_region_labels == "ventromedial hypothalamus"] = "VM HT"
all_region_labels[all_region_labels == "ventromedial thalamus"] = "VM T"
all_region_labels[all_region_labels == "zona incerta"] = "Zona Inc"



png(file = paste0("/home/antonio/Documents/tim_brown/GLMCC/aggregate_circos_all.png"), units = "in",  width = 24, height = 24, res = 1200)
ref_conn = conn_ongoing_day

layout(matrix(1:9, 3, 3))

sectors = factor(all_region_labels, levels = all_region_labels)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 10))
set_track_gap(cm_h(0.1))
circos.track(sectors, ylim = c(0, 4), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.1))
circos.track(sectors, ylim = c(0, 4), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.1))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 3, col = "black", sector.index = si, track.index = 2)
}

values = array(ref_conn);
for(ii in 1:dim(ref_conn)[1]){
  for(jj in 1:dim(ref_conn)[2]){
    if (ref_conn[ii,jj] > )
    {  
      circos.link(all_region_labels[ii], c(3,5), all_region_labels[jj] , c(3,5), directional = 0,  col = col_fun_edges_pos(ref_conn[ii,jj]*3),
                  lwd = (ref_conn[ii,jj]*2)/(max_values_pos) )
    }
  }
}

col_fun_edges = colorRamp2(  seq(0, 1, 0.051), carto.pal(pal1 = "blue.pal" , n1 = depth_values+1), transparency = 0)
# leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
#                col_fun = col_fun_edges, title = "Strength", grid_width = unit(0.4, "in"),
#                labels_gp = gpar(col = "black", fontsize = 16),
#                title_gp = gpar(col = "black", fontsize = 20),
#                title_gap = unit(4, "mm"))
# draw(leg1, x = unit(1, "in"), y = unit(1, "in")) #just = c("right", "top")

ref_conn = conn_ongoing_night
circos.clear()
sectors = factor(all_region_labels, levels = all_region_labels)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 10))
set_track_gap(cm_h(0.1))
circos.track(sectors, ylim = c(0, 4), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.1))
circos.track(sectors, ylim = c(0, 4), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.1))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 3, col = "black", sector.index = si, track.index = 2)
}

values = array(ref_conn);
for(ii in 1:dim(ref_conn)[1]){
  for(jj in 1:dim(ref_conn)[2]){
    if (ref_conn[ii,jj] > )
    {  
      circos.link(all_region_labels[ii], c(3,5), all_region_labels[jj] , c(3,5), directional = 0,  col = col_fun_edges_pos(ref_conn[ii,jj]*3),
                  lwd = (ref_conn[ii,jj]*2)/(max_values_pos) )
    }
  }
}

# col_fun_edges = colorRamp2(  seq(0, 1, 0.051), carto.pal(pal1 = "blue.pal" , n1 = depth_values+1), transparency = 0)
# leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
#                col_fun = col_fun_edges, title = "Strength", grid_width = unit(0.4, "in"),
#                labels_gp = gpar(col = "black", fontsize = 16),
#                title_gp = gpar(col = "black", fontsize = 20),
#                title_gap = unit(4, "mm"))
# draw(leg1, x = unit(1, "in"), y = unit(1, "in")) #just = c("right", "top")


ref_conn = conn_ongoing_day_vs_night
circos.clear()
sectors = factor(all_region_labels, levels = all_region_labels)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 10))
set_track_gap(cm_h(0.1))
circos.track(sectors, ylim = c(0, 4), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.1))
circos.track(sectors, ylim = c(0, 4), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.1))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 3, col = "black", sector.index = si, track.index = 2)
}

values = array(ref_conn);

for(ii in 1:dim(ref_conn)[1]){
  for(jj in 1:dim(ref_conn)[2]){
    if (abs(ref_conn[ii,jj]) > median(abs(values)))
    {  
      circos.link(all_region_labels[ii], c(3,5), all_region_labels[jj] , c(3,5), directional = 0,  col = col_fun_edges_vs(ref_conn[ii,jj]*3),
                  lwd = abs((ref_conn[ii,jj]*2)/(max_values_pos)) )
    }
  }
}

# col_fun_edges = colorRamp2(seq(-1, 1, 0.1), rev(coolwarm(length(seq(-1, 1, 0.1)))), transparency = 0)
# leg1 <- Legend(at = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1), labels = c("-1", "-0.75", "-0.5", "-0.25", "0", "0.25", "0.5", "0.75", "1"), 
#                col_fun = col_fun_edges, title = "Strength", grid_width = unit(0.4, "in"),
#                labels_gp = gpar(col = "black", fontsize = 16),
#                title_gp = gpar(col = "black", fontsize = 20),
#                title_gap = unit(4, "mm"))
# draw(leg1, x = unit(1, "in"), y = unit(1.5, "in")) #just = c("right", "top")


ref_conn = conn_evoked_day
circos.clear()
sectors = factor(all_region_labels, levels = all_region_labels)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 10))
set_track_gap(cm_h(0.1))
circos.track(sectors, ylim = c(0, 4), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.1))
circos.track(sectors, ylim = c(0, 4), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.1))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 3, col = "black", sector.index = si, track.index = 2)
}

values = array(ref_conn);
for(ii in 1:dim(ref_conn)[1]){
  for(jj in 1:dim(ref_conn)[2]){
    if (ref_conn[ii,jj] > median(abs(values)))
    {  
      circos.link(all_region_labels[ii], c(3,5), all_region_labels[jj] , c(3,5), directional = 0,  col = col_fun_edges_pos(ref_conn[ii,jj]*3),
                  lwd = (ref_conn[ii,jj]*2)/(max_values_pos) )
    }
  }
}

# col_fun_edges = colorRamp2(  seq(0, 1, 0.051), carto.pal(pal1 = "blue.pal" , n1 = depth_values+1), transparency = 0)
# leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
#                col_fun = col_fun_edges, title = "Strength", grid_width = unit(0.4, "in"),
#                labels_gp = gpar(col = "black", fontsize = 16),
#                title_gp = gpar(col = "black", fontsize = 20),
#                title_gap = unit(4, "mm"))
# draw(leg1, x = unit(1, "in"), y = unit(1, "in")) #just = c("right", "top")

ref_conn = conn_evoked_night
circos.clear()
sectors = factor(all_region_labels, levels = all_region_labels)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 10))
set_track_gap(cm_h(0.1))
circos.track(sectors, ylim = c(0, 4), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.1))
circos.track(sectors, ylim = c(0, 4), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.1))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 3, col = "black", sector.index = si, track.index = 2)
}

values = array(ref_conn);
for(ii in 1:dim(ref_conn)[1]){
  for(jj in 1:dim(ref_conn)[2]){
    if (ref_conn[ii,jj] > median(abs(values)))
    {  
      circos.link(all_region_labels[ii], c(3,5), all_region_labels[jj] , c(3,5), directional = 0,  col = col_fun_edges_pos(ref_conn[ii,jj]*3),
                  lwd = (ref_conn[ii,jj]*2)/(max_values_pos) )
    }
  }
}

# col_fun_edges = colorRamp2(  seq(0, 1, 0.051), carto.pal(pal1 = "blue.pal" , n1 = depth_values+1), transparency = 0)
# leg1 <- Legend(at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1"), 
#                col_fun = col_fun_edges, title = "Strength", grid_width = unit(0.4, "in"),
#                labels_gp = gpar(col = "black", fontsize = 16),
#                title_gp = gpar(col = "black", fontsize = 20),
#                title_gap = unit(4, "mm"))
# draw(leg1, x = unit(1, "in"), y = unit(1, "in")) #just = c("right", "top")

ref_conn = conn_evoked_day_vs_night
circos.clear()
sectors = factor(all_region_labels, levels = all_region_labels)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 10))
set_track_gap(cm_h(0.1))
circos.track(sectors, ylim = c(0, 4), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.1))
circos.track(sectors, ylim = c(0, 4), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.1))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 3, col = "black", sector.index = si, track.index = 2)
}

values = array(ref_conn);
for(ii in 1:dim(ref_conn)[1]){
  for(jj in 1:dim(ref_conn)[2]){
    if (abs(ref_conn[ii,jj]) > median(abs(values)))
    {  
      circos.link(all_region_labels[ii], c(3,5), all_region_labels[jj] , c(3,5), directional = 0,  col = col_fun_edges_vs(ref_conn[ii,jj]*3),
                  lwd = abs((ref_conn[ii,jj]*2)/(max_values_vs)) )
    }
  }
}

# col_fun_edges = colorRamp2(seq(-1, 1, 0.1), rev(coolwarm(length(seq(-1, 1, 0.1)))), transparency = 0)
# leg1 <- Legend(at = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1), labels = c("-1", "-0.75", "-0.5", "-0.25", "0", "0.25", "0.5", "0.75", "1"), 
#                col_fun = col_fun_edges, title = "Strength", grid_width = unit(0.4, "in"),
#                labels_gp = gpar(col = "black", fontsize = 16),
#                title_gp = gpar(col = "black", fontsize = 20),
#                title_gap = unit(4, "mm"))
# draw(leg1, x = unit(1, "in"), y = unit(1.5, "in")) #just = c("right", "top")

ref_conn = conn_evoked_day - conn_ongoing_day
circos.clear()
sectors = factor(all_region_labels, levels = all_region_labels)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 10))
set_track_gap(cm_h(0.1))
circos.track(sectors, ylim = c(0, 4), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.1))
circos.track(sectors, ylim = c(0, 4), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.1))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 3, col = "black", sector.index = si, track.index = 2)
}

values = array(ref_conn);
for(ii in 1:dim(ref_conn)[1]){
  for(jj in 1:dim(ref_conn)[2]){
    if (abs(ref_conn[ii,jj]) > median(abs(values)))
    {  
      circos.link(all_region_labels[ii], c(3,5), all_region_labels[jj] , c(3,5), directional = 0,  col = col_fun_edges_vs(ref_conn[ii,jj]*3),
                  lwd = abs((ref_conn[ii,jj]*2)/(max_values_vs)) )
    }
  }
}

# col_fun_edges = colorRamp2(seq(-1, 1, 0.1), rev(coolwarm(length(seq(-1, 1, 0.1)))), transparency = 0)
# leg1 <- Legend(at = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1), labels = c("-1", "-0.75", "-0.5", "-0.25", "0", "0.25", "0.5", "0.75", "1"), 
#                col_fun = col_fun_edges, title = "Strength", grid_width = unit(0.4, "in"),
#                labels_gp = gpar(col = "black", fontsize = 16),
#                title_gp = gpar(col = "black", fontsize = 20),
#                title_gap = unit(4, "mm"))
# draw(leg1, x = unit(1, "in"), y = unit(1.5, "in")) #just = c("right", "top")

ref_conn = conn_evoked_night - conn_ongoing_night
circos.clear()
sectors = factor(all_region_labels, levels = all_region_labels)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(sectors, xlim = c(0, 10))
set_track_gap(cm_h(0.1))
circos.track(sectors, ylim = c(0, 4), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.1))
circos.track(sectors, ylim = c(0, 4), bg.border = "white", track.height = 0.2)
set_track_gap(cm_h(0.1))

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
  circos.text(mean(xlim), ylim[1], si, facing = "clockwise", adj = c(0, 0.1),
              niceFacing = TRUE, cex = 3, col = "black", sector.index = si, track.index = 2)
}

values = array(ref_conn);

for(ii in 1:dim(ref_conn)[1]){
  for(jj in 1:dim(ref_conn)[2]){
    if (abs(ref_conn[ii,jj]) > median(abs(values)))
    {  
      circos.link(all_region_labels[ii], c(3, 5), all_region_labels[jj] , c(3, 5), directional = 0,  col = col_fun_edges_vs(ref_conn[ii,jj]*3),
                  lwd = abs((ref_conn[ii,jj]*2)/(max_values_vs)) )
    }
  }
}

# col_fun_edges = colorRamp2(seq(-1, 1, 0.1), rev(coolwarm(length(seq(-1, 1, 0.1)))), transparency = 0)
# leg1 <- Legend(at = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1), labels = c("-1", "-0.75", "-0.5", "-0.25", "0", "0.25", "0.5", "0.75", "1"), 
#                col_fun = col_fun_edges, title = "Strength", grid_width = unit(0.4, "in"),
#                labels_gp = gpar(col = "black", fontsize = 16),
#                title_gp = gpar(col = "black", fontsize = 20),
#                title_gap = unit(4, "mm"))
# draw(leg1, x = unit(1, "in"), y = unit(1.5, "in")) #just = c("right", "top")

circos.clear()
dev.off()
