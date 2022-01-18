# This code loads all necessary libraries, seed column neurons, and some connectivity info.
# One can either load neurons and connectivity data from FAFB, which requires access permission, or from saved .Rdata files.
# To reproduce a plot, see "plot_centraleye.R" and "plot_DRA.R" for examples

# load library ----------------------------------------------------------
# usually you can install missing library using: install.packages("name of the library")

library(natverse)
library(tidyverse)

library(rjson)
library(ggpubr)
library(alphashape3d)
library(RColorBrewer) # palette
# library(webshot2) #snapshot3d, remotes::install_github("rstudio/chromote"),remotes::install_github("rstudio/webshot2")

# - for making the gallery
library(png)
library(grid)
library(gridExtra)

# clean everythign up.
rm(list=ls())

#close any open rgl windows 
while (rgl.cur() > 0) { rgl.close() }

# load some useful functions
source("R7R8_func.R")

#  palette ----------------------------------------------------------------

# pal_n <- c("#f781bf", "#1f78b4", "#984ea3", "#b2df8a", brewer.pal(9, "Paired")[c(2,6,9)] )

pal_syn <- brewer.pal(9,"Set1")[c(4,3)]
names(pal_syn) <- c("R7","R8")

# FAFB space axes --------------------------------------------------------------

axis_ori <- c(300e3, 150e3, 250e3)
axis_lat <- c(-10000, 0, 0) #green
axis_dor <- c(0, -10000, 0) #magenta
axis_post <- c(0, 0, 10000) #blue

# # neuron data from CATMAID server ------------------------------------------------------------
# 
# # - seed column R7 R8
# skid_pR7 <- c(10082582, 10538510)# pR7
# skid_yR7 <- c(10585940, 10653593) #yR7
# skid_pR8 <- c(10086691, 11466408) #pR8
# skid_yR8 <- c(10629254, 11468318) #yR8
# 
# skid_pyR7R8 <- c(skid_pR7, skid_yR7, skid_pR8, skid_yR8)
# 
# pR7 = read.neurons.catmaid(skid_pR7, .progress='text')
# yR7 = read.neurons.catmaid(skid_yR7, .progress='text')
# pR8 = read.neurons.catmaid(skid_pR8, .progress='text')
# yR8 = read.neurons.catmaid(skid_yR8, .progress='text')
# 
# pR7 <- fix_jump(pR7)
# yR7 <- fix_jump(yR7)
# pR8 <- fix_jump(pR8)
# yR8 <- fix_jump(yR8)
# 
# # - all p and y R7s
# pR7_all_info <- catmaid_query_by_annotation("pR7") #all pR7 annotation
# yR7_all_info <- catmaid_query_by_annotation("yR7")
# R7_p <-  read.neurons.catmaid(pR7_all_info$skid, .progress='text')
# R7_y <-  read.neurons.catmaid(yR7_all_info$skid, .progress='text')
# 
# nopen3d()
# plot3d(pR7, col = 'black')
# plot3d(yR7, col = 'brown')
# plot3d(pR8, col = 'grey')
# plot3d(yR8, col = 'gold')
# plot3d(R7_p, col='gray')
# plot3d(R7_y, col='yellow')
# 
# # seed DRAR7/R8
# skid_DRAR7 <- c(10191735, 10300949, 11728779)
# skid_DRAR8 <- c(10190508, 10300963, 11728827)
# skid_DRAR7R8 <- c(skid_DRAR7, skid_DRAR8)
# DRAR8 = read.neurons.catmaid(skid_DRAR8, .progress='text')
# DRAR7 = read.neurons.catmaid(skid_DRAR7, .progress='text')
# 
# # - Tm20 for LO mesh
# skid <- c(10423774, 11444392, 11450552, 11473668)
# Tm20 <- read.neurons.catmaid(skid, .progress='text')
# 
# # - complete DmDRA
# # Putative Dm-DRA1 10247371 TO	10247370
# # Putative Dm-DRA1 10440161 TO	10440160
# # Putative Dm-DRA1 12106450 TO	12106449
# # Putative Dm-DRA2 10411789 TO	10411788
# # Putative Dm-DRA2 10449077 TO	10449076
# # Putative Dm-DRA2 11710980 TO	11710979
# 
# # DmDRA1_skid <- c(10247370, 10440160, 12106449)
# # DmDRA1 <- read.neurons.catmaid(DmDRA1_skid, .progress='text')
# # DmDRA2_skid <- c(10411788, 10449076, 11710979)
# # DmDRA2 <- read.neurons.catmaid(DmDRA1_skid, .progress='text')
# #
# # - aMe12
# skid_aMe12 <- c(28841,164544,7038035)
# aMe12 <-  read.neurons.catmaid(skid= skid_aMe12, .progress='text' )
# aMe12 <- fix_jump(aMe12)
# 
# # # SAVE
# save(pR7, pR8, yR7, yR8, R7_p, R7_y, DRAR7, DRAR8, Tm20, aMe12, pR7_all_info, yR7_all_info, file = "data/neu_R7R8.RData")

# neuron data from saved data --------------------------------------------

# seed column R7 R8
skid_pR7 <- c(10082582, 10538510)# pR7
skid_pR8 <- c(10086691, 11466408) #pR8
skid_yR7 <- c(10585940, 10653593)
skid_yR8 <- c(10629254, 11468318) #yR8
skid_pyR7R8 <- c(skid_pR7, skid_yR7, skid_pR8, skid_yR8)

# seed col DRA
skid_DRAR7 <- c(10191735, 10300949, 11728779)
skid_DRAR8 <- c(10190508, 10300963, 11728827)
skid_DRAR7R8 <- c(skid_DRAR7, skid_DRAR8)

skid_aMe12 <- c(28841,164544,7038035)

# R7/8
load("data/neu_R7R8.RData")

# - Mi1
# DRA Mi1, access FAFB server
# anno_Mi1_DRA <- catmaid_query_by_annotation('DRA_column') # DRA col
# anno_Mi1_DRA_excl <- catmaid_query_by_annotation('non_PR_column') # DRA col excluded

load("data/Mi1_info.RData") #load locally

ind_Mi1_DRA <- which(anno_Mi1$skid %in% anno_Mi1_DRA$skid)
ind_Mi1_DRA_hcol <- c(133, 268, 42) #seed col
ind_Mi1_DRA <- ind_Mi1_DRA[!(ind_Mi1_DRA %in% ind_Mi1_DRA_hcol)] #exclude seed column
ind_Mi1_DRA_excl <- which(anno_Mi1$skid %in% anno_Mi1_DRA_excl$skid)

# save(anno_Mi1, anno_Mi1_DRA, anno_Mi1_DRA_excl, xyz_M5_avg, Mi1_xyz, Mi1_16_me, DRA_Mi1, file = "data/Mi1_info.RData")

# meshes ------------------------------------------------------------------

# # - Medulla layer 5
# M5_msh <- catmaid_get_volume("v14.ME_R_M5_Mi1")

# - medulla layer, old
# ME_msh <- as.mesh3d(FAFBNP.surf, 'ME_R')

# JFRC2010
load("data/JFRC2NP.surf.fafb.rda") # data from Greg Jefferis
ME_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="ME_R")
ME_L_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="ME_L")
LO_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="LO_R")
AOTU_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="AOTU_R")

LH_R_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="LH_R")
LOP_R_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="LOP_R")
LOP_L_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="LOP_L")
FLA_R_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="FLA_R")
IPS_L_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="IPS_L")
IPS_R_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="IPS_R")
SPS_R_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="SPS_R")
PLP_R_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="PLP_R")
PLP_L_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="PLP_L")
AME_R_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="AME_R")
AME_L_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="AME_L")
MB_CA_R_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="MB_CA_R")
MB_CA_L_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="MB_CA_L")
SCL_R_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="SCL_R")

FB_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="FB")
EB_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="EB")
MB_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="MB")
PB_msh <- as.mesh3d(JFRC2NP.surf.fafb, Regions="PB")

# # connectivity data from CATMAID ----------------------------------------------------------------------------------------------
# 
# # - central eye
# pyR7R8_conn <- list()
# pyR7R8_conn_in <- list()
# pyR7R8_conn_out <- list()
# for (j in 1:length(skid_pyR7R8)) {
#   pyR7R8_conn[[j]] <- catmaid_query_connected(skid_pyR7R8[j])
#   # group incoming / outgoing
#   pyR7R8_conn_in[[j]] <- pyR7R8_conn[[j]]$incoming %>%
#     group_by(partner) %>%
#     summarise(syn.count=sum(syn.count)) %>%
#     data.frame()
#   pyR7R8_conn_out[[j]] <- pyR7R8_conn[[j]]$outgoing %>%
#     group_by(partner) %>%
#     summarise(syn.count=sum(syn.count)) %>%
#     data.frame()
# }
# 
# # -- incoming
# R7R8_in <- pyR7R8_conn_in %>%
#   reduce(full_join, by="partner") %>%
#   replace(is.na(.), 0)
# colnames(R7R8_in) <- c("partner", "pR7a","pR7b","yR7a","yR7b","pR8a","pR8b","yR8a","yR8b")
# R7R8_in %<>% as_tibble() %>%
#   mutate(pR7 = rowSums(select(., starts_with("pR7")))) %>%
#   mutate(yR7 = rowSums(select(., starts_with("yR7")))) %>%
#   mutate(pR8 = rowSums(select(., starts_with("pR8")))) %>%
#   mutate(yR8 = rowSums(select(., starts_with("yR8"))))
# R7R8_in_ori <- R7R8_in
# 
# # - outgoing
# R7R8_out <- pyR7R8_conn_out %>%
#   reduce(full_join, by="partner") %>%
#   replace(is.na(.), 0)
# colnames(R7R8_out) <- c("partner", "pR7a","pR7b","yR7a","yR7b","pR8a","pR8b","yR8a","yR8b")
# R7R8_out %<>% as_tibble() %>%
#   mutate(pR7 = rowSums(select(., starts_with("pR7")))) %>%
#   mutate(yR7 = rowSums(select(., starts_with("yR7")))) %>%
#   mutate(pR8 = rowSums(select(., starts_with("pR8")))) %>%
#   mutate(yR8 = rowSums(select(., starts_with("yR8"))))
# R7R8_out_ori <- R7R8_out
# 
# 
# 
# # - DRA
# DRAR7R8_conn <- list()
# DRAR7R8_conn_in <- list()
# DRAR7R8_conn_out <- list()
# for (j in 1:length(skid_DRAR7R8)) {
#   DRAR7R8_conn[[j]] <- catmaid_query_connected(skid_DRAR7R8[j])
#   # group incoming / outgoing
#   DRAR7R8_conn_in[[j]] <- DRAR7R8_conn[[j]]$incoming %>%
#     group_by(partner) %>%
#     summarise(syn.count=sum(syn.count)) %>%
#     data.frame()
#   DRAR7R8_conn_out[[j]] <- DRAR7R8_conn[[j]]$outgoing %>%
#     group_by(partner) %>%
#     summarise(syn.count=sum(syn.count)) %>%
#     data.frame()
# }
# 
# # -- DRA incoming
# DRAR7R8_in <- DRAR7R8_conn_in %>%
#   reduce(full_join, by="partner") %>%
#   replace(is.na(.), 0)
# colnames(DRAR7R8_in) <- c("partner", "DRAR7a","DRAR7b","DRAR7c","DRAR8a","DRAR8b","DRAR8c")
# DRAR7R8_in %<>% as_tibble() %>%
#   mutate(DRAR7 = rowSums(select(., starts_with("DRAR7")))) %>%
#   mutate(DRAR8 = rowSums(select(., starts_with("DRAR8"))))
# DRAR7R8_in_ori <- DRAR7R8_in
# 
# # -- DRA outgoing
# DRAR7R8_out <- DRAR7R8_conn_out %>%
#   reduce(full_join, by="partner") %>%
#   replace(is.na(.), 0)
# colnames(DRAR7R8_out) <- c("partner", "DRAR7a","DRAR7b","DRAR7c","DRAR8a","DRAR8b","DRAR8c")
# DRAR7R8_out %<>% as_tibble() %>%
#   mutate(DRAR7 = rowSums(select(., starts_with("DRAR7")))) %>%
#   mutate(DRAR8 = rowSums(select(., starts_with("DRAR8"))))
# DRAR7R8_out_ori <- DRAR7R8_out
# 
# # SAVE
# save(R7R8_in_ori, R7R8_out_ori, DRAR7R8_in_ori, DRAR7R8_out_ori, file = "data/conn_R7R8.RData")

# connectivity data from saved ----------------------------------------------------------------------------------------------

load("data/conn_R7R8.RData")

# central eye, define a coord transformation for plotting --------------------
# use seed columns and surrounding columns to define a plane (via pca principal component analysis)
# transform (translation + rotation using pca) all relavent neurons/meshes 
# such that they are centered at the origin of the coordinate system and aligned with the xyz axes

Mi1_ind_p <- c(38, 40) # Mi1 index
Mi1_ind_y <- c(39, 41)
Mi1_ind_Nnb_p <- c(500, 657, 364, 742, 599, 209)
Mi1_ind_Nnb_y <- c(365, 35, 208, 210, 207, 605) 
Mi1_ind_pca <- c(Mi1_ind_p, Mi1_ind_y, Mi1_ind_Nnb_p, Mi1_ind_Nnb_y)

# Mi1_16_me <- nlapply(Mi1[Mi1_ind_pca], function(x) subset(x, pointsinside(x, ME_msh))) #me portion

# for Dm8 plot only
Mi1_ind_R7_p <- c(334, 457, 211, 205, 720, 248, 707, 342, 343) 
Mi1_ind_R7_y <- c(441,  10, 583, 718, 206, 484, 606, 667) 

# # DEBUG index of pale and yellow col
# nopen3d()
# plot3d(R7_p, lwd =2)
# points3d(xyz_M5_avg, size  = 10, col='gray')
# # identify3d(xyz_M5_avg)
# points3d(xyz_M5_avg[Mi1_ind_Nnb_p,], size = 13)
# points3d(xyz_M5_avg[Mi1_ind_Nnb_y,], size = 13, col='cyan')

# use pca to obtain transformation (translation + rotation matrix)
node_xyz <- xyzmatrix(Mi1_16_me)
me_pca <- prcomp(node_xyz)
if (me_pca$rotation[,1] %*% c(-0.84,  0.20, -0.49) < 0) {
  me_pca$rotation <- - me_pca$rotation
}
if (t(cross3D(me_pca$rotation[,1],me_pca$rotation[,2])) %*% me_pca$rotation[,3] < 0 ) {
  me_pca$rotation[,3] <- - me_pca$rotation[,3]
}

# transform neuron(s)
# Mi1_xform <-  xEucl_neu(Mi1, me_pca$rotation, me_pca$center) 
pR7_xform <-  xEucl_neu(pR7, me_pca$rotation, me_pca$center)
pR8_xform <-  xEucl_neu(pR8, me_pca$rotation, me_pca$center)
yR7_xform <-  xEucl_neu(yR7, me_pca$rotation, me_pca$center)
yR8_xform <-  xEucl_neu(yR8, me_pca$rotation, me_pca$center)

# axes for plotting
axis_ori_me <- c(-1.05e5, 30000, 35000)
axis_lat_me <- axis_lat %*% me_pca$rotation
axis_dor_me <- axis_dor %*% me_pca$rotation
axis_post_me <- axis_post %*% me_pca$rotation

# DRA, define a coord transformation -------------------------------

# use the medulla portion only
DRAR7_me <- nlapply(DRAR7, function(x) subset(x, pointsinside(x, ME_msh))) 
DRAR8_me <- nlapply(DRAR8, function(x) subset(x, pointsinside(x, ME_msh))) 
DRA_ref_com <- rbind(xyzmatrix(DRAR7_me[[1]]), xyzmatrix(DRAR8_me[[1]])) %>% colMeans() # use polar col

ii <- sweep(xyz_M5_avg, 2, DRA_ref_com)^2 %>% rowSums() %>% sqrt() %>% order() %>% head(.,7)
DRA_Mi1 <- nlapply(Mi1[ii], function(x) subset(x, pointsinside(x,ME_msh,rval='distance') > 0))

# pca
node_xyz <- xyzmatrix(DRA_Mi1)
DRA_me_pca <- prcomp(node_xyz)
if (DRA_me_pca$rotation[,1] %*% c(-0.84,  0.20, -0.49) < 0) {
  DRA_me_pca$rotation <- - DRA_me_pca$rotation
}
if (t(cross3D(DRA_me_pca$rotation[,1],DRA_me_pca$rotation[,2])) %*% DRA_me_pca$rotation[,3] < 0 ) {
  DRA_me_pca$rotation[,3] <- - DRA_me_pca$rotation[,3]
}

# Mi1_xform_DRA <-  xEucl_neu(Mi1, DRA_me_pca$rotation, DRA_me_pca$center)
# Mi1_me_xform_DRA <-  xEucl_neu(Mi1_me, DRA_me_pca$rotation, DRA_me_pca$center)

DRAR7_xform <- xEucl_neu(DRAR7, DRA_me_pca$rotation, DRA_me_pca$center)
DRAR8_xform <- xEucl_neu(DRAR8, DRA_me_pca$rotation, DRA_me_pca$center)

# DmDRA1_xform <- xEucl_neu(DmDRA1, DRA_me_pca$rotation, DRA_me_pca$center)
# DmDRA2_xform <- xEucl_neu(DmDRA2, DRA_me_pca$rotation, DRA_me_pca$center)

# # DEBUG, id pale and yellow col Mi1 index
# nopen3d()
# plot3d(DRAR7, lwd =2)
# points3d(xyz_M5_avg, size  = 10, col='gray')
# # identify3d(xyz_M5_avg)

# - transform meshes
DRA_ME_msh_xform <- ME_msh
DRA_ME_msh_xform$vb[1:3,] <- sweep(t(ME_msh$vb[1:3,]), 2, DRA_me_pca$center) %*% DRA_me_pca$rotation %>%  t()

DRA_AOTU_msh_xform <- AOTU_msh
DRA_AOTU_msh_xform$vb[1:3,] <- sweep(t(AOTU_msh$vb[1:3,]), 2, DRA_me_pca$center) %*% DRA_me_pca$rotation %>%  t()

DRA_ME_L_msh_xform <- ME_L_msh
DRA_ME_L_msh_xform$vb[1:3,] <- sweep(t(ME_L_msh$vb[1:3,]), 2, DRA_me_pca$center) %*% DRA_me_pca$rotation %>%  t()

DRA_PLP_R_msh_xform <- PLP_R_msh
DRA_PLP_R_msh_xform$vb[1:3,] <- sweep(t(PLP_R_msh$vb[1:3,]), 2, DRA_me_pca$center) %*% DRA_me_pca$rotation %>%  t()

# - DRA axes 
DRA_axis_ori <- c(-110000, -50000, -40000)
DRA_axis_lat <- -DRA_me_pca$rotation[1,] * 15000
DRA_axis_dor <- DRA_me_pca$rotation[2,] * -15000
DRA_axis_post <- DRA_me_pca$rotation[3,] * 15000

# - scale bar 
DRA_xAng <- -30 #rotation angle around x-axis, align with long axis of ME
DRA_xAngRot <- (DRA_xAng - 90)/180*pi
DRA_xRot <- matrix(c(1,0,0,
                     0, cos(DRA_xAngRot), sin(DRA_xAngRot),
                     0, -sin(DRA_xAngRot), cos(DRA_xAngRot)), ncol = 3, byrow = T)
DRA_xAngRot2 <- DRA_xAng/180*pi
DRA_xRot2 <- matrix(c(1,0,0,
                      0, cos(DRA_xAngRot2), sin(DRA_xAngRot2),
                      0, -sin(DRA_xAngRot2), cos(DRA_xAngRot2)), ncol = 3, byrow = T)

# top view
DRA_xAng_top <- 30
DRA_xAngRot_top <- (DRA_xAng_top)/180*pi
DRA_xRot_top <- matrix(c(1,0,0,
                     0, cos(DRA_xAngRot_top), sin(DRA_xAngRot_top),
                     0, -sin(DRA_xAngRot_top), cos(DRA_xAngRot_top)), ncol = 3, byrow = T)
DRA_yAng_top <- 90
DRA_yAngRot_top <- DRA_yAng_top/180*pi
DRA_yRot_top <- matrix(c(cos(DRA_yAngRot_top), 0, sin(DRA_yAngRot_top),
                     0, 1, 0,
                     -sin(DRA_yAngRot_top), 0, cos(DRA_yAngRot_top)), ncol = 3, byrow = T)
DRA_zAng_top <- 80 
DRA_zAngRot_top <- (DRA_zAng_top)/180*pi
DRA_zRot_top <- matrix(c(cos(DRA_zAngRot_top), sin(DRA_zAngRot_top), 0,
                         -sin(DRA_zAngRot_top), cos(DRA_zAngRot_top), 0,
                         0, 0, 1), ncol = 3, byrow = T)
DRA_sbar <- matrix(c(0,0,0, 0,1e4,0), ncol=3,byrow=T) %*%
  DRA_zRot_top %*%
  DRA_yRot_top %*%
  DRA_xRot_top

# side
DRA_sbar_side <- matrix(c(30000,5000,10000, 30000,5000,20000), ncol=3,byrow=T) %*% 
  t(matrix(c(1,0,0,
             0,0,1,
             0,-1,0), ncol=3, byrow=T)) 

#rotation angle around y-axis 
DRA_yAng <- 0
DRA_yAngRot <- -DRA_yAng/180*pi
DRA_yRot <- matrix(c(cos(DRA_yAngRot), 0, sin(DRA_yAngRot),
                     0,1,0,
                     -sin(DRA_yAngRot), 0, cos(DRA_yAngRot)), ncol = 3, byrow = T)

DRA_sbar_rot90 <- matrix(c(30000,5000,10000, 30000,5000,20000), ncol=3,byrow=T)

# a coord transformation for aMe12 plot -------------------------------------------------------

# for aMe12
M5_pca <- prcomp(xyz_M5_avg)
if (M5_pca$rotation[,1] %*% c(0,1,0) < 0) {
  M5_pca$rotation <- - M5_pca$rotation
}
if (t(cross3D(M5_pca$rotation[,1],M5_pca$rotation[,2])) %*% M5_pca$rotation[,3] < 0 ) {
  M5_pca$rotation[,3] <- - M5_pca$rotation[,3]
}

ME_msh_xform_M5 <- ME_msh
ME_msh_xform_M5$vb[1:3,] <- sweep(t(ME_msh_xform_M5$vb[1:3,]), 2, M5_pca$center) %*% M5_pca$rotation %>%  t()
AME_R_msh_xform_M5 <- AME_R_msh
AME_R_msh_xform_M5$vb[1:3,] <- sweep(t(AME_R_msh_xform_M5$vb[1:3,]), 2, M5_pca$center) %*% M5_pca$rotation %>%  t()

axis_ori_M5 <- c(1e5, -1e5, 0)
axis_lat_M5 <- axis_lat %*% M5_pca$rotation
axis_dor_M5 <- axis_dor %*% M5_pca$rotation
axis_post_M5 <- axis_post %*% M5_pca$rotation


aMe12_xform <-  xEucl_neu(aMe12, M5_pca$rotation, M5_pca$center)
aMe12_ME_xform <-  nlapply(aMe12, function(x) subset(x, pointsinside(x,ME_msh,rval='distance') > -80000)) %>%
  xEucl_neu(., M5_pca$rotation, M5_pca$center)


# exemplary cells -----------------------------------

eg_name <- c(
  "Putative Tm_OTu_other 11455157 CL",
  "Putative Mi15 11445262 MF",
  "Putative Dm2 11448823 HL",
  "Putative ML1 11472158 CL",
  "Putative Tm5a 11447184 MF",
  "Putative Tm5b 10356413 MF",
  "Putative ML2b 11458373 CL",
  "Putative Tm5c 11485095 MF",
  "Putative Tm5c 11450506 CL",
  "Putative Dm11 11450454 CL",
  "Putative Dm9 11452428 MF",
  "Putative Dm8 10208776 MF"
)

