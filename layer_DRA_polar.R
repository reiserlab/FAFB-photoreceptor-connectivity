# DRA use Mi1 to define a coord trans  -----------------------------------------

# - 2 col of DRA R7/8
DRAR7_me <- nlapply(DRAR7, function(x) subset(x, pointsinside(x, ME_msh))) #me portion
DRAR8_me <- nlapply(DRAR8, function(x) subset(x, pointsinside(x, ME_msh)))
# DRA_ref_com <- rbind(xyzmatrix(DRAR7_me), xyzmatrix(DRAR8_me)) %>% colMeans()
DRA_ref_com <- rbind(xyzmatrix(DRAR7_me[[1]]), xyzmatrix(DRAR8_me[[1]])) %>% colMeans() #change back to polar col

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

# DRA Mi1
anno_Mi1_DRA <- catmaid_query_by_annotation('DRA_column') # DRA col, light blue, or darker gray
anno_Mi1_DRA_excl <- catmaid_query_by_annotation('non_PR_column') # DRA col excluded

ind_Mi1_DRA <- which(anno_Mi1$skid %in% anno_Mi1_DRA$skid)
ind_Mi1_DRA_hcol <- c(133, 268, 42)
ind_Mi1_DRA <- ind_Mi1_DRA[!(ind_Mi1_DRA %in% ind_Mi1_DRA_hcol)] #exclude home column
ind_Mi1_DRA_excl <- which(anno_Mi1$skid %in% anno_Mi1_DRA_excl$skid)

# # DEBUG, id pale and yellow col Mi1 index
# nopen3d()
# plot3d(DRAR7, lwd =2)
# points3d(xyz_M5_avg, size  = 10, col='gray')
# # identify3d(xyz_M5_avg)

# - meshes
# - transfomr medulla mesh
DRA_ME_msh_xform <- ME_msh
DRA_ME_msh_xform$vb[1:3,] <- sweep(t(ME_msh$vb[1:3,]), 2, DRA_me_pca$center) %*% DRA_me_pca$rotation %>%  t()

DRA_AOTU_msh_xform <- AOTU_msh
DRA_AOTU_msh_xform$vb[1:3,] <- sweep(t(AOTU_msh$vb[1:3,]), 2, DRA_me_pca$center) %*% DRA_me_pca$rotation %>%  t()

DRA_ME_L_msh_xform <- ME_L_msh
DRA_ME_L_msh_xform$vb[1:3,] <- sweep(t(ME_L_msh$vb[1:3,]), 2, DRA_me_pca$center) %*% DRA_me_pca$rotation %>%  t()

DRA_PLP_R_msh_xform <- PLP_R_msh
DRA_PLP_R_msh_xform$vb[1:3,] <- sweep(t(PLP_R_msh$vb[1:3,]), 2, DRA_me_pca$center) %*% DRA_me_pca$rotation %>%  t()


# - axes
DRA_axis_ori <- c(-110000, -50000, -40000)
DRA_axis_lat <- -DRA_me_pca$rotation[1,] * 15000
DRA_axis_dor <- DRA_me_pca$rotation[2,] * -15000
DRA_axis_post <- DRA_me_pca$rotation[3,] * 15000

# - scale bar

DRA_xAng <- -30 #rotation angle around x-axis, align with long axis of ME
# DRA_xAng <- 20 # align col
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
DRA_yRot_top <- matrix(c(1,0,0,
                         0, cos(DRA_yAngRot_top), sin(DRA_yAngRot_top),
                         0, -sin(DRA_yAngRot_top), cos(DRA_yAngRot_top)), ncol = 3, byrow = T)
DRA_zAng_top <- 80
DRA_zAngRot_top <- (DRA_zAng_top)/180*pi
DRA_zRot_top <- matrix(c(1,0,0,
                         0, cos(DRA_zAngRot_top), sin(DRA_zAngRot_top),
                         0, -sin(DRA_zAngRot_top), cos(DRA_zAngRot_top)), ncol = 3, byrow = T)
DRA_sbar      <- matrix(c(30000,-5000,40000, 30000,5000,40000), ncol=3,byrow=T) %*% 
  t(matrix(c(1,0,0,
             0,0,1,
             0,-1,0), ncol=3, byrow=T))  %*% t(DRA_xRot2)

# # -- side, anterior view
# DRA_sbar_side      <- matrix(c(40000,-5000,20000, 40000,5000,20000), ncol=3,byrow=T) %*% 
#   t(matrix(c(1,0,0,
#              0,0,1,
#              0,-1,0), ncol=3, byrow=T)) %*% t(DRA_xRot2)
# 
# #rotation angle around y-axis to level layer 5
# DRA_yAng <- 0
# DRA_yAngRot <- -DRA_yAng/180*pi
# DRA_yRot <- matrix(c(cos(DRA_yAngRot), 0, sin(DRA_yAngRot),
#                      0,1,0,
#                      -sin(DRA_yAngRot), 0, cos(DRA_yAngRot)), ncol = 3, byrow = T)
# 
# # DRA_sbar_rot90 <- DRA_sbar %*% t(matrix(c(1,0,0,
# #                                           0,0,1,
# #                                           0,-1,0), ncol=3, byrow=T)) %*% t(DRA_xRot) %*% t(DRA_yRot)
# 
# DRA_sbar_rot90 <- matrix(c(30000,5000,10000, 30000,5000,20000), ncol=3,byrow=T)


# -- side, polar view
DRA_sbar_side = matrix(c(30000,5000,10000, 30000,5000,20000), ncol=3,byrow=T) %*% t(DRA_xRot2)

#rotation angle around y-axis to level layer 5
DRA_yAng <- -20
DRA_yAngRot <- -DRA_yAng/180*pi
DRA_yRot <- matrix(c(cos(DRA_yAngRot), 0, sin(DRA_yAngRot),
                     0,1,0,
                     -sin(DRA_yAngRot), 0, cos(DRA_yAngRot)), ncol = 3, byrow = T)

DRA_sbar_rot90 <- DRA_sbar %*% t(DRA_xRot) %*% t(DRA_yRot)

DRA_sbar_side_5 <- DRA_sbar_side %*% t(DRA_xRot) %*% t(DRA_yRot)


# col for top view --------------------------------------------------------

aMe12_ind_Mi1 <- list()
Mi1_ind_aMe12 <- list()
for (ii_n in 1:3) {
  tar <- aMe12[[ii_n]]
  ii_col = match(tar$tags$`columnar branch`, tar$d$PointNo)
  xyz_col = xyzmatrix(tar$d[ii_col,])
  aMe12_ind_Mi1[[ii_n]] <- apply(xyz_col, 1, function(x) order(rowSums(sweep(Mi1_M5_xyz,2,x)^2))[1] )
  Mi1_ind_aMe12[[ii_n]] <- apply(Mi1_M5_xyz, 1, function(x) order(rowSums(sweep(xyz_col,2,x)^2))[1] )
}

i1 <- sqrt(rowSums(sweep(xyz_M5_avg, 2, xyz_M5_avg[ind_Mi1_DRA_hcol[1],], '-')^2)) <60000
i2 <- sqrt(rowSums(sweep(xyz_M5_avg, 2, xyz_M5_avg[ind_Mi1_DRA_hcol[2],], '-')^2)) <30000
i3 <- sqrt(rowSums(sweep(xyz_M5_avg, 2, xyz_M5_avg[ind_Mi1_DRA_hcol[3],], '-')^2)) <30000
ind_Mi1_top <- which(i1|i2|i3)

ind_Mi1_top_grey <- ind_Mi1_top[!(ind_Mi1_top %in% c(ind_Mi1_DRA_excl, ind_Mi1_DRA))]
ind_Mi1_top_red <- ind_Mi1_top[ind_Mi1_top %in% c(ind_Mi1_DRA)]
ind_Mi1_top_y <- ind_Mi1_top_grey[!(ind_Mi1_top_grey %in% unlist(aMe12_ind_Mi1))]
ind_Mi1_top_p <- ind_Mi1_top_grey[(ind_Mi1_top_grey  %in% unlist(aMe12_ind_Mi1))]

# meshes ------------------------------------------------------------------


# - Mi1 M5 nodes
xyz_M5_avg_xform_DRA <- sweep(xyz_M5_avg, 2, DRA_me_pca$center) %*% DRA_me_pca$rotation

xyz_M5_avg_xform_DRA_yz <- xyz_M5_avg_xform_DRA
xyz_M5_avg_xform_DRA_yz[,1] <- 0

# Nnb <- 4
# ind_nb4 <- xyz_M5_avg_xform_DRA_yz^2 %>% rowSums() %>% order() %>% .[1:Nnb]

Nnb <- 16
ind_nb16 <- xyz_M5_avg_xform_DRA_yz^2 %>% rowSums() %>% order() %>% .[1:Nnb]

# Nnb <- 25
# ind_nb25 <- xyz_M5_avg_xform_DRA_yz^2 %>% rowSums() %>% order() %>% .[1:Nnb]


# nopen3d()
# points3d(xyz_M5_avg_xform_DRA_yz, col='grey', size = 20)
# points3d(xyz_M5_avg_xform_DRA, col='grey', size = 10)
# # points3d(ref_com4,col='red',size=10)
# points3d(xyz_M5_avg_xform_DRA_yz[ind_nb25,], size=25, col='blue')
#
# # # ind_16 <- identify3d(xyz_M5_avg_xform_yz)
# # Mi1_16 <- xEucl_neu(Mi1[ind_16], me_pca$rotation, me_pca$center)
#
# plot3d(DRAR7_xform, col='green', lwd=2)
# plot3d(DRAR8_xform, col='gold2', lwd=2)

# make new local mesh 
# ref_com_avg <- colMeans(ref_com4)
xyz <- t(DRA_ME_msh_xform$vb[1:3,])
dd <- sweep(xyz[,2:3], 2, colMeans(xyz_M5_avg_xform_DRA_yz[ind_nb16,])[-1]) %>% .^2 %>% rowSums() %>% sqrt()
xyz_msh <- xyz[dd < 1.2*max(dist(xyz_M5_avg_xform_DRA_yz[ind_nb16,])), ]
DRA_ME_msh_local_large <- ashape3d(xyz_msh, alpha = 60000) %>% as.mesh3d()

xyz_msh <- xyz[dd < 0.8*max(dist(xyz_M5_avg_xform_DRA_yz[ind_nb16,])), ]
DRA_ME_msh_local <- ashape3d(xyz_msh, alpha = 60000) %>% as.mesh3d()

# xyz_msh <- xyz[dd < max(dist(ref_com4)), ]
# ME_msh_local_2 <- ashape3d(xyz_msh, alpha = 50000) %>% as.mesh3d()

# - transfomr LO mesh
DRA_LO_msh_xform <- LO_msh
DRA_LO_msh_xform$vb[1:3,] <- sweep(t(LO_msh$vb[1:3,]), 2, DRA_me_pca$center) %*% DRA_me_pca$rotation %>%  t()


# ME layer middle DRAR7, polar view ----------------------------------------------------------------

# -- ME boundary
yz_msh <- t(DRA_ME_msh_local$vb)
# yz_dorsal <- data.frame(y = seq(20000, -20000, length.out = 20),
#                         z = seq(0, -0, length.out = 20))
# 
# yz_posterior <- data.frame(y = seq(0, 0, length.out = 16),
#                            z = seq(-20000, 20000, length.out = 16))

yz_dorsal <- data.frame(y = seq(5000, -45000, length.out = 20),
                        z = seq(0, -0, length.out = 20))

yz_posterior <- data.frame(y = seq(0, 0, length.out = 16),
                           z = seq(-45000, 5000, length.out = 16))

# choose yz
nopen3d()
points3d(yz_dorsal)
points3d(yz_posterior, col='red')
# points3d(yz_msh[,1:3])
shade3d(DRA_ME_msh_local, alpha=0.1, col='gray')
# shade3d(DRA_ME_msh_xform, alpha=0.1, col='gold')
axes3d(c('x','y','z')); title3d('','','x','y','z')

# points on mesh
xyz <- t(DRA_ME_msh_xform$vb[1:3,])
dd <- sweep(xyz[,2:3], 2, colMeans(xyz_M5_avg_xform_DRA_yz[ind_nb16,])[-1]) %>% .^2 %>% rowSums() %>% sqrt()
xyz_msh <- xyz[dd < max(dist(xyz_M5_avg_xform_DRA_yz[ind_nb16,])), ]
# top
xyz_top <- xyz_msh[xyz_msh[,1] > 17000, ]
x <- xyz_top[,1]; y <- xyz_top[,2]; z <- xyz_top[,3]
fitlm <- lm(x ~ poly(y, z, degree = 2, raw = T))
valfit <- predict(fitlm, yz_dorsal) #generate values from the fit
top_bd_dorsal <- cbind(valfit, yz_dorsal)
valfit <- predict(fitlm, yz_posterior) #generate values from the fit
top_bd_posterior <- cbind(valfit, yz_posterior)
# bottom
xyz_bot <- xyz_msh[xyz_msh[,1] < -24000, ]
x <- xyz_bot[,1]; y <- xyz_bot[,2]; z <- xyz_bot[,3]
fitlm <- lm(x ~ poly(y, z, degree = 2, raw = T))
valfit <- predict(fitlm, yz_dorsal) #generate values from the fit
bot_bd_dorsal <- cbind(valfit, yz_dorsal)
valfit <- predict(fitlm, yz_posterior) #generate values from the fit
bot_bd_posterior <- cbind(valfit, yz_posterior)

# # 2D mesh as outline
# DRA_ME_msh_local_bd <- rbind(top_bd_dorsal, apply(bot_bd_dorsal, MARGIN = 2, rev)) 
# DRA_ME_msh_local_bd <- rbind(DRA_ME_msh_local_bd[-1,], DRA_ME_msh_local_bd[1,]) %>%
#   cbind(DRA_ME_msh_local_bd, .) %>%
#   t() %>%
#   matrix(., ncol = 3, byrow = T)


# # -- layer Fischbach
# bp <- c(0, 0.55,0.4,0.5,0.2,0.35,0.4,0.6,0.3,0.7,0.4) %>% rev()
# bp_prob <- cumsum(bp)/sum(bp)
# x_qua <- apply(cbind(top_bd$valfit, bot_bd$valfit), MARGIN = 1,
#             function(x) quantile(x, probs = bp_prob) )
# layers <- matrix(ncol = 3, nrow = 0)
# for (j in 1:length(bp_prob)) {
#   cc <- cbind(x_qua[j,], yz)
#   c2 <- cbind(cc[-nrow(cc),], cc[-1,]) %>%
#     t() %>%
#     matrix(., ncol = 3, byrow = T)
#   layers <- rbind(layers, c2)
# }
# layers[,1] <- layers[,1] + 3000

# --layers from 7-column
# bp7c <- c(0, 5.5,14.5,20,23,28,36,42,48,56,61) #layer 1 to 10
bp7c <- c(0, 5.5, 14.5, 20, 23, 28, 33, 39, 45, 56, 61) #layer 1 to 10, modify 6-7-8-9 by Mi1, C2, Tm5 and Tm20
bp7c <- max(bp7c) - bp7c %>% rev() #layer 10 to 1
bp7c_prob <- bp7c/max(bp7c)

x_qua <- apply(cbind(top_bd_dorsal$valfit, bot_bd_dorsal$valfit), MARGIN = 1,
               function(x) quantile(x, probs = bp7c_prob) )
DRA_layers_ME_dorsal <- matrix(ncol = 3, nrow = 0)
for (j in 1:length(bp7c_prob)) {
  cc <- cbind(x_qua[j,], yz_dorsal)
  c2 <- cbind(cc[-nrow(cc),], cc[-1,]) %>%
    t() %>%
    matrix(., ncol = 3, byrow = T)
  DRA_layers_ME_dorsal <- rbind(DRA_layers_ME_dorsal, c2)
}
# DRA_layers_ME_dorsal[,1] <- DRA_layers_ME_dorsal[,1] * 0.97 + 1100

x_qua <- apply(cbind(top_bd_posterior$valfit, bot_bd_posterior$valfit), MARGIN = 1,
               function(x) quantile(x, probs = bp7c_prob) )
DRA_layers_ME_posterior <- matrix(ncol = 3, nrow = 0)
for (j in 1:length(bp7c_prob)) {
  cc <- cbind(x_qua[j,], yz_posterior)
  c2 <- cbind(cc[-nrow(cc),], cc[-1,]) %>%
    t() %>%
    matrix(., ncol = 3, byrow = T)
  DRA_layers_ME_posterior <- rbind(DRA_layers_ME_posterior, c2)
}
bot2 <- DRA_layers_ME_posterior[(nrow(yz_posterior)*2-3):(nrow(yz_posterior)*2-2),]
top2 <- DRA_layers_ME_posterior[(nrow(DRA_layers_ME_posterior)-1):(nrow(DRA_layers_ME_posterior)), ]
DRA_layers_ME_posterior <- DRA_layers_ME_posterior[-c(
  (nrow(yz_posterior)*2-3):(nrow(yz_posterior)*2-2),
  (nrow(DRA_layers_ME_posterior)-1):(nrow(DRA_layers_ME_posterior))), ]

cc <- cbind(x_qua[2:10,nrow(yz_posterior)], yz_posterior[nrow(yz_posterior),])
c2 <- cbind(cc[-nrow(cc),], cc[-1,]) %>%
  t() %>%
  matrix(., ncol = 3, byrow = T)
DRA_layers_ME_posterior <- rbind(DRA_layers_ME_posterior, c2)

# add round corner
pt_corner <- rbind(bot2, 
                   as.numeric(c(x_qua[2,nrow(yz_posterior)], yz_posterior[nrow(yz_posterior),])))
ptArc <- round_corner(pt_corner)
cc <- rbind(pt_corner[3,], ptArc, pt_corner[1,])
c2 <- cbind(cc[-nrow(cc),], cc[-1,]) %>%
  t() %>%
  matrix(., ncol = 3, byrow = T)
DRA_layers_ME_posterior <- rbind(DRA_layers_ME_posterior, c2)

pt_corner <- rbind(top2, 
                   as.numeric(c(x_qua[10,nrow(yz_posterior)], yz_posterior[nrow(yz_posterior),])))
ptArc <- round_corner(pt_corner)
cc <- rbind(pt_corner[3,], ptArc, pt_corner[1,])
c2 <- cbind(cc[-nrow(cc),], cc[-1,]) %>%
  t() %>%
  matrix(., ncol = 3, byrow = T)
DRA_layers_ME_posterior <- rbind(DRA_layers_ME_posterior, c2)

DRA_layers_ME_posterior_rot90 <- DRA_layers_ME_posterior %*% t(DRA_xRot2)


# bd_Mx <- x_qua[,ncol(x_qua)] %>% rev()
# bd_Mx <- bd_Mx * 0.97 + 1100
# bd_Mx <- bd_Mx[-length(bd_Mx)] + diff(bd_Mx)/2
# DRA_layers_ME_anno <- cbind(bd_Mx, yz_dorsal[nrow(yz_dorsal),])



# # ME layer side DRAR7, 2 - posterior  ---------------------------------------------------
# 
# # Nnb <- 9
# # ind_nb <- sweep(xyz_M5_avg_xform_DRA, 2,
# #             xyz_M5_avg_xform_DRA[ind_Mi1_DRA_hcol[2],])^2 %>%
# #   rowSums() %>%
# #   order() %>%
# #   .[1:Nnb]
# # ind_nb <- c(704, 644, 419, 217, 323, 537, 642,221)
# 
# ind_nb <- c(521, 325, 656, 770, 268, 518, 516, 517, 269, 774, 322, 494, 425)
# 
# # nopen3d()
# # points3d(xyz_M5_avg_xform_DRA_yz, col='grey', size = 20)
# # # points3d(xyz_M5_avg_xform_DRA, col='grey', size = 10)
# # # points3d(ref_com4,col='red',size=10)
# # points3d(xyz_M5_avg_xform_DRA_yz[ind_nb,], size=25, col='blue')
# # plot3d(DRAR7_xform, col='green', lwd=2)
# # plot3d(DRAR8_xform, col='gold2', lwd=2)
# 
# # -- make local mesh 
# xyz <- t(DRA_ME_msh_xform$vb[1:3,])
# dd <- sweep(xyz[,2:3], 2, colMeans(xyz_M5_avg_xform_DRA_yz[ind_nb,])[-1]) %>% .^2 %>% rowSums() %>% sqrt()
# xyz_msh <- xyz[dd < 1.2*max(dist(xyz_M5_avg_xform_DRA_yz[ind_nb,])), ]
# DRA_ME_msh_local_2 <- ashape3d(xyz_msh, alpha = 60000) %>% as.mesh3d()
# 
# # -- ME boundary
# yz_msh <- t(DRA_ME_msh_local_2$vb)
# yz_dorsal <- data.frame(y = seq(-20000, -50000, length.out = 20),
#                         z = seq(-0000, -60000, length.out = 20))
# yz_posterior <- data.frame(y = seq(-10000, -45000, length.out = 16),
#                            z = seq(-40000, -25000, length.out = 16))
# 
# # choose yz
# nopen3d()
# plot3d(DRAR7_xform[[2]], lwd=2)
# # points3d(yz_msh[,1:3])
# shade3d(DRA_ME_msh_local_2, alpha=0.1, col='gray')
# points3d(yz_dorsal)
# points3d(yz_posterior, col='red')
# # shade3d(DRA_ME_msh_xform, alpha=0.1, col='gold')
# axes3d(c('x','y','z')); title3d('','','x','y','z')
# 
# points3d(xyz_bot)
# # points3d(top_bd_dorsal_2)
# # points3d(bot_bd_dorsal_2)
# 
# # points on mesh
# xyz <- t(DRA_ME_msh_xform$vb[1:3,])
# dd <- sweep(xyz[,2:3], 2, colMeans(xyz_M5_avg_xform_DRA_yz[ind_nb,])[-1]) %>% .^2 %>% rowSums() %>% sqrt()
# xyz_msh <- xyz[dd < max(dist(xyz_M5_avg_xform_DRA_yz[ind_nb,])), ]
# # top
# xyz_top <- xyz_msh[xyz_msh[,1] > 30000, ]
# x <- xyz_top[,1]; y <- xyz_top[,2]; z <- xyz_top[,3]
# fitlm <- lm(x ~ poly(y, z, degree = 2, raw = T))
# valfit <- predict(fitlm, yz_dorsal) #generate values from the fit
# top_bd_dorsal_2 <- cbind(valfit, yz_dorsal)
# valfit <- predict(fitlm, yz_posterior) #generate values from the fit
# top_bd_posterior_2 <- cbind(valfit, yz_posterior)
# # bottom
# xyz_bot <- xyz_msh[xyz_msh[,1] < -18000 & xyz_msh[,1] > -28000, ]
# x <- xyz_bot[,1]; y <- xyz_bot[,2]; z <- xyz_bot[,3]
# fitlm <- lm(x ~ poly(y, z, degree = 2, raw = T))
# valfit <- predict(fitlm, yz_dorsal) #generate values from the fit
# bot_bd_dorsal_2 <- cbind(valfit, yz_dorsal)
# valfit <- predict(fitlm, yz_posterior) #generate values from the fit
# bot_bd_posterior_2 <- cbind(valfit, yz_posterior)
# 
# x_qua <- apply(cbind(top_bd_dorsal_2$valfit, bot_bd_dorsal_2$valfit), MARGIN = 1,
#                function(x) quantile(x, probs = bp7c_prob) )
# DRA_layers_ME_dorsal_2 <- matrix(ncol = 3, nrow = 0)
# for (j in 1:length(bp7c_prob)) {
#   cc <- cbind(x_qua[j,], yz_dorsal)
#   c2 <- cbind(cc[-nrow(cc),], cc[-1,]) %>%
#     t() %>%
#     matrix(., ncol = 3, byrow = T)
#   DRA_layers_ME_dorsal_2 <- rbind(DRA_layers_ME_dorsal_2, c2)
# }
# # manual translation
# DRA_layers_ME_dorsal_2[,1] <- DRA_layers_ME_dorsal_2[,1] * 1.3 - 12000
# 
# x_qua <- apply(cbind(top_bd_posterior_2$valfit, bot_bd_posterior_2$valfit), MARGIN = 1,
#                function(x) quantile(x, probs = bp7c_prob) )
# DRA_layers_ME_posterior_2 <- matrix(ncol = 3, nrow = 0)
# for (j in 1:length(bp7c_prob)) {
#   cc <- cbind(x_qua[j,], yz_posterior)
#   c2 <- cbind(cc[-nrow(cc),], cc[-1,]) %>%
#     t() %>%
#     matrix(., ncol = 3, byrow = T)
#   DRA_layers_ME_posterior_2 <- rbind(DRA_layers_ME_posterior_2, c2)
# }
# bot2 <- DRA_layers_ME_posterior_2[(nrow(yz_posterior)*2-3):(nrow(yz_posterior)*2-2),]
# top2 <- DRA_layers_ME_posterior_2[(nrow(DRA_layers_ME_posterior_2)-1):(nrow(DRA_layers_ME_posterior_2)), ]
# DRA_layers_ME_posterior_2 <- DRA_layers_ME_posterior_2[-c(
#   (nrow(yz_posterior)*2-3):(nrow(yz_posterior)*2-2),
#   (nrow(DRA_layers_ME_posterior_2)-1):(nrow(DRA_layers_ME_posterior_2))), ]
# 
# cc <- cbind(x_qua[2:10,nrow(yz_posterior)], yz_posterior[nrow(yz_posterior),])
# c2 <- cbind(cc[-nrow(cc),], cc[-1,]) %>%
#   t() %>%
#   matrix(., ncol = 3, byrow = T)
# DRA_layers_ME_posterior_2 <- rbind(DRA_layers_ME_posterior_2, c2)
# 
# # add round corner
# pt_corner <- rbind(bot2, 
#                    as.numeric(c(x_qua[2,nrow(yz_posterior)], yz_posterior[nrow(yz_posterior),])))
# ptArc <- round_corner(pt_corner)
# cc <- rbind(pt_corner[3,], ptArc, pt_corner[1,])
# c2 <- cbind(cc[-nrow(cc),], cc[-1,]) %>%
#   t() %>%
#   matrix(., ncol = 3, byrow = T)
# DRA_layers_ME_posterior_2 <- rbind(DRA_layers_ME_posterior_2, c2)
# 
# pt_corner <- rbind(top2, 
#                    as.numeric(c(x_qua[10,nrow(yz_posterior)], yz_posterior[nrow(yz_posterior),])))
# ptArc <- round_corner(pt_corner)
# cc <- rbind(pt_corner[3,], ptArc, pt_corner[1,])
# c2 <- cbind(cc[-nrow(cc),], cc[-1,]) %>%
#   t() %>%
#   matrix(., ncol = 3, byrow = T)
# DRA_layers_ME_posterior_2 <- rbind(DRA_layers_ME_posterior_2, c2)
# 
# # manual adjust
# DRA_layers_ME_posterior_2[,1] <- DRA_layers_ME_posterior_2[,1]*1.1 - 9000
# 
# # DRA_layers_ME_posterior_2_rot90 <- DRA_layers_ME_posterior_2 %*% t(DRA_xRot2)
# 
# 
# 
# 
# # # for DRAR7 - 3
# # # -- ME boundary
# # yz_msh <- t(DRA_ME_msh_local_2$vb)
# # yz_dorsal <- data.frame(y = seq(-20000, 30000, length.out = 20),
# #                         z = seq(100000, 50000, length.out = 20))
# # 
# # # bottom
# # xyz_bot <- xyz_msh[xyz_msh[,1] < -45000 & xyz_msh[,1] > -70000, ]
# # # top
# # xyz_top <- xyz_msh[xyz_msh[,1] > 21000, ]

# # ME layer side DRAR7, 3 - anterior ---------------------------------------
# 
# Nnb <- 16
# ind_nb <- sweep(xyz_M5_avg_xform_DRA, 2,
#             xyz_M5_avg_xform_DRA[ind_Mi1_DRA_hcol[3],])^2 %>%
#   rowSums() %>%
#   order() %>%
#   .[1:Nnb]
# 
# nopen3d()
# points3d(xyz_M5_avg_xform_DRA_yz, col='grey', size = 20)
# # points3d(xyz_M5_avg_xform_DRA, col='grey', size = 10)
# # points3d(ref_com4,col='red',size=10)
# points3d(xyz_M5_avg_xform_DRA_yz[ind_nb,], size=25, col='blue')
# plot3d(DRAR7_xform, col='green', lwd=2)
# plot3d(DRAR8_xform, col='gold2', lwd=2)
# 
# # -- make local mesh
# xyz <- t(DRA_ME_msh_xform$vb[1:3,])
# dd <- sweep(xyz[,2:3], 2, colMeans(xyz_M5_avg_xform_DRA_yz[ind_nb,])[-1]) %>% .^2 %>% rowSums() %>% sqrt()
# xyz_msh <- xyz[dd < 1.2*max(dist(xyz_M5_avg_xform_DRA_yz[ind_nb,])), ]
# DRA_ME_msh_local_3 <- ashape3d(xyz_msh, alpha = 60000) %>% as.mesh3d()
# 
# # -- ME boundary
# # yz_msh <- t(DRA_ME_msh_local_3$vb)
# # yz_dorsal <- data.frame(y = seq(80000, 35000, length.out = 20),
# #                         z = seq(-40000, -0000, length.out = 20))
# # 
# # yz_posterior <- data.frame(y = seq(30000, 70000, length.out = 16),
# #                            z = seq(-40000, -8000, length.out = 16))
# 
# yz_msh <- t(DRA_ME_msh_local_3$vb)
# yz_dorsal <- data.frame(y = seq(-35000, 30000, length.out = 16),
#                         z = seq(-30000, 10000, length.out = 16))
# 
# yz_posterior  <- data.frame(y = seq( -15000,7500, length.out = 20),
#                             z = seq( 21300,-12000, length.out = 20))
# 
# # choose yz
# nopen3d()
# shade3d(DRA_ME_msh_local_3, alpha=0.1, col='gray')
# plot3d(DRAR7_xform[[3]], lwd=2)
# points3d(yz_dorsal)
# points3d(yz_posterior, col='red')
# # points3d(yz_msh[,1:3])
# # shade3d(DRA_ME_msh_xform, alpha=0.1, col='gold')
# axes3d(c('x','y','z')); title3d('','','x','y','z')
# # points3d(xyz_bot)
# # points3d(top_bd_dorsal_3)
# # points3d(bot_bd_dorsal_3)
# 
# # points on mesh
# xyz <- t(DRA_ME_msh_xform$vb[1:3,])
# dd <- sweep(xyz[,2:3], 2, colMeans(xyz_M5_avg_xform_DRA_yz[ind_nb,])[-1]) %>% .^2 %>% rowSums() %>% sqrt()
# xyz_msh <- xyz[dd < max(dist(xyz_M5_avg_xform_DRA_yz[ind_nb,])), ]
# # top
# xyz_top <- xyz_msh[xyz_msh[,1] > 0, ]
# x <- xyz_top[,1]; y <- xyz_top[,2]; z <- xyz_top[,3]
# fitlm <- lm(x ~ poly(y, z, degree = 2, raw = T))
# valfit <- predict(fitlm, yz_dorsal) #generate values from the fit
# top_bd_dorsal_3 <- cbind(valfit, yz_dorsal)
# valfit <- predict(fitlm, yz_posterior) #generate values from the fit
# top_bd_posterior_3 <- cbind(valfit, yz_posterior)
# # bottom
# # xyz_bot <- xyz_msh[xyz_msh[,1] < -24000, ]
# xyz_bot <- xyz_msh[xyz_msh[,1] < -24000 & xyz_msh[,2] < 68000, ]
# x <- xyz_bot[,1]; y <- xyz_bot[,2]; z <- xyz_bot[,3]
# fitlm <- lm(x ~ poly(y, z, degree = 2, raw = T))
# valfit <- predict(fitlm, yz_dorsal) #generate values from the fit
# bot_bd_dorsal_3 <- cbind(valfit, yz_dorsal)
# valfit <- predict(fitlm, yz_posterior) #generate values from the fit
# bot_bd_posterior_3 <- cbind(valfit, yz_posterior)
# 
# x_qua <- apply(cbind(top_bd_dorsal_3$valfit, bot_bd_dorsal_3$valfit), MARGIN = 1,
#                function(x) quantile(x, probs = bp7c_prob) )
# DRA_layers_ME_dorsal_3 <- matrix(ncol = 3, nrow = 0)
# for (j in 1:length(bp7c_prob)) {
#   cc <- cbind(x_qua[j,], yz_dorsal)
#   c2 <- cbind(cc[-nrow(cc),], cc[-1,]) %>%
#     t() %>%
#     matrix(., ncol = 3, byrow = T)
#   DRA_layers_ME_dorsal_3 <- rbind(DRA_layers_ME_dorsal_3, c2)
# }
# 
# # manual translation
# DRA_layers_ME_dorsal_3[,1] <- DRA_layers_ME_dorsal_3[,1] - 4000
# 
# 
# x_qua <- apply(cbind(top_bd_posterior_3$valfit, bot_bd_posterior_3$valfit), MARGIN = 1,
#                function(x) quantile(x, probs = bp7c_prob) )
# DRA_layers_ME_posterior_3 <- matrix(ncol = 3, nrow = 0)
# for (j in 1:length(bp7c_prob)) {
#   cc <- cbind(x_qua[j,], yz_posterior)
#   c2 <- cbind(cc[-nrow(cc),], cc[-1,]) %>%
#     t() %>%
#     matrix(., ncol = 3, byrow = T)
#   DRA_layers_ME_posterior_3 <- rbind(DRA_layers_ME_posterior_3, c2)
# }
# bot2 <- DRA_layers_ME_posterior_3[(nrow(yz_posterior)*2-3):(nrow(yz_posterior)*2-2),]
# top2 <- DRA_layers_ME_posterior_3[(nrow(DRA_layers_ME_posterior_3)-1):(nrow(DRA_layers_ME_posterior_3)), ]
# DRA_layers_ME_posterior_3 <- DRA_layers_ME_posterior_3[-c(
#   (nrow(yz_posterior)*2-3):(nrow(yz_posterior)*2-2),
#   (nrow(DRA_layers_ME_posterior_3)-1):(nrow(DRA_layers_ME_posterior_3))), ]
# 
# cc <- cbind(x_qua[2:10,nrow(yz_posterior)], yz_posterior[nrow(yz_posterior),])
# c2 <- cbind(cc[-nrow(cc),], cc[-1,]) %>%
#   t() %>%
#   matrix(., ncol = 3, byrow = T)
# DRA_layers_ME_posterior_3 <- rbind(DRA_layers_ME_posterior_3, c2)
# 
# # add round corner
# pt_corner <- rbind(bot2,
#                    as.numeric(c(x_qua[2,nrow(yz_posterior)], yz_posterior[nrow(yz_posterior),])))
# ptArc <- round_corner(pt_corner)
# cc <- rbind(pt_corner[3,], ptArc, pt_corner[1,])
# c2 <- cbind(cc[-nrow(cc),], cc[-1,]) %>%
#   t() %>%
#   matrix(., ncol = 3, byrow = T)
# DRA_layers_ME_posterior_3 <- rbind(DRA_layers_ME_posterior_3, c2)
# 
# pt_corner <- rbind(top2,
#                    as.numeric(c(x_qua[10,nrow(yz_posterior)], yz_posterior[nrow(yz_posterior),])))
# ptArc <- round_corner(pt_corner)
# cc <- rbind(pt_corner[3,], ptArc, pt_corner[1,])
# c2 <- cbind(cc[-nrow(cc),], cc[-1,]) %>%
#   t() %>%
#   matrix(., ncol = 3, byrow = T)
# DRA_layers_ME_posterior_3 <- rbind(DRA_layers_ME_posterior_3, c2)
# 
# # manual adjust
# DRA_layers_ME_posterior_3[,1] <- DRA_layers_ME_posterior_3[,1] - 5500
# 
# # DRA_layers_ME_posterior_3_rot90 <- DRA_layers_ME_posterior_3 %*% t(DRA_xRot2)



# # SAVE --------------------------------------------------------------------
# 
# save(DRA_layers_ME_dorsal, DRA_layers_ME_dorsal_2, DRA_layers_ME_dorsal_3,
#      DRA_layers_ME_posterior, DRA_layers_ME_posterior_2, DRA_layers_ME_posterior_3,
#      DRA_layers_ME_posterior_rot90, file = "data/DRA_layers.RData")

# # Dm-DRA1 test view ---------------------------------------------------------------------
# 
# type <- "Dm-DRA1"
# tb <- read.csv(paste("table_by_type_out_DRA/", "DRAR7R8_outgoing_", type, ".csv", sep = ''))
# 
# skid <- tb$skid %>% na.omit()
# neu <-  read.neurons.catmaid(skid, .progress='text' )
# neu_xform <- xEucl_neu(neu, DRA_me_pca$rotation, DRA_me_pca$center)
# 
# ii <- 4 #4, 2, 1 for 3 seed col positions
# 
# 
# # - top
# nopen3d()
# par3d('windowRect' = c(100,100,1300,1300))
# plot3d(neu_xform[[ii]], col='gray40', lwd=2, soma = T, WithNodes = F, lit=F)
# pch3d(xyz_M5_avg_xform_DRA_yz[ind_Mi1_DRA_hcol,,drop=F], radius=5000,col="#b7252a",pch=16,alpha=0.7)
# pch3d(xyz_M5_avg_xform_DRA_yz[ind_Mi1_top_grey,], radius=5000,col="grey",pch=16,alpha=0.2)
# pch3d(xyz_M5_avg_xform_DRA_yz[ind_Mi1_top_red,], radius=5000,col="#b7252a",pch=16,alpha=0.3)
# segments3d(sweep(DRA_sbar,2,c(0,0,-20000)), lwd=2)
# text3d(colMeans(DRA_sbar)+c(0,0,23000), texts = "10 um", cex = 1.5)
# rgl.viewpoint(fov=0,zoom=1, userMatrix= t(translationMatrix(1e4,0,0)) %*%
#                 rotationMatrix(-80/180*pi,0,0,1) %*%
#                 rotationMatrix(-90/180*pi,0,1,0) %*%
#                 rotationMatrix(DRA_xAng/180*pi,1,0,0) )
# 
# # add syn from R7
# conn_in_R7 <- catmaid_get_connectors_between(pre_skids = skid_R7s, post_skids = skid[ii])
# if (!is.null(conn_in_R7)) {
#   conn_in_R7 <- conn_in_R7[, c("post_node_x", "post_node_y", "post_node_z")]
#   conn_in_R7_xform <- sweep(as.matrix(conn_in_R7), 2, DRA_me_pca$center) %*% DRA_me_pca$rotation
#   pch3d(conn_in_R7_xform, pch=16, radius=1000, alpha=0.7, col=pal_syn['R7'])
# }
# 
# # # save
# # rgl.snapshot(filename = paste("F5_", type, "_", skid[ii], "_top.png", sep = ''))
# 
# 
# # # - side
# # nopen3d()
# # par3d('windowRect' = c(100,100,1300,1300))
# # plot3d(neu_xform[[ii]], col='gray40', lwd=2, soma = T, WithNodes = F, lit=F)
# # segments3d(DRA_sbar_side, lwd=2)
# # plot3d(DRAR7_xform[[1]], col=pal_syn["R7"], lwd=2, soma = T, WithNodes = F, lit=F)
# # segments3d(DRA_layers_ME_dorsal, lwd=1)
# # rgl.viewpoint(fov=0,zoom=0.75, userMatrix= rotationMatrix(90/180*pi,0,0,1) %*%
# #                 rotationMatrix(0/180*pi,1,0,0) %*%
# #                 rotationMatrix(DRA_yAng/180*pi,0,1,0) ) # add 15 deg for the tilt
# # conn_in_R7 <- catmaid_get_connectors_between(pre_skids = skid_R7s, post_skids = skid[ii])
# # if (!is.null(conn_in_R7)) {
# #   conn_in_R7 <- conn_in_R7[, c("post_node_x", "post_node_y", "post_node_z")]
# #   conn_in_R7_xform <- sweep(as.matrix(conn_in_R7), 2, DRA_me_pca$center) %*% DRA_me_pca$rotation
# #   pch3d(conn_in_R7_xform, pch=16, radius=1000, alpha=0.7, col=pal_syn['R7'])
# # }
# # 
# # # save
# # rgl.snapshot(filename = paste("F5_", type, "_", skid[ii], "_side_dorsal.png", sep = ''))
# # 
# 
# 
# # rot 90
# ii_hc <- sweep(xyz_M5_avg_xform_DRA[ind_Mi1_DRA_hcol,], 2, colMeans(xyzmatrix(neu_xform[[ii]]$d)), '-') %>%
#   .^2 %>% rowSums() %>% order()
# ii_hc[1]
# 
# nopen3d()
# par3d('windowRect' = c(100,100,1300,1300))
# plot3d(neu_xform[[ii]], col='gray40', lwd=2, soma = T, WithNodes = F, lit=F)
# if (ii_hc[1] == 1) {
#   segments3d(DRA_layers_ME_posterior_rot90, lwd=1)
#   rgl.viewpoint(fov=0,zoom=1, userMatrix= rotationMatrix(90/180*pi,0,0,1) %*%
#                   rotationMatrix((DRA_xAng-90)/180*pi,1,0,0) %*%
#                   rotationMatrix(DRA_yAng/180*pi,0,1,0))
# } else if (ii_hc[1] == 2) {
#   segments3d(DRA_layers_ME_posterior_2, lwd=1)
#   rgl.viewpoint(fov=0,zoom=1,
#                 userMatrix= matrix(c(-0.084,  0.870, -0.485,  -1044.4827,
#                                      0.996,  0.077, -0.034,    392.4445,
#                                      0.0077, -0.486, -0.8738, -13105.9336,
#                                      0,0,0,1), ncol = 4, byrow = T))
# } else if (ii_hc[1] == 3) {
#   segments3d(DRA_layers_ME_posterior_3, lwd=1)
#   rgl.viewpoint(fov=0,zoom=1,
#                 userMatrix= matrix(c(0.164, -0.752, -0.6374, -11506.360,
#                                      0.984,  0.1626,  0.062,   5351.796,
#                                      0.0569, -0.6379,  0.7679,      0.000,
#                                      0,0,0,1), ncol = 4, byrow = T))
# }
# 
# conn_in_R7 <- catmaid_get_connectors_between(pre_skids = skid_R7s, post_skids = skid[ii])
# if (!is.null(conn_in_R7)) {
#   conn_in_R7 <- conn_in_R7[, c("post_node_x", "post_node_y", "post_node_z")]
#   conn_in_R7_xform <- sweep(as.matrix(conn_in_R7), 2, DRA_me_pca$center) %*% DRA_me_pca$rotation
#   pch3d(conn_in_R7_xform, pch=16, radius=1000, alpha=0.7, col=pal_syn['R7'])
# }
# 
# # # save
# # rgl.snapshot(filename = paste("F5_", type, "_", skid[ii], "_side_posterial.png", sep = ''))

