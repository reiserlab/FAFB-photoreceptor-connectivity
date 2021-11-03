# For DRA plots, there are two views
# For polar view, 
# Run "startup.R" and "layer_DRA_polar.R"
# scale bar is 10 um

# polar DmDRA1 in Figure 7 as an example ---------------------------------------------------------------------


type <- "Dm-DRA1"

tb <- read.csv(paste("table_by_type_out_DRA/", "DRAR7R8_outgoing_", type, ".csv", sep = ''))

skid <- tb$skid %>% na.omit()
neu <-  read.neurons.catmaid(skid, .progress='text' )
neu_xform <- xEucl_neu(neu, DRA_me_pca$rotation, DRA_me_pca$center)

ii_match <- skid %in% 12106449 %>% which()
ii <- ii_match[1]

# - topview - seed column input
nopen3d()
par3d('windowRect' = c(10, 10, 2110, 2110))
plot3d(neu_xform[[ii]], col='black', lwd=2, soma = T, WithNodes = F, lit=F)
pch3d(xyz_M5_avg_xform_DRA_yz[ind_Mi1_DRA_hcol,,drop=F], radius=5000,col="#b7252a",pch=16,alpha=0.7)
pch3d(xyz_M5_avg_xform_DRA_yz[ind_Mi1_top_grey,], radius=5000,col="grey",pch=16,alpha=0.2)
pch3d(xyz_M5_avg_xform_DRA_yz[ind_Mi1_top_red,], radius=5000,col="#b7252a",pch=16,alpha=0.3)
segments3d(sweep(rbind(c(0, 0, 0),
                       rotate3d(rotate3d(rotate3d(c(10000, 0, 0), -80/180*pi, 0,0,1), -90/180*pi,0,1,0), DRA_xAng/180*pi,1,0,0)), 2, c(-35000,10000,-10000)), lwd=2)

rgl.viewpoint(fov=0,zoom=0.8, userMatrix= t(translationMatrix(5.9e3,-1.1e4,-6.1e-13)) %*%
                rotationMatrix(-80/180*pi,0,0,1) %*%
                rotationMatrix(-90/180*pi,0,1,0) %*%
                rotationMatrix(DRA_xAng/180*pi,1,0,0) )

conn_in_R7 <- catmaid_get_connectors_between(pre_skids = skid_DRAR7[1], post_skids = skid[ii])
if (!is.null(conn_in_R7)) {
  conn_in_R7 <- conn_in_R7[, c("connector_x", "connector_y", "connector_z")]
  conn_in_R7_xform <- sweep(as.matrix(conn_in_R7), 2, DRA_me_pca$center) %*% DRA_me_pca$rotation
  pch3d(conn_in_R7_xform, pch=16, radius=1000, alpha=0.7, col=pal_syn['R7'])
}

# rgl.snapshot(filename = paste("F7_", type, "_", skid[ii], "_seed_column_input_top.png", sep = ''))


# - sideview
nopen3d()
par3d('windowRect' = c(10, 10, 1510, 2110))
plot3d(neu_xform[[ii]], col='black', lwd=2, soma = T, WithNodes = F, lit=F)
plot3d(DRAR7_xform[[1]], col=pal_syn["R7"], lwd=2, soma = T, WithNodes = F, lit=F)
segments3d(DRA_layers_ME_posterior_rot90, lwd=1)
# segments3d(sweep(DRA_sbar_side,2,c(0,0,-20000)), lwd=2)
segments3d(sweep(rbind(c(0, 0, 0),
                       rotate3d(rotate3d(rotate3d(c(10000, 0, 0), 95/180*pi, 0,0,1), -1*120/180*pi, 1,0,0), 1*20/180*pi, 0,1,0)
), 2, c(-35000,10000,0)), 
lwd=2)
# text3d(colMeans(DRA_sbar_side)+c(5000,-5000,25000), texts = "10 um", cex = 1.5)
rgl.viewpoint(fov=0,zoom=1, userMatrix= rotationMatrix(95/180*pi,0,0,1) %*%
                rotationMatrix(1*(-30-90)/180*pi,1,0,0) %*%
                rotationMatrix(1*20/180*pi,0,1,0))

conn_in_R7 <- catmaid_get_connectors_between(pre_skids = skid_DRAR7[1], post_skids = skid[ii])
if (!is.null(conn_in_R7)) {
  conn_in_R7 <- conn_in_R7[, c("connector_x", "connector_y", "connector_z")]
  conn_in_R7_xform <- sweep(as.matrix(conn_in_R7), 2, DRA_me_pca$center) %*% DRA_me_pca$rotation
  pch3d(conn_in_R7_xform, pch=16, radius=1000, alpha=0.7, col=pal_syn['R7'])
}

# rgl.snapshot(filename = paste("F7_", type, "_", skid[ii], "_side.png", sep = ''))
