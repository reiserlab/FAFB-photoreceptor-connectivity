# find Mi1's that identifies the seed columns -------------------------------------------

# - 4 col of 8 R7/8
pR78_me <- nlapply(c(pR7,pR8), function(x) subset(x, pointsinside(x, ME_msh))) #me portion 
yR78_me <- nlapply(c(yR7,yR8), function(x) subset(x, pointsinside(x, ME_msh)))
R78_4col_me <- c(pR78_me, yR78_me) #combine all 8
R78_4col_me_xform <- xEucl_neu(R78_4col_me, me_pca$rotation, me_pca$center)
ref_com <- sapply(R78_4col_me_xform, function(x) colMeans(xyzmatrix(x$d)))
ref_com <- t(ref_com)
ref_com[,1] <- 0
# group into 4
hc <- hclust(dist(ref_com))
hc_ind <- cutree(hc, k=4)
ref_com4 <- matrix(ncol = 3, nrow = 4)
for (j in 1:4) {
  xyz <- ref_com[hc_ind == j, ,drop=F]
  ref_com4[j,] <- colMeans(xyz)
}

# - R7 p and y
R7_p_xform <- xEucl_neu(R7_p, me_pca$rotation, me_pca$center)
R7_y_xform <- xEucl_neu(R7_y, me_pca$rotation, me_pca$center)

# -- Mi1 M5 nodes
xyz_M5_avg_xform <- sweep(xyz_M5_avg, 2, me_pca$center) %*% me_pca$rotation
xyz_M5_avg_xform_yz <- xyz_M5_avg_xform
xyz_M5_avg_xform_yz[,1] <- 0

xyz_M5_avg_xform_yz_Nnb <- xyz_M5_avg_xform_yz[c(Mi1_ind_Nnb_p,Mi1_ind_Nnb_y), ] #neighbor columns

# DEBUG
# nopen3d()
# points3d(xyz_M5_avg_xform_yz, col='grey', size = 20)
# # points3d(ref_com4,col='red',size=10)
# points3d(xyz_M5_avg_xform_yz[Mi1_ind_p,], size=25, col='blue')
# points3d(xyz_M5_avg_xform_yz[Mi1_ind_y,], size=25, col='gold')
# plot3d(R7_p_xform, col='green', lwd=2)
# plot3d(R7_y_xform, col='gold2', lwd=2)



# transform and generate meshes in the new coord given me_pca -------------

# - make a medulla mesh 
xyz <- t(ME_msh$vb[1:3,])
xyz <- sweep(xyz, 2, me_pca$center) %*% me_pca$rotation
ME_msh_xform <- ashape3d(xyz, alpha = 50000) %>% as.mesh3d()

# make new local mesh 
xyz <- t(ME_msh_xform$vb[1:3,])
dd <- sweep(xyz[,2:3], 2, colMeans(xyz_M5_avg_xform_yz_Nnb)[-1]) %>% .^2 %>% rowSums() %>% sqrt()
xyz_msh <- xyz[dd < 1.2*max(dist(xyz_M5_avg_xform_yz_Nnb)), ]
ME_msh_local_large <- ashape3d(xyz_msh, alpha = 60000) %>% as.mesh3d()

# a smaller version
xyz <- t(ME_msh_xform$vb[1:3,])
dd <- sweep(xyz[,2:3], 2, colMeans(xyz_M5_avg_xform_yz_Nnb)[-1]) %>% .^2 %>% rowSums() %>% sqrt()
xyz_msh <- xyz[dd < 0.8*max(dist(xyz_M5_avg_xform_yz_Nnb)), ]
ME_msh_local <- ashape3d(xyz_msh, alpha = 60000) %>% as.mesh3d()


# - transfomr LO mesh
LO_msh_xform <- LO_msh
LO_msh_xform$vb[1:3,] <- sweep(t(LO_msh$vb[1:3,]), 2, me_pca$center) %*% me_pca$rotation %>%  t()

# use Tm20 to make new local mesh 
Tm20_LO <- nlapply(Tm20, subset, function(x) pointsinside(x, LO_msh))
Tm20_LO_xform <- xEucl_neu(Tm20_LO, me_pca$rotation, me_pca$center)

LO_pca <- prcomp(xyzmatrix(Tm20_LO))

# make local mesh
xyz <- sweep(t(LO_msh$vb[1:3,]), 2, LO_pca$center) %*% LO_pca$rotation 
dd <- xyz[,2:3]^2 %>% rowSums() %>% sqrt() #distance to origin defined by Tm20

xyz <- sweep(t(LO_msh$vb[1:3,]), 2, me_pca$center) %*% me_pca$rotation
xyz_msh <-  xyz[dd < 1.2*max(dist(xyz_M5_avg_xform_yz_Nnb)), ]
LO_msh_local_large <- ashape3d(xyz_msh, alpha = 50000) %>% as.mesh3d()

xyz <- sweep(t(LO_msh$vb[1:3,]), 2, me_pca$center) %*% me_pca$rotation
xyz_msh <-  xyz[dd < 0.8*max(dist(xyz_M5_avg_xform_yz_Nnb)), ]
LO_msh_local <- ashape3d(xyz_msh, alpha = 50000) %>% as.mesh3d()


# - transfomr other mesh
AOTU_msh_xform <- AOTU_msh
AOTU_msh_xform$vb[1:3,] <- sweep(t(AOTU_msh$vb[1:3,]), 2, me_pca$center) %*% me_pca$rotation %>%  t()

LOP_L_msh_xform <- LOP_L_msh
LOP_L_msh_xform$vb[1:3,] <- sweep(t(LOP_L_msh$vb[1:3,]), 2, me_pca$center) %*% me_pca$rotation %>%  t()

FLA_R_msh_xform <- FLA_R_msh
FLA_R_msh_xform$vb[1:3,] <- sweep(t(FLA_R_msh$vb[1:3,]), 2, me_pca$center) %*% me_pca$rotation %>%  t()

IPS_L_msh_xform <- IPS_L_msh
IPS_L_msh_xform$vb[1:3,] <- sweep(t(IPS_L_msh$vb[1:3,]), 2, me_pca$center) %*% me_pca$rotation %>%  t()

IPS_R_msh_xform <- IPS_R_msh
IPS_R_msh_xform$vb[1:3,] <- sweep(t(IPS_R_msh$vb[1:3,]), 2, me_pca$center) %*% me_pca$rotation %>%  t()

SPS_R_msh_xform <- SPS_R_msh
SPS_R_msh_xform$vb[1:3,] <- sweep(t(SPS_R_msh$vb[1:3,]), 2, me_pca$center) %*% me_pca$rotation %>%  t()

PLP_R_msh_xform <- PLP_R_msh
PLP_R_msh_xform$vb[1:3,] <- sweep(t(PLP_R_msh$vb[1:3,]), 2, me_pca$center) %*% me_pca$rotation %>%  t()

AME_R_msh_xform <- AME_R_msh
AME_R_msh_xform$vb[1:3,] <- sweep(t(AME_R_msh$vb[1:3,]), 2, me_pca$center) %*% me_pca$rotation %>%  t()

# medulla layers ----------------------------------------------------------------
# first find boundaries of the neuropil, then draw layers in between them using 7-col paper and some known neurons

# - this transformation is related to later plots, such that the layers are mostly parallel to the plane of the plots
xAng <- -105 #rotation angle around x-axis
xAngRot <- (xAng - 90)/180*pi
xRot <- matrix(c(1,0,0,
                 0, cos(xAngRot), sin(xAngRot),
                 0, -sin(xAngRot), cos(xAngRot)), ncol = 3, byrow = T)
xyz <- t(ME_msh_local$vb[1:3,])
xz <- xyz %*% xRot %>% .[,c(1,3)]
xyz_chull <- cbind(xz[chull(xz),1], 0, xz[chull(xz),2])

# scale bar
sbar <- matrix(c(-60000,0,-40000, -60000,0,-30000), ncol=3,byrow=T) %*% t(xRot)

# - ME boundary
# fit a surface to the mesh boundary and them sample the fit in the desired positions
vbd <- xyz_chull %*% t(xRot)
yz <- data.frame(y = seq(min(vbd[,2]), max(vbd[,2]), length.out = 20),
                 z = seq(min(vbd[,3]), max(vbd[,3]), length.out = 20))

# points on mesh
xyz <- t(ME_msh_xform$vb[1:3,])
dd <- sweep(xyz[,2:3], 2, colMeans(xyz_M5_avg_xform_yz_Nnb)[-1]) %>% .^2 %>% rowSums() %>% sqrt()
xyz_msh <- xyz[dd < max(dist(xyz_M5_avg_xform_yz_Nnb)), ]
# top
xyz_top <- xyz_msh[xyz_msh[,1] > 0, ]
x <- xyz_top[,1]; y <- xyz_top[,2]; z <- xyz_top[,3]
fitlm <- lm(x ~ poly(y, z, degree = 2, raw = T))
valfit <- predict(fitlm, yz) #generate values from the fit
top_bd <- cbind(valfit, yz)
# bottom
xyz_bot <- xyz_msh[xyz_msh[,1] < 0, ]
x <- xyz_bot[,1]; y <- xyz_bot[,2]; z <- xyz_bot[,3]
fitlm <- lm(x ~ poly(y, z, degree = 2, raw = T))
valfit <- predict(fitlm, yz) #generate values from the fit
bot_bd <- cbind(valfit, yz)

# 2D mesh as outline
ME_msh_local_bd <- rbind(top_bd, apply(bot_bd, MARGIN = 2, rev)) 
ME_msh_local_bd <- rbind(ME_msh_local_bd[-1,], ME_msh_local_bd[1,]) %>%
  cbind(ME_msh_local_bd, .) %>%
  t() %>%
  matrix(., ncol = 3, byrow = T)


# --layer definition from 7-column paper
# https://www.pnas.org/content/112/44/13711?ijkey=682c42b2318fb1efb9fe6d7ad6da91752ce3fa22&keytype2=tf_ipsecsha
# bp7c <- c(0, 5.5,14.5,20,23,28,36,42,48,56,61) #layer 1 to 10
bp7c <- c(0, 5.5, 14.5, 20, 23, 28, 33, 39, 45, 56, 61) #layer 1 to 10, modify 6-7-8-9 by Mi1, C2, Tm5 and Tm20
bp7c <- max(bp7c) - bp7c %>% rev() #layer 10 to 1
bp7c_prob <- bp7c/max(bp7c)
x_qua <- apply(cbind(top_bd$valfit, bot_bd$valfit), MARGIN = 1,
               function(x) quantile(x, probs = bp7c_prob) )
layers <- matrix(ncol = 3, nrow = 0)
for (j in 1:length(bp7c_prob)) {
  cc <- cbind(x_qua[j,], yz)
  c2 <- cbind(cc[-nrow(cc),], cc[-1,]) %>%
    t() %>%
    matrix(., ncol = 3, byrow = T)
  layers <- rbind(layers, c2)
}
layers[,1] <- layers[,1] * 0.97 + 1100

bd_Mx <- x_qua[,ncol(x_qua)] %>% rev()
bd_Mx <- bd_Mx * 0.97 + 1100
bd_Mx <- bd_Mx[-length(bd_Mx)] + diff(bd_Mx)/2
layer_M_anno <- cbind(bd_Mx, yz[nrow(yz),])


# layers side view
xAngRot_layer <- 30/180*pi
xRot_layer <- matrix(c(1,0,0,
                       0, cos(xAngRot_layer), sin(xAngRot_layer),
                       0, -sin(xAngRot_layer), cos(xAngRot_layer)), ncol = 3, byrow = T)
layers_side <- layers %*% (xRot_layer)

# LObula layer -------------------------------------------------------------

xyz_msh <- t(LO_msh_local$vb[1:3,])

# right / top
xyz_top <- xyz_msh[xyz_msh[,3] > - 40000, ]
i1 <- which.min(xyz_top[,1])
i2 <- which.max(xyz_top[,1])
xy <- data.frame(x = seq(xyz_top[i1,1], xyz_top[i2,1], length.out = 20),
                 y = seq(mean(xyz_top[c(i1,i2),2]), mean(xyz_top[c(i1,i2),2]), length.out = 20) )
x <- xyz_top[,1]; y <- xyz_top[,2]; z <- xyz_top[,3]
fitlm <- lm(z ~ poly(x, y, degree = 2, raw = T))
valfit <- predict(fitlm, xy) #generate values from the fit
top_bd <- data.frame(x=xy$x, y=xy$y, z=valfit)

# bottom
xyz_bot <- xyz_msh[xyz_msh[,3] < - 40000, ]
i1 <- which.min(xyz_bot[,1])
i2 <- which.max(xyz_bot[,1])
xy <- data.frame(x = seq(xyz_bot[i1,1], xyz_bot[i2,1], length.out = 20),
                 y = seq(mean(xyz_bot[c(i1,i2),2]), mean(xyz_bot[c(i1,i2),2]), length.out = 20) )
x <- xyz_bot[,1]; y <- xyz_bot[,2]; z <- xyz_bot[,3]
fitlm <- lm(z ~ poly(x, y, degree = 2, raw = T))
valfit <- predict(fitlm, xy) #generate values from the fit
bot_bd <- data.frame(x=xy$x, y=xy$y, z=valfit)

# -- layer Wu 2016
# https://elifesciences.org/articles/21022
# bp <- c(0, 0.8, 1.35, 2.1, 2.85, 3.8, 4.85, 6.7)
bp <- c(0, 0.8, 1.35, 2.1, 2.85, 3.8*0.6, 4.85*0.6, 6.7) # mod by Tm5a
bp_prob <- cumsum(bp)/sum(bp)
# bp_prob <- 1 - bp_prob
x_qua <- matrix(ncol = nrow(top_bd), nrow = length(bp_prob))
y_qua <- matrix(ncol = nrow(top_bd), nrow = length(bp_prob))
z_qua <- matrix(ncol = nrow(top_bd), nrow = length(bp_prob))
for (j in 1:length(bp_prob)) {
  x_qua[j,] <- top_bd$x + (bot_bd$x - top_bd$x)*bp_prob[j]
  y_qua[j,] <- top_bd$y + (bot_bd$y - top_bd$y)*bp_prob[j]
  z_qua[j,] <- top_bd$z + (bot_bd$z - top_bd$z)*bp_prob[j]
}
layers_LO <- matrix(ncol = 3, nrow = 0)
for (j in 1:length(bp_prob)) {
  cc <- cbind(x_qua[j,], y_qua[j,], z_qua[j,])
  c2 <- cbind(cc[-nrow(cc),], cc[-1,]) %>%
    t() %>%
    matrix(., ncol = 3, byrow = T)
  layers_LO <- rbind(layers_LO, c2)
}
layers_LO[,3] <- layers_LO[,3] - 2000 # manual correction

bd_Mx <- x_qua[,1]
bd_Mx <- bd_Mx[-length(bd_Mx)] + diff(bd_Mx)/2
bd_My <- y_qua[,1] 
bd_My <- bd_My[-length(bd_My)] + diff(bd_My)/2
bd_Mz <- z_qua[,1] 
bd_Mz <- bd_Mz[-length(bd_Mz)] + diff(bd_Mz)/2

layer_LO_anno <- cbind(bd_Mx + c(0,-1,-2,-2,-1.5,-1, -1)*2000, bd_My, bd_Mz)

# generate cylinder shades for plots -----------------------------------------------------

t(ME_msh_local$vb)[,1] %>% range()

x_top <- 38000
x_bot <- -48000

cp1 <- cylinder3d(cbind(seq(x_bot, x_top, by=5000),
                        xyz_M5_avg_xform_yz[Mi1_ind_p[1],2],
                        xyz_M5_avg_xform_yz[Mi1_ind_p[1],3]), radius = 2200)
cp2 <- cylinder3d(cbind(seq(x_bot, x_top, by=5000),
                        xyz_M5_avg_xform_yz[Mi1_ind_p[2],2],
                        xyz_M5_avg_xform_yz[Mi1_ind_p[2],3]), radius = 2200)
cy1 <- cylinder3d(cbind(seq(x_bot, x_top, by=5000),
                        xyz_M5_avg_xform_yz[Mi1_ind_y[1],2],
                        xyz_M5_avg_xform_yz[Mi1_ind_y[1],3]), radius = 2200)
cy2 <- cylinder3d(cbind(seq(x_bot, x_top, by=5000),
                        xyz_M5_avg_xform_yz[Mi1_ind_y[2],2],
                        xyz_M5_avg_xform_yz[Mi1_ind_y[2],3]), radius = 2200)