# some useful functions ---------------------------------------------------


# make conn table for hist ------------------------------------------------

mktb_in_py78 <- function(skid){
  conn_x <- list()
  for (j in 1:4) {
    conn <- catmaid_get_connectors_between(pre_skids = skid_pyR7R8[(1:2)+2*(j-1)], post_skids = skid)
    if (is.null(conn)) {
      conn_x[[j]] <- NaN
    } else {
      conn <- conn[, c("post_node_x", "post_node_y", "post_node_z")]
      conn <- sweep(as.matrix(conn), 2, me_pca$center) %*% me_pca$rotation
      conn_x[[j]] <- conn[,1]
    }
  } 
  tb <- matrix(NaN, ncol = 4, nrow = max(sapply(conn_x, length)))
  for (j in 1:4) {
    tb[1:length(conn_x[[j]]),j] <- conn_x[[j]]
  }
  tb <- round(tb,2)
  tb <- as.data.frame(tb)
  colnames(tb) <- c('in_pR7', 'in_yR7', 'in_pR8', 'in_yR8')
  return(tb)
}

mktb_out_py78 <- function(skid){
  conn_x <- list()
  for (j in 1:4) {
    conn <- catmaid_get_connectors_between(pre_skids = skid, post_skids = skid_pyR7R8[(1:2)+2*(j-1)])
    if (is.null(conn)) {
      conn_x[[j]] <- NaN
    } else {
      conn <- conn[, c("post_node_x", "post_node_y", "post_node_z")]
      conn <- sweep(as.matrix(conn), 2, me_pca$center) %*% me_pca$rotation
      conn_x[[j]] <- conn[,1]
    }
  } 
  tb <- matrix(NaN, ncol = 4, nrow = max(sapply(conn_x, length)))
  for (j in 1:4) {
    tb[1:length(conn_x[[j]]),j] <- conn_x[[j]]
  }
  tb <- round(tb,2)
  tb <- as.data.frame(tb)
  colnames(tb) <- c('out_pR7', 'out_yR7', 'out_pR8', 'out_yR8')
  return(tb)
}


#  smooth 1-node jump -----------------------------------------------------
# assuming the first node of a segment is ok

# # - set threshold
# tar <- neu[[ii]]
# nopen3d()
# plot3d(tar)
# identify3d(xyzmatrix(tar$d))
# iim <- which(!is.na(sapply(tar$SegList, function(x) match(98, x))))
# aa <- xyzmatrix(tar$d)[tar$SegList[[iim]],] %>% diff() %>% diff() %>% .^2 %>% rowSums() %>% sqrt()
# hist(aa, 100)

# - func
fix_jump <- function(xneu){
  thr <- 2800
  # thr <- 1500
  for (jn in 1:length(xneu)) {
    xyz <- xyzmatrix(xneu[[jn]]$d) # node xyz
    # maxSeg <- which.max(sapply(xneu[[jn]]$SegList, length)) #longest seg
    # xyzmatrix(tar$d)[tar$SegList[[8]],] %>% diff() %>% diff() %>% .^2 %>% rowSums() %>% sqrt()
    for (js in 1:length(xneu[[jn]]$SegList)) {
      seg <- xneu[[jn]]$SegList[[js]]
      if (length(seg) > 2) {
        v <- diff(xyz[seg,])
        a <- diff(v)
        ii_jump <- which(sqrt(rowSums(a^2)) > thr)
        ii_0 <- seg[ii_jump]
        ii_1 <- seg[ii_jump + 1]
        ii_2 <- seg[ii_jump + 2]
        if (!identical(ii_0, integer(0))) {
          xneu[[jn]]$d[ii_1, c("X","Y","Z")] <- 
            (xneu[[jn]]$d[ii_0, c("X","Y","Z")] + xneu[[jn]]$d[ii_2, c("X","Y","Z")])/2
        }
      }
    }
  }
  return(xneu)
}

#  --(13)-- rigid/Euclidean transform neurons by pc  ---------------------------------------------------------------------------
xEucl_neu <- function(n, Rm, Tm){
  n_xform <- n
  for (j in 1:length(n)) {
    xyz <- n[[j]]$d[, c("X","Y","Z")]
    xyz_xform <- sweep(as.matrix(xyz), 2, Tm, '-') %*% Rm
    n_xform[[j]]$d[, c("X","Y","Z")] <- xyz_xform
  }
  return(n_xform)
}

# -- (6) -- rotate a point set ------------------------------------------------------------------------------------
# rotate a point set (eg T4 dendrite nodes) assuming a plane norm as z-axis (eg. plane normal to pc3) and +x direction
# input: Nx3 point set - pts, 1X6 norm vec - zv, 1x6 +x dir - xv, both zv and xv are [tail x, y ,z, head x, y, z]
# zv[1:3] is the origin
# output: dataframe [x,y,z] 

rot_zx <- function(pts, zv, xv){
  
  zv <- as.numeric(zv)
  xv <- as.numeric(xv)
  
  # first normalize zv and xv
  zvu <- zv[4:6] - zv[1:3]
  zvu <- zvu/sqrt(sum(zvu^2))
  xvu <- xv[4:6] - xv[1:3]
  xvu <- xvu - c(xvu %*% zvu)*zvu 
  xvu <- xvu/sqrt(sum(xvu^2))
  
  pts0 <- sweep(pts, 2, zv[1:3]) # vec of pts wrt the origin, zv[1:3]
  
  # find the rotation that makes zv3 z-axis and nx3 x-axis
  # first rotation about z to align zv3 to xz plane
  t1 <- 2*pi*(zvu[2]<0) + (-1)^(zvu[2]<0)*acos(zvu[1]/sqrt(sum(zvu[1:2]^2)))
  # t1 <- -t1 # cf wiki ???
  R1 <- matrix(c(cos(t1), -sin(t1), 0, sin(t1), cos(t1), 0, 0, 0, 1), ncol = 3)
  zvu1 <- R1 %*% zvu
  
  # second about y to align zvu to z
  t2 <- 2*pi*(zvu1[1]<0) + (-1)^(zvu1[1]<0)*acos(zvu1[3]/sqrt(sum(zvu1[c(1,3)]^2)))
  R2 <- matrix(c(cos(t2), 0, sin(t2), 0, 1, 0, -sin(t2), 0, cos(t2)), ncol = 3)
  zvu2 <- R2 %*% zvu1
  xvu2 <- R2 %*% R1 %*% xvu
  
  # third rotation about z to align central point on x
  t3 <- 2*pi*(xvu2[2]<0) + (-1)^(xvu2[2]<0)*acos(xvu2[1]/sqrt(sum(xvu2[1:2]^2)))
  R3 <- matrix(c(cos(t3), -sin(t3), 0, sin(t3), cos(t3), 0, 0, 0, 1), ncol = 3)
  xvu3 <- R3 %*% xvu2
  
  # apply to all points
  pts_rot <- t(R3 %*% R2 %*% R1 %*% t(pts0))
  pts__rot <- as.data.frame(pts_rot)
  colnames(pts_rot) <- c("X","Y","Z")
  
  return(pts_rot)
}

# -- (7) --  cross product ----------------------------------------------------------------------------------------

cross3D <- function(a, b){
  prod <- matrix(c(a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1]), ncol = 1)
  return(prod)
}


# -- () -- round corner
# an angle (3 points) on a plane, find bisect, 
# and draw inscribed circle of a chosen radius
# input pt_corner -- 3 pts define the corner. R -- radius
# output ptArc_rot -- pts sampling the arc

round_corner <- function(pt_corner, R = NA){
  pt_corner <- as.matrix(pt_corner)
  colnames(pt_corner) <- c('x','y','z')
  
  p3 <- sweep(pt_corner, 2, pt_corner[2,], '-')
  
  # rotate to z-x coord
  zv <- p3[1,]
  xv <- p3[3,] - c(p3[3,] %*% p3[1,])*p3[1,]
  p3r <- rot_zx(p3, c(0,0,0, zv), c(0,0,0,xv) )
  
  # bisect
  ang <- acos(c(p3r[1,] %*% p3r[3,]) / sqrt(sum(p3r[1,]^2)) / sqrt(sum(p3r[3,]^2))) / 2
  # radius
  if (is.na(R)) {
    R <- sqrt(sum(p3r[1,]^2))/2 / cos(ang) * sin(ang)
    R <- as.numeric(R)
  }
  
  # on z-x plane [z,x]
  cc <- c(sqrt(sum(p3r[1,]^2))/2, R) #center of circle
  angArc <- c(pi/2*3 - (pi - ang*2), pi/2*3) #angle range of the arc
  aa <- seq(angArc[1], angArc[2], length.out = 10) # sample angArc
  
  # in rotated frame [x,y,z]
  ptArc <- cbind(R * sin(aa) + cc[2],
                 0,
                 R * cos(aa) + cc[1])
  
  # inverse rotate
  axisxyz <- rot_zx(matrix(c(1,0,0,
                             0,1,0,
                             0,0,1), ncol = 3, byrow = T), c(0,0,0, zv), c(0,0,0,xv) )
  ptArc_rot <- rot_zx(ptArc, c(0,0,0,axisxyz[3,]), c(0,0,0,axisxyz[1,]))
  ptArc_rot <- sweep(ptArc_rot, 2, pt_corner[2,], '+')
  
  return(ptArc_rot)
}