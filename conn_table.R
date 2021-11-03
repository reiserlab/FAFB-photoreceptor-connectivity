# This code query CATMAID to construct various connectivity tables
# "table by cell" sections save .csv files for each cell type, which will be used later for plotting

# WARNING: saving is by default enabled. Search for "save" and "write.csv" to see how data is saved

# load conn data ----------------------------------------------------------

load("data/conn_R7R8.RData")

# conn table cutoff value
Ncutoff <- 1 # >= 2 or 3

# type as rows ----------------------------------------------------------------
# func to make table by type 

mktb_by_type <- function(tb_by_name){
  tb_by_type <- tb_by_name %>% 
    as_tibble() %>%
    transmute(pR7, yR7, pR8, yR8, type) %>%
    group_by(type) %>%
    mutate(N = n()) %>%
    summarise(pR7 = sum(pR7), yR7 = sum(yR7), pR8 = sum(pR8), yR8 = sum(yR8), N = mean(N)) %>%
    # transmute(type = paste(type, ' (', N, ')', sep = ''), pR7, yR7, pR8, yR8, total= pR7+yR7+pR8+yR8, of_target = NA)
    transmute(type = type, num = N, pR7, yR7, pR8, yR8, total= pR7+yR7+pR8+yR8, of_target = NA)
  cs <- colSums(tb_by_type[,-1]) #col sum
  tb_by_type %<>% 
    mutate(perc_Total = round(total/cs['total']*100,3),
           perc_Total_R7 = round((pR7+yR7)/(cs['pR7']+cs['yR7'])*100,3),
           perc_Total_R8 = round((pR8+yR8)/(cs['pR8']+cs['yR8'])*100,3),
           perc_Total_p = round((pR7+pR8)/(cs['pR7']+cs['pR8'])*100,3),
           perc_Total_y = round((yR7+yR8)/(cs['yR7']+cs['yR8'])*100,3)) %>%
    mutate(of_R7R8 = NA) %>%
    data.frame()
  # order by total
  tb_by_type <- tb_by_type[rev(order(tb_by_type$total)),]
  tb_by_type <- rbind(tb_by_type, c(0, colSums(tb_by_type[,-1])) )
  tb_by_type[dim(tb_by_type)[1],1] <-  'Total'
  tb_by_type %<>% as_tibble() %>%
    mutate(perc_Total = round(perc_Total,1),
           perc_Total_R7 = round(perc_Total_R7,1),
           perc_Total_R8 = round(perc_Total_R8,1),
           perc_Total_p = round(perc_Total_p,1),
           perc_Total_y = round(perc_Total_y,1)) %>%
    mutate(perc_R7 = round((pR7+yR7)/total*100,1),
           perc_R8 = round((pR8+yR8)/total*100,1),
           perc_p = round((pR7+pR8)/total*100,1),
           perc_y = round((yR7+yR8)/total*100,1) ) %>%
    transmute(type, num, pR7, yR7, pR8, yR8, sum=total,
              perc_R7, perc_R8, perc_p, perc_y,
              perc_Total, perc_Total_R7, perc_Total_R8, perc_Total_p, perc_Total_y) %>%
    data.frame()
  
  ii <- c()
  tmp <- matrix(ncol = ncol(tb_by_type), nrow= 0)
  colnames(tmp) <- colnames(tb_by_type)
  tmp2 <- matrix(0, ncol = ncol(tb_by_type)-1, nrow= 1)
  colnames(tmp2) <- colnames(tb_by_type)[-1]
  for (ii_type in c("Identified_<3*", "Unidentified_>=3*", "Unidentified_<3*")) {
    i3 <- grepl(ii_type, tb_by_type$type) %>% which()
    ii <- c(ii, i3)
    if (length(i3) == 0) {
      tmp <- rbind(tmp, data.frame(type = gsub("\\*", "", ii_type), tmp2))
    } else {
      tmp <- rbind(tmp, tb_by_type[i3,])
    }
  }
  
  rr <- nrow(tb_by_type)
  tb_by_type <- rbind(tb_by_type[-c(ii,rr),], tmp, tb_by_type[rr,])
  
  colnames(tb_by_type) <- gsub("perc_", "%", colnames(tb_by_type))
  
  return(as.data.frame(tb_by_type))
}

mktb_by_type_DRA <- function(tb_by_name){
  tb_by_type_DRA <- tb_by_name %>%
    as_tibble() %>% 
    # tb_by_type_DRA <- DRAR7R8_in_name %>%
    transmute(DRAR7,DRAR8, type) %>%
    group_by(type) %>%
    mutate(N = n()) %>%
    summarise(DRAR7 = sum(DRAR7), DRAR8 = sum(DRAR8), N = mean(N)) %>%
    mutate(of_R = NA) %>%
    # transmute(type = paste(type, ' (', N, ')', sep = ''), DRAR7, DRAR8, total=DRAR7+DRAR8)
    transmute(type = type, num = N, DRAR7, DRAR8, total=DRAR7+DRAR8)
  cs <- colSums(tb_by_type_DRA[,-1]) #col sum
  tb_by_type_DRA %<>% 
    mutate(perc_Total = round(total/cs['total']*100,3),
           perc_Total_R7 = round(DRAR7/cs['DRAR7']*100,3),
           perc_Total_R8 = round(DRAR8/cs['DRAR8']*100,3)) %>%
    data.frame()
  
  # order by total
  tb_by_type_DRA <- tb_by_type_DRA[rev(order(tb_by_type_DRA$total)),]
  tb_by_type_DRA <- rbind(tb_by_type_DRA, c(0, colSums(tb_by_type_DRA[,-1])) )
  tb_by_type_DRA[dim(tb_by_type_DRA)[1],1] <-  'Total'
  tb_by_type_DRA %<>% as_tibble() %>%
    mutate(perc_Total = round(perc_Total,1),
           perc_Total_R7 = round(perc_Total_R7,1),
           perc_Total_R8 = round(perc_Total_R8,1)) %>%
    mutate(perc_R7 = round(DRAR7/total*100,1),
           perc_R8 = round(DRAR8/total*100,1) ) %>%
    transmute(type, num, DRAR7, DRAR8, sum=total,
              perc_R7, perc_R8,
              perc_Total, perc_Total_R7, perc_Total_R8) %>%
    data.frame()

  ii <- c()
  tmp <- matrix(ncol = ncol(tb_by_type_DRA), nrow= 0)
  colnames(tmp) <- colnames(tb_by_type_DRA)
  tmp2 <- matrix(0, ncol = ncol(tb_by_type_DRA)-1, nrow= 1)
  colnames(tmp2) <- colnames(tb_by_type_DRA)[-1]
  for (ii_type in c("Identified_<3*", "Unidentified_>=3*", "Unidentified_<3*")) {
    i3 <- grepl(ii_type, tb_by_type_DRA$type) %>% which()
    ii <- c(ii, i3)
    if (length(i3) == 0) {
      tmp <- rbind(tmp, data.frame(type = gsub("\\*", "", ii_type), tmp2))
    } else {
      tmp <- rbind(tmp, tb_by_type_DRA[i3,])
    }
  }
  
  rr <- nrow(tb_by_type_DRA)
  tb_by_type_DRA <- rbind(tb_by_type_DRA[-c(ii,rr),], tmp, tb_by_type_DRA[rr,])
  
  colnames(tb_by_type_DRA) <- gsub("perc_", "%", colnames(tb_by_type_DRA))
  colnames(tb_by_type_DRA) <- gsub("total", "sum", colnames(tb_by_type_DRA))

  return(as.data.frame(tb_by_type_DRA))
}



# conn table by type  -----------------------------------------------
# this code block will query CATMAID

# - central out
R7R8_out_name <- R7R8_out_ori %>%
  mutate(name = catmaid_get_neuronnames(partner)) %>%
  mutate(type = if_else(grepl("(P|p)utative", name), word(name, 2),"Unidentified")) %>%
  mutate(total = pR7+yR7+pR8+yR8) %>%
  arrange(desc(name)) %>%
  # filter(total >= Ncutoff) %>%
  mutate(type = if_else((total < Ncutoff & type != "Unidentified"), "Identified_<3", type)) %>%
  mutate(type = if_else((total < Ncutoff & type == "Unidentified"), "Unidentified_<3", type)) %>%
  mutate(type = if_else((total >= Ncutoff & type == "Unidentified"), "Unidentified_>=3", type)) %>%
  # transmute(pR7, yR7, pR8, yR8, type, name) %>%
  data.frame()
R7R8_out_type <- mktb_by_type(R7R8_out_name)

# - central in
R7R8_in_name <- R7R8_in_ori %>%
  mutate(name = catmaid_get_neuronnames(partner)) %>%
  mutate(type = if_else(grepl("(P|p)utative", name), word(name, 2),"Unidentified")) %>%
  mutate(total = pR7+yR7+pR8+yR8) %>%
  mutate(type = if_else((total < Ncutoff & type != "Unidentified"), "Identified_<3", type)) %>%
  mutate(type = if_else((total < Ncutoff & type == "Unidentified"), "Unidentified_<3", type)) %>%
  mutate(type = if_else((total >= Ncutoff & type == "Unidentified"), "Unidentified_>=3", type)) %>%
  data.frame()
R7R8_in_type <- mktb_by_type(R7R8_in_name)


# - DRA in
DRAR7R8_in_name <- DRAR7R8_in_ori %>%
  mutate(name = catmaid_get_neuronnames(partner)) %>%
  mutate(type = if_else(grepl("(P|p)utative", name), word(name, 2),"Unidentified")) %>%
  mutate(total = rowSums(.[2:7])) %>%
  mutate(type = if_else((total < Ncutoff & type != "Unidentified"), "Identified_<3", type)) %>%
  mutate(type = if_else((total < Ncutoff & type == "Unidentified"), "Unidentified_<3", type)) %>%
  mutate(type = if_else((total >= Ncutoff & type == "Unidentified"), "Unidentified_>=3", type)) %>%
  data.frame()
DRAR7R8_in_type <- mktb_by_type_DRA(DRAR7R8_in_name)

# - DRA out
DRAR7R8_out_name <- DRAR7R8_out_ori %>%
  mutate(name = catmaid_get_neuronnames(partner)) %>%
  mutate(type = if_else(grepl("^(P|p)utative", name), word(name, 2),"Unidentified")) %>%
  mutate(total = rowSums(.[2:7])) %>%
  arrange(desc(name)) %>%
  mutate(type = if_else((total < Ncutoff & type != "Unidentified"), "Identified_<3", type)) %>%
  mutate(type = if_else((total < Ncutoff & type == "Unidentified"), "Unidentified_<3", type)) %>%
  mutate(type = if_else((total >= Ncutoff & type == "Unidentified"), "Unidentified_>=3", type)) %>%
  # transmute(DRAR7, DRAR8, type, name) %>%
  data.frame()
DRAR7R8_out_type <- mktb_by_type_DRA(DRAR7R8_out_name)


# SAVE 
# conn_data.RData file is used by "*.Rmd" files for making tables and galleries
# *.csv files is used by "Figure_x.R" files

save(R7R8_in_type, R7R8_out_type, DRAR7R8_in_type, DRAR7R8_out_type,
     R7R8_in_name, R7R8_out_name, DRAR7R8_in_name, DRAR7R8_out_name, file = "data/conn_data.RData")

write.csv(R7R8_in_type, "R7R8_incoming_type.csv", fileEncoding = "UTF-8", na="")
write.csv(R7R8_out_type, "R7R8_outgoing_type.csv", fileEncoding = "UTF-8", na="")
write.csv(DRAR7R8_in_type, "DRAR7R8_incoming_type.csv", fileEncoding = "UTF-8")
write.csv(DRAR7R8_out_type, "DRAR7R8_outgoing_type.csv", fileEncoding = "UTF-8")


# table by cell ---------------------------------------------------------------------------------------------------
# The next two code blocks will save a conn table for each cell type as a .csv file


# -incoming, create a folder if not exist
foldername <- "table_by_type_in"
if (dir.exists(foldername)) {
  setwd(foldername)
} else {
  dir.create(foldername)
  setwd(foldername)
}

# loop over cell types
for (type_ind in R7R8_in_type$type) {
  type_name <- word(type_ind, 1)
  num <- word(type_ind, 2)
  if (type_name %in% R7R8_in_name$type) {
    # get partner info
    R7R8_in <- R7R8_in_name %>%
      filter(type == type_name) %>%
      transmute(pR7a, pR7b, yR7a, yR7b, pR8a, pR8b, yR8a, yR8b, total, name, skid = partner) %>%
      arrange(desc(total)) %>%
      rbind(colSums(.[,1:9])) %>%
      relocate(name, skid) %>%
      mutate(perc_R7 = round((pR7a+pR7b+yR7a+yR7b)/total*100,1),
             perc_R8 = round((pR8a+pR8b+yR8a+yR8b)/total*100,1),
             perc_p = round((pR7a+pR7b+pR8a+pR8b)/total*100,1),
             perc_y = round((yR7a+yR7b+yR8a+yR8b)/total*100,1)) %>%
      data.frame()
    R7R8_in[nrow(R7R8_in),1] <-  'Total'
    R7R8_in[nrow(R7R8_in),2] <-  ''
    colnames(R7R8_in)[12:15] <- c("%R7", "%R8", "%p", "%y")
    
    if (grepl(3, type_name)) {
      new_name <- gsub("<3", "2-", type_name)
      new_name <- gsub(">=3", "3+", new_name)
      write.csv(R7R8_in, paste("R7R8_incoming_", new_name, ".csv",sep = ''), row.names=F, fileEncoding = "UTF-8")
    } else {
      write.csv(R7R8_in, paste("R7R8_incoming_", type_name, ".csv",sep = ''), row.names=F, fileEncoding = "UTF-8") 
    }
  }
}
setwd("../")


# - outgoing, loop over cell types
foldername <- "table_by_type_out"
if (dir.exists(foldername)) {
  setwd(foldername)
} else {
  dir.create(foldername)
  setwd(foldername)
}

for (type_ind in R7R8_out_type$type) {
  type_name <- word(type_ind, 1)
  num <- word(type_ind, 2)
  if (type_name %in% R7R8_out_name$type) {
    # get partner info
    R7R8_out <- R7R8_out_name %>%
      filter(type == type_name) %>%
      transmute(pR7a, pR7b, yR7a, yR7b, pR8a, pR8b, yR8a, yR8b, total, name, skid = partner) %>%
      arrange(desc(total)) %>%
      rbind(colSums(.[,1:9])) %>%
      relocate(name, skid) %>%
      mutate(perc_R7 = round((pR7a+pR7b+yR7a+yR7b)/total*100,1),
             perc_R8 = round((pR8a+pR8b+yR8a+yR8b)/total*100,1),
             perc_p = round((pR7a+pR7b+pR8a+pR8b)/total*100,1),
             perc_y = round((yR7a+yR7b+yR8a+yR8b)/total*100,1)) %>%
      data.frame()
    R7R8_out[nrow(R7R8_out),1] <-  'Total'
    R7R8_out[nrow(R7R8_out),2] <-  ''
    colnames(R7R8_out)[12:15] <- c("%R7", "%R8", "%p", "%y")
    
    if (grepl(3, type_name)) {
      new_name <- gsub("<3", "2-", type_name)
      new_name <- gsub(">=3", "3+", new_name)
      write.csv(R7R8_out, paste("R7R8_outgoing_", new_name, ".csv",sep = ''), row.names=F, fileEncoding = "UTF-8")
    } else {
      write.csv(R7R8_out, paste("R7R8_outgoing_", type_name, ".csv",sep = ''), row.names=F, fileEncoding = "UTF-8") 
    }
  }
}

setwd("../")



# table by cell DRA -------------------------------------------------------

# - incoming, loop over cell types
foldername <- "table_by_type_in_DRA"
if (dir.exists(foldername)) {
  setwd(foldername)
} else {
  dir.create(foldername)
  setwd(foldername)
}

for (type_ind in DRAR7R8_in_type$type) {
  type_name <- word(type_ind, 1)
  num <- word(type_ind, 2)
  if (type_name %in% DRAR7R8_in_name$type) {
    # get partner info
    DRAR7R8_in <- DRAR7R8_in_name %>%
      filter(type == type_name) %>%
      transmute(DRAR7a, DRAR7b, DRAR7c, DRAR8a, DRAR8b, DRAR8c, total, name, skid = partner) %>%
      arrange(desc(total)) %>%
      rbind(colSums(.[,1:7])) %>%
      relocate(name, skid) %>%
      mutate(perc_R7 = round((DRAR7a+DRAR7b+DRAR7c)/total*100,1),
             perc_R8 = round((DRAR8a+DRAR8b+DRAR8c)/total*100,1)) %>%
      data.frame()
    DRAR7R8_in[nrow(DRAR7R8_in),1] <-  'Total'
    DRAR7R8_in[nrow(DRAR7R8_in),2] <-  ''
    colnames(DRAR7R8_in)[10:11] <- c("%R7", "%R8")
    
    if (grepl(3, type_name)) {
      new_name <- gsub("<3", "2-", type_name)
      new_name <- gsub(">=3", "3+", new_name)
      write.csv(DRAR7R8_in, paste("DRAR7R8_incoming_", new_name, ".csv",sep = ''), row.names=F, fileEncoding = "UTF-8")
    } else {
      write.csv(DRAR7R8_in, paste("DRAR7R8_incoming_", type_name, ".csv",sep = ''), row.names=F, fileEncoding = "UTF-8") 
    }
  }
}
setwd("../")



# - outgoing, loop over cell types
foldername <- "table_by_type_out_DRA"
if (dir.exists(foldername)) {
  setwd(foldername)
} else {
  dir.create(foldername)
  setwd(foldername)
}

for (type_ind in DRAR7R8_out_type$type) {
  type_name <- word(type_ind, 1)
  num <- word(type_ind, 2)
  if (type_name %in% DRAR7R8_out_name$type) {
    # get partner info
    DRAR7R8_out <- DRAR7R8_out_name %>%
      filter(type == type_name) %>%
      transmute(DRAR7a, DRAR7b, DRAR7c, DRAR8a, DRAR8b, DRAR8c, total, name, skid = partner) %>%
      arrange(desc(total)) %>%
      rbind(colSums(.[,1:7])) %>%
      relocate(name, skid) %>%
      mutate(perc_R7 = round((DRAR7a+DRAR7b+DRAR7c)/total*100,1),
             perc_R8 = round((DRAR8a+DRAR8b+DRAR8c)/total*100,1)) %>%
      data.frame()
    DRAR7R8_out[nrow(DRAR7R8_out),1] <-  'Total'
    DRAR7R8_out[nrow(DRAR7R8_out),2] <-  ''
    colnames(DRAR7R8_out)[10:11] <- c("%R7", "%R8")
    
    if (grepl(3, type_name)) {
      new_name <- gsub("<3", "2-", type_name)
      new_name <- gsub(">=3", "3+", new_name)
      write.csv(DRAR7R8_out, paste("DRAR7R8_outgoing_", new_name, ".csv",sep = ''), row.names=F, fileEncoding = "UTF-8")
    } else {
      write.csv(DRAR7R8_out, paste("DRAR7R8_outgoing_", type_name, ".csv",sep = ''), row.names=F, fileEncoding = "UTF-8") 
    }
  }
}
setwd("../")
