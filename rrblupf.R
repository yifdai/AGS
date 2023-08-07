##################################
# This is the module for rrblup. #
##################################

## origin little change
c.z.hglm <- function(kin){
  relmat <- kin
  relmat[upper.tri(relmat)] <- t(relmat)[upper.tri(relmat)]
  cat("relmat")
  print(dim(relmat))
  svd <- svd(relmat)
  Z <- svd$u %*% diag(sqrt(svd$d))
  print(dim(Z))
  return(Z)
}


prep_input <- function(input1, input2, sample_ID, trait){
  # na <- which(is.na(input1[,trait]))
  ph2 <- input1[!is.na(input1[,trait]),c(sample_ID, trait)]
  interid <- intersect(colnames(input2), input1[,sample_ID])
  if (length(interid) == 0) {
      message("No id corresponding. Stopping execution.")
      return()
  } else {
      cat(paste0(length(interid), " intersect ids"))
      cat("The following is the interids:")
      cat(interid)
  }
  # interid <- intersect(colnames(input2),as.character(ph2[[sample_ID]]))
  y <- input1[input1[,sample_ID] %in% interid, trait]
  # y <- input1[interid,trait]
  # cat(head(y))
  input2[upper.tri(input2)] <- t(input2)[upper.tri(input2)]
  mat <- input2[rownames(input2) %in% interid, colnames(input2) %in% interid]
  # mat <- input2[interid,interid]
  # cat(dim(mat))
  Z <- c.z.hglm(mat)
  return(list("y"=y,"Z" = Z))
}


# five fold cross validation
require(rrBLUP)
cv_rrblup <- function(y,Z,t_now,dir){
  cv <- sample(1:5,size = length(y),prob = c(0.2,0.2,0.2,0.2,0.2),replace = T)
  require(rrBLUP)
  acc <- c()
  for(i in 1:5){
    #train <- interid[!(cv ==i)]
    #test <- interid[cv ==i]
    # est g first
    est <- mixed.solve(y=y,Z = Z,X = as.matrix(rep(1,length(y))))
    gt <- Z %*% est$u
    
    ## test g
    y2 <- y
    y2[cv ==i] <- NA
    est2 <- mixed.solve(y=y2,Z = Z,X = as.matrix(rep(1,length(y2))))
    gt2 <- Z %*% est2$u
    
    r_now <- cor(gt[cv ==i], gt2[cv ==i])
    
    acc <- c(acc,r_now)
    cat(i,"\n")
  }
  pre1 <- gt[cv ==i]
  pre2 <- gt2[cv ==i]
  save(y,est,pre1,pre2,acc,Z,cv,file = paste0(dir ,t_now,".RData"))
  return(acc)
}

loop_cv<-function(input1, input2, sample_ID, out, acc_result, dir){
  traits <- out$traits
  cat("loop_cv")
  for (i in 1:length(traits)){
    t_now <- traits[i]
    input_now <- prep_input(input1 = input1, input2 = input2, sample_ID = sample_ID, trait = t_now)
    out_now <- cv_rrblup(y = input_now$y,Z = input_now$Z,t_now=t_now, dir = dir)
    r <- mean(out_now,na.rm = T)
    sd <- sd(out_now,na.rm = T)
    out[i,2:3] <- c(r,sd)
    acc_result <- rbind(acc_result, data.frame(trait = t_now, acc1=out_now[1],acc2=out_now[2],acc3=out_now[3],acc4=out_now[4],acc5=out_now[5],mean_acc = r, sd_acc = sd))
  }
  return(acc_result)
}
