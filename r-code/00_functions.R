options(dplyr.print_max = 1e9)
options(max.print = 1e9)

# install packages if not present
using<-function(...) {
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  if(length(need)>0){ 
    install.packages(need)
    lapply(need,require,character.only=TRUE)
  }
}
# usage: using("RCurl","ggplot2","jsonlite","magrittr")

# update and append objects to a .RData file
resave <- function(..., list = character(), file) {
  previous  <- load(file)
  var.names <- c(list, as.character(substitute(list(...)))[-1L])
  for (var in var.names) assign(var, get(var, envir = parent.frame()))
  save(list = unique(c(previous, var.names)), file = file)
}



#function to extract the legend of a ggplot; source:
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs

get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


#
split_tibble <- function(tibble, col = 'col') tibble %>% split(., .[, col]) # function to spilt data by group


# convert factors to numeric
# https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-integer-numeric-without-loss-of-information
as.charnum <- function(x) {as.numeric(levels(x))[x]}

# usage:
# cols = c("temperature")
# dfcols] <- lapply(df[cols], as.charnum)


PredictorSelection<-function(predictors.cor){
  df<-cbind(row.names(predictors.cor),predictors.cor)
  colnames(df)<-c("pred","predcor")
  rownames(df)<-1:nrow(df)
  
  
  
  # SELECTION starts here:
  
  # get rid of factors!, sonst Wanrmeldung beim letzten loop:
  # Warnmeldung:
  # In `[<-.factor`(`*tmp*`, iseq, value = 
  # "rs13,rs265,bgl,fk_klasse,nfk_klasse,kak_klasse,kfa_klasse,natveg,kultpfla,fipu_a_ln,fipu_o_ln,fipu_s_ln") :
  #  invalid factor level, NAs generated
  
  xy<-data.frame(lapply(df,as.character),stringsAsFactors=FALSE)
  
  
  # SELECT variables with 0 correlation
  y<-as.vector(xy[xy[,2]==0,1])
  if(length(which(xy[,2]==0))!=0){xy<-xy[-which(xy[,2]==0),]} # delete selected variables
  
  
  gg<-NULL
  while(nrow(xy)>0){
    cat("\nLeft variables:",nrow(xy),"\n")
    print(as.vector(xy[,1]))
    #print(xy)
    
    cat("\nSelect an additional environmental variable! (Type 'z' to quit) \n")
    x<-readline(prompt ="Input: ")

    

    
    while(!(x %in% xy[,1])){
      if (x == "z") {
      print("Quitting variable selection!");
      break
    } else {
      cat("\nWrong spelling - please try again!\n")
      x<-readline(prompt ="Input: ")
      }
    }
    
    
    # STEP 1: add variable to selection
    y<-c(y,x)
    
    # STEP 2: delete correlated variables from dataframe
    x.cor1<-strsplit(as.vector(xy[which(xy[,1]==x),2]),"\\,") # select corr variables
    x.cor2<-which(xy[,1] %in% x.cor1[1][[1]]) # which rows?
    if(length(x.cor2)!=0){
      xy<-xy[-x.cor2,] # delete selected variables (their rows)
    }
    # STEP 3: delete variable x
    xy<-xy[-which(xy[,1] %in% x),]
    
    
    
    # Step 4: delete correlated variables from the "correlated variables column" of other variables
    if(nrow(xy)>0){
      for (i in 1:nrow(xy)){
        g<-strsplit(as.vector(xy[i,2]),"\\,")[[1]] # splitted string of correlated variables(column 2) belonging to variable i
        xy[i,2]<-paste(g[!(g %in% x.cor1[1][[1]])],collapse=",") # get 
        if(length(g[!(g %in% x.cor1[1][[1]])])==0){ # if variable has 0 correlated variables -> add it to y
          y<-c(y,as.vector(xy[i,1]))
          gg<-c(gg,which(xy[,1] %in% xy[i,1])) # record rownumber of variable to be dropped
        }
        
      }
      if (length(gg)>0){
        xy<-xy[-gg,]
        gg<-NULL}
    }
    
  }
  
  
  
  pred.uncor<-y

  if (x == "z") {
    return()
  } else{
    print(pred.uncor);
    return(pred.uncor);
  } 
   
}



# speciesMin <- function(data,cols,species,genus,family,order,class){
#   
#   m     = data[,cols]>0 # read >0
#   zu    = data.frame(species=species,genus=genus,family=family,order=order,class=class)
#   
#   
#   #species
#   ms     <- m[zu$species!="",]  # reads with identified species
#   as     <- zu[zu$species!="",] # names of identified species
#   
#   #genus
#   ag      = zu[zu$species==""& zu$genus!="" & zu$family!="" & zu$order!="",]     # names of genus
#   mg      = m[zu$species=="" & zu$genus!="" & zu$family!="" & zu$order!="",]     # presence of genus
#   
#   if(is.vector(mg)==T){mg <- t(mg)}
#   
#   if(nrow(ag)>0){
#     for(i in 1:nrow(ag)){
#       mn <- m[zu$genus==ag$genus[i],]
#       if(is.vector(mn)==F){
#         r <- apply(mn,2,sum) # we have species und genus
#         rr <- r==1&mg[i,]==T
#         mg[i,] <- rr
#       }
#       else {
#         mg[i,] <- mn
#       }
#     }
#   }
#   
#   
#   
#   #family
#   af <- zu[zu$species==""&zu$genus==""&zu$family!=""&zu$order!="",]
#   mf <- m[zu$species==""&zu$genus==""&zu$family!=""&zu$order!="",]
#   
#   if(is.vector(mf)==T){mf <- t(mf)}
#   if(nrow(af)>0){
#     for(i in 1:nrow(af)){
#       mn <- m[zu$family==af$family[i],]
#       if(is.vector(mn)==F){
#         r <- apply(mn,2,sum)
#         r <- r==1&&mf[i,]==T
#         mf[i,] <- r
#       }
#       else {
#         mf[i,] <- mn
#       }
#     }
#   }
#   
#   #order
#   ao <- zu[zu$species==""&zu$genus==""&zu$family==""&zu$order!="",]
#   mo <- m[zu$species==""&zu$genus==""&zu$family==""&zu$order!="",]
#   if(is.vector(mo)==T){mo <- t(mo)}
#   
#   if(nrow(ao)>0){
#     for(i in 1:nrow(ao)){
#       mn <- m[zu$order==af$order[i],]
#       if(is.vector(mn)==F){
#         r <- apply(mn,2,sum)
#         r <- r==1&&mo[i,]==T
#         mo[i,] <- r
#       }
#       else {
#         mo[i,] <- mn
#       }
#     }
#   }
#   
#   #class
#   ac <- zu[zu$species==""&zu$genus==""&zu$family==""&zu$order==""&zu$class!="",]
#   mc <- m[zu$species==""&zu$genus==""&zu$family==""&zu$order==""&zu$class!="",]
#   if(is.vector(mc)==T){mc <- t(mc)}
#   
#   if(nrow(ac)>0){
#     for(i in 1:nrow(ac)){
#       mn <- m[zu$class==ac$class[i],]
#       if(is.vector(mn)==F){
#         r <- apply(mn,2,sum)
#         r <- r==1&&mc[i,]==T
#         mc[i,] <- r
#       }
#       else {
#         mc[i,] <- mn
#       }
#     }
#   }
#   
#   sp <- rbind(as,ag,af,ao,ac)
#   co <- rbind(ms,mg,mf,mo,mc)
#   
#   result_occ <- cbind(sp,co)
#   
#   cols <- sapply(result_occ, is.logical)
#   result_occ[,cols] <- lapply(result_occ[,cols], as.numeric)
#   
#   
#   return(result_occ)
# }
# 
# 
# 
# # example usage -presence/absence
# taxon_id = 1:11
# class = c("Ellipura","Arachnida",rep("Insecta",9))
# order = c("Collembola","","Coleoptera","Coleoptera","Coleoptera","Hymenoptera","Hymenoptera","Hemiptera","Hemiptera","Hemiptera","Hemiptera")
# family = c("Sminthuridae","","Carabidae","Carabidae","Carabidae","Ichneumonidae","Apidae","Cicadellidae","Cicadellidae","Cicadellidae","Cicadellidae")
# genus = c("Allacma","", "Carabus","Carabus","Carabus","","","Macrosteles","Macrosteles","Psammotettix","")
# species = c("Allacma fusca","", "Carabus auratus","","Carabus nemoralis","","","Macrosteles sexnotatus","","Psammotettix alienus","")
# site1 = c(0,rep(1,10))
# site2 = c(0,0,2,1,0,0,0,0,0,0,0)
# site3 = c(0,0,2,1,3,0,0,0,0,0,0)
# site4 = c(0,1,0,0,0,1,1,1,0,0,0)
# site5 = c(0,0,0,2,0,0,0,1,1,0,0)
# site6 = c(0,0,0,0,0,0,0,1,0,1,0)
# site7 = c(0,0,0,0,0,0,0,1,0,0,1)
# site8 = c(0,0,0,0,0,0,0,0,1,1,0)
# 
# df = data.frame(taxon_id=taxon_id,class=class,order=order,family = family,genus=genus,species=species,site1=site1,site2=site2,site3=site3,site4=site4,site5=site5,site6=site6,site7=site7)
# 
# q <- speciesMin(data=df,cols=7:ncol(df),species=df$species,genus=df$genus,family=df$family,order=df$order,class=df$class)
# 
# 
# 
# 

# data frame to list based on grouping column
# https://stackoverflow.com/questions/39638233/grouped-data-frame-to-list
split_tibble <- function(tibble, column = 'col') {
  tibble %>% split(., .[,column]) %>% lapply(., function(x) x[,setdiff(names(x),column)])
}




mytheme <-  ggplot2::theme(axis.title = element_text(size =18),
                  axis.text = element_text(size =14),
                  strip.text= element_text(size = 18),
                  plot.title = element_text(hjust = 0.5,size =18),
                  legend.text=element_text(size=14))



mytheme2 <-  ggplot2::theme(axis.title = element_text(size =18),
                   axis.text = element_text(size =10),
                   strip.text= element_text(size = 14),
                   legend.position = "none",
                   plot.title = element_text(hjust = 0.5))

# Custom ggplot theme to make pretty plots
theme_clean <- function() {
  ggplot2::theme_minimal(base_size=18) +
    ggplot2::theme(
      panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = rel(1)),
          strip.background = element_rect(fill = "grey80", color = NA),
          legend.title = element_text(face = "bold"),
      axis.line = element_line())
}




### modified spCodes function
# omits error message
spCodes2 <- function(species, nchar.gen = 3, nchar.sp = 3, nchar.ssp = 0, sep.species = " ", sep.spcode = "") {
  # version 1.8 (7 Mai 2013)
  
  species <- as.character(species)
  splits <- strsplit(species, split = sep.species)
  gen <- sp <- vector("character", length(splits))
  for (i in 1:length(splits)) {
    gen[i] <- splits[[i]][1]
    sp[i] <- splits[[i]][2]
    if (is.na(sp[i]))  sp[i] <- ""
  }  # end for i
  
  abbrev <- function (string, nchar) {
    abv <- substr(string, start = 1, stop = nchar)
    return(abv)
  }  # end abbrev function
  gen.code <- abbrev(gen, nchar.gen)
  sp.code <- abbrev(sp, nchar.sp)
  spcode <- paste(gen.code, sp.code, sep = sep.spcode)
  
  if (nchar.ssp > 0) {
    ssp <- vector("character", length(splits))
    for (i in 1:length(splits)) {
      ssp[i] <- splits[[i]][3]
      if (is.na(ssp[i]))  ssp[i] <- ""
    }
    ssp.code <- abbrev(ssp, nchar.ssp)
    spcode <- paste(spcode, ssp.code, sep = sep.spcode)
  }
  return(spcode)
}






# iNext.link data preparation function
# reguirements
# PD: all species from coocc dataset must be present in phylogenetic trees
# PD: trees in newick format

# FD: all species from coocc must be present in trait dataset and vice versa
# FD: trait table as long format and species in 1. column

#cooccDF=bee.coocc;nTreatments=3;iTraitDF=iTraits;pTraitDF=pTraits;iTree=iTree;pTree=pTree


iNext.link.dataprep = function(cooccDF=NULL,cooccDF_PD =NULL, cooccDF_FD =NULL,
                               nTreatments=1,
                               iTraitDF=NULL,pTraitDF=NULL,
                               iTree=NULL,pTree=NULL){
  
  out <- list()
  
  # data cleaning: remove species and plants with no occurrences
  
  for(i in 1:nTreatments){
    if(length(which(colSums(cooccDF[[i]])==0))>0){
      cooccDF[[i]] = cooccDF[[i]][,-which(colSums(cooccDF[[i]])==0)]}
    
    if(length(which(rowSums(cooccDF[[i]])==0))>0){
      cooccDF[[i]] = cooccDF[[i]][-which(rowSums(cooccDF[[i]])==0),]}
  } 
  
  # rename species names
  for(i in 1:nTreatments){
    colnames(cooccDF[[i]]) = sub(" ", "_", colnames(cooccDF[[i]]))
    colnames(cooccDF[[i]]) = sub("\\.", "",colnames(cooccDF[[i]]))
    
    rownames(cooccDF[[i]]) = sub(" ", "_", rownames(cooccDF[[i]]))
    rownames(cooccDF[[i]]) = sub("\\.", "",rownames(cooccDF[[i]]))
  } 
  
  ### ***  TD  *** ###
  out$cooccTD = cooccDF
  
  cat("Information on TD input data\n") 
  x1 = length(unique(unlist(lapply(out$cooccTD, function(x) colnames(x)))))
  x2 = length(unique(unlist(lapply(out$cooccTD, function(x) rownames(x)))))
  
  cat(x1, " unique insect species and ",x2," unique plant species\n\n\n") 

  
  ### ***  PD  *** ###
  if(!is.null(iTree) & !is.null(pTree)){
    if(!is.null(cooccDF_PD)){coocc = cooccDF_PD}else{coocc = cooccDF}
    
    iTree$node.label =  paste("match",rep(1:length(iTree$node.label)))
    pTree$node.label =  paste("match",rep(1:length(pTree$node.label)))
    
    # remove species in treatment datasets which are not in trees
    for(i in 1:nTreatments){
      # insect in colnames
      if(length(which(!colnames(coocc[[i]]) %in% iTree$tip.label ))>0){
        coocc[[i]] = coocc[[i]][,-which(!colnames(coocc[[i]]) %in% iTree$tip.label)]}
      
      # plants in rows
      if(length(which(!rownames(coocc[[i]]) %in% pTree$tip.label))>0){
        coocc[[i]] = coocc[[i]][-which(!rownames(coocc[[i]]) %in% pTree$tip.label),]}
    } # end remove species in treatments of coocc dataset

    for(i in 1:nTreatments){
      if(length(which(colSums(coocc[[i]])==0))>0){
        coocc[[i]] = coocc[[i]][,-which(colSums(coocc[[i]])==0)]}

      if(length(which(rowSums(cooccDF[[i]])==0))>0){
        coocc[[i]] = coocc[[i]][-which(rowSums(coocc[[i]])==0),]}
    }
    

    out$cooccPD = coocc
    out$iTree = iTree
    out$pTree = pTree
    
    
  }  # close PD

  cat("Information on PD input data\n") 
  x1 = length(which(unique(unlist(lapply(coocc, function(x) colnames(x)))) %in% iTree$tip.label))
  x2 = length(unique(unlist(lapply(coocc, function(x) colnames(x)))))      
  
  cat("Are all insect species of coocc dataset in insect tree: ", identical(x1,x2),"\n")
  
  x1 = length(which(unique(unlist(lapply(coocc, function(x) rownames(x)))) %in% pTree$tip.label))
  x2 = length(unique(unlist(lapply(coocc, function(x) rownames(x)))))    
  cat("Are all plant species of coocc dataset in plant tree: ", identical(x1,x2),"\n\n\n")
  
  
  ### ***  FD  *** ###
  if(!is.null(iTraitDF) & !is.null(pTraitDF)){
    coocc = cooccDF
    
    # remove species with less than 2 traits
    iTraitNA = which(rowSums(!is.na(iTraitDF[,-1]))<2)
    if(length(iTraitNA)>0){iDF = iTraitDF[-iTraitNA,]}else{iDF = iTraitDF}
    
    pTraitNA = which(rowSums(!is.na(pTraitDF[,-1]))<2)
    if(length(pTraitNA)>0){pDF = pTraitDF[-pTraitNA,]}else{pDF = pTraitDF}
    
    
    iDF[,1] = sub(" ", "_", iDF[,1])
    iDF[,1] = sub("\\.", "",iDF[,1])
    
    pDF[,1] = sub(" ", "_", pDF[,1])
    pDF[,1] = sub("\\.", "",pDF[,1])
    
    
    for(i in 1:nTreatments){
      # insect in colnames
      if(length(which(!colnames(coocc[[i]]) %in% iDF[,1]))>0){
        coocc[[i]] = coocc[[i]][,-which(!colnames(coocc[[i]]) %in% iDF[,1])]}
      
      # plants in rows
      if(length(which(!rownames(coocc[[i]]) %in% pDF[,1]))>0){
        coocc[[i]] = coocc[[i]][-which(!rownames(coocc[[i]]) %in% pDF[,1]),]}
    }
    
    
    # remove species in traits dataset that are not in interaction dataset
i=0
    while(any(unlist(lapply(coocc, function(x) rowSums(x)))==0) || # any plant with no interaction
          any(unlist(lapply(coocc, function(x) colSums(x)))==0) ){ # any insect with no interaction
i=i+1
print(i)

      if(any(unlist(lapply(coocc, function(x) rowSums(x)))==0)){
      for(i in 1:nTreatments){
        # plants with no interactions
        if(any(rowSums(coocc[[i]])==0)){
          coocc[[i]] = coocc[[i]][-which(rowSums(coocc[[i]])==0),]}
      } # end remove species in treatments of coocc dataset
      }
      
      if(any(unlist(lapply(coocc, function(x) colSums(x)))==0)){
        for(i in 1:nTreatments){
          # insects with no interactions
          if(any(colSums(coocc[[i]])==0)){
            coocc[[i]] = coocc[[i]][,-which(colSums(coocc[[i]])==0)]}
        } # end remove species in treatments of coocc dataset
      }
    }



      
    if(length(which(!iDF[,1] %in% sort(unique(unlist(lapply(coocc, function(x) colnames(x)))))))>0){
      iDF = iDF[-which(!iDF[,1] %in% sort(unique(unlist(lapply(coocc, function(x) colnames(x)))))),]} #  insect traits rows in col 1 - cooc colnames
    
    if(length(which(!pDF[,1] %in% sort(unique(unlist(lapply(coocc, function(x) rownames(x)))))))>0){
      pDF = pDF[-which(!pDF[,1] %in% sort(unique(unlist(lapply(coocc, function(x) rownames(x)))))),]} # plants traits rows in col 1 - cooc rownames
    

    
    # calculation of distance matrices for traits
    iTraits.gwd = StatMatch::gower.dist(iDF[,-1])
    colnames(iTraits.gwd) =  iDF[,1]
    rownames(iTraits.gwd) =  iDF[,1]
    
    pTraits.gwd = StatMatch::gower.dist(pDF[,-1])
    colnames(pTraits.gwd) =  pDF[,1]
    rownames(pTraits.gwd) =  pDF[,1]
    
    out$cooccFD = coocc
    out$iTraits.gwd = iTraits.gwd
    out$pTraits.gwd = pTraits.gwd
    
    cat("Information on FD input data\n")  
    cat("Number of insect species in coocc dataset: ", length(unique(unlist(lapply(out$cooccFD, function(x) colnames(x))))),"\n")
    cat("Number of insect species in trait matrix: ", dim(out$iTraits.gwd)[1],"\n")
    
    
    x1 = sort(unique(unlist(lapply(out$cooccFD, function(x) colnames(x)))))
    x2 = sort(colnames(out$iTraits.gwd))
    cat("Are insect species in coocc and trait matrix identical? ", identical(x1,x2),"\n\n")
    
    
    cat("Number of plant species in coocc dataset: ", length(unique(unlist(lapply(out$cooccFD, function(x) rownames(x))))),"\n")
    cat("Number of plant species in trait matrix: ", dim(out$pTraits.gwd)[1],"\n")
    
    x1 = sort(unique(unlist(lapply(out$cooccFD, function(x) rownames(x)))))
    x2 = rownames(out$pTraits.gwd)
    cat("Are plant species in coocc and trait matrix identical? ", identical(x1,x2),"\n")

    
  } # close FD
  
  return(out)
}







# TAPNET

tapnet.dataprep = function(cooccDF=NULL,iTraitDF=NULL,pTraitDF=NULL,iTree=NULL,pTree=NULL,nTraits = 1){
  
  out <- list()
  
  
  
  # data cleaning: remove species and plants with no occurrences
  
  
  if(length(which(colSums(cooccDF)==0))>0){
    cooccDF = cooccDF[,-which(colSums(cooccDF)==0)]}
  
  if(length(which(rowSums(cooccDF)==0))>0){
    cooccDF = cooccDF[-which(rowSums(cooccDF)==0),]}
  
  # rename species names
  colnames(cooccDF) = sub(" ", "_", colnames(cooccDF))
  colnames(cooccDF) = sub("\\.", "",colnames(cooccDF))
  
  rownames(cooccDF) = sub(" ", "_", rownames(cooccDF))
  rownames(cooccDF) = sub("\\.", "",rownames(cooccDF))
  
  
  coocc = cooccDF
  
  
  
  ### ***  PD  *** ###
  
  iTree$node.label =  paste("match",rep(1:length(iTree$node.label)))
  pTree$node.label =  paste("match",rep(1:length(pTree$node.label)))
  
  
  # remove species in treatment datasets which are not in trees
  
  # insect in colnames
  if(length(which(!colnames(coocc) %in% iTree$tip.label ))>0){
    coocc = coocc[,-which(!colnames(coocc) %in% iTree$tip.label)]}
  
  # plants in rows
  if(length(which(!rownames(coocc) %in% pTree$tip.label))>0){
    coocc = coocc[-which(!rownames(coocc) %in% pTree$tip.label),]}
  
  
  out$iTree = iTree
  out$pTree = pTree
  
  
  
  
  
  ### ***  FD  *** ###
  if(!is.null(iTraitDF) & !is.null(pTraitDF)){
    
    # remove species with less than 2 traits
    if(nTraits >1){
    iTraitNA = which(rowSums(!is.na(iTraitDF[,-1]))<1)
    if(length(iTraitNA>0)){iDF = iTraitDF[-iTraitNA,]}else{iDF = iTraitDF}
    
    pTraitNA = which(rowSums(!is.na(pTraitDF[,-1]))<1)
    if(length(pTraitNA>0)){pDF = pTraitDF[-pTraitNA,]}else{pDF = pTraitDF}
    }
    
    if(nTraits ==1){
      iTraitNA = which(is.na(iTraitDF[,-1]))
      if(length(iTraitNA>0)){iDF = iTraitDF[-iTraitNA,]}else{iDF = iTraitDF}

     pTraitNA = which(is.na(pTraitDF[,-1]))
     if(length(pTraitNA>0)){pDF = pTraitDF[-pTraitNA,]}else{pDF = pTraitDF}
    }
    
    iDF[,1] = sub(" ", "_", iDF[,1])
    iDF[,1] = sub("\\.", "",iDF[,1])
    
    pDF[,1] = sub(" ", "_", pDF[,1])
    pDF[,1] = sub("\\.", "",pDF[,1])
    
    # remove species in traits dataset that are not in interaction dataset
    round1 = 1 # go through the while fct the 1. time
    while(any(rowSums(coocc)==0) || round1 == 1){
      round1 = 2
      
      
      # plants with no interactions
      if(any(rowSums(coocc)==0)){
        coocc = coocc[-which(rowSums(coocc)==0),]}
      
      
      
      if(length(which(!iDF[,1] %in% colnames(coocc)))>0){
        iDF = iDF[-which(!iDF[,1] %in% colnames(coocc)),]} #  insect traits rows in col 1 - cooc colnames
      
      if(length(which(!pDF[,1] %in% rownames(coocc)))>0){
        pDF = pDF[-which(!pDF[,1] %in% rownames(coocc)),]} # plants traits rows in col 1 - cooc rownames
      
      # remove species in treatment datasets which are not in traits datasets
      
      # insect in colnames
      if(length(which(!colnames(coocc) %in% iDF[,1]))>0){
        coocc = coocc[,-which(!colnames(coocc) %in% iDF[,1])]}
      
      # plants in rows
      if(length(which(!rownames(coocc) %in% pDF[,1]))>0){
        coocc = coocc[-which(!rownames(coocc) %in% pDF[,1]),]}
      
    } # end remove species in treatments of coocc dataset
    
    
    # calculation of distance matrices for traits
    iTraits.gwd = StatMatch::gower.dist(iDF[,-1])
    colnames(iTraits.gwd) =  iDF[,1]
    rownames(iTraits.gwd) =  iDF[,1]
    
    pTraits.gwd = StatMatch::gower.dist(pDF[,-1])
    colnames(pTraits.gwd) =  pDF[,1]
    rownames(pTraits.gwd) =  pDF[,1]
    
    iMat           = as.matrix(iDF[,-1])
    rownames(iMat) = iDF[,1]
    
    pMat              = as.matrix(pDF[,-1])
    rownames(pMat) = pDF[,1]
    
    
    out$coocc = coocc
    out$iTraits = iMat
    out$pTraits = pMat
    out$iTraits.gwd = iTraits.gwd
    out$pTraits.gwd = pTraits.gwd
    
    
    
    
  } # close FD
  
  return(out)
}



generateWebs = function(Nnetworks =100, seed=2023, maxNspecies=10){
  # minimum number of species is set to 5
  set.seed(NULL); set.seed = seed
  networks = list()
  
  for (i in 1:Nnetworks){
    
    # random number for number of species
    Nrow = sample(5:maxNspecies,1)
    Ncol = sample(5:maxNspecies,1)
    Ncell = Nrow*Ncol
    
    # generate abundances and zeros with prob 0.7
    dat_abu = matrix(sample(0:10,Ncell,replace=TRUE),ncol=Ncol)
    dat_bin = matrix(rbinom(n=Ncell,size=1,prob=0.2),ncol=Ncol)
    
    dat_abu[which(dat_bin==0)] <- 0
    
    colnames(dat_abu) = LETTERS[sort(sample(1:(maxNspecies+2), Ncol))]
    rownames(dat_abu) = letters[sort(sample(1:(maxNspecies+2), Nrow))]
    
    
    dat_abu = data.frame(dat_abu)
    web = data.frame(network=i,species2 = row.names(dat_abu),dat_abu)
    web_long=tidyr::pivot_longer(web,  cols= -c(network,species2),names_to = "species1",
                                 values_to = "value", values_drop_na = TRUE)
    
    
    networks[[i]] = web_long
  }
  
  networkList = networks %>% bind_rows
  
  # Convert dataframe to 3D array
  webs_array <- bipartite::frame2webs(networkList, 
                                      varnames=c("species1", "species2","network","value"), 
                                      type.out="array")
  
  return(webs_array)

}


# webs_array = generateWebs()


rarebeta = function(webs_array,iter=50, seed=2023){
  set.seed(NULL); set.seed = seed
  Nnetworks = dim(webs_array)[3]
  result_samplingEffort = list()
  result_iter           = list()  
  
  metaweb <- as.array(rowSums(webs_array, dims=2), dim=3, 
                      dimnames=list("species1", "species2", "network"))
  
  for(i in c(1:iter)) {
    for(Nw in c(1:Nnetworks)){
      Nets = sample(c(1:Nnetworks),Nw,replace=FALSE)
      if(length(Nets) == 1){local_metaweb = webs_array[,,Nets]}
      if(length(Nets) > 1){
        local_metaweb <- as.array(rowSums(webs_array[,,Nets], dims=2), dim=3, 
                                  dimnames=list("species1", "species2", "network"))}
      network_array = abind::abind(webs_array,metaweb,local_metaweb)
      dimnames(network_array)[[3]][Nnetworks+1] = "metaweb"
      dimnames(network_array)[[3]][Nnetworks+2] = "local_metaweb"
      
      #NuniqueSpecies = sum(dim(local_metaweb))
      #NuniqueLinks   = sum(local_metaweb)
      #Connectance    = sum(local_metaweb)/prod(dim(local_metaweb))
      #Specificity    = mean(ESM::getspe(local_metaweb))
      #Nestedness     = vegan::nestednodf(local_metaweb)$statistic[3]/100
      
      
      result_samplingEffort[[Nw]] = c(Nets = length(Nets), 
                                      iter = i,
                                      bipartite::networklevel(local_metaweb,index=c("number of species","connectance","NODF")),
                                      NuniqueLinks   = sum(local_metaweb),
                                      bipartite::betalinkr(network_array[,,c("metaweb","local_metaweb")], partitioning="poisot"))
    }
    result_iter[[i]] = result_samplingEffort
  }
  
  result = result_iter %>% bind_rows
  return(result)
}

# result = rarebeta(webs_array)



# https://stackoverflow.com/questions/5620885/how-does-one-reorder-columns-in-a-data-frame
##arrange df vars by position
##'vars' must be a named vector, e.g. c("var.name"=1)
arrange.vars <- function(data, vars){
  ##stop if not a data.frame (but should work for matrices as well)
  stopifnot(is.data.frame(data))
  
  ##sort out inputs
  data.nms <- names(data)
  var.nr <- length(data.nms)
  var.nms <- names(vars)
  var.pos <- vars
  ##sanity checks
  stopifnot( !any(duplicated(var.nms)), 
             !any(duplicated(var.pos)) )
  stopifnot( is.character(var.nms), 
             is.numeric(var.pos) )
  stopifnot( all(var.nms %in% data.nms) )
  stopifnot( all(var.pos > 0), 
             all(var.pos <= var.nr) )
  
  ##prepare output
  out.vec <- character(var.nr)
  out.vec[var.pos] <- var.nms
  out.vec[-var.pos] <- data.nms[ !(data.nms %in% var.nms) ]
  stopifnot( length(out.vec)==var.nr )
  
  ##re-arrange vars by position
  data <- data[ , out.vec]
  return(data)
}



# round to 3 significant digit
my_signif <- function(x, digits = 3) {
  
  round(x, pmax(digits - ceiling(log10(abs(x))), 0))
  
} 




my_upper <- function(data, mapping, ...) {
  
  # get the x and y data to use the other code
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  r  <- unname(cor.test(x, y)$estimate)
  rt <- format(r, digits = 2)[1]
  tt <- as.character(rt)
  
  # plot the cor value
  ggally_text(
    label = tt, 
    mapping = aes(),
    size = 4,
    color = "grey20") +
    theme_void() 
}

my_diag <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) + 
    geom_density(fill = "grey50") +
    theme_ipsum()
}

my_lower <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) + 
    geom_point(shape = 1, size = 1/2, alpha = 1/6) +
    theme_ipsum()
}



# https://stackoverflow.com/questions/5620885/how-does-one-reorder-columns-in-a-data-frame

##arrange df vars by position
##'vars' must be a named vector, e.g. c("var.name"=1)
arrange.vars <- function(data, vars){
  ##stop if not a data.frame (but should work for matrices as well)
  stopifnot(is.data.frame(data))
  
  ##sort out inputs
  data.nms <- names(data)
  var.nr <- length(data.nms)
  var.nms <- names(vars)
  var.pos <- vars
  ##sanity checks
  stopifnot( !any(duplicated(var.nms)), 
             !any(duplicated(var.pos)) )
  stopifnot( is.character(var.nms), 
             is.numeric(var.pos) )
  stopifnot( all(var.nms %in% data.nms) )
  stopifnot( all(var.pos > 0), 
             all(var.pos <= var.nr) )
  
  ##prepare output
  out.vec <- character(var.nr)
  out.vec[var.pos] <- var.nms
  out.vec[-var.pos] <- data.nms[ !(data.nms %in% var.nms) ]
  stopifnot( length(out.vec)==var.nr )
  
  ##re-arrange vars by position
  data <- data[ , out.vec]
  return(data)
}

# usage: arrange.vars(table, c("Out"=2))