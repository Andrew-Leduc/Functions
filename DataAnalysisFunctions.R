## Andrew Leduc Function Files 

### Functions for Max Quant proteomic Data analysis
library(dplyr)
library(reshape2)
library(tibble)
library(ComplexHeatmap)

## Filtering 

getDat <- function(evidence){
  evidence <- evidence %>% filter(PEP < .02, Potential.contaminant != '+',Reverse != '+',PIF > .85)
  evidence <- evidence %>% select('Sequence','Intensity','Retention.time','Leading.razor.protein',
                                 'Missed.cleavages',contains('Reporter.intensity.corrected'))
  
  return(evidence)
}
getRI <- function(evidence){
  b <- NULL
  for (i in 1:ncol(evidence)){
     if (median(evidence[,i]==0)){
        a <- colnames(evidence)[i]
        b <- c(b,a)
    }
  }
  evidence <- evidence %>% select(-all_of(b))
  return(evidence)
}
p2p <- function(data,keep_genes){
  if (keep_genes == TRUE){
    genes <- data %>% select('Leading.razor.protein','Gene.names')
    genes$Gene.names <- as.character(genes$Gene.names)
    genes <- genes[order(nchar(genes$Gene.names)),]
    genes <- genes %>% distinct(Leading.razor.protein,.keep_all = TRUE)
  }
 
  data <- data %>% select('Leading.razor.protein', contains('Reporter'))
  data.m <- melt(data,id.vars = 'Leading.razor.protein')
  data.m <- data.m %>% group_by(Leading.razor.protein,variable) %>% summarise(medRI = median(value,na.rm=TRUE))
  data <- dcast(data.m, Leading.razor.protein ~ variable, value.var = 'medRI')
  
  if (keep_genes == TRUE){
    data <- data[order(data$Leading.razor.protein),]
    genes <- genes[order(genes$Leading.razor.protein),]
    data <- add_column(data,genes$Gene.names,.after = 1)
  }
  return(data)
}

## Normalization
row_norm <- function(data,start,norm){
  a <- colnames(data)[1]
  data_cut <- as.matrix(data[,start:ncol(data)])
  
  if(norm == 0 | missing(norm)){
    for (i in 1:nrow(data_cut)){
      data_cut[i,] <- data_cut[i,]/mean(data_cut[i,],na.rm = TRUE)
    }
  } else {
    for (i in 1:nrow(data_cut)){
      data_cut[i,] <- data_cut[i,]-mean(data_cut[i,],na.rm = TRUE)
    }
  }
  data_cut <- as.data.frame(data_cut)
  data<- cbind(data[,1:(start-1)],data_cut)
  if (start == 2){
    colnames(data)[1] <- a
  }
  return(data)
}
col_norm <- function(data,start,norm){
  if(missing(norm) | norm == 0){
    for (i in start:ncol(data)){
      data[,i] <- data[,i]/median(data[,i],na.rm = TRUE)
    }
  } else {
    for (i in start:ncol(data)){
      data[,i] <- data[,i]-median(data[,i],na.rm = TRUE)
    }
  }
  return(data)
}
row_norm_ref <- function(data,ref,start){
  data_cut <- as.matrix(data[,start:ncol(data)])
  for (i in 1:nrow(data_cut)){
    data_cut[i,] <- data_cut[i,]/data_cut[i,(ref-1)]
  }
  data_cut <- as.data.frame(data_cut)
  data<- cbind(data[,1:(start-1)],data_cut)
  if (start == 2){
    colnames(data)[1] <- "Proteins"
  }
  return(data)
}
naZeros <- function(data){
  data_cut <- as.matrix(data)
  data_cut[data_cut == 0] = NA
  data <- as.data.frame(data_cut)
  return(data)
}
norm_FC <- function(data,norm,start){
  data_cut <- as.matrix(data[,start:ncol(data)])
  for (i in 1:nrow(data)){
    data_cut[i,] <- data_cut[i,]/mean(data_cut[i,],na.rm = TRUE)
  }
  data_cut <- as.data.frame(data_cut)
  data<- cbind(data[,1:(start-1)],data_cut)
}
f2num <- function(data,start,stop){
  
  for (i in start:stop){
    data[,i] <- as.numeric(as.character(data[,i]))
  }
  
  return(data)
}
f2char <- function(data,start,stop){
  data <- as.matrix(data)
  for (i in start:stop){
    data[,i] <- as.character(data[,i])
  }
  data <- as.data.frame(data)
  return(data)
}

## Plotting

Heater <- function(data,title){
  x <- cor(data,use="pairwise.complete.obs")
  h <- Heatmap(x, column_title = title,
               row_names_gp = gpar(fontsize = 7),row_order = rownames(x),column_order = colnames(x),show_row_dend = FALSE,name = "Correlation")
  ht = draw(h)
  
}

## Math Functions

# PLS

# TLS 0 Centered Data



# PCA
PCA <- function(data){
  Q <- cor(data, use="pairwise.complete.obs")
  a<-svd(Q)
  eig_val <- a$d
  eig_vect <- as.data.frame(a$v)
  tot <- sum(eig_val)
  for (i in 1:ncol(eig_vect)){
    colnames(eig_vect)[i] <- paste0("PC", i ,', ',round(eig_val[i]/tot, digits=2),' % of varience')
  }
 
  
  return(eig_vect)
  
}






