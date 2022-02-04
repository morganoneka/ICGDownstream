library(ggplot2)
library(ggalluvial)
library(reshape2)
library(stringr)
library(colorspace)

# function for reading in a CSV file containing the Giotto ICG output
read_in_data <-function(fname){
  input_data <- read.table(fname, header=TRUE, sep=",", fill=TRUE, stringsAsFactors=FALSE)
  input_data[which(input_data$CPGscores.p.adj <= 0.05),]
}

# creates sankey plot
create_sankey <- function(input_df){
  input_df = input_df[,c("CPGscores.cell_type", "CPGscores.int_cell_type")] %>% table() %>% as.data.frame()
  ggplot(input_df, aes(y=Freq, axis1=CPGscores.cell_type, axis2=CPGscores.int_cell_type)) + 
    geom_alluvium(aes(fill=CPGscores.int_cell_type), width=1/12) +
    geom_stratum(width = 1/12, fill = "black", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("Cell Type", "Interacting Cell"), expand = c(.05, .05)) 
}

# identifies the n genes that are most perturbed by cell interactions
get_top_genes <- function(df, n=10){
  df$CPGscores.genes %>% table() %>% sort(decreasing=TRUE) %>% head(n)
}

# identifies genes that are under-expressed in cells of type ct1 with type ct2 neighbors
underexpressed_genes_for_inx <- function(df, ct1, ct2){
  brain_sub = df[which(df$CPGscores.cell_type == ct1 & df$CPGscores.int_cell_type == ct2), c("CPGscores.genes","CPGscores.diff") ] 
  brain_sub[which(brain_sub$CPGscores.diff<0),] 
}

# identifies genes that are under-expressed in cells of type ct1 with type ct2 neighbors
overexpressed_genes_for_inx <- function(df, ct1, ct2){
  brain_sub = df[which(df$CPGscores.cell_type == ct1 & df$CPGscores.int_cell_type == ct2), c("CPGscores.genes","CPGscores.diff") ] 
  brain_sub[which(brain_sub$CPGscores.diff>0),] 
}

# compiles a df of underexpressed genes for all interactions
all_underexpressed_genes <- function(df){
  a = lapply(unique(df$CPGscores.cell_type), function(x){
    ug = lapply(unique(df$CPGscores.int_cell_type), function(y){
      tmp=underexpressed_genes_for_inx(df,x,y)
      
      if(nrow(tmp) >0){
        tmp$Cell1 = x
        tmp$Cell2=y
        return(tmp)
      } else{
        return(NA)
      }
      
    }  )
    
    # print(ug)
    do.call(rbind, ug)
    
    
  })
  
  remove_from_a = lapply(1:length(a), function(x){
    if (ncol(a[[x]]) == 1){
      return (x)
    }
  }) %>% unlist()
  a[remove_from_a] <- NULL
  
  return(do.call(rbind, a))  
}

# compiles a df of overexpressed genes for all interactions
all_overexpressed_genes <- function(df){
  a = lapply(unique(df$CPGscores.cell_type), function(x){
    ug = lapply(unique(df$CPGscores.int_cell_type), function(y){
      tmp=overexpressed_genes_for_inx(df,x,y)
      
      if(nrow(tmp) >0){
        tmp$Cell1 = x
        tmp$Cell2=y
        return(tmp)
      } else{
        return(NA)
      }
      
    }  )
    
    # print(ug)
    do.call(rbind, ug)
    
    
  })
  
  remove_from_a = lapply(1:length(a), function(x){
    if (ncol(a[[x]]) == 1){
      return (x)
    }
  }) %>% unlist()
  a[remove_from_a] <- NULL
  
  return(do.call(rbind, a))  
}

# identify all interactions that influence the expression of a given gene 
subset_by_gene <- function(df, gene){
  df[which(brain_1a$CPGscores.genes == gene), c("CPGscores.unif_int", "CPGscores.diff")]
}

# compare and contrast differentially expressed genes in an interaction differ between two samples
compare_diff_exp <- function(sample1, sample2, ct1, ct2){
  o_1 = all_overexpressed_genes(sample1)
  u_1 = all_underexpressed_genes(sample1)
  
  o_1_spef = o_1[which(o_1$Cell1==ct1 & o_1$Cell2==ct2),] 
  u_1_spef = u_1[which(u_1$Cell1==ct1 & u_1$Cell2==ct2),]
  
  o_2 = all_overexpressed_genes(sample2)
  u_2 = all_underexpressed_genes(sample2)
  
  o_2_spef = o_2[which(o_2$Cell1==ct1 & o_2$Cell2==ct2),]
  u_2_spef = u_2[which(u_2$Cell1==ct1 & u_2$Cell2==ct2),]
  
  o1_o2 = intersect(o_1_spef$CPGscores.genes, o_2_spef$CPGscores.genes)  %>% paste(sep = " ", collapse=", ")
  o1_u2 = intersect(o_1_spef$CPGscores.genes, u_2_spef$CPGscores.genes)  %>% paste(sep = " ", collapse=", ")
  o1_n2 = setdiff(o_1_spef$CPGscores.genes, c(o_2_spef$CPGscores.genes, u_2_spef$CPGscores.genes))  %>% paste(sep = " ", collapse=", ") 
  
  u1_o2 = intersect(u_1_spef$CPGscores.genes, o_2_spef$CPGscores.genes)  %>% paste(sep = " ", collapse=", ")
  u1_u2 = intersect(u_1_spef$CPGscores.genes, u_2_spef$CPGscores.genes)  %>% paste(sep = " ", collapse=", ")
  u1_n2 = setdiff(u_1_spef$CPGscores.genes, c(o_2_spef$CPGscores.genes, u_2_spef$CPGscores.genes))  %>% paste(sep = " ", collapse=", ")
  
  n1_o2 = setdiff(o_2_spef$CPGscores.genes, c(o_1_spef$CPGscores.genes, u_1_spef$CPGscores.genes))  %>% paste(sep = " ", collapse=", ")
  n1_u2 = setdiff(u_2_spef$CPGscores.genes, c(o_1_spef$CPGscores.genes, u_1_spef$CPGscores.genes))  %>% paste(sep = " ", collapse=", ")
  n1_n2 = ""
  
  l = list(
    o1_o2, o1_u2, o1_n2,
    u1_o2, u1_u2, u1_n2,
    n1_o2, n1_u2, n1_n2
  )
  
  names(l) = c("O1_O2", "O1_U2", "O1_N2",
               "U1_O2", "U1_U2", "U1_N2",
               "N1_O2", "N1_U2", "N1_N2")
  
  return(l)
}

# heatmap for the output of compare_diff_exp
plot_comparison <- function(l){
  l_matrix = matrix(unlist(l), nrow=3, ncol=3)
  colnames(l_matrix) = c("Overexpressed in 1st", "Underexpressed in 1st", "Neither in 1st")
  rownames(l_matrix) = c("Overexpressed in 2nd", "Underexpressed in 2nd", "Neither in 2nd")
  
  l_df = l_matrix %>% melt()
  l_df$num = lapply(l_df$value, function(x){
    strsplit(x, ",") %>% unlist %>% length
  }) %>% unlist()
  
  ggplot(l_df, aes(Var1, Var2, fill = num, label=value)) + 
    geom_tile() +
    scale_fill_continuous_sequential(palette = "PinkYl") +
    theme(panel.grid.major.x=element_blank(), #no gridlines
          panel.grid.minor.x=element_blank(), 
          panel.grid.major.y=element_blank(), 
          panel.grid.minor.y=element_blank(),
          panel.background=element_rect(fill="white"), # background=white
          axis.text.x = element_text(angle=45, hjust = 1,vjust=1,size = 12,face = "bold"),
          plot.title = element_text(size=20,face="bold"),
          axis.text.y = element_text(size = 12,face = "bold")) + 
    guides(fill=guide_colorbar(title="number of genes"), label="none") +
    scale_x_discrete(name="") +
    scale_y_discrete(name="")  +  coord_fixed() +  geom_fit_text(reflow = TRUE)
}