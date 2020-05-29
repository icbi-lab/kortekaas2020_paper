# Packages
library(Seurat)
library(sctransform)
library(dplyr)
library(SingleR)
library(ggplot2)

# Create directories for results 
dir.create("../Results", showWarnings = FALSE)
dir.create("../Results/RData", showWarnings = FALSE)
dir.create("./Results/Figures", showWarnings = FALSE)
dir.create("../Results/Tables", showWarnings = FALSE)

# Define Path to directory
datadir <- "../GEX/"
subjects <- list.files(datadir)

sobjs<-vector(mode="list", length=length(subjects))
names(sobjs)<-subjects

for (i in 1:length(subjects)){
  
  cat(i, "/", length(subjects), "\n")
  
  subject<-subjects[i]
  
  cdir<-file.path(datadir, subject, "filtered_feature_bc_matrix")

  seurat_obj<-CreateSeuratObject(counts = Read10X(data.dir=cdir), 
                             project = subject, 
                             min.cells = 3, 
                             min.features = 200)
  
  seurat_obj@meta.data$subject<-subject
  
  sobjs[[subject]]<-seurat_obj

}

# Merge
hncGEX<-merge(sobjs[[1]], 
                   y=c(sobjs[[subjects[2]]], sobjs[[subjects[3]]], sobjs[[subjects[4]]], 
                       sobjs[[subjects[5]]], sobjs[[subjects[6]]], sobjs[[subjects[7]]],
                       sobjs[[subjects[8]]], sobjs[[subjects[9]]], sobjs[[subjects[10]]], 
                       sobjs[[subjects[11]]], sobjs[[subjects[12]]], sobjs[[subjects[13]]]),
                   add.cell.ids=subjects, 
                   project = "HNC_HPV")

# Save
fileout<-"../Results/RData/hncGEX_Seurat_obj.RData"
save(hncGEX, file=fileout)
