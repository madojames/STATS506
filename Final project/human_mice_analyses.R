# Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)
library(ComplexHeatmap)

# load series and platform data from GEO human data
gset <- getGEO("GSE18565", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "000111"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("low","hi"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)

# 提取显著性基因表
tT <- topTable(fit2, adjust="none", sort.by="logFC", number=Inf)  # 提取所有基因
# 设定阈值，筛选差异表达基因
deg <- tT[tT$adj.P.Val < 0.05, ]  # 筛选调整后的 p 值小于 0.05 的基因
deg <- deg[abs(deg$logFC) > 2, ]  # 筛选调整后的 p 值小于 0.05 的基因
n_deg <- nrow(deg)  # 差异表达基因的数量
print(n_deg)

human_gene_deg <- unique(deg$Gene.symbol)
length(human_gene_deg)

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
# Please note that the code provided to generate graphs serves as a guidance to
# the users. It does not replicate the exact GEO2R web display due to multitude
# of graphical options.
# 
# The following will produce basic volcano plot using limma function:
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))



# load series and platform data from GEO mice data

gset1 <- getGEO("GSE17256", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset1) > 1) idx <- grep("GPL1261", attr(gset1, "names")) else idx <- 1
gset1 <- gset1[[idx]]

# make proper column names to match toptable 
fvarLabels(gset1) <- make.names(fvarLabels(gset1))

# group membership for all samples
gsms <- "11110000"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset1)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset1) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("low","hi"))
levels(gs) <- groups
gset1$group <- gs
design <- model.matrix(~group + 0, gset1)
colnames(design) <- levels(gs)

gset1 <- gset1[complete.cases(exprs(gset1)), ] # skip missing values

fit <- lmFit(gset1, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)


# getting differential expression genes
tT <- topTable(fit2, adjust="none", sort.by="logFC", number=Inf) 
deg <- tT[tT$adj.P.Val < 0.05, ]  # p<0.05
deg <- deg[abs(deg$logFC) > 2, ]  # logFC>2
n_deg <- nrow(deg)  
print(n_deg)

mice_gene_deg <- unique(deg$Gene.symbol)
length(mice_gene_deg)

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
# Please note that the code provided to generate graphs serves as a guidance to
# the users. It does not replicate the exact GEO2R web display due to multitude
# of graphical options.
# 
# The following will produce basic volcano plot using limma function:
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

#the union of two set of genes is
intersect(tolower(mice_gene_deg), tolower(human_gene_deg))



#Gene markers identified by the paper (only part of 130 genes)
human_genes <- c("F13A1", "CAPG", "VCAN", "HPR", "F5", "ETHE1", "CCR2", "CLEC5A", "PADI4", "BST1", "TFEC", "HPSE", "APRT", "TGFB1","FCGR1B")
mice_genes <- c("F13a1", "Capg", "Vcan", "Hp", "F5", "Ethe1", "Ccr2", "Clec5a", "Padi4", "Bst1", "Tfec", "Hpse", "Aprt", "Tgfb1","Fcgr1")

#get all gene symbol list
human_gene_name <- gset@featureData@data$Gene.symbol
mice_gene_name <- gset1@featureData@data$Gene.symbol


#which(mice_gene_name %in% "Tgfb1")

#find the corresponding id
position1 <- which(human_gene_name %in% human_genes)
position2 <- which(mice_gene_name %in% mice_genes)

#查看找到那些基因
#在这个地方会存在多个基因名的情况，这是由于微阵列探针本身的问题
#可以建议按照后满的图取想要的行作为展示
vec_hum <- gset@featureData@data$Gene.symbol[position1]
vec_mic <- gset1@featureData@data$Gene.symbol[position2]

#remove the none unique ones
position1 <- position1[match(unique(vec_hum), vec_hum)]
position2 <- position2[match(unique(vec_mic), vec_mic)]

vec_hum <- gset@featureData@data$Gene.symbol[position1]
vec_mic <- gset1@featureData@data$Gene.symbol[position2]

position1 <- position1[order(vec_hum)]
position2 <- position2[order(vec_mic)]

vec_hum <- gset@featureData@data$Gene.symbol[position1]
vec_mic <- gset1@featureData@data$Gene.symbol[position2]

#get the probs ID
ID1 <- gset@featureData@data$ID[position1]
ID2 <- gset1@featureData@data$ID[position2]

# now for human data
temp <- gset@assayData$exprs[ID1,]

#normalize the data by row
temp <- t(scale(t(temp)))


#plot the differential expression genes in human 
#change the rowname from probs id to gene names
row.names(temp) <- vec_hum
Heatmap(temp, show_row_names = TRUE, show_column_names = T,row_dend_side = "left", show_heatmap_legend = FALSE,
        colorRampPalette(c("#7BAFDE", "white", "firebrick3"))(100),border = 'black',rect_gp = gpar(col = "black", lwd = 1),
        row_names_gp = gpar(fontsize = 15),column_names_gp = gpar(fontsize = 15), cluster_rows = F, 
        cluster_columns=F, row_names_side = "left", column_names_side = "top",  column_dend_side = "bottom")

# now for mice data
temp <- gset1@assayData$exprs[ID2,]

#normalize the data by row
temp <- t(scale(t(temp)))


#plot the differential expression genes in mice 
#change the rowname from probs id to gene names
row.names(temp) <- vec_mic
Heatmap(temp, show_row_names = TRUE, show_column_names = T,row_dend_side = "left", show_heatmap_legend = FALSE,
        colorRampPalette(c("#7BAFDE", "white", "firebrick3"))(100),border = 'black',rect_gp = gpar(col = "black", lwd = 1),
        row_names_gp = gpar(fontsize = 15),column_names_gp = gpar(fontsize = 15), cluster_rows = F, 
        cluster_columns=F, row_names_side = "left", column_names_side = "top",  column_dend_side = "bottom")
