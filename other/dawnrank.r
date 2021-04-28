load("~/Documents/nempi_backup/backup3/TCGA-BRCA_nempi.rda")
M <- Pmut
patients.mut <- unlist(lapply(colnames(M), function(x) {
    y <- paste(unlist(strsplit(x, "\\-"))[1:3], collapse = "-")
    return(y)
}))
snv_matrix <- M

D <- readRDS("~/Documents/mutclust/TCGA-BRCA_counts.rds")
DT <- D$T
DN <- D$N
patients.exprs <- unlist(lapply(colnames(DT), function(x) {
    y <- paste(unlist(strsplit(x, "\\-"))[1:3], collapse = "-")
    return(y)
}))
# colnames(DN) <- unlist(lapply(colnames(DN), function(x) {
#     y <- paste(unlist(strsplit(x, "\\-"))[1:3], collapse = "-")
#     return(y)
# }))

DTN <- cbind(DT, DN)
samples <- colnames(DT)[which(patients.exprs %in% patients.mut)]
expression_matrix <- DTN
sample_origins <- c(rep("tumor",ncol(DT)),rep("normal",ncol(DN)))

colnames(expression_matrix) <- unlist(lapply(colnames(expression_matrix), function(x) {
    y <- paste(unlist(strsplit(x, "\\-"))[1:4], collapse = "-")
    y <- unlist(strsplit(y, ""))
    y <- paste(y[-length(y)], collapse = "")
    return(y)
}))

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genesymbols <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), values=rownames(expression_matrix), mart= mart)
rownames(expression_matrix) <- genesymbols[match(rownames(expression_matrix), genesymbols[, 1]), 2]

library(PRODIGY)
# Load STRING network data 
data(STRING_network)
network = STRING_network

# Get differentially expressed genes (DEGs) for all samples
expression_matrix = expression_matrix[which(rownames(expression_matrix) %in% unique(c(network[,1],network[,2]))),]

library(DESeq2)
DEGs = get_DEGs(expression_matrix,samples,sample_origins=sample_origins,beta=2,gamma=0.05)

names(DEGs) <- unlist(lapply(names(DEGs), function(x) {
    y <- paste(unlist(strsplit(x, "\\-"))[1:4], collapse = "-")
    y <- unlist(strsplit(y, ""))
    y <- paste(y[-length(y)], collapse = "")
    return(y)
}))

# Run PRODIGY
samples <- colnames(snv_matrix)
all_patients_scores = PRODIGY_cohort(snv_matrix,expression_matrix,network=network,samples=samples,DEGs=DEGs,alpha=0.05,pathwayDB="reactome",num_of_cores=1,sample_origins=sample_origins,write_results = TRUE, results_folder = "prodigy/",beta=2,gamma=0.05,delta=0.05)

tmp <- PRODIGY(snv_matrix,expression_matrix,network,
               sample=names(DEGs)[1],DEGs[[1]],
               sample_origins=sample_origins,
               write_results = TRUE,
               results_folder = "prodigy/")

# Get driver gene rankings for all samples 
results = analyze_PRODIGY_results(all_patients_scores)

## try dawnrank:

library(DawnRank)

DT2 <- DT
colnames(DT2) <- unlist(lapply(colnames(DT2), function(x) {
    y <- paste(unlist(strsplit(x, "\\-"))[1:4], collapse = "-")
    y <- unlist(strsplit(y, ""))
    y <- paste(y[-length(y)], collapse = "")
    return(y)
}))

DT2 <- DT2[,which(colnames(DT2) %in% colnames(M))]
M <- M[,which(colnames(M) %in% colnames(DT2))]
DT2 <- DT2[,colnames(M)]

DT2.log2 <- log2(DT2 + 1)
DN.log2 <- log2(DN + 1)
file <- "~/Documents/testing/nempi/other/nempi_dawnrank_data.rds"
if (file.exists(file)) {
    normalizedDawn <- readRDS(file)
} else {
    normalizedDawn <- DawnNormalize(tumorMat = DT2.log2, normalMat = DN.log2)
    saveRDS(normalizedDawn, file = file)
}

n <- length(unique(as.vector(network[,1:2])))
pathway <- matrix(0, n, n)
rownames(pathway) <- colnames(pathway) <- sort(unique(as.vector(network[,1:2])))
a <- match(network[,1],rownames(pathway))
b <- match(network[,2],rownames(pathway))
pathway[cbind(a,b)] <- 1

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genesymbols <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), values=rownames(normalizedDawn), mart= mart)
rownames(normalizedDawn) <- genesymbols[match(rownames(normalizedDawn), genesymbols[, 1]), 2]

tmp <- which(rownames(M) %in% rownames(normalizedDawn))
M <- M[tmp,]
tmp <- which(rownames(normalizedDawn) %in% rownames(M))
normalizedDawn <- normalizedDawn[tmp,]
normalizedDawn <- rowsum(normalizedDawn, rownames(normalizedDawn))/as.numeric(table(rownames(normalizedDawn)))

tmp <- rownames(M)[which(!(rownames(M) %in% rownames(pathway)))]
pathway2 <- matrix(0, length(tmp)+nrow(pathway),length(tmp)+ncol(pathway))
pathway2[1:nrow(pathway),1:ncol(pathway)] <- pathway
rownames(pathway2) <- colnames(pathway2) <- c(rownames(pathway),tmp)

tmp <- which(colnames(M) %in% colnames(normalizedDawn))
M <- M[,tmp]
M <- M[rownames(M),colnames(M)]
normalizedDawn <- normalizedDawn[rownames(M),colnames(M)]
pathway2 <- pathway2[rownames(M),rownames(M)]

my.genes <- c("TBX3","MAP2K4","GATA3","CBFB","NCOR1","PTPRD","GPS2","CDKN1B")

## load("~/Documents/testing/nempi/other/Dawnrank_temp.rda")

file <- "~/Documents/testing/nempi/other/nempi_dawnrank.rds"
if (!file.exists(file)) {
    # get the DawnRank Score Get some coffee, this might take a while!
    dawnRankScore <- DawnRank(adjMatrix = pathway2, mutationMatrix = M,
                          expressionMatrix = normalizedDawn, mu = 3,
                          parallel = NULL)
    saveRDS(dawnRankScore,file=file)
} else {
    dawnRankScore <- readRDS(file)
}

rank.mat <- matrix(NA, length(my.genes), ncol(M))
rownames(rank.mat) <- my.genes
colnames(rank.mat) <- colnames(M)

for (gene in 1:length(my.genes)) {
    idx <- which(rownames(dawnRankScore[[2]]) == my.genes[gene])
    if (length(idx) > 0) {
        rank.mat[gene,] <- dawnRankScore[[2]][idx,]
    }
}

Msub <- M[which(rownames(M) %in% my.genes),]
Msub <- Msub[rownames(rank.mat),]
rank.mat.sub <- rank.mat[,which(apply(Msub, 2, sum) != 0)]
Msub <- Msub[,which(apply(Msub, 2, sum) != 0)]
cut <- 90
TP <- sum(rank.mat.sub > cut & Msub == 1)
FP <- sum(rank.mat.sub > cut & Msub == 0)
TN <- sum(rank.mat.sub <= cut & Msub == 0)
FN <- sum(rank.mat.sub <= cut & Msub == 1)

TP/(TP+FP)

TP/(TP+FN)

gamsave <- rank.mat.sub/100
Gamma <- Msub
auc <- roc <- 0
ppv <- rec <- spec <- NULL
for (cut in c(2,seq(1,0, length.out = 100),-1)) {
    gamtmp <- apply(gamsave, 2, function(x) {
        y <- x*0
        y[which(x > cut)] <- 1
        return(y)
    })
    gamtmp2 <- Gamma
    tp <- sum(gamtmp == 1 & gamtmp2 == 1)
    fp <- sum(gamtmp == 1 & gamtmp2 == 0)
    tn <- sum(gamtmp == 0 & gamtmp2 == 0)
    fn <- sum(gamtmp == 0 & gamtmp2 == 1)
    ppvtmp <-  tp/(tp+fp)
    if (is.na(ppvtmp)) { ppvtmp <- 0.5 }
    rectmp <- tp/(tp+fn)
    spectmp <- 1-tn/(tn+fp)
    if (length(ppv) > 0) {
        auc <- auc +
            (rectmp-rec[length(rec)])*(ppvtmp+ppv[length(ppv)])/2
        roc <- roc +
            (spectmp-spec[length(spec)])*(rectmp+rec[length(rec)])/2
    }
    ppv <- c(ppv, ppvtmp)
    rec <- c(rec, rectmp)
    spec <- c(spec, spectmp)
}

par(mfrow=c(1,3))
plot(ppv)
plot(rec)
plot(spec)


pdf("~/Documents/temp.pdf", width = 16, height = 8)
tmp <- rank.mat.sub
colnames(tmp) <- NULL
csc <- numeric(ncol(tmp))
csc[which(apply(Msub, 2, sum) != 0)] <- 4
epiNEM::HeatmapOP(tmp/100,
                  breaks=seq(0,1,length.out=100),
                  col="RdBu",
                  bordercol="transparent",
                  colSideColors=csc)
dev.off()

## check with mutation matrix:

if (file.exists("nempi_dawnrank_agg.rds")) {
    aggregateDawnRankScore <- condorcetRanking(scoreMatrix =dawnRankScore[[2]],mutationMatrix = M,parallel = NULL)
    saveRDS(aggregateDawnRankScore,file="nempi_dawnrank_agg.rds")
} else {
    aggregateDawnRankScore <- readRDS("nempi_dawnrank_agg.rds")
}

tmp <- aggregateDawnRankScore[[2]]
tmp[which(names(tmp) %in% my.genes)]

## process mutations for 2020 plus (makes no sense?):

mutSave <- mut
for (i in 1:4) {
    mut[[i]] <- as.data.frame(mut[[i]])
}

m2020 <- do.call("rbind",mut)
write.table(m2020, file = "~/Documents/2020plus/data/brca.txt")

blad <- read.delim("~/Documents/2020plus/data/bladder.txt")

brca <- read.table("~/Documents/2020plus/data/brca.txt")

colnames(brca)[which(colnames(brca) == "Gene")] <- "EnsembleId"
colnames(brca)[which(colnames(brca) == "SYMBOL")] <- "Gene"
colnames(brca)[which(colnames(brca) == "Tumor_Sample_Barcode")] <- "Tumor_Sample"
colnames(brca)[which(colnames(brca) == "Tumor_Seq_Allele1")] <- "Tumor_Allele"
colnames(brca)[which(colnames(brca) == "HGVSp")] <- "Protein_Change"
colnames(brca)[which(colnames(brca) == "HGVSc")] <- "DNA_Change"
brca <- data.frame(brca, Tumor_Type = rep("Breast Cancer", nrow(brca)))

# "Symbol" "Tumor_Sample_Barcode" "?" "Chromosome" "Start_Position" "End_Position" "Variant_Classification" "Reference_Allele" "Allele/Tumor_Seq_Allele1/2!!!" "HGVSp/HGVSp_Short" "HGVSc"

brca <- brca[,which(colnames(brca) %in% colnames(blad))]
brca <- brca[,colnames(blad)]
brca[,"Tumor_Sample"] <- unlist(lapply(brca[,"Tumor_Sample"], function(x) {
    y <- paste(unlist(strsplit(x,"-"))[1:3],collapse="-")
    return(y)
}))
brca$Start_Position <- as.numeric(brca$Start_Position)
brca$End_Position <- as.numeric(brca$End_Position)
## brca_save <- brca
tmp <- which(table(brca$Start_Position) >= 3)
brca <- brca[which(brca$Start_Position %in% names(tmp)),]
brca <- brca[!duplicated(brca$Start_Position),]

write.table(brca, file = "~/Documents/2020plus/data/brca_reduced.txt", quote = FALSE, sep = "\t")

## 2020+

rm error.txt
rm output.txt

ram=32000 # data: 64gb, wilcox sc: 8gb for 32 nodes, edger: 64gb, voom: 64gb
queue=24
cores=1

bsub -M ${ram} -q normal.${queue}h -n ${cores} -e error.txt -o output.txt -R "rusage[mem=${ram}]" 'snakemake -s Snakefile pretrained_predict -p --cores 1 \
--config mutations="data/brca.txt" output_dir="brca/" trained_classifier="data/2020plus_100k.Rdata"'

## read results:

res2020 <- read.delim("~/Mount/Eulershare/2020plus/2020plus/brca/pretrained_output/results/r_random_forest_prediction.txt")









