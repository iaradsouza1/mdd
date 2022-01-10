
# Read and summarise transcript to gene level -----------------------------

library(tximport)
library(GenomicFeatures)

# Load annotation 
load("results/important_variables/ann.rda")

# Load counts from kallisto -----------------------------------------------
gtf <- "data/genome/Homo_sapiens.GRCh38.97.gtf.gz"
txdb.filename <- "data/genome/Homo_sapiens.GRCh38.97.gtf.sqlite"

if(!("Homo_sapiens.GRCh38.97.gtf.sqlite" %in% list.files("data/genome"))) {
  txdb <- makeTxDbFromGFF(gtf, format = "gtf")
  saveDb(txdb, txdb.filename)
}

# Load db
txdb <- loadDb(txdb.filename)
txdf <- AnnotationDbi::select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]
tx2gene <- data.frame(tx = txdf$TXNAME, gene = txdf$GENEID, stringsAsFactors = F)

# Read files
files <- list.files(path = "data/kallisto", pattern="tsv", recursive = TRUE, full.names = TRUE)
files <- files[sapply(rownames(ann), function(x) grep(x, files))]
names(files) <- rownames(ann)
txi <- tximport(files = files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = T)

# Save
if(!dir.exists("results/txi")) {
  dir.create("results/txi")
}
save(txi, file = "results/txi/txi_gene.rda")
