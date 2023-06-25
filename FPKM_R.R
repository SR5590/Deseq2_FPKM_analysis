# Read the feature count result file
featureCounts <- read.table("FeatureCount_GRCm39.tsv", header = TRUE, stringsAsFactors = FALSE)
colnames(featureCounts) <- gsub("Aligned.sortedByCoord.out.bam", "", colnames(featureCounts))

# Extract necessary columns for FPKM
geneNames <- featureCounts$GeneID
geneLengths <- featureCounts$Length
expressionData <- featureCounts[, -(1:2)]

# Calculate total mapped reads per sample for FPKM
totalMappedReads <- colSums(expressionData)

# Replicate gene lengths to match the number of samples
geneLengths_replicated <- rep(geneLengths, each = ncol(expressionData))

# Calculate FPKM values
fpkmData <- data.frame(Gene = geneNames)
expressionValues <- expressionData
fpkm <- 1e9 * expressionValues / (geneLengths_replicated * totalMappedReads)

# Assign FPKM values to the corresponding sample names
colNames <- colnames(expressionValues)
fpkmData[colNames] <- fpkm

# Print FPKM values
write.table(fpkmData, file = "Fpkm_results_allsamples.tsv", sep = "\t", quote = FALSE, col.names = NA)
