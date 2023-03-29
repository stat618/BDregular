metadata <- read.table("./data/lajosJCO2006_genomicMatrix", header=TRUE, sep="\t")

new_df <- metadata[, -1, drop=FALSE]
