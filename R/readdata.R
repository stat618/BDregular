metadata <- read.table("./data/testdata", header=TRUE, sep="\t")

new_df <- metadata[, -1, drop=FALSE]
