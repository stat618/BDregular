# metadata <- read.table("./data/testdata", header=TRUE, sep="\t")
#
# new_df <- metadata[, -1, drop=FALSE]
#
# new_df_transposed <- t(new_df)
#
# data6 =t((read.table("./data/testdata", header=TRUE, sep="\t"))[, -1, drop=FALSE])
#
# result = data6[1:30, 1:900]
#
# AAA = BDregular(result)
# AAA = BDregular::BDregular(result)
#
#
# library(glasso)
# BBB <- glasso(result, rho = 0.01)
# # 这个运算的特别慢，差的极其多，甚至能差大概十倍
#
# # 定义一个函数来计算两个向量之间的余弦相似度
# cosine_similarity <- function(x, y) {
#   sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
# }
#
# # 计算两个矩阵之间的余弦相似度
# similarity <- cosine_similarity(AAA, BBB)
