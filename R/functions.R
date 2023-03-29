# 这段代码用于实现矩阵的中心化处理，
# 即将每一列减去该列的均值。
# 通过这个函数，可以避免数据集中的变量之间
# 存在较大的差异，使得数据更容易进行比较和分析。
dmean <- function(A){
  # 计算矩阵 A 中每一列的均值，结果保存在 a 变量中
  a = colMeans(A)
  # 对矩阵 A 进行中心化处理，即矩阵 A 的转置减去均值
  # 再次转置
  A = t(t(A) - a)
  return(A)
}


# 标准化
standard=function(A){
  # 获取矩阵A每一列的均值
  a=colMeans(A)
  # 使用apply函数计算矩阵A每一列的标准差
  # 并赋值给变量b，其中2表示按列计算标准差
  b=apply(A,2,sd)
  # 将矩阵A转置，减去均值a
  # 再将结果再次转置，赋值给变量B
  B=t(t(A)-a)
  # 将变量B转置后除以每一列的标准差b
  # 再将结果转置回来，得到标准化后的矩阵
  B=t(t(B)/b)
  return(B)
}



# 将两个矩阵合并为一个块对角矩阵
blockdiag=function(A,B){
  # 如果矩阵A有多列，则执行以下操作
  if(length(A)>1){
    # 定义p1和p2分别为矩阵A和矩阵B的列数
    p1=dim(A)[2];p2=dim(B)[2]
    # 创建一个大小为(p1+p2) x (p1+p2)的零矩阵C
    # 然后将其转化为对角矩阵
    C=diag(p1+p2)
    # 将矩阵A放在块对角矩阵的左上角
    C[1:p1,1:p1]=A
    # 将矩阵B放在块对角矩阵的右下角
    C[c((p1+1):(p1+p2)),c((p1+1):(p1+p2))]=B
  # 否则（即，如果矩阵A只有一列），执行以下操作
  }else{
    p1=1;
    # 定义p2为矩阵B的列数
    p2=dim(B)[2]
    # 创建一个大小为
    # (p1+p2) x (p1+p2)的零矩阵C
    # 然后将其转化为对角矩阵
    C=diag(p1+p2)
    # 将矩阵A的唯一值放在块对角矩阵的左上角
    C[1,1]=A
    # 将矩阵B放在块对角矩阵的右下角
    C[c(2:(1+p2)),c(2:(1+p2))]=B
  }
  return(C)
}


# SCAD的阈值函数
# 接受三个参数：a，lam和gamma
# 其中a是一个数字或向量，lam是一个正数
# gamma是一个可选参数，默认值为3.7
dSCAD=function(a,lam,gamma=3.7){
  # 将输入向量a的每个元素转换成它的绝对值
  # 并将其赋值给变量a，以确保所有值都是非负数
  a=abs(a)
  z=a
  # 对于z中小于阈值lam的元素，将它们设置为lam
  z[a<lam]=lam
  # 对于z中大于阈值lam的元素，将它们替换为如下
  # 其中i是当前元素的索引
  # 这个公式可以想象成是将原来的值线性地削减到比lam低的值
  # 这是SCAD方法的关键
  z[a>lam]=(gamma*lam-z[a>lam])/(gamma-1)
  # 将z中大于gamma * lam的元素设置为0
  z[a>(gamma*lam)]=0
  return(z)
}



# 这段代码实现的是一种变种的软阈值函数
# 它将每个值减去一个量`a/gamma`
# 如果值小于`gamma*lam`则将其保留，否则设为0
# 该函数可以在信号处理、压缩感知、机器学习等领域中使用
dMCP=function(a,lam,gamma=3){
  a=abs(a)
  z=lam-a/gamma
  z[a>(gamma*lam)]=0
  return(z)
}



# 这段代码实现的是一种标准的软阈值函数
# 它将每个值减去一个量b
# 如果值小于0则将其设为0
soft=function(a,b){
  c=abs(a)-b
  c[c<0]=0
  c=c*sign(a)
  return(c)
}




shrinking=function(A,eps){
  # 通过eigen()函数计算矩阵A的特征值，将其保存在向量d中
  d=eigen(A)$values
  d1=min(d)
  dp=max(d)
  us=max(eps,(d1+dp)/2)
  # 计算uf变量，该变量是与对角线上的元素有关的一个值
  # 具体地，它等于(d-d1)^2的和除以d-d1的和
  uf=sum((d-d1)^2)/sum(d-d1)
  u=max(us,uf)
  # 计算alpha变量，该变量是一个缩放因子
  # 用于调整矩阵A的大小
  # 它等于(eps-d1)/(u-d1)的差异
  # 其中eps是输入参数
  # 该计算基于一个假设
  # 即矩阵A的特征值分布在[d1, dp]范围内
  alpha=1-(eps-d1)/(u-d1)
  B=alpha*as.matrix(A)+(1-alpha)*u*diag(ncol(A))
  B=cov2cor(B)
  S=list(est=B,alpha=alpha,u=u)
  return(S)
}




SCADthreshold=function(S,lam,k=3){
  a=as.vector(S)
  b=abs(a)
  for(i in 1:k){
    c=dSCAD(b,lam)
    b=soft(b,c)
  }
  d=sign(a)*b
  d=matrix(d,ncol(S),ncol(S))
  diag(d)=1
  d=shrinking(d,eps=0.001)
  d=cov2cor(d$est)
  return(d)
}




MCPthreshold=function(S,lam,k=3){
  a=as.vector(S)
  b=abs(a)
  for(i in 1:k){
    c=dMCP(b,lam)
    b=soft(b,c)
  }
  d=sign(a)*b
  d=matrix(d,ncol(S),ncol(S))

  diag(d)=1
  d=shrinking(d,eps=0.001)
  d=cov2cor(d$est)
  return(d)
}


