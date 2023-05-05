#include <RcppEigen.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

using namespace Rcpp;

// 写一段Rcpp的代码，定义一个名为dmean的函数，它接受一个矩阵A作为输入。函数首先计算矩阵A每列的平均值并将结果存储在向量a中。然后，使用矩阵减法将矩阵A的每一列减去对应的平均值。最后，返回处理后的矩阵A。
// [[Rcpp::export]]
NumericMatrix dmean(NumericMatrix A){
  int row = A.nrow();
  int col = A.ncol();
  NumericMatrix B(row, col);
  NumericVector a(col);
  for(int i = 0; i < col; i++){
    double sum = 0;
    for(int j = 0; j < row; j++){
      sum += A(j, i);
    }
    double average = sum / row;
    a[i] = average;
  }
  for(int i = 0; i < col; i++){
    for(int j = 0; j < row; j++){
      B(j, i) = A(j, i) - a[i];
    }
  }
  return B;
}


//写一段Rcpp的代码，定义一个名为standard的函数，它接受一个矩阵A作为输入。函数首先计算矩阵A每列的平均值并将结果存储在向量a中。然后，使用apply函数计算矩阵A每列的标准差并将结果存储在向量b中。接下来，使用矩阵减法将矩阵A的每一列减去对应的平均值，然后再除以对应的标准差。最后，返回处理后的矩阵B。
// [[Rcpp::export]]
NumericMatrix standard(NumericMatrix A){
  int row = A.nrow();
  int col = A.ncol();
  NumericMatrix B(row, col);
  NumericVector a(col);
  NumericVector b(col);
  for(int i = 0; i < col; i++){
    double sum = 0;
    for(int j = 0; j < row; j++){
      sum += A(j, i);
    }
    double average = sum / row;
    a[i] = average;
  }
  for(int i = 0; i < col; i++){
    double sum = 0;
    for(int j = 0; j < row; j++){
      sum += pow(A(j, i) - a[i], 2);
    }
    double standard = sqrt(sum / (row - 1));
    b[i] = standard;
  }
  for(int i = 0; i < col; i++){
    for(int j = 0; j < row; j++){
      B(j, i) = (A(j, i) - a[i]) / b[i];
    }
  }
  return B;
}

//写一段Rcpp的代码，定义一个名为blockdiag的函数，它接受两个矩阵A和B作为输入。函数首先检查矩阵A的长度是否大于1，如果是，则执行以下操作：获取矩阵A和B的列数分别为p1和p2，然后创建一个大小为(p1+p2)*(p1+p2)的对角矩阵C。接着，将矩阵A赋值给矩阵C的左上角，将矩阵B赋值给矩阵C的右下角。如果矩阵A的长度不大于1，则执行以下操作：将变量p1设为1，获取矩阵B的列数为p2，然后创建一个大小为(p1+p2)*(p1+p2)的对角矩阵C。接着，将矩阵A(此时应为一个标量)赋值给矩阵C的左上角，将矩阵B赋值给矩阵C的右下角。最后，返回处理后的矩阵C.
// [[Rcpp::export]]
NumericMatrix blockdiag(NumericMatrix A, NumericMatrix B){
  int p1 = A.ncol();
  int p2 = B.ncol();
  NumericMatrix C(p1 + p2, p1 + p2);
  if(p1 > 1){
    for(int i = 0; i < p1; i++){
      for(int j = 0; j < p1; j++){
        C(i, j) = A(i, j);
      }
    }
    for(int i = 0; i < p2; i++){
      for(int j = 0; j < p2; j++){
        C(i + p1, j + p1) = B(i, j);
      }
    }
  }else{
    for(int i = 0; i < p2; i++){
      for(int j = 0; j < p2; j++){
        C(i + p1, j + p1) = B(i, j);
      }
    }
  }
  return C;
}


//写一段Rcpp的代码，定义了一个名为 dSCAD 的函数。该函数接受三个参数：a（a是一个矩阵），lam 和 gamma。其中 gamma 的默认值为 3.7。函数首先对输入的 a 取绝对值，然后创建一个新变量 z 并将其赋值为 a。接下来，函数根据 a 的值对 z 进行修改。如果 a < lam，则将 z 设为 lam；如果 a > lam，则将 z 设为 (gamma*lam-z[a>lam])/(gamma-1)；如果 a > gamma*lam，则将 z 设为 0。最后，函数返回 z。
// [[Rcpp::export]]
NumericMatrix dSCAD(NumericMatrix a, double lam, double gamma = 3.7){
  int row = a.nrow();
  int col = a.ncol();
  NumericMatrix z(row, col);
  for(int i = 0; i < row; i++){
    for(int j = 0; j < col; j++){
      z(i, j) = a(i, j);
    }
  }
  for(int i = 0; i < row; i++){
    for(int j = 0; j < col; j++){
      if(a(i, j) < lam){
        z(i, j) = lam;
      }else if(a(i, j) > lam && a(i, j) < gamma * lam){
        z(i, j) = (gamma * lam - a(i, j)) / (gamma - 1);
      }else if(a(i, j) > gamma * lam){
        z(i, j) = 0;
      }
    }
  }
  return z;
}


//写一段Rcpp的代码，定义了一个名为 dMCP 的函数。该函数接受三个参数：a（a是一个矩阵），lam 和 gamma。其中 gamma 的默认值为 3。函数首先对输入的 a 取绝对值，然后创建一个新变量 z 并将其赋值为 lam-a/gamma。接下来，如果 a > gamma*lam，则将 z 设为 0。最后，函数返回 z。这个函数是一个阈值函数，用于在统计建模中进行变量选择和正则化。
// [[Rcpp::export]]
NumericMatrix dMCP(NumericMatrix a, double lam, double gamma = 3){
  int row = a.nrow();
  int col = a.ncol();
  NumericMatrix z(row, col);
  for(int i = 0; i < row; i++){
    for(int j = 0; j < col; j++){
      z(i, j) = lam - a(i, j) / gamma;
    }
  }
  for(int i = 0; i < row; i++){
    for(int j = 0; j < col; j++){
      if(a(i, j) > gamma * lam){
        z(i, j) = 0;
      }
    }
  }
  return z;
}

//写一段Rcpp的代码，定义了一个名为 sign 的函数，函数的输入为一个数字x，如果x>0，则返回1；如果x<0，则返回-1；如果x=0，则返回0。
// [[Rcpp::export]]
int sign(double x){
  if(x > 0){
    return 1;
  }else if(x < 0){
    return -1;
  }else{
    return 0;
  }
}

//写一段Rcpp的代码，定义了一个名为 soft 的函数。
// 该函数接受两个参数：a 和 b(a和b都是矩阵)。
// 函数首先创建一个新变量 c 并将其赋值为 abs(a)-b。
// 接下来，如果 c < 0，则将 c 设为 0。
// 然后，将 c 乘以 a 的符号。最后，函数返回 c。
// 这个函数实现了软阈值算法，用于在统计建模中进行变量选择和正则化。
// [[Rcpp::export]]
NumericMatrix soft(NumericMatrix a, NumericMatrix b){
  int row = a.nrow();
  int col = a.ncol();
  NumericMatrix c(row, col);
  for(int i = 0; i < row; i++){
    for(int j = 0; j < col; j++){
      c(i, j) = abs(a(i, j)) - b(i, j);
    }
  }
  for(int i = 0; i < row; i++){
    for(int j = 0; j < col; j++){
      if(c(i, j) < 0){
        c(i, j) = 0;
      }
    }
  }
  for(int i = 0; i < row; i++){
    for(int j = 0; j < col; j++){
      c(i, j) = c(i, j) * sign(a(i, j));
    }
  }
  return c;
}


