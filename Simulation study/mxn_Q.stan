functions {
  matrix kronecker_prod(matrix A, matrix B) {
    int m = rows(A);
    int n = cols(A);
    int p = rows(B);
    int q = cols(B);
    matrix[m*p, n*q] C;
    for (i in 1:m) {
      for (j in 1:n) {
        int row_start = (i - 1) * p + 1;
        int row_end = i*p;
        int col_start = (j - 1) * q + 1;
        int col_end = j*q;
        C[row_start:row_end, col_start:col_end] = A[i, j] * B;
      }
    }
    return C;
  }
  vector kp_diag(vector A, vector B) {
    int m = rows(A);
    int p = rows(B);
    vector[m*p] C;
    for (i in 1:m) {
      int start = (i-1)*p + 1;
      int end = i*p;
      C[start:end] = A[i]*B;
    }
    return C;
  }
}

data {
  int<lower=0> n_node;
  int<lower=0> n_time;
  
  vector[n_node*n_time] X; //Target to evaluate the density on
  vector[n_node*n_time] Y; //For mod=2, this is X1
  
  array[2] int<lower=0, upper=n_node> zero;
  
  int<lower=0, upper=3> mod;
  
  matrix[n_time, n_time] evecQAR;
  vector[n_time] seq;
  
  vector[2] QX; //First column of Q, QX[1] = q1 and QX[2] = q2
}

parameters {
  real<lower=-1, upper=1> rho;
  corr_matrix[n_node] cor;
}

transformed parameters {
  vector[n_time] evalQAR = 1 + rho^2 + 2*rho*seq;
  matrix[n_node, n_node] evecNode;
  vector[n_node] evalNode;
  
  (evecNode, evalNode) = eigendecompose_sym(cor);
  for (i in 1:n_node) {
    if (evecNode[1,i] < 0) {
      evecNode[:,i] = -1*evecNode[:,i];
    }
  }
}

model {
  if (zero[1] != 0) {
    cor[zero[1],zero[2]] ~ normal(0, 0.001);
  } 
  {
    matrix[n_node*n_time, n_node*n_time] QR = kronecker_prod(evecQAR, evecNode);
    vector[n_node*n_time] AB = kp_diag(inv(evalQAR), evalNode);
    if (mod == 0) {
      QR'*X ~ normal(0, sqrt(AB));
    } else if (mod == 1) {
      QR'*X ~ normal(0, sqrt(QX[1]*QX[1]*(AB-1)+1));
    } else if (mod == 2) {
      QR'*X ~ normal(0, sqrt(QX[2]*QX[2]*(AB-1)+1));
    } else {
      QR'*X ~ normal((QX[1]*QX[2]*(AB-1) ./ (QX[1]*QX[1]*(AB-1)+1)) .* (QR'*Y), sqrt(AB ./ (QX[1]*QX[1]*(AB-1)+1)));
    }
  }
}

generated quantities {
  real ll; 
  {
    matrix[n_node*n_time, n_node*n_time] QR = kronecker_prod(evecQAR, evecNode);
    vector[n_node*n_time] AB = kp_diag(inv(evalQAR), evalNode);
    if (mod == 0) {
        ll = normal_lpdf(QR'*X | 0, sqrt(AB));
      } else if (mod == 1) {
        ll = normal_lpdf(QR'*X | 0, sqrt(QX[1]*QX[1]*(AB-1)+1));
      } else if (mod == 2) {
        ll = normal_lpdf(QR'*X | 0, sqrt(QX[2]*QX[2]*(AB-1)+1));
      } else {
        ll = normal_lpdf(QR'*X | (QX[1]*QX[2]*(AB-1) ./ (QX[1]*QX[1]*(AB-1)+1)) .* (QR'*Y), sqrt(AB ./ (QX[1]*QX[1]*(AB-1)+1)));
      }
  }
}
