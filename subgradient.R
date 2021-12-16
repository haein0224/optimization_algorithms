# Subgradient# Dataset ( 적당한 parameter값 이용해서 lasso solution 구하기)
# 참고자료 : https://hydragon-cv.info/entry/Lasso-Regression-3
library(glmnet)
data=read.table("C:/R_directory/conv_opt/conv/data/prostate.txt")

# train 데이터만 사용할 예정
tr.X = data[data[, "train"], ][,1:8]
tr.y = data[data[, "train"], ][,9]

tr.X = as.matrix(scale(tr.X,T,T)) # 정규화 & 매트릭스화해서 사용
tr.y = scale(tr.y, T, F)
n = nrow(tr.X)

fb = function(pp, ll) { # minimize 하고자하는 식
  if (ll <0) print("lambda value는 0 이상이어야함")
  1/2*sum((tr.y-tr.X%*%pp)^2)+ll*sum(abs(pp))
}


subfb = function(pp, ll) {
  front = -(t(tr.X)%*%(tr.y-tr.X%*%pp))
  subg = NULL
  for (i in 1:length(pp)) {
    if ( pp[i] > 0 ) {
      subg[i] = front[i]  + ll
      } else if ( pp[i] < 0 ) {
      subg[i] = front[i]  - ll
      } else if ( pp[i] == 0) {
      subg[i] = front[i]  + (ll*runif(1,-1,1))
      }
  }
  return(subg)
}

# 하나만 할때
# library(lm.beta)
# comp = lm.beta(lm(tr.y~tr.X))$standardized.coefficients[-1]
par(mfrow=c(2,4))

bestbeta = pbeta1 = c(0.5, -0.5,-0.5,-0.5,-0.5,0.5,-0.5,-0.5)
mat = t(matrix(pbeta1))
t = 0.05 # step size
iter = 0
lamb = 0.3
fit <- glmnet(x = tr.X, tr.y, alpha=1, lambda=lamb)

while ( iter < 10000 ) {
  iter = iter+1 # update iteration
  nbeta = pbeta1 - t/iter*subfb(pbeta1, n*lamb) # next beta  ## 문제 되던곳 찾음!!!
  if ( fb(nbeta,lamb) < fb(pbeta1,lamb) ) { 
    bestbeta = nbeta
  }
  pbeta1 = nbeta # update new beta (pbeta : present beta)
  mat = rbind(mat, pbeta1) # beta의 변화 저장
  if ( iter == 10000 ) {
    for (i in 1:ncol(mat)) {
      plot(mat[,i])
    }
    print(cbind(bestbeta, fit$beta))
    break
  }
}


