# Gradient Discent
getwd()
# problem 1.1 f(x) = (10*x1^2 + x2^2)/2
# 수렴시킬 함수 정의
fx <- function(x1, x2) {
  (10*x1^2 + x2^2)/2 # 미분가능함수!
}

dfx <- function(x1, x2) {
  c(10*x1, x2)
}


# 등고선 그리기
# 참고한 자료 : https://m.blog.naver.com/PostView.naver?isHttpsRedirect=true&blogId=pmw9440&logNo=221597742238

x1 <- x2 <- seq(-20,20,length = 100)

values <- matrix(0,length(x1),length(x2))
for ( i in 1:length(x1)) {
  for ( j in 1:length(x2)) {
    values[i,j] <- fx(x1[i],x2[j])
  }
}

contour(x1, x2, values)
points(0, 0, pch=8, cex=1)


# gradient discect algorithm
t = c(0.05, 0.1, 0.15, 0.2) # 시도해 볼 stepsize 
ep = 1e-10 # 기준이 될 epsilon

par(mfrow = c(2,2))
for ( i in 1:length(t)) {
  contour(x1, x2, values)
  points(0, 0, pch=8, cex=1, col='red')
  st = 100
  mat = px = c(15,18) # 초기값 지정
  iter = 0
  
  while ( st > ep) {
    nxtx = px - t[i]*dfx(px[1],px[2]) # next x
    # nxty = fx(nxtx[1], nxtx[2])
    st = sum(dfx(px[1],px[2])^2) # stop rule에 대조
    # py = nxty
    px = nxtx # update new x (px : present x)
    mat = rbind(mat, px) # x의 이동 경로 남기기
    iter = iter+1 # update iteration
    if (iter >= 10000) {
      print("수렴실패")
      points(mat[,1], mat[,2], type='b')
      break
      }
    if ( st <= ep ) {
      points(mat[,1], mat[,2], type='b')
      print(paste("total iteration is", iter)) # 수렴한 경우, 총 update 횟수 프린트
      print(paste("final X values", px[1], px[2]))
      }
  }9
}

# problem 1.2 backtracking line search # 근데 왜 프린트에 있는 거랑 다르지,,p.14
par(mfrow = c(1,1))
beta = alpha = 1/2 # 두 값 모두 (0,1)

contour(x1, x2, values)
points(0, 0, pch=8, cex=1, col='red')
st = 100 # stoping값 큰값으로 초기 지정(실제로는 이 값은 사용되지 않음)
mat = px = c(4, 20) # 초기값 지정
iter = 0
t = 0.6 # 초기 stepsize (크게잡아줌)

while ( st > ep) {
  nxtx = px - t*dfx(px[1], px[2]) # next x
  if ( fx(nxtx[1], nxtx[2]) > fx(px[1], px[2]) - alpha*t*sum(dfx(px[1], px[2])^2)) {
    t = beta*t 
  }
  st = sum(dfx(px[1], px[2])^2) # stop rule에 대조
  px = nxtx # update new x (px : present x)
  mat = rbind(mat, px) # x의 이동 경로 남기기
  iter = iter+1 # update iteration
  if (iter >= 10000) {
    points(mat[,1], mat[,2], type='b')
    print("수렴실패")
    break
  }
  if ( st <= ep ) {
    points(mat[,1], mat[,2], type='b')
    print(paste("total iteration is", iter)) # 수렴한 경우, 총 update 횟수 프린트
    print(paste("final X values", px[1], px[2])) 
  }
}




# subgradient discent
# problem 2. 

# 심장병 데이터 (?) : logistic 모형에서 beta추정지로 수렴하고싶음
Heart <- read.table("http://www.stat.ufl.edu/~aa/cat/data/Heart.dat", header=T)

Heart$score = c(0,2,4,5)

# 우선 찾고 싶은 값(beta 추정치를 확인하자)
fit <- glm(cbind(Heart$yes, Heart$no)~c(0,2,4,5), family=binomial(link="logit"))
summary(fit) # 우리의 타겟 값(베타값) :  0.39734

# 함수 정의
xx = c(0,2,4,5)
y = Heart$yes
n = Heart$yes + Heart$no

L = function(b0, b1) {
  sum(y*(b0+b1*xx)-n*log(1+exp(b0+b1*xx)))
}

dfL = function(b0, b1) {
  a = sum(y-n*(exp(b0+b1*xx)/(1+exp(b0+b1*xx))))
  b = sum(xx*y-n*xx*(exp(b0+b1*xx)/(1+exp(b0+b1*xx))))
  return(c(a,b))
}

# function visualization
b0 <- seq(-10,5, length = 100)
b1 <- seq(-5,3, length = 100)
yvalue <- matrix(0,length(b0),length(b0))
for ( i in 1:length(b0)) {
  for ( j in 1:length(b1)) {
    yvalue[i,j] <- L(b0[i],b1[j])
  }
}

contour(b0, b1, yvalue)
points(-3.8662481, 0.3973366, pch=8, cex=1, col='red') # fit을 통해서 구한 값 표시

# backtracking gradient discent 적용 => not working....
alpha = beta = 1/2
st = 100 # stoping값 큰값으로 초기 지정(실제로는 이 값은 사용되지 않음)
mat = px = c(0,0) # 초기값 지정
iter = 0
contour(b0, b1, yvalue)
points(-3.8662481, 0.3973366, pch=8, cex=1, col='red')
points(px[1], px[2], col='red')
t = 0.01 # 초기 stepsize (크게잡아줌)

while ( st > ep) {
  nxtx = px + t*dfL(px[1], px[2]) # next x
  if ( L(nxtx[1], nxtx[2]) > L(px[1], px[2]) - alpha*t*sum(dfL(px[1], px[2])^2)) {
    t = beta*t # update step size
  }
  st = sum(dfL(px[1], px[2])^2) # stopping rule에 대조
  px = nxtx # update new x (px : present x)
  mat = rbind(mat, px) # x의 이동 경로 남기기
  iter = iter+1 # update iteration
  if (iter >= 10000) {
    points(mat[,1], mat[,2], type='b')
    print("수렴실패")
    break
  }
  if ( st <= ep ) {
    points(mat[,1], mat[,2], type='b')
    print(paste("total iteration is", iter)) # 수렴한 경우, 총 update 횟수 프린트
    print(paste("final beta values", px[1], px[2])) 
  }
}


# gradient discent
t = c(0.0001, 0.001, 0.0015, 0.002) # 시도해 볼 stepsize 
ep = 1e-10 # 기준이 될 epsilon

par(mfrow = c(2,2))
for ( i in 1:length(t)) {
  contour(b0, b1, yvalue)
  points(-3.8662481, 0.3973366, pch=8, cex=1, col='red')
  st = 100
  mat = px = c(0,0) # 초기값 지정
  iter = 0
  
  while ( st > ep) {
    nxtx = px + t[i]*dfL(px[1],px[2]) # next x # 함수가 concave기 때문에 부호가 +가 되어야함!!!!!
    st = sum(dfL(px[1],px[2])^2) # stopping rule에 대조
    px = nxtx # update new x (px : present x)
    mat = rbind(mat, px) # x의 이동 경로 남기기
    iter = iter+1 # update iteration
    if (iter >= 10000) {
      print("수렴실패")
      points(mat[,1], mat[,2], type='b')
      break
    }
    if ( st <= ep ) {
      points(mat[,1], mat[,2], type='b')
      print(paste("total iteration is", iter)) # 수렴한 경우, 총 update 횟수 프린트
      print(paste("final beta values", px[1], px[2]))
    }
  }
}

