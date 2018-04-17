# 8.3 

# c. Gibbs Sampler

set.seed(920804)
alpha = 10
beta =5
n = 10

pt.x = function (x){
  return(rbeta(1, alpha+x, n - x + beta))
}
px.t = function (t){
  return(rbinom(1, 10, t))
}

iter    = 2000
burn.in = 100
xint    = 0
x0      = xint
x       = rep(0, iter)
theta   = rep(0, iter)

for (i in 1:iter){
  theta[i] = pt.x(x0)
  x[i]     = px.t(theta[i])
  x0       = x[i]
}

plot(x[-burn.in], theta[-burn.in])
hist(x[-burn.in])

# d. CFTP

set.seed(920804)
alpha = 10
beta  = 5
n     = 10
q.xUV = function(U,V, x){
  res = NULL
  for (t in x){
    theta = qbeta(U, alpha+t, beta +n -t)
    res   = c(res, qbinom(V, n, theta))
  }
  return(res)
}
CFTP = function(){
  U = NULL
  V = NULL
  map2zero = rep(0, n+1)
  tau      = -1
  repeat{
    U0 = runif(1,0,1)
    V0 = runif(1,0,1)
    U = c(U0, U)
    V = c(V0, V)
    XtauPlus1 = q.xUV(U = U0,V = V0, 0:n)
    if (tau==-1){
      map2zero = XtauPlus1
    }else{
      map2zero = map2zero[XtauPlus1+1]
    }
    if (length(unique(map2zero)) == 1){
      break
    }else{
      tau = tau-1
    }
  }
  return(list(tau = tau, x = map2zero, U= U, V = V))
}
x.CFTP = NULL
tau.CFTP = NULL
for (i in 1:100){
  a        = CFTP()
  x.CFTP   = c(x.CFTP, a$x)
  tau.CFTP = c(tau.CFTP, a$tau)
}
hist(tau.CFTP)
hist(x.CFTP)

# e.

set.seed(920804)
alpha = 1.001
beta  = 1
n     = 10
repeat{
  a = CFTP()
  if (a$tau <= -15){
    break
  }
}
x.axis = rep(a$tau:0,11)
y.axis = NULL
for (t in 0:10){
  y.axis = c(y.axis, rep(t, abs(a$tau)+1))
}
plot(x.axis,y.axis)
x.start = 0:n
for(i in 1:abs(a$tau)){
  x.end = q.xUV(U = a$U[i], V =a$V[i], x.start)
  for (t in 1:(n+1)){
    segments(y0 = x.start[t], y1 = x.end[t], x0=a$tau +i -1, x1 = a$tau+i)
  }
  x.start = x.end
}

# f.

a = 10
b = 5
n = 10
k=11
set.seed(k)
xx = ss
U = runif(2)
u = U[1]; v = U[2]
for(i in 1:(n+1)) xx[i] = q(ss[i], u, v)

while( length(table(xx)) != 1){
  U = rbind(U, runif(2))
  path = xx = ss
  for(tau in (dim(U)[1]:1) ){
    u = U[tau, 1]; v = U[tau, 2]
    for(i in 1:length(ss) ) xx[i] = q(xx[i], u, v)
    path = rbind(xx, path)
  }
  #path
}
#path

theta = gibbs = perfect = rep(0, 20)
gibbs[1] = perfect[1] = 0
theta[1] = rbeta(1, a+gibbs[1], b+n-gibbs[1])

for(i in 2:20){
  ## perfect
  u = runif(1)
  v = runif(1)
  perfect[i] = q(perfect[i-1], u, v)
  
  ## gibbs
  gibbs[i] = rbinom(1, n, prob = theta[i-1])
  theta[i] = rbeta(1, a+gibbs[i-1], b+n-gibbs[i-1])
}


plot(rep((0:19), 11), rep(rep(0:10), each = 20), cex = 0.3, xlab = "tau", ylab = "sample space" )
points(0:19, perfect, "l")
points(0:19, gibbs, lty = 2, "l")
legend("topleft", legend = c("CFTP", "Gibbs"), lty=1:2)

