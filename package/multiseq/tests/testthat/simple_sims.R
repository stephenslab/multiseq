set.seed(100)
lambda = c(1,10,100)
x = list()
for(i in 1:length(lambda)){
  x[[i]] = rpois(1024,rep(lambda[i],1024))
  x.m[[i]] = multiseq(x[[i]])
}

plot(exp(x.m[[1]]$baseline.mean))
plot(exp(x.m[[2]]$baseline.mean))
plot(exp(x.m[[3]]$baseline.mean))

y.m = multiseq(c(x[[1]],x[[2]]))
