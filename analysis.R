source("fun.R")

#---set seed for reproducibility---

set.seed(1)

#---import data and define binary exposure---

library(ivtools)
data(VitD)
n <- nrow(VitD)
VitD$vitdcat <- as.numeric(VitD$vitd>30)

#---compute observed probabilities---

#pzxy
p000 <- mean(VitD$filaggrin==0 & VitD$vitdcat==0 & VitD$death==0)
p001 <- mean(VitD$filaggrin==0 & VitD$vitdcat==0 & VitD$death==1)
p010 <- mean(VitD$filaggrin==0 & VitD$vitdcat==1 & VitD$death==0)
p011 <- mean(VitD$filaggrin==0 & VitD$vitdcat==1 & VitD$death==1)
p100 <- mean(VitD$filaggrin==1 & VitD$vitdcat==0 & VitD$death==0)
p101 <- mean(VitD$filaggrin==1 & VitD$vitdcat==0 & VitD$death==1)
p110 <- mean(VitD$filaggrin==1 & VitD$vitdcat==1 & VitD$death==0)
p111 <- mean(VitD$filaggrin==1 & VitD$vitdcat==1 & VitD$death==1)

#---estimate bounds for a and b---

p00 <- p000+p100
p01 <- p001+p101
p10 <- p010+p110
p11 <- p011+p111
pa <- c(p00, p01, p10, p11)

p0 <- p000+p001+p010+p011
p1 <- p100+p101+p110+p111
p00.0 <- p000/p0
p00.1 <- p100/p1
p01.0 <- p001/p0
p01.1 <- p101/p1
p10.0 <- p010/p0
p10.1 <- p110/p1
p11.0 <- p011/p0
p11.1 <- p111/p1
pb <- c(p00.0, p00.1, p01.0, p01.1, p10.0, p10.1, p11.0, p11.1)

ba <- bounds(pa, "a")
print("a")
print(round(c(ba$l, ba$u), 2))
bb <- bounds(pb, "b")
print("b")
print(round(c(bb$l, bb$u), 2))

#---estimate bounds and 95% confidence intervals for c,d,e,f,g,h---

pY1 <- mean(VitD$death)
bS <- 0.5
pS1 <- c(0.1, 0.5, 0.9)
#using the full (S=0 + S=1) sample proportion of Z=1 as the true probability p(Z=1)
pZ1 <- mean(VitD$filaggrin)
B <- 100

scenarios <- c("c", "d", "e", "f", "g", "h")
ls <- length(scenarios)
l <- matrix(nrow=length(pS1), ncol=ls)
u <- matrix(nrow=length(pS1), ncol=ls)
ll <- matrix(nrow=length(pS1), ncol=ls)
ul <- matrix(nrow=length(pS1), ncol=ls)
lu <- matrix(nrow=length(pS1), ncol=ls)
uu <- matrix(nrow=length(pS1), ncol=ls)
for(i in 1:length(pS1)){
  a <- uniroot(f=function(a) plogis(a+bS)*pY1+plogis(a)*(1-pY1)-pS1[i],
    lower=-10, upper=10)$root
  #simulate outcome-dependent selection
  VitD$S <- rbinom(n, 1, plogis(a+bS*VitD$death))
  ss <- vector(length=8)
  j <- 0
  for(z in 0:1){
    for(x in 0:1){
      for(y in 0:1){
        j <- j+1
        ss[j] <- sum(VitD$filaggrin==z & VitD$vitdcat==x & VitD$death==y & VitD$S==1)
      }
    }
  }
  nss <- sum(ss)
  b <- sample.to.bounds(ss, pS1[i], pZ1, scenarios)
  l[i, ] <- b$l
  u[i, ] <- b$u
  pss <- ss/nss
  lB <- matrix(nrow=B, ncol=ls)
  uB <- matrix(nrow=B, ncol=ls)
  for(k in 1:B){
    ssB <- rmultinom(n=1, size=nss, prob=pss)
    b <- sample.to.bounds(ssB, pS1[i], pZ1, scenarios)
    lB[k, ] <- b$l
    uB[k, ] <- b$u
  }
  cil <- apply(X=lB, MARGIN=2, FUN=quantile, probs=c(0.025, 0.975), na.rm=TRUE)
  ciu <- apply(X=uB, MARGIN=2, FUN=quantile, probs=c(0.025, 0.975), na.rm=TRUE)
  ll[i, ] <- cil[1, ]
  ul[i, ] <- cil[2, ]
  lu[i, ] <- ciu[1, ]
  uu[i, ] <- ciu[2, ]
}

ml <- matrix(paste0(as.character(round(l, 2)), " (", as.character(round(ll, 2)),
  ",", as.character(round(ul, 2)), ")"), nrow=length(pS1), ncol=ls)
mu <- matrix(paste0(as.character(round(u, 2)), " (", as.character(round(lu, 2)),
  ",", as.character(round(uu, 2)), ")"), nrow=length(pS1), ncol=ls)
m <- rbind(ml, mu)
rownames(m) <- rep(paste0("p{S=1}=", pS1), 2)
colnames(m) <- scenarios

print("estimated bounds with CIs")
print(m)








