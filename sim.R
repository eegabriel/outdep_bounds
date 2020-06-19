source("fun.R")

set.seed(1)

#----------Numeric example for Section 4.3----------

pU1 <- 0.9
pZ1 <- 0.7
pX1.U0Z0 <- 0.3
pX1.U0Z1 <- 0.3
pX1.U1Z0 <- 0.9
pX1.U1Z1 <- 0.1
pY1.U0X0 <- 0.3
pY1.U0X1 <- 0.4
pY1.U1X0 <- 0.1
pY1.U1X1 <- 0.8
pS1.U0X0Y0 <- 0.9
pS1.U0X0Y1 <- 0.2
pS1.U0X1Y0 <- 0.9
pS1.U0X1Y1 <- 0.2
pS1.U1X0Y0 <- 0.9
pS1.U1X0Y1 <- 0.2
pS1.U1X1Y0 <- 0.9
pS1.U1X1Y1 <- 0.2
p1 <- c(pX1.U0Z0, pX1.U0Z1, pX1.U1Z0, pX1.U1Z1,
  pY1.U0X0, pY1.U0X1, pY1.U1X0, pY1.U1X1,
  pS1.U0X0Y0, pS1.U0X0Y1, pS1.U0X1Y0, pS1.U0X1Y1,
  pS1.U1X0Y0, pS1.U1X0Y1, pS1.U1X1Y0, pS1.U1X1Y1)
RD <- p1.to.RD(pU1, p1)
p2 <- p1.to.p2(pZ1, pU1, p1)
p3 <- p2.to.p3(p2)
bd <- bounds.inner(p=p3$pd, scenario="d")
bf <- bounds.inner(p=p3$pf, scenario="f")

print("---Numeric example---")
print("true causal risk difference")
print(RD)
print("p(Z,X,Y|S=1)")
print(p2)
print("bounds under scenario d")
print(bd)
print("bounds under scenario f")
print(bf)

#----------Figure 2----------

scenarios <- c("a", "b", "c", "d", "e", "f", "g", "h")
ls <- length(scenarios)
N <- 100000
tmp <- simfun1(N, scenarios, sdX=0, sdU=0)
width <- as.vector(tmp$u-tmp$l)
ss <- rep(scenarios, each=N)

pdf("width0.pdf", height = 3.5, width = 3.5)
par(mar = c(4, 4, .25, .5), family = "serif")
boxplot(width~ss, xlab="scenario", ylab="width", outline=FALSE)
dev.off()

RDcut <- seq(0, 1, 0.1)
p <- matrix(nrow=length(RDcut)-1, ncol=ls)
for(i in 1:(length(RDcut)-1)){
  incl <- abs(tmp$RD)>RDcut[i] & abs(tmp$RD)<RDcut[i+1]
  p[i, ] <- colMeans(tmp$l[incl, ]>0 | tmp$u[incl, ]<0)
}
pch <- c(0,22,1,21,2,24,6,25)
bg <- rep(c("white", "black"), 4)

pdf("exclude.pdf", height = 3.5, width = 3.5)
par(mar = c(4, 4, .25, .5), family = "serif")
matplot(RDcut[2:length(RDcut)], p, type="l", col="black", lty="solid",
  xlab=expression(paste("|"*theta*"|")), ylab=expression(paste("proportion excluding  ",theta,"=0")),
  xaxt="n")
axis(1, at=RDcut[-1], labels=paste0(RDcut[1:(length(RDcut)-1)],"-",RDcut[2:length(RDcut)]),
     cex.axis=.75, las = 3)
matpoints(RDcut[2:length(RDcut)], p, pch=pch, col="black", bg=bg)
legend("topleft", pch=pch, legend=scenarios, col="black", pt.bg=bg)
dev.off()

#----------Figure 3----------


lty <- "solid"
pch <- c(0,22,1,21,2,24,6,25)
bg <- rep(c("white", "black"), 4)
scenarios <- c("a", "b", "c", "d", "e", "f", "g", "h")
ls <- length(scenarios)
N <- 100000

#Varying the effect of U on S.
sdU <- seq(0, 10, 1)
lsd <- length(sdU)
pwrong <- matrix(nrow=lsd, ncol=ls)
width <- matrix(nrow=lsd, ncol=ls)
rownames(pwrong) <- sdU
colnames(pwrong) <- scenarios
for(k in 1:lsd){
  tmp <- simfun1(N, scenarios, sdX=0, sdU=sdU[k])
  pwrong[k, ] <- colMeans(tmp$u<tmp$RD | tmp$l>tmp$RD)
  width[k, ] <- apply(X=tmp$u-tmp$l, MARGIN=2, FUN=median)
}

pdf("width.pdf", width = 7, height = 5.5)
layout(matrix(c(1, 2, 5, 3, 4, 5), nrow = 2, ncol = 3, byrow = TRUE),
       widths = c(.45, .45, .1), heights = c(1, 1))
par(mar = c(4, 4, .5, .5), family = "serif", mgp = c(2, 1, 0))
matplot(sdU, width, type="l", col="black", lty=lty,
  xlab=expression(sigma[U]), ylab="median width")
matpoints(sdU, width, pch=pch, col="black", bg=bg)
matplot(sdU, pwrong, type="l", col="black", lty=lty,
  xlab=expression(sigma[U]), ylab="proportion violated")
matpoints(sdU, pwrong, pch=pch, col="black", bg=bg)
#legend("topleft", pch=pch, legend=scenarios, col="black", pt.bg=bg)

#Varying the effect of X on S.
sdX <- seq(0, 10, 1)
lsd <- length(sdX)
pwrong <- matrix(nrow=lsd, ncol=ls)
width <- matrix(nrow=lsd, ncol=ls)
rownames(pwrong) <- sdX
colnames(pwrong) <- scenarios
for(k in 1:lsd){
  tmp <- simfun1(N, scenarios, sdX=sdX[k], sdU=0)
  pwrong[k, ] <- colMeans(tmp$u<tmp$RD | tmp$l>tmp$RD)
  width[k, ] <- apply(X=tmp$u-tmp$l, MARGIN=2, FUN=median)
}

matplot(sdU, width, type="l", col="black", lty=lty,
  xlab=expression(sigma[X]), ylab="median width")
matpoints(sdU, width, pch=pch, col="black", bg=bg)
matplot(sdU, pwrong, type="l", col="black", lty=lty,
  xlab=expression(sigma[X]), ylab="proportion violated")
matpoints(sdU, pwrong, pch=pch, col="black", bg=bg)

par(mar = c(1, 1, 1, 1))
plot(1~1, type = "n", axes = FALSE)
legend("center", pch=pch, legend=scenarios, col="black", pt.bg=bg)
dev.off()

#----------Tables 3 and 4----------

#Generate true bounds and true causal risk difference.
pU1 <- 0.5
pZ1 <- 0.5
a0 <- -1
aZ <- 0.5
aU <- 0.5
aUZ <- 0
b0 <- -1
bX <- 0.5
bU <- 0.5
bUX <- 0
c0 <- -1
cY <- 0.5
cU <- 0
cX <- 0
p1 <- par.to.p1(c(a0, aZ, aU, aUZ, b0, bX, bU, bUX, c0, cY, cU, cX))
p2 <- p1.to.p2(pZ1, pU1, p1)
p3 <- p2.to.p3(p2)
RD <- p1.to.RD(pU1, p1)
scenarios <- c("c", "d", "e", "f", "g", "h")
ls <- length(scenarios)
l <- vector(length=ls)
names(l) <- scenarios
u <- vector(length=ls)
names(u) <- scenarios
for(j in 1:ls){
  b <- bounds(p=p3[[paste0("p",scenarios[j])]], scenario=scenarios[j])
  l[j] <- b$l
  u[j] <- b$u
}

#Simulate samples and bootstrap when pS1 is known and unknown from
#true distribution.
p <- p3$pd[1:8]
pS1 <- p3$pd[9]
B <- 1000
N <- 1000
n <- c(100, 1000)
sim.pS1known <- simfun2(p, l, u, n, B, N, scenarios, pS1, pZ1, pS1known=TRUE)
sim.pS1unknown <- simfun2(p, l, u, n, B, N, scenarios, pS1, pZ1, pS1known=FALSE)

bias <- rbind(sim.pS1known$biasl, sim.pS1known$biasu,
  sim.pS1unknown$biasl, sim.pS1unknown$biasu)
colnames(bias) <- scenarios
rownames(bias) <- rep(paste0("n=", n), 4)

sd.cover <- rbind(cbind(sim.pS1known$sdl, sim.pS1known$coverprobl),
  cbind(sim.pS1known$sdu, sim.pS1known$coverprobu),
  cbind(sim.pS1unknown$sdl, sim.pS1unknown$coverprobl),
  cbind(sim.pS1unknown$sdu, sim.pS1unknown$coverprobu))
colnames(sd.cover) <- rep(scenarios, 2)
rownames(sd.cover) <- rep(paste0("n=", n), 4)



print("---Simulation for estimated bounds---")
print("true causal risk difference")
print(round(RD, 2))
print("true bounds")
print(round(rbind(l,u), 2))
print("bias in estimated bounds")
print(round(bias, 2))
print("sd of estimated bounds and coverage probability of bootstrap CIs")
print(round(sd.cover, 2))











