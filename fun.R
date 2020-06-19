
#-------------------------------------------------------------------------

#The functions bounds and bounds.inner are useful for users who want to
#compute the bounds, the remaining functions are mainly for the simulations
#in the paper.

#-------------------------------------------------------------------------

#The code contains three probability representations; p1, p2 and p3.
#p1 is vector of p(X|U,Z), p(Y|U,X), p(S|U,X,Y)
#p2 is vector of p(Z,X,Y|S=1)
#p3 is vector of probabilities of the form required for the particular bound,
#that is, p3 is vector of p(X,Y) for scenario a,
#a vector of p(X,Y|Z) for scenario b,
#a vector of p(X,Y|S=1) concatenated with p(S=1) for scenarios c, e, and g,
#and a vector of p(Z,X,Y|S=1) concatenated with {P(S=1), p(Z=1)} for
#scenarios d, f and h

#-------------------------------------------------------------------------

#Other common abbreviations:
#par = model parameters for the DAG
#pS1 = p(S=1)
#pU1 = p(U=1)
#pZ1 = p(Z=1)
#RD = causal risk difference
#sdU = standard deviation for effect of U on S
#sdX = standard deviation for effect of X on S

#-------------------------------------------------------------------------

#----------bootfun-------------
#Description: generates bootstrap replicates of bounds.
#Arguments:
# sample numeric vector containing the number of subjects for each level of (Z,X,Y), all with S=1.
# pS1 numeric scalar
# pZ1 numeric scalar
# scenarios string vector containg scenarios for which to produce bounds for
# B numeric scalar; number of bootstrap replicates
# pS1known logical indicating whether pS1 is condisdered known or not
#Value: list containing
# l numeric matrix; lower bounds for bootstrap replicates (rows) and scenarios (columns)
# u numeric matrix; upper bounds for bootstrap replicates (rows) and scenarios (columns)

bootfun <- function(sample, pS1, pZ1, scenarios, B, pS1known){
  ls <- length(scenarios)
  n <- sum(sample)
  psample <- sample/n
  l <- matrix(nrow=B, ncol=ls)
  u <- matrix(nrow=B, ncol=ls)
  for(k in 1:B){
    if(pS1known){
        sampleB <- rmultinom(n=1, size=n, prob=psample)
        pss <- pS1
    }
    else{
      fullsampleB <- rmultinom(n=1, size=n/pS1, prob=c(psample*pS1, 1-pS1))
      sampleB <- fullsampleB[1:8]
      pS1est <- sum(fullsampleB[1:8])/sum(fullsampleB)
      pss <- pS1est
    }
    b <- sample.to.bounds(sampleB, pss, pZ1, scenarios)
    l[k, ] <- b$l
    u[k, ] <- b$u
  }
  return(list(l=l, u=u))
}

#----------bounds-------------
#Description: computes the bounds. If scenario = "d", then switches to
#the bounds for scenario f if these are tighter.
#Arguments:
# p numeric vector containing p3
# scenario scalar string indicating which scenario to compute bounds for
#Value: list containing
# l numeric scalar; the lower bound
# il numeric scalar; the index for lower bound element
# u numeric scalar; the upper bound
# iu numeric scalar; the index for upper bound element
# switch.lower logical scalar indicating if switched lower d to lower f
# switch.upper logical scalar indicating if switched upper d to upper f

bounds <- function(p, scenario){
  b <- bounds.inner(p, scenario)
  b$switch.lower <- FALSE
  b$switch.upper <- FALSE
  if(scenario=="d"){
    bf <- bounds.inner(p, scenario="f")
    if(!is.na(bf$l) & !is.na(b$l))
      if(bf$l>b$l){
        b$l <- bf$l
        b$il <- bf$il
        b$switch.lower <- TRUE
      }
    if(!is.na(bf$u) & !is.na(b$u))
      if(bf$u<b$u){
        b$u <- bf$u
        b$iu <- bf$iu
        b$switch.upper <- TRUE
      }
  }
  return(b)
}

#----------bounds.inner-------------
#Description: computes the bounds, but does not switch d to f.
#Arguments:
# p numeric vector containing p3
# scenario scalar string indicating which scenario to compute bounds for
#Value: list containing
# l numeric scalar; the lower bound
# il numeric scalar; the index for lower bound element
# u numeric scalar; the upper bound
# iu numeric scalar; the index for upper bound element

bounds.inner <- function(p, scenario){

  if(scenario=="a"){
    p00 <- p[1]
    p01 <- p[2]
    p10 <- p[3]
    p11 <- p[4]

    l <- -(p01+p10)
    u <- 1-(p01+p10)
    il <- 1
    iu <- 1
  }

  if(scenario=="b"){
    p00.0 <- p[1]
    p00.1 <- p[2]
    p01.0 <- p[3]
    p01.1 <- p[4]
    p10.0 <- p[5]
    p10.1 <- p[6]
    p11.0 <- p[7]
    p11.1 <- p[8]

    l1 <- p11.1+p00.0-1
    l2 <- p11.0+p00.1-1
    l3 <- p11.0-p11.1-p01.1-p10.0-p01.0
    l4 <- p11.1-p11.0-p01.0-p10.1-p01.1
    l5 <- -p10.1-p01.1
    l6 <- -p10.0-p01.0
    l7 <- p00.1-p10.1-p01.1-p10.0-p00.0
    l8 <- p00.0-p10.0-p01.0-p10.1-p00.1
    lall <- c(l1, l2, l3, l4, l5, l6, l7, l8)
    l <- max(lall)
    il <- which.max(lall)

    u1 <- 1-p10.1-p01.0
    u2 <- 1-p10.0-p01.1
    u3 <- -p10.0+p10.1+p00.1+p11.0+p00.0
    u4 <- -p10.1+p11.1+p00.1+p10.0+p00.0
    u5 <- p11.1+p00.1
    u6 <- p11.0+p00.0
    u7 <- -p01.1+p11.1+p00.1+p11.0+p01.0
    u8 <- -p01.0+p11.0+p00.0+p11.1+p01.1
    uall <- c(u1, u2, u3, u4, u5, u6, u7, u8)
    u <- min(uall)
    iu <- which.min(uall)
  }

  if(scenario=="c"){
    p00.1 <- p[1]
    p01.1 <- p[2]
    p10.1 <- p[3]
    p11.1 <- p[4]
    pS1 <- p[5]
    A01 <- p10.1/p00.1
    A11 <- p11.1/p01.1

    r <- pS1
    l <- -(p01.1+p10.1)*r-max(1/(1+A11),A01/(1+A01))*(1-r)
    il <- 1
    u <- 1-(p01.1+p10.1)*r-min(1/(1+A11),A01/(1+A01))*(1-r)
    iu <- 1
  }

  if(scenario=="d"){
    p000.1 <- p[1]
    p001.1 <- p[2]
    p010.1 <- p[3]
    p011.1 <- p[4]
    p100.1 <- p[5]
    p101.1 <- p[6]
    p110.1 <- p[7]
    p111.1 <- p[8]
    pS1 <- p[9]
    pZ1<- p[10]

    pZ1.S1 <- p100.1+p101.1+p110.1+p111.1
    pS1.Z1 <- pZ1.S1*pS1/pZ1
    pS1.Z0 <- (1-pZ1.S1)*pS1/(1-pZ1)
    r1 <- pS1.Z1
    r0 <- pS1.Z0

    p00.01 <- p000.1/(1-pZ1.S1)
    p00.11 <- p100.1/pZ1.S1
    p01.01 <- p001.1/(1-pZ1.S1)
    p01.11 <- p101.1/pZ1.S1
    p10.01 <- p010.1/(1-pZ1.S1)
    p10.11 <- p110.1/pZ1.S1
    p11.01 <- p011.1/(1-pZ1.S1)
    p11.11 <- p111.1/pZ1.S1

    B001 <- p10.01/p00.01
    B011 <- p10.11/p00.11
    B101 <- p11.01/p01.01
    B111 <- p11.11/p01.11

    l1 <- p11.11*r1+p00.01*r0-1
    l2 <- p11.01*r0+p00.11*r1-1
    l3 <- (p11.01-p10.01-p01.01)*r0-(p11.11+p01.11)*r1-B001/(1+B001)*(1-r0)-(1-r1)
    l4 <- (p11.11-p10.11-p01.11)*r1-(p11.01+p01.01)*r0-B011/(1+B011)*(1-r1)-(1-r0)
    l5 <- -(p10.11+p01.11)*r1-max(1/(1+B111), B011/(1+B011))*(1-r1)
    l6 <- -(p10.01+p01.01)*r0-max(1/(1+B101), B001/(1+B001))*(1-r0)
    l7 <- (p00.11-p10.11-p01.11)*r1-(p10.01+p00.01)*r0-
      max(1/(1+B111), (B011-1)/(1+B011))*(1-r1)-(1-r0)
    l8 <- (p00.01-p10.01-p01.01)*r0-(p10.11+p00.11)*r1-
      max(1/(1+B101),-(1-B001)/(1+B001))*(1-r0)-(1-r1)
    lall <- c(l1, l2, l3, l4, l5, l6, l7, l8)
    l <- max(lall)
    il <- which.max(lall)

    u1 <- 1-p10.11*r1-p01.01*r0
    u2 <- 1-p10.01*r0-p01.11*r1
    u3 <- (p10.11+p00.11)*r1+(p11.01+p00.01-p10.01)*r0+
      +(1-r1)+max((1-B001)/(1+B001), B101/(1+B101))*(1-r0)
    u4 <- (p11.11+p00.11-p10.11)*r1+(p10.01+p00.01)*r0+
      max((1-B011)/(1+B011), B111/(1+B111))*(1-r1)+(1-r0)
    u5 <- (p11.11+p00.11)*r1+max(1/(1+B011), B111/(1+B111))*(1-r1)
    u6 <- (p11.01+p00.01)*r0+max(1/(1+B001), B101/(1+B101))*(1-r0)
    u7 <- (p11.11+p00.11-p01.11)*r1+(p11.01+p01.01)*r0+
      max(1/(1+B011), -(1-B111)/(1+B111))*(1-r1)+(1-r0)
    u8 <- (p11.11+p01.11)*r1+(p11.01+p00.01-p01.01)*r0+
      max(1/(1+B001), -(1-B101)/(1+B101))*(1-r0)+(1-r1)
    uall <- c(u1, u2, u3, u4, u5, u6, u7, u8)
    u <- min(uall)
    iu <- which.min(uall)
  }

  if(scenario=="e" | scenario=="g"){
    p00.1 <- p[1]
    p01.1 <- p[2]
    p10.1 <- p[3]
    p11.1 <- p[4]
    pS1 <- p[5]

    p001 <- p00.1*pS1
    p011 <- p01.1*pS1
    p101 <- p10.1*pS1
    p111 <- p11.1*pS1

  l <- (p111+p001)-1
  il <- 1
  u <- 1-(p011+p101)
  iu <- 1
  }

  if(scenario=="f"){
    pp000.1 <- p[1]
    pp001.1 <- p[2]
    pp010.1 <- p[3]
    pp011.1 <- p[4]
    pp100.1 <- p[5]
    pp101.1 <- p[6]
    pp110.1 <- p[7]
    pp111.1 <- p[8]
    pS1 <- p[9]
    pZ1<- p[10]

    p001.0 <- pp000.1*pS1/(1-pZ1)
    p001.1 <- pp100.1*pS1/pZ1
    p011.0 <- pp001.1*pS1/(1-pZ1)
    p011.1 <- pp101.1*pS1/pZ1
    p101.0 <- pp010.1*pS1/(1-pZ1)
    p101.1 <- pp110.1*pS1/pZ1
    p111.0 <- pp011.1*pS1/(1-pZ1)
    p111.1 <- pp111.1*pS1/pZ1

    l1 <- p001.1-p011.0-p111.0+2*p111.1-1
    l2 <- p001.1+p111.0-1
    l3 <- p001.0+p111.1-1
    l4 <- p001.1+p111.1-1
    l5 <- p001.0-p011.1+2*p111.0-p111.1-1
    l6 <- -p001.0+2*p001.1-p101.0+p111.1-1
    l7 <- p001.0+p111.0-1
    l8 <- 2*p001.0-p001.1-p101.1+p111.0-1
    l9 <- 2*p001.0-p001.1-p101.1-p011.1+2*p111.0-p111.1-1
    l10 <- -p001.0+2*p001.1-p101.0-p011.0-p111.0+2*p111.1-1
    lall <- c(l1, l2, l3, l4, l5, l6, l7, l8, l9, l10)
    l <- max(lall)
    il <- which.max(lall)

    u1 <- -p101.0-2*p011.0+p011.1+p111.1+1
    u2 <- -p101.1+p011.0-2*p011.1+p111.0+1
    u3 <- p001.1-2*p101.0+p101.1-p011.0+1
    u4 <- -p101.0-p011.1+1
    u5 <- -p101.1-p011.0+1
    u6 <- -p101.0-p011.0+1
    u7 <- p001.1-2*p101.0+p101.1-2*p011.0+p011.1+p111.1+1
    u8 <- -p101.1-p011.1+1
    u9 <- p001.0+p101.0-2*p101.1-p011.1+1
    u10 <- p001.0+p101.0-2*p101.1+p011.0-2*p011.1+p111.0+1
    uall <- c(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10)
    u <- min(uall)
    iu <- which.min(uall)
  }

  if(scenario=="h"){
    pp000.1 <- p[1]
    pp001.1 <- p[2]
    pp010.1 <- p[3]
    pp011.1 <- p[4]
    pp100.1 <- p[5]
    pp101.1 <- p[6]
    pp110.1 <- p[7]
    pp111.1 <- p[8]
    pS1 <- p[9]
    pZ1<- p[10]

    p001.0 <- pp000.1*pS1/(1-pZ1)
    p001.1 <- pp100.1*pS1/pZ1
    p011.0 <- pp001.1*pS1/(1-pZ1)
    p011.1 <- pp101.1*pS1/pZ1
    p101.0 <- pp010.1*pS1/(1-pZ1)
    p101.1 <- pp110.1*pS1/pZ1
    p111.0 <- pp011.1*pS1/(1-pZ1)
    p111.1 <- pp111.1*pS1/pZ1

    l1 <- p001.1+p111.1-1
    l2 <- p001.0+p111.1-1
    l3 <- p001.1+p111.0-1
    l4 <- p001.0+p111.0-1
    l5 <- 2*p001.1+p011.0+p111.0+p111.1-2
    l6 <- 2*p001.0+p011.1+p111.0+p111.1-2
    l7 <- p001.0+p001.1+p101.0+2*p111.1-2
    l8 <- p001.0+p001.1+p101.1+2*p111.0-2
    lall <- c(l1, l2, l3, l4, l5, l6, l7, l8)
    l <- max(lall)
    il <- which.max(lall)

    u1 <- -p101.0-p011.0+1
    u2 <- -p101.1-p011.0+1
    u3 <- -p101.0-p011.1+1
    u4 <- -p101.1-p011.1+1
    u5 <- -p001.0-p101.0-p101.1-2*p011.1+2
    u6 <- -p001.1-p101.0-p101.1-2*p011.0+2
    u7 <- -2*p101.0-p011.0-p011.1-p111.1+2
    u8 <- -2*p101.1-p011.0-p011.1-p111.0+2
    uall <- c(u1, u2, u3, u4, u5, u6, u7, u8)
    u <- min(uall)
    iu <- which.min(uall)
  }
  return(list(l=l, il=il, u=u, iu=iu))
}

#----------par.to.p1-------------
#Description: computes p1 from par
#Arguments:
# vector par
#Value:
# vector p1

par.to.p1 <- function(par){

  a0 <- par[1]
  aZ <- par[2]
  aU <- par[3]
  aUZ <- par[4]
  b0 <- par[5]
  bX <- par[6]
  bU <- par[7]
  bUX <- par[8]
  c0 <- par[9]
  cY <- par[10]
  cU <- par[11]
  cX <- par[12]

  pX1.U0Z0 <- plogis(a0)
  pX1.U0Z1 <- plogis(a0+aZ)
  pX1.U1Z0 <- plogis(a0+aU)
  pX1.U1Z1 <- plogis(a0+aZ+aU+aUZ)
  pY1.U0X0 <- plogis(b0)
  pY1.U0X1 <- plogis(b0+bX)
  pY1.U1X0 <- plogis(b0+bU)
  pY1.U1X1 <- plogis(b0+bX+bU+bUX)
  pS1.U0X0Y0 <- plogis(c0)
  pS1.U0X0Y1 <- plogis(c0+cY)
  pS1.U0X1Y0 <- plogis(c0+cX)
  pS1.U0X1Y1 <- plogis(c0+cY+cX)
  pS1.U1X0Y0 <- plogis(c0+cU)
  pS1.U1X0Y1 <- plogis(c0+cY+cU)
  pS1.U1X1Y0 <- plogis(c0+cX+cU)
  pS1.U1X1Y1 <- plogis(c0+cY+cX+cU)
  p1 <- c(pX1.U0Z0, pX1.U0Z1, pX1.U1Z0, pX1.U1Z1,
    pY1.U0X0, pY1.U0X1, pY1.U1X0, pY1.U1X1,
    pS1.U0X0Y0, pS1.U0X0Y1, pS1.U0X1Y0, pS1.U0X1Y1,
    pS1.U1X0Y0, pS1.U1X0Y1, pS1.U1X1Y0, pS1.U1X1Y1)

  return(p1)

}

#----------p1.to.p2-------------
#Description: computes p2 from p1
#Arguments:
# p1 numeric vector
#Value:
# p2 numeric vector

p1.to.p2 <- function(pZ1, pU1, p1){

  pX1.U0Z0 <- p1[1]
  pX1.U0Z1 <- p1[2]
  pX1.U1Z0 <- p1[3]
  pX1.U1Z1 <- p1[4]
  pY1.U0X0 <- p1[5]
  pY1.U0X1 <- p1[6]
  pY1.U1X0 <- p1[7]
  pY1.U1X1 <- p1[8]
  pS1.U0X0Y0 <- p1[9]
  pS1.U0X0Y1 <- p1[10]
  pS1.U0X1Y0 <- p1[11]
  pS1.U0X1Y1 <- p1[12]
  pS1.U1X0Y0 <- p1[13]
  pS1.U1X0Y1 <- p1[14]
  pS1.U1X1Y0 <- p1[15]
  pS1.U1X1Y1 <- p1[16]

  p0000 <- (1-pZ1)*((1-pU1)*(1-pX1.U0Z0)*(1-pY1.U0X0)*(1-pS1.U0X0Y0)+pU1*(1-pX1.U1Z0)*(1-pY1.U1X0)*(1-pS1.U1X0Y0))
  p0001 <- (1-pZ1)*((1-pU1)*(1-pX1.U0Z0)*(1-pY1.U0X0)*pS1.U0X0Y0+pU1*(1-pX1.U1Z0)*(1-pY1.U1X0)*pS1.U1X0Y0)
  p0010 <- (1-pZ1)*((1-pU1)*(1-pX1.U0Z0)*pY1.U0X0*(1-pS1.U0X0Y1)+pU1*(1-pX1.U1Z0)*pY1.U1X0*(1-pS1.U1X0Y1))
  p0011 <- (1-pZ1)*((1-pU1)*(1-pX1.U0Z0)*pY1.U0X0*pS1.U0X0Y1+pU1*(1-pX1.U1Z0)*pY1.U1X0*pS1.U1X0Y1)

  p0100 <- (1-pZ1)*((1-pU1)*pX1.U0Z0*(1-pY1.U0X1)*(1-pS1.U0X1Y0)+pU1*pX1.U1Z0*(1-pY1.U1X1)*(1-pS1.U1X1Y0))
  p0101 <- (1-pZ1)*((1-pU1)*pX1.U0Z0*(1-pY1.U0X1)*pS1.U0X1Y0+pU1*pX1.U1Z0*(1-pY1.U1X1)*pS1.U1X1Y0)
  p0110 <- (1-pZ1)*((1-pU1)*pX1.U0Z0*pY1.U0X1*(1-pS1.U0X1Y1)+pU1*pX1.U1Z0*pY1.U1X1*(1-pS1.U1X1Y1))
  p0111 <- (1-pZ1)*((1-pU1)*pX1.U0Z0*pY1.U0X1*pS1.U0X1Y1+pU1*pX1.U1Z0*pY1.U1X1*pS1.U1X1Y1)

  p1000 <- pZ1*((1-pU1)*(1-pX1.U0Z1)*(1-pY1.U0X0)*(1-pS1.U0X0Y0)+pU1*(1-pX1.U1Z1)*(1-pY1.U1X0)*(1-pS1.U1X0Y0))
  p1001 <- pZ1*((1-pU1)*(1-pX1.U0Z1)*(1-pY1.U0X0)*pS1.U0X0Y0+pU1*(1-pX1.U1Z1)*(1-pY1.U1X0)*pS1.U1X0Y0)
  p1010 <- pZ1*((1-pU1)*(1-pX1.U0Z1)*pY1.U0X0*(1-pS1.U0X0Y1)+pU1*(1-pX1.U1Z1)*pY1.U1X0*(1-pS1.U1X0Y1))
  p1011 <- pZ1*((1-pU1)*(1-pX1.U0Z1)*pY1.U0X0*pS1.U0X0Y1+pU1*(1-pX1.U1Z1)*pY1.U1X0*pS1.U1X0Y1)

  p1100 <- pZ1*((1-pU1)*pX1.U0Z1*(1-pY1.U0X1)*(1-pS1.U0X1Y0)+pU1*pX1.U1Z1*(1-pY1.U1X1)*(1-pS1.U1X1Y0))
  p1101 <- pZ1*((1-pU1)*pX1.U0Z1*(1-pY1.U0X1)*pS1.U0X1Y0+pU1*pX1.U1Z1*(1-pY1.U1X1)*pS1.U1X1Y0)
  p1110 <- pZ1*((1-pU1)*pX1.U0Z1*pY1.U0X1*(1-pS1.U0X1Y1)+pU1*pX1.U1Z1*pY1.U1X1*(1-pS1.U1X1Y1))
  p1111 <- pZ1*((1-pU1)*pX1.U0Z1*pY1.U0X1*pS1.U0X1Y1+pU1*pX1.U1Z1*pY1.U1X1*pS1.U1X1Y1)

  p2 <- c(p0000, p0001, p0010, p0011,
    p0100, p0101, p0110, p0111,
    p1000, p1001, p1010, p1011,
    p1100, p1101, p1110, p1111)

  return(p2)

}

#----------p2.to.p3-------------
#Description: computes p3 from p2
#Arguments:
# p2 numeric vector
#Value: list containing elements pa,...,ph corresponding to p3 for each scenario

p2.to.p3 <- function(p2){

  p0000 <- p2[1]
  p0001 <- p2[2]
  p0010 <- p2[3]
  p0011 <- p2[4]
  p0100 <- p2[5]
  p0101 <- p2[6]
  p0110 <- p2[7]
  p0111 <- p2[8]
  p1000 <- p2[9]
  p1001 <- p2[10]
  p1010 <- p2[11]
  p1011 <- p2[12]
  p1100 <- p2[13]
  p1101 <- p2[14]
  p1110 <- p2[15]
  p1111 <- p2[16]

  p00 <- p0000+p0001+p1000+p1001
  p01 <- p0010+p0011+p1010+p1011
  p10 <- p0100+p0101+p1100+p1101
  p11 <- p0110+p0111+p1110+p1111
  pa <- c(p00, p01, p10, p11)

  p0 <- p0000+p0001+p0010+p0011+p0100+p0101+p0110+p0111
  p1 <- p1000+p1001+p1010+p1011+p1100+p1101+p1110+p1111
  p00.0 <- (p0000+p0001)/p0
  p00.1 <- (p1000+p1001)/p1
  p01.0 <- (p0010+p0011)/p0
  p01.1 <- (p1010+p1011)/p1
  p10.0 <- (p0100+p0101)/p0
  p10.1 <- (p1100+p1101)/p1
  p11.0 <- (p0110+p0111)/p0
  p11.1 <- (p1110+p1111)/p1
  pb <- c(p00.0, p00.1, p01.0, p01.1, p10.0, p10.1, p11.0, p11.1)

  pS1 <- p0001+p0011+p0101+p0111+p1001+p1011+p1101+p1111
  p00.1 <- (p0001+p1001)/pS1
  p01.1 <- (p0011+p1011)/pS1
  p10.1 <- (p0101+p1101)/pS1
  p11.1 <- (p0111+p1111)/pS1
  pc <- c(p00.1, p01.1, p10.1, p11.1, pS1)

  pS1 <- p0001+p0011+p0101+p0111+p1001+p1011+p1101+p1111
  pZ1 <- p1000+p1001+p1010+p1011+p1100+p1101+p1110+p1111
  p000.1 <- p0001/pS1
  p001.1 <- p0011/pS1
  p010.1 <- p0101/pS1
  p011.1 <- p0111/pS1
  p100.1 <- p1001/pS1
  p101.1 <- p1011/pS1
  p110.1 <- p1101/pS1
  p111.1 <- p1111/pS1
  pd <- c(p000.1, p001.1, p010.1, p011.1, p100.1, p101.1, p110.1, p111.1, pS1, pZ1)

  return(list(pa=pa, pb=pb, pc=pc, pd=pd, pe=pc, pf=pd, pg=pc, ph=pd))

}

#----------p1.to.RD-------------
#Description: computes RD from p1
#Arguments:
# pU1 numeric scalar
# p1 numeric vector
#Value:
# RD numeric scalar

p1.to.RD <- function(pU1, p1){

  pY1.U0X0 <- p1[5]
  pY1.U0X1 <- p1[6]
  pY1.U1X0 <- p1[7]
  pY1.U1X1 <- p1[8]
  (1-pU1)*(pY1.U0X1-pY1.U0X0)+pU1*(pY1.U1X1-pY1.U1X0)

}

#----------p2RDfun-------------
#Description: generates par randomly, then computes p2 and RD from par
#Arguments:
# sdU numeric scalar
# sdX numeric scalar
#Value: list containing
# p2 numeric vector
# RD numeric scalar

p2RDfun <- function(sdU, sdX){

  pZ1 <- runif(1)
  pU1 <- runif(1)
  a0 <- rnorm(1, sd=5)
  aZ <- rnorm(1, sd=5)
  aU <- rnorm(1, sd=5)
  aUZ <- rnorm(1, sd=5)
  b0 <- rnorm(1, sd=5)
  bX <- rnorm(1, sd=5)
  bU <- rnorm(1, sd=5)
  bUX <- rnorm(1, sd=5)
  c0 <- rnorm(1, sd=5)
  cY <- rnorm(1, sd=5)
  cU <- rnorm(1, sd=sdU)
  cX <- rnorm(1, sd=sdX)

  p1 <- par.to.p1(c(a0, aZ, aU, aUZ, b0, bX, bU, bUX, c0, cY, cU, cX))
  p2 <- p1.to.p2(pZ1, pU1, p1)
  RD <- p1.to.RD(pU1, p1)

  return(list(p2=p2, RD=RD))

}

#----------sample.to.bounds-------------
#Description: converts a sample to bounds.
#Arguments:
# sample numeric vector containing the number of subjects for each level of (Z,X,Y), all with S=1.
# pS1 numeric scalar
# pZ1 numeric scalar
# scenarios string vector containg scenarios for which to produce bounds for
#Value: list containing
# l numeric vector containing lower bounds
# u numeric vector containing upper bounds

sample.to.bounds <- function(sample, pS1, pZ1, scenarios){
  ls <- length(scenarios)
  p <- sample/sum(sample)
  p000.1 <- p[1]
  p001.1 <- p[2]
  p010.1 <- p[3]
  p011.1 <- p[4]
  p100.1 <- p[5]
  p101.1 <- p[6]
  p110.1 <- p[7]
  p111.1 <- p[8]
  p00.1 <- p000.1+p100.1
  p01.1 <- p001.1+p101.1
  p10.1 <- p010.1+p110.1
  p11.1 <- p011.1+p111.1
  pc <- pe <- pg <- c(p00.1, p01.1, p10.1, p11.1, pS1)
  pd <- pf <- ph <- c(p, pS1, pZ1)
  l <- vector(length=ls)
  u <- vector(length=ls)
  for(i in 1:ls){
    b <- bounds(p=get(paste0("p", scenarios[i])), scenario=scenarios[i])
    l[i] <- b$l
    u[i] <- b$u
  }
  return(list(l=l, u=u))
}

#----------simfun1-------------
#Description: generates distributions and computes bounds and causal risk
#differences from these for the simulation in section 6.1.
#Arguments:
# N numeric scalar, number of bounds
# scenarios string vector containg scenarios for which to produce bounds for
# sdX numeric scalar
# sdU numeric scalar
#Value: list containing
# l, u numeric matricies containing lower and upper bounds, respectively,
# for each distribution (rows) and scenario (columns)
# RD numeric vector containing causal risk difference for each distribution.

simfun1 <- function(N, scenarios, sdX, sdU){
  ls <- length(scenarios)
  RD <- vector(length=N)
  l <- matrix(nrow=N, ncol=ls)
  colnames(l) <- scenarios
  u <- matrix(nrow=N, ncol=ls)
  colnames(u) <- scenarios
  for(i in 1:N){
    p2RD <- p2RDfun(sdU=sdU, sdX=sdX)
    RD[i] <- p2RD$RD
    p2 <- p2RD$p2
    p3 <- p2.to.p3(p2)
    for(j in 1:ls){
      b <- bounds(p=p3[[paste0("p",scenarios[j])]], scenario=scenarios[j])
      l[i, j] <- b$l
      u[i, j] <- b$u
    }
  }
  return(list(l=l, u=u, RD=RD))
}

#----------simfun2-------------
#Description: simulates samples and bootstraps for the simulation in section 6.2.
#Arguments:
# p numeric vector; true p3
# l numeric vector; true lower bounds
# u numeric vector; true upper bounds
# n numeric vector; desired sample sizes, i.e. observed subjects with S=1.
# If pS1 is unknown, then the number of subjects with S=1 varies from
# sample to sample in simulation; the "full" (S=1+S=0) sample size is chosen
# to give n number of subjects with S=1 in expectation.
# B numeric scalar; number of bootstrap replicates
# N numeric scalar; number of samples
# scenarios string vector containg scenarios for which to produce bounds for
# pS1 numeric scalar
# pZ1 numeric scalar
# pS1known logical indicating whether pS1 is condisdered known or not
#Value: list containing
# biasl, biasu numeric matricies containing mean (over samples) bias in estimated
#lower and upper bounds, respectively, for each n (rows) and scenario (columns)
# sdl, sdu numeric matricies containing sd (over samples) of estimated lower
# and upper bounds, respectively, for each n (rows) and scenario (columns)
# coverprobl, coverprobu numeric matricies containing coverage probabilities
# (over samples) of bootstgrap CI for estimated lower and upper bounds,
# respectively, for each n (rows) and scenario (columns)

simfun2 <- function(p, l, u, n, B, N, scenarios, pS1, pZ1, pS1known){
  ls <- length(scenarios)
  biasl <- matrix(nrow=length(n), ncol=ls)
  biasu <- matrix(nrow=length(n), ncol=ls)
  sdl <- matrix(nrow=length(n), ncol=ls)
  sdu <- matrix(nrow=length(n), ncol=ls)
  coverprobl <- matrix(nrow=length(n), ncol=ls)
  coverprobu <- matrix(nrow=length(n), ncol=ls)
  colnames(biasl) <- scenarios
  colnames(biasu) <- scenarios
  colnames(sdl) <- scenarios
  colnames(sdu) <- scenarios
  colnames(coverprobl) <- scenarios
  colnames(coverprobu) <- scenarios
  for(i in 1:length(n)){
    coverl <- matrix(nrow=N, ncol=ls)
    coveru <- matrix(nrow=N, ncol=ls)
    bl <- matrix(nrow=N, ncol=ls)
    bu <- matrix(nrow=N, ncol=ls)
    for(j in 1:N){
      if(pS1known){
        sample <- rmultinom(n=1, size=n[i], prob=p)
        pss <- pS1
      }
      else{
        fullsample <- rmultinom(n=1, size=n[i]/pS1, prob=c(p*pS1, 1-pS1))
        sample <- fullsample[1:8]
        pS1est <- sum(fullsample[1:8])/sum(fullsample)
        pss <- pS1est
      }
      b <- sample.to.bounds(sample, pss, pZ1, scenarios)
      boot <- bootfun(sample, pss, pZ1, scenarios, B, pS1known)
      bl[j, ] <- b$l
      bu[j, ] <- b$u
      lB <- boot$l
      uB <- boot$u
      cil <- apply(X=lB, MARGIN=2, FUN=quantile, probs=c(0.025, 0.975),
        na.rm=TRUE)
      coverl[j, ] <- l>cil[1, ] & l<cil[2, ]
      ciu <- apply(X=uB, MARGIN=2, FUN=quantile, probs=c(0.025, 0.975),
        na.rm=TRUE)
      coveru[j, ] <- u>ciu[1, ] & u<ciu[2, ]
    }
    biasl[i, ] <- l-colMeans(bl)
    biasu[i, ] <- u-colMeans(bu)
    sdl[i, ] <- apply(X=bl, MARGIN=2, FUN=sd, na.rm=TRUE)
    sdu[i, ] <- apply(X=bu, MARGIN=2, FUN=sd, na.rm=TRUE)
    coverprobl[i, ] <- colMeans(coverl)
    coverprobu[i, ] <- colMeans(coveru)
  }
  return(list(biasl=biasl, biasu=biasu, sdl=sdl, sdu=sdu,
    coverprobl=coverprobl, coverprobu=coverprobu))
}

