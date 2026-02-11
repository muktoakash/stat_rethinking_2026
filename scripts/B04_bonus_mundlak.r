# endogenous group confound example

library(rethinking)

set.seed(8672)

sim_dat <- function(
    N_groups = 30,
    N_id = 200,
    a0 = (-2),
    bZY = (-1),
    bXY = 0
) {
    g <- sample(1:N_groups,size=N_id,replace=TRUE) # sample into groups
    Ug <- rnorm(N_groups,1.5) # group confounds
    X <- rnorm(N_id, Ug[g] ) # individual varying trait
    Z <- rnorm(N_groups) # group varying trait (observed)
    Y <- rbern(N_id, p=inv_logit( a0 + bXY*X + Ug[g] + bZY*Z[g] ) )
    dat <- list(g=g,Ug=Ug,X=X,Z=Z,Y=Y,Ng=N_groups,N_id=N_id,a0=a0,bZY=bZY,bXY=bXY)
    return(dat)
}

dat <- sim_dat( N_groups=30 , N_id=2000 , bZY=1 , bXY=0 )

xbar <- sapply( 1:dat$Ng , function(j) mean(dat$X[dat$g==j]) )
dat$Xbar <- xbar
dat$Xbar_g <- xbar[dat$g]
dat$Z_g <- dat$Z[dat$g]

table(dat$g)

if (FALSE) {

plot( dat$g , dat$X , xlab="group" , ylab="X" , col=2 )

plot( dat$g , dat$X - dat$Xbar_g , xlab="group" , ylab="X - Xbar[g]" , col=2 )

plot( dat$X , jitter(dat$Y) , col=2 , xlab="X" , ylab="Y" )

plot( dat$X - dat$Xbar_g , jitter(dat$Y) , col=2 , xlab="X-Xbar[g]" , ylab="Y" )

}

# confounded by correlation
#precis(glm(Y~X,family=binomial,data=list(Y=dat$Y,X=dat$X)),2)

# fixed effects
# X deconfounded, but Z unidentified now!
#precis(glm(Y~X+as.factor(g),family=binomial,data=list(Y=dat$Y,X=dat$X,g=dat$g)),pars=c("X","Z"),2)
#precis(glm(Y~X+Z+as.factor(g),family=binomial,data=list(Y=dat$Y,X=dat$X,g=dat$g,Z=dat$Z[dat$g])),pars=c("X","Z"),2)

# naive model
m0 <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a + bxy*X + bzy*Z_g,
        a ~ dnorm(0,10),
        c(bxy,bzy) ~ dnorm(0,1)
    ) , data=dat , chains=4 , cores=4 , refresh=0 )

# fixed effects
mf <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a[g] + bxy*X + bzy*Z_g,
        a[g] ~ dnorm(0,10),
        c(bxy,bzy) ~ dnorm(0,1)
    ) , data=dat , chains=4 , cores=4 , refresh=0 )

# fixed effects by subtracting Xbar
mf2 <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a + bxy*(X-Xbar_g) + bzy*Z_g,
        a ~ dnorm(0,10),
        c(bxy,bzy) ~ dnorm(0,1)
    ) , data=dat , chains=4 , cores=4 , refresh=0 )

# random effects
mr <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a[g] + bxy*X + bzy*Z_g,
        transpars> vector[Ng]:a <<- abar + z*tau,
        z[g] ~ dnorm(0,1),
        c(bxy,bzy) ~ dnorm(0,1),
        abar ~ dnorm(0,1),
        tau ~ dexp(1)
    ) , data=dat , chains=4 , cores=4 , sample=TRUE , refresh=0 )

# random effects + Xbar
# The Mundlak Machine

mrx <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- a[g] + bxy*X + buy*Xbar_g + bzy*Z_g ,
        transpars> vector[Ng]:a <<- abar + z*tau,
        z[g] ~ dnorm(0,1),
        c(bxy,buy,bzy) ~ dnorm(0,1),
        abar ~ dnorm(0,1),
        tau ~ dexp(1)
    ) , data=dat , chains=4 , cores=4 , sample=TRUE , refresh=0 )

# random effects + latent U
# The Latent Mundlak Machine
mru <- ulam(
    alist(
        # Y model
        Y ~ bernoulli(p),
        logit(p) <- a[g] + bxy*X + bzy*Z_g + buy*u[g],
        transpars> vector[Ng]:a <<- abar + z*tau,
        # X model
        X ~ normal(mu,sigma),
        mu <- aX + bux*u[g],
        vector[Ng]:u ~ normal(0,1),
        # priors
        z[g] ~ dnorm(0,1),
        c(aX,bxy,buy,bzy) ~ dnorm(0,1),
        bux ~ dexp(1),
        abar ~ dnorm(0,1),
        tau ~ dexp(1),
        sigma ~ dexp(1)
    ) , data=dat , chains=4 , cores=4 , sample=TRUE , refresh=0 )

precis(m0,pars=c("bxy","bzy"))
precis(mf,pars=c("bxy","bzy"))
precis(mf2,pars=c("bxy","bzy"))
precis(mr,pars=c("bxy","bzy"))
precis(mrx,pars=c("bxy","bzy"))
precis(mru,pars=c("bxy","bzy"))

# compute accuracy of intervention - the real estimand
acc_doX <- function(model) {
    post <- extract.samples(model)
    bXY_true <- model@data$bXY
    abserr <- abs( bXY_true - post$bxy )
    # predict intervention in each group
    dat1 <- model@data
    dat1$X <- dat1$X+1
    l0 <- link(model)
    l1 <- link(model,data=dat1)
    if ( length(l0)==2 )
        diff_pred <- l1$p - l0$p
    else 
        diff_pred <- l1-l0
    # true effect
    #logit_p0 <- with( dat1 , a0 + bXY*X + Ug[g] + bZY*Z[g] )
    #logit_p1 <- with( dat1 , a0 + bXY*(X+1) + Ug[g] + bZY*Z[g] )
    #diff_true <- inv_logit(logit_p1) - inv_logit(logit_p0)
    # if true bXY is zero just need mean(abs(diff_pred))
    mean(abs(diff_pred))
}

acc_doX(m0)
acc_doX(mf)
acc_doX(mr)
acc_doX(mrx)
acc_doX(mru)

blank2(w=1.35,h=0.7)

# density plots
# bxy
post <- extract.samples(mf)
#dens(post$bxy,lwd=3,col=1,xlab="b_XY",ylim=c(0,11),xlim=c(-0.2,0.6))
dens(post$bxy,lwd=3,col=1,xlab="b_XY",ylim=c(0,2.5),xlim=c(-1.2,0.7))
abline(v=dat$bXY,lty=2)

post <- extract.samples(m0)
dens(post$bxy,lwd=3,col=grau(),add=TRUE)

post <- extract.samples(mr)
dens(post$bxy,lwd=3,col=2,add=TRUE)

post <- extract.samples(mrx)
dens(post$bxy,lwd=3,col=4,add=TRUE)

post <- extract.samples(mru)
dens(post$bxy,lwd=8,col="white",add=TRUE)
dens(post$bxy,lwd=4,col=3,add=TRUE)

# bzy
post <- extract.samples(mf)
dens(post$bzy,lwd=3,col=1,xlab="b_ZY",ylim=c(0,7))
abline(v=dat$bZY,lty=2)

post <- extract.samples(m0)
dens(post$bzy,lwd=3,col=grau(),add=TRUE)

post <- extract.samples(mr)
dens(post$bzy,lwd=3,col=2,add=TRUE)

post <- extract.samples(mrx)
dens(post$bzy,lwd=3,col=4,add=TRUE)

post <- extract.samples(mru)
dens(post$bzy,lwd=8,col="white",add=TRUE)
dens(post$bzy,lwd=4,col=3,add=TRUE)


##########
# show better estimates of intercepts

af <- coef(mf)[1:N_groups]
ar <- coef(mr)[1:N_groups]

plot( af , col=4 )
points( 1:N_groups , ar , col=2 )
points( 1:N_groups , a0+Ug , col=1 )


# treatment effect in each group now
# counterfactual increase of X at individual level, stratified by each group

# fixed estimates
pf0 <- link(mf,data=list(g=1:N_groups,X=rep(0,N_groups)))
pf1 <- link(mf,data=list(g=1:N_groups,X=rep(1,N_groups)))
cf <- apply( pf1 - pf0 , 2 , mean )

# random estimates
pr0 <- link(mr,data=list(g=1:N_groups,X=rep(0,N_groups)))
pr1 <- link(mr,data=list(g=1:N_groups,X=rep(1,N_groups)))
cr <- apply( pr1 - pr0 , 2 , mean )

# true
ctrue <- inv_logit( a0 + Ug + 1 ) - inv_logit( a0 + Ug )

plot( ctrue , ylim=c(0,0.3) )
points( 1:N_groups , cf , col=4 )
points( 1:N_groups , cr , col=2 )

plot( cf - ctrue , col=4 , ylim=c(-0.1,0.1) )
points( 1:N_groups , cr - ctrue , col=2 )
abline(h=0,lty=2)

mean((cf - ctrue)^2)
mean((cr - ctrue)^2)
