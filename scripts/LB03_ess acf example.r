library(rethinking)
data(bangladesh)
d <- bangladesh

dat <- list(
    C = d$use.contraception,
    D = as.integer(d$district),
    U = ifelse(d$urban==1,1,0) )

# total U
mCDU <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- a[D] + b[D]*U,
        vector[61]:a ~ normal(abar,sigma),
        vector[61]:b ~ normal(bbar,tau),
        c(abar,bbar) ~ normal(0,1),
        c(sigma,tau) ~ exponential(1)
    ) , data=dat , chains=4 , cores=4 )





# ess example
D <- 1000
mess <- ulam(
    alist(
        vector[D]:theta ~ normal(0,1),
        save> vector[D]:theta_sq <<- theta^2
    ), 
    data=list(D=D) , chains=4 , cores=4, refresh=0 )

trankplot(mess,pars=c("theta[1]","theta_sq[1]"))

draws <- attr(mess,"cstanfit")$draws()

i <- 2
acf(draws[,1,i+1],main="",lwd=5,col=2)
mtext("theta[1]")

acf(draws[,1,i+1+D],main="",lwd=5,col=2)
mtext("theta_sq[1]")

pp <- precis(mess,2)

essb <- pp['ess_bulk'][,1]
esst <- pp['ess_tail'][,1]

blank2()

plot( essb[1:D] , esst[1:D] , xlab="theta ess_bulk" , ylab="theta ess_tail" , col=2 )

plot( essb[(D+1):(2*D)] , esst[(D+1):(2*D)] , xlab="theta^2 ess_bulk" , ylab="theta^2 ess_tail" , col=2 )

plot( essb[1:D] , essb[(D+1):(2*D)] , xlab="theta ess_bulk" , ylab="theta^2 ess_bulk" , col=2  )

plot( esst[1:D] , esst[(D+1):(2*D)] , xlab="theta ess_tail" , ylab="theta^2 ess_tail" , col=2 )



