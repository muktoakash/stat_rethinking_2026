sim_weight <- function(H,b,sd) {
    U <- rnorm( length(H) , 0 , sd )
    W <- b*H + U
    return(W)
}

# simulate a sample of 10 people
set.seed(93)
H <- runif(10,130,170)
W <- sim_weight(H,b=0.5,sd=5)

# run the model
library(rethinking)
m3.1 <- quap(
    alist(
        W ~ dnorm(mu,sigma),
        mu <- a + b*H,
        a ~ dnorm(0,10),
        b ~ dunif(0,1),
        sigma ~ dunif(0,10)
    ) , data=list(W=W,H=H) )

# summary
precis( m3.1 )

f <- function(n=10,b=0.5) {
    H <- runif(n,130,170)
    W <- sim_weight(H,b=b,sd=5)
    m <- quap(
        alist(
            W ~ dnorm(mu,sigma),
            mu <- a + b*H,
            a ~ dnorm(0,10),
            b ~ dunif(0,1),
            sigma ~ dunif(0,10)
        ) , data=list(W=W,H=H) )
    return(precis(m)['b',])
}

blank2(h=2)

b_true <- 0.5
b_ests <- replicate( 100 , f(n=100,b=b_true) )

plot(NULL,xlim=c(0,1),ylim=c(1,dim(b_ests)[3]),xlab="b",ylab="")
for ( i in 1:dim(b_ests)[3] ) {
    points( b_ests[1,1,i] , i , col=2 , pch=16 )
    lines( b_ests[1,3:4,i] , c(i,i) , lwd=3 , col=2 )
}
abline(v=b_true,lty=2,lwd=3)
