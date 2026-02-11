library(rethinking)
data(Howell1)
d <- Howell1[Howell1$age >= 18,]

blank2(w=1.6)

h <- round( d$height/2 )*2
hm <- h[d$male==1]
hf <- h[d$male==0]

pth <- function(y,col=2,lwd=4,ox=0.25) {
    thf <- table(y)
    for ( i in 1:length(thf) ) {
        x <- as.numeric(names(thf)[i])
        lines(c(x,x)+ox,c(0,thf[i]),lwd=lwd,col=col)
    }
}

plot( table(h) , lwd=3 , col="white" , xlab="height (cm)" , ylab="frequency" )
pth(hm,col=4,ox=-0.25)
pth(hf,col=2)

w <- round( d$weight/2 )*2
wm <- w[d$male==1]
wf <- w[d$male==0]

plot( table(w) , lwd=3 , col="white" , xlab="weight (kg)" , ylab="frequency" )
pth(wm,col=4,ox=-0.25)
pth(wf,col=2)

