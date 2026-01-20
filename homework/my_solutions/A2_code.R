# Statistical Rethinking 2026
## A2 Solution (Mukto Akash)
library(rethinking)

set.seed(89)

# define grid
p_grid <- seq( from=0 , to=1 , length.out=1000 )
# define prior Uniform(0, 1)
prior <- rep( 1 , 1000 )
# compute likelihood at each value in grid
likelihood <- dbinom( 3 , size=14 , prob=p_grid )
# compute product of likelihood and prior
unstd.posterior <- likelihood * prior
# standardize the posterior, so it sums to 1
posterior <- unstd.posterior / sum(unstd.posterior)

plot( p_grid , posterior , type="b" ,
      xlab="probability of water" , ylab="posterior probability" )
mtext( "Posterior Distribution" )

probabilities <- sample( p_grid , prob=posterior , size=1e4 , replace=TRUE )

w <- rbinom( 1e4 , size=5 , prob=probabilities )
simplehist( w , xlab="water count from 5 globe tosses" )
