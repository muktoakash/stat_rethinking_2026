possib <- c("H", "T")
count <- 0

x_heads_from_n_flips <- function(num_flips, num_heads){
  if (num_heads < 0 || num_heads > num_flips) {
    return(0)
  } else if (num_heads == num_flips){
    return(1)
  }
  
  for (flip in possib) {
    if (flip == "H"){
      count <- count + x_heads_from_n_flips(num_flips - 1, num_heads - 1)
    } else {
      count <- count + x_heads_from_n_flips(num_flips -1, num_heads)
    }
  }
  
  return(count)
}

print(x_heads_from_n_flips(10, 0))
print(count)
print(x_heads_from_n_flips(10, 10))
print(count)
print(x_heads_from_n_flips(3, 2))
print(count)
print(x_heads_from_n_flips(4, 2))
print(count)
print(x_heads_from_n_flips(10, 8))