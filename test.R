library(ALDEx2)
library(profvis)

data(selex)
# 14 samples, 1600 features

## devtools::load_all('~/Documents/0_git/projects/ALDEx3')
condition <- c(rep(0,7,), rep(1,7))
## X <- formula(~condition) # not needed, can pass formula directly to aldex
data <- data.frame(condition=condition)


profvis({
foo <- aldex(selex, ~condition, data, nsample=1024, scale=default, streamsize=1000)
})
