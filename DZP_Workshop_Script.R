# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3718544/#R26

# DZP Workshop ------------------------------------------------------------

n <- 1 #number of individuals
t <- 1:10   #number of time periods
dat <- expand.grid(t=1:t,id=1:n)

df = data.frame(t = rep(t, n))

df$y = 24 + 10*t - 1.2*t^2 #+ rnorm(nrow(df), 0, 3)
df$y = 24 + 1*t - 1*t^2 + 0.5*t^3#+ rnorm(nrow(df), 0, 3)

df$y = 24 - 4*log(t+1)  #+ rnorm(nrow(df), 0, 3)

ggplot(data = df, aes(x = t, y = y)) +
  geom_point() +
  coord_cartesian(ylim = c(0,50)) +
  theme_classic()

cor(hama_scores)

psych::describe(hama_scores)

data = round(
  faux::rnorm_multi(n = n, 
                    mu = psych::describe(hama_scores)$mean, 
                    sd = psych::describe(hama_scores)$sd, 
                    r = cor(hama_scores), 
                    varnames = colnames(hama_scores), 
                    empirical = F), 2)
