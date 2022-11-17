# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3718544/#R26

# DZP Workshop ------------------------------------------------------------


# Simulate Practice data ---------------------------------------------------

n <- 400 #number of individuals
t <- 1:10   #number of time periods

hama_names = c("Anxious Mood",
               "Tension",
               "Fears",
               "Insomnia",
               "Concentration and Memory",
               "Depressed Mood",
               "General somatic symptoms: muscular",
               "General somatic symptoms: sensory",
               "Cardiovascular symptoms",
               "Respiratory symptoms",
               "Gastro-intestinal symptoms",
               "Genito-urinary symptoms",
               "Other autonomic symptoms")

hama_means = c(2.9669421,
               2.8044077,
               2.4559229,
               2.4297521,
               1.3815427,
               1.3071625,
               1.6129477,
               1.5633609,
               1.4531680,
               1.0330579,
               1.7190083,
               0.6694215,
               1.6198347)

hama_sds = c(0.7136179,
             0.8339568,
             1.1314254,
             1.2283532,
             1.2313786,
             1.2454069,
             1.2715661,
             1.2264810,
             1.2583008,
             1.2200659,
             1.2815280,
             1.0691029,
             1.1585213)



dat <- expand.grid(t=1:t,id=1:n)

df = data.frame(t = rep(t, n))

df$y = 24 + 10*t - 1.2*t^2 #+ rnorm(nrow(df), 0, 3)
df$y = 24 + 1*t - 1*t^2 + 0.5*t^3#+ rnorm(nrow(df), 0, 3)

df$y = 24 - 4*log(t+1)  #+ rnorm(nrow(df), 0, 3)

ggplot(data = df, aes(x = t, y = y)) +
  geom_point() +
  coord_cartesian(ylim = c(0,50)) +
  theme_classic()


psych::describe(hama_scores)
cor(hama_scores)

data = round(
  faux::rnorm_multi(n = n,
                    mu = psych::describe(hama_scores)$mean,
                    sd = psych::describe(hama_scores)$sd,
                    r = cor(hama_scores),
                    varnames = colnames(hama_scores),
                    empirical = F), 2)

hama_imp = scale(data)

hamatransposed = t(na.omit(hama_imp))

d = dist(hamatransposed, method = "euclidean")
clust = hclust(d, method = "ward.D2")

dend <- as.dendrogram(clust, hang = -1)
pruned = dynamicTreeCut::cutreeDynamic(clust, distM = as.matrix(d), method = "hybrid", minClusterSize = 1)

names(pruned) = colnames(hama_scores)

marg = c(4, 4, 10, 35)
marg2 = c(4, 1, 10, 40)

labels_colors(dend) = pruned[c(clust$order)]
labels_cex(dend) = 2
par(mar = marg, font = 1, cex = 0.4, cex.axis = 1.7, cex.lab = 2)
plot(rev(dend), horiz = T, main = "Dendrogram", edgePar = list(lwd = 2))
axis(1, font=2, cex = 5)
p1 <- recordPlot()
