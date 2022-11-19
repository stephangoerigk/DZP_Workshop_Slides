# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3718544/#R26

# DZP Workshop ------------------------------------------------------------

library(lme4)
library(lmerTest)
library(dendextend)


# Simulate Practice data ---------------------------------------------------

n = 350 # number of individuals
t = 1:10   # number of time periods

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

hama_cor = read.csv("https://raw.githubusercontent.com/stephangoerigk/DZP_Workshop_Slides/master/hamacor.csv")
hama_cor = as.matrix(hama_cor[,-1])

data_cluster = round(
  faux::rnorm_multi(n = n,
                    mu = hama_means,
                    sd = hama_sds,
                    r = hama_cor,
                    varnames = hama_names,
                    empirical = F), 2)

psych::describe(data_cluster)

data_cluster = scale(data_cluster)

data_transposed = t(na.omit(data_cluster))

d = dist(data_transposed, method = "euclidean")
clust = hclust(d, method = "ward.D2")

dend <- as.dendrogram(clust, hang = -1)
pruned = dynamicTreeCut::cutreeDynamic(clust, distM = as.matrix(d), method = "hybrid", minClusterSize = 1)

names(pruned) = hama_names

labels_colors(dend) = pruned[c(clust$order)]
labels_cex(dend) = 2
marg = c(4, 4, 10, 35)
par(mar = marg, font = 1, cex = 0.4, cex.axis = 1.7, cex.lab = 2)
plot(rev(dend), horiz = T, edgePar = list(lwd = 2))
p1 <- recordPlot()

df = expand.grid(t = 1:max(t),
                 id = 1:n)
df$group = c(rep("active", nrow(df)/2), rep("placebo", nrow(df)/2))

trajectory = c("Linear response",
               "Deteriorate",
               "Rev. U-shape",
               "Rapid response",
               "No change")

set.seed(123)
for(ch in unique(df$id)){

  if(df$group[df$id == ch][1] == "active"){
    df$trajectory[df$id == ch] = rep(sample(trajectory, size = 1, replace = T, prob = c(.5, .05, .2, .2, .05)), max(t))
  }
  if(df$group[df$id == ch][1] == "placebo"){
    df$trajectory[df$id == ch] = rep(sample(trajectory, size = 1, replace = T, prob = c(.2, .2, .1, .05, .45)), max(t))
  }

  if(df$trajectory[df$id == ch][1] == "No change"){
    df$y[df$id == ch] = 24 + 0*t  + rnorm(nrow(df[df$id == ch,]), 0, 3)
  }
  if(df$trajectory[df$id == ch][1] == "Rev. U-shape"){
    df$y[df$id == ch] = 24 + 8*t - 0.9*t^2 + rnorm(nrow(df[df$id == ch,]), 0, 3)
  }
  if(df$trajectory[df$id == ch][1] == "Linear response"){
    df$y[df$id == ch] = 24 - 1*t  + rnorm(nrow(df[df$id == ch,]), 0, 3)
  }
  if(df$trajectory[df$id == ch][1] == "Deteriorate"){
    df$y[df$id == ch] = 24 + 2*t  + rnorm(nrow(df[df$id == ch,]), 0, 3)
  }
  if(df$trajectory[df$id == ch][1] == "Rapid response"){
    df$y[df$id == ch] = 24 - 10 * log(t) +  rnorm(nrow(df[df$id == ch,]), 0, 3)
  }
}

ggplot(data = df, aes(x = t, y = y, colour = group)) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean") +
  scale_x_continuous(breaks = t) +
  coord_cartesian(ylim = c(0,50)) +
  labs(x = "Time", colour = "Group") +
  theme_classic()

summary(lmer(y ~ t * group + (1|id), data = df))

ggplot(data = df, aes(x = t, y = y, colour = trajectory)) +
  geom_point() +
  scale_x_continuous(breaks = t) +
  coord_cartesian(ylim = c(0,50)) +
  facet_grid(.~trajectory) +
  labs(colour = "Trajectory", x = "Time") +
  theme_classic()

ggplot(data = df, aes(x = t, y = y, colour = trajectory)) +
  geom_point() +
  scale_x_continuous(breaks = t) +
  coord_cartesian(ylim = c(0,50)) +
  facet_grid(cols = vars(trajectory), rows = vars(group)) +
  labs(colour = "Trajectory", x = "Time") +
  theme_classic()

set.seed(222)
results_total = data.frame(ng = NA,
                           Polynomial = NA,
                           Random = NA,
                           BIC = NA,
                           AIC = NA,
                           loglik = NA)

for(ng in 2:5){
  for(random in c("~ 1", "~ 1 + t")){

    mi_sq <- lcmm::hlme(fixed = y ~ 1 + t + I(t^2),
                        mixture = ~ 1 + t + I(t^2),
                        random = as.formula(random),
                        ng = ng,
                        nwg = FALSE,
                        idiag = FALSE,
                        data = df,
                        subject = "id")
    mi_cub <- lcmm::hlme(fixed = y ~ 1 + t + I(t^3),
                         mixture = ~ 1 + t + I(t^3),
                         random = as.formula(random),
                         ng = ng,
                         nwg = FALSE,
                         idiag = FALSE,
                         data = df,
                         subject = "id")

    lin <- c(mi_lin$ng, 1, random, mi_lin$BIC, mi_lin$AIC, mi_lin$loglik)
    sq <- c(mi_sq$ng, 2, random, mi_sq$BIC, mi_sq$AIC, mi_sq$loglik)
    cub <- c(mi_cub$ng, 3, random, mi_cub$BIC, mi_cub$AIC, mi_cub$loglik)
    results_total = rbind(results_total, lin)
    results_total = rbind(results_total, sq)
    results_total = rbind(results_total, cub)
  }
}

results_total = results_total[order(results_total$BIC, decreasing = T),]

mi = lcmm::hlme(fixed = y ~ 1 + t + I(t^2),
           mixture = ~ 1 + t + I(t^2),
           random = ~ 1,
           ng = 5,
           nwg = FALSE,
           idiag = FALSE,
           data = df,
           subject = "id")

LCTMtoolkit_total = LCTMtoolkit(mi)
postprob_total = lcmm::postprob(mi)

transfer_class = function(data, mi){
  data$class = NA
  for(ch in unique(data$id)){
    data$class[data$id == ch] = mi$pprob$class[mi$pprob$id == ch]
  }
  return(data)
}

df = transfer_class(data = df, mi = mi)

plot_traj = function(mi, data, var.time){
  datnew   <- data.frame(t = seq(0, max(data[, var.time]), length = 100))
  plotpred <- lcmm::predictY(mi, datnew, var.time = var.time, draws = TRUE)

  frame_traj = as.data.frame(expand.grid(Time = plotpred$times$t,
                                         trajectory = unique(mi$pprob$class),
                                         pred = NA,
                                         upper = NA,
                                         lower = NA))

  for(traj in unique(frame_traj$trajectory)){
    for(i in 1:100){
      frame_traj$pred[frame_traj$trajectory == traj][i] = plotpred$pred[,which(grepl(paste0("^Ypred_class", as.character(traj)), colnames(plotpred$pred)))][i]
      frame_traj$upper[frame_traj$trajectory == traj][i] = plotpred$pred[,which(grepl(paste0("^lower.Ypred_class", as.character(traj)), colnames(plotpred$pred)))][i]
      frame_traj$lower[frame_traj$trajectory == traj][i] = plotpred$pred[,which(grepl(paste0("^upper.Ypred_class", as.character(traj)), colnames(plotpred$pred)))][i]
    }
  }
  frame_traj$trajectory = factor(frame_traj$trajectory)
  ggplot(data = frame_traj, aes(x = Time, y = pred, ymin = lower, ymax = upper)) +
    # geom_vline(xintercept = c(5, 10, 14, 17), linetype = "dotted") +
    geom_line(aes(colour = trajectory)) +
    geom_ribbon(aes(fill = trajectory), alpha = .2, linetype = "dotted") +
    theme_classic()
}



# Dynamic Time Warping ----------------------------------------------------


nodep <- 5
dep <- 20

# no depressive episode
ts_sim_nodep <- abs(arima.sim(n = 420, mean = 0.001, model = list(order = c(1,0,0), ar = 0.9))) + nodep

ts.plot(ts_sim_nodep, xlab = "Time", ylab = "Negative Mood", main = "Example no depressive episode", ylim=c(0,25))

# depressive episode
ts_sim_part1 <- abs(arima.sim(n = 210, model = list(order = c(1,0,0), ar = 0.9))) + nodep
ts_sim_part2 <- ts(arima.sim(n = 210, model = list(order = c(1,0,0), ar = 0.9)) + dep, start = 211,end =420)

ts_sim_dep <- ts.union(ts_sim_part1,ts_sim_part2)
ts_sim_dep<- pmin(ts_sim_dep[,1], ts_sim_dep[,2], na.rm = TRUE)

ts.plot(ts_sim_dep, xlab = "Time", ylab = "Negative Mood", main = "Example depressive episode", ylim=c(0,25))

ts_sim_boot_classic <- ts_sim_nodep %>%
  tseries::tsbootstrap(., nb=5, b=200, type = "block") %>%
  as.data.frame(.) %>%
  dplyr::rename_all(funs(c(paste0("nodep_",.))))

ts_sim_boot_dep <- ts_sim_dep %>%
  tseries::tsbootstrap(., nb=5, b=350, type = "block") %>%
  as.data.frame(.) %>%
  dplyr::rename_all(funs(c(paste0("dep_",.))))

ts_sim_df <- cbind(ts_sim_boot_classic,ts_sim_boot_dep)

dtw_dist <- function(x){dist(x, method="DTW")}

ts_sim_df %>%
  as.matrix() %>%
  gplots::heatmap.2 (
    # dendrogram control
    distfun = dtw_dist,
    hclustfun = hclust,
    dendrogram = "column",
    Rowv = FALSE,
    labRow = FALSE
  )

cluster_dtw_h2 <- dtwclust::tsclust(t(ts_sim_df),
                                    type = "h",
                                    k = 2,
                                    distance = "dtw",
                                    control = hierarchical_control(method = "complete"),
                                    preproc = NULL,
                                    args = tsclust_args(dist = list(window.size = 5L)))

hclus <- stats::cutree(cluster_dtw_h2, k = 2) %>%
  as.data.frame(.) %>%
  dplyr::rename(.,cluster_group = .) %>%
  tibble::rownames_to_column("type_col")

hcdata <- ggdendro::dendro_data(cluster_dtw_h2)
names_order <- hcdata$labels$label

p1 <- hcdata %>%
  ggdendro::ggdendrogram(., rotate=TRUE, leaf_labels=FALSE)

p2 <- ts_sim_df %>%
  dplyr::mutate(index = 1:420) %>%
  tidyr::gather(key = type_col,value = value, -index) %>%
  dplyr::full_join(., hclus, by = "type_col") %>%
  mutate(type_col = factor(type_col, levels = rev(as.character(names_order)))) %>%
  ggplot(aes(x = index, y = value, colour = cluster_group)) +
  geom_line() +
  facet_wrap(~type_col, ncol = 1, strip.position="left") +
  guides(color=FALSE) +
  theme_bw() +
  theme(strip.background = element_blank(), strip.text = element_blank())


gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2+ labs(y = "Negative Mood", x = "Time"))

gridExtra::grid.arrange(gp2, gp1, ncol=2, widths=c(4,2))

