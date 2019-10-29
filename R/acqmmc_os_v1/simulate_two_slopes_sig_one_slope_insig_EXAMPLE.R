 ##
##
## Simulate 
##
##
library(ggplot2)
library(xtable)
library(lmtest)

yf <- function(a, b, x, mu, sig) {
  return(a + b*x + rnorm(length(x), mu, sig))
}

printPval <- function(fit, param=2, digits=8) {
  eff <- summary(fit)
  p <- eff$coefficients[param,'Pr(>|t|)']
  if (p < 0.001) return('(p<0.001)***')
  if (p < 0.01) return(sprintf('(p=%.3f)**',p))
  if (p < 0.05) return(sprintf('(p=%.3f)*',p))
  return(sprintf('(p=%.3f)',p))
}

## RANDOM SAMPLE
for (. in 1)
{
  # seed <- sample(1:9e5, 1, 1)
  seed <- 549897
  set.seed(seed)
  n <- 30 
  #
  mu1  <- 0
  sig1 <- 13
  mu2  <- 0
  sig2 <- 7
  #
  a1 <- 55
  b1 <- -12
  x1 <- runif(n, min=0.5, max=5.5)
  y1 <- yf(a1, b1, x1, mu1, sig1)
  df1 <- data.frame(x=x1, y=y1, Group='A')
  #
  a2 <- 39
  b2 <- -1.5
  x2 <- runif(n, min=0.5, max=5.5)
  y2 <- yf(a2, b2, x2, mu2, sig2)
  df2 <- data.frame(x=x2, y=y2, Group='B')
  #
  idx.sample <- sample(seq_len(n), size = n/2, replace = F)
  #
  df1 <- df1[idx.sample, ]
  df2 <- df2[idx.sample, ]
  df12 <- rbind(df1, df2)
  #
  m1 <- lm(y ~ x, data = df1)
  m2 <- lm(y ~ x, data = df2)
  m3 <- lm(y ~ x, data = df12)
  
  dfl <- rbind(within(df12,{Functions='Separate'}), within(df12,{Group='A or B'; Functions='Combined'}))
}

# PLOTTING
ann_text <- data.frame(x = c(2.5,3.8,4),
                       y = c(15,58,55),
                       lab = c(sprintf('b1 = %.1f %s',m1$coefficients[2],printPval(m1)),
                               sprintf('b2 = %.1f %s',m2$coefficients[2],printPval(m2)),
                               sprintf('b = %.1f %s', m3$coefficients[2],printPval(m3))),
                       Group = factor(c('A','B','A or B'), levels = c('A', 'B', 'A or B')),
                       Functions = factor(c('Separate','Separate','Combined'), levels = c('Combined', 'Separate')))
ggplot(data=dfl, aes(x=x, y=y, colour=Group, shape=Group, text=NULL)) +
  geom_point(size=3) +  geom_smooth( method='lm',se=F) +
  facet_wrap(.~Functions)  + theme_bw() + theme(legend.position='top') +
  geom_text(data = ann_text, label=ann_text$lab, lwd=4, show.legend = F) +
  scale_color_manual(values = c('darkblue','red','black'))
ggsave('simulate_separate_slopes_signif_combined_insig_EXAMPLE.png', width = 6.5, height = 3.8, units='in', dpi = 200)

## Table
screenreg(list(m1, m2, Combined=m3))
cat(sprintf('\n seed = %s\n', seed))


## NO Heteroskedasticity in full model 
bptest(m3)



# 
# 
# ggplot(data=df12, aes(x=x, y=y, colour=Group, shape=Group)) + 
#   geom_point(size=3) +  geom_smooth(method='lm',se=F) +
#   theme_bw() + ggtitle('Separate Group Trends')
# 
# 
# df12<- data.frame(x=x12, y=y12)
# 
# dfl <- rbind(within(),within())
# 
# ggplot(data=df1, aes(x=x1, y=y1)) + 
#   geom_point(col='black') +  geom_smooth(method='lm', col='black') + theme_bw()
# ggplot(data=df2, aes(x=x2, y=y2)) + 
#   geom_point(col='red') +  geom_smooth(method='lm', col='red') + theme_bw()
# ggplot(data=df3, aes(x=x2, y=y2)) + 
#   geom_point(col='red') +  geom_smooth(method='lm', col='red') + theme_bw()
# 

