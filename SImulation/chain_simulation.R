#ZR

set.seed(1956)

library(MASS)
library(brms)
library(lme4)
library(lmerTest)
library(ggplot2)

theme_set(theme_grey())


##Data Gen

N2 <- 30 # Number of L2 Groups
N  <- 30  # N per group 

cor_thresh <- c(0.0,.9)


# Makes 100, covariance matricies 
# Random value for cov, by group
index <- 1:N2
sig <- c()

for(i in index){
  x <- runif(n = 1,min = cor_thresh[[1]],max = cor_thresh[[2]]) #correlation boundaries
  sig[[i]] <- matrix(c(1,x, 
                       x,1),2)  
}


# Generate continuos Data
dat1 <- c()
for(i in 1:length(sig)){
  dat1[[i]]<- MASS::mvrnorm(n = N,mu = runif(n = 2,min = -1,max = 1),Sigma = sig[[i]])
}

dat2 <- data.frame(do.call("rbind",dat1))


# Generate Catagorical Data
cat_dat <- c()

for(i in 1:N2){
  cat_dat[[i]] <- rep(i,N)
}

cat_dat2<- unlist(cat_dat)

#put it together
dat2$g <- cat_dat2

dat2$id <- seq(1:nrow(dat2))

colnames(dat2) <- c("y","x","g_id","id")

head(dat2)
cor(dat2[,"x"],dat2[,"y"])

#Plot of data by G_id

#pdf(file = "plotty",height = 50,width = 6)
ggplot(data = dat2) +
  geom_point(aes(x = x, y = y)) +
  geom_smooth(aes(x = x, y = y),method = "lm",se = FALSE)+
  facet_grid(g_id~.)


#ICC
mod0<- lmer(formula = y ~ 1 + (1 | g_id), data = dat2)

summary(mod0)

covrand <- as.data.frame(VarCorr(mod0))
sigmau0 <- covrand[1,4]
sigmaeps <- covrand[2,4]

ICC <- sigmau0/(sigmau0+sigmaeps)

ICC



#Random intercept only model
mod1 <- lmer(formula = y ~ x + (1 | g_id),data = dat2)

summary(mod1)

mod2 <- lmer(formula = y ~ 1 + x + (1 + id | g_id),data = dat2)

summary(mod2)

###################################################3
#Write out data and frequentist model

setwd("/home/z458r456/Documents/Simulation/data/")

# saveRDS(dat2,file = "sim_data.RDS")
# 
# saveRDS(mod1,file = "sim_freq_mod1.RDS")
# 
# saveRDS(mod2,file = "sim_freq_mod2.RDS")
# 
# saveRDS(ICC,file = "sim_ICC.RDS")

dat2 <- readRDS("sim_data.RDS")

# Random intercept and slope

holder <- c()
chainz <- rep(2:10,10)

index <- 1:length(chainz)


form   <- bf(y ~ x + (1 + id | g_id)) 

#good priors
priors <- c(prior(normal(0,1),class = "b", coef = "x"),
            prior(lkj_corr_cholesky(1), class = "L"),
            prior(student_t(3,0,10), class = "sd", coef = "id", group = "g_id"),
            prior(student_t(3,0,10), class = "sigma"))

#diffuse priors 
priors <- c(prior(uniform(-100,100),class = "b", coef = "x"),
            prior(lkj_corr_cholesky(1), class = "L"),
            prior(uniform(0,100), class = "sd", coef = "id", group = "g_id"),
            prior(uniform(0,100), class = "sigma"))

for(i in index){
  
fit <- brm(formula = form, 
            data = dat2, prior = priors,family = gaussian(),
            chains = chainz[[i]],iter = 2000,warmup = 1000)

fit$chainz <- chainz[[i]]

holder[[i]] <- fit 

}
saveRDS(holder,file = "fits_diffuse2_priors.RDS")

###############
#fit 2
library(pander)

fit4$fit

fit3$fit

###############

plot(fit2)

stancode(fit2)
sort(brms::rhat(fit2))

pairs(fit2)

