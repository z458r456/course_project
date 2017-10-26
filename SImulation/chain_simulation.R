#ZR

set.seed(1856)

##Data Gen
library(MASS)

N2 <- 100 # Number of L2 Groups
N  <- 30  # N per group 

# Makes 100, covariance matricies 
# Random value for cov, by group
index <- 1:N2
sig <- c()

for(i in index){
x <- runif(n = 1,min = -0.4,max = 0.90) #correlation boundaries
sig[[i]] <- matrix(c(1,x, 
                 x,1),2)  
}

# Generate continuos Data
dat1 <- c()
for(i in 1:length(sig)){
dat1[[i]]<- MASS::mvrnorm(n = N,mu = c(0,0),Sigma = sig[[i]])
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

dat2$id <- seq(1:length(dat2))

colnames(dat2) <- c("y","x","g_id","id")

head(dat2)
cor(dat2[,1],dat2[,2])



library(brms)

fit6 <- brm(bf(y ~ x, sigma ~ 0 + x), data = dat2)

plot(fit6)



fit2 <- brm(y ~ x + (1+ x|g_id), 
          data = dat2, prior = set_prior("normal(0,1)"), chains = 4,iter = 2000,warmup = 1000)

plot(fit2)

brms::rhat(fit2)[5]



###Model and Priors

fit1 <- brm(formula =
                y ~ x + (1 +|g_id),
            data = dat2, family = lognormal(),
            prior = c(set_prior("normal(0,0.01)", class = "b"),
                      set_prior("cauchy(0,2)", class = "sd"),
                      set_prior("lkj(2)", class = "cor")), warmup = 1000, 
            iter = 5000, chains = 10, control = list(adapt_delta = 0.95))


#####Output

brms::stancode(fit1)
