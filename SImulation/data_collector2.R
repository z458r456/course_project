#ZR
#Dec 3 2017
#Chain aggregator



#Collector

library(asbio)
library(reshape2)

setwd("/home/z458r456/Documents/Simulation")

index <- list.files("data/fits/")
folder <- c()

for(i in 1:length(index)) {
  
folder[[i]]<- readRDS(paste0("data/fits/",index[i]))

}

fit_dif   <- folder[[3]][[9]]
fit_str_inc <- folder[[4]][[9]]
fit_str_cor <- folder[[5]][[9]]

rm(folder) 

fit_dif_samples      <- fit_dif$fit@sim$samples
fit_str_inc_samples  <- fit_str_inc$fit@sim$samples
fit_str_cor_samples  <- fit_str_cor$fit@sim$samples

# # # # # # # # # # 
#Population values
# # # # # # # # # # 

# pop_b_Intercept <- 0.00
# pop_b_x <- 0.35
# pop_sd_g_id__Intercept <- 0.95
# pop_sd_g_id__id <- 0.00
# pop_cor_g_id__Intercept__id <- 0.95
# pop_sigma <- 1



pop_vals <- list(0.00,0.35,0.95,0.00,0.95,1)
names(pop_vals) <-   c('b_Intercept','b_x','sd_g_id__Intercept','sd_g_id__id','cor_g_id__Intercept__id',
  'sigma')





# # # # # # # # # # 

#constructing data frames for processing Rhats

index <- c('b_Intercept','b_x','sd_g_id__Intercept','sd_g_id__id','cor_g_id__Intercept__id',
      'sigma','r_g_id[1,Intercept]','r_g_id[2,Intercept]','r_g_id[3,Intercept]',
      'r_g_id[4,Intercept]','r_g_id[5,Intercept]','r_g_id[6,Intercept]',
      'r_g_id[7,Intercept]','r_g_id[8,Intercept]','r_g_id[9,Intercept]',
      'r_g_id[10,Intercept]','r_g_id[11,Intercept]','r_g_id[12,Intercept]',
      'r_g_id[13,Intercept]','r_g_id[14,Intercept]','r_g_id[15,Intercept]',
      'r_g_id[16,Intercept]','r_g_id[17,Intercept]','r_g_id[18,Intercept]',
      'r_g_id[19,Intercept]','r_g_id[20,Intercept]','r_g_id[21,Intercept]',
      'r_g_id[22,Intercept]','r_g_id[23,Intercept]','r_g_id[24,Intercept]',
      'r_g_id[25,Intercept]','r_g_id[26,Intercept]','r_g_id[27,Intercept]',
      'r_g_id[28,Intercept]','r_g_id[29,Intercept]','r_g_id[30,Intercept]',
      'r_g_id[1,id]','r_g_id[2,id]','r_g_id[3,id]','r_g_id[4,id]','r_g_id[5,id]',
      'r_g_id[6,id]','r_g_id[7,id]','r_g_id[8,id]','r_g_id[9,id]','r_g_id[10,id]',
      'r_g_id[11,id]','r_g_id[12,id]','r_g_id[13,id]','r_g_id[14,id]','r_g_id[15,id]',
      'r_g_id[16,id]','r_g_id[17,id]','r_g_id[18,id]','r_g_id[19,id]','r_g_id[20,id]',
      'r_g_id[21,id]','r_g_id[22,id]','r_g_id[23,id]','r_g_id[24,id]','r_g_id[25,id]',
      'r_g_id[26,id]','r_g_id[27,id]','r_g_id[28,id]','r_g_id[29,id]','r_g_id[30,id]',
      'lp__')
if_index<-  c('b_Intercept','b_x','sd_g_id__Intercept','sd_g_id__id','cor_g_id__Intercept__id',
  'sigma')
    

holder <- list()
#Makes list of each parameter estimate by chain.

   
 
# DOING IT 3 TIMES!!! 1 for each condition
dat_sets <-list( fit_dif_samples,    
              fit_str_inc_samples,  
              fit_str_cor_samples) 

final_list <- list()

#THIS DETERMINES THE CONDITION
dat <- dat_sets[[3]]

for(j in index){# j = 'b_Intercept'

  if(j %in% if_index){
    
    pop_value  <- rep(pop_vals[[j]],2000)
    
    holder[[j]]   <- data.frame(dat[[1]][[j]],
                                dat[[2]][[j]],
                                dat[[3]][[j]],
                                dat[[4]][[j]],
                                dat[[5]][[j]],
                                dat[[6]][[j]],
                                dat[[7]][[j]],
                                dat[[8]][[j]],
                                dat[[9]][[j]],
                                dat[[10]][[j]],
                                pop_value)
    
    
  }else{

holder[[j]]   <- data.frame(dat[[1]][[j]],
                            dat[[2]][[j]],
                            dat[[3]][[j]],
                            dat[[4]][[j]],
                            dat[[5]][[j]],
                            dat[[6]][[j]],
                            dat[[7]][[j]],
                            dat[[8]][[j]],
                            dat[[9]][[j]],
                            dat[[10]][[j]])
  }
}

out_list <- list()

itter <- 100

hats<- data.frame(index)
r_bias <- data.frame(if_index)
full_r_bias <- list()
rm(bias)

for(k in 1:itter){#k = 1
  
  for(j in 2:10){# j = 2
    
  cols <- sample(x = 1:10,replace = FALSE,size = j)

    for(i in 1:length(index)){# i = 1
      
    pre_mat<- as.matrix(holder[[i]])
    sample_mat<- pre_mat[,c(cols)]
    
    hats[i,j] <- R.hat(M = sample_mat,burn.in = 0.5)
    
      if(i < 7){ #relative bias
      nam<- if_index[i]
        
      pop <- holder[[i]][,11][1001:2000]  
      bias <- pop - rowSums(sample_mat[1001:2000,])/length(cols)
      r_bias[i,j] <- mean(sqrt((bias)^2))
      
      }

    }

  }
  colnames(r_bias) <- c("parm","2_c","3_c","4_c","5_c","6_c","7_c","8_c","9_c","10_c")  
  full_r_bias[[k]] <- r_bias

colnames(hats) <- c("parm","2_c","3_c","4_c","5_c","6_c","7_c","8_c","9_c","10_c")  
out_list[[k]]  <- hats

}

##rhats
out_dat <- do.call("rbind",out_list)

out_parm <- list()

for(j in index){# j = 'b_Intercept'

id <- 1:length(out_dat[which(out_dat$parm== j ),])
  
out_parm[[j]]  <- data.frame(out_dat[which(out_dat$parm== j ),],id)
 
}

fin_list <- list()
for(j in 1:length(out_parm)){
  
t_dat <- melt(data = out_parm[[j]],id.vars = id,measure.vars = 
            c("X2_c","X3_c","X4_c","X5_c","X6_c","X7_c","X8_c","X9_c","X10_c"),
            variable.name = "chains")

fin_list[[j]]<- t_dat[c("parm","chains","value")]

}

plot_dat_rhat<- do.call("rbind",fin_list)

##Bias

out_dat <- do.call("rbind",full_r_bias)

out_parm <- list()

for(j in if_index){# j = 'b_Intercept'
  
  id <- 1:length(out_dat[which(out_dat$parm== j ),])
  
  out_parm[[j]]  <- data.frame(out_dat[which(out_dat$parm== j ),],id)
  
}

fin_list <- list()
for(j in 1:length(out_parm)){
  
  t_dat <- melt(data = out_parm[[j]],id.vars = id,measure.vars = 
                  c("X2_c","X3_c","X4_c","X5_c","X6_c","X7_c","X8_c","X9_c","X10_c"),
                variable.name = "chains")
  
  fin_list[[j]]<- t_dat[c("parm","chains","value")]
  
}

plot_dat_bias<- do.call("rbind",fin_list)

#Save

setwd("/home/z458r456/Dropbox/Graduate/Courses/Stat_Computing/GIT/course_project/SImulation/")

plot_dat_rhat2 <-  plot_dat_rhat[which(plot_dat_rhat$parm %in% c("sd_g_id__id","sd_g_id__Intercept","sigma","b_Intercept","b_x", "cor_g_id__Intercept__id")),]
plot_dat_bias2 <-  plot_dat_bias[which(plot_dat_bias$parm %in% c("sd_g_id__id","sd_g_id__Intercept","sigma","b_Intercept","b_x", "cor_g_id__Intercept__id")),]

saveRDS(plot_dat_rhat2, file = "rhat_strong_cor.RDS")

saveRDS(plot_dat_bias2, file = "bias_strong_cor.RDS")

############Results Plots
library(ggplot2)



ggplot(data = plot_dat) +
  geom_density(aes(x = value,color = chains,
                   fill = chains), alpha = 0.35) +
  facet_grid(parm ~ .)
  






