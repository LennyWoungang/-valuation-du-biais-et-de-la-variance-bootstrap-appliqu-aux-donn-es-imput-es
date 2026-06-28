###STT 6220 projet
#Check how long the simulation takes 
#Consider n = 100 and missing rate is 40%

#install.packages("mice")
library("mice")
#check the document behind mice.impute.pmm() by 
#https://cran.r-project.org/web/packages/mice/mice.pdf

#And the logic behind was taken from flexible imputation of missing data p73
#https://pzs.dstu.dp.ua/DataMining/preprocessing/bibl/fimd.pdf


#Create dataset 
#Attention: our case X is fixed, so we use the same X for the iteration
#but only the epsilon should change ()
#B is the number of bootstrap
#M is the number of simulation


##############################################################################
#Create function

final_project <- function(n, missing_rate, data_list){
  
  B <- 1000
  #M : number of the simulations 
  
  
  s_sq_imp <- numeric()
  var_s_sq_imp_star <- numeric()
  bias_s_sq_imp_star <- numeric()
  
  for(i in 1:M){
    #Create data set
    
    df <- data_list[[i]]
    
    indicator <- runif(n, min = 0, max = 1)
    df$missing_y <- ifelse(indicator < missing_rate, NA, df$y)
    
    #Impute by PMM
    imp <- mice.impute.pmm(y = df$missing_y,
                           ry = !is.na(df$missing_y),
                           x = as.matrix(df[, c("X1","X2", "X3")]))
    df$imputed_y <- ifelse(is.na(df$missing_y), imp, df$missing_y)
    
    
    s_sq_imp[i] <- var(df$imputed_y)
    
    s_sq_imp_star <- numeric()
    
    #Bootstrap
    for(j in 1:B){
      BS_df <- df[sample(1:n, replace = T), ]
      BS_imp <- mice.impute.pmm(y = BS_df$missing_y,
                                ry = !is.na(BS_df$missing_y),
                                x = as.matrix(BS_df[, c("X1","X2", "X3")]))
      BS_df$BS_imputed_y <- ifelse(is.na(BS_df$missing_y), BS_imp, BS_df$missing_y)
      
      #sample variance of the each bootstrap data 
      s_sq_imp_star[j] <- var(BS_df$BS_imputed_y)
      
      
      
    }
    var_s_sq_imp_star[i] <- var(s_sq_imp_star)
    bias_s_sq_imp_star[i] <- mean(s_sq_imp_star) - s_sq_imp[i]
    
  }
  
  
  #M=1000 10 minutes 
  
  
  
  #Get the true variance
  total_data <- do.call(rbind, data_list)
  sigma_sq <- sum((total_data$y - mean(total_data$y))^2)/ (n*M) 
  
  
  #For Variance 
  var_s_sq_imp <- var(s_sq_imp)
  
  var_abs_bias <- abs(mean(var_s_sq_imp_star) - var_s_sq_imp)
  var_relative_bias <- var_abs_bias / var_s_sq_imp
  
  #For Bias
  bias_s_sq_imp <- mean(s_sq_imp) - sigma_sq
  bias_abs_bias <- abs(mean(bias_s_sq_imp_star) - bias_s_sq_imp)
  bias_relative_bias <- bias_abs_bias / bias_s_sq_imp
  
  return(
    list(
      sigma_sq = sigma_sq,
      var_s_sq_imp = var_s_sq_imp,
      mean_of_var_s_sq_imp_star = mean(var_s_sq_imp_star),
      var_abs_bias = var_abs_bias,
      var_relative_bias = var_relative_bias,
      bias_s_sq_imp = bias_s_sq_imp,
      mean_of_bias_s_sq_imp_star = mean(bias_s_sq_imp_star),
      bias_abs_bias = bias_abs_bias,
      bias_relative_bias = bias_relative_bias
    )
  )
}



#################################################################################
#Simulation 

set.seed(221)
simulation_data40 <- list()
simulation_data100 <- list()


M <- 1000
n <- 40 
X1 <- rnorm(n)
X2 <- rnorm(n)
X3 <- rnorm(n)

b1 <- 1
b2 <- 13
b3 <- -2

#Create dataset n = 40
for(i in 1:M){
   epsilon <- rnorm(n, mean = 0, sd = 4)
  y <- b1*X1 + b2*X2 + b3*X3 + epsilon
  df <- data.frame(X1=X1, X2=X2, X3=X3, epsilon, y=y)
  
  simulation_data40[[i]] <- df
}

#Create datasets n = 100
n <- 100
X1 <- rnorm(n)
X2 <- rnorm(n)
X3 <- rnorm(n)

for(i in 1:M){
  epsilon <- rnorm(n, mean = 0, sd = 4)
  y <- b1*X1 + b2*X2 + b3*X3 + epsilon
  df <- data.frame(X1=X1, X2=X2, X3=X3, epsilon, y=y)
  
  simulation_data100[[i]] <- df
}


n40p10 <- final_project(40, 0.1, simulation_data40)
print(1)
n40p20 <- final_project(40, 0.2, simulation_data40)
print(2)
n40p30 <- final_project(40, 0.3, simulation_data40)
print(3)
n100p10 <- final_project(100, 0.1, simulation_data100)
print(4)
n100p20 <- final_project(100, 0.2, simulation_data100)
print(5)
n100p30 <- final_project(100, 0.3, simulation_data100)
print(6)

