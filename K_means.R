require(data.table)



# Simulate group data
# Functions:
# (1) generate data
# (2) Lloyds algorithm

# (1)
gen_data <-function(dgp,prob_G){
  Ndim = dgp@Ndim
  Tdim = dgp@Tdim
  G = dgp@G
  siga = dgp@sig_alpha
  sigerr = dgp@sig_err
  
  alpha <- sort(rnorm(G,0,siga))
  error <- matrix(rnorm(Ndim*Tdim,0,sigerr),Ndim,Tdim)
  
  assign <- sample(sample(1:G,Ndim,replace=TRUE,prob=prob_G))
  
  alpha_mat <- matrix(rep(alpha[assign],Tdim),Ndim,Tdim)
  y = alpha_mat + error
  return(list(y,alpha,assign))
}

# (2) Lloyd's algorithm
Lloyds <- function(dgp,y,alpha_ini,itermax){
  Ndim <- dgp@Ndim
  Tdim <- dgp@Tdim
  G <- dgp@G
  tol <- 0.0001
  dif <- 1
  data <- CJ(i=1:Ndim,t=1:Tdim)
  data[,y:=as.vector(t(y))]
  
  res = matrix(0,Ndim,G)
  alpha_new = rep(0,G)
  iter = 1
  while(dif > tol &&  iter < 50){
    # assignment
    for(gg in 1:G){
      resaux = y -  matrix(rep(alpha_ini[gg],Ndim*Tdim),Ndim,Tdim)
      resaux2 = resaux^2
      res[,gg] = rowSums(resaux2)
    }
    new_assign = apply(res, 1, FUN = which.min)
    data[,assign:=as.character(new_assign[i])]
    # update alpha
    model <- glm(formula = y ~ assign + 0, data = data)
    dif = 0
    for(gg in 1:G){
      alpha_new[gg] = model$coefficients[[gg]]
      dif = dif +  (alpha_new[gg] - alpha_ini[gg])^2
    }
    #print(dif)
    alpha_ini = alpha_new
    iter = iter+1
    #print(iter)
  }
  data[,alpha_out:=alpha_ini[as.numeric(assign)]]
  return(data)
}


# MAIN

setClass("DGP", slots = list(Ndim = "numeric",Tdim = "numeric",G = "numeric",sig_alpha = "numeric", sig_err = "numeric"))


dgp <- new("DGP",Ndim=1000,Tdim=10,G=3,sig_alpha = 10,sig_err=1)
prob_G = c(0.3,0.3,0.4)
itermax = 50

# Generate data:
data_out <- gen_data(dgp,prob_G)
y <- data_out[[1]]
alpha0 <- data_out[[2]]
assign0 <- data_out[[3]]

# Estimate the model
ini_cond = 10
for(ii in 1:ini_cond){
  print(ii)
  alpha_ini = quantile(y,seq(1/(dgp@G+1),1-1/(dgp@G+1),1/(dgp@G+1))) + rnorm(dgp@G)
  dataoutput <- Lloyds(dgp,y,alpha_ini)
  loss = sum((dataoutput$y-dataoutput$alpha_out)^2)
  
  missclass <- 1-sum(dataoutput$assign[seq(1, dgp@Ndim*dgp@Tdim, by = dgp@Tdim)] == assign0)/dgp@Ndim
  print(list(missclass,loss))
  
  if(ii ==1){
    min_loss = loss
    best_data = dataout
    missclass <- 1-sum(best_data$assign[seq(1, dgp@Ndim*dgp@Tdim, by = dgp@Tdim)] == assign0)/dgp@Ndim
    print(list(missclass,loss))
    best = ii
  }else{ 
    if(loss < min_loss){
      min_loss = loss
      best_data = data_out
      # Compute missclassification probability:
      best = ii
    }
  }
}

missclass <- 1-sum(best_data$assign[seq(1, dgp@Ndim*dgp@Tdim, by = dgp@Tdim)] == assign0)/dgp@Ndim
print(list(missclass,loss))



