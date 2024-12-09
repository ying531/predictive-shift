## =========================================== ##
###  compute KL divergence across two sites   
# generalize from df1 to df2 to be consistent with other estimators
# compute KL divergence E_P[Q/P * log(Q/P)] 
# where P is the weighted df1 distribution
## =========================================== ##

compute_KL_weighted <- function(df1, df2, alg = 'logistic',
                                covariates, 
                                analysis_formula, 
                                treatment_variable
){
  lm.1 = lm(analysis_formula, data = df1)
  lm.2 = lm(analysis_formula, data = df2)
  X.1 = model.matrix(formula(analysis_formula), data=df1)
  y.1 = model.response(model.frame(formula(analysis_formula), data=df1))
  X.2 = model.matrix(formula(analysis_formula), data=df2)
  y.2 = model.response(model.frame(formula(analysis_formula), data=df2))
  df1$infl = t(solve(t(X.1) %*% X.1 / nrow(X.1)) %*%  t(diag(y.1) %*% X.1))[,treatment_variable]
  df2$infl = t(solve(t(X.2) %*% X.2 / nrow(X.2)) %*%  t(diag(y.2) %*% X.2))[,treatment_variable]
  
  # compute KL divergence between joint distributions of (X, influence function)
  
  df1_KL = df1[,c(covariates, "infl")]
  df2_KL = df2[,c(covariates, "infl")] 
  
  df_KL = rbind(df1_KL, df2_KL)
  df_KL$y = c(rep(0, nrow(df1_KL)), rep(1, nrow(df2_KL)))
  
  if (alg == 'logistic'){
    mdl = glm(y~., family='binomial', data = df_KL)
    drs = predict(mdl, type = 'response')
    drs = drs / (1-drs) * nrow(df1_KL) / nrow(df2_KL)
    
    cov.mdl = glm(y~., family = 'binomial', data=df_KL[colnames(df_KL)!='infl'])
    cov.drs = predict(cov.mdl, type= 'response')
    cov.drs = cov.drs / (1-cov.drs) * nrow(df1_KL) / nrow(df2_KL)
  }
  
  if (alg == 'grf'){
    mdl = regression_forest(rbind(df1_KL, df2_KL),
                            c(rep(0, nrow(df1_KL)), rep(1, nrow(df2_KL))))
    drs = predict(mdl)$predictions
    drs = drs / (1-drs) * nrow(df1_KL) / nrow(df2_KL)
    
    cov.mdl = regression_forest(rbind(df1[,covariates], df2[,covariates]),
                                c(rep(0, nrow(df1_KL)), rep(1, nrow(df2_KL))))
    cov.drs = predict(cov.mdl)$predictions
    cov.drs = cov.drs / (1-cov.drs) * nrow(df1_KL) / nrow(df2_KL)
  }
  
  drs1 = drs[1:nrow(df1_KL)]
  drs1 = drs1 / mean(drs1)
  cov.drs1 = cov.drs[1:nrow(df1_KL)]
  cov.drs1 = cov.drs1 / mean(cov.drs1)
  
  KL_est = as.numeric(mean(drs1 * log(drs1 / cov.drs1)))
  
  return(KL_est)
  
}
 
## =========================================== ##
# utilitiy functions for optimizing KL-bounded DRO

grad.kl.w <- function(ifl, lam, delta, w){
  return( - log(mean( w*exp(-ifl/lam ))) - delta - mean(w * ifl * exp(-ifl / lam)) / lam / mean(w *exp(-ifl/lam)) )
}

obj.kl.w <- function(ifl, lam, delta, w){
  return( - lam  * ( log(mean(w*exp(- ifl/lam ))) + delta) )
}

KL_single.w <- function(ifl, delta, w){ 
  lam.left = 0.00001
  while (obj.kl.w(ifl, lam.left, delta, w) %in% c(Inf, - Inf)){
    lam.left = lam.left * 2
  }
  lam.right = 10
  while (grad.kl.w(ifl, lam.right, delta, w) >0){
    lam.right = lam.right * 2
  }
  lam.mid = (lam.left + lam.right) / 2
  grad.left = grad.kl.w(ifl, lam.left, delta, w)
  grad.right = grad.kl.w(ifl, lam.right, delta, w)
  grad.mid = grad.kl.w(ifl, lam.mid, delta, w)
  
  while (abs(grad.mid) > 0.0001){
    if (grad.left * grad.mid > 0){
      lam.left = lam.mid 
      lam.mid = (lam.left + lam.right) / 2
      grad.left = grad.mid 
      grad.mid = grad.kl.w(ifl, lam.mid, delta, w)
    }else{
      lam.right = lam.mid
      lam.mid = (lam.left + lam.right) / 2
      grad.right = grad.mid 
      grad.mid = grad.kl.w(ifl, lam.mid, delta, w)
    }
  }
  
  lam.final = lam.mid
  obj.final = obj.kl.w(ifl, lam.final, delta, w)
  
  return(list("opt.lam" = lam.final, 
              "opt.obj" = obj.final))
}


## ================================================================ ##
###  compute worst-case KL bounds under constraint  
#    E_P[Q/P * log(Q/P)]  <= delta
#    where P is the weighted df1 distribution
#    generalize from df1 to df2 to be consistent with other estimators 
## ================================================================ ##

compute_KL_bound_weighted <- function(df1, df2, delta, alg,
                             covariates, 
                             analysis_formula, 
                             treatment_variable){
  lm.1 = lm(analysis_formula, data = df1)
  lm.2 = lm(analysis_formula, data = df2)
  X.1 = model.matrix(formula(analysis_formula), data=df1)
  y.1 = model.response(model.frame(formula(analysis_formula), data=df1))
  X.2 = model.matrix(formula(analysis_formula), data=df2)
  y.2 = model.response(model.frame(formula(analysis_formula), data=df2))
  df1$infl = t(solve(t(X.1) %*% X.1 / nrow(X.1)) %*%  t(diag(y.1) %*% X.1))[,treatment_variable]
  df2$infl = t(solve(t(X.2) %*% X.2 / nrow(X.2)) %*%  t(diag(y.2) %*% X.2))[,treatment_variable]
  
  df1_KL = df1[,c(covariates, "infl")]
  df2_KL = df2[,c(covariates, "infl")] 
  
  df_KL = rbind(df1_KL, df2_KL)
  df_KL$y = c(rep(0, nrow(df1_KL)), rep(1, nrow(df2_KL)))
  
  if (alg == 'logistic'){ 
    cov.mdl = glm(y~., family = 'binomial', data=df_KL[colnames(df_KL)!='infl'])
    cov.drs = predict(cov.mdl, type= 'response')
    cov.drs = cov.drs / (1-cov.drs) * nrow(df1_KL) / nrow(df2_KL)
  }
  
  if (alg == 'grf'){ 
    
    cov.mdl = regression_forest(rbind(df1[,covariates], df2[,covariates]),
                                c(rep(0, nrow(df1_KL)), rep(1, nrow(df2_KL))))
    cov.drs = predict(cov.mdl)$predictions
    cov.drs = cov.drs / (1-cov.drs) * nrow(df1_KL) / nrow(df2_KL)
  }
  
  cov.drs1 = cov.drs[1:nrow(df1_KL)]
  cov.drs1 = cov.drs1 / mean(cov.drs1)
  
  ifl = df1$infl
  kl.res = NULL
  
  kl.res <- withTimeout({
    KL_single.w(ifl, delta, cov.drs1) 
  }, timeout = 20, onTimeout = "warning") 
  
  if (is.null(kl.res)){
    est.lower = min(ifl)
  }else{
    est.lower = kl.res$opt.obj
  }
  
  kl.res.upp <- withTimeout({
    KL_single.w(-ifl, delta, cov.drs1) 
  }, timeout = 20, onTimeout = "warning") 
  
  if (is.null(kl.res.upp)){
    est.upper = max(ifl)
  }else{
    est.upper = - kl.res.upp$opt.obj
  } 
  
  return(list("lower" = est.lower, 
              "upper" = est.upper,
              "theta1" = mean(df1$infl),
              "theta2" = mean(df2$infl)
  ))
  
}
 
###################################################
# generalize from df1 (rep study in decomp analysis function) to df2 (org study in decomp analysis function)
# to be consistent with standard deviation calculation and doubly robust estimator
# use stable diff-in-mean estimator for Delta_X / s_X
###################################################

compute_basic <- function(df1, df2, K=5,  
                          covariates = NULL,
                          analysis_formula = NULL, 
                          treatment_variable = NULL){ 
  
  lm.1 = lm(analysis_formula, data = df1)
  lm.2 = lm(analysis_formula, data = df2)
  X.1 = model.matrix(formula(analysis_formula), data=df1)
  y.1 = model.response(model.frame(formula(analysis_formula), data=df1))
  X.2 = model.matrix(formula(analysis_formula), data=df2)
  y.2 = model.response(model.frame(formula(analysis_formula), data=df2))
  df1$infl = t(solve(t(X.1) %*% X.1 / nrow(X.1)) %*%  t(diag(y.1) %*% X.1))[,treatment_variable]
  df2$infl = t(solve(t(X.2) %*% X.2 / nrow(X.2)) %*%  t(diag(y.2) %*% X.2))[,treatment_variable]
  df1$infl_centered = df1$infl - mean(df1$infl)
  df2$infl_centered = df2$infl - mean(df2$infl)
  
  # stablize conditional shift 
  X1 = df1[,covariates]
  X2 = df2[,covariates]
  
  theta1 = mean(df1$infl)
  theta2 = mean(df2$infl)
  
  stab.x.shift = 0
  cov.count = 0
  if (length(covariates)==1){
    if (sd(X1)>0){
      stab.x.shift = abs((mean(X1)-mean(X2)) / sd(X1)) 
      cov.count = cov.count + 1
    }else{
      stab.x.shift = 0
    }
    
  }else{
    for (j in 1:ncol(X1)){
      if (sd(X1[,j])>0){
        x.rat.j = (mean(X1[,j]) - mean(X2[,j]))/ sd(X1[,j])
        cov.count = cov.count + 1
      }else{ 
        x.rat.j = 0
      }  
      stab.x.shift = stab.x.shift + x.rat.j^2
    }
    stab.x.shift = sqrt(stab.x.shift / cov.count)
  }
  
  
  
  return(list("theta1" = theta1, 
              "theta2" = theta2,
              "delta_total" = theta2-theta1, 
              "stab_x_shift" = stab.x.shift,
              "infl_ratio" = (mean(df1$infl) - mean(df2$infl))/sd(df1$infl)))
  
}


## =========================================================================== ##
# compute distribution shift measures when generalizing from df1 (P) to df2 (Q) #
# returns all intermediate quantities 
# for computing ratios and distribution shift measures
## =========================================================================== ##

compute_stable <- function(df1, df2, K=5, algorithm = 'grf', estimator = 'ebal',
                           covariates = NULL,
                           analysis_formula = NULL, 
                           treatment_variable = NULL, 
                           focal_variable = NULL,
                           covariate_formula = NULL, 
                           center = FALSE, eps_w=0.05){
  if (is.null(focal_variable)){
    focal_variable = treatment_variable
  }
  if (is.null(covariate_formula)){
    covariate_formula = paste("~(", paste(covariates, collapse = "+"), ")")
  }
  
  lm.1 = lm(analysis_formula, data = df1)
  lm.2 = lm(analysis_formula, data = df2)
  X.1 = model.matrix(formula(analysis_formula), data=df1)
  y.1 = model.response(model.frame(formula(analysis_formula), data=df1))
  X.2 = model.matrix(formula(analysis_formula), data=df2)
  y.2 = model.response(model.frame(formula(analysis_formula), data=df2))
  df1$infl = t(solve(t(X.1) %*% X.1 / nrow(X.1)) %*%  t(diag(y.1) %*% X.1))[,treatment_variable]
  df2$infl = t(solve(t(X.2) %*% X.2 / nrow(X.2)) %*%  t(diag(y.2) %*% X.2))[,treatment_variable]
  df1$infl_centered = df1$infl - mean(df1$infl)
  df2$infl_centered = df2$infl - mean(df2$infl)
  
  K = max(K, 2) # avoid trivial K
  n1 = nrow(df1)
  n2 = nrow(df2) 
  split.1 = split.index(n1, K)
  split.2 = split.index(n2, K)
  idx.1 = split.1$fold
  idx.1.comp = split.1$comp
  idx.2 = split.2$fold
  idx.2.comp = split.2$comp
  
  if (estimator == 'double'){ 
    # compute parameter shift components (doubly robust estimator)
    phi.pred.1 = c()
    phi.pred.2 = c()
    w.pred = c()
    for (k in 1:K){
      if (length(covariates) == 1){
        df.w = data.frame("x" = c(df1[idx.1.comp[[k]],covariates], df2[idx.2.comp[[k]],covariates]),
                          "y" = c(rep(0, length(idx.1.comp[[k]])), rep(1, length(idx.2.comp[[k]]))))
        
        if (center){
          df.phi = data.frame("x" = df1[idx.1.comp[[k]], covariates], "infl" = df1$infl_centered[idx.1.comp[[k]]]) 
        }else{
          df.phi = data.frame("x" = df1[idx.1.comp[[k]], covariates], "infl" = df1$infl[idx.1.comp[[k]]]) 
        }
        
        
        if (algorithm == 'grf'){
          
          if (center){
            c.phi = regression_forest(data.frame("x" = df.phi$x), df1$infl_centered[idx.1.comp[[k]]])
          }else{
            c.phi = regression_forest(data.frame("x" = df.phi$x), df1$infl[idx.1.comp[[k]]])
          }
          
          c.w = regression_forest(df.w$x, df.w$y) 
          phi.pred.1 = c(phi.pred.1, predict(c.phi, data.frame("x" = df1[idx.1[[k]], covariates]))$predictions)
          phi.pred.2 = c(phi.pred.2, predict(c.phi, data.frame("x" = df2[idx.2[[k]], covariates]))$predictions)
          w.pred.k = predict(c.w, data.frame("x" = df1[idx.1[[k]],covariates]))$predictions
          w.pred.k = pmax(pmin(1-eps_w, w.pred.k / (1-w.pred.k)), eps_w)
          w.pred = c(w.pred, w.pred.k) 
          
        }else if (algorithm == 'loess'){
          
          c.phi = suppressWarnings(loess(infl~., data = df.phi, span=1))
          phi.pred.1 = c(phi.pred.1, predict(c.phi, df1[idx.1[[k]],covariates]))
          phi.pred.2 = c(phi.pred.2, predict(c.phi, df2[idx.2[[k]],covariates]))
          
          c.w = suppressWarnings(loess(y~., data = df.w, span=1))
          w.pred.k = predict(c.w, newdata = df1[idx.1[[k]], covariates])
          w.pred.k = pmax(pmin(1-eps_w, w.pred.k / (1-w.pred.k)), eps_w)
          w.pred = c(w.pred, w.pred.k) 
          
        } 
      }else{
        df.w = rbind(df1[idx.1.comp[[k]],covariates], df2[idx.2.comp[[k]],covariates])
        df.w$y = c(rep(0, length(idx.1.comp[[k]])), rep(1, length(idx.2.comp[[k]])))
        
        df.phi = df1[idx.1.comp[[k]], covariates]
        if (center){
          df.phi$infl = df1$infl_centered[idx.1.comp[[k]]]
        }else{
          df.phi$infl = df1$infl[idx.1.comp[[k]]]
        }
        
        
        if (algorithm == 'grf'){
          
          if (center){
            c.phi = regression_forest(df1[idx.1.comp[[k]],covariates], df1$infl_centered[idx.1.comp[[k]]]) 
          }else{
            c.phi = regression_forest(df1[idx.1.comp[[k]],covariates], df1$infl[idx.1.comp[[k]]]) 
          }
          
          phi.pred.1 = c(phi.pred.1, predict(c.phi, df1[idx.1[[k]], covariates])$predictions)
          phi.pred.2 = c(phi.pred.2, predict(c.phi, df2[idx.2[[k]], covariates])$predictions)
          
          c.w = regression_forest(df.w[,covariates], df.w$y) 
          w.pred.k = predict(c.w, df1[idx.1[[k]],covariates])$predictions
          w.pred.k = pmax(pmin(1-eps_w, w.pred.k / (1-w.pred.k)), eps_w)
          w.pred = c(w.pred, w.pred.k)  
          
        }else if (algorithm == 'loess'){
          
          c.phi = suppressWarnings(loess(infl~., data = df.phi, span=1))
          phi.pred.1 = c(phi.pred.1, predict(c.phi, df1[idx.1[[k]],covariates]))
          phi.pred.2 = c(phi.pred.2, predict(c.phi, df2[idx.2[[k]],covariates]))
          
          c.w = suppressWarnings(loess(y~., data = df.w, span=1))
          w.pred.k = predict(c.w, newdata = df1[idx.1[[k]], covariates])
          w.pred.k = pmax(pmin(1-eps_w, w.pred.k / (1-w.pred.k)), eps_w)
          w.pred = c(w.pred, w.pred.k)
          
        } 
      }
      
    } 
    
    w.pred = w.pred * nrow(df1) / nrow(df2)
     
    theta.x = mean(phi.pred.2) - mean(phi.pred.1 * w.pred) + mean(w.pred * df1$infl)  
    
    delta.yx = mean(df2$infl) - theta.x
    delta.x = theta.x - mean(df1$infl)
    theta1 = mean(df1$infl)
    theta2 = mean(df2$infl)
  } 
  
  # covariate balancing estimator
  if (estimator == 'ebal'){
    df1$myid = 1:nrow(df1)
    df2$myid = 1:nrow(df2)
    
    fit = decomposition_helper(
      data1 = df2, # "orginal study, transfer to"
      data2 = df1, # "replication study, transfer from"
      analysis_formula = analysis_formula,
      focal_variable = focal_variable,
      covariate_formula = covariate_formula,  
      cluster_id = "myid",
      selection_variable = focal_variable
    )
    delta.x = as.numeric(fit$decomp['Covariates'])
    delta.yx = as.numeric(fit$decomp['Residual'])
    theta2 = as.numeric(fit$decomp['Theta_org'])
    theta1 = as.numeric(fit$decomp['Theta_rep'])
    
    # prepare for estimating variance
    phi.pred.1 = c() 
    for (k in 1:K){  
      if (length(covariates)==1){
        
        if (center){
          df.phi = data.frame("x" = df1[idx.1.comp[[k]], covariates], "infl" = df1$infl_centered[idx.1.comp[[k]]]) 
        }else{
          df.phi = data.frame("x" = df1[idx.1.comp[[k]], covariates], "infl" = df1$infl[idx.1.comp[[k]]]) 
        }
        
        if (algorithm == 'grf'){
          c.phi = regression_forest(data.frame("x" = df.phi$x), df.phi$infl)  
          phi.pred.1 = c(phi.pred.1, predict(c.phi, data.frame("x"=df1[idx.1[[k]], covariates]))$predictions) 
        }else if (algorithm == 'loess'){
          c.phi = suppressWarnings(loess(infl~., data = df.phi, span=1)) 
          phi.pred.1 = c(phi.pred.1, predict(c.phi, df1[idx.1[[k]],covariates])) 
        } 
        
      }else{
        df.phi = df1[idx.1.comp[[k]], covariates]
        
        if (center){
          df.phi$infl = df1$infl_centered[idx.1.comp[[k]]]
        }else{
          df.phi$infl = df1$infl[idx.1.comp[[k]]]
        }
        
        if (algorithm == 'grf'){
          if (center){
            c.phi = regression_forest(df1[idx.1.comp[[k]], covariates], df1$infl_centered[idx.1.comp[[k]]])  
          }else{
            c.phi = regression_forest(df1[idx.1.comp[[k]], covariates], df1$infl[idx.1.comp[[k]]])  
          }
          
          phi.pred.1 = c(phi.pred.1, predict(c.phi, df1[idx.1[[k]], covariates])$predictions) 
        }else if (algorithm == 'loess'){
          c.phi = suppressWarnings(loess(infl~., data = df.phi, span=1)) 
          phi.pred.1 = c(phi.pred.1, predict(c.phi, df1[idx.1[[k]],covariates])) 
        } 
      }
      
      
      
    } 
  } 
  
  # compute sensitivity measures
  if (center){
    cal.df = data.frame("x" = phi.pred.1 + mean(df1$infl), "y" = df1$infl)
  }else{
    cal.df = data.frame("x" = phi.pred.1, "y" = df1$infl)
  }
  
  cal.lm = lm(y~., data = cal.df)
  phi.pred.1.cal = predict(cal.lm)
  # s.yx.squared = min(mean((df1$infl - phi.pred.1.cal)^2), sd(df1$infl)^2)
  s.x.squared = mean(phi.pred.1.cal * (2*df1$infl - phi.pred.1.cal)) - mean(df1$infl)^2
  s.yx.squared = sd(df1$infl)^2 - s.x.squared
  s.yx = sqrt(s.yx.squared)
  s.x = sqrt(s.x.squared) 
  
  s.x.naive = sd(phi.pred.1.cal)
  s.yx.naive = sqrt(sd(df1$infl)^2 - s.x.naive^2)
  
  # compute shift strength
  cond.shift = delta.yx / s.yx
  x.shift = delta.x / s.x
  
  # stablize conditional shift 
  X1 = df1[,covariates]
  X2 = df2[,covariates]
  
  stab.x.shift = 0
  cov.count = 0
  if (length(covariates)==1){
    if (sd(X1)>0){
      stab.x.shift = abs((mean(X1)-mean(X2)) / sd(X1)) 
      cov.count = cov.count + 1
    }else{
      stab.x.shift = 0
    }
    
  }else{
    for (j in 1:ncol(X1)){
      if (sd(X1[,j])>0){
        x.rat.j = (mean(X1[,j]) - mean(X2[,j]))/ sd(X1[,j])
        cov.count = cov.count + 1
      }else{ 
        x.rat.j = 0
      }   
      stab.x.shift = stab.x.shift + x.rat.j^2
    }
    stab.x.shift = sqrt(stab.x.shift / cov.count)
  }
  
  
  
  return(list("theta1" = theta1, 
              "theta2" = theta2,
              "delta_total" = delta.yx + delta.x, 
              "delta_yx" = delta.yx, "delta_x" = delta.x, 
              "cond_sensitivity" = s.yx, "x_sensitivity" = s.x, 
              "cond_sens_naive" = s.yx.naive, "x_sens_naive" = s.x.naive,
              "cond_shift" = cond.shift, "x_shift" = x.shift,
              "stab_x_shift" = stab.x.shift,
              "infl_ratio" = (mean(df1$infl) - mean(df2$infl))/sd(df1$infl)))
}


 
## ============================================================== ##
# transfer estimator + PI uncer covariate shift (from df1 to df2)
# method = 'double': doubly robust estimator 
# method = 'ebal': entropy balancing estimator
# K is cross-fitting folder number
## ============================================================== ##

transfer_dr <- function(df1, df2, K=5, algorithm = 'grf', 
                        method = 'ebal',
                        covariates = NULL,
                        analysis_formula = NULL, 
                        treatment_variable = NULL, 
                        focal_variable = NULL,
                        covariate_formula = NULL, 
                        center=TRUE, eps_w=0.05){
  if (is.null(focal_variable)){
    focal_variable = treatment_variable
  }
  if (is.null(covariate_formula)){
    covariate_formula = paste("~(", paste(covariates, collapse = "+"), ")")
  }
  
  lm.1 = lm(analysis_formula, data = df1)
  lm.2 = lm(analysis_formula, data = df2)
  X.1 = model.matrix(formula(analysis_formula), data=df1)
  y.1 = model.response(model.frame(formula(analysis_formula), data=df1))
  X.2 = model.matrix(formula(analysis_formula), data=df2)
  y.2 = model.response(model.frame(formula(analysis_formula), data=df2))
  df1$infl = t(solve(t(X.1) %*% X.1 / nrow(X.1)) %*%  t(diag(y.1) %*% X.1))[,treatment_variable]
  df2$infl = t(solve(t(X.2) %*% X.2 / nrow(X.2)) %*%  t(diag(y.2) %*% X.2))[,treatment_variable]
  df1$infl_centered = df1$infl - mean(df1$infl)
  df2$infl_centered = df2$infl - mean(df2$infl)
  
  K = max(K, 2) # avoid trivial K
  n1 = nrow(df1)
  n2 = nrow(df2) 
  split.1 = split.index(n1, K)
  split.2 = split.index(n2, K)
  idx.1 = split.1$fold
  idx.1.comp = split.1$comp
  idx.2 = split.2$fold
  idx.2.comp = split.2$comp
  
  # compute parameter shift components (doubly robust estimator)
  phi.pred.1 = c()
  phi.pred.2 = c()
  w.pred = c()
  for (k in 1:K){
    if (length(covariates) == 1){
      # preparation
      df.w = data.frame("x" = c(df1[idx.1.comp[[k]],covariates], df2[idx.2.comp[[k]],covariates]),
                        "y" = c(rep(0, length(idx.1.comp[[k]])), rep(1, length(idx.2.comp[[k]]))))
      if (center){
        df.phi = data.frame("x" = df1[idx.1.comp[[k]], covariates], "infl" = df1$infl_centered[idx.1.comp[[k]]]) 
      }else{
        df.phi = data.frame("x" = df1[idx.1.comp[[k]], covariates], "infl" = df1$infl[idx.1.comp[[k]]]) 
      } 
      
      if (algorithm == 'grf'){
        if (center){
          c.phi = regression_forest(data.frame("x" = df.phi$x), df1$infl_centered[idx.1.comp[[k]]])
        }else{
          c.phi = regression_forest(data.frame("x" = df.phi$x), df1$infl[idx.1.comp[[k]]])
        }
        
        c.w = regression_forest(df.w$x, df.w$y)
        
        phi.pred.1 = c(phi.pred.1, predict(c.phi, data.frame("x" = df1[idx.1[[k]], covariates]))$predictions)
        phi.pred.2 = c(phi.pred.2, predict(c.phi, data.frame("x" = df2[idx.2[[k]], covariates]))$predictions)
        w.pred = c(w.pred, predict(c.w, data.frame("x" = df1[idx.1[[k]],covariates]))$predictions) 
        
        
      }else if (algorithm == 'loess'){
        
        c.phi = suppressWarnings(loess(infl~., data = df.phi, span=1))
        c.w = suppressWarnings(loess(y~., data = df.w, span=1))
        phi.pred.1 = c(phi.pred.1, predict(c.phi, df1[idx.1[[k]],covariates]))
        phi.pred.2 = c(phi.pred.2, predict(c.phi, df2[idx.2[[k]],covariates]))
        w.pred = c(w.pred, predict(c.w, newdata = df1[idx.1[[k]], covariates]))
        
      } 
      
    }else{
      
      df.w = rbind(df1[idx.1.comp[[k]],covariates], df2[idx.2.comp[[k]],covariates])
      df.w$y = c(rep(0, length(idx.1.comp[[k]])), rep(1, length(idx.2.comp[[k]])))
      
      df.phi = df1[idx.1.comp[[k]], covariates]
      
      if (center){
        df.phi$infl = df1$infl_centered[idx.1.comp[[k]]]
      }else{
        df.phi$infl = df1$infl[idx.1.comp[[k]]]
      } 
      
      if (algorithm == 'grf'){
        if (center){
          c.phi = regression_forest(df1[idx.1.comp[[k]],covariates], df1$infl_centered[idx.1.comp[[k]]]) 
        }else{
          c.phi = regression_forest(df1[idx.1.comp[[k]],covariates], df1$infl[idx.1.comp[[k]]]) 
        }
        phi.pred.1 = c(phi.pred.1, predict(c.phi, df1[idx.1[[k]], covariates])$predictions)
        phi.pred.2 = c(phi.pred.2, predict(c.phi, df2[idx.2[[k]], covariates])$predictions)
        
        c.w = regression_forest(df.w[,covariates], df.w$y) 
        w.pred.k = predict(c.w, df1[idx.1[[k]],covariates])$predictions
        w.pred.k = pmax(pmin(1-eps_w, w.pred.k / (1-w.pred.k)), eps_w)
        w.pred = c(w.pred, w.pred.k) 
        
      }else if (algorithm == 'loess'){
        
        c.phi = suppressWarnings(loess(infl~., data = df.phi, span=1))
        phi.pred.1 = c(phi.pred.1, predict(c.phi, df1[idx.1[[k]],covariates]))
        phi.pred.2 = c(phi.pred.2, predict(c.phi, df2[idx.2[[k]],covariates]))
        
        c.w = suppressWarnings(loess(y~., data = df.w, span=1))
        w.pred.k = predict(c.w, newdata = df1[idx.1[[k]], covariates])
        w.pred.k = pmax(pmin(1-eps_w, w.pred.k / (1-w.pred.k)), eps_w)
        w.pred = c(w.pred, w.pred.k)
        
      } 
    }
    
  } 
  
  
  w.pred = w.pred * nrow(df1) / nrow(df2)
  
  if (center){
    theta.x = mean(df1$infl) - mean(phi.pred.1 * w.pred) + mean(w.pred * df1$infl_centered)
  }else{
    theta.x = mean(df1$infl) - mean(phi.pred.1 * w.pred) + mean(w.pred * df1$infl)
  }
  
  
  delta.yx = mean(df2$infl) - theta.x
  delta.x = theta.x - mean(df1$infl)
  theta1 = mean(df1$infl)
  theta2 = mean(df2$infl) 
  
  # compute sensitivity measures
  if (center){
    cal.df = data.frame("x" = phi.pred.1 + mean(df1$infl), "y" = df1$infl)
  }else{
    cal.df = data.frame("x" = phi.pred.1, "y" = df1$infl)
  }
  
  cal.lm = lm(y~., data = cal.df)
  phi.pred.1.cal = predict(cal.lm) # linear calibration
  s.yx.squared = min(mean((df1$infl - phi.pred.1.cal)^2), sd(df1$infl)^2) # conditional variance
  s.x.squared = mean(phi.pred.1.cal * (2*df1$infl - phi.pred.1.cal)) - mean(df1$infl)^2 # covariate variance
  s.yx = sqrt(s.yx.squared)
  s.x = sqrt(s.x.squared) 
  
  # transfer variance 
  trans.var.2 = min(mean(w.pred * (df1$infl - phi.pred.1.cal)^2), mean(w.pred * (df1$infl - mean(df1$infl))^2)) / nrow(df2) 
  trans.var.1 = sd(w.pred * (df1$infl - phi.pred.1.cal))^2 / nrow(df1)
  trans.var = trans.var.2 + trans.var.1 
  
  return(list("theta1" = theta1, 
              "theta2" = theta2,
              # "hat_theta_x" = theta.x, 
              "hat_var" = trans.var,
              "alg" = method
  ))
}

## ============================================ ##
# transfer estimator + PI uncer iid assumption 
# transfer from df1 to df2 
## ============================================ ##  
  
transfer_plain <- function(df1, df2, 
                           analysis_formula = NULL, 
                           treatment_variable = NULL){
  
  lm.1 = lm(analysis_formula, data = df1)
  lm.2 = lm(analysis_formula, data = df2)
  X.1 = model.matrix(formula(analysis_formula), data=df1)
  y.1 = model.response(model.frame(formula(analysis_formula), data=df1))
  X.2 = model.matrix(formula(analysis_formula), data=df2)
  y.2 = model.response(model.frame(formula(analysis_formula), data=df2))
  infl.1 = t(solve(t(X.1) %*% X.1 / nrow(X.1)) %*%  t(diag(y.1) %*% X.1))
  infl.2 = t(solve(t(X.2) %*% X.2 / nrow(X.2)) %*%  t(diag(y.2) %*% X.2)) 
  
  resid.1 = (as.numeric(y.1) - as.numeric(X.1 %*% as.matrix(as.numeric(colMeans(infl.1)), ncol=1)))
  resid.2 = (as.numeric(y.2) - as.numeric(X.2 %*% as.matrix(as.numeric(colMeans(infl.2)), ncol=1)))
  c.infl.1 = t(solve(t(X.1) %*% X.1 / nrow(X.1)) %*% t(X.1) %*% diag(resid.1))
  c.infl.2 = t(solve(t(X.2) %*% X.2 / nrow(X.2)) %*% t(X.2) %*% diag(resid.2))
  
  df1$infl = c.infl.1[,treatment_variable]
  df2$infl = c.infl.2[,treatment_variable]
  
  sd.1 = diag(sqrt(cov(c.infl.1)))[treatment_variable]
  
  theta1 = mean(infl.1[,treatment_variable]) 
  theta2 = mean(infl.2[,treatment_variable])
  # transfer variance   
  trans.var = sd.1^2 * (1/nrow(df1) + 1/nrow(df2))
  
  return(list("hat_theta" = theta1, 
              "theta1" = theta1,
              "theta2" = theta2,
              "hat_var" = trans.var
  ))
}

## ============================================ ##  
# other utility functions

 
split.index <- function(n, K=2, seed=NULL){
  if (!(is.null(seed))){
    set.seed(seed)
  }
  idx = sample(1:n)
  fold.idx = list()
  comp.idx = list()
  for (k in 1:K){
    if (k==1){
      fold.idx[[k]] = idx[1:floor(n/K)]
      comp.idx[[k]] = setdiff(1:n, fold.idx[[k]])
    }else if (k==K){
      fold.idx[[k]] = idx[(floor(n*(K-1)/K)+1):n]
      comp.idx[[k]] = setdiff(1:n, fold.idx[[k]])
    }else{
      fold.idx[[k]] = idx[(floor(n*(k-1)/K)+1): floor(n*k/K)]
      comp.idx[[k]] = setdiff(1:n, fold.idx[[k]])
    }
    
  }
  return(list("fold" = fold.idx, "comp" = comp.idx))
}


run_diagnosis <- function(
    org_df,
    rep_df,
    analysis_formula,
    treatment_variable, # column name for the treatment indicator (1 = treatment)
    focal_variable = NULL, # focal variable
    covariates = NULL,
    covariate_formula = NULL,
    mediators = NULL,
    mediation_formula = NULL,
    cluster_id = NULL,
    selection_variable = NULL,
    alpha = 0.1,
    verbose = TRUE,
    if_selective = TRUE,
    pub_pvalue_threshold = 0.05
){
  
  df1 = data.frame(org_df)
  df2 = data.frame(rep_df)
  
  
  
  #####################################################
  ######### sanity check for covariate inputs #########
  #####################################################
  
  # if covariate formula is not provided, check whether the covariates argument is meaningful
  # otherwise just use the covariate formula
  if (is.null(covariate_formula)){
    if (length(covariates) > 0){
      for (i.par in 1:length(covariates)){
        if (!is.na(suppressWarnings(as.integer(covariates[i.par])))){
          covariates[i.par] = colnames(org_df)[as.integer(covariates[i.par])]
        }
      }
      # if (is.null(treatment_variable)){
      #   stop("Please provide the treatment variable argument!")
      # }
      covariate_formula = paste("~(", paste(covariates, collapse = "+"), ") *", treatment_variable)
    }
  }
  
  # extract covariates and check whether it behaves normally
  if (!is.null(covariate_formula)){
    X1 = tryCatch({model.matrix(formula(covariate_formula), data=org_df)},
                  error = function(e) { return(NA) })
    X2 = tryCatch({model.matrix(formula(covariate_formula), data=rep_df)},
                  error = function(e) { return(NA) })
    if (is.na(X1) || is.na(X2)){ # any of the data fails to be extracted
      stop("Covariate input does not work! \n")
    }else{# extract covariate names into a vector
      covariates = setdiff(intersect(colnames(org_df), colnames(X1)[2:ncol(X1)]), c(treatment_variable))
      if (is.null(dim(X1[,2:ncol(X1)]))){
        if (sd(X1[,2:ncol(X1)])==0){stop("Constant covariate error! \n")}
      }
    }
  }
  
  #####################################################
  ######### sanity check for mediation inputs #########
  #####################################################
  
  # if mediation formula is not provided, check whether the mediators argument is meaningful
  if (is.null(mediation_formula)){
    if (length(mediators) > 0){
      # if (is.null(treatment_variable)){
      #   stop("Please provide the treatment variable argument!")
      # }
      
      for (i.par in 1:length(mediators)){
        if (!is.na(suppressWarnings(as.integer(mediators[i.par])))){
          mediators[i.par] = colnames(org_df)[as.integer(mediators[i.par])]
        }
      }
      
      cov_and_med = union(covariates, mediators)
      # generate mediation formula
      mediation_formula = paste("~(", paste(cov_and_med, collapse="+"), ")*", treatment_variable)
      
    }
  }
  
  # extract mediators and check whether it behaves normally
  if (!is.null(mediation_formula)){
    X1 = tryCatch({model.matrix(formula(mediation_formula), data=df1)},
                  error = function(e) { return(NA) })
    X2 = tryCatch({model.matrix(formula(mediation_formula), data=df2)},
                  error = function(e) { return(NA) })
    if (is.na(X1) || is.na(X2)){ # any of the data fails to be extracted
      stop("Mediation input does not work! \n")
    }else{# extract covariate names into a vector
      mediators = setdiff(intersect(colnames(org_df), colnames(X1)[2:ncol(X1)]),
                          c(treatment_variable, covariates))
      if (is.null(dim(X1[,2:ncol(X1)]))){
        if (sd(X1[,2:ncol(X1)])==0){stop("Constant mediator error! \n")}
      }
    }
  }
  
  #####################################################
  ######### sanity check for other inputs #########
  #####################################################
  
  if (is.null(cluster_id)) {
    df1$idformyuse = 1:nrow(df1)
    df2$idformyuse = 1:nrow(df2)
    cluster_id = "idformyuse"
  }
  
  no_covariate_formula = is.null(covariate_formula) || covariate_formula==""
  no_mediation_formula = is.null(mediation_formula) || mediation_formula==""
  
  if (no_covariate_formula & no_mediation_formula) {
    warning("covariate_formula and mediation_formula cannot both be NULL.")
    return()
  }
  
  if (is.null(focal_variable)){
    focal_variable = treatment_variable
  }
  
  # selection variable equals the focal if not specified
  if (is.null(selection_variable)){
    selection_variable = focal_variable
  }
  
  fit = decomposition_helper(
    data1 = df1,
    data2 = df2,
    analysis_formula = analysis_formula,
    focal_variable = focal_variable,
    covariate_formula = covariate_formula,
    mediation_formula = mediation_formula,
    cluster_id = cluster_id,
    selection_variable = selection_variable
  )
  
  results = fit$decomp
  clusters1 = unique(df1[[cluster_id]])
  clusters2 = unique(df2[[cluster_id]])
  
  
  N1 = length(clusters1)
  N2 = length(clusters2)
  phi1 = phi2 = data.frame()
  pb = progress_bar$new("Jackknifing [:bar] :percent", width=80, total=N1 + N2+1)
  
  # Add progress bar
  total_steps = N1+N2 + 1
  
  
  for (i in 1:N1) {
    if (verbose) { pb$tick() }
    
    tryCatch({phi1 = bind_rows(phi1, decomposition_helper(
      data1 = df1[df1[[cluster_id]] != clusters1[i],],
      data2 = df2,
      analysis_formula = analysis_formula,
      focal_variable = focal_variable,
      covariate_formula = covariate_formula,
      mediation_formula = mediation_formula,
      cluster_id = cluster_id,
      selection_variable = selection_variable
    )$decomp)}, error = function(e) {})
  }
  for (i in 1:N2) {
    if (verbose) { pb$tick() }
    
    tryCatch({phi2 = bind_rows(phi2, decomposition_helper(
      data1 = df1,
      data2 = df2[df2[[cluster_id]] != clusters2[i],],
      analysis_formula = analysis_formula,
      focal_variable = focal_variable,
      covariate_formula = covariate_formula,
      mediation_formula = mediation_formula,
      cluster_id = cluster_id,
      selection_variable = selection_variable
    )$decomp)}, error = function(e) {})
  }
  
  SEs = sqrt(apply(phi1, 2, var)*(nrow(phi1)-1)^2/N1 + apply(phi2, 2, var)*(nrow(phi2)-1)^2/N2)
  Sigma = var(phi1)*(nrow(phi1)-1)^2/N1 + var(phi2)*(nrow(phi2)-1)^2/N2
  
  # ==================================================
  # Plotting results
  
  # non-selective results
  pvals = 2 - 2 * pnorm(abs(as.numeric(fit$decomp))/SEs)
  ret_table = cbind(as.numeric(fit$decomp), SEs, as.numeric(fit$decomp)/SEs, pvals)
  colnames(ret_table) = c("Estimate", "Std. Error", "t-stat", "Pr(>|z|)")
  ret_table = ret_table[rownames(ret_table) != "Selected",]
  rownames(ret_table)[1] = "Observed discrepancy"
  ret_table_df <- as.data.frame(ret_table)
  ret_table_df$Component <- rownames(ret_table)
  ret_table_df <- ret_table_df[, c("Component", "Estimate", "Std. Error", "t-stat", "Pr(>|z|)")]
  
  clean_decomp = select(fit$decomp, -Selected)
  clean_SE = SEs[2:length(SEs)]
  clean_low = clean_decomp
  clean_low[2:length(clean_low)] = clean_low[2:length(clean_low)] - qnorm(1-alpha/2) * clean_SE[2:length(clean_SE)]
  clean_high = clean_decomp
  clean_high[2:length(clean_high)] = clean_high[2:length(clean_high)] + qnorm(1-alpha/2) * clean_SE[2:length(clean_SE)]
  
  component = c("Observed discrepancy")
  if (!no_covariate_formula){
    component = c(component, "+ Covariate shift")
  }
  if (!no_mediation_formula){
    component = c(component, "+ Mediation shift")
  }
  component = c(component, "+ Residual", "= Sampling variability")
  bounds_for_plot = data.frame("low" = c(as.numeric(clean_low),
                                         -qnorm(1-alpha/2) * (as.numeric(clean_SE[1]))),
                               "high" = c(as.numeric(clean_high),
                                          qnorm(1-alpha/2) * (as.numeric(clean_SE[1]))),
                               "estimate" = c(as.numeric(clean_decomp),0),
                               "component" = component)
  suppressWarnings({
    decomp.plot = discrepancy_plot(bounds_for_plot, alpha=0.1, base_size=28, digital=2)
  })
  
  # selective part
  decomp.plot.sel = NULL
  sel.message = NULL
  ret_table_sel_df = NULL
  if (if_selective){
    sel.message = "NOTE:"
    
    if ((as.numeric(pvals['Selected']) > pub_pvalue_threshold)){
      sel.message = paste(sel.message,
                          "Original p-value > publication bias threshold, stop running selective inference! \n")
      pb$tick()
      return(list("plot" = decomp.plot, "table" = ret_table_df,
                  "sel.plot" = NULL, "sel.table" = NULL,
                  "sel.message" = sel.message))
    }else{
      sel.message = paste(sel.message,
                          paste("Running selective inference conditional on p <= ", pub_pvalue_threshold,"! \n",sep=''))
      bounds_sel = selective_decomposition(fit$decomp, Sigma,
                                           qnorm(1-pub_pvalue_threshold/2) * sqrt(Sigma["Selected", "Selected"]),
                                           alpha=alpha)
    }
    
    
    # generate table
    ret_table_sel = bounds_sel[1:(nrow(bounds_sel)-1),]
    rownames(ret_table_sel) = ret_table_sel$component
    rownames(ret_table_sel)[2:nrow(ret_table_sel)] = sapply(rownames(ret_table_sel)[2:nrow(ret_table_sel)],
                                                            function(x) substr(x, 3, nchar(x)))
    ret_table_sel_df = data.frame(ret_table_sel)
    
    # selective table
    ret_table_sel_df$Component <- rownames(ret_table_sel_df)
    ret_table_sel_df <- ret_table_sel_df[, c("Component", "low", "high", "estimate")]
    colnames(ret_table_sel_df) = c("Component", "CI.Low", "CI.High", "Estimate")
    ret_table_sel_df$Estimate = as.numeric(ret_table_sel_df$Estimate)
    
    ret_table_sel_df$Estimate[ret_table_sel_df$Estimate==-Inf] = NULL
    ret_table_sel_df$Estimate[ret_table_sel_df$Estimate==Inf] = NULL
    
    # generate error message
    components = rownames(ret_table_sel_df)
    
    for (v in 1:length(components)){
      
      if ((ret_table_sel_df$CI.Low[v] == -Inf) || (ret_table_sel_df$CI.High[v] == Inf)){
        sel.message = paste(sel.message,
                            "No meaningful selective CI within est +- 256*se for ",
                            components[v], "component! \n")
      }
      
    }
    
    
    # generate selective decomposition plot
    bounds_for_plot_sel = data.frame(bounds_sel) %>% drop_na()
    if (nrow(bounds_for_plot_sel) >0 ){
      if (("Observed discrepancy" %in% bounds_for_plot_sel$component)){
        bounds_for_plot_sel$low[1] = bounds_for_plot_sel$estimate[1]
        bounds_for_plot_sel$high[1] = bounds_for_plot_sel$estimate[1]
      }
      suppressWarnings({
        decomp.plot.sel = discrepancy_plot(bounds_for_plot_sel, alpha=alpha, base_size=28, digital=2)
      })
      
    }else{
      decomp.plot.sel = NULL
    }
  }
  
  if (verbose){
    cat("\n")
    cat("\n")
    cat("Message for selective inferece:\n")
    cat(sel.message)
  }
  
  
  pb$tick()
  
  ret_table_df = select(ret_table_df, -Component) 
  ret_table_df = ret_table_df %>% mutate_if(is.numeric, function(x) round(x,4))
  
  if (if_selective){
    ret_table_sel_df = select(ret_table_sel_df, -Component)
    ret_table_sel_df = ret_table_sel_df %>% mutate_if(is.numeric, function(x) round(x,4))
  }
  
  return(list("plot" = decomp.plot, "table" = ret_table_df,
              "sel.plot" = decomp.plot.sel, "sel.table" = ret_table_sel_df,
              "sel.message" = sel.message))
  
}


decomposition_helper = function(
    data1, # original data
    data2, # replication data
    analysis_formula, # a formula for OLS
    focal_variable, # the name of the focal variable
    covariate_formula = NULL, # covariates to be balanced
    # mediation_formula = NULL, # covariates + mediators to be balanced
    cluster_id, # variable to cluster on. well-defined before calling this function
    selection_variable = NULL # variable on which the selection happens
) {
  
  
  group1 = data1 %>% group_by_at(cluster_id) %>% slice(1)
  group2 = data2 %>% group_by_at(cluster_id) %>% slice(1)
  
  
  no_covariate_formula = is.null(covariate_formula) || covariate_formula==""
  # # no_mediation_formula = is.null(mediation_formula) || mediation_formula==""
  # 
  # 
  # # Create design matrices.  Care to drop unused factor levels, for the jackknife
  # if (!no_covariate_formula) { 
    
  X1 = model.matrix.lm(formula(covariate_formula), data=group1, na.action = "na.pass")[,-1]
  if (is.null(dim(X1))){
    X1 = data.frame("X" = as.numeric(X1))
    colnames(X1)[1] = colnames(model.matrix(formula(covariate_formula), data=group1))[2]
    X2 = model.matrix(formula(covariate_formula), data=group2)[,colnames(X1)] 
    X2 = data.frame("X" = as.numeric(X2))
    colnames(X2) = colnames(X1)
  }else{
    X1 = X1[,apply(X1, 2, function(x) { !all(x==0) })]
    X2 = model.matrix.lm(formula(covariate_formula), data=group2, na.action = "na.pass")[,colnames(X1)] 
  }
  
  group1 = group1[rowSums(is.na(X1)) ==0, ]
  group2 = group2[rowSums(is.na(X2)) ==0, ]
  X1 = X1[rowSums(is.na(X1)) ==0, ]
  X2 = X2[rowSums(is.na(X2)) ==0, ]
  
  # standardize 
  mean.1 = colMeans(X1)
  std.1 = sapply(1:ncol(X1), function(j) sd(X1[,j]))
  X1 = (X1 - matrix(1, nrow = nrow(X1), ncol=1) %*% matrix(mean.1,nrow=1)) %*% diag(1/std.1)
  X2 = (X2 - matrix(1, nrow = nrow(X2), ncol=1) %*% matrix(mean.1,nrow=1)) %*% diag(1/std.1)
  
  X = rbind(X1, X2)
   
  
  # Fitting covariate weights
  S = c(rep(1, nrow(X1)), rep(0, nrow(X2)))
  if (no_covariate_formula) {
    group2$w = rep(1, nrow(X2))
  } else {
    group2$w = ebalance(Treatment=S, X=X, print.level=-1)$w
  } 
  data2_joined = left_join(data2, group2[,c(cluster_id, "w" )], by = c(cluster_id))  
  # Fitting mediation weights
  # if (no_mediation_formula) {
  #   group2$o = group2$w
  # } else {
  #   group2$o = ebalance(Treatment=S, X=M, print.level=-1)$w
  # } 
  # data2_joined = left_join(data2, group2[,c(cluster_id, "w", "o")], by = c(cluster_id)) 
  
  # Re-performing analysis
  selection = coef(lm(formula(analysis_formula), data=data1))[selection_variable]
  theta1 = coef(lm(formula(analysis_formula), data=data1))[focal_variable] 
  theta2 = coef(lm(formula(analysis_formula), data=data2_joined))[focal_variable]  
  theta2X= coef(lm(formula(analysis_formula), data=data2_joined, weights=w))[focal_variable]  
  # theta2M= coef(lm(formula(analysis_formula), data=data2_joined, weights=data2_joined$o))[focal_variable] 
  
  # Computing decomposition
  decomp = data.frame(
    "Theta_org" = theta1, 
    "Theta_rep" = theta2,
    "Selected" = selection,
    "Observed" = theta1 - theta2,
    "Covariates" = theta2X - theta2,
    "Residual" = theta1 - theta2X
    # "Mediators" = theta2M - theta2X,
    # "Residual" = theta1 - theta2M
  )
  if (no_covariate_formula) {
    decomp = select(decomp, -Covariates)
  }
  # if (no_mediation_formula) {
  #   decomp = select(decomp, -Mediators)
  # }
  
  return(list("decomp" = decomp, "w" = group2$w#, "o" = group2$o
              ))
  
}




selective_decomposition = function(decomp, Sigma, threshold, alpha=0.1) {
  # Selective confidence intervals and conditional MLEs
  x = decomp#$decomp
  S = Sigma#$Sigma
  observed = x[["Observed"]]
  
  # Compute selection-adjusted confidence intervals for each component
  components = setdiff(names(x), "Selected")
  bounds = lapply(
    X = components,
    FUN = function(v) { vars = c(v, "Selected"); selective_ci(x[vars], S[vars, vars], threshold, alpha)}
  ) %>% bind_rows()
  rownames(bounds) = components
  
  # Compute selection-adjusted point estimates for the non-residual components 
  vars = setdiff(components[which(bounds$low != -Inf)], "Residual") 
  SX = S[vars, vars]
  SXY = S[vars, "Selected"]
  Omega = solve(SX)
  beta = as.numeric(Omega %*% SXY)
  sigma = as.numeric(sqrt(t(beta) %*% SX %*% beta))
  R = as.numeric(x["Selected"] - sum(beta*x[vars]))
  
  suppressWarnings({
    estimate = optim(
      par = as.numeric(unlist((bounds["high"] - bounds["low"])/2))[setdiff(which(bounds$low!=-Inf), length(components))], # Midpoint initialization
      fn = function(par) {
        0.5*t(as.numeric(x[vars]-par)) %*% Omega %*% as.numeric(x[vars]-par) +
          truncprob(mean = sum(par*beta), sd = sigma, R = R, threshold = threshold, log.prob = T)
      }
    )$par # Maximum conditional likelihood
  })
  
  if (length(estimate >1)){
    names(estimate) = vars
  }else{
    estimate = c(estimate,0)
    names(estimate) = c(vars, "placeholder")
  }
  names(estimate) = vars
  bounds$estimate = as.numeric(decomp[components]) 
  
  bounds$estimate[which(bounds$low!=-Inf)] = sapply(
    X = components[which(bounds$low!=-Inf)],
    FUN = function(comp) {
      if (comp == "Observed") {
        return(observed)
      } else if (comp == "Residual") {
        return(as.numeric(2*estimate["Observed"] - sum(estimate)))
      } else {
        return(estimate[[comp]])
      }
    }
  )
  
  # Compute the estimate of sampling variability
  if (bounds$low[1] != -Inf){
    bounds = rbind(bounds, observed - c(bounds["Observed", "high"], bounds["Observed", "low"], estimate["Observed"]))
  }else{
    bounds = rbind(bounds, observed - c(bounds["Observed", "high"], bounds["Observed", "low"], 0))
  }
  
  
  # Return results
  bounds$component0 = c(components, "Sampling")
  clean_names = data.frame(
    component0 = c("Observed", "Covariates", "Mediators", "Residual", "Sampling"),
    component = c("Observed discrepancy", "+ Covariate shift", "+ Mediation shift", "+ Residual", "= Sampling variability")
  )
  bounds = left_join(data.frame(bounds), clean_names, by = c("component0")) %>% select(-component0)
  
  bounds$low[is.na(bounds$low)] = -Inf
  bounds$high[is.na(bounds$high)] = Inf
  return(bounds)
}


discrepancy_plot = function(bounds, alpha=0.1, base_size=28, digital=2) {
  # Bounds should be a data frame of the type returned by "selective_decomposition"
  df_plot = bounds %>% mutate(
    level = factor(component, levels=c("+ Residual", "+ Mediation shift", "+ Covariate shift", "= Sampling variability", "Observed discrepancy"))
  )
  # Location of the numbers and line
  xnums = min(df_plot$low) - 0.2*(max(df_plot$high) - min(df_plot$low))
  xline = min(df_plot$low) - 0.6*(max(df_plot$high) - min(df_plot$low)) - 0.008*(digital-2)
  balance_sheet_plot = df_plot %>% ggplot() +
    geom_crossbar(aes(y = level, x = estimate, col = level, fill = level,
                      xmin = low, xmax = high),
                  width = 0.5, alpha = 0.5) +
    geom_text(aes(x = xnums, y = level, label=paste(ifelse(estimate>=0, "+", ""),
                                                    as.character(round(estimate, digital)), sep="")),
              hjust = 1, size = base_size*0.27) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_vline(xintercept = df_plot$estimate[1], lty=2) +
    ylab("") +
    xlab("") +
    theme_bw(base_size = base_size) +
    theme(
      legend.position = "None",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_text(face = c(rep("italic", 4), "bold")),
      axis.text.x = element_text(size=base_size*0.5),
      plot.margin = unit(c(0.7, 0.7, -0.5, -0.5), "cm")
    ) +
    annotate(geom = "segment", x = xline, xend = xline, y = -Inf, yend = Inf) +
    scale_x_continuous(breaks = c(0, round(df_plot$estimate[1], 2)))
  return(balance_sheet_plot)
}


truncprob = function(mean, sd, R, threshold, log.prob=TRUE) {
  # Probability that |N(mean, sd) + R| > threshold > 0
  upper = 1-pnorm((threshold-R-mean)/sd) # P{N(mu, sigma^2) > threshold - R}
  lower = pnorm(-(threshold+R+mean)/sd) # P{N(mu, sigma^2) < -threshold - R}
  return(ifelse(log.prob, log(upper+lower), upper+lower))
}

selective_pval = function(mu, x, Sigma, threshold) {
  x = as.numeric(x)
  threshold = as.numeric(threshold)
  beta = as.numeric(Sigma[1,2]/Sigma[1,1])
  s1 = as.numeric(sqrt(Sigma[1,1]))
  R = as.numeric(x[2] - beta*x[1])
  if (beta >= 0) {
    num =
      pnorm((pmin(-(threshold+R)/beta, x[1]) - mu)/s1) +
      pmax(0, pnorm((x[1]-mu)/s1) - pnorm(((threshold-R)/beta - mu)/s1))
    denom =
      pnorm((-(threshold+R)/beta - mu)/s1) + pnorm(((threshold-R)/beta - mu)/s1, lower.tail=F)
    if (denom==0){return(1)}else{return(num/denom)}
    # return(num/denom)
  } else {
    num =
      pnorm((pmin((threshold-R)/beta, x[1]) - mu)/s1) +
      pmax(0, pnorm((x[1]-mu)/s1) - pnorm((-(threshold+R)/beta - mu)/s1))
    denom =
      pnorm(((threshold-R)/beta - mu)/s1) + pnorm((-(threshold+R)/beta - mu)/s1, lower.tail=F)
    if (denom==0){return(1)}else{return(num/denom)}
  }
}

selective_ci = function(xx, SSigma, threshold, alpha = 0.1) {
  # Confidence interval for mu[1], valid conditional on |x[2]| > threshold
  base.factor = 10
  low = tryCatch({
    uniroot(
      f = function(mu) { selective_pval(mu, xx, SSigma, threshold) - alpha/2 },
      lower = as.numeric(xx[1] - 10*sqrt(SSigma[1,1])),
      upper = as.numeric(xx[1] + 10*sqrt(SSigma[1,1]))
    )$root
  }, error = function(e) {NA}
  )
  while (is.na(low) & base.factor < 2^8){
    base.factor = base.factor * 2
    low = tryCatch({
      uniroot(
        f = function(mu) { selective_pval(mu, xx, SSigma, threshold) - alpha/2 },
        lower = as.numeric(xx[1] - base.factor*sqrt(SSigma[1,1])),
        upper = as.numeric(xx[1] + base.factor*sqrt(SSigma[1,1]))
      )$root
    }, error = function(e) {NA}
    )
  }
  base.factor = 10
  high = tryCatch({
    uniroot(
      f = function(mu) { selective_pval(mu, xx, SSigma, threshold) - (1-alpha/2) },
      lower = as.numeric(xx[1] - 10*sqrt(SSigma[1,1])),
      upper = as.numeric(xx[1] + 10*sqrt(SSigma[1,1]))
    )$root
  }, error = function(e) {NA}
  )
  while (is.na(high) & base.factor < 2^8){
    base.factor = base.factor * 2
    high = tryCatch({
      uniroot(
        f = function(mu) { selective_pval(mu, xx, SSigma, threshold) - (1-alpha/2) },
        lower = as.numeric(xx[1] - 10*sqrt(SSigma[1,1])),
        upper = as.numeric(xx[1] + 10*sqrt(SSigma[1,1]))
      )$root
    }, error = function(e) {NA}
    )
  }
  if (!is.na(low) & !is.na(high)){
    return(data.frame("low" = min(low, high), "high" = max(low, high)))
  }else{
    return(data.frame("low"=-Inf, "high"=Inf))
  }
  
}





clean_vars <- function(data, turn_numerics, othername=FALSE){
  data_all = data %>% mutate(age = clean_age(yearbirth)) %>% filter(age < 100)
  if (othername){
    data_all$pincome = suppressWarnings(as.numeric(gsub(">", "", gsub("\\+", "", gsub(",", "", gsub("\\$", "", data_all$familyinc)))))) / 10^4 
    # data_all$pincome[is.na(data_all$pincome)] = median(data_all$pincome, na.rm=TRUE)
  }else{
    data_all$pincome = suppressWarnings(as.numeric(gsub(">", "", gsub("\\+", "", gsub(",", "", gsub("\\$", "", data_all$faminc)))))) / 10^4
    # data_all$pincome[is.na(data_all$pincome)] = median(data_all$pincome, na.rm=TRUE)
  }
  # clean up extreme values of pincome 
  data_all$pincome = pmin(data_all$pincome, 1000)
  
  data_all = data_all %>% mutate(gender = gender - 1, pincome2 = pincome^2, parented2 = parented^2, age2 = age^2)
  
  data_all[turn_numerics] = sapply(data_all[turn_numerics], as.numeric)
  
  return(data_all)
}




between_study_regression <- function(data, all_sites, all_grid, covariates, if_center = FALSE, if_scale = FALSE){
  all_res = matrix(0, nrow=0, ncol=length(covariates))
  # data[,covariates] = data.frame(data.frame(data[,covariates]) %>% scale(center = TRUE, scale = TRUE))
  data[,covariates] = whiten(as.matrix(data[,covariates]))
  for (ii in 1:(nrow(all_grid)/2)){
    # covariates = c("pltclideo", "gender", "age", "parented")
    data_org = data %>% filter(datacollection == all_sites[all_grid[ii,1]]) %>% 
      mutate(across(covariates, ~replace_na(., mean(., na.rm=TRUE))))  
    data_rep = data %>% filter(datacollection == all_sites[all_grid[ii,2]]) %>% 
      mutate(across(covariates, ~replace_na(., mean(., na.rm=TRUE))))   
    
    # analysis   
    data_org = data.frame(data_org %>% select(covariates) %>% scale(center = if_center, scale = if_scale))
    data_rep = data.frame(data_rep %>% select(covariates) %>% scale(center = if_center, scale = if_scale))
    glm.data = rbind(data_org %>% mutate(study = 1),
                     data_rep %>% mutate(study = -1))
    # glm.mdl = tryCatch({glm(study ~., family = 'binomial', data = glm.data)}, 
    #                    error = function(e){NULL})
    glm.mdl = tryCatch({lm(study ~., data = glm.data)}, 
                       error = function(e){NULL})
    if (!is.null(glm.mdl)){
      glm.coef = as.numeric(coef(glm.mdl)[covariates])
      all_res = rbind(all_res, glm.coef)
    }
    
  }
  return(all_res)
}

pooled_study_regression <- function(data, all_sites, covariates, if_center = FALSE, if_scale = FALSE){
  pooled_res = matrix(0, nrow = 0, ncol = length(covariates))
  for (ii in 1:length(all_sites)){
    # covariates = c("pltclideo", "gender", "age", "parented")
    data_single = data %>% filter(datacollection == all_sites[ii]) %>% 
      mutate(across(covariates, ~replace_na(., mean(., na.rm=TRUE))))  
    data_pooled = data %>% 
      mutate(across(covariates, ~replace_na(., mean(., na.rm=TRUE))))   
    
    # analysis   
    data_single = data.frame(data_single %>% select(covariates) %>% scale(center = if_center, scale = if_scale))
    data_pooled = data.frame(data_pooled %>% select(covariates) %>% scale(center = if_center, scale = if_scale))
    glm.data = rbind(data_single %>% mutate(study = 1),
                     data_pooled %>% mutate(study = 0))
    glm.mdl = tryCatch({glm(study ~., family = 'binomial', data = glm.data)}, 
                       error = function(e){NULL})
    if (!is.null(glm.mdl)){
      glm.coef = as.numeric(coef(glm.mdl)[covariates])
      pooled_res = rbind(pooled_res, glm.coef)
    }
  }
  return(pooled_res)
}


clean_age <- function(yearbirth){ 
  age = yearbirth
  age[yearbirth>=1000 & yearbirth <= 2014 & !is.na(yearbirth)] = 2014 - yearbirth[yearbirth>=1000 & yearbirth <= 2014 & !is.na(yearbirth)]
  age[age >= 80 & age <= 99 & !is.na(age)] = 114 - age[age >= 80 & age <= 99 & !is.na(age)]
  return(age)
} 

compute_covariate <- function(df1, df2, covariates){
  ratio = c()
  for (name in covariates){
    ratio = c(ratio, (mean(df1[,name], na.rm=T) - mean(df2[,name], na.rm=T))/ (sd(df1[,name], na.rm=T) ))
  }
  return(ratio)
}

compute_covariate_normalize <- function(df1, df2, covariates){
  ratio = c()
  for (name in covariates){
    ratio = c(ratio, 
              (mean(df1[,name], na.rm=T) - mean(df2[,name], na.rm=T))/ 
                (sqrt(sd(df1[,name], na.rm=T)^2 / sum(!is.na(df1[,name])) + sd(df2[,name], na.rm=T)^2 / sum(!is.na(df2[,name])) )))
  }
  return(ratio)
}



 

