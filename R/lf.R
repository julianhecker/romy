

#' ROMY-LF
#'
#' This function performs the conditional latent factor model computations
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param inds_X List of indices for each column of Y that links covariates in X to the specific Y_m component.
#' @param epsilon Convergence threshold, default value is 0.01.
#' @param K_min Minimum number of latent factors, default value is 1.
#' @param K_max Maximum number of latent factors, default value is 5.
#' @param max_iter Maximum number of iterations, default value is 20.
#' @param eta Eta parameter, default is value is 0.5.
#'
#' @export
romy_lf=function(Y, X, inds_X, epsilon=0.01, K_min=1, K_max=5, max_iter=20, eta=0.5)
{
  
  
  results=list()
  ############################################################################################
  c_lambda=1
  c_delta=0.2
  convergence=FALSE
  iter=1
  while(!convergence & iter<=5)
  {
    res=.romy_lf_algorithm(Y=Y, X=X, inds_X=inds_X, epsilon=epsilon, K_min=K_min, K_max=K_max, 
                          c_lambda=c_lambda, c_delta=c_delta, max_iter=max_iter, eta=eta)
    convergence=res$convergence
    iter=iter+1
  }
  tmp_result=list(convergence=convergence, K_est=res$K_est, K_hat=res$K_hat,
                  Z_hat=as.matrix(res$Z_hat), Gamma_hat=res$Gamma_hat)
  results[[1]]=tmp_result
  ############################################################################################
  c_lambda=2
  c_delta=0.2
  convergence=FALSE
  iter=1
  while(!convergence & iter<=5)
  {
    res=.romy_lf_algorithm(Y=Y, X=X, inds_X=inds_X, epsilon=epsilon, K_min=K_min, K_max=K_max, 
                          c_lambda=c_lambda, c_delta=c_delta, max_iter=max_iter, eta=eta)
    convergence=res$convergence
    iter=iter+1
  }
  tmp_result=list(convergence=convergence, K_est=res$K_est, K_hat=res$K_hat,
                  Z_hat=as.matrix(res$Z_hat), Gamma_hat=res$Gamma_hat)
  results[[2]]=tmp_result
  ############################################################################################
  c_lambda=2
  c_delta=0.4
  convergence=FALSE
  iter=1
  while(!convergence & iter<=5)
  {
    res=.romy_lf_algorithm(Y=Y, X=X, inds_X=inds_X, epsilon=epsilon, K_min=K_min, K_max=K_max, 
                          c_lambda=c_lambda, c_delta=c_delta, max_iter=max_iter, eta=eta)
    convergence=res$convergence
    iter=iter+1
  }
  tmp_result=list(convergence=convergence, K_est=res$K_est, K_hat=res$K_hat,
                  Z_hat=as.matrix(res$Z_hat), Gamma_hat=res$Gamma_hat)
  results[[3]]=tmp_result
  
  ############################################################################################
  convergence=FALSE
  Z_hat=as.matrix(0)
  Gamma_hat=as.matrix(0)
  K_hat=0
  K_est=0
  
  if(results[[3]]$convergence)
  {
    convergence=results[[3]]$convergence
    Z_hat=as.matrix(results[[3]]$Z_hat)
    Gamma_hat=as.matrix(results[[3]]$Gamma_hat)
    K_hat=results[[3]]$K_hat
	K_est=results[[3]]$K_est
  }
  if(results[[2]]$convergence & results[[2]]$K_est<K_max)
  {
    convergence=results[[2]]$convergence
    Z_hat=as.matrix(results[[2]]$Z_hat)
    Gamma_hat=as.matrix(results[[2]]$Gamma_hat)
    K_hat=results[[2]]$K_hat
	K_est=results[[2]]$K_est
  }
  if(results[[1]]$convergence & results[[1]]$K_est<K_max)
  {
    convergence=results[[1]]$convergence
    Z_hat=as.matrix(results[[1]]$Z_hat)
    Gamma_hat=as.matrix(results[[1]]$Gamma_hat)
    K_hat=results[[1]]$K_hat
	K_est=results[[1]]$K_est
  }
  
  return(list(convergence=convergence, Z_hat=Z_hat, K_est=K_est, K_hat=K_hat, Gamma_hat=Gamma_hat))
  
}

#' Performs the proximal gradient descent algorithm as described by Chen 2022 (Appendix B).
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param inds_X List of indices for each column of Y that links covariates in X to the specific Y_m component.
#' @param epsilon Convergence threshold, default value is 0.01.
#' @param K_min Minimum number of latent factors, default value is 1.
#' @param K_max Maximum number of latent factors, default value is 5.
#' @param c_lambda factor for Lambda.
#' @param c_delta factor for Delta.
#' @param max_iter Maximum number of iterations, default value is 20.
#' @param eta Eta parameter, default is value is 0.5.
#'
.romy_lf_algorithm=function(Y, X, inds_X, epsilon, K_min, K_max, c_lambda, c_delta, max_iter=20, eta=0.5)
{
  N=nrow(Y)
  M=ncol(Y)
  P=sum(lengths(inds_X))
  
  
  w_km1=1
  w_k=1
  
  Gamma_s_km1=matrix(rnorm(N*P), nrow=N, ncol=P)
  Gamma_s_k=matrix(rnorm(N*P), nrow=N, ncol=P)
  
  #################
  L=-1
  for(m in 1:M)
  {
    tmp_X=as.matrix(X[,inds_X[[m]]])
    tmp_X_sqsum=rowSums(tmp_X**2)
    if(L<max(tmp_X_sqsum))
    {
      L=max(tmp_X_sqsum)
    }
  }
  #################
  
  lambda_NP=c_lambda*sqrt((N+P)*log(N))
  delta_NP=c_delta*(N+P)*log(N)
  
  
  tau_0=L
  tau_km1=tau_0
  
  
  k=1
  
  convergence=FALSE
  
  while(!convergence & k<=max_iter){
		# step 1
		
		Gamma_k = Gamma_s_k + (w_km1-1)/w_k *(Gamma_s_k - Gamma_s_km1)
		
		# step 2
		tau_hat_j = eta * tau_km1
		run=TRUE
		iter_in=0
		while(run & iter_in<=max_iter)
		{
		  iter_in=iter_in+1
		  #print(tau_hat_j)
		  x = 1/tau_hat_j * lambda_NP
		  A = Gamma_k - 1/tau_hat_j*.gradient_f_func(Y, X, inds_X, Gamma_k)
		  A_j = .S_operator(A, x)
		  Ftmp=.cF_func(Y, X, inds_X, A_j)
		  Qtmp=.Q_func(Y, X, inds_X, A_j, Gamma_k, tau_hat_j, lambda_NP)
		  if(Ftmp<=Qtmp)
		  {
			tau_k = tau_hat_j
			run = FALSE
			break
		  }
		  else{
			tau_hat_j = min(1/eta*tau_hat_j, tau_0)
		  }
		}
		tau_k=tau_0
		# step 3
		x = 1/tau_k * lambda_NP
		A = Gamma_k - 1/tau_k * .gradient_f_func(Y, X, inds_X, Gamma_k)
		
		Gamma_s_kp1 = .S_operator(A, x)
		
		
		# step 4
		w_kp1 = (1 + sqrt(1 + 4*w_k**2))/2
		
		# step 5
		
		D_kp1= tau_k *(Gamma_k - Gamma_s_kp1) + .gradient_f_func(Y, X, inds_X, Gamma_s_kp1) - .gradient_f_func(Y, X, inds_X, Gamma_k)
		norm_D_kp1=sum(D_kp1**2)
		factor_value=tau_k*max(1, sum(Gamma_s_kp1**2))
		#cat("norm/factor:", norm_D_kp1/factor, "\n")
		if(norm_D_kp1/factor_value<=epsilon)
		{
		  Gamma_out=Gamma_s_kp1
		  convergence=TRUE
		  break
		}
		else{
		  w_km1=w_k
		  w_k=w_kp1
		  tau_km1=tau_k
		  Gamma_s_km1=Gamma_s_k
		  Gamma_s_k=Gamma_s_kp1
		  k=k+1
		}
    
  }
  
  if(convergence)
  {
		Gamma_hat=Gamma_out
		
		
		M_N=(diag(N)-1/N*matrix(1, nrow=N, ncol=N))
		
		M_Gamma=M_N %*%Gamma_hat
		
		M_GammaGamma_M=M_Gamma %*% t(M_Gamma)
		
		evd=eigen(M_GammaGamma_M)
		inds=1:N
		K_est=0
		if(sum(evd$values>=delta_NP)>=1)
		{
		  K_est=max(inds[evd$values>=delta_NP])
		}
		K_hat=K_est
		if(K_hat<=K_min){K_hat=K_min}
		if(K_hat>=K_max){K_hat=K_max}
		Z_hat=evd$vectors[,1:K_hat]*sqrt(N)
  }else{
		K_hat=0
		Z_hat=0
		K_est=0
		Gamma_hat=0
  }
  
  return(list(convergence=convergence, Z_hat=as.matrix(Z_hat), K_hat=K_hat, Gamma_hat=as.matrix(Gamma_hat), K_est=K_est))
}

#' Computes nuclear norm of Matrix Gamma
#'
#' @param Gamma Matrix
#'
.nuclear_norm=function(Gamma)
{
    svd_tmp=svd(Gamma)
    return(sum(svd_tmp$d))
}

#' f function as described in Appendix B of Chen 2022.
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param inds_X List of indices for each column of Y that links covariates in X to the specific Y_m component.
#' @param Gamma Matrix Gamma
#'
.f_func=function(Y, X, inds_X, Gamma)
{
	M=ncol(Y)
    Yp=matrix(0, nrow=nrow(Y), ncol=ncol(Y))
    for(m in 1:M)
    {
      Yp[,m]=rowSums(as.matrix(X[,inds_X[[m]]]) * as.matrix(Gamma[,inds_X[[m]]]))
    }
    diff=Y-Yp
    f=sum(diff**2)
    return(f)
}

#' gradient f function as described in Appendix B of Chen 2022.
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param inds_X List of indices for each column of Y that links covariates in X to the specific Y_m component.
#' @param Gamma Matrix Gamma
#'
.gradient_f_func=function(Y, X, inds_X, Gamma)
{
    N=nrow(Y)
	M=ncol(Y)
    P=ncol(Gamma)
    grad=matrix(0, nrow=N, ncol=P)
    for(m in 1:M)
    {
      pred=rowSums(as.matrix(X[,inds_X[[m]]]) * as.matrix(Gamma[,inds_X[[m]]]))
      grad[, inds_X[[m]]]=-(Y[,m]-pred)*X[,inds_X[[m]]]
    }
    return(grad)
}

#' F function as described in Appendix B of Chen 2022.
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param inds_X List of indices for each column of Y that links covariates in X to the specific Y_m component.
#' @param Gamma Matrix Gamma
#'
.cF_func=function(Y, X, inds_X, Gamma)
{
    f=.f_func(Y, X, inds_X, Gamma)
    nn=.nuclear_norm(Gamma)
    return(f+nn)
}

#' Q function as described in Appendix B of Chen 2022.
#'
#' @param Y A matrix containing the measurements for Y.
#' @param X A matrix containing the measurements for X.
#' @param inds_X List of indices for each column of Y that links covariates in X to the specific Y_m component.
#' @param Gamma Matrix Gamma
#' @param Gamma_k Matrix Gamma_k
#' @param tau_k Parameter tau_k
#' @param lambda Parameter lambda
#'
.Q_func=function(Y, X, inds_X, Gamma, Gamma_k, tau_k, lambda)
{
    f=.f_func(Y, X, inds_X, Gamma_k)
    grad=.gradient_f_func(Y, X, inds_X, Gamma_k)
    trace= sum((Gamma-Gamma_k) * grad)
    diff=tau_k/2 * sum((Gamma-Gamma_k)**2)
    penalty=lambda*.nuclear_norm(Gamma)
    return(f+trace+diff+penalty)
  
}
#' S operator as described in Appendix B of Chen 2022.
#'
#' @param A Matrix A.
#' @param x Threshold x.
#'
.S_operator=function(A, x)
{
    svd_tmp=svd(A)
    diag_mat=diag(ifelse(svd_tmp$d-x>0, svd_tmp$d-x,0))
    A_new=svd_tmp$u %*% diag_mat %*% t(svd_tmp$v)
    return(A_new)
}



