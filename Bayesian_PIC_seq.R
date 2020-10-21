library(Ternary)
library(rstan)
library(Matrix)
library(parallel)
library(foreach)
library(doParallel)

Core_demultiplexing_function = function(Q,G,K,X) {
  
  Data = list(Q=Q,G=G,K=K,X=X)
  
  Rstan_fit = stan(refresh=0, #Make Rstan silent
               file = "PIC_seq_Stan_script.stan",  #Path to the Rstan script
               data = Data,    
               chains = 4,            
               warmup = 50,        
               iter = 250,          
               cores = 1,           
  )
  Sumamry_fit = summary(Rstan_fit)
  
  #If the fitting worked : extracting the median value of the sampled parameter
  if (!is.null(Sumamry_fit)) {
    p = Sumamry_fit$summary
    Alpha_vector = p[grepl(rownames(p),pattern = "alpha",fixed = T),5]
  }
  
  if (is.null(Sumamry_fit)) {
    Alpha_vector = rep(NA,K)
  }
  return(Alpha_vector)
} 

PIC_demultiplexing = function(PIC_Expression_matrix,Singlet_mean_profile,N_cores = 6) {
  
  #Managing parallelization
  N_detected_cores = detectCores()
  
  if (N_cores>=N_detected_cores) {
    warning("Not enough cores are available on this computer. Reducing the number of used threads !")
    N_cores = N_detected_cores - 1
  }
  
  registerDoParallel(N_cores)
  
  #Checking the data
  if (nrow(PIC_Expression_matrix)!=nrow(Singlet_mean_profile)) {
    stop("The number of genes in the PIC and singlet expression matrix are not identical. Please check your input data !")
  }
  
  #Defining variables
  G = nrow(PIC_Expression_matrix) # Number of genes measured
  K = ncol(Singlet_mean_profile) # Number of singlet clusters
  N_PIC = ncol(PIC_Expression_matrix)
  
  #Demutiplexing per se : done in parallel using a foreach loop
  #Alpha table : rows correspond to the contribution of each cell type and columns to PIC/doublets
  
  cat("Starting the fitting of individual doublets. \n")
  
  Alpha_table =foreach (i=1:N_PIC, .combine=cbind) %dopar% {
    X = PIC_Expression_matrix[,i]
    Alpha_vector = Core_demultiplexing_function(Q = Singlet_mean_profile,G = G,K = K,X = X)
  }
  
  cat("Fitting of the doublets has been performed. \n")
  
  #Renaming row and columns of the Alpha table 
  colnames(Alpha_table) = colnames(PIC_Expression_matrix)
  rownames(Alpha_table) = colnames(Singlet_mean_profile)
  
  #Counting then number of cells for which the fitting did not worked
  N_failure_fitting = colSums(is.na(Alpha_table))
  N_failure_fitting = sum(N_failure_fitting!=0)
  
  cat(paste("Out of ",N_PIC, " PICs, ", N_PIC-N_failure_fitting," were successfully fitted ! \n",sep=""))

  return(Alpha_table)
}

Quality_control_PIC_demultiplexing = function(Alpha_table,PIC_Expression_matrix,Singlet_mean_profile,Non_PIC_clusters=NULL) {
  
  N_PIC = ncol(PIC_Expression_matrix)
  L = colSums(PIC_Expression_matrix)
  #Checking the input data 
  
  if (nrow(Alpha_table)!=ncol(Singlet_mean_profile)){
    stop("Number of single-cell clusters is not coherent between the alpha table and the single-cell mean profiles")
  }
  
  if (nrow(PIC_Expression_matrix)!=nrow(Singlet_mean_profile)) {
    stop("The number of genes in the PIC and singlet expression matrix are not identical. Please check your input data !")
  }
  
  #Computing the expected PIC profile and the correlation between expected and observed mean profiles
  
  Expected_profile = Singlet_mean_profile%*%Alpha_table
  Observed_profile = t(t(PIC_Expression_matrix)/colSums(PIC_Expression_matrix))
  
  List_correlation = c()
  
  for (k in 1:N_PIC){
    
    x = asinh(Expected_profile[,k]*10^3)
    y = asinh(Observed_profile[,k]*10^3)
    z = !is.infinite(x) & !is.infinite(y)
    List_correlation = c(List_correlation,cor(x[z],y[z],use = "pairwise.complete"))
    
  }
  names(List_correlation) = colnames(PIC_Expression_matrix)
  
  #Main covariate explaining the quality of the fitting : total number of UMIs 
  
  par(las=1,bty="l")
  plot(L,List_correlation,log="x",pch=21,bg="orange",
       xlab="Number of PIC UMIs used for fitting",ylim=c(0,1),yaxs='i',
       ylab="Correlation between observed and predicted profile")
  m = lm(List_correlation~log10(L))
  abline(coef(m),lwd=2,lty=2,col="red")
  Corrected_quality = m$residuals
  Correlation = summary(m)$r.squared
  legend("topleft",legend =  paste("R=",round(sqrt(Correlation),digits = 2),sep = ""),bty="n",cex = 1.4)
  
  Corrected_quality = Corrected_quality[colnames(PIC_Expression_matrix)]
  names(Corrected_quality) = colnames(PIC_Expression_matrix)
  
  if (!is.null(Non_PIC_clusters)) {
    Non_PIC_proportion = colSums(Alpha_table[Non_PIC_clusters,],na.rm = T)*100
    par(las=1,bty="l")
    hist(Non_PIC_proportion,xaxs="i",yaxs="i",n=40,col="grey70",main=NULL,xlab="Contribution of un-expected cell in the PIC (%)",cex.lab=1.3)
    abline(v=median(Non_PIC_proportion,na.rm = T),lwd=2,lty=2,col="red")
    CI_95 = quantile(Non_PIC_proportion,na.rm = T,probs = c(0.05,0.95))
    abline(v=CI_95,lwd=1,lty=2,col="black")
    legend("topright",legend = paste("Median contribution : ",round(median(Non_PIC_proportion,na.rm = T),2),"% \n 95% IQR : ",round(CI_95[1],2),"% -",round(CI_95[2],2),"%"),bty="n",cex=1.3)
  }

  
  
  return(Corrected_quality)
}

##Visualization of the fitting

Ternary_QC_plot = function(Alpha_table,Classification_clusters,Experimental_condition=NULL) {
  if (length(Classification_clusters)!=3){
    stop("The clusters have not been devided into three as required, please check the Classification_clusters object (must be a list of length 3)")
  }
  Ternary_tables =  cbind(colSums(Alpha_table[Classification_clusters[[1]],],),
                          colSums(Alpha_table[Classification_clusters[[2]],],),
                          colSums(Alpha_table[Classification_clusters[[3]],],))
  colnames(Ternary_tables) = names(Classification_clusters)
  
  if (is.null(Experimental_condition)) {
    TernaryPlot(colnames(Ternary_tables)[1],colnames(Ternary_tables)[2],colnames(Ternary_tables)[3],axis.cex = 0.5)
    TernaryPoints(Ternary_tables,pch=21,bg="orange")
    TernaryDensityContour((Ternary_tables))
  }
  
  if (!is.null(Experimental_condition)) {
    
    n = length(unique(Experimental_condition))
    order_condition =unique(as.character(Experimental_condition))
    
    if (is.factor(Experimental_condition)){
      order_condition = levels(Experimental_condition)
    }
    
    par(mfrow=c(1,n))
    
    for (k in order_condition){
      TernaryPlot(colnames(Ternary_tables)[1],colnames(Ternary_tables)[2],colnames(Ternary_tables)[3],axis.cex = 0.5)
      title(main = paste(k," (n=",sum(Experimental_condition==k)," PICs)",sep = ""))
      TernaryPoints(Ternary_tables[Experimental_condition==k,],pch=21,bg="orange")
      TernaryDensityContour((Ternary_tables[Experimental_condition==k,]))
      
    }
    
  }
  
}

##How to detect genes specifically expressed in PIC-seq

gene = "Foxp3"
X = Data_to_fit[gene,]
p = Q%*%Alpha_table
L = colSums(Data_to_fit)

NB_log_likelihood = function(X,Expected_lambda){
  mean_gene_expression = mean(X)
  var_gene_expres = var(X)
  phi = mean_gene_expression^2/(var_gene_expres-mean_gene_expression)
  LL = sum(lgamma(X+phi)-lgamma(X+1)-lgamma(phi)+X*log(mean_gene_expression/(mean_gene_expression+phi))+phi*log(phi/(mean_gene_expression+phi)))
  return(LL)
}


Compute_gene_divergence = function(Alpha_table,PIC_Expression_matrix,Singlet_mean_profile) {
  Expected_expression_profile = Singlet_mean_profile%*%Alpha_table
  L = colSums(PIC_Expression_matrix)
  Gene_LL = c()
  
  for (k in rownames(PIC_Expression_matrix)) {
    LL = NB_log_likelihood(PIC_Expression_matrix[k,],L*Expected_expression_profile[k,])
    
    Gene_LL = c(Gene_LL,LL)
  }
  names(Gene_LL) = rownames(PIC_Expression_matrix)
}




Likelihood_ratio_test = function(X,p,L) {
  function_to_optimize = function(alpha) {
    S = -sum(-L*(p+alpha)+X*log(L*(p+alpha))-lgamma(X+1))
    return(S)
  }
  
  
  alpha_fitting = optimize(function_to_optimize,interval = c(-1,1))
  alpha = alpha_fitting$minimum
  p_prime = p+alpha
  LR = -2*log(sum(-L*p+X*log(L*p)-lgamma(X+1))/sum(-L*p_prime+X*log(L*p_prime)-lgamma(X+1)))
  P_value = pchisq(q = LR,df = 1)
}
