// Input data : number of genes (G),  Number of clusters (K),
// Expression vector (X) and the list of probability vectors (Q)

data {
  int G;
  int K;
  matrix<lower=0>[G, K] Q; 
  int X[G];
}

// The parameters accepted by the model. 
//Our model accepts one parameter 
parameters {
  simplex[K] alpha;
}

// The transformed parameters : the final probability vector

transformed parameters {
  vector[G] p;
  p = Q*alpha ;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  X ~ multinomial(p);
}

