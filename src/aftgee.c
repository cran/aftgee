#include <R.h>
#include <Rmath.h>
#include <math.h>

/* uFun    : smoothed Gehan */
/* unsFun  : non smoothed Gehan */
/* ulogFun : smoothed log rank */

double get_rikjl(double *X, double *sigma, int *clsize, int *p, int *N, int *n, int ik_idx, int jl_idx) {
  double *xdif = Calloc(*p, double);
  double rikjl = 0.0;
  int m = 0, q = 0;

 for (m = 0; m < *p; m++) {
    xdif[m] = 0.0;
    xdif[m] = X[ik_idx + m * *N] - X[jl_idx + m * *N];
  }

  for (m = 0; m < *p; m++) {
    for (q = 0; q < *p; q++) {
      rikjl += xdif[m] * sigma[m * *p + q] * xdif[q];
    }
  }

  rikjl = sqrt(rikjl);
  Free(xdif);
  return(rikjl);
}

void abargehanfun(double *beta, double *Y, double *X, double *delta, int *clsize,
	     double *sigma,
	     int *n,
	     int *p,
	     int *N,
	     double *weight,
	     //output
	     double *abar) {
  int i, j, k, l, ik_idx = 0, jl_idx, r, s, xdif_idx, abar_idx;
  double *e = Calloc(*N, double), *xdif = Calloc(*p, double), *nu1 = Calloc(*p, double), *nu2 = Calloc(*p, double), *nu3 = Calloc(*p, double);
  double rikjl, edif, H, h, u, z, sqrtn = sqrt(*n), coef, de, coef1;

  /* compute the e vector */
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
  }

  /*The main abar part */
  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
      if (delta[ik_idx] != 0) {
	/* reset */
	for ( r = 0; r < *p; r++) {
	  nu1[r] = 0.0;
	  nu2[r] = 0.0;
	  nu3[r] = 0.0;
	}
	de = 0.0;
	coef1 = 0.0;
	jl_idx = 0;
	for (j = 0; j < *n; j++) {
	  for (l = 0; l < clsize[j]; l++) {
	    rikjl = get_rikjl(X, sigma, clsize, p, N, n, ik_idx, jl_idx);
	    if (rikjl != 0) {
	      edif = e[jl_idx] - e[ik_idx];
	      z = sqrtn * edif / rikjl;
	      h = dnorm(z, 0.0, 1.0, 0);
	      coef = weight[ik_idx] * weight[jl_idx] * h * sqrtn / rikjl;
	      for (r = 0; r < *p; r++) {
		xdif[r] = 0.0;
		xdif[r] = X[ik_idx + r * *N] - X[jl_idx + r * *N];
	      }
	      xdif_idx = 0;
	      for (r = 0; r < *p; r++) {
		for (s = 0; s < *p; s++) {
		  abar[xdif_idx] += coef * xdif[r] * xdif[s];  // p by p matrix
		  xdif_idx++;
		} // end for s
	      } // end for r
	    } // end if rikjl != 0
	    jl_idx++;
	  } // end for l
	} // end for j
      } // end if delta[ik_idx] != 0
      ik_idx++;
    } // end for k
  } // end for i

  Free(xdif);
  Free(e);
  Free(nu1);
  Free(nu2);
  Free(nu3);
  /* for (r = 0; r < (*p * *p); r++) { */
  /*   abar[r] /= (*n * (*n - 1)); */
  /* } */
  abar;
}


void abarlogfun(double *beta, double *Y, double *X, double *delta, int *clsize,
	     double *sigma,
	     int *n,
	     int *p,
	     int *N,
	     double *weight,
		/* output */
	     double *abar) {
  int i, j, k, l, ik_idx = 0, jl_idx, r, s, xdif_idx, xdif_idx2, abar_idx;
  double *e = Calloc(*N, double), *xdif = Calloc(*p, double), *nu1 = Calloc(*p * *p, double), *nu2 = Calloc(*p, double), *nu3 = Calloc(*p, double);
  double rikjl, edif, H, h, u, z, sqrtn = sqrt(*n), coef, de, coef1;

  /* compute the e vector */
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
  }

  /*The main abar part */
  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
      if (delta[ik_idx] != 0) {
	/* reset */
	xdif_idx = 0;
	for ( r = 0; r < *p; r++) {
	  for ( s = 0; s < *p; s++) {
	    nu1[xdif_idx] = 0.0;
	    xdif_idx++;
	  }
	  nu2[r] = 0.0;
	  nu3[r] = 0.0;
	}
	de = 0.0;
	coef = 0.0;
	jl_idx = 0;
	for (j = 0; j < *n; j++) {
	  for (l = 0; l < clsize[j]; l++) {
	    rikjl = get_rikjl(X, sigma, clsize, p, N, n, ik_idx, jl_idx);
	    if (rikjl != 0) {
	      edif = e[jl_idx] - e[ik_idx];
	      z = sqrtn * edif / rikjl;
	      h = dnorm(z, 0.0, 1.0, 0);
	      H = pnorm(z, 0.0, 1.0, 1, 0);
	      coef = weight[ik_idx] * weight[jl_idx] * h * sqrtn / rikjl;
	      for (r = 0; r < *p; r++) {
		xdif[r] = 0.0;
		xdif[r] = X[ik_idx + r * *N] - X[jl_idx + r * *N];
	      }
	      xdif_idx = 0;
	      for (r = 0; r < *p; r++) {
		for (s = 0; s < *p; s++) {
		  nu1[xdif_idx] += coef * X[jl_idx + s * *N] * xdif[r];  // s, r
		  xdif_idx++;
		}
		nu2[r] += weight[jl_idx] * H * X[jl_idx + r * *N];
		nu3[r] += coef * xdif[r];
	      }
	      de += weight[jl_idx] * H;
	    } // end if rikjl != 0
	    jl_idx++;
	  } // end for l
	} // end for j
	xdif_idx2 = 0;
	for (r = 0; r < *p; r++) {
	  for (s = 0; s < *p; s++) {
	    abar[xdif_idx2] += -1 * (de * nu1[xdif_idx2] - nu2[s] * nu3[r]) / (de * de); //
	    xdif_idx2++;
	  }
	}
      } // end if delta[ik_idx] != 0
      ik_idx++;
    } // end for k
  } // end for i

  Free(xdif);
  Free(e);
  Free(nu1);
  Free(nu2);
  Free(nu3);
  /* for (r = 0; r < (*p * *p); r++) { */
  /*   abar[r] /= (*n * (*n - 1)); */
  /* } */
  abar;
}


void ufun(double *beta,
	  double *Y, double *X, double *delta, int *clsize,
	  double *sigma,
	  int *n, // numbers of clusters
	  int *p, // ncol(X), or dimension
	  int *N, // nrow(X), or numbers of dataset. Also equals to sum(clsize)
	  //output
	  double *Z,
	  double *weight,
	  double *sn) {
  int i, j, k, l, ik_idx = 0, jl_idx, m, r;
  double *e = Calloc(*N, double), *xdif = Calloc(*p, double);
  double rikjl, z, H, sqrtn = sqrt(*n), edif;

  /* compute the e vector */
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
  }

  /* e is in here, cant compute it else where*/
  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
      if (delta[ik_idx] != 0) {
	jl_idx = 0;
	for (j = 0; j < *n; j++) {
	  for (l = 0; l < clsize[j]; l++) {
	    // if (e[ik_idx] - e[jl_idx] <= 0) {
	    rikjl = get_rikjl(X, sigma, clsize, p, N, n, ik_idx, jl_idx);
	    if (rikjl != 0) {
	      edif = e[jl_idx] - e[ik_idx];
	      z = sqrtn * edif / rikjl;
	      H = pnorm(z, 0.0, 1.0, 1, 0);
	      for (r = 0; r < *p; r++) {
		xdif[r] = 0.0;
		xdif[r] = X[ik_idx + r * *N] - X[jl_idx + r * *N];
		sn[r] += weight[ik_idx] * weight[jl_idx] * Z[ik_idx] * Z[jl_idx] * xdif[r] * H;
	      } // end for r
	    }
	    jl_idx++;
	  }
	}
      }
      ik_idx++;
    }
  }
  //  Free(edif);
  Free(xdif);
  Free(e);
  /* for (r = 0; r < *p; r++) { */
  /*   sn[r] /= (*n * (*n - 1)); */
  /* } */
  sn;
}


void unsfun(double *beta,
	    double *Y, double *X, double *delta, int *clsize,
	    double *sigma,
	    int *n, // numbers of clusters
	    int *p, // ncol(X), or dimension
	    int *N, // nrow(X), or numbers of dataset. Also equals to sum(clsize)
	    //output
	    double *Z,
	    double *weight,
	    double *sn) {
  int i, j, k, l, ik_idx = 0, jl_idx, m, r;
  double *e = Calloc(*N, double), *xdif = Calloc(*p, double);
  double rikjl, z, H, sqrtn = sqrt(*n), edif;

  /* compute the e vector */
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
  }

  /* e is in here, cant compute it else where*/
  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
      if (delta[ik_idx] != 0) {
	jl_idx = 0;
	for (j = 0; j < *n; j++) {
	  for (l = 0; l < clsize[j]; l++) {
	    if (e[ik_idx] - e[jl_idx] <= 0) {
	      for (r = 0; r < *p; r++) {
		xdif[r] = 0.0;
		xdif[r] = X[ik_idx + r * *N] - X[jl_idx + r * *N];
		sn[r] += weight[ik_idx] * weight[jl_idx] * Z[ik_idx] * Z[jl_idx] * xdif[r];
	      } // end e[ik] - e[jl]
	    }
	    jl_idx++;
	  }
	}
      }
      ik_idx++;
    }
  }
  //  Free(edif);
  Free(xdif);
  Free(e);
  /* for (r = 0; r < *p; r++) { */
  /*   sn[r] /= (*n * (*n - 1)); */
  /* } */
  sn;
}


void ulogfun(double *beta,
	     double *Y, double *X, double *delta, int *clsize,
	     double *sigma,
	     int *n, // numbers of clusters
	     int *p, // ncol(X), or dimension
	     int *N, // nrow(X), or numbers of dataset. Also equals to sum(clsize)
	     //output
	     double *Z,
	     double *weight,
	     double *pw,
	     double *sn) {
  int i, j, k, l, ik_idx = 0, jl_idx, m, r;
  double *e = Calloc(*N, double), *nu = Calloc(*p, double);
  double rikjl, z, H, sqrtn = sqrt(*n), edif, de;

  /* compute the e vector */
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
  }

  /* e is in here, cant compute it else where*/
  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
      if (delta[ik_idx] != 0) {
	/* reset nu */
	for ( r = 0; r < *p; r++) {
	  nu[r] = 0.0;
	}
	de = 0.0;
	jl_idx = 0;
	for (j = 0; j < *n; j++) {
	  for (l = 0; l < clsize[j]; l++) {
	    rikjl = get_rikjl(X, sigma, clsize, p, N, n, ik_idx, jl_idx);
	    if (rikjl != 0) {
	      edif = e[jl_idx] - e[ik_idx];
	      z = sqrtn * edif / rikjl;
	      H = pnorm(z, 0.0, 1.0, 1, 0);
	      for (r = 0; r < *p; r++) {
		nu[r] += Z[jl_idx] * weight[jl_idx] * X[jl_idx + r * *N] * H;
	      }
	      de += weight[jl_idx] * Z[jl_idx] * H;
	    }
	    jl_idx++;
	  }
	}
	for ( r = 0; r < *p; r ++) {
	  sn[r] += weight[ik_idx] * pw[ik_idx] * Z[ik_idx] * (X[ik_idx + r * *N] - nu[r] / de);
	}
      }
      ik_idx++;
    }
  }
  //  Free(edif);
  Free(nu);
  //  Free(de);
  Free(e);
  /* for (r = 0; r < *p; r++) { */
  /*   sn[r] /= (*n * (*n - 1)); */
  /* } */
  sn;
}


void ulognsfun(double *beta,
	       double *Y, double *X, double *delta, int *clsize,
	       double *sigma,
	       int *n, // numbers of clusters
	       int *p, // ncol(X), or dimension
	       int *N, // nrow(X), or numbers of dataset. Also equals to sum(clsize)
	       //output
	       double *Z,
	       double *weight,
	       double *pw,
	       double *sn) {
  int i, j, k, l, ik_idx = 0, jl_idx, m, r;
  double *e = Calloc(*N, double), *nu = Calloc(*p, double);
  double rikjl, z, H, sqrtn = sqrt(*n), edif, de = 0.0;

  /* compute the e vector */
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
  }

  /* e is in here, cant compute it else where*/
  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
      if (delta[ik_idx] != 0) {
	/* reset nu */
	for ( r = 0; r < *p; r++) {
	  nu[r] = 0.0;
	}
	de = 0.0;
	jl_idx = 0;
	for (j = 0; j < *n; j++) {
	  for (l = 0; l < clsize[j]; l++) {
	    if (e[ik_idx] - e[jl_idx] <= 0) {
	      for ( r = 0; r < *p; r++) {
		nu[r] += X[jl_idx + r * *N] * weight[jl_idx] * Z[jl_idx];
	      }
	      de += Z[jl_idx] * weight[jl_idx];
	    }
	    jl_idx++;
	  }
	}  // end jl
	for (r =  0; r < *p; r++) {
	  sn[r] += weight[ik_idx] * pw[ik_idx] * Z[ik_idx] * (X[ik_idx + r * *N] - nu[r] / de);
	}
      }
      ik_idx++;
    }
  }
  Free(nu);
  Free(e);
  /* for (r = 0; r < *p; r++) { */
  /*   sn[r] /= (*n * (*n - 1)); */
  /* } */
  sn;
}


void lfun(double *beta,
	  double *Y, double *X, double *delta, int *clsize,
	  double *sigma,
	  int *n, // numbers of clusters
	  int *p, // ncol(X), or dimension
	  int *N, // nrow(X), or numbers of dataset. Also equals to sum(clsize)
	  double *weight,
	  double *Z,
	  //output
	  double *ln) {
  int i, j, k, l, ik_idx = 0, jl_idx, m;
  double eik, ejl, rikjl, edif, z, H, h, sqrtn = sqrt(*n);
  double *e = Calloc(*N, double);

  /* compute the e vector */
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
  }

  /* e is in here, cant compute it else where*/

  *ln = 0.0;
  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
      if (delta[ik_idx] != 0) {
	jl_idx = 0;
	for (j = 0; j < *n; j++) {
	  for (l = 0; l < clsize[j]; l++) {
	    rikjl = get_rikjl(X, sigma, clsize, p, N, n, ik_idx, jl_idx);
	    if (rikjl != 0) {
	      edif = e[jl_idx] - e[ik_idx];
	      z = sqrtn * edif / rikjl;
	      H = pnorm(z, 0.0, 1.0, 1, 0);
	      h = dnorm(z, 0.0, 1.0, 0);
	      *ln += weight[ik_idx] * weight[jl_idx] * Z[ik_idx] * Z[jl_idx] * (edif * H + rikjl * h / sqrtn);
	    }
	    jl_idx++;
	  }
	}
      }
      ik_idx++;
    }
  }
  Free(e);
  /* *ln /= (*n * (*n - 1)); */
  *ln;
}

void omegafun(double *beta,
	      double *Y,
	      double *X,
	      double *delta,
	      int *clsize,
	      int *n,
	      int *p,
	      int *N,
	      double *weight,
	      //output
	      double *omega) {
  int i, k, j, l, m, r, s, ik_idx = 0, jl_idx, il_idx, rs_idx, emk = 0, ind = 0, omega_idx;
  double *xdif = Calloc(*p, double), *e = Calloc(*N, double), *ksi = Calloc(*p * *N, double);

  // compute e
  for (i = 0; i < *N; i++) {
    e[i] = 0.0;
    for (j = 0; j < *p; j++) {
      e[i] += X[j * *N + i] * beta[j];
    }
    e[i] = Y[i] - e[i];
   }

  // compute ksi
  for (i = 0; i < *n; i++) {
    for (k = 0; k < clsize[i]; k++) {
      jl_idx = 0;
      for (j = 0; j < *n; j++) {
	for (l = 0; l < clsize[j]; l++) {

	  if (delta[ik_idx] != 0) {
	    if (e[ik_idx] - e[jl_idx] < 0) {
	      for (m = 0; m < *p; m++) {
		ksi[ik_idx + m * *N] += delta[ik_idx] * (X[ik_idx + m * *N] - X[jl_idx + m * *N]) * weight[jl_idx] / *n;  //changed
	      } // end for m
	    } // end if e - e <0
	  } // end if delta[ik_idx] != 0

	  if (delta[jl_idx] != 0) {
	    if (e[ik_idx] - e[jl_idx] >= 0) {
	      rs_idx = 0;
	      emk = 0;
	      for (r = 0; r < *n; r++) {
		for (s = 0; s < clsize[r]; s++) {
		  if (e[rs_idx] - e[jl_idx] >= 0) {
		    for (m = 0; m < *p; m++) {
		      xdif[m] += weight[rs_idx] * ( X[ik_idx + m * *N] - X[rs_idx + m * *N] );  //changed
		    } // end for m
		    emk++;
		  } // end if
		  rs_idx++;
		} // end for s
	      } // end for r

	      for (m = 0; m < *p; m++) {
		ksi[ik_idx + m * *N] -= xdif[m] / (*n * emk);
		xdif[m] = 0;
	      } // end for m
	    } // end if e -e >= 0
	  }// end if delta jl != 0
	  jl_idx++;
	} // end for l
      } // end for j
      ik_idx++;
    } // end for k1
  } // end for i

  for (i = 0; i < *n; i++) {
    ind += clsize[i];
    for (k = 0; k < clsize[i]; k++) {
      for (l = 0; l < clsize[i]; l++) {
	omega_idx = 0;
	for (r = 0; r < *p; r++) {
	  for (s = 0; s < *p; s++) {
	    omega[omega_idx] += ksi[ind + k - clsize[0] + r * *N] * ksi[ind + l - clsize[0] + s * *N];
	    omega_idx++;
	  } // end for r
	} // end for s
      } // end for l
    } // end for k
  } // end for i

  /* for (r = 0; r < *p * *p; r++) { */
  /* omega[r] /= *n; */
  /* } */

  omega;

  Free(e);
  Free(xdif);
  Free(ksi);
}
