#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include <errno.h>
#include <sys/stat.h>
#include <vector>
#include <map>
#include <algorithm>
#include <math.h>
#include <float.h>
#include <utility>
#include "moment.estimator.h"

using namespace std;

#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
#define SIGMA_TINY 1.0e-20

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	/* log(sqrt(2*pi)) */
#endif

#ifndef M_LN_2PI
#define M_LN_2PI	1.837877066409345483560659472811	/* log(2*pi) */
#endif

double R_D_exp (double x, bool log_p) {
	return (log_p	? (x)	: exp(x));
}

/* This file contains the p-value calculation functions of LARVA */

// Sort a pair vector based on the value of its first value
bool reverseSortFirst (vector<double> a, vector<double> b) {
	return a[0] > b[0];
}

// Sort a pair vector based on the value of its second value
bool sortVectorsSecond (vector<double> a, vector<double> b) {
	return a[1] < b[1];
}

// Returns the value ln[gamma(xx)] for xx > 0.
double gammln(double xx) {
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
												24.01409824083091,-1.231739572450155,
												0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

double stirlerr(double n) {

#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508/* 1/1188 */

/*
  error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
*/
    const static double sferr_halves[31] = {
	0.0, /* n=0 - wrong, place holder only */
	0.1534264097200273452913848,  /* 0.5 */
	0.0810614667953272582196702,  /* 1.0 */
	0.0548141210519176538961390,  /* 1.5 */
	0.0413406959554092940938221,  /* 2.0 */
	0.03316287351993628748511048, /* 2.5 */
	0.02767792568499833914878929, /* 3.0 */
	0.02374616365629749597132920, /* 3.5 */
	0.02079067210376509311152277, /* 4.0 */
	0.01848845053267318523077934, /* 4.5 */
	0.01664469118982119216319487, /* 5.0 */
	0.01513497322191737887351255, /* 5.5 */
	0.01387612882307074799874573, /* 6.0 */
	0.01281046524292022692424986, /* 6.5 */
	0.01189670994589177009505572, /* 7.0 */
	0.01110455975820691732662991, /* 7.5 */
	0.010411265261972096497478567, /* 8.0 */
	0.009799416126158803298389475, /* 8.5 */
	0.009255462182712732917728637, /* 9.0 */
	0.008768700134139385462952823, /* 9.5 */
	0.008330563433362871256469318, /* 10.0 */
	0.007934114564314020547248100, /* 10.5 */
	0.007573675487951840794972024, /* 11.0 */
	0.007244554301320383179543912, /* 11.5 */
	0.006942840107209529865664152, /* 12.0 */
	0.006665247032707682442354394, /* 12.5 */
	0.006408994188004207068439631, /* 13.0 */
	0.006171712263039457647532867, /* 13.5 */
	0.005951370112758847735624416, /* 14.0 */
	0.005746216513010115682023589, /* 14.5 */
	0.005554733551962801371038690  /* 15.0 */
    };
    double nn;

    if (n <= 15.0) {
	nn = n + n;
	if (nn == (int)nn) return(sferr_halves[(int)nn]);
	return(gammln(n + 1.) - (n + 0.5)*log(n) + n - M_LN_SQRT_2PI);
    }

    nn = n*n;
    if (n>500) return((S0-S1/nn)/n);
    if (n> 80) return((S0-(S1-S2/nn)/nn)/n);
    if (n> 35) return((S0-(S1-(S2-S3/nn)/nn)/nn)/n);
    /* 15 < n <= 35 : */
    return((S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n);
}

double bd0(double x, double np) {
    double ej, s, s1, v;
    int j;

    if (isnan(x) || isnan(np) || np == 0.0) {
    	perror("Error: bd0 given an argument that is not a number.\n");
			exit(1);
		}

    if (fabs(x-np) < 0.1*(x+np)) {
	v = (x-np)/(x+np);  // might underflow to 0
	s = (x-np)*v;/* s using v -- change by MM */
	if(fabs(s) < DBL_MIN) return s;
	ej = 2*x*v;
	v = v*v;
	for (j = 1; j < 1000; j++) { /* Taylor series; 1000: no infinite loop
					as |v| < .1,  v^2000 is "zero" */
	    ej *= v;// = v^(2j+1)
	    s1 = s+ej/((j<<1)+1);
	    if (s1 == s) /* last term was effectively 0 */
		return s1 ;
	    s = s1;
	}
    }
    /* else:  | x - np |  is not too small */
    return(x*log(x/np)+np-x);
}

// Our own dbinom function
double dbinom (double x, double n, double p, bool log_p) {
	
	// Error checking
	if ((p < 0) || (p > 1)) {
		perror("p must be between 0 and 1\n");
		exit(1);
	}
	if (x < 0) {
		perror("x must be >=0\n");
		exit(1);
	}
	if (n < x) {
		perror("x must be <= than the binomial denominator\n");
		exit(1);
	}
	double q = 1 - p;
	double lf, lc;
	
	if (p == 0) return((x == 0) ? (log_p ? 0. : 1.) : (log_p ? -DBL_MAX : 0.));
  if (q == 0) return((x == n) ? (log_p ? 0. : 1.) : (log_p ? -DBL_MAX : 0.));
  
  if (x == 0) {
		if(n == 0) return (log_p ? 0. : 1.);
		lc = (p < 0.1) ? -bd0(n,n*q) - n*p : n*log(q);
		return( (log_p	?  (lc)	 : exp(lc)) );
  }
  if (x == n) {
		lc = (q < 0.1) ? -bd0(n,n*p) - n*q : n*log(p);
		return( (log_p	?  (lc)	 : exp(lc)) );
  }
  if (x < 0 || x > n) return( (log_p ? -DBL_MAX : 0.) );
  
  /* n*p or n*q can underflow to zero if n and p or q are small.  This
		 used to occur in dbeta, and gives NaN as from R 2.3.0.  */
	lc = stirlerr(n) - stirlerr(x) - stirlerr(n-x) - bd0(x,n*p) - bd0(n-x,n*q);

	/* f = (M_2PI*x*(n-x))/n; could overflow or underflow */
	/* Upto R 2.7.1:
	 * lf = log(M_2PI) + log(x) + log(n-x) - log(n);
	 * -- following is much better for  x << n : */
	lf = M_LN_2PI + log(x) + log1p(- x/n);

	return R_D_exp((lc - 0.5*lf), log_p);
}

// Used by betai: Evaluates continued fraction for incomplete beta function by modified Lentz's method.
double betacf(double a, double b, double x) {
	int m,m2;
	double aa,c,d,del,h,qab,qam,qap;
	
	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) {
		d=FPMIN;
	}
	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (m > MAXIT) {
		perror("a or b too big, or MAXIT too small in betacf");
		exit(1);
	}
	return h;
}

// Returns the incomplete beta function for I_x(a,b)
// To evaluate pbinom(q, n, p), use I_p(k, n-k+1) or betai(k, n-k+1, p)
double betai(double a, double b, double x) {

	double bt;

	if (x < 0.0 || x > 1.0) {
		perror("Bad x in routine betai");
		exit(1);
	}
	
	if (x == 0.0 || x == 1.0) {
		bt=0.0;
	} else {
		bt = exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
	}
	
	if (x < (a+1.0)/(a+b+2.0)) {
		return bt*betacf(a,b,x)/a;
	} else {
		return 1.0-bt*betacf(b,a,1.0-x)/b;
	}
}

double do_search(double y, double *z, double p, double n, double pr, double incr) {
	if (*z >= p) {
		/* search to the left */
		for(;;) {
			double newz;
	    // if (y == 0 || (newz = pbinom(y - incr, n, pr)) < p) {
	    if (y == 0 || (newz = betai(y - incr, n-(y-incr)+1, pr)) < p) {
	    	return y;
	    }
	    y = fmax(0, y - incr);
	    *z = newz;
		}
	} else {
		/* search to the right */
		for(;;) {
			y = fmin(y + incr, n);
			// if (y == n || (*z = pbinom(y, n, pr, /*l._t.*/TRUE, /*log_p*/FALSE)) >= p) {
			if (y == n || (*z = betai(y, n-y+1, pr)) >= p) {
				return y;
			}
		}
	}
}

// Helper function for rounding doubles
double round_to_digits(double value, int digits)
{
	if (value == 0.0) {
		return 0.0;
	} else if (value > 0.0 && value < 0.05) {
		return 0.0;
	} else {
    double factor = pow(10.0, digits - ceil(log10(fabs(value))));
    return round(value * factor) / factor;   
  }
}

// Helper function for rounding doubles so that there are always 3 s.f.
double final_rounding(double value) {
	if (isnan(value)) {
		return 0.0;
	} else {
		double factor = pow(10.0, 3);
		double res = round(value * factor) / factor;
  	return res;
  }
}

double R_D_Lval(double p, bool lower_tail) {
	return (lower_tail ? (p) : (0.5 - (p) + 0.5));
}

double R_DT_qIv(double p, bool lower_tail, bool log_p) {
	return (log_p ? (lower_tail ? exp(p) : - expm1(p)) : R_D_Lval(p, lower_tail));
}

double R_D_Cval(double p, bool lower_tail) {
	return (lower_tail ? (0.5 - (p) + 0.5) : (p));
}

double R_DT_CIv(double p, bool lower_tail, bool log_p) {
	return (log_p ? (lower_tail ? -expm1(p) : exp(p)) : R_D_Cval(p, lower_tail));
}

double qnorm(double p, double mu, double sigma, bool lower_tail, bool log_p) {
	double p_, q, r, val;
	
	// Error checking
	if (isnan(p) || isnan(mu) || isnan(sigma)) {
		perror("Error: qnorm given an argument that is not a number.\n");
		exit(1);
	}
	
	if (sigma < SIGMA_TINY) {
		sigma = SIGMA_TINY;
	}
	
	if (sigma < 0) {
		perror("Error: qnorm given a negative sigma\n");
		exit(1);
	} else if (sigma == 0) {
		return mu;
	}
	
	p_ = R_DT_qIv(p, true, false);/* real lower_tail prob. p */
  q = p_ - 0.5;
  
  if (fabs(q) <= .425) {/* 0.075 <= p <= 0.925 */
  	r = .180625 - q * q;
  	val =
            q * (((((((r * 2509.0809287301226727 +
                       33430.575583588128105) * r + 67265.770927008700853) * r +
                     45921.953931549871457) * r + 13731.693765509461125) * r +
                   1971.5909503065514427) * r + 133.14166789178437745) * r +
                 3.387132872796366608)
            / (((((((r * 5226.495278852854561 +
                     28729.085735721942674) * r + 39307.89580009271061) * r +
                   21213.794301586595867) * r + 5394.1960214247511077) * r +
                 687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
  } else { /* closer than 0.075 from {0,1} boundary */
  	/* r = min(p, 1-p) < 0.075 */
		if (q > 0) {
	    r = R_DT_CIv(p, lower_tail, log_p);/* 1-p */
		} else {
	    r = p_;/* = R_DT_Iv(p) ^=  p */
	  }
	  
	  r = sqrt(- ((log_p &&
		     ((lower_tail && q <= 0) || (!lower_tail && q > 0))) ?
		    p : /* else */ log(r)));
        /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */
        
  	if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
				r += -1.6;
				val = (((((((r * 7.7454501427834140764e-4 +
									 .0227238449892691845833) * r + .24178072517745061177) *
								 r + 1.27045825245236838258) * r +
								3.64784832476320460504) * r + 5.7694972214606914055) *
							r + 4.6303378461565452959) * r +
						 1.42343711074968357734)
						/ (((((((r *
										 1.05075007164441684324e-9 + 5.475938084995344946e-4) *
										r + .0151986665636164571966) * r +
									 .14810397642748007459) * r + .68976733498510000455) *
								 r + 1.6763848301838038494) * r +
								2.05319162663775882187) * r + 1.);
		} else { /* very close to  0 or 1 */
				r += -5.;
				val = (((((((r * 2.01033439929228813265e-7 +
									 2.71155556874348757815e-5) * r +
									.0012426609473880784386) * r + .026532189526576123093) *
								r + .29656057182850489123) * r +
							 1.7848265399172913358) * r + 5.4637849111641143699) *
						 r + 6.6579046435011037772)
						/ (((((((r *
										 2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
										r + 1.8463183175100546818e-5) * r +
									 7.868691311456132591e-4) * r + .0148753612908506148525)
								 * r + .13692988092273580531) * r +
								.59983220655588793769) * r + 1.);
		}
		if(q < 0.0) {
	    val = -val;
	  }
    /* return (q >= 0.)? r : -r ;*/
  }
  return mu + sigma * val;
}

// pbinom implementation, upper tail only
double pbinom (int x, int n, double p) {
	double result = 0;
	for (int i = x; i <= n; i++) {
		result += dbinom((double)i, (double)n, p, false);
	}
	return result;
}

// Benjamini-Hochberg implementation
// Assumes we always use n = length(p)
vector<double> bh_adjust (vector<double> p) {
	// Save only the non-NaNs
	// vector<double> p_new;
	vector<vector<double> > p_ind;
	int p_ind_index = 0;
	for (unsigned int i = 0; i < p.size(); i++) {
		if (!isnan(p[i])) {
			vector<double> temp;
			temp.push_back(p[i]);
			temp.push_back((double)p_ind_index);
			p_ind_index++;
			p_ind.push_back(temp);
		}
	}
	
	// Sort p_ind by decreasing order
	sort(p_ind.begin(), p_ind.end(), reverseSortFirst);
	
	// Do the n/i calculation
	// Dot product between seq and p_ind
	for (unsigned int i = 0; i < p_ind.size(); i++) {
		
		unsigned int i_inv = p_ind.size() - i;
		double temp = (double)p.size()/(double)i_inv;
		p_ind[i][0] = temp*(p_ind[i][0]);
	}
	
	// Cumulative minimum
	double curmin;
	for (unsigned int i = 0; i < p_ind.size(); i++) {
	
		// Initialization
		if (i == 0) {
			curmin = p_ind[i][0];
		}
	
		double this_iter = fmin(curmin, p_ind[i][0]);
		
		// Update curmin if necessary
		if (p_ind[i][0] < curmin) {
			curmin = p_ind[i][0];
		}
		
		// Don't let this value exceed 1
		if (this_iter > 1) {
			this_iter = 1.0;
		}
		
		p_ind[i][0] = this_iter;
	}
	
	sort(p_ind.begin(), p_ind.end(), sortVectorsSecond);
	
	// Save the firsts
	vector<double> retval;
	for (unsigned int i = 0; i < p_ind.size(); i++) {
		retval.push_back(p_ind[i][0]);
	}
	return retval;
}

// qbinom implementation
double qbinom(double p, double n, double pr, bool lower_tail, bool log_p) {
	double q, mu, sigma, gamma, z, y;
	
	// Error checking
	if (isnan(p) || isnan(n) || isnan(pr)) {
		perror("Error: Input to qbinom is NaN\n");
		exit(1);
	}
	if (pr < 0 || pr > 1 || n < 0) {
		perror("Error: Input to qbinom is out of bounds\n");
		exit(1);
	}
	
	if (pr == 0. || n == 0) return 0.;
	
	q = 1 - pr;
	if(q == 0.) return n; /* covers the full range of the distribution */
	mu = n * pr;
	sigma = sqrt(n * pr * q);
	gamma = (q - pr) / sigma;
	
	if (!lower_tail || log_p) {
		p = R_DT_qIv(p, true, false);
		if (p == 0.) return 0.;
		if (p == 1.) return n;
  }
  if (p + 1.01*DBL_EPSILON >= 1.) {
  	return n;
  }
  z = qnorm(p, 0., 1., /*lower_tail*/true, /*log_p*/false);
  y = floor(mu + sigma * (z + gamma * (z*z - 1) / 6) + 0.5);
  
  if (y > n) { /* way off */
  	y = n;
  }
  
  z = betai(y, n-y+1, pr);
  
  /* fuzz to ensure left continuity: */
  p *= 1 - 64*DBL_EPSILON;
  
  if (n < 1e5) {
  	return do_search(y, &z, p, n, pr, 1);
  } else { /* Otherwise be a bit cleverer in the search */
  	double incr = floor(n * 0.001), oldincr;
  	do {
	    oldincr = incr;
	    y = do_search(y, &z, p, n, pr, incr);
	    incr = fmax(1, floor(incr/100));
		} while(oldincr > 1 && incr > n*1e-15);
		return y;
	}
}

// Density of BBD distribution
double dBB (double x, double mu, double sigma, double bd, bool log_p) {

	if (sigma < SIGMA_TINY) {
		sigma = SIGMA_TINY;
	}

	if (mu <= 0 || mu >= 1) {
		perror("dBB: mu must be between 0 and 1\n");
		exit(1);
	}
	if (sigma <= 0) {
		perror("sigma must be greater than 0\n");
		exit(1);
	}
	if (x < 0) {
		perror("x must be >=0\n");
		exit(1);
	}
	if (bd < x) {
		perror("x must be <= than the binomial denominator\n");
		exit(1);
	}
	if (sigma < SIGMA_TINY) {
		sigma = SIGMA_TINY;
	}
	double logfy = (gammln(bd + 1) - gammln(x + 1) - gammln(bd - x + 
        					1) + gammln((1/sigma)) + gammln(x + mu * (1/sigma)) + 
        					gammln(bd + ((1 - mu)/sigma) - x) - gammln(mu * (1/sigma)) - 
        					gammln((1 - mu)/sigma) - gammln(bd + (1/sigma)));
	if (sigma < 0.0001) {
		logfy = dbinom(x, bd, mu, true /*log_p*/);
	}
	if (!log_p) {
		logfy = exp(logfy);
	}
	return logfy;
}

// CDF of BBD distribution
double pBB (double q, double mu, double sigma, double bd, bool lower_tail, bool log_p) {
	
	if (sigma < SIGMA_TINY) {
		sigma = SIGMA_TINY;
	}

	if (mu <= 0 || mu >= 1) {
		perror("pBB: mu must be between 0 and 1\n");
		exit(1);
	}
	if (sigma <= 0) {
		perror("sigma must be greater than 0\n");
		exit(1);
	}
	if (q < 0) {
		perror("y must be >=0\n");
		exit(1);
	}
	if (bd < q) {
		perror("y must be <= than the binomial denominator\n");
		exit(1);
	}
	double FFF = 0.0;
	double pdfall_sum = 0;
	double allval_max = floor(q);
	for (int i = 0; i < (int)allval_max; i++) {
		double pdfall = dBB((double)i, mu, sigma, bd, false);
		pdfall_sum += pdfall;
	}
	FFF = pdfall_sum;
	double cdf = FFF;
	if (!(lower_tail)) {
		cdf = 1 - cdf;
	}
	if (log_p) {
		cdf = log(cdf);
	}
	if (sigma <= 0.0001) {
		cdf = pbinom(q, bd, mu);
	}
	return cdf;
}

// P-value estimation for equal length bins
vector<vector<double> > p_value_equal_len (vector<int> x, vector<int> len) {

	// Truncate any x greater than len
	for (unsigned int i = 0; i < x.size(); i++) {
		if (x[i] > len[i]) {
			x[i] = len[i];
		}
	}

	// Moment estimation
	vector<double> p;
	for (unsigned int i = 0; i < x.size(); i++) {
		double temp = (double)x[i]/(double)len[i];
		p.push_back(temp);
	}
	vector<double> re = bbd_moment_estimator_equal_length(p, len);
	
	// MLE calculations
	double p_mle = 0;
	for (unsigned int i = 0; i < p.size(); i++) {
		p_mle += p[i];
	}
	p_mle = p_mle/(double)p.size();
	
	// pbinom calculations
	vector<double> p_binomial;
	for (unsigned int i = 0; i < x.size(); i++) {
		// double b = (double)len[i] - (double)x[i] + 1;
		// double this_p_binomial = betai((double)x[i], b, p_mle);
		double this_p_binomial;
		if (p_mle == 0) {
			this_p_binomial = 1;
		} else {
			this_p_binomial = pbinom(x[i], len[i], p_mle);
		}
		p_binomial.push_back(this_p_binomial);
	}
	
	// pBB calculations
	vector<double> p_bb;
	for (unsigned int i = 0; i < x.size(); i++) {
		int this_x = x[i];
// 		if (this_x == 0) {
// 			this_x = 1;
// 		}
// 		this_x--;
		double this_p_bb;
		if (re[0] == 0) {
			this_p_bb = 1;
		} else {
			this_p_bb = pBB((double)this_x, re[0], re[4], (double)len[i], false, false);
		}
		p_bb.push_back(this_p_bb);
	}
	
	// Prepare return values
	vector<vector<double> > retval;
	retval.push_back(p_binomial);
	retval.push_back(p_bb);
	return retval;
}

vector<vector<double> > pval_fixed_length_rep_time_correction (vector<int> x, vector<int> len, vector<double> rep_time) {

	// Truncate any x greater than len
	for (unsigned int i = 0; i < x.size(); i++) {
		if (x[i] > len[i]) {
			x[i] = len[i];
		}
	}

	vector<double> rep_index;
	for (unsigned int i = 0; i < rep_time.size(); i++) {
		double transformed_reptime = round_to_digits(rep_time[i], 1);
		transformed_reptime = round(transformed_reptime*10.0);
		rep_index.push_back(transformed_reptime);
	}
	
	// fit the local model
	vector<double> miu_val;
	vector<double> sigma_val;
	vector<double> p_mle;
	
	for (unsigned int i = 0; i <= 10; i++) {
		
		// Step through the rep_index vector, and find those with reptiming rank i
		// Then get the corresponding x (mut count) and len, save to vector for moment estimation
		vector<int> rank_x;
		vector<int> rank_len;
		// Divide x by len
		vector<double> rank_p;
		int rank_x_sum = 0;
		int rank_len_sum = 0;
		
		for (unsigned int j = 0; j < rep_index.size(); j++) {
			if (rep_index[j] == i) {
				rank_x.push_back(x[j]);
				rank_len.push_back(len[j]);
				rank_x_sum += x[j];
				double this_rank_p = (double)x[j]/(double)len[j];
				rank_p.push_back(this_rank_p);
				rank_len_sum += len[j];
			}
		}
		
		if (rank_x_sum > 0) {
		
			vector<double> re = bbd_moment_estimator_equal_length(rank_p, rank_len);
			miu_val.push_back(re[0]);
			sigma_val.push_back(re[4]);
			double this_mle = (double)rank_x_sum/(double)rank_len_sum;
			p_mle.push_back(this_mle);
		} else {
			miu_val.push_back(0.0);
			sigma_val.push_back(0.0);
			p_mle.push_back(0.0);
		}
	}
	
	// p-value calculation
	vector<double> p_binomial_correction;
	for (unsigned int i = 0; i < x.size(); i++) {
		double this_p_binomial;
		if (p_mle[rep_index[i]] == 0) {
			this_p_binomial = 1;
		} else {
			this_p_binomial = pbinom(x[i], len[i], p_mle[rep_index[i]]);
		}
		p_binomial_correction.push_back(this_p_binomial);
	}
	
	vector<double> p_bbd_correction;
	for (unsigned int i = 0; i < x.size(); i++) {
		int this_x = x[i];
// 		if (this_x == 0) {
// 			this_x = 1;
// 		}
// 		this_x--;
		double this_p_bb;
		if (miu_val[rep_index[i]] == 0) {
			this_p_bb = 1;
		} else {
			this_p_bb = pBB((double)this_x, miu_val[rep_index[i]], sigma_val[rep_index[i]], (double)len[i], false, false);
		}
		p_bbd_correction.push_back(this_p_bb);
	}
	
	vector<vector<double> > retval;
	retval.push_back(p_binomial_correction);
	retval.push_back(p_bbd_correction);
	retval.push_back(p_mle);
	return retval;
}

vector<vector<double> > pval_varying_length (vector<int> x, vector<int> len) {

	// Truncate any x greater than len
	for (unsigned int i = 0; i < x.size(); i++) {
		if (x[i] > len[i]) {
			x[i] = len[i];
		}
	}

	// Moment estimation
	vector<double> p;
	double sum_x = 0;
	double sum_len = 0;
	for (unsigned int i = 0; i < x.size(); i++) {
		double temp = (double)x[i]/(double)len[i];
		p.push_back(temp);
		
		sum_x += (double)x[i];
		sum_len += (double)len[i];
	}
	
	vector<double> re = bbd_moment_estimator_inequal_length(p, len);
	
	double p_mle = sum_x/sum_len;
	vector<double> p_binomial;
	for (unsigned int i = 0; i < x.size(); i++) {
		double this_p_binomial;
		if (p_mle == 0) {
			this_p_binomial = 1;
		} else {
			this_p_binomial = pbinom(x[i], len[i], p_mle);
		}
		p_binomial.push_back(this_p_binomial);
	}
	
	vector<double> p_bbd;
	for (unsigned int i = 0; i < x.size(); i++) {
		int this_x = x[i];
// 		if (this_x == 0) {
// 			this_x = 1;
// 		}
// 		this_x--;
		double this_p_bb;
		if (re[0] == 0) {
			this_p_bb = 1;
		} else {
			this_p_bb = pBB((double)this_x, re[0], re[4], (double)len[i], false, false);
		}
		p_bbd.push_back(this_p_bb);
	}
	
	// Prepare return values
	vector<vector<double> > retval;
	retval.push_back(p_binomial);
	retval.push_back(p_bbd);
	return retval;
}

vector<vector<double> > pval_varying_length_rep_time_correction (vector<int> x, vector<int> len, vector<double> rep_time) {

	// Truncate any x greater than len
	for (unsigned int i = 0; i < x.size(); i++) {
		if (x[i] > len[i]) {
			x[i] = len[i];
		}
	}

	vector<double> rep_index;
	for (unsigned int i = 0; i < rep_time.size(); i++) {
	
		if (rep_time[i] == 0) {
			rep_index.push_back(0.0);
		} else {
	
			double transformed_reptime = round_to_digits(rep_time[i], 1);
			transformed_reptime = round(transformed_reptime*10.0);
		
			rep_index.push_back(transformed_reptime);
		}
	}
	
	// fit the local model
	vector<double> miu_val;
	vector<double> sigma_val;
	vector<double> p_mle;
	
	for (unsigned int i = 0; i <= 10; i++) {
	
		// Step through the rep_index vector, and find those with reptiming rank i
		// Then get the corresponding x (mut count) and len, save to vector for moment estimation
		vector<int> rank_x;
		vector<int> rank_len;
		// Divide x by len
		vector<double> rank_p;
		int rank_x_sum = 0;
		int rank_len_sum = 0;
		
		for (unsigned int j = 0; j < rep_index.size(); j++) {
			if (rep_index[j] == i) {
				rank_x.push_back(x[j]);
				rank_len.push_back(len[j]);
				rank_x_sum += x[j];
				double this_rank_p = (double)x[j]/(double)len[j];
				rank_p.push_back(this_rank_p);
				rank_len_sum += len[j];
			}
		}
		
		if (rank_x_sum > 0) {
			vector<double> re = bbd_moment_estimator_inequal_length(rank_p, rank_len);
			miu_val.push_back(re[0]);
			sigma_val.push_back(re[4]);
			double this_mle = (double)rank_x_sum/(double)rank_len_sum;
			p_mle.push_back(this_mle);
		} else {
			miu_val.push_back(0.0);
			sigma_val.push_back(0.0);
			p_mle.push_back(0.0);
		}
	}
	
	// p-value calculation
	vector<double> p_binomial_correction;
	for (unsigned int i = 0; i < x.size(); i++) {
		double this_p_binomial;
		if (p_mle[rep_index[i]] == 0) {
			this_p_binomial = 1;
		} else {
			this_p_binomial = pbinom(x[i], len[i], p_mle[rep_index[i]]);
		}
		p_binomial_correction.push_back(this_p_binomial);
	}
	
	vector<double> p_bbd_correction;
	for (unsigned int i = 0; i < x.size(); i++) {
		
		int this_x = x[i];
// 		if (this_x == 0) {
// 			this_x = 1;
// 		}
// 		this_x--;
		double this_p_bb;
		if (miu_val[rep_index[i]] == 0) {
			this_p_bb = 1;
		} else {
			this_p_bb = pBB((double)this_x, miu_val[rep_index[i]], sigma_val[rep_index[i]], (double)len[i], false, false);
		}
		p_bbd_correction.push_back(this_p_bb);
	}
	
	vector<vector<double> > retval;
	retval.push_back(p_binomial_correction);
	retval.push_back(p_bbd_correction);
	retval.push_back(p_mle);
	return retval;
}
