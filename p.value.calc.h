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

using namespace std;

double round_to_digits (double value, int digits);
double final_rounding (double value);
vector<double> bh_adjust (vector< double > p);
vector<vector<double> > p_value_equal_len (vector<int> x, vector<int> len);
vector<vector<double> > pval_fixed_length_rep_time_correction (vector<int> x, vector<int> len, vector<double> rep_time);
vector<vector<double> > pval_varying_length (vector<int> x, vector<int> len);
vector<vector<double> > pval_varying_length_rep_time_correction (vector<int> x, vector<int> len, vector<double> rep_time);
double pBB (double q, double mu, double sigma, double bd, bool lower_tail, bool log_p);
double pbinom (int x, int n, double p);
