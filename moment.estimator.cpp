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

/* This file contains the moment estimation functions of LARVA */

/* Parameter change from miu gamma to alfa and beta */
vector<double> miu2alfa (double miu_val, double gamma_val) {
	double alfa_val = miu_val*(1.0-gamma_val)/gamma_val;
	double beta_val = (1.0-gamma_val)*(1.0-miu_val)/gamma_val;
	vector<double> retval;
	retval.push_back(alfa_val);
	retval.push_back(beta_val);
	return retval;
}

vector<double> bbd_moment_estimator_equal_length (vector<double> p, vector<int> n) {
	
	// Mean of p
	double p_mean = 0;
	for (unsigned int i = 0; i < p.size(); i++) {
		p_mean += p[i];
	}
	p_mean = p_mean/(double)p.size();
	
	int k = p.size();
	
	// Sum of squares calculation
	double s = 0;
	for (unsigned int i = 0; i < p.size(); i++) {
		double residual = p[i] - p_mean;
		residual *= residual;
		s += residual;
	}
	
	double miu_val = p_mean;
	double gamma_val;
	for (unsigned int i = 0; i < n.size(); i++) {
		double cand = n[i]/(n[i]-1.0)*s/(miu_val)/(1.0-miu_val)/((double)k-1.0) - 1.0/(n[i]-1.0);
		if (i == 0) {
			gamma_val = cand;
		} else {
			if (cand > gamma_val) {
				gamma_val = cand;
			}
		}
	}
	
	if (miu_val < 0) {
		miu_val = 0;
	}
	
	if (miu_val > 1) {
		miu_val = 1;
	}
	
	vector<double> re = miu2alfa(miu_val, gamma_val);
	double alfa_val = re[0];
	double beta_val = re[1];
	double sigma_val = gamma_val/(1.0-gamma_val);
	
	// Prepare output
	vector<double> retval;
	retval.push_back(miu_val);
	retval.push_back(gamma_val);
	retval.push_back(alfa_val);
	retval.push_back(beta_val);
	retval.push_back(sigma_val);
	return retval;
}

vector<double> bbd_moment_estimator_inequal_length (vector<double> p, vector<int> n) {
	double miu_val = 0;
	double gamma_val = 0;
	
	// first round: wi equal to ni
	vector<double> wi;
	for (unsigned int i = 0; i < n.size(); i++) {
		wi.push_back((double)n[i]);
	}
	double sum_wi = 0;
	double p_bar = 0;
	
	// Summations
	for (unsigned int i = 0; i < wi.size(); i++) {
		sum_wi += wi[i];
		p_bar += wi[i]*p[i];
	}
	
	p_bar = p_bar/sum_wi;
	double q_bar = 1.0 - p_bar;
	
	// Sum of squares calculation
	double ss = 0;
	for (unsigned int i = 0; i < p.size(); i++) {
		double value = p[i]-p_bar;
		value *= value;
		value *= wi[i];
		ss += value;
	}
	
	double first_sum = 0;
	double second_sum = 0;
	for (unsigned int i = 0; i < wi.size(); i++) {
		double temp = wi[i]/sum_wi;
		double temp2 = 1.0 - temp;
		first_sum += temp2;
		
		double temp3 = wi[i]*temp2;
		second_sum += temp3;
	}
	
	double gamma_bar_up = ss - (p_bar * q_bar * first_sum);
	double gamma_bar_down = p_bar * q_bar * (second_sum - first_sum);
	double gamma_bar = gamma_bar_up/gamma_bar_down;
	
	if (isnan(gamma_bar) || gamma_bar < 0) {
		gamma_bar = 0;
	}
	
	// empirical weight
	for (unsigned int i = 0; i < wi.size(); i++) {
		wi[i] = n[i]/(1+gamma_bar*(n[i]-1));
	}
	
	// Summations
	sum_wi = 0;
	p_bar = 0;
	for (unsigned int i = 0; i < wi.size(); i++) {
		sum_wi += wi[i];
		p_bar += wi[i]*p[i];
	}
	
	p_bar = p_bar/sum_wi;
	q_bar = 1.0 - p_bar;
	
	// Sum of squares calculation
	ss = 0;
	for (unsigned int i = 0; i < p.size(); i++) {
		double value = p[i]-p_bar;
		value *= value;
		value *= wi[i];
		ss += value;
	}
	
	miu_val = p_bar;
	
	first_sum = 0;
	second_sum = 0;
	for (unsigned int i = 0; i < wi.size(); i++) {
		double temp = wi[i]/sum_wi;
		double temp2 = 1.0 - temp;
		double temp3 = (wi[i]/n[i])*temp2;
		first_sum += temp3;
		
		double temp4 = wi[i]*temp2;
		second_sum += temp4;
	}
	
	gamma_bar_up = ss - (p_bar * q_bar * first_sum);
	gamma_bar_down = p_bar * q_bar * (second_sum - first_sum);
	gamma_val = gamma_bar_up/gamma_bar_down;
	
	if (isnan(gamma_val) || gamma_val < 0) {
		gamma_val = 0;
	}
	
	vector<double> re = miu2alfa(miu_val, gamma_val);
	double alfa_val = re[0];
	double beta_val = re[1];
	double sigma_val = gamma_val/(1-gamma_val);
	
	// Prepare return values
	vector<double> retval;
	retval.push_back(miu_val);
	retval.push_back(gamma_val);
	retval.push_back(alfa_val);
	retval.push_back(beta_val);
	retval.push_back(sigma_val);
	return retval;
}
