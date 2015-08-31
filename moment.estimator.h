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

using namespace std;

vector<double> bbd_moment_estimator_equal_length (vector<double> p, vector<int> n);
vector<double> bbd_moment_estimator_inequal_length (vector<double> p, vector<int> n);
