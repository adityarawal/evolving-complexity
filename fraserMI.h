#include <stdio.h>
#include <iostream>
#include <math.h>
#include <algorithm>    // std::sort
#include <cmath>        // std::abs
#include <stdlib.h>
#include <memory.h>
#include <cstring>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <list>
#include <vector>
#include <string>
#include "minfomodr.h"

using namespace std;

void read_file_data(char *fname, std::vector <double> &X);

double fraser_entropy(const vector<double> &X, int feature_count); //Ranges between 0-max_entropy

double fraser_mutual_information(const vector<double> &X, const vector<double> &Y);
