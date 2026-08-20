#ifndef __utils_h__
#define __utils_h__

#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>
#include <memory>
#include <cmath>
#include <cassert>
#include <errors.h>



#include <chrono>


static constexpr double INF     = std :: numeric_limits < double > :: infinity(); ///< shortcut to infinity value
static constexpr double epsilon = std :: numeric_limits < double > :: epsilon(); ///< shortcut to epsilon value


/**
* @brief split string to token
*
* @param txt input string
* @param del delimiter as string
*
* @returns std::vector<std::string> vector of token
*
*/
std :: vector < std :: string > split (const std :: string & txt, const std :: string & del);

/**
* @brief check if a given file exists
*
* @param filename filename/path to check
*
* @returns bool true if the filename is found else false
*
*/
bool file_exists (const std :: string & filename);

bool containsNaN(const std::vector<std::vector<double>>& matrix);

double* finesample(double** & scant,
                   std :: vector < std :: vector <double> > & blood,
                   long int N,
                   long int & Nint,
                   long int res = 1L,
                   const std::string & itype = "linear");

double* finesample2(const std::vector<double*>& scant,
                   const std::vector<double>& blood,
                   long int& Nint,
                   long int res,
                   const std::string& itype);

#endif // __utils_h__
