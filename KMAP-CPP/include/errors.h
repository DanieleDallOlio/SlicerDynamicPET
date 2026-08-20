#ifndef __errors_h__
#define __errors_h__

#include <iostream>
#include <string>

static constexpr int ERROR_POSITIVE         = 1;  ///< Return Error
static constexpr int ERROR_PATTERN          = 2;  ///< Return Error
static constexpr int ERROR_MODEL            = 3; ///< Return Error
static constexpr int ERROR_NAN              = 4; ///< Return Error
static constexpr int ERROR_SIZE             = 5; ///< Return Error
static constexpr int ERROR_INTERP           = 6; ///< Return Error

/**
* @brief Raise error because max_iters is not positive.
*
* @param max_iters The given value of max_iters.
*/
inline void error_maxit (const long int & maxit)
{
  std :: cerr << "maximum iteration must be non-negative; given: " << maxit << std :: endl;
  std :: exit(ERROR_POSITIVE);
}

/**
* @brief Raise error because input filename is not found.
*
* @param filename The given value of filename.
*/
inline void error_file (const std :: string & filename)
{
  std :: cerr << "Input file not found! Given: " << filename << std :: endl;
  std :: exit(ERROR_PATTERN);
}


/**
* @brief Raise error because compartment model type is invalid
*
* @param model The given value of compartment model
*/
inline void error_model (const std :: string & model)
{
  std :: cerr << "Invalid model name. Given: " << model
              << ". Possible values are \"onetcm\", \"twotcm\"."
              << std :: endl;
  std :: exit(ERROR_MODEL);
}


/**
* @brief Raise error when a NaN is found in an object
*
* @param obj The object storing an NaN
*/
inline void error_nan (const std :: string & obj)
{
  std :: cerr << "NaN found. Variable: " << obj
              << ". Only non-NaN values accepted."
              << std :: endl;
  std :: exit(ERROR_NAN);
}

/**
* @brief Raise error when a NaN is found in an object
*
* @param obj The object storing an NaN
*/
inline void error_size (const std :: string & obj1,
                        const std :: string & obj2,
                        const int s1,
                        const int s2
                      )
{
  std :: cerr << "Different sizes (" << s1 << "," << s2
              << ") between: " << obj1
              << " and " << obj2
              << std :: endl;
  std :: exit(ERROR_SIZE);
}

/**
* @brief Raise error when an interpolation method is not available
*
* @param method The method
*/
inline void error_interpolation (const std :: string & method)
{
  std :: cerr << "Method: " << method
              << ". Not implemented yet."
              << std :: endl;
  std :: exit(ERROR_INTERP);
}

#endif // __errors_h__
