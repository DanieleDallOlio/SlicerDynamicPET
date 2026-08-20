#ifndef __utils_hpp__
#define __utils_hpp__

#include <utils.h>
#include <errors.h>
#include <fstream>
#include <sstream>
#include <algorithm>

template < typename T >
void read_matrix (std::vector< std::vector<T> > & matrix,
                  long int & Nrow,
                  long int & Ncol,
                  const std :: string & filename,
                  const std :: string & del) {
  #ifdef DEBUG
  long int ncol;
  #endif
  std :: vector < std :: string > row_;
  std :: ifstream is(filename);
  if ( !is ) error_file(filename);

  std :: stringstream buff;
  buff << is.rdbuf();
  is.close();

  std :: string line;
  std::getline(buff, line);
  row_ = split(line, del);
  Ncol = row_.size();
  std::vector<T> row_values(Ncol);
  std :: transform( row_.begin(), row_.end(), row_values.begin(), [](std :: string & i){return std :: stod(i);} );
  matrix.push_back(row_values);
  long int rowi = 1L;
  while (std::getline(buff, line)) {
    row_ = split(line, del);
    #ifdef DEBUG
    ncol = row_.size();
    assert (ncol == Ncol);
    #endif
    std :: transform( row_.begin(), row_.end(), row_values.begin(), [](std :: string & i){return std :: stod(i);} );
    matrix.push_back(row_values);
    ++rowi;
  }
  Nrow = rowi;
}

template < typename T >
void print_matrix (T matrix,
                   long int Nrow,
                   long int Ncol) {
  for (long int i = 0L; i < Nrow; ++i) {
    for (long int j = 0L; j < Ncol; ++j) {
      std :: cout << matrix[i][j] << "\t";
    }
    std :: cout << std :: endl;
  }
}

template < typename T >
void print_vec (T vec,
                long int N) {
  for (long int i = 0L; i < N; ++i) {
      std :: cout << vec[i] << "\t";
  }
  std :: cout << std :: endl;
}

#endif // __utils_hpp__
