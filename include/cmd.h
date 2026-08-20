#ifndef __cmd_h__
#define __cmd_h__

#include <errors.h>
#include <parse_args.hpp>
#include <string>
#include <utils.h>

void parse_test_kmap_levmar(int argc, char *argv[],
                            std :: string & cfile,
                            std :: string & wfile,
                            std :: string & sfile,
                            std :: string & bfile,
                            std :: string & wbfile,
                            double & dk,
                            std :: string & pinit,
                            std :: string & lb,
                            std :: string & ub,
                            std :: string & psens,
                            int & maxit,
                            long int & td,
                            int & nth,
                            std :: string & del
                          );

// void parse_read_matrix(int argc, char *argv[],
//                        std :: string & cfile,
//                        std :: string & del
//                       );

void parse_demo_tcm_realdata(int argc, char *argv[],
                             std :: string & cfile,
                             std :: string & bfile,
                             std :: string & tdfile,
                             int & nth,
                             std :: string & del
                            );

#endif // __cmd_h__
