#include <cmd.h>

/**
* @brief Parse command line for running kmap_levmar procedure using parse_args library.
*
* @note This function is just an utility for the usage of the library from command line interface.
* It can be substituted by any equivalent command-line parser library.
* The variable are set into the function and thus are all passed by reference.
*
* @param[in] argc number of arguments in command line.
* @param[in] argv list of arguments in command line.
* @param[out] cfile filename storing time activity curve (TAC) data.
* @param[out] wfile filename storing weights for the TAC data.
* @param[out] sfile filename storing scan time data.
* @param[out] bfile filename storing blood data.
* @param[out] wbfile filename storing whole blood data.
* @param[out] dk decay constant.
* @param[out] pinit filename storing initial parameters for the model.
* @param[out] lb filename storing lower bounds for the parameters.
* @param[out] ub filename storing upper bounds for the parameters.
* @param[out] psens filename storing sensitivity matrix for the parameters.
* @param[out] maxit maximum number of iterations.
* @param[out] td time duration for the scan.
*
*/
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
                          )
{
  ArgumentParser argparse("Kinetic model fit");// + std :: to_string(MAJOR) + "." + std :: to_string(MINOR) + "." + std :: to_string(REVISION));

// #ifdef _OPENMP
//   nth  = omp_get_max_threads();
//   nth -= nth % 2;
//   argparse.add_argument < int >(          "tArg",  "t",  "threads", "Max number of threads exploitable",           false,  nth);
//
// #else
  nth = 1;

// #endif

  argparse.add_argument < std :: string >("cArg",  "c",  "cfile",   "TAC filename (with extension)",                true,  "");
  argparse.add_argument < std :: string >("wArg",  "w",  "wfile",   "Weights filename (with extension)",            true,  "");
  argparse.add_argument < std :: string >("sArg",  "s",  "sfile",   "Scan time filename (with extension)",          true,  "");
  argparse.add_argument < std :: string >("bArg",  "b",  "bfile",   "Blood data filename (with extension)",         true,  "");
  argparse.add_argument < std :: string >("wbArg", "wb", "wbfile",  "Whole Blood data filename (with extension)",   true,  "");
  argparse.add_argument < double >(       "dkArg", "dk", "decay",   "Decay constant(default: 0.1)",                 false, 0.);
  argparse.add_argument < std :: string >("pArg",  "p",  "pinit",   "Initial parameters filename (with extension)", true,  "");
  argparse.add_argument < std :: string >("lbArg", "lb", "lower",   "Lower bounds filename (with extension)",       true,  "");
  argparse.add_argument < std :: string >("ubArg", "ub", "upper",   "Upper bounds filename (with extension)",       true,  "");
  argparse.add_argument < std :: string >("psArg", "ps", "psens",   "Sensitivity matrix filename (with extension)", true,  "");
  argparse.add_argument < int >(          "iArg",  "i",  "iter",    "Maximum number of iterations (default:10)",    false, 10);
  argparse.add_argument < long int >(     "tdArg",  "td", "delay",  "Time delay (default:0)",                       false, 0L);
  argparse.add_argument < std :: string >("dmArg",  "dm", "delim",  "Delimiter (default: \"\\t\")",                 false, "\t");

  argparse.parse_args(argc, argv);

// #ifdef _OPENMP
//
//   argparse.get < int >(          "tArg",  nth);
//
// #endif
  argparse.get < std :: string >("cArg",  cfile);
  argparse.get < std :: string >("wArg",  wfile);
  argparse.get < std :: string >("sArg",  sfile);
  argparse.get < std :: string >("bArg",  bfile);
  argparse.get < std :: string >("wbArg", wbfile);
  argparse.get < double >(       "dkArg", dk);
  argparse.get < std :: string >("pArg",  pinit);
  argparse.get < std :: string >("lbArg", lb);
  argparse.get < std :: string >("ubArg", ub);
  argparse.get < std :: string >("psArg", psens);
  argparse.get < int >(          "iArg",  maxit);
  argparse.get < long int >(     "tdArg", td);
  argparse.get < std :: string >("dmArg", del);

  if( !file_exists(cfile) )  error_file(cfile);
  if( !file_exists(wfile) )  error_file(wfile);
  if( !file_exists(sfile) )  error_file(sfile);
  if( !file_exists(bfile) )  error_file(bfile);
  if( !file_exists(wbfile) ) error_file(wbfile);
  if( !file_exists(pinit) )  error_file(pinit);
  if( !file_exists(lb) )     error_file(lb);
  if( !file_exists(ub) )     error_file(ub);
  if( !file_exists(psens) )  error_file(psens);

  return;
}


// void parse_read_matrix(int argc, char *argv[],
//                        std :: string & cfile,
//                        std :: string & del
//                       )
// {
//   ArgumentParser argparse("Parse read matrix");
//   argparse.add_argument < std :: string >("cArg",  "c",  "cfile",   "TAC filename (with extension)",                true,  "");
//   argparse.add_argument < std :: string >("dmArg",  "dm", "delim",  "Delimiter (default: \"\\t\")",                 false, "\t");
//   argparse.parse_args(argc, argv);
//   argparse.get < std :: string >("cArg",  cfile);
//   argparse.get < std :: string >("dmArg", del);
//   if( !file_exists(cfile) )  error_file(cfile);
//   return;
// }


void parse_demo_tcm_realdata(int argc, char *argv[],
                             std :: string & cfile,
                             std :: string & bfile,
                             std :: string & tdfile,
                             int & nth,
                             std :: string & del
                            )
{
  ArgumentParser argparse("Demo " + std :: to_string(KMAP_MAJOR) + "." + std :: to_string(KMAP_MINOR) + "." + std :: to_string(KMAP_REVISION) + " for TCM models using real data ");

// #ifdef _OPENMP
//   nth  = omp_get_max_threads();
//   nth -= nth % 2;
//   argparse.add_argument < int >(          "tArg",  "t",  "threads", "Max number of threads exploitable",           false,  nth);
//
// #else
  nth = 1;

// #endif

  argparse.add_argument < std :: string >("cArg",  "c",  "cfile",   "TAC filename (with extension)",                true,  "");
  argparse.add_argument < std :: string >("bArg",  "b",  "bfile",   "Blood data filename (with extension)",         true,  "");
  argparse.add_argument < std :: string >("fArg", "f", "ffile",  "Framing protocol (with extension)",            true,  "");
  argparse.add_argument < std :: string >("dmArg", "dm", "delim",   "Delimiter (default: \"\\t\")",                 false, "\t");

  argparse.parse_args(argc, argv);

// #ifdef _OPENMP
//
//   argparse.get < int >(          "tArg",  nth);
//
// #endif
  argparse.get < std :: string >("cArg",  cfile);
  argparse.get < std :: string >("bArg",  bfile);
  argparse.get < std :: string >("fArg", tdfile);
  argparse.get < std :: string >("dmArg", del);

  if( !file_exists(cfile) )  error_file(cfile);
  if( !file_exists(bfile) )  error_file(bfile);
  if( !file_exists(tdfile) ) error_file(tdfile);

  return;
}
