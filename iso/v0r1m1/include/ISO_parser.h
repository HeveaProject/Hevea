//************************************************************************
//                                                    ISO_parser.h
//                                                  -------------------
//    written by:        : Vincent Borrelli, Saïd Jabrane, Francis Lazarus and Boris Thibert
//    email                : francis.lazarus@gipsa-lab.fr, boris.thibert@imag.fr
//************************************************************************
#include <iostream>
#include<fstream> 
#include "argstream.h"
#include "UTI_utils.h"
#include <vector>

using namespace std;

template<typename NT> int parse_command_line_iso( int argc, 
				      char* argv[],  
				      int & output_level,
				      int & debug_level,
				      TOR_embed_type & embedding_type,
				      int & nx, 
				      int * n_oscil,
				      NT & r1,
				      NT & r2 )
{
  bool is_standard_bis;
  bool is_no_output;

  argstream as(argc,argv);  
  as >> parameter('n', "nb-points", nx, "Number of x/y samples ", false)
     >> option('b', "standard-bis", is_standard_bis, 
	       "If true use yx parameterization\n                            of standard torus." )
     >> option("no-output", is_no_output, "No output file")
     >> parameter('r', "r1", r1, "Small torus radius. ", false)
     >> parameter('R', "r2", r2, "Large torus radius. ", false)
     >> parameter('o', "Output-level", output_level,
		  "-- 0: subsampled embeddings and local piece in VTK\n -- 1: + N*N torus in VTK and VRML with N=min(n, 3500)\n -- 2: + Torus before gluing and before reparametrization\n -- 3: + Integral curves in the flat square\n -- 4: + W Field in the flat square and its embedding. ", false )
     >> parameter('d', "Debug-level", debug_level,
		  "From 0: no debug output to 3: more outputs. ", false )
     >> parameter('L', "N0", n_oscil[0], "N_oscil0. ", false)
     >> parameter('M', "N1", n_oscil[1], "N_oscil1. ", false)
     >> parameter('N', "N2", n_oscil[2], "N_oscil2. ", false)
     >> help();

  if (as.helpRequested())
    {
      cout << as.usage() << endl;
      exit(-1);
    }

  if (nx < 0 ) {
    cout << "You must specify a positive nx.\n Type -h for help" << endl;
    exit(-1);
  }
    if (is_no_output) output_level = 0;
  embedding_type = TOR_STANDARD;
  if (is_standard_bis) embedding_type = TOR_STANDARD_BIS;

  // Write the parameters on std::cout.
  cout << "----------------------------------"
       << "-----------------------------------" << endl;  
  if (embedding_type == TOR_STANDARD) cout << "Uses TOR_STANDARD ";
  else cout << "Uses TOR_STANDARD_BIS ";
  cout << " nx = " << nx 
       << " Oscillation Numbers = " << n_oscil[0] << " " << n_oscil[1] << " " << n_oscil[2] 
      << endl << "r1 = " << r1 << " r2 = " << r2 << " output level = " << output_level
       << " debug level = " << debug_level << endl << endl;

  return 0;
}
