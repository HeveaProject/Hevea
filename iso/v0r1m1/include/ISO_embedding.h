//************************************************************************
//                                                 ISO_embedding.h
//                                                  -------------------
//    written by : Vincent Borrelli, Saïd Jabrane, Francis Lazarus and Boris Thibert
//    email        : francis.lazarus@gipsa-lab.fr, boris.thibert@imag.fr
//************************************************************************
#ifndef _ISO_EMBEDDING_H_
#define _ISO_EMBEDDING_H_

#include "NonStiffIntegratorT.h"
#include "IntegratorT.h"
#include "NonStiffIntegratorT.cpp"
#include "IntegratorT.cpp"

//#include "LOC_param.h"


//----------------------------------------------------------//
//                  class ISO_embedding
//----------------------------------------------------------//
template <typename Num_type > class ISO_embedding 
{
  TOR_embedding<Num_type> f_0_;   // contains the initial embedding.
  TOR_embedding<Num_type> f_;     // contains the current embedding.
  TOR_embedding<Num_type> f_old_; // contains the previous embedding.
  UTI_2vector<Num_type> V_[3], U_[3], bezout_perp_[3];
  Num_type lambda_[3];
  TOR_primitive<Num_type> prim_metric_[3]; // the 3 primitive metrics for 
                                           // the decomposition of the metrics.
  UTI_3vector<Num_type> rho_matrix_[3];  // The matrix for the calculus of 
                                           // these decompositions.
  int stage_;            // the k in delta_k.
  int n_x_, n_y_;        // Dimensions of matrices. 
  Num_type const r1_;    // Small radius of the torus.
  Num_type const r2_;    // Large radius of the torus.
  int output_level_;     // 0 : subsampled embeddings and local piece in VTK, 
                         // 1 : + N*N torus in VTK and VRML with N=min(n, 3500)
                         // 2 : + Torus before gluing and before reparametrization
                         // 3 : + Integral curves in the flat square
                         // 4 : + W Field in the flat square and its embedding. 
  int debug_level_;      // 0 : From 0: no debug output to 3: more outputs. 

  friend class ISO_h_s; // So that ISO_h_s has an access to the fonctions of ISO_..
  friend class ISO_last_h_s; // So that ISO_last_h_s has an access to the fonctions of ISO_..

  void init_();

  UTI_sym_2form<Num_type> Delta(Num_type x,
				Num_type y,
				const int k, 
				TOR_embedding<Num_type> const & f ) const;
    
  Num_type zeta(Num_type x,
		 Num_type y,
		 UTI_2vector<Num_type> const & V,
		 UTI_2vector<Num_type> const & U,
		 TOR_embedding<Num_type> const & f ) const;

  UTI_2vector<Num_type> W(Num_type x,
			    Num_type y,
			    UTI_2vector<Num_type> const & V,
			    UTI_2vector<Num_type> const & U,
			    TOR_embedding<Num_type> const & f ) const;

  UTI_2vector<Num_type> W(UTI_2vector<Num_type> const & param,
			    UTI_2vector<Num_type> const & V,
			    UTI_2vector<Num_type> const & U,
			    TOR_embedding<Num_type> const & f ) const
  { return W(param.x(), param.y(), V, U, f); }

  Num_type radius(Num_type x,
		 Num_type y,
		 int const stage,
		 int const dir,
		 UTI_2vector<Num_type> const & V,
		 UTI_2vector<Num_type> const & U,
		 TOR_embedding<Num_type> const & f ) const;

  Num_type radius(UTI_2vector<Num_type> const & param,
		 int const stage,
		 int const dir,
		 UTI_2vector<Num_type> const & V,
		 UTI_2vector<Num_type> const & U,
		 TOR_embedding<Num_type> const & f ) const
  { return radius(param.x(), param.y(), stage, dir, V, U, f); }

  CYL_function_S1<Num_type> flow(TOR_embedding<Num_type> const & f,
				 UTI_2vector<Num_type> const & V,
				 UTI_2vector<Num_type> const & U
				 ) const; 

  void integrate_func( HAI_function & f,
		       int const dim, // dimension of problem
		       Num_type *init_cond,    
		       std::vector< std::vector<Num_type> > & flot,
		       int Ny ) const;

  CYL_embedding<Num_type> integrate_h_hai(
				 CYL_function_S1< Num_type > const & phi_X,
				 TOR_embedding<Num_type> const & f,
				 UTI_2vector<Num_type> const & V,
				 UTI_2vector<Num_type> const & U,
				 int const stage,
				 int const dir,
				 int const N_oscil) const;
 
  CYL_embedding<Num_type> integrate_h_trap(
				 CYL_function_S1< Num_type > const & phi_X,
				 TOR_embedding<Num_type> const & f,
				 UTI_2vector<Num_type> const & V,
				 UTI_2vector<Num_type> const & U,
				 int const stage,
				 int const dir,
				 int const N_oscil
				 ) const;

  string generate_file_name( const char* nom_fichier,
			       const char* extension,
			       int stage,
			       int dir,
			       int n_oscil ) const;

  string generate_file_name( const char* nom_fichier,
			       const char* extension,
			       int stage,
			       int dir,
			       int n_oscil,
			       int nx ) const;

  void check_dt_apply_f_tilde(CYL_embedding<Num_type> const & f_tilde,
			      TOR_embedding<Num_type> const & f_initial,
			      int const stage,
			      int const dir,
			      CYL_function_S1<Num_type> const & phi_X,
			      UTI_2vector<Num_type> const & V,
			      UTI_2vector<Num_type> const & U) const;
  
  void check_pull_back_o_phi( TOR_embedding<Num_type> & f,
			      CYL_embedding<Num_type> & f_o_phi, 
			      CYL_function_S1<Num_type> & phi_X,
			      int dir ) const;

  UTI_sym_2form<Num_type> pull_back_o_phi( Num_type s,
					   Num_type t,
					   CYL_embedding<Num_type> & f_o_phi, 
					   CYL_function_S1<Num_type> & phi_X,
					   int dir ) const;

/*   UTI_sym_2form<Num_type> pull_back_o_phi( Num_type s, */
/* 					   Num_type t, */
/* 					   LOC_embedding<Num_type> & f_o_phi,  */
/* 					   LOC_function<Num_type> & phi_X, */
/* 					   int dir ) const; */



  void print_lengths(std::vector<std::vector<UTI_3point<Num_type> > > & curves,
		     Num_type param_length ) const;
public:
  ISO_embedding( TOR_embed_type e_type, 
		  int matrix_size,
		  Num_type r1 = 0,
		  Num_type r2 = 0,
		  int output_level = 0, 
		  int debug_level = 0 );

  ISO_embedding( string const & fichier_f0,
		  string const & fichier_f,
		  int stage,
		  int output_level = 0, 
		  int debug_level = 0 );

  Num_type coef_delta(int k=0) const {    //calculus of delta_k
    if (k == 0) return 4* (r2_*r2_ - r1_*r1_)/ (3 + r2_*r2_ - 3*r1_*r1_);
    else if (k == 1) return 0.8; 
    else if (k == 2) return 0.95; 
    else if (k > 2) return (Num_type) (k*k - 1) / (Num_type) (k*k);
    else return 1;
  }  

  void corrugation(int stage, 
		   int dir, 
		   int n_oscil);


  void resample( Num_type ratio_x, Num_type ratio_y);
  void reset_f0() { f_0_ = f_; }

  Num_type print_min_max_rho(int stage = 0); 
  Num_type print_iso_default(int k, int dir) const;

  void read_binary_f0( string const & nomFichier )  {
    f_0_ = TOR_embedding<Num_type>::read_binary( nomFichier );
  }
  void read_binary_f( string const & nomFichier )  {
    f_ = TOR_embedding<Num_type>::read_binary( nomFichier );
  }
  void write_binary_f0( const char* nomFichier )  {
    f_0_.write_binary( generate_file_name( nomFichier, 
					   "bin", 
					   0, 
					   0, 
					   0 ) );
  }
  void write_binary_f( const char* nomFichier, int dir, int n_oscil )  {
    f_.write_binary( generate_file_name( nomFichier, 
					   "bin", 
					   stage_, 
					   dir, 
					   n_oscil ) );
  }
};


#include "ISO_integrate.h"
#include "ISO_implementation.h"
#include "ISO_analyse.h"

#endif // _ISO_EMBEDDING_H_
