//************************************************************************
//                                               ISO_implementation.h
//                                                  -------------------
//    written by:        : Vincent Borrelli, Sa�d Jabrane, Francis Lazarus and Boris Thibert
//    email                : francis.lazarus@gipsa-lab.fr, boris.thibert@imag.fr
//************************************************************************
#ifndef _ISO_IMPLEMENTATION_H_
#define _ISO_IMPLEMENTATION_H_

#include <sstream>

//-------------------------------------------------------------
//          ISO_embedding<Num_type>:: Delta()
//
// Computation of the metric Delta at the point (x,y) i.e.
// (1 - delta_k)*f0^*<,> + delta_k*Id - f^*<,> 
// 
//-------------------------------------------------------------
template <typename Num_type > UTI_sym_2form<Num_type>
ISO_embedding <Num_type>::
Delta(Num_type x,
        Num_type y,
        const int k, 
        TOR_embedding<Num_type> const & f ) const
{
  UTI_sym_2form<Num_type> Id(1.,0.,1.); // Euclidean metric

  Num_type delta_k = coef_delta(k);

  return (1. - delta_k)*f_0_.pull_back(x,y) + delta_k*Id - f.pull_back(x,y);
}

//-------------------------------------------------------------
//          ISO_embedding<Num_type>:: zeta()
//
// zeta = - (V_perp(f) * V(f)) / (V(f) * V(f))
//-------------------------------------------------------------
template <typename Num_type > Num_type
ISO_embedding <Num_type >::
zeta(Num_type x,
      Num_type y,
      UTI_2vector<Num_type> const & V,
      UTI_2vector<Num_type> const & U,
      TOR_embedding<Num_type> const & f ) const
{
  UTI_3vector<Num_type> V_f = f.deriv(V,x,y);

  return - (f.deriv(U,x,y) * V_f) / (V_f*V_f);
}

//-------------------------------------------------------------
//          ISO_embedding<Num_type>:: W()
//
// W = U +  zeta *  V
//-------------------------------------------------------------
template <typename Num_type > UTI_2vector<Num_type>
ISO_embedding <Num_type >::
W(Num_type x,
    Num_type y,
    UTI_2vector<Num_type> const & V,
    UTI_2vector<Num_type> const & U,
    TOR_embedding<Num_type> const & f ) const
{
  return U +  zeta(x,y,V,U,f) *  V;
}

//-------------------------------------------------------------
//          ISO_embedding<Num_type>:: radius()
//
// 
//-------------------------------------------------------------
template <typename Num_type > Num_type
ISO_embedding<Num_type>::
radius(Num_type x,
      Num_type y,
      int const stage,
      int const dir,
      UTI_2vector<Num_type> const & V,
      UTI_2vector<Num_type> const & U,
      TOR_embedding<Num_type> const & f ) const
{
  UTI_3vector<Num_type> z = f.deriv(W(x, y, V, U, f), x, y);
 
  // radius = sqrt( r = z*z + rho_i/ V*V )
  return sqrt( z*z + (rho_matrix_[dir] * 
		    (UTI_3vector<Num_type>)Delta(x,y,stage,f)) / (V * V) );
}

//-------------------------------------------------------------
//  ISO_embedding<Num_type>::ISO_embedding()
//
// Constructor.
//-------------------------------------------------------------
template < typename Num_type >  ISO_embedding<Num_type>::
ISO_embedding( TOR_embed_type e_type, 
		int matrix_size,
		Num_type r1,
		Num_type r2,
		int ol,
		int dl ) : 
  f_0_(matrix_size, matrix_size, e_type, r1, r2),
  f_(matrix_size, matrix_size, e_type, r1, r2), 
  f_old_(f_),
  n_x_(matrix_size),
  n_y_(matrix_size),
  r1_(r1),
  r2_(r2),
  stage_(0),
  output_level_(ol),
  debug_level_(dl)
{
  init_();
}

//-------------------------------------------------------------
//  ISO_embedding<Num_type>::ISO_embedding()
//
// Constructor.
//-------------------------------------------------------------
template < typename Num_type >  ISO_embedding<Num_type>::
ISO_embedding( string const & fichier_f0,
		string const & fichier_f,
		int stage,
		int ol,
		int dl ) : 
  f_0_( TOR_embedding<Num_type>::read_binary( fichier_f0 ) ),
  f_( TOR_embedding<Num_type>::read_binary( fichier_f ) ),
  f_old_(f_),
  n_x_(f_.n_x()),
  n_y_(f_.n_y()),
  r1_(0),
  r2_(0),
  stage_(stage),
  output_level_(ol),
  debug_level_(dl)
{
  init_();
}

//-------------------------------------------------------------
//  ISO_embedding<Num_type>::init_()
//
//-------------------------------------------------------------
template < typename Num_type >  void ISO_embedding<Num_type>::
init_()
{
  UTI_2vector<Num_type> l[3];   // The 3 constant 1-forms l1, l2, l3.
  l[0] = UTI_2vector<Num_type>(1,0);
  l[1] = UTI_2vector<Num_type>(1/sqrt(5), 2/sqrt(5));
  l[2] = UTI_2vector<Num_type>(1/sqrt(5), -2/sqrt(5));

  V_[0] = UTI_2vector<Num_type>(0,1); 
  V_[1] = UTI_2vector<Num_type>(-2,1);
  V_[2] = UTI_2vector<Num_type>(2,1);
  
  UTI_2vector<Num_type> bezout[3]; // associated Bezout coefficients.
  bezout[0] = UTI_2vector<Num_type>(0,1);
  bezout[1] = UTI_2vector<Num_type>(-1,-1);
  bezout[2] = UTI_2vector<Num_type>(1,-1);

  for (int dir = 0; dir < 3; dir++) {
    assert( l[dir] * V_[dir] == 0 ); // V_[dir] is in Ker l_dir
    assert( bezout[dir] * V_[dir] == 1 ); 
  }

  // Intitialize the primitive metrics for the metric decompositions
  for ( int i=0; i < 3; i++ )
    prim_metric_[i].set_from_1form(l[i]); 

  // Computation of the "inversion matrix" for the decompositions
  Num_type det = prim_metric_[0][0]*( prim_metric_[1][1]*prim_metric_[2][2] - 
				      prim_metric_[1][2]*prim_metric_[2][1] ) -
                 prim_metric_[0][1]*( prim_metric_[1][0]*prim_metric_[2][2] - 
				      prim_metric_[1][2]*prim_metric_[2][0] ) +
                 prim_metric_[0][2]*( prim_metric_[1][0]*prim_metric_[2][1] - 
				      prim_metric_[1][1]*prim_metric_[2][0] );

  for (int i=0; i < 3; i++) {
    for (int j=0; j < 3; j++) {
      rho_matrix_[i][j] = 
	prim_metric_[(i+2)%3][(j+2)%3] * prim_metric_[(i+1)%3][(j+1)%3] -
	prim_metric_[(i+2)%3][(j+1)%3] * prim_metric_[(i+1)%3][(j+2)%3];
    }
    rho_matrix_[i] /= det;
  }

  // parameters in the direction Ker(l_dir)
  for (int dir = 0; dir < 3; dir++) {
    U_[dir] = UTI_2vector<Num_type>(-V_[dir].y(), V_[dir].x());
    U_[dir] /=  V_[dir]*V_[dir];
    bezout_perp_[dir] = UTI_2vector<Num_type>(-bezout[dir].y(), bezout[dir].x());
    lambda_[dir] = bezout[dir] * U_[dir];
  }
}

//-------------------------------------------------------------
//  ISO_embedding<Num_type>::corrugation()
//
//-------------------------------------------------------------
template <typename Num_type> void ISO_embedding <Num_type>::
corrugation(int stage, int dir, int n_oscil) 
{
  stage_ = stage;

  f_old_ = f_;

  if (output_level_ > 3) { // Display of W and of the torus before the corrugation with the image field
    TOR_2vector_field<Num_type> W_field(n_x_, n_y_);
    TOR_3vector_field<Num_type> z(n_x_, n_y_);
    for (int i = 0; i< n_x_; i++) {
      Num_type x_i = (Num_type)i/(Num_type)n_x_;
      for (int j = 0; j< n_y_; j++) {
	Num_type y_j = (Num_type)j/(Num_type)n_y_;
	W_field.mat(i,j) = W( x_i, y_j, V_[dir], U_[dir], f_ );
	z.mat(i,j) = f_.deriv(W_field.mat(i,j), x_i, y_j);
      }
    }
    W_field.write_VTK(generate_file_name("W_field", "vtk", 
				       stage, dir, n_oscil));
    f_.write_VTK(generate_file_name("torus_with_W", "vtk", 
					stage, dir, n_oscil), z);
  }
  // Computation of the flow (in fact its first coordinate in the basis 
  // (V[dir],U[dir])).
  CYL_function_S1< Num_type > phi_X = flow( f_, V_[dir], U_[dir] ); 

  if (output_level_ > 2) 
    phi_X.write_VTK(generate_file_name("W_Flow", "vtk", stage, dir, n_oscil), 
		    V_[dir], U_[dir], 100);

  CYL_embedding< Num_type > f_tilde = integrate_h_hai( phi_X, // _hai or _trap
						       f_, 
						       V_[dir], 
						       U_[dir],
						       stage,
						       dir,
						       n_oscil ); 

  // FOR DEBUG
  if (debug_level_ > 1)
    check_dt_apply_f_tilde(f_tilde, f_, stage, dir, phi_X,
			   V_[dir], U_[dir]);
  // END FOR DEBUG

  if (output_level_ > 1) {
    f_tilde.write_VTK(generate_file_name("unglued_torus", "vtk", // change 10 to other 
			      stage, dir, n_oscil, n_x_/10), 10, 10);  // value to modify the sampling.
                                                                                   
    f_tilde.write_VTK_local(generate_file_name("local_unglued_torus", "vtk", 
						 stage, dir, n_oscil), 200, 1000); // Subgrid = 200*1000.
  }

  f_tilde.gluing( f_, phi_X, V_[dir], U_[dir] );

  f_ = f_tilde.cyl_to_torus( phi_X, lambda_[dir], V_[dir], bezout_perp_[dir], f_ );

  if (debug_level_ > 0)
    check_pull_back_o_phi( f_, f_tilde, phi_X, dir );

    f_.write_VTK(generate_file_name("sampled_torus","vtk", 
				    stage, dir, n_oscil, n_x_/10), 10, 10);
    f_.write_VTK_local(generate_file_name("local_torus","vtk", 
					    stage, dir, n_oscil), 500, 500);
  if (output_level_ > 0) {
    f_.write_VTK(generate_file_name("torus","vtk", stage, dir, n_oscil));
    // The next command outputs a k*k grid VRML file with k = min(n_x_, 3500). 
    // If you want to use a 3D printer for the third corrugated torus it seems accurate to  
    // use the default grid size (10,000*10,1000) for the computations and an  
    // undersampled grid of size 3,500*3,500 for the ouput is sufficient (larger grids may 
    // not be accepted by 3D printers). 
    f_.write_VRML(generate_file_name("torus","wrl", stage, dir, n_oscil, min(n_x_, 3500 ) ), 
		         3500 );
    //    f_.write_OFF(generate_file_name("torus","off", stage, dir, n_oscil));
    //    f_.write_Povray(generate_file_name("torus","pov", stage, dir, n_oscil));
  }
}

//-------------------------------------------------------------
//  ISO_embedding<Num_type>::resample()
//
// Re-samples f_ with nx' = floor(nx*ratio)
// ny' = floor(ny*ratio). The x and y dimension of the 
// grid may differ!
//-------------------------------------------------------------
template <typename Num_type> void ISO_embedding <Num_type>::
resample( Num_type ratio_x, Num_type ratio_y )
{
  f_.resample(ratio_x, ratio_y);
  n_x_ = f_.n_x();
  n_y_ = f_.n_y();
}

//-------------------------------------------------------------
//  ISO_embedding<Num_type>::generate_file_name()
//
// Generates an indexed file name with parameters
//-------------------------------------------------------------
template <typename Num_type > string  ISO_embedding <Num_type>::
generate_file_name( const char* nom_fichier,
		             const char* extension,
		             int stage,
		             int dir,
		             int n_theta,
		             int nx ) const
{
  ostringstream nomFinal;
  nomFinal << nom_fichier << "_r=" << r1_ << "_" << r2_ << "_k=" << stage
	   << "_dir=" << dir << "_Npts=" << nx << "_Nosc=" 
	   << n_theta << '.' << extension; 
  return nomFinal.str();
}

//-------------------------------------------------------------
//  ISO_embedding<Num_type>::generate_file_name()
//
// Generates an indexed file name with parameters
//-------------------------------------------------------------
template <typename Num_type > string  ISO_embedding <Num_type>::
generate_file_name( const char* nom_fichier,
		             const char* extension,
		             int stage,
		             int dir,
		             int n_theta ) const
{
  return generate_file_name( nom_fichier, extension, stage, dir, n_theta, n_x_ );
}

#endif // _ISO_IMPLEMENTATION_H_
