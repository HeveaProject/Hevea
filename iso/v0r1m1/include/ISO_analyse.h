//************************************************************************
//                                                    ISO_analyse.h
//                                                  -------------------
//    written by:        : Vincent Borrelli, Saïd Jabrane, Francis Lazarus and Boris Thibert
//    email                : francis.lazarus@gipsa-lab.fr, boris.thibert@imag.fr
//************************************************************************
#ifndef _ISO_ANALYSE_H_
#define _ISO_ANALYSE_H_

#include <sstream>

//-------------------------------------------------------------
//  ISO_embedding<Num_type>::print_min_max_rho()
//
// Return the min and max values of rho_0,1,2 of
// Delta = (1 - delta_k)*f0^*<,> + delta_k*Id - f^*<,> 
//         = Sum_i rho_i li*li
//-------------------------------------------------------------
template <typename Num_type > 
Num_type ISO_embedding<Num_type>::print_min_max_rho(int k)
{
  int n = n_x_;
  
  Num_type max_rho[3];
  Num_type min_rho[3];

  UTI_sym_2form<Num_type> D0 = Delta(0., 0., k, f_);
  for (int l=0; l<3; l++ ) {
    max_rho[l] = min_rho[l] = rho_matrix_[l] * (UTI_3vector<Num_type>)D0;
  }
  int i;
#pragma omp parallel shared(max_rho, min_rho) private(i)
  {
    Num_type my_min_rho[3], my_max_rho[3];
    for (int l=0; l<3; l++ ) {
      my_min_rho[l] = min_rho[l];
      my_max_rho[l] = max_rho[l];
    }
#pragma omp for nowait
    for (i=1; i < n; i++) {
      Num_type x = (Num_type)i/(Num_type)n;
      for (int j=0; j < n; j++) {
	Num_type y = (Num_type)j/(Num_type)n;
      
	UTI_sym_2form<Num_type> D = Delta(x, y, k, f_);

	for (int l=0; l<3; l++ ) {
	  Num_type rho = rho_matrix_[l] * (UTI_3vector<Num_type>)D;
	  my_max_rho[l] = max( my_max_rho[l], rho );
	  my_min_rho[l] = min( my_min_rho[l], rho );
	}
      }
    }
#pragma omp critical
    {
      for (int l=0; l<3; l++ ) {
	max_rho[l] = max( my_max_rho[l], max_rho[l] );
	min_rho[l] = min( my_min_rho[l], min_rho[l] );
      }
    }
  }

  for (int l=0; l<3; l++ ) {
    cout << "rho_" << l << " min = " << min_rho[l]
	 << " max = " << max_rho[l] << endl; 
  }
  return min( min(min_rho[0],min_rho[1]), min_rho[2] );
}

//-------------------------------------------------------------
//  ISO_embedding<Num_type>::pull_back_o_phi()
//
// Computation of the metric F^*<,>_3 at the point phi(s,t). We have an access to 
// F o phi instead of F, thus we need to perform an adequate change of variables.
//-------------------------------------------------------------
template <typename Num_type > UTI_sym_2form<Num_type> 
ISO_embedding <Num_type>::
pull_back_o_phi( Num_type s,
		 Num_type t,
		 CYL_embedding<Num_type> & f_o_phi, 
		 CYL_function_S1<Num_type> & phi_X,
		 int dir ) const
{
  UTI_2vector<Num_type> const & V = V_[dir];
  UTI_2vector<Num_type> const & U = U_[dir];
  assert( (V*V)*(U*U)<1.000001 && (V*V)*(U*U)>.999999);

  Num_type dPhi_X_ds = phi_X.deriv_x(s,t);
  Num_type zeta = phi_X.deriv_y(s,t);

  UTI_sym_2form<Num_type> f_o_phi_pull_back = f_o_phi.pull_back(s,t);

  Num_type e = f_o_phi_pull_back.E();
  Num_type f = dPhi_X_ds*f_o_phi_pull_back.F() - zeta*e;
  Num_type g = zeta*zeta*e - 2*dPhi_X_ds*zeta*f_o_phi_pull_back.F() +
               dPhi_X_ds*dPhi_X_ds*f_o_phi_pull_back.G();

  Num_type V_2 = V*V;
  Num_type a_1x1x = V.x()*V.x() / (V_2*V_2);
  Num_type a_1x2x = V.x()*U.x();
  Num_type a_2x2x = U.x()*U.x()*V_2*V_2;

  UTI_sym_2form<Num_type> f_pull_back_o_phi;
  f_pull_back_o_phi.E() = e*a_1x1x + 2*f*a_1x2x + g*a_2x2x;

  Num_type a_1x1y = V.x()*V.y() / (V_2*V_2);
  Num_type a_2x1y = U.x()*V.y();
  Num_type a_1x2y = V.x()*U.y();
  Num_type a_2x2y = U.x()*U.y()*V_2*V_2;

  f_pull_back_o_phi.F() = e*a_1x1y + f*(a_2x1y + a_1x2y) + g*a_2x2y;

  Num_type a_1y1y = V.y()*V.y() / (V_2*V_2);
  Num_type a_1y2y = V.y()*U.y();
  Num_type a_2y2y = U.y()*U.y()*V_2*V_2;

  f_pull_back_o_phi.G() = e*a_1y1y + 2*f*a_1y2y + g*a_2y2y;

  return f_pull_back_o_phi / (dPhi_X_ds * dPhi_X_ds);
}

//-------------------------------------------------------------
//  ISO_embedding<Num_type>::print_iso_default()
//
// Print the isometric default = |mu_i - F_i^*<,>| = 
// |f_old_^*<,> + rho_dir(g_k - f_old_^*<,>)*prim[dir] - f_^*<,>|
// où |D|^2 = D_xx^2 + 2*D_xy^2 + D_yy^2.
//-------------------------------------------------------------
template <typename Num_type > Num_type ISO_embedding <Num_type>::
print_iso_default(int k, int dir ) const
{
  int nx = n_x_;
  int ny = n_y_;
  
  Num_type max_default = 0;
  Num_type mean_default = 0;
  Num_type max_default_iso = 0;
  Num_type mean_default_iso = 0;

  int i;
#pragma omp parallel shared(max_default, max_default_iso) private(i) 
  {
    Num_type my_max_default = 0, my_max_default_iso = 0;

#pragma omp for nowait reduction(+:mean_default, mean_default_iso)
    for (i=0; i < nx; i++) {
      Num_type x = (Num_type)i/(Num_type)nx;
      for (int j=0; j < ny; j++) {
	Num_type y = (Num_type)j/(Num_type)ny;

	Num_type rho = rho_matrix_[dir] * 
	  (UTI_3vector<Num_type>)Delta(x, y, k, f_old_);

	UTI_sym_2form<Num_type> prim_2_form = 
	  (UTI_3vector<Num_type>)prim_metric_[dir];


	UTI_sym_2form<Num_type> default_metric = f_old_.pull_back(x, y) + 
	  rho * prim_2_form - f_.pull_back(x, y);

	Num_type norm_default = sqrt( norm2_2form(default_metric) );
      
	my_max_default = max(norm_default, my_max_default);
	mean_default += norm_default;

	UTI_sym_2form<Num_type> default_metric_iso(1,0,1);
	default_metric_iso -= f_.pull_back(x, y);

	Num_type norm_default_iso = sqrt( norm2_2form(default_metric_iso) );
      
	my_max_default_iso = max(norm_default_iso, my_max_default_iso);
	mean_default_iso += norm_default_iso;
      }
    }
#pragma omp critical
    {
      max_default_iso = max( max_default_iso, my_max_default_iso );
      max_default = max( max_default, my_max_default );
    }
  }
  cout << "Max |mu_" << dir << " - F_" << dir << "^*<,>| = " << max_default 
       << " Mean |mu_" << dir << " - F_" << dir << "^*<,>| = " 
       << mean_default / (Num_type)(nx*ny) << endl;

  cout << "Max |I - F_" << dir << "^*<,>| = " << max_default_iso 
       << " Mean |I - F_" << dir << "^*<,>| = " 
       << mean_default_iso / (Num_type)(nx*ny)  << endl;
 
  return max_default;
}
  
//------------------------------------------------------------------------
//  ISO_embedding< Num_type > check_dt_apply_f_tilde()
//
// Check that || df_tilde/dt || = r o phi.
//-----------------------------------------------------------------------
template <typename Num_type>  void ISO_embedding<Num_type>::
check_dt_apply_f_tilde(CYL_embedding<Num_type> const & f_tilde,
		       TOR_embedding<Num_type> const & f_initial,
		       int const stage,
		       int const dir,                               
		       CYL_function_S1<Num_type> const & phi_X,
		       UTI_2vector<Num_type> const & V,
		       UTI_2vector<Num_type> const & U) const
{
  Num_type delta_max = 0;
  Num_type delta_mean = 0;
  
  int i;
#pragma omp parallel shared(delta_max) private(i) 
  {
    Num_type my_delta_max = 0;

#pragma omp for nowait reduction(+:delta_mean)
    for (i = 0; i < n_x_; i++) {
      for(int j=2; j < n_y_-1; j++) {
	UTI_3vector<Num_type> df_tilde_dt =
	  (Num_type)8.*(f_tilde.mat(i,j+1) - f_tilde.mat(i,j-1)) +
	  (f_tilde.mat(i,j-2) - f_tilde.mat(i,j+2));
	df_tilde_dt *= (Num_type) n_y_ / (Num_type)12.;
	
	UTI_2vector<Num_type> phi_ij = phi_X.mat(i,j)*V + 
	  ((Num_type)j/(Num_type)n_y_)*U;
	Num_type ray = radius(phi_ij, stage, dir, V, U, f_initial);
	Num_type delta = abs(df_tilde_dt.norm() - ray);
	my_delta_max = max( delta, my_delta_max );
	delta_mean += delta;
      }
    }
#pragma omp critical
    {
      delta_max = max( delta_max, my_delta_max );
    }
  }

  cout << "Max  | || df_tilde_dt || - r o phi | = " << delta_max << endl
       << "Mean | || df_tilde_dt || - r o phi | = "
       << delta_mean/ (Num_type) (n_x_*(n_y_ - 3))<< endl;
}


//-------------------------------------------------------------
//  ISO_embedding<Num_type>::check_pull_back_o_phi()
//
// Computation of the metric F^*<,>_3 at the point phi(s,t) by two 
// different ways. Either from F, or from F o phi, 
// in the second case, we need a change of variable... 
//-------------------------------------------------------------
template <typename Num_type > void 
ISO_embedding <Num_type>::
check_pull_back_o_phi( TOR_embedding<Num_type> & f,
		       CYL_embedding<Num_type> & f_o_phi, 
		       CYL_function_S1<Num_type> & phi_X,
		       int dir ) const
{
  for (int i = -1; i <= 1; i++) {
    Num_type s = i*.04;
    cout << "phi_X(" << s << ", 0) = " << phi_X(s, 0) 
	 << " phi_X(" << s << ", 1) = " << phi_X(s, 1) << endl;

    for (int j = 0; j < 3; j++) {
      Num_type t = (Num_type) (j*10 + 1) / 30.;
      UTI_2vector<Num_type> phi_s_t = phi_X(s,t)*V_[dir] + t*U_[dir];
      
      cout << "f*<>" << phi_s_t << " = " << f.pull_back(phi_s_t) << endl
	   << "f*<> en phi(" << s << ", " << t << ")  = " 
	   << pull_back_o_phi(s, t, f_o_phi, phi_X, dir) << endl << endl;
    }
  }
}

#endif // _ISO_ANALYSE_H_
