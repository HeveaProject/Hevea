//************************************************************************
//                                                    CYL_output.h
//                                                  -------------------
//    written by : Vincent Borrelli, Saïd Jabrane, Francis Lazarus and Boris Thibert
//    date         : Aug 2009
//    email        : francis.lazarus@gipsa-lab.fr, boris.thibert@imag.fr
//************************************************************************

#ifndef _CYL_OUTPUT_H_
#define _CYL_OUTPUT_H_

//-------------------------------------------------------------
//          CYL_embedding<Num_type>::write_VTK_local()
//
// Writes a piece [0,nbre_x]*[0,nbre_y] of the embedding in a 
// VTK file format
//-------------------------------------------------------------
template < typename Num_type > void CYL_embedding <Num_type>::
write_VTK_local(string const & nomFichier, int nbre_x, int nbre_y) const
{
  if (nbre_x > n_x()) nbre_x = n_x();
 if (nbre_y > n_y() + 1) nbre_y = n_y() + 1;
  int npointsx(nbre_x);
  if ( nbre_x == n_x() ) npointsx++;
  int npointsy(nbre_y);

  std::ofstream fd(nomFichier.c_str());
  fd << "# vtk DataFile Version 3.0" << std::endl;
  fd << "plot Xuv" << std::endl;
  fd << "ASCII" << std::endl;
  fd << "DATASET STRUCTURED_GRID" <<  std::endl;
  fd << "DIMENSIONS " << npointsy << " " << npointsx << " "
     << 1 << std::endl;
  fd << "POINTS "<< npointsx * npointsy << " double" << std::endl;
  for (int i=0; i < nbre_x; i++) {
    for (int j=0; j < nbre_y; j++) {
      UTI_3point<Num_type> point = m_(i,j);
      fd << point.x() << "  " <<  point.y() << "  " <<  point.z()
	 <<  std::endl;
    }
  }
  // Case i=0
  if (nbre_x == n_x()) {
    for (int j=0;j < nbre_y;j++) {
      UTI_3point<Num_type> point = m_(0,j);
      fd << point.x() << "  " <<  point.y() << "  " <<  point.z()
	 << std::endl;
    }
  }
}

//-------------------------------------------------------------
//          CYL_embedding<Num_type>::write_VTK()
//
// Writes the embedding in a VTK file format
// We display one i-index out of step_x,and one j-index out of step_y
// so as to reduce the size.
//-------------------------------------------------------------
template < typename Num_type > void CYL_embedding <Num_type>::
write_VTK(string const & nomFichier, int step_x, int step_y) const
{
  int nx(n_x()), ny(n_y());
  // We add 1 to glue
  int npointsx = floor( ((Num_type)nx - 1.) / (Num_type)step_x ) + 2;
  int npointsy = floor( (Num_type)ny / (Num_type)step_y ) + 1;

  std::ofstream fd(nomFichier.c_str());
  fd << "# vtk DataFile Version 3.0" << std::endl;
  fd << "plot Xuv" << std::endl;
  fd << "ASCII" << std::endl;
  fd << "DATASET STRUCTURED_GRID" <<  std::endl;
  fd << "DIMENSIONS " << npointsy << " " << npointsx << " " << 1 << std::endl;
  fd << "POINTS "<< npointsy * npointsx << " double"
     << std::endl;
  for (int i=0; i < nx; i += step_x) {
    for (int j=0; j < ny + 1; j += step_y) {
      UTI_3point<Num_type> point = m_(i,j);
      fd << point.x() << "  " <<  point.y() << "  " <<  point.z()
	 << std::endl;
    }
  }
  // Case i=0
  for (int j=0; j < ny + 1; j += step_y) {
    UTI_3point<Num_type> point = m_(0,j);
    fd << point.x() << "  " <<  point.y() << "  " <<  point.z()
       << std::endl;
  }
}

//-------------------------------------------------------------
//          CYL_function_S1<Num_type>::write_VTK()
//
// Writes the integral curves in a VTK file
// The function is given by function*d1 + j*d2.
// We display one curve out of step_x, so as to reduce the size.
//-------------------------------------------------------------
template < typename Num_type > void CYL_function_S1 <Num_type>::
write_VTK(string const & nomFichier, 
	  UTI_2vector<Num_type> const & d1, 
	  UTI_2vector<Num_type> const & d2,
	  int step_x) const 
{
  int nx(n_x()), ny(n_y());

  int nb_lines = floor( ((Num_type)nx - 1.) / (Num_type)step_x ) + 1;
  int npointsy = ny + 1;

  std::ofstream fd(nomFichier.c_str());
  fd << "# vtk DataFile Version 3.0" << std::endl;
  fd << "Courbes integrales" << std::endl;
  fd << "ASCII" << std::endl;
  fd << "DATASET POLYDATA" <<  std::endl;
  fd << "POINTS "<< nb_lines * npointsy  << "  " << "float" << std::endl;
  for (int i=0; i< nx; i+= step_x) {
    for (int j=0; j <= ny; j++) {
      UTI_2vector<Num_type> const & vecteur =  
	              m_(i,j) * d1 + (j/((Num_type) ny)) *d2 ;
      fd << vecteur.x() <<"   "<<  vecteur.y()  << "  " << "0" << std::endl;
    }
  }
  fd << "LINES " << nb_lines  << " " << nb_lines * (npointsy + 1) << std::endl;
  for (int i=0; i < nb_lines; i++) {  
    fd << npointsy;
    for (int j = 0; j < npointsy; j++) {
      fd << " " << i * npointsy + j;
    }
    fd << std::endl;
  }
  fd.close();
}

#endif // _CYL_OUTPUT_H_
