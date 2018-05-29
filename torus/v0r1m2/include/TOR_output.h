//************************************************************************
//                                                    TOR_output.h
//                                                  -------------------
//    written by : Vincent Borrelli, Saïd Jabrane, Francis Lazarus and Boris Thibert
//    date         : Oct 2009
//    email        : francis.lazarus@gipsa-lab.fr, boris.thibert@imag.fr
//************************************************************************

// =======================================================
//
// This file contains the implementation of the class members of 
// various TOR_ embedding, vector fields,... for ouput such as VTK files.
//
// =======================================================

#ifndef _TOR_OUTPUT_H_
#define _TOR_OUTPUT_H_

#include <TOR_torus.h>


//-------------------------------------------------------------
//          TOR_2vector_field<Num_type>::write_VTK()
//
// Writes the vector field encoded at each vertex
//  in a VTK file format
//-------------------------------------------------------------
template < typename Num_type > void TOR_2vector_field<Num_type>::
write_VTK(string const & nomFichier) const
{
  std::ofstream fd(nomFichier.c_str());
  fd << "# vtk DataFile Version 3.0" << std::endl;
  fd << "Champ de vecteur dans le carre" << std::endl;
  fd << "ASCII" << std::endl;
  fd << "DATASET STRUCTURED_GRID" <<  std::endl;
  fd << "DIMENSIONS " << n_y()+1 << " " << n_x()+1 << " "
     << 1 << std::endl;
  fd << "POINTS "<< (n_x()+1) * (n_y()+1) << "  " << "float"
     << std::endl;
  for (int i=0; i< n_x(); i++) {
    for (int j=0; j< n_y(); j++) {
      fd << i/((double) n_x()) <<"   "<<  j/((double) n_x()) 
	 << "   " <<  0 << std::endl;
    }
    fd << i/((double) n_x()) <<"   "<<  0 << "   " << 0 <<  std::endl;
  }
  //Case i=0
  for (int j=0; j< n_y(); j++) {
    fd << 0 <<"   "<<  j/((double) n_x()) << "   " <<  0 << std::endl;
  }
  fd << 0 <<"   "<<  0 << "   " << 0 <<  std::endl;

  //The normals
  fd << "POINT_DATA "<< (n_x()+1) * (n_y()+1) <<  std::endl;
  fd << "VECTORS ChampDeVecteur float " << std::endl;
  for (int i=0; i< n_x(); i++) {
    for (int j=0; j< n_y(); j++) {
      fd << m_(i,j).x() <<"   "<<  m_(i,j).y() << "   " <<  0 
	 <<  std::endl;
    }
    fd << m_(i,0).x() <<"   "<<  m_(i,0).y() << "   " <<  0 <<  
      std::endl;
  }
  //Case i=0
  for (int j=0; j< n_y(); j++) {
    fd << m_(0,j).x() <<"   "<<  m_(0,j).y() << "   " <<  0
       <<  std::endl;
  }
  fd << m_(0,0).x() <<"   "<< m_(0,0).y() << "   " <<  0 
     <<  std::endl;
  fd.close();
}

//-------------------------------------------------------------
//          TOR_embedding<Num_type>::write_VTK()
//
// Writes the embedding in a VTK file format. We write
// one point over step_x_y so as to reduce the file's size.
//-------------------------------------------------------------
template < typename Num_type > void TOR_embedding <Num_type>::
write_VTK(string const & nomFichier, int step_x, int step_y) const
{
  int nx(n_x()), ny(n_y());
  // We add 1 to glue
  int npointsx = floor( ((Num_type)nx - 1.) / (Num_type)step_x ) + 2;
  int npointsy = floor( ((Num_type)ny - 1.) / (Num_type)step_y ) + 2;

  std::ofstream fd(nomFichier.c_str());//,ofstream::binary);
  fd << "# vtk DataFile Version 3.0" << std::endl;
  fd << "Torus" << std::endl;
  fd << "ASCII" << std::endl;
  fd << "DATASET STRUCTURED_GRID" <<  std::endl;
  fd << "DIMENSIONS " << npointsy << " " <<  npointsx << " " << 1 << std::endl;
  fd << "POINTS "<<  npointsy * npointsx << " double" << std::endl;
  //FILE *fd; //for binary case
  //fprintf(fd, "# vtk DataFile Version 3.0\n");
  //fprintf(fd, "Torus\n");
  //fprintf(fd, "BINARY\n");
  //fprintf(fd, "DATASET STRUCTURED_GRID\n");
  //fprintf(fd, "DIMENSIONS %d %d %d\n", npointsy, npointsx, 1);
  //fprintf(fd, "POINTS %d double\n", npointsy*npointsx);
  for (int i=0; i< nx; i += step_x) {
    for (int j=0; j< ny; j += step_y) {
      UTI_3point<Num_type> const & point = m_(i,j);
      fd << point.x() <<"   "<<  point.y() << "   " <<  point.z() << std::endl;
      //double x=point.x(),y=point.y(),z=point.z();
      //fwrite(&x, sizeof(double), 1, fd);
      //fwrite(&y, sizeof(double), 1, fd);
      //fwrite(&z, sizeof(double), 1, fd);
    }
    // Case j=0
    UTI_3point<Num_type> const & point = m_(i,0);
    fd << point.x() << "  " <<  point.y() << "  " <<  point.z() << std::endl;
  }
  // Case i=0
  for (int j=0; j < ny; j += step_y) {
    UTI_3point<Num_type> const & point = m_(0,j);
    fd << point.x() << "  " <<  point.y() << "  " <<  point.z() << std::endl;
  }
  // Casee i=j=0
  UTI_3point<Num_type> const & point = m_(0,0);
  fd << point.x() << "  " <<  point.y() << "  " <<  point.z() << std::endl;

  //fclose(fd);
}

//-------------------------------------------------------------
//          TOR_embedding<Num_type>::write_VTK_local()
//
// Ecrit le plongement dans un fichier au format VTK
//-------------------------------------------------------------
template < typename Num_type > void TOR_embedding <Num_type>::
write_VTK_local(string const & nomFichier, int nbre_x, int nbre_y) const
{
  // When we glue, we add 1
  if (nbre_x > n_x()) nbre_x = n_x();
  if (nbre_y > n_y()) nbre_y = n_y();
  int npointsx(nbre_x);
  if ( nbre_x == n_x() ) npointsx++;
  int npointsy(nbre_y);
  if ( nbre_y == n_y() ) npointsy++;

  std::ofstream fd(nomFichier.c_str());//,ofstream::binary);
  fd << "# vtk DataFile Version 3.0" << std::endl;
  fd << "Torus" << std::endl;
  fd << "ASCII" << std::endl;
  fd << "DATASET STRUCTURED_GRID" <<  std::endl;
  fd << "DIMENSIONS " << npointsy << " " <<  npointsx << " " << 1 << std::endl;
  fd << "POINTS "<<  npointsy * npointsx << " double" << std::endl;
  for (int i=0; i< nbre_x; i++) {
    for (int j=0; j< nbre_y; j++) {
      UTI_3point<Num_type> const & point = m_(i,j);
      fd << point.x() <<"   "<<  point.y() << "   " <<  point.z() << std::endl;
    }
    if ( nbre_y == n_y() ) {//Eventually add j=0
      UTI_3point<Num_type> const & point = m_(i,0);
      fd << point.x() <<"   "<<  point.y() << "   " <<  point.z() << std::endl;
    }
  }
  //Eventually add i=0
  if ( nbre_x == n_x() ) {
    for (int j=0; j< nbre_y; j++) {
      UTI_3point<Num_type> const & point = m_(0,j);
      fd << point.x() <<"   "<<  point.y() << "   " <<  point.z() << std::endl;
    }
    if ( nbre_y == n_y() ) { // Eventually add j=0
      UTI_3point<Num_type> const & point = m_(0,0);
      fd << point.x() <<"   "<<  point.y() << "   " <<  point.z() << std::endl;
    }
  }
  // fclose(fd);
}

//-------------------------------------------------------------
//          TOR_embedding<Num_type>::write_VTK()
//
// Writes the embedding in a VTK file
// Also writes a vector field encoded at each vertex
//-------------------------------------------------------------
template < typename Num_type > void TOR_embedding <Num_type>::
write_VTK(string const & nomFichier, 
	  TOR_3vector_field<Num_type> const & champ, 
	  bool reduit, int nbre_x, int nbre_y) const
{
  // When we glue, we add 1
  int npointsx, npointsy;
  if (! reduit) {
    nbre_x=n_x();
    nbre_y=n_y();
  }
  if ( nbre_x == n_x() ) npointsx = nbre_x +1; else npointsx = nbre_x;
  if ( nbre_y == n_y() ) npointsy = nbre_y +1; else npointsy = nbre_y;

  std::ofstream fd(nomFichier.c_str()); //,ofstream::binary);
  fd << "# vtk DataFile Version 3.0" << std::endl;
  fd << "Torus" << std::endl;
  fd << "ASCII" << std::endl;//MODIFIER
  fd << "DATASET STRUCTURED_GRID" <<  std::endl;
  fd << "DIMENSIONS " << npointsy << " " <<  npointsx << " " << 1 << std::endl;
  fd << "POINTS "<<  npointsy * npointsx << "  " << "double" << std::endl;
  for (int i=0; i< nbre_x; i++) {
    for (int j=0; j< nbre_y; j++) {
      UTI_3point<Num_type> const & point = m_(i,j);
      fd << point.x() <<"   "<<  point.y() << "   " << point.z() << std::endl;
    }
    if ( nbre_y == n_y() ) {
      UTI_3point<Num_type> const & point = m_(i,0);
      fd << point.x() << "   " << point.y() << "   " << point.z() << std::endl;
    }
  }
  if ( nbre_x == n_x() ) {
    for (int j=0; j< nbre_y; j++) {
      UTI_3point<Num_type> const & point = m_(0,j);
      fd << point.x() << "   " << point.y() << "   " << point.z() << std::endl;
    }
    if ( nbre_y == n_y() ) {
      UTI_3point<Num_type> const & point = m_(0,0);
      fd << point.x() << "   " << point.y() << "   " << point.z() << std::endl;
    }
  }

  //THE VECTOR FIELD
  fd << "POINT_DATA "<<  npointsy *  npointsx <<  std::endl;
  fd << "VECTORS ChampDeVecteur double " << std::endl;
  for (int i=0; i< nbre_x; i++) {
    for (int j=0; j< nbre_y; j++) {
      UTI_3vector<Num_type> const & vecteur = champ.mat(i,j);
      fd << vecteur.x() <<"   "<<  vecteur.y() << "   " <<  vecteur.z() 
	 <<  std::endl;
    }
    if ( nbre_y == n_y() ) {
      UTI_3vector<Num_type> const & vecteur = champ.mat(i,0);
      fd << vecteur.x() <<"   "<<  vecteur.y() << "   " <<  vecteur.z() <<  
	std::endl;
    }
  }
  if ( nbre_x == n_x() ) {
    for (int j=0; j < nbre_y; j++) {
      UTI_3vector<Num_type> const & vecteur = champ.mat(0,j);
      fd << vecteur.x() <<"   "<<  vecteur.y() << "   " <<  vecteur.z() 
	 <<  std::endl;
    }
    if ( nbre_y == n_y() ) {
      UTI_3vector<Num_type> const & vecteur = champ.mat(0,0);
      fd << vecteur.x() <<"   "<<  vecteur.y() << "   " <<  vecteur.z() 
	 <<  std::endl;
    }
  }
  fd.close();
}

//-------------------------------------------------------------
//          TOR_embedding<Num_type>::write_Povray()
//
// Writes the embedding in a Povray file
// Also writes a vector field encoded at each vertex
//-------------------------------------------------------------
template < typename Num_type > void TOR_embedding <Num_type>::
write_Povray(string const & nomFichier, 
	     bool reduit, int nbre_x, int nbre_y) const
{
  std::ofstream fd(nomFichier.c_str()); 
  //header
  fd << "// POVRay file" << endl;
  fd << "//" << endl;
  fd << "// +W1348 +H914" << endl;

  fd << "global_settings {" << endl;
  fd << "\t ambient_light color rgb <1.0, 1.0, 1.0>" << endl;
  fd << "\t assumed_gamma 2" << endl;
  fd << "}" << endl;

  fd << "background { color rgb <0.319997, 0.340002, 0.429999>}" << endl;

  fd << "camera {" << endl;
  fd << "\t perspective" << endl;
  fd << "\t location <-0.360250, 0.033578, -0.207749>" << endl;
  fd << "\t sky <-0.145303, -0.984896, 0.094162>" << endl;
  fd << "\t right <-1, 0, 0>" << endl;
  fd << "\t angle 30.000000" << endl;
  fd << "\t look_at <0.204216, -0.018577, 0.117776>" << endl;
  fd << "}" << endl;

  fd << "light_source {" << endl;
  fd << "\t <-0.360250, 0.033578, -0.207749>" << endl;
  fd << "\t color <0.999800, 0.999800, 0.999800>*0.250000" << endl;
  fd << "\t parallel" << endl;
  fd << "\t point_at <0.204216, -0.018577, 0.117776>" << endl;
  fd << "}" << endl;

  fd << "light_source {" << endl;
  fd << "\t <0.111619, 0.766044, 0.633022>" << endl;
  fd << "\t color <1.000000, 0.972320, 0.902220>*0.750000" << endl;
  fd << "\t parallel" << endl;
  fd << "\t point_at <0.000000, 0.000000, 0.000000>" << endl;
  fd << "}" << endl;

  fd << "light_source {" << endl;
  fd << "\t <-0.044943, -0.965926, 0.254887>" << endl;
  fd << "\t color <0.908240, 0.933140, 1.000000>*0.250000" << endl;
  fd << "\t parallel" << endl;
  fd << "\t point_at <0.000000, 0.000000, 0.000000>" << endl;
  fd << "}" << endl;

  fd << "light_source {" << endl;
  fd << "\t <0.939693, 0.000000, -0.342020>" << endl;
  fd << "\t color <0.999800, 0.999800, 0.999800>*0.214286" << endl;
  fd << "\t parallel" << endl;
  fd << "\t point_at <0.000000, 0.000000, 0.000000>" << endl;
  fd << "}" << endl;
  // end header

  // When we glue, we add 1
  int npointsx,npointsy;
  if (! reduit) {
    nbre_x=n_x();
    nbre_y=n_y();
  }
  if ( nbre_x == n_x() ) npointsx = nbre_x +1; else npointsx = nbre_x;
  if ( nbre_y == n_y() ) npointsy = nbre_y +1; else npointsy = nbre_y;

  // Display of points
  fd << "mesh2 {" << endl;
  fd << "\t vertex_vectors {" << endl;
  fd << "\t\t" << npointsx*npointsy << "," << endl;
  for (int i=0; i< nbre_x; i++) {
    for (int j=0; j< nbre_y; j++) {
      UTI_3point<Num_type> const & point = m_(i,j);
      fd << "\t\t<" << point.x() <<", "<<  point.y() << ", " 
	 <<  point.z() << ">," << endl;
    }
    if ( nbre_y == n_y() ) {
      UTI_3point<Num_type> const & point = m_(i,0);
      fd << "\t\t<" << point.x() <<", "<<  point.y() << ", " 
	 <<  point.z() << ">," << endl;
    }
  }
  if ( nbre_x == n_x() ) {
    for (int j=0; j< nbre_y; j++) {
      UTI_3point<Num_type> const & point = m_(0,j);
      fd << "\t\t<" << point.x() <<", "<<  point.y() << ", " 
	 <<  point.z() << ">," << endl;
    }
    if ( nbre_y == n_y() ) {
      UTI_3point<Num_type> const & point = m_(0,0);
      fd << "\t\t<" << point.x() <<", "<<  point.y() << ", " 
	 <<  point.z() << ">," << endl;
    }
  }
  fd << "\t }" << endl;

  
  // The normals
  TOR_3vector_field<Num_type> champ(n_x(), n_y());
  champ = (this->normal_field)(champ);


  fd << "\t normal_vectors {" << endl;
  fd << "\t\t" << npointsx*npointsy << "," << endl;
  for (int i=0; i< nbre_x; i++) {
    for (int j=0; j< nbre_y; j++) {
      UTI_3vector<Num_type> const & vecteur = champ.mat(i,j);
      fd << "\t\t<" << vecteur.x() <<", "<<  vecteur.y() << ", " 
	 <<  vecteur.z() << ">," << endl;
    }
    if ( nbre_y == n_y() ) {
      UTI_3vector<Num_type> const & vecteur = champ.mat(i,0);
      fd << "\t\t<" << vecteur.x() <<", "<<  vecteur.y() << ", " 
	 <<  vecteur.z() << ">," << endl;
    }
  }
  if ( nbre_x == n_x() ) {
    for (int j=0; j < nbre_y; j++) {
      UTI_3vector<Num_type> const & vecteur = champ.mat(0,j);
      fd << "\t\t<" << vecteur.x() <<", "<<  vecteur.y() << ", " 
	 <<  vecteur.z() << ">," << endl;
    }
    if ( nbre_y == n_y() ) {
      UTI_3vector<Num_type> const & vecteur = champ.mat(0,0);
      fd << "\t\t<" << vecteur.x() <<", "<<  vecteur.y() << ", " 
	 <<  vecteur.z() << ">," << endl;
    }
  }
  fd <<  "\t }" << endl;

  // Indices
  fd << "\t face_indices {" << endl;
  fd << "\t\t" << 2*(npointsx-1)*(npointsy-1) << "," << endl;
  for (int i=0; i< npointsx-1; i++) {
    for (int j=0; j< npointsy-1; j++) {
      fd << "\t\t<" <<  j+i*npointsy <<", "<<  j+i*npointsy +1 << ", " 
	 << j+i*npointsy +1 +npointsy << ">," << endl;
      fd << "\t\t<" <<  j+i*npointsy <<", "<<  j+i*npointsy+npointsy << ", " 
	 << j+i*npointsy +1 +npointsy << ">," << endl;
    }
  }
  fd << "\t }" << endl;

  // end of file
  fd << "\t matrix < 1.000000, 0.000000, 0.000000," << endl;
  fd << "\t\t 0.000000, 1.000000, 0.000000," << endl;
  fd << "\t\t 0.000000, 0.000000, 1.000000," << endl;
  fd << "\t\t 0.000000, 0.000000, 0.000000 >" << endl;
  fd << "\t texture {" << endl;
  fd << "\t\t pigment {" << endl;
  fd << "\t\t\t color rgbf <1.000000, 1.000000, 1.000000 0.000000>" << endl;
  fd << "\t\t }" << endl;
  fd << "\t\t finish {" << endl;
  fd << "\t\t\t ambient 0.000000  diffuse 1.000000  phong 0.100000  phong_size 100.000000  " << endl;
  fd << "\t\t }" << endl;
  fd << "\t }" << endl;

  fd << "}" << endl;

  fd.close();
}


//-------------------------------------------------------------
//          TOR_embedding<Num_type>::SortieOFF()
//
// Writes the embedding in an OFF file format
//-------------------------------------------------------------
template < typename Num_type > void TOR_embedding <Num_type>::
write_OFF(string const & nomFichier, 
	  bool reduit, int nbre_x, int nbre_y) const
{
  std::ofstream fd(nomFichier.c_str()); 

  // When we glue, we add 1
  int npointsx,npointsy;
  if (! reduit) {
    nbre_x=n_x();
    nbre_y=n_y();
  }
  int nquadx,nquady;
  if ( nbre_x >= n_x() ) {npointsx = n_x(); nquadx = npointsx;} 
  else {npointsx = nbre_x; nquadx=npointsx-1;}
  if ( nbre_y >= n_y() ) {npointsy = n_y(); nquady = npointsy;} 
  else {npointsy = nbre_y; nquady=npointsy-1;}

  fd << "OFF " <<std::endl;
  fd << npointsx*npointsy << " " << nquadx*nquady << " " << 0 << std::endl;

  // coordinates of points
  for (int i=0; i< nbre_x; i++) {
    for (int j=0; j< nbre_y; j++) {
      UTI_3point<Num_type> const & point = m_(i,j);
      fd << point.x() << "   " << point.y() << "   " << point.z() << std::endl;
    }
  }

  //indices
  for (int i=0; i< nquadx; i++) {
    int iplus1 = (i+1) % npointsx;
    for (int j=0; j< nquady; j++) {
      int jplus1 = (j+1) % npointsy;
      fd << 4 << " " << j+i*npointsy << " " << i*npointsy + jplus1 << " " 
	 << iplus1*npointsy + jplus1 << " " << iplus1*npointsy+j << endl;
    }
  }

  fd.close();
}

//-------------------------------------------------------------
//          TOR_embedding<Num_type>::write_VMRL()
//
// Write a VRML file of the torus with square faces rather
// than triangles. Add an interior standard torus so that 
// 3D printers can spare some consumable material.
//-------------------------------------------------------------
template < typename Num_type > void TOR_embedding <Num_type>::
write_VRML(string const & nomFichier, int n_max) const
{
  int nbre_x = min(n_x(), n_max);
  int nbre_y = min(n_y(), n_max);

  std::ofstream fd(nomFichier.c_str());

  // coef of multiplication -- so that the torus lies in a box of side 250mm -- 
  float coef = 25.4*41.145;

  //header
  fd <<"#VRML V2.0 utf8"<<endl;
fd <<"Shape {"<<endl;
fd <<"\t"<<"appearance Appearance {"<<endl;
fd <<"\t\t"<<"material Material { diffuseColor 0.67 0.77 0.96}"<<endl;
fd <<"\t"<<"}"<<endl;
fd <<"\t"<<"geometry IndexedFaceSet {"<<endl;
fd <<"\t\t"<<"coord Coordinate {"<<endl;
fd <<"\t\t\t"<<"point [ "<<endl;

  // coordinates of points
  for (int i=0; i< nbre_x; i++) {
    for (int j=0; j< nbre_y; j++) {
      UTI_3point<Num_type> const & 
	point = (*this)((Num_type)i/(Num_type)nbre_x, (Num_type)j/(Num_type)nbre_y);
      fd <<"\t\t\t\t"<< fixed << setprecision(2) << coef*point.x() << "  " << coef*point.y() << "  " << coef*point.z() <<"," <<std::endl;
    }
  }
  //coordinates of points of an interior small torus of revolution
  int n_small_torus = 100;
  double r_1=0.12; double r_2=0.5;
  TOR_embedding< Num_type > torus(n_small_torus, n_small_torus, TOR_STANDARD, r_1, r_2);
  for (int i=0; i< n_small_torus; i++) {
    for (int j=0; j< n_small_torus; j++) {
      UTI_3point<Num_type> const & point = torus.mat(i,j);
      fd <<"\t\t\t\t"<< fixed << setprecision(2) << coef*point.x() << "  " << coef*point.y() << "  " << coef*point.z() <<"," <<std::endl;
    }
  }

  //
fd <<"\t\t\t"<<"]"<<endl;
fd <<"\t\t"<<"}"<<endl;
fd <<"\t\t"<<"coordIndex ["<<endl;

  // indices of triangles
  for (int i=0; i< nbre_x; i++) {
    int iplus1 = (i+1) % nbre_x;
    for (int j=0; j< nbre_y; j++) {
      int jplus1 = (j+1) % nbre_y;
      fd << "\t\t" <<  j + i * nbre_y <<", "<<  jplus1 +i * nbre_y << ", " << jplus1 + iplus1 * nbre_y << ", " << j+iplus1*nbre_y  << ", -1," << endl;
    }
  }
  // indices of triangles (small interior torus)
  int N_dec = nbre_x*nbre_y;
  for (int i=0; i< n_small_torus; i++) {
    int iplus1 = (i+1) % n_small_torus;
    for (int j=0; j< n_small_torus; j++) {
      int jplus1 = (j+1) % n_small_torus;
      fd << "\t\t" <<  j + i * n_small_torus + N_dec <<", "<< j + iplus1 * n_small_torus  + N_dec  <<", "<<   jplus1 + iplus1 * n_small_torus  + N_dec  << ", " << jplus1 +i * n_small_torus + N_dec  << ", -1," << endl;
    }
  }

fd <<"\t\t"<<"]"<<endl;
fd <<"\t"<<"}"<<endl;
fd <<"}"<<endl;

  fd.close();
} 

//-------------------------------------------------------------
//          TOR_embedding<Num_type>::write_binary()
//
// Writes the embedding in a binary format
//-------------------------------------------------------------
template < typename Num_type > void TOR_embedding <Num_type>::
write_binary(string const & nomFichier) const
{
  int nx(n_x()), ny(n_y());

  std::ofstream fd(nomFichier.c_str(), ofstream::binary);
  fd.write((char*)&nx, sizeof(int));
  fd.write((char*)&ny, sizeof(int));
  
  Num_type pts[3*ny];
  for (int i=0; i< nx; i++) {
    for (int j=0; j< ny; j++) {
      UTI_3point<Num_type> const & point = m_(i,j);
      pts[3*j] = point.x(); 
      pts[3*j+1] = point.y();
      pts[3*j+2] = point.z();
    }
    fd.write ((char*)pts, 3*ny*sizeof(Num_type));
  }
  fd.close();
}

//-------------------------------------------------------------
//          TOR_embedding<Num_type>::read_binary()
//
// Reads the embedding in a binary format
//-------------------------------------------------------------
template < typename Num_type > TOR_embedding <Num_type>
TOR_embedding <Num_type>::
read_binary(string const & nomFichier)
{
  int nx, ny;

  std::ifstream fd(nomFichier.c_str(), ifstream::binary);
  fd.read((char*)&nx, sizeof(int));
  fd.read((char*)&ny, sizeof(int));

  TOR_embedding<Num_type> torus(nx, ny, TOR_EMPTY);
  
  Num_type pts[3*ny];
  for (int i=0; i< nx; i++) {
    fd.read ((char*)pts, 3*ny*sizeof(Num_type));
    for (int j=0; j< ny; j++) {
      UTI_3point<Num_type> & point = torus.mat(i,j);      
      point.x() = pts[3*j];
      point.y() = pts[3*j+1];
      point.z() = pts[3*j+2];
    }
  }
  if ( fd.good() )
    cout << "File " << nomFichier << " loaded." << endl;
  else
    cout << "Problem while loading the file " << nomFichier << endl;

  fd.close();

  return torus;
}

#endif // _TOR_OUTPUT_H_
