/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
/**
 * @file implicitToPCL.cpp
 * @ingroup Tools
 * @author Jeremy Levallois (\c jeremy.levallois@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), Universite de Lyon, France
 * LAboratoire de MAthematiques - LAMA (CNRS, UMR 5807), Universite de Savoie, France
 *
 * @date 2012/06/20
 *
 * DGtal 3D curvature shape comparator
 *
 * This file is part of the DGtalTools library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>

#include "DGtal/base/Common.h"
#include "DGtal/base/Clock.h"

//shapes
#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/helpers/Surfaces.h"

#include "DGtal/topology/SetOfSurfels.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/kernel/sets/SetPredicate.h"

//Digitizer
#include "DGtal/shapes/GaussDigitizer.h"


#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/images/imagesSetsUtils/SimpleThresholdForegroundPredicate.h"
#include "DGtal/math/MPolynomial.h"
#include "DGtal/io/readers/MPolynomialReader.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"


//Vol Export
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/colormaps/GrayscaleColorMap.h"
#include "DGtal/io/writers/VolWriter.h"
#include "DGtal/images/ImageHelper.h"

using namespace DGtal;


template <typename Space, typename Shape>
bool
processShape( const std::string & name,
              Shape & aShape,
              double border_min[],
              double border_max[],
              double h,
              const std::string & namePCLFile = "" )
{
  typedef typename Space::RealPoint RealPoint;
  typedef typename Space::Integer Integer;
  typedef GaussDigitizer< Space, Shape > Digitizer;
  typedef KhalimskySpaceND< Space::dimension, Integer > KSpace;
  typedef LightImplicitDigitalSurface< KSpace, Digitizer > LightImplicitDigSurface;
  typedef HyperRectDomain< Space > Domain;
  
  typedef DigitalSurface< LightImplicitDigSurface > MyDigitalSurface;
  typedef typename MyDigitalSurface::ConstIterator ConstIterator;
  typedef typename KSpace::Surfel Surfel;
  
  
  Digitizer dig;
  dig.attach( aShape );
  dig.init( RealPoint( border_min ), RealPoint( border_max ), h );
  Domain domain = dig.getDomain();
  
  typedef typename ImageSelector< Domain, unsigned int >::Type Image;
  Image image( domain );
  DGtal::imageFromRangeAndValue( domain.begin(), domain.end(), image );
  
  KSpace K;
  bool ok = K.init( domain.lowerBound(), domain.upperBound(), true );
  if ( ! ok )
  {
    std::cerr << "[compareShapeEstimators]" << " error in creating KSpace." << std::endl;
    return false;
  }
  
  try
  {
    // Extracts shape boundary
    SurfelAdjacency< KSpace::dimension > SAdj ( true );
    Surfel bel = Surfaces<KSpace>::findABel ( K, dig, 10000 );
    
    LightImplicitDigSurface LightImplDigSurf ( K, dig, SAdj, bel );
    MyDigitalSurface surf ( LightImplDigSurf );
    typedef typename MyDigitalSurface::ConstIterator SurfelConstIterator;

    std::ofstream PCL;
    trace.info() << "Filename = "<<namePCLFile<<std::endl;
    PCL.open( namePCLFile.c_str() );


      
    //Count the number of points
    long int cpt=0;
    for(SurfelConstIterator it = surf.begin(), itend=surf.end(); it != itend; ++it)
      cpt++;
    
    trace.info() << "Surface size = "<<cpt<<std::endl;
    
    PCL << "# range size = " << surf.size() << std::endl;
    PCL << "# h = " << h << std::endl;
    PCL << "# .PCD v.7 - Point Cloud Data file format"<< std::endl;
    PCL <<" VERSION .7"<<std::endl;
    PCL <<" FIELDS x y z"<<std::endl;
    PCL <<" SIZE 4 4 4 4"<<std::endl;
    PCL <<" TYPE I I I I"<<std::endl;
    PCL <<" COUNT 1 1 1 1"<<std::endl;
    PCL <<" WIDTH "<< cpt <<std::endl;
    PCL <<" HEIGHT 1"<<std::endl;
    PCL <<" VIEWPOINT 0 0 0 1 0 0 0"<<std::endl;
    PCL <<" POINTS "<< cpt <<std::endl;
    PCL <<" DATA ascii"<<std::endl;
    PCL <<std::endl;
   
    for(SurfelConstIterator it = surf.begin(), itend=surf.end(); it != itend; ++it)
      PCL << K.sCoord(*it , 0)  << " " << K.sCoord(*it , 1) <<" "<<  K.sCoord(*it , 2) <<std::endl;
    
    PCL.close();
   	
  }
  catch ( InputException e )
  {
    std::cerr << "[estimatorCurvatureComparator3D]"
    << " error."
    << e.what() << std::endl;
    return false;
  }
  return true;
}

void usage( int /*argc*/, char** argv )
{
  std::cerr << "Usage: " << argv[ 0 ] << " <Polynomial> <Px> <Py> <Pz> <Qx> <Qy> <Qz> <step>   <name_of_PCL_file>" << std::endl;
  std::cerr << "\t - displays the boundary of a shape defined implicitly by a 3-polynomial <Polynomial>." << std::endl;
  std::cerr << "\t - P and Q defines the bounding box." << std::endl;
  std::cerr << "\t - step is the grid step." << std::endl;
  std::cerr << "\t - path is optional. It's the path where you want to generate a .vol file of the polynomial shape (for external computation). If no path was set, we don't export as a .vol file." << std::endl;
  std::cerr << "\t - name is optional. It's the name of your .vol file you want to generate (for external computation). If no name was set, we don't export as a .vol file." << std::endl;
}

int main( int argc, char** argv )
{
  if ( argc < 9 )
  {
    usage( argc, argv );
    return 1;
  }
  double border_min[ 3 ];
  double border_max[ 3 ];
  for ( unsigned int i = 0; i < 3; ++i )
  {
    border_min[ i ] = atof( argv[ 2+i ] );
    border_max[ i ] = atof( argv[ 5+i ] );
  }
  double h = atof( argv[ 8 ] );
  
  typedef Z3i::Space::RealPoint RealPoint;
  typedef Z3i::Space::RealPoint::Coordinate Ring;
  typedef MPolynomial< 3, Ring > Polynomial3;
  typedef MPolynomialReader<3, Ring> Polynomial3Reader;
  typedef ImplicitPolynomial3Shape<Z3i::Space> ImplicitShape;
  
  /// Construction of the polynomial shape
  Polynomial3 poly;
  Polynomial3Reader reader;
  std::string poly_str = argv[ 1 ];
  std::string::const_iterator iter = reader.read( poly, poly_str.begin(), poly_str.end() );
  if ( iter != poly_str.end() )
  {
    std::cerr << "ERROR: I read only <"
    << poly_str.substr( 0, iter - poly_str.begin() )
    << ">, and I built P=" << poly << std::endl;
    return 1;
  }
  
  
  bool export_vol = false;
  std::string pathToSavePCLFile = "";
  std::string namePCLFile = "";
  if( argc >= 10 )
  {
    export_vol = true;
    namePCLFile = argv[ 9 ];
  }
  
  ImplicitShape shape( poly );
  
  /// Computation of 3D curvature estimators (mean & Gaussian) on the implicit shape
  processShape< Z3i::Space, ImplicitShape >( poly_str,
     shape,
     border_min, border_max,
     h,
      namePCLFile
     );
}
