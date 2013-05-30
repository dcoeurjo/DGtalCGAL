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

#pragma once

/**
 * @file MongeJetFittingCurvatureEstimator.h
 * @brief Computes the true quantity to each element of a range associated to a parametric shape.
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2011/06/27
 *
 * Header file for module MongeJetFittingCurvatureEstimator.cpp
 *
 * This file is part of the DGtal library.
 *
 * @see testLengthEstimators.cpp, testTrueLocalEstimator.cpp
 */

#if defined(MongeJetFittingCurvatureEstimator_RECURSES)
#error Recursive header files inclusion detected in MongeJetFittingCurvatureEstimator.h
#else // defined(MongeJetFittingCurvatureEstimator_RECURSES)
/** Prevents recursive inclusion of headers. */
#define MongeJetFittingCurvatureEstimator_RECURSES

#if !defined MongeJetFittingCurvatureEstimator_h
/** Prevents repeated inclusion of headers. */
#define MongeJetFittingCurvatureEstimator_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/topology/SCellsFunctors.h>

//CGAL
#include <CGAL/Cartesian.h>
#include <CGAL/Monge_via_jet_fitting.h>
#include <vector>

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{
  /////////////////////////////////////////////////////////////////////////////
  // template class MongeJetFittingCurvatureEstimator
  /**
   * Description of template class 'MongeJetFittingCurvatureEstimator' <p>
   * \brief Aim: Estimates curvature using CGAL Jet Fitting and Monge Form.
   *
   * model of CLocal
   */
  template <typename TSurfel, typename TEmbedder>
  class MongeJetFittingCurvatureEstimator
  {
  public:

    typedef TSurfel Surfel;
    typedef TEmbedder SCellEmbedder;  
    typedef double Quantity;
    typedef typename SCellEmbedder::RealPoint RealPoint;

    typedef CGAL::Cartesian<double> CGALKernel;
    typedef CGALKernel::Point_3  CGALPoint;
    typedef CGAL::Monge_via_jet_fitting<CGALKernel>  CGALMongeViaJet;
    typedef CGALMongeViaJet::Monge_form CGALMongeForm;


    MongeJetFittingCurvatureEstimator(ConstAlias<SCellEmbedder> anEmbedder):
      myEmbedder(anEmbedder) {};


    //We explicitely store the surfel embedding
    void pushSurfel(const Surfel & aSurf)
    {
      RealPoint p = myEmbedder->operator()(aSurf);
      trace.info()<<" Got = "<<p<<std::endl;
      CGALPoint pp(p[0],p[1],p[2]);
      myPoints.push_back(pp);
    }
    
    //Jet fitting and Monge form evaluation
    Quantity eval(const double h)
    {
      CGALMongeForm monge_form;
      CGALMongeViaJet monge_fit;
      
      trace.info() << "Fitting..."<<std::endl;
      
      monge_form = monge_fit(myPoints.begin() , myPoints.end(), 4, 4); 
      
      //OUTPUT on std::cout
      CGAL::set_pretty_mode(std::cout);
      std::cout << "number of points used : " << myPoints.size() << std::endl
                << monge_form;
      std::cout  << "condition_number : " << monge_fit.condition_number() << std::endl
                 << "pca_eigen_vals and associated pca_eigen_vecs :"  << std::endl;
      for (int i=0; i<3; i++)
        std::cout << monge_fit.pca_basis(i).first << std::endl
                  << monge_fit.pca_basis(i).second  << std::endl;
      
      double k1 = monge_form.principal_curvatures ( 0 );
      double k2 = monge_form.principal_curvatures ( 1 );
      //Gaussian curvature 
      return k1*k2;
    }
    
    
    void reset()
    {
      myPoints.clear();
    }
    

  private:
    
    const SCellEmbedder * myEmbedder;

    std::vector<CGALPoint> myPoints;

    
    

  }; // end of class MongeJetFittingCurvatureEstimator

} // namespace DGtal


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined MongeJetFittingCurvatureEstimator_h

#undef MongeJetFittingCurvatureEstimator_RECURSES
#endif // else defined(MongeJetFittingCurvatureEstimator_RECURSES)
