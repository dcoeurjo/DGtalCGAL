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
   *
   * @tparam TConstIteratorOnPoints type of iterator on points used as
   * query points.
   * @tparam TParametricShape type of the parametric shape.
   * @tparam TParametricShapeFunctor type of Functor used to evaluate
   * the quantity.
   */
  template <typename TKSpace, typename TShapeFunctor>
  class MongeJetFittingCurvatureEstimator
  {
  public:
    typedef TKSpace KSpace;
    typedef typename Z3i::Domain Domain;
    typedef typename KSpace::Space::RealPoint RealPoint;
    typedef typename KSpace::SCell Cell;

    typedef TShapeFunctor ShapeCellFunctor;

    
    typedef double Quantity;
    
    typedef CGAL::Cartesian<Quantity> CGALKernel;
    typedef CGALKernel::Point_3  CGALPoint;
    typedef CGAL::Monge_via_jet_fitting<CGALKernel>  CGALMongeViaJet;
    typedef CGALMongeViaJet::Monge_form CGALMongeForm;
    
    

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Constructor.
     *
     * @param space space in which the shape is defined.
     * @param f functor on cell of the shape.
     */
    MongeJetFittingCurvatureEstimator (ConstAlias<KSpace>  space,
                                       ConstAlias<ShapeCellFunctor> aShape );
    
    /**
     * Destructor.
     */
    ~MongeJetFittingCurvatureEstimator()
    {}
    
    // ----------------------- Interface --------------------------------------
  public:
    
    /**
     * Initialise the MongeJetFittingCurvatureEstimator with a specific Euclidean kernel radius re, and grid step h.
     *
     * @param _h precision of the grid
     *
     * @bug known bug with radius of kernel. Small hack for the moment.
     */
    template< typename ConstIteratorOnCells >
    void init ( const double _h, ConstIteratorOnCells &itb, ConstIteratorOnCells &ite);
    
    /**
     * Compute the integral invariant mean curvature to cell *it of a shape.
     *
     * @tparam ConstIteratorOnCells iterator on a Cell
     *
     * @param it iterator of a cell (from a shape) we want compute the integral invariant curvature.
     *
     * @return quantity of the result of Integral Invariant estimator at position *it
     */
    template< typename ConstIteratorOnCells >
    Quantity eval ( const ConstIteratorOnCells & it );
    
    /**
     * Compute the integral invariant mean curvature from two cells (from *itb to *ite (exclude) ) of a shape.
     * Return the result on an OutputIterator (param).
     *
     * @tparam ConstIteratorOnCells iterator on a Cell
     * @tparam OutputIterator Output iterator type
     *
     * @param ite iterator of the begin position on the shape where we compute the integral invariant curvature.
     * @param itb iterator of the end position (excluded) on the shape where we compute the integral invariant curvature.
     * @param result iterator of the result of the computation.
     */
    template< typename ConstIteratorOnCells, typename OutputIterator >
    void eval ( const ConstIteratorOnCells & itb,
                const ConstIteratorOnCells & ite,
                OutputIterator & result );
    
    
    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    // ------------------------- Protected Datas ------------------------------
  protected:

    // ------------------------- Private Datas --------------------------------
  private:

   
    //Copy of the KSpace
    KSpace *mySpace;
    
    ///Grid size
    double myH; 
    
    /// origin spel of the kernel support.
    Cell myOrigin;
    

    ///Aliased Shape
    ShapeCellFunctor *myShape;
    
    ///Explicit point set surface
    std::vector<CGALPoint> myPointSet;
    
    
    // ------------------------- Hidden services ------------------------------
  private:
    
    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    MongeJetFittingCurvatureEstimator<KSpace,ShapeCellFunctor> ( const MongeJetFittingCurvatureEstimator<KSpace,ShapeCellFunctor> & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    MongeJetFittingCurvatureEstimator<KSpace,ShapeCellFunctor> & operator= ( const MongeJetFittingCurvatureEstimator<KSpace,ShapeCellFunctor> & other );


  }; // end of class MongeJetFittingCurvatureEstimator

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "MongeJetFittingCurvatureEstimator.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined MongeJetFittingCurvatureEstimator_h

#undef MongeJetFittingCurvatureEstimator_RECURSES
#endif // else defined(MongeJetFittingCurvatureEstimator_RECURSES)
