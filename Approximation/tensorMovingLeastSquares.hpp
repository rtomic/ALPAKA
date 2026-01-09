/*
 *  tensorMovingLeastSquares.hpp
 *
 *  Created on: 03.06.2024
 *  Last changed: 25.02.2025
 *      Author: kempf
 *
 *      Class TensorMovingLeastSquares that represents the tensorized moving least squares for given multi-index
 */


#pragma once

#include <vector>
#include <functional>


#include "IndexSet.hpp"
#include "MultiIndex.hpp"
#include "KernelFactory.h"
#include "WendlandKernel.h"
#include "Function.h"
#include "hwtree.h"
#include "hwvector.h"
#include "hwfunction.h"
#include "HWstdInterfaces.hpp"
#include "movingLeastSquares.hpp"

#include <Eigen/Core>
#include <Eigen/QR>

//#include "Approximation.h"

class TensorMovingLeastSquares //: public Approximation
{
public:
    TensorMovingLeastSquares(std::vector<std::vector<std::vector<double>>> const &TPcenters,
							 std::vector<std::size_t> const &kernelTypes,
							 std::vector<double> const &deltas,
							 std::vector<std::size_t> const &polDegrees,
							 std::vector<std::vector<double>> const &shifts,
							 std::vector<double> const &scalings);
    ~TensorMovingLeastSquares();
    
    double at(std::vector<std::vector<double>> const &x,
			  std::function<double(std::vector<std::vector<double>> const &x)> func);
    
private:
    std::size_t m_nDirections;
    std::vector<std::size_t> m_dimensions;
    std::vector<std::vector<std::vector<double>>> m_tpPoints;
    //std::vector<HWpointArray> m_pointsHW;
    std::vector<std::size_t> m_nPts;
    std::vector<double> m_deltas;
    std::vector<std::size_t> m_polDegrees;
    std::vector<std::vector<double>> m_shifts;
    std::vector<double> m_scalings;
    std::vector<std::size_t> m_kernelTypes;
    std::vector<std::size_t> m_dimPolSpaces;
    //std::vector<double> m_functionValues;
    std::vector<MLSApproximation> MLS;
    
    std::vector<HWtree*> m_trees;
    
    std::vector<HWidxArray> m_indices;
    
    
    std::size_t binom(std::size_t const &n, std::size_t const &k);
    double evalPolynomial(std::vector<double> const &evalPoint, std::vector<double> const &shift, double const &scaling, MultiIndex const &degree, std::size_t dimension);
    Eigen::VectorXd forwardSubstitution(Eigen::MatrixXd const &R, Eigen::VectorXd const &RHS);
    Eigen::VectorXd backwardSubstitution(Eigen::MatrixXd const &R, Eigen::VectorXd const &RHS);
};
