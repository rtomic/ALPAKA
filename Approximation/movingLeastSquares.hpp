/*
 * movingLeastSquares.hpp
 *
 *  Created on: 29.04.2024
 *  Last edited on: 24.02.2025
 *      Author: kempf
 *
 *      Class that computes the moving least squares approximation
 *      to a target function.
 *      For theoretical results, see Holger Wendland - Scattered Data
 *      Approximation
 */


#pragma once

#include "KernelFactory.h"
#include "WendlandKernel.h"
#include "Function.h"
#include "hwtree.h"
#include "hwvector.h"
#include "hwfunction.h"
#include "CParam.h"
#include "IndexSet.hpp"
#include "HWstdInterfaces.hpp"

#include <Eigen/Core>
#include <Eigen/QR>

#include <vector>
#include <iostream>
#include <functional>

#include <stdio.h>

//#include "Approximation.h"



class MLSApproximation //: public Approximation
{
    public:
//    MLSApproximation(std::vector<std::vector<double>> const &centers, std::size_t const &kernelType, double const &delta, HWfunction const *f, std::size_t const &polDegree, std::vector<double> const &shift, double const &scaling);
    MLSApproximation(std::vector<std::vector<double>> const &centers,
					 std::size_t const &kernelType,
                     double const &delta,
					 std::size_t const &polDegree,
					 std::vector<double> const &shift,
					 double const &scaling);

    ~MLSApproximation();
    
    double at(std::vector<double> const &x, HWfunction const *f);
    double at(std::vector<double> const &x, std::function<double(std::vector<double> const &x)> f);

    std::vector<double> shapeFuncValues(std::vector<double> const &x);
    std::vector<double> at(std::vector<std::vector<double>> const &multipleX, HWfunction const *func);
    
    private:
    std::vector<std::vector<double>> m_points;
    HWpointArray m_pointsHW;
    HWtree *m_tree;
    HWidxArray m_index;
    std::size_t m_dimension;
    std::size_t m_nPts;
    std::size_t m_kernelType;
    std::size_t m_polDegree;
    int m_dimPolSpace;
    WendlandKernel* m_kernel;
    HWdouble m_delta;
    
    std::vector<double> m_shift;
    double m_scaling;

    //Eigen::VectorXd m_rightHandSide;
    
    Eigen::VectorXd computeShapeFunctions(std::vector<double> const &x, HWidxArray const &indices,
                                          int const &nLocal);
    Eigen::VectorXd computeShapeFunctionsQR(std::vector<double> const &x, HWidxArray const &indices,
                                            int const &nLocal);
    
    std::size_t binom(std::size_t const &n, std::size_t const &k);
    
    double evalPolynomial(std::vector<double> const &evalPoint, std::vector<double> const &shift, double const &scaling, MultiIndex const &degree);
    
    Eigen::VectorXd forwardSubstitution(Eigen::MatrixXd const &R, Eigen::VectorXd const &RHS);
    Eigen::VectorXd backwardSubstitution(Eigen::MatrixXd const &R, Eigen::VectorXd const &RHS);

};
