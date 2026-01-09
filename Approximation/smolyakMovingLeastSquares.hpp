/*
 * smolyakMovingLeastSquares.hpp
 *
 *  Created on: 20.08.2024
 *  Last edited on: 17.02.2025
 *      Author: kempf
 *
 *      Class that computes the Smolyak moving least squares approximation
 *      to a target function.
 */


#pragma once

#include <iostream>
#include <vector>
#include <functional>
#include <list>

//#include "Approximation.h"

#include "tensorMovingLeastSquares.hpp"
#include "IndexSet.hpp"
#include "MultiIndex.hpp"
#include "SurfaceIndexSet.hpp"

#include "CParam.h"


class SmolyakMLS //: public Approximation
{
public:
    SmolyakMLS(std::vector<std::vector<std::vector<std::vector<double>>>> const &centers,
			   std::vector<std::vector<std::size_t>> const &kernelTypes,
			   std::vector<std::vector<double>> const &deltas,
			   std::vector<std::vector<std::size_t>> const &polDegrees,
			   std::vector<std::vector<std::vector<double>>> const &shifts,
			   std::vector<std::vector<double>> const &scalings);
    SmolyakMLS(std::vector<std::vector<std::vector<std::vector<double>>>> const &centers,
    		   std::vector<std::vector<std::size_t>> const &kernelTypes,
    		   std::vector<std::vector<double>> const &deltas,
    		   std::vector<std::vector<std::size_t>> const &polDegrees,
    		   std::vector<std::vector<std::vector<double>>> const &shifts,
    	       std::vector<std::vector<double>> const &scalings,
    		   std::vector<double> const &weights,
    		   double const &threshold);
    ~SmolyakMLS();
    
    double at(std::vector<std::vector<double>> const &x, std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &targetFunc);
    
private:
    std::vector<double> m_weights;
    double m_threshold;
    std::size_t m_nDirections;
    std::vector<std::size_t> m_nDimensions;
    std::list<MultiIndex> m_surfaceIndexSet;
    
    std::vector<std::vector<std::vector<std::vector<double>>>> m_centers;
    std::vector<std::vector<std::size_t>> m_kernelTypes;
    std::vector<std::vector<double>> m_deltas;
    std::vector<std::vector<std::size_t>> m_polDegrees;
    std::vector<std::vector<std::vector<double>>> m_shifts;
    std::vector<std::vector<double>> m_scalings;
    
    double eval(std::vector<std::vector<double>> const &x, std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &targetFunc,MultiIndex const &lambda);
    
};
