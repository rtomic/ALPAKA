/*
 * CError.h
 *
 *  Created on:   3.09.2019
 *  Last changed: 10.11.2025
 *
 *      Author: kempf
 */


#pragma once

//General includes
#include <vector>
#include <list>
#include <iostream>
#include <cmath>
#include <thread>
#include <iterator>
#include <functional>
#include <memory>

#include "Approximation.h"
//Tensor Problem includes
#include "IndexSet.hpp"
//#include "SparseGrid.h"
//#include "smolyakMovingLeastSquares.hpp"
#include "CombinationTechniqueApproximation.h"

/*
//Multilevel includes
#include "multilevel/hwmultilevel.h"
#include "multilevel/hwnestedmultilevel.h"
#include "general/hwfunction.h"
#include "general/hwcenters.h"
#include "general/hwutility.h"
#include "general/hwbasics.h"
#include "general/CParam.h"

*/


class CError
{
    public:

		virtual ~CError(){};

		virtual void computeError(std::function<double(std::vector<double> const &evalPoint)> const &f,
									std::unique_ptr<Approximation> const &approximation)  = 0;
		virtual void computeError(std::function<double(std::vector<double> const &evalPoint)> &f,
									std::unique_ptr<Approximation> const &approximation)  = 0;
		virtual void computeError
			(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
			 CombinationTechniqueApproximation const &approximation) = 0;
        virtual double getError() = 0;
        virtual double getRelativeError() = 0;
        //virtual double computeError(std::vector<HWfunction*> const &f, TensorInterpolantVariable &interpolant) = 0;
        //virtual double computeError(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
		//						SmolyakMLS &smolyakMLSApprox) = 0;
        //virtual double computeRelativeError(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
		//								SmolyakMLS &smolyakMLSApprox) = 0;

    //virtual double computeError(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
    //							Approximation& approx) = 0;
    //virtual double computeRelativeError(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
	//										Approximation& approx) = 0;

};
