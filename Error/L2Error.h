/*
 * L2Error.h
 *
 *  Created on:   3.09.2019
 *  Last changed: 11.01.2022
 *
 *      Author: kempf
 */



#pragma once

#include "CError.h"
#include "SparseGrid.h"


class L2Error: public CError
{
    public:
        L2Error();
        double computeError(std::vector<HWfunction*> const &f, Interpolation &interpolant);
        double computeError(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f, SmolyakMLS &smolyakMLSApprox)
    {
            std::cout << "Not yet implemented!" << std::endl;
            
            return -100000000;
        };
        //double computeError(std::vector<HWfunction*> const &f, TensorInterpolantVariable &interpolant);

        
    private:
        std::size_t                          m_dimension;
        std::size_t                          m_threshold;
        std::vector<HWorthRect>              m_domain;
        std::vector<std::vector<HWpoint>>    m_evalGrid;
        std::size_t							 m_choiceGrid;
        
        double                               m_error       = 0.0;
        
        void fillSparseGrid();
		void fillUniformGrid();
		void fillRandomPoints();


		double uniformGridQuadrature(std::vector<HWfunction*> const &f, Interpolation &interpolant);
		//double uniformGridQuadrature(std::vector<HWfunction*> const &f, TensorInterpolantVariable &interpolant);
		double sparseGridQuadrature(std::vector<HWfunction*> const &f, Interpolation &interpolant);
		//double sparseGridQuadrature(std::vector<HWfunction*> const &f, TensorInterpolantVariable &interpolant);
		double monteCarloQuadrature(std::vector<HWfunction*> const &f, Interpolation &interpolant);
		//double monteCarloQuadrature(std::vector<HWfunction*> const &f, TensorInterpolantVariable &interpolant);
};



