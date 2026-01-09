/*
 * LooError.h
 *
 *  Created on:   4.09.2019
 *  Last changed: 19.08.2025
 *
 *      Author: kempf
 */


#pragma once


#include <cmath>
#include <mutex>

#include "CError.h"
#include "SitesGenerator.h"

class LooError: public CError
{    
    public:
        LooError(std::vector<std::vector<double>> errorSites);
        LooError(std::shared_ptr<std::vector<std::vector<std::vector<double>>>> tensorErrorSites);
        ~LooError(){};
        void computeError(std::function<double(std::vector<double> const &evalPoint)> const &f,
        							std::unique_ptr<Approximation> const &approximation) override;
        void computeError(std::function<double(std::vector<double> const &evalPoint)> &f,
        							std::unique_ptr<Approximation> const &approximation) override;
        void computeError(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
        		CombinationTechniqueApproximation const &approximation) override;
        double getError() override;
        double getRelativeError() override;
/*
        double computeRelativeError(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f, SmolyakMLS &smolyakMLSApprox)
        {
        	std::cout << "Not yet implemented!" << std::endl;

        	return -1000000;
		};

  */
    private:
        //std::size_t                          m_dimension;
        //OrthRect			                 m_domain;
        //std::size_t 						 m_n1DPoints;
        //double 								 m_stepSize;
        //std::vector<HWpoint>			     m_evalGridHW;
        //std::vector<std::size_t>			 m_nPoints;
        std::vector<std::vector<double>> 	 m_grid;
        std::shared_ptr<std::vector<std::vector<std::vector<double>>>> m_tensorGrid_ptr;
        double                               m_error       = 0.0;
        double 								 m_relativeError = 0.0;
        
        //void fillGrid();
        //void create();
        //std::vector<std::vector<double>> cartesian_product(const std::vector<std::vector<double>>& vector) const;
};
