/*
 * RandomPointsError.h
 *
 *  Created on:   4.09.2019
 *  Last changed: 19.08.2025
 *
 *      Author: kempf
 *
 *      Changed m_dimension to vector
 */


#pragma once

#include <mutex>


#include "CError.h"
#include "SitesGenerator.h"



class RandomPointsError: public CError
{
    public:
        RandomPointsError(std::vector<std::vector<double>> errorSites);
        RandomPointsError(std::shared_ptr<std::vector<std::vector<std::vector<double>>>> tensorErrorSites);

        ~RandomPointsError(){};
        void computeError (std::function<double(std::vector<double> const &evalPoint)> const &f,
							 std::unique_ptr<Approximation> const &approximation) override;
        void computeError (std::function<double(std::vector<double> const &evalPoint)> &f,
							 std::unique_ptr<Approximation> const &approximation) override;
        void computeError(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
        		CombinationTechniqueApproximation const &approximation) override;
        double getError() override;
        double getRelativeError() override;

    /*
        double computeError(std::vector<HWfunction*> const &f, TensorInterpolant &interpolant);
        double computeError(std::vector<HWfunction*> const &f, TensorInterpolantVariable &interpolant);
    double computeError(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f, SmolyakMLS &smolyakMLSApprox);

        double computeError_parallel(std::vector<HWfunction*> const &f, TensorInterpolant &interpolant);
        double computeError_parallel(std::vector<HWfunction*> const &f, TensorInterpolantVariable &interpolant);

        double computeError_sequential(std::vector<HWfunction*> const &f, TensorInterpolant &interpolant);
        double computeError_sequential(std::vector<HWfunction*> const &f, TensorInterpolantVariable &interpolant);
   */
    
 //   double computeError(HWfunction* const &f, Interpolation &interpolant)
 //   {
 //       std::cout << "Not yet implemented! " << std::endl;
 //
 //       return -10000000;
 //   };
    
//    double computeError(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f, SmolyakMLS &smolyakMLSApprox);
//    double computeRelativeError(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f, SmolyakMLS &smolyakMLSApprox);
/*
    double computeError(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
						Approximation& approx);
    double computeRelativeError(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
								Approximation& approx);
*/
//    double computeError_sequential(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
//								   SmolyakMLS &smolyakMLSApprox);
//    double computeRelativeError_sequential(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
//										   SmolyakMLS &smolyakMLSApprox);
/*
    double computeError_sequential(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
								   Approximation& approx);
    double computeRelativeError_sequential(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
										   Approximation& approx);
*/
//    double computeError_parallel(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
//								 SmolyakMLS &smolyakMLSApprox);
//    double computeRelativeError_parallel(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
//										 SmolyakMLS &smolyakMLSApprox);
/*
    double computeError_parallel(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
								 Approximation& approx);
    double computeRelativeError_parallel(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
										 Approximation& approx);
*/
        

    private:
        std::vector<std::vector<double>> 	 m_grid;
        double                               m_error       = 0.0;
        double 								 m_relativeError = 0.0;
        
        std::shared_ptr<std::vector<std::vector<std::vector<double>>>> m_tensorGrid_ptr;
        std::vector<std::vector<std::vector<double>>> m_evalPoints;
        std::vector<std::size_t> unravelIndex(std::size_t flatIndex,
                                              std::vector<std::size_t> const &sizes);

        std::vector<std::vector<std::vector<double>>> fillTensorGrid();

        //std::vector<std::size_t>             m_dimension;
        //std::size_t                          m_nDirections;
        //std::size_t                          m_nRandomPoints;
        //std::vector<OrthRect>              	m_domain;
        //std::vector<std::vector<std::vector<double>>>        m_evalPoints;
        //double                               m_error       = 0.0;

        //void fillGrid();
};
