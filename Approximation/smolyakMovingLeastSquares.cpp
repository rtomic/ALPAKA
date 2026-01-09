/*
 *  smolyakMovingLeastSquares.cpp
 *
 *  Created on: 20.06.2024
 *  Last changed: 17.02.2025
 *      Author: kempf
 *
 *      Implementation of class SmolyakMLS
 *
 *
 *      centers = Direction LevelInDirection IndexInVector ComponentOfPoint
 */


#include "smolyakMovingLeastSquares.hpp"

SmolyakMLS::SmolyakMLS(std::vector<std::vector<std::vector<std::vector<double>>>> const &centers,
					   std::vector<std::vector<std::size_t>> const &kernelTypes,
					   std::vector<std::vector<double>> const &deltas,
					   std::vector<std::vector<std::size_t>> const &polDegrees,
					   std::vector<std::vector<std::vector<double>>> const &shifts,
					   std::vector<std::vector<double>> const &scalings)
{
    m_weights = CParam::getValue<double>("SparseGridWeights");
    m_threshold = CParam::getFirstValue<std::size_t>("threshold");
    m_nDirections = CParam::getFirstValue<std::size_t>("nDirections");
    
    m_nDimensions.resize(m_nDirections);
    for(std::size_t i = 0; i < m_nDirections; i++)
    {
        m_nDimensions[i] = centers[i][0][0].size();
    }
    
    m_kernelTypes = kernelTypes;
    m_deltas      = deltas;
    m_polDegrees  = polDegrees;
    m_shifts      = shifts;
    m_scalings    = scalings;
    
    SurfaceIndexSet surfaceIndexSet(m_weights,m_nDirections,m_threshold);
    m_surfaceIndexSet = surfaceIndexSet.getSurfaceIndexSet();
    
    m_centers = centers;
}

SmolyakMLS::SmolyakMLS(std::vector<std::vector<std::vector<std::vector<double>>>> const &centers,
					   std::vector<std::vector<std::size_t>> const &kernelTypes,
					   std::vector<std::vector<double>> const &deltas,
					   std::vector<std::vector<std::size_t>> const &polDegrees,
					   std::vector<std::vector<std::vector<double>>> const &shifts,
					   std::vector<std::vector<double>> const &scalings,
					   std::vector<double> const &weights,
					   double const &threshold)
{
    m_weights = weights;
    m_threshold = threshold;
    m_nDirections = CParam::getFirstValue<std::size_t>("nDirections");

    m_nDimensions.resize(m_nDirections);
    for(std::size_t i = 0; i < m_nDirections; i++)
    {
        m_nDimensions[i] = centers[i][0][0].size();
    }

    m_kernelTypes = kernelTypes;
    m_deltas      = deltas;
    m_polDegrees  = polDegrees;
    m_shifts      = shifts;
    m_scalings    = scalings;

    SurfaceIndexSet surfaceIndexSet(m_weights,m_nDirections,m_threshold);
    m_surfaceIndexSet = surfaceIndexSet.getSurfaceIndexSet();

    m_centers = centers;
}

SmolyakMLS::~SmolyakMLS()
{
    
}

double SmolyakMLS::at(std::vector<std::vector<double>> const &x,
					  std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &targetFunc)
{
    
    double result = 0.0;
        
    double weightL1Norm = 0.0;
    for(std::size_t i = 0; i < m_nDirections; i++)
    {
        weightL1Norm += m_weights[i];
    }
    
    std::vector<std::vector<std::vector<double>>> 	productGrid(m_nDirections);
    std::vector<std::size_t> 						kernelTypes (m_nDirections);
    std::vector<double> 							deltas(m_nDirections);
    std::vector<std::size_t> 						polDegrees(m_nDirections);
    std::vector<std::vector<double>> 				shifts(m_nDirections);
    std::vector<double> 							scalings(m_nDirections);
    
    double value;
    
    auto innerSum = [&](std::list<MultiIndex>::iterator const &indexSetIter,
    					TensorMovingLeastSquares  &TMLS)
    {
        double sum = 0.0;
        
        MultiIndex beta(m_nDirections);
        beta.initZeroes();

        while(true)
        {
        	if((*indexSetIter).getWeightedL1Norm(m_weights) + beta.getWeightedL1Norm(m_weights)
            		<= m_threshold * m_weights[0] + weightL1Norm + 1e-10)
            {
                value = TMLS.at(x,targetFunc);
                // odd case
                if(beta.getL1Norm() % 2 == 1)
                {
                    sum -= value;
                }
                // even case
                else
                {
                    sum += value;
                }
            }
            
            for(std::size_t j = 0; j < m_nDirections; j++)
            {
                beta.setValues(j,beta.getMultiIndex()[j]+1);
                if(beta.getMultiIndex()[j] > 1)
                {
                    if(j == m_nDirections-1)
                    {
                        return sum;
                    }
                    beta.setValues(j,0);
                }
                else
                {
                    break;
                }
            }
        }
    };
    
    for(std::list<MultiIndex>::iterator indexSetIter = m_surfaceIndexSet.begin();
        indexSetIter != m_surfaceIndexSet.end(); )
    {
        // Getting parameters for current multiIndex
        for(std::size_t i = 0; i < m_nDirections; i++)
        {
            productGrid[i] = m_centers[i][(*indexSetIter).getMultiIndex()[i]-1];
            kernelTypes[i] = m_kernelTypes[i][(*indexSetIter).getMultiIndex()[i]-1];
            deltas[i] = m_deltas[i][(*indexSetIter).getMultiIndex()[i]-1];
            polDegrees[i] = m_polDegrees[i][(*indexSetIter).getMultiIndex()[i]-1];
            shifts[i] = m_shifts[i][(*indexSetIter).getMultiIndex()[i]-1];
            scalings[i] = m_scalings[i][(*indexSetIter).getMultiIndex()[i]-1];
        }

        TensorMovingLeastSquares TMLS(productGrid,kernelTypes,deltas,polDegrees,shifts,scalings);

        result += innerSum(indexSetIter, TMLS);
        
        /*
        while(true)
        {
            if((*indexSetIter).getWeightedL1Norm(m_weights) + beta.getWeightedL1Norm(m_weights) <= m_threshold * m_weights[0] + weightL1Norm)
            {
                value = TMLS.at(x,targetFunc);
                
                std::cout << "Values: " << value << std::endl;
                // odd case
                if(beta.getL1Norm() % 2 == 1)
                {
                    sum -= value;
                }
                // even case
                else
                {
                    sum += value;
                }
            }
            
            for(std::size_t j = 0; j < m_nDirections; j++)
            {
                beta.setValues(j,beta.getMultiIndex()[j]+1);
                if(beta.getMultiIndex()[j] > 1)
                {
                    if(j == m_nDirections-1)
                    {
                        return sum;
                    }
                    beta.setValues(j,0);
                }
                else
                {
                    break;
                }
            }
        }
         */
        ++indexSetIter;
    }
    
    return result;
}


