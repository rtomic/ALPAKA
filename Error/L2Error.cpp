/*
 * L2Error.cpp
 *
 *  Created on:   4.09.2019
 *  Last changed: 11.01.2022
 *
 *      Author: kempf
 */


#include "L2Error.h"
#include "CFullGrid.h"
#include "SparseGrid.h"

#include <random>
#include <limits>


L2Error::L2Error()
{
    m_dimension       =  CParam::getFirstValue<std::size_t>("dimension");
	m_threshold       =  CParam::getFirstValue<std::size_t>("errorSparseGridThreshold");

    
    std::vector<HWdouble> low  = CParam::getValue<double>("boundingBox.lo");
    std::vector<HWdouble> high = CParam::getValue<double>("boundingBox.hi");
    
    //Setting Tensor Domain by 1-dim Rectangles

    m_domain.reserve(m_dimension);
    for(std::size_t i = 0; i < m_dimension; i++)
    {
        m_domain.emplace_back(1);
        m_domain[i].lo[0] = low[i];
        m_domain[i].hi[0] = high[i];
    }

    m_choiceGrid = CParam::getFirstValue<std::size_t>("l2ErrorGridChoice");

    switch(m_choiceGrid)
    {
    	case 1:
    	{
    		this->fillUniformGrid();
    		break;
    	}
    	case 2:
    	{
    		this->fillSparseGrid();
    		break;
    	}
    	case 3:
    	{
    		this->fillRandomPoints();
    		break;
    	}
    	default:
    	{
    		std::cout << "Give l2ErrorGridChoice parameter in parameter file" << std::endl;
    		break;
    	}
    }
/*
    std::size_t choiceGrid = CParam::getFirstValue<std::size_t>("l2ErrorGridChoice");
	
	switch(choiceGrid) {
		case 1: this->fillUniformGrid();
			break;
		
		case 2: this->fillSparseGrid();
			break;
		case 3:

		default: std::cout << "Give l2ErrorGridChoice parameter in parameter file" << std::endl;
	}
	*/
}

double L2Error::computeError(std::vector<HWfunction*> const &f, Interpolation &interpolant)
{
    std::size_t choiceGrid = CParam::getFirstValue<std::size_t>("l2ErrorGridChoice");
    double result;

	switch(choiceGrid) {
		case 1: this->fillUniformGrid();
				result = this->uniformGridQuadrature(f, interpolant);
			break;

		case 2: this->fillSparseGrid();
				result = this->sparseGridQuadrature(f, interpolant);
			break;

		case 3: this->fillRandomPoints();
				result = this->monteCarloQuadrature(f, interpolant);
			break;
		
		default: std::cout << "Give l2ErrorGridChoice parameter in parameter file" << std::endl;
	} 
	return result;
}
/*
double L2Error::computeError(std::vector<HWfunction*> const &f, TensorInterpolantVariable &interpolant)
{
    //std::size_t choiceGrid = CParam::getFirstValue<std::size_t>("l2ErrorGridChoice");
    double result;

	switch (m_choiceGrid)
	{
		case 1:
			{
				//this->fillUniformGrid ();
				result = this->uniformGridQuadrature (f, interpolant);
				break;
			}
		case 2:
			{
				//this->fillSparseGrid ();
				result = this->sparseGridQuadrature (f, interpolant);
				break;
			}
		case 3:
			{
				//this->fillRandomPoints ();
				result = this->monteCarloQuadrature (f, interpolant);
				break;
			}
		default:
			{
				std::cout
						<< "Give l2ErrorGridChoice parameter in parameter file"
						<< std::endl;
				break;
			}
	}
	return result;
}
*/
/*
double L2Error::computeError(std::vector<HWfunction*> const &f, TensorInterpolant &interpolant)
{

    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "Evaluating L2-Error on " << m_evalGrid.size() << " points" << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
    
    m_error = 0.0;
    
    std::size_t nThreads = thread::hardware_concurrency();
    std::vector<std::thread> threads;
    threads.resize(nThreads);
    
    std::vector<double> partialResults;
    partialResults.resize(nThreads);
    
    std::size_t nPointEvalPerThread = m_evalGrid.size() / nThreads;
    
    auto func = [&](std::size_t const &startIndex, std::size_t const &endIndex, std::size_t const &threadIndex)
    {
         double tmpProduct;
         double tmpValue;
    
         for(std::size_t i = startIndex; i < endIndex; i++)
         {
            tmpProduct = 1.0;
            for(std::size_t j = 0; j < m_dimension; j++)
            {
                tmpProduct *= f[j]->at(m_evalGrid[i][j]);
            }
            tmpValue = interpolant.at(m_evalGrid[i]);

            partialResults[threadIndex] += (tmpProduct - tmpValue) * (tmpProduct - tmpValue);
        }
    };
    
    
    std::size_t startIdx = 0;
    std::size_t endIdx   = nPointEvalPerThread;
    // Threads work
    for(std::size_t i = 0; i < nThreads - 1; i++)
    {
        threads[i] = std::thread(func, startIdx, endIdx, i);
    
        startIdx += nPointEvalPerThread;
        endIdx   += nPointEvalPerThread;
    }
    
    startIdx = (nThreads - 1) * nPointEvalPerThread;
    endIdx   = m_evalGrid.size();
    threads[nThreads - 1] = std::thread(func, startIdx, endIdx, nThreads - 1);

    for(std::size_t i = 0; i < nThreads; i++)
    {
        threads[i].join();
    }
    
    
    for(std::size_t i = 0; i < nThreads; i++)
    {
        m_error += partialResults[i];
    }
    
    return sqrt(m_error / m_evalGrid.size());
}


double L2Error::computeError(std::vector<HWfunction*> const &f, TensorInterpolantVariable &interpolant)
{

    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "Evaluating L2-Error on " << m_evalGrid.size() << " points" << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;

    m_error = 0.0;

    std::size_t nThreads = thread::hardware_concurrency();
    std::vector<std::thread> threads;
    threads.resize(nThreads);

    std::vector<double> partialResults;
    partialResults.resize(nThreads);

    std::size_t nPointEvalPerThread = m_evalGrid.size() / nThreads;

    auto func = [&](std::size_t const &startIndex, std::size_t const &endIndex, std::size_t const &threadIndex)
    {
         double tmpProduct;
         double tmpValue;

         for(std::size_t i = startIndex; i < endIndex; i++)
         {
            tmpProduct = 1.0;
            for(std::size_t j = 0; j < m_dimension; j++)
            {
                tmpProduct *= f[j]->at(m_evalGrid[i][j]);
            }
            tmpValue = interpolant.at(m_evalGrid[i]);

            partialResults[threadIndex] += (tmpProduct - tmpValue) * (tmpProduct - tmpValue);
        }
    };


    std::size_t startIdx = 0;
    std::size_t endIdx   = nPointEvalPerThread;
    // Threads work
    for(std::size_t i = 0; i < nThreads - 1; i++)
    {
        threads[i] = std::thread(func, startIdx, endIdx, i);

        startIdx += nPointEvalPerThread;
        endIdx   += nPointEvalPerThread;
    }

    startIdx = (nThreads - 1) * nPointEvalPerThread;
    endIdx   = m_evalGrid.size();
    threads[nThreads - 1] = std::thread(func, startIdx, endIdx, nThreads - 1);

    for(std::size_t i = 0; i < nThreads; i++)
    {
        threads[i].join();
    }


    for(std::size_t i = 0; i < nThreads; i++)
    {
        m_error += partialResults[i];
    }

    return sqrt(m_error / m_evalGrid.size());
}
*/
void L2Error::fillSparseGrid()
{
	std::vector<double> weights;
	weights.resize(m_dimension);
	for(std::size_t i = 0; i < m_dimension; i++)
	{
		weights[i] = 1.0;
	}

    std::size_t offset = CParam::getFirstValue<std::size_t>("initialgrid");

    RKIndexSet indexSet(weights, m_dimension, m_threshold);
    std::list<MultiIndex> multiIndices = indexSet.getIndexSet();
    std::vector<double> maxLevelPerDim = indexSet.getMaxLevelPerDim();

    std::vector<RKnestedGrid*> grids;
	grids.resize(m_dimension);
    for(std::size_t i = 0; i < m_dimension; i++)
    {
    	grids[i] = new RKnestedGrid(1,offset,0.5, m_domain[i], maxLevelPerDim[i]);
    }

    SparseGrid sparseGrid(multiIndices, grids);
    m_evalGrid = sparseGrid.getSparseGrid();
}   



void L2Error::fillUniformGrid()
{
	CFullGrid fullGrid(m_domain); // @suppress("Type cannot be resolved")
	
	m_evalGrid = fullGrid.provide_grid_HWpoint(); // @suppress("Method cannot be resolved")
}



void L2Error::fillRandomPoints()
{
	std::size_t nRandomPoints = CParam::getFirstValue<std::size_t>("nRandomPoints");
    m_evalGrid.resize(nRandomPoints);

    float epsilon = std::numeric_limits<float>::denorm_min();

    for(std::size_t i = 0; i < m_evalGrid.size(); i++)
    {
        m_evalGrid[i].resize(m_dimension);
        for(std::size_t j = 0; j < m_dimension; j++)
        {
            std::random_device rd;  //Will be used to obtain a seed for the random number engine
            std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
            std::uniform_real_distribution<> dis(m_domain[j].lo[0], m_domain[j].hi[0] + epsilon);

            m_evalGrid[i][j] = new double [1];
            m_evalGrid[i][j][0] = dis(gen);
        }
    }
}
    

double L2Error::uniformGridQuadrature(std::vector<HWfunction*> const &f, Interpolation &interpolant)
{
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "Evaluating L2-Error on " << m_evalGrid.size() << " points" << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;

    m_error = 0.0;

    std::size_t nThreads = thread::hardware_concurrency();
    std::vector<std::thread> threads;
    threads.resize(nThreads);

    std::vector<double> partialResults;
    partialResults.resize(nThreads);

    std::size_t nPointEvalPerThread = m_evalGrid.size() / nThreads;

    auto func = [&](std::size_t const &startIndex, std::size_t const &endIndex, std::size_t const &threadIndex)
    {
         double tmpProduct;
         double tmpValue;

         for(std::size_t i = startIndex; i < endIndex; i++)
         {
            tmpProduct = 1.0;
            for(std::size_t j = 0; j < m_dimension; j++)
            {
                tmpProduct *= f[j]->at(m_evalGrid[i][j]);
            }
            tmpValue = interpolant.at(m_evalGrid[i]);

            partialResults[threadIndex] += (tmpProduct - tmpValue) * (tmpProduct - tmpValue);
        }
    };


    std::size_t startIdx = 0;
    std::size_t endIdx   = nPointEvalPerThread;
    // Threads work
    for(std::size_t i = 0; i < nThreads - 1; i++)
    {
        threads[i] = std::thread(func, startIdx, endIdx, i);

        startIdx += nPointEvalPerThread;
        endIdx   += nPointEvalPerThread;
    }

    startIdx = (nThreads - 1) * nPointEvalPerThread;
    endIdx   = m_evalGrid.size();
    threads[nThreads - 1] = std::thread(func, startIdx, endIdx, nThreads - 1);

    for(std::size_t i = 0; i < nThreads; i++)
    {
        threads[i].join();
    }


    for(std::size_t i = 0; i < nThreads; i++)
    {
        m_error += partialResults[i];
    }

    return sqrt(m_error / m_evalGrid.size());
}
/*
double L2Error::uniformGridQuadrature(std::vector<HWfunction*> const &f, TensorInterpolantVariable &interpolant)
{

    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "Evaluating L2-Error on " << m_evalGrid.size() << " points" << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;

    m_error = 0.0;

    std::size_t nThreads = thread::hardware_concurrency();
    std::vector<std::thread> threads;
    threads.resize(nThreads);

    std::vector<double> partialResults;
    partialResults.resize(nThreads);

    std::size_t nPointEvalPerThread = m_evalGrid.size() / nThreads;

    auto func = [&](std::size_t const &startIndex, std::size_t const &endIndex, std::size_t const &threadIndex)
    {
         double tmpProduct;
         double tmpValue;

         for(std::size_t i = startIndex; i < endIndex; i++)
         {
            tmpProduct = 1.0;
            for(std::size_t j = 0; j < m_dimension; j++)
            {
                tmpProduct *= f[j]->at(m_evalGrid[i][j]);
            }
            tmpValue = interpolant.at(m_evalGrid[i]);

            partialResults[threadIndex] += (tmpProduct - tmpValue) * (tmpProduct - tmpValue);
        }
    };


    std::size_t startIdx = 0;
    std::size_t endIdx   = nPointEvalPerThread;
    // Threads work
    for(std::size_t i = 0; i < nThreads - 1; i++)
    {
        threads[i] = std::thread(func, startIdx, endIdx, i);

        startIdx += nPointEvalPerThread;
        endIdx   += nPointEvalPerThread;
    }

    startIdx = (nThreads - 1) * nPointEvalPerThread;
    endIdx   = m_evalGrid.size();
    threads[nThreads - 1] = std::thread(func, startIdx, endIdx, nThreads - 1);

    for(std::size_t i = 0; i < nThreads; i++)
    {
        threads[i].join();
    }


    for(std::size_t i = 0; i < nThreads; i++)
    {
        m_error += partialResults[i];
    }

    return sqrt(m_error / m_evalGrid.size());
}
*/
double L2Error::sparseGridQuadrature(std::vector<HWfunction*> const &f, Interpolation &interpolant)
{
	std::cout << "Not yet implemented!" << std::endl;
	return 1;
}
/*
double L2Error::sparseGridQuadrature(std::vector<HWfunction*> const &f, TensorInterpolantVariable &interpolant)
{
	std::cout << "Not yet implemented!" << std::endl;
	return 1;
}
*/

double L2Error::monteCarloQuadrature(std::vector<HWfunction*> const &f, Interpolation &interpolant)
{
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "Evaluating L2-Error on " << m_evalGrid.size() << " points" << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;

    m_error = 0.0;

    std::size_t nThreads = thread::hardware_concurrency();
    std::vector<std::thread> threads;
    threads.resize(nThreads);

    std::vector<double> partialResults;
    partialResults.resize(nThreads);

    std::size_t nPointEvalPerThread = m_evalGrid.size() / nThreads;

    auto func = [&](std::size_t const &startIndex, std::size_t const &endIndex, std::size_t const &threadIndex)
    {
         double tmpProduct;
         double tmpValue;

         for(std::size_t i = startIndex; i < endIndex; i++)
         {
            tmpProduct = 1.0;
            for(std::size_t j = 0; j < m_dimension; j++)
            {
                tmpProduct *= f[j]->at(m_evalGrid[i][j]);
            }
            tmpValue = interpolant.at(m_evalGrid[i]);

            partialResults[threadIndex] += (tmpProduct - tmpValue) * (tmpProduct - tmpValue);
        }
    };


    std::size_t startIdx = 0;
    std::size_t endIdx   = nPointEvalPerThread;
    // Threads work
    for(std::size_t i = 0; i < nThreads - 1; i++)
    {
        threads[i] = std::thread(func, startIdx, endIdx, i);

        startIdx += nPointEvalPerThread;
        endIdx   += nPointEvalPerThread;
    }

    startIdx = (nThreads - 1) * nPointEvalPerThread;
    endIdx   = m_evalGrid.size();
    threads[nThreads - 1] = std::thread(func, startIdx, endIdx, nThreads - 1);

    for(std::size_t i = 0; i < nThreads; i++)
    {
        threads[i].join();
    }


    for(std::size_t i = 0; i < nThreads; i++)
    {
        m_error += partialResults[i];
    }

    return sqrt(m_error / m_evalGrid.size());
}
/*
double L2Error::monteCarloQuadrature(std::vector<HWfunction*> const &f, TensorInterpolantVariable &interpolant)
{
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "Evaluating L2-Error on " << m_evalGrid.size() << " points" << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;

    m_error = 0.0;

    std::size_t nThreads = thread::hardware_concurrency();
    std::vector<std::thread> threads;
    threads.resize(nThreads);

    std::vector<double> partialResults;
    partialResults.resize(nThreads);

    std::size_t nPointEvalPerThread = m_evalGrid.size() / nThreads;

    auto func = [&](std::size_t const &startIndex, std::size_t const &endIndex, std::size_t const &threadIndex)
    {
         double tmpProduct;
         double tmpValue;

         for(std::size_t i = startIndex; i < endIndex; i++)
         {
            tmpProduct = 1.0;
            for(std::size_t j = 0; j < m_dimension; j++)
            {
                tmpProduct *= f[j]->at(m_evalGrid[i][j]);
            }
            tmpValue = interpolant.at(m_evalGrid[i]);

            partialResults[threadIndex] += (tmpProduct - tmpValue) * (tmpProduct - tmpValue);
        }
    };


    std::size_t startIdx = 0;
    std::size_t endIdx   = nPointEvalPerThread;
    // Threads work
    for(std::size_t i = 0; i < nThreads - 1; i++)
    {
        threads[i] = std::thread(func, startIdx, endIdx, i);

        startIdx += nPointEvalPerThread;
        endIdx   += nPointEvalPerThread;
    }

    startIdx = (nThreads - 1) * nPointEvalPerThread;
    endIdx   = m_evalGrid.size();
    threads[nThreads - 1] = std::thread(func, startIdx, endIdx, nThreads - 1);

    for(std::size_t i = 0; i < nThreads; i++)
    {
        threads[i].join();
    }


    for(std::size_t i = 0; i < nThreads; i++)
    {
        m_error += partialResults[i];
    }

    return sqrt(m_error / m_evalGrid.size());
}
*/

