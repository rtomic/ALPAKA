/*
 * RandomPointsError.cpp
 *
 *  Created on:  28.08.2019
 *  Last changed:04.09.2025
 *
 *      Author: kempf
 */

#include "RandomPointsError.h"

#include <random>
#include <limits>


RandomPointsError::RandomPointsError(std::vector<std::vector<double>> errorSites)
{
	m_grid = errorSites;
	/*
    m_dimension       =  CParam::getValue<std::size_t>("dimension");
    m_nRandomPoints   =  CParam::getFirstValue<std::size_t>("nRandomPoints");
    m_nDirections = CParam::getFirstValue<std::size_t>("nDirections");
    
    std::vector<double> low  = CParam::getValue<double>("boundingBox.lo");
    std::vector<double> high = CParam::getValue<double>("boundingBox.hi");
    
    //Setting Tensor Domain by 1-dim Rectangles

    m_domain.reserve(m_nDirections);
    for(std::size_t i = 0; i < m_nDirections; i++)
    {
        m_domain.emplace_back(1);
        m_domain[i].lo[0] = low[i];
        m_domain[i].hi[0] = high[i];
    }
    this->fillGrid();
	*/
}


RandomPointsError::RandomPointsError
		(std::shared_ptr<std::vector<std::vector<std::vector<double>>>> tensorErrorSites)
{
	m_tensorGrid_ptr = tensorErrorSites;
	m_evalPoints = fillTensorGrid();
}

double RandomPointsError::getError()
{
	return m_error;
}

double RandomPointsError::getRelativeError()
{
	return m_relativeError;
}

void RandomPointsError::computeError (std::function<double(std::vector<double> const &evalPoint)> const &f,
										 std::unique_ptr<Approximation> const &approximation)
{

		std::size_t nPoints = m_grid.size();
		std::size_t nThreads = std::thread::hardware_concurrency();
        std::mutex mutex;
		std::vector<std::thread> threads;
		std::size_t chunkSize = (nPoints + nThreads - 1) / nThreads;

		double globalMaxTargetValue = 0.0;
		m_error = 0.0;
		m_relativeError = 0.0;

		for (std::size_t t = 0; t < nThreads; ++t)
		{
			std::size_t begin = t * chunkSize;
			std::size_t end = std::min (begin + chunkSize, nPoints);

			threads.emplace_back ( [&, begin, end] ()
			{
				double localMaxAbsError = 0.0;
				double localMaxTargetValue = 0.0;

				for (std::size_t i = begin; i < end; ++i)
				{
					const auto& pt = m_grid[i];
					double approxVal = approximation->at(pt);
					double targetVal = f(pt);
					double absError = std::abs(approxVal - targetVal);
					localMaxAbsError = std::max(localMaxAbsError, absError);
					localMaxTargetValue = std::max(localMaxTargetValue, std::abs(targetVal));
				}

				std::lock_guard<std::mutex> lock(mutex);
				m_error = std::max(m_error, localMaxAbsError);
				globalMaxTargetValue = std::max(globalMaxTargetValue, localMaxTargetValue);
			});
		}

		for (auto &th : threads)
		{
			if (th.joinable ())
				th.join ();
		}

		if (globalMaxTargetValue == 0.0)
		{
			std::cout << "Approximating zero Function, no relative Error" << std::endl;
		}

		m_relativeError = m_error / globalMaxTargetValue;
}

void RandomPointsError::computeError (std::function<double(std::vector<double> const &evalPoint)> &f,
										 std::unique_ptr<Approximation> const &approximation)
{

		std::size_t nPoints = m_grid.size();
		std::size_t nThreads = std::thread::hardware_concurrency();
        std::mutex mutex;
		std::vector<std::thread> threads;
		std::size_t chunkSize = (nPoints + nThreads - 1) / nThreads;

		double globalMaxTargetValue = 0.0;
		m_error = 0.0;
		m_relativeError = 0.0;

		for (std::size_t t = 0; t < nThreads; ++t)
		{
			std::size_t begin = t * chunkSize;
			std::size_t end = std::min (begin + chunkSize, nPoints);

			threads.emplace_back ( [&, begin, end] ()
			{
				double localMaxAbsError = 0.0;
				double localMaxTargetValue = 0.0;

				for (std::size_t i = begin; i < end; ++i)
				{
					const auto& pt = m_grid[i];
					double approxVal = approximation->at(pt);
					double targetVal = f(pt);
					double absError = std::abs(approxVal - targetVal);
					localMaxAbsError = std::max(localMaxAbsError, absError);
					localMaxTargetValue = std::max(localMaxTargetValue, std::abs(targetVal));
				}

				std::lock_guard<std::mutex> lock(mutex);
				m_error = std::max(m_error, localMaxAbsError);
				globalMaxTargetValue = std::max(globalMaxTargetValue, localMaxTargetValue);
			});
		}

		for (auto &th : threads)
		{
			if (th.joinable ())
				th.join ();
		}

		if (globalMaxTargetValue == 0.0)
		{
			std::cout << "Approximating zero Function, no relative Error" << std::endl;
		}

		m_relativeError = m_error / globalMaxTargetValue;
}



void RandomPointsError::computeError
(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
	CombinationTechniqueApproximation const &approximation)
{
	//std::size_t nDirections = m_tensorGrid_ptr->size();

	//std::vector<std::vector<std::vector<double>>> evalPoints = fillTensorGrid();


	std::size_t nPoints = m_evalPoints.size();

	std::cout << "COMPUTING RANDOM POINT ERROR ON " << nPoints << " POINTS" << std::endl;

		std::size_t nThreads;
		if(CParam::getFirstValue<std::size_t>("ErrorParallel") == 0)
		{
			nThreads = std::thread::hardware_concurrency();
		}
		else
		{
			nThreads = 1;
		}

		std::vector<std::thread> threads;
		std::size_t chunkSize = (nPoints + nThreads - 1) / nThreads;

		m_error = 0.0;
		m_relativeError = 0.0;

		double globalMaxTargetValue = 0.0;

	    std::mutex mutex;

		for (std::size_t t = 0; t < nThreads; ++t)
		{
			std::size_t begin = t * chunkSize;
			std::size_t end = std::min (begin + chunkSize, nPoints);

			threads.emplace_back ( [&, begin, end] ()
			{
				double localMaxAbsError = 0.0;
				double localMaxTargetValue = 0.0;
/*
				std::vector<std::size_t> sizes(nDirections);
				for (std::size_t d = 0; d < nDirections; d++)
				{
					sizes[d] = (*m_tensorGrid_ptr)[d].size();
				}
*/
				for (std::size_t i = begin; i < end; ++i)
				{
/*
				    auto multiIndex = unravelIndex(i, sizes);

				    std::vector<std::vector<double>> pt(nDirections);
				    for (std::size_t j = 0; j < nDirections; j++)
				    {
				        pt[j] = (*m_tensorGrid_ptr)[j][multiIndex[j]];
				    }
*/
				    double approxVal = approximation.at(m_evalPoints[i], f);

				    //std::cout << "approxVal: " << approxVal << std::endl;
				    double targetVal = f(m_evalPoints[i]);

				    double absError  = std::abs(approxVal - targetVal);
				    localMaxAbsError = std::max(localMaxAbsError, absError);
				    localMaxTargetValue = std::max(localMaxTargetValue, std::abs(targetVal));
				}

				{
					std::lock_guard<std::mutex> lock(mutex);
					m_error = std::max(m_error, localMaxAbsError);
					globalMaxTargetValue = std::max(globalMaxTargetValue, localMaxTargetValue);
				}
				});
		}

		for (auto &th : threads)
		{
			if (th.joinable ())
				th.join ();
		}

		if (globalMaxTargetValue == 0.0)
		{
			std::cout << "Approximating zero Function, no relative Error" << std::endl;
			m_relativeError = -10000;
		}
		else
		{
			m_relativeError = m_error / globalMaxTargetValue;
		}

}


std::vector<std::vector<std::vector<double>>> RandomPointsError::fillTensorGrid()
{
	std::size_t nDirections = m_tensorGrid_ptr->size();
	std::size_t nPoints = 1;
	for(std::size_t i = 0; i < nDirections; i++)
	{
		nPoints *= (*m_tensorGrid_ptr)[i].size();
	}
	std::vector<std::vector<std::vector<double>>> evalPoints(nPoints);

	std::vector<std::size_t> pointCounter(nDirections,0);
	std::size_t counter = 0;

	while(true)
	{
		evalPoints[counter].resize(nDirections);
		for(std::size_t dir = 0; dir < nDirections; dir++)
		{
			evalPoints[counter][dir] = (*m_tensorGrid_ptr)[dir][pointCounter[dir]];
		}
		counter++;

		for(std::size_t j = 0; j < nDirections; j++)
		{
			pointCounter[j] += 1;
			if(pointCounter[j] > (*m_tensorGrid_ptr)[j].size()-1)
			{
				if(j == nDirections-1)
				{
					return evalPoints;
				}
				pointCounter[j] = 0;
			}
			else
			{
				break;
			}
		}
	}
}

/*
void RandomPointsError::computeError
	(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
		CombinationTechniqueApproximation const &approximation)
{
	std::size_t nDirections = m_tensorGrid_ptr->size();
	std::size_t nPoints = 1;
	for(std::size_t i = 0; i < nDirections; i++)
	{
		nPoints *= (*m_tensorGrid_ptr)[i].size();
	}

	std::cout << "COMPUTING RANDOM POINT ERROR ON " << nPoints << " POINTS" << std::endl;

	std::size_t nThreads;
	if(CParam::getFirstValue<std::size_t>("ErrorParallel") == 0)
	{
		nThreads = std::thread::hardware_concurrency();
	}
	else
	{
		nThreads = 1;
	}

	std::vector<std::thread> threads;
	std::size_t chunkSize = (nPoints + nThreads - 1) / nThreads;

	m_error = 0.0;
	m_relativeError = 0.0;

	double globalMaxTargetValue = 0.0;

    std::mutex mutex;

	for (std::size_t t = 0; t < nThreads; ++t)
	{
		std::size_t begin = t * chunkSize;
		std::size_t end = std::min (begin + chunkSize, nPoints);

		threads.emplace_back ( [&, begin, end] ()
		{
			double localMaxAbsError = 0.0;
			double localMaxTargetValue = 0.0;

			std::vector<std::size_t> sizes(nDirections);
			for (std::size_t d = 0; d < nDirections; d++)
			{
				sizes[d] = (*m_tensorGrid_ptr)[d].size();
			}

			for (std::size_t i = begin; i < end; ++i)
			{
			    auto multiIndex = unravelIndex(i, sizes);

			    std::vector<std::vector<double>> pt(nDirections);
			    for (std::size_t j = 0; j < nDirections; j++)
			    {
			        pt[j] = (*m_tensorGrid_ptr)[j][multiIndex[j]];
			    }

			    double approxVal = approximation.at(pt, f);
			    double targetVal = f(pt);
			    double absError  = std::abs(approxVal - targetVal);
			    localMaxAbsError = std::max(localMaxAbsError, absError);
			    localMaxTargetValue = std::max(localMaxTargetValue, std::abs(targetVal));
			}

			{
				std::lock_guard<std::mutex> lock(mutex);
				m_error = std::max(m_error, localMaxAbsError);
				globalMaxTargetValue = std::max(globalMaxTargetValue, localMaxTargetValue);
			}
			});
	}

	for (auto &th : threads)
	{
		if (th.joinable ())
			th.join ();
	}

	if (globalMaxTargetValue == 0.0)
	{
		std::cout << "Approximating zero Function, no relative Error" << std::endl;
		m_relativeError = -10000;
	}
	else
	{
		m_relativeError = m_error / globalMaxTargetValue;
	}
}

std::vector<std::size_t> RandomPointsError::unravelIndex(std::size_t flatIndex,
                                      std::vector<std::size_t> const &sizes)
{
    std::vector<std::size_t> multiIndex(sizes.size());
    for (std::size_t dim = sizes.size(); dim-- > 0; )
    {
        multiIndex[dim] = flatIndex % sizes[dim];
        flatIndex /= sizes[dim];
    }
    return multiIndex;
}
*/





/*
void RandomPointsError::fillGrid()
{
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


    m_evalPoints.resize(m_nDirections);
    for(std::size_t i = 0; i < m_nDirections; i++)
    {
        m_evalPoints[i].resize(m_nRandomPoints);
        for(std::size_t j = 0; j < m_nRandomPoints; j++)
        {
            m_evalPoints[i][j].resize(m_dimension[i]);
            for(std::size_t k = 0; k < m_dimension[i]; k++)
            {
                std::random_device rd;  //Will be used to obtain a seed for the random number engine
                std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
                std::uniform_real_distribution<> dis(m_domain[i].lo[0], m_domain[i].hi[0] + epsilon);
                
                m_evalPoints[i][j][k] = dis(gen);
            }
        }
    }
    
}
 */
/*

double RandomPointsError::computeError(std::vector<HWfunction*> const &f, TensorInterpolant &interpolant)
{
	double result;

	if(m_nRandomPoints > 1000)
	{
		result = computeError_parallel(f, interpolant);
	}
	else
	{
		result = computeError_sequential(f, interpolant);
	}

	return result;
}

double RandomPointsError::computeError(std::vector<HWfunction*> const &f, TensorInterpolantVariable &interpolant)
{
	double result;

	if(m_nRandomPoints > 1000)
	{
		result = computeError_parallel(f, interpolant);
	}
	else
	{
		result = computeError_sequential(f, interpolant);
	}

	return result;
}


double RandomPointsError::computeError_parallel(std::vector<HWfunction*> const &f, TensorInterpolant &interpolant)
{
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
            tmpValue = std::abs(tmpProduct - interpolant.at(m_evalGrid[i]));

            if(partialResults[threadIndex] < tmpValue)
            {
                partialResults[threadIndex] = tmpValue;
            }
        }
    };


    std::size_t startIdx = 0;
    std::size_t endIdx   = nPointEvalPerThread;
    // Threads work
    for(std::size_t i = 0; i < nThreads - 1; i++)
    {
        partialResults[i] = 0.0;

        threads[i] = std::thread(func, startIdx, endIdx, i);

        startIdx += nPointEvalPerThread;
        endIdx   += nPointEvalPerThread;
    }

    startIdx = (nThreads - 1) * nPointEvalPerThread;
    endIdx   = m_evalGrid.size();

    partialResults[nThreads - 1] = 0.0;

    threads[nThreads - 1] = std::thread(func, startIdx, endIdx, nThreads - 1);



    for(std::size_t i = 0; i < nThreads; i++)
    {
        threads[i].join();
    }


    for(std::size_t i = 0; i < nThreads; i++)
    {
        if(m_error < partialResults[i])
        {
            m_error = partialResults[i];
        }
    }

//     double tmp = 0.0;
//     double tmpProduct;
//     for(std::size_t i = 0; i < m_nRandomPoints; i++)
//     {
//         tmpProduct = 1.0;
//         for(std::size_t j = 0; j < m_dimension; j++)
//         {
//             tmpProduct *= f[j]->at(m_evalGrid[i][j]);
//         }
//         tmp = std::abs(tmpProduct - interpolant.at(m_evalGrid[i]));
//
//         if(m_error < tmp)
//         {
//             m_error = tmp;
//         }
//     }

    return m_error;
}


double RandomPointsError::computeError_parallel(std::vector<HWfunction*> const &f, TensorInterpolantVariable &interpolant)
{
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
            tmpValue = std::abs(tmpProduct - interpolant.at(m_evalGrid[i]));

            if(partialResults[threadIndex] < tmpValue)
            {
                partialResults[threadIndex] = tmpValue;
            }
        }
    };
    
    
    std::size_t startIdx = 0;
    std::size_t endIdx   = nPointEvalPerThread;
    // Threads work
    for(std::size_t i = 0; i < nThreads - 1; i++)
    {
        partialResults[i] = 0.0;
        
        threads[i] = std::thread(func, startIdx, endIdx, i);
    
        startIdx += nPointEvalPerThread;
        endIdx   += nPointEvalPerThread;
    }
    
    startIdx = (nThreads - 1) * nPointEvalPerThread;
    endIdx   = m_evalGrid.size();
    
    partialResults[nThreads - 1] = 0.0;

    threads[nThreads - 1] = std::thread(func, startIdx, endIdx, nThreads - 1);
    
    

    for(std::size_t i = 0; i < nThreads; i++)
    {
        threads[i].join();
    }
    
    
    for(std::size_t i = 0; i < nThreads; i++)
    {
        if(m_error < partialResults[i])
        {
            m_error = partialResults[i];
        }
    }
     
//     double tmp = 0.0;
//     double tmpProduct;   
//     for(std::size_t i = 0; i < m_nRandomPoints; i++)
//     {        
//         tmpProduct = 1.0;
//         for(std::size_t j = 0; j < m_dimension; j++)
//         {
//             tmpProduct *= f[j]->at(m_evalGrid[i][j]);
//         }
//         tmp = std::abs(tmpProduct - interpolant.at(m_evalGrid[i]));
//             
//         if(m_error < tmp)
//         {
//             m_error = tmp;
//         }
//     }
    
    return m_error;
}


double RandomPointsError::computeError_sequential(std::vector<HWfunction*> const &f, TensorInterpolant &interpolant)
{
    m_error = 0.0;

    double tmp = 0.0;
    double tmpProduct;   
    
    for(std::size_t i = 0; i < m_nRandomPoints; i++)
    {        
        tmpProduct = 1.0;
        for(std::size_t j = 0; j < m_dimension; j++)
        {
            tmpProduct *= f[j]->at(m_evalGrid[i][j]);
        }
        tmp = std::abs(tmpProduct - interpolant.at(m_evalGrid[i]));
            
        if(m_error < tmp)
        {
            m_error = tmp;
        }
    }
    
    return m_error;
}

double RandomPointsError::computeError_sequential(std::vector<HWfunction*> const &f, TensorInterpolantVariable &interpolant)
{
    m_error = 0.0;

    double tmp = 0.0;
    double tmpProduct;

    for(std::size_t i = 0; i < m_nRandomPoints; i++)
    {
        tmpProduct = 1.0;
        for(std::size_t j = 0; j < m_dimension; j++)
        {
            tmpProduct *= f[j]->at(m_evalGrid[i][j]);
        }
        tmp = std::abs(tmpProduct - interpolant.at(m_evalGrid[i]));

        if(m_error < tmp)
        {
            m_error = tmp;
        }
    }

    return m_error;
}

*/
/*
double RandomPointsError::computeError(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
									   SmolyakMLS &smolyakMLSApprox)
{
    double result;

    if(m_nRandomPoints > 100)
    {
        result = computeError_parallel(f, smolyakMLSApprox);
    }
    else
    {
        result = computeError_sequential(f, smolyakMLSApprox);
    }

    return result;
}


double RandomPointsError::computeRelativeError(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
											   SmolyakMLS &smolyakMLSApprox)
{
    double result;

    if(m_nRandomPoints > 100)
    {
        result = computeRelativeError_parallel(f, smolyakMLSApprox);
    }
    else
    {
        result = computeRelativeError_sequential(f, smolyakMLSApprox);
    }

    return result;
}

/*

double RandomPointsError::computeError(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
									   Approximation &approx)
{
    double result;

    if(m_nRandomPoints > 100)
    {
        result = computeError_parallel(f, approx);
    }
    else
    {
        result = computeError_sequential(f, approx);
    }

    return result;
}


double RandomPointsError::computeRelativeError(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
											   Approximation &approx)
{
    double result;

    if(m_nRandomPoints > 100)
    {
        result = computeRelativeError_parallel(f, approx);
    }
    else
    {
        result = computeRelativeError_sequential(f, approx);
    }

    return result;
}

*/
/*
double RandomPointsError::computeError_parallel(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
												SmolyakMLS &smolyakMLSApprox)
{
    m_error = 0.0;

    std::size_t nThreads = thread::hardware_concurrency();
    std::vector<std::thread> threads;
    threads.resize(nThreads);

    std::vector<double> partialResults;
    partialResults.resize(nThreads);
    
    std::size_t nPointEvalPerThread = m_nRandomPoints / nThreads;

    auto func = [&](std::size_t const &startIndex, std::size_t const &endIndex, std::size_t const &threadIndex)
    {
         double tmpValue;
         std::vector<std::vector<double>> evalPoint(m_nDirections);

         for(std::size_t i = startIndex; i < endIndex; i++)
         {
             for(std::size_t direction = 0; direction < m_nDirections; direction++)
             {
            	 evalPoint[direction].resize(m_dimension[direction]);
            	 for(std::size_t dimension = 0; dimension < m_dimension[direction]; dimension++)
            	 {
            		 evalPoint[direction][dimension] = m_evalPoints[direction][i][dimension];
            	 }
             }
            //tmpValue = std::abs(f(m_evalPoints[i]) - smolyakMLSApprox.at(m_evalPoints[i],f));
             tmpValue = std::abs(f(evalPoint) - smolyakMLSApprox.at(evalPoint,f));


            if(partialResults[threadIndex] < tmpValue)
            {
                partialResults[threadIndex] = tmpValue;
            }
        }
    };


    std::size_t startIdx = 0;
    std::size_t endIdx   = nPointEvalPerThread;
    // Threads work
    for(std::size_t i = 0; i < nThreads - 1; i++)
    {
        partialResults[i] = 0.0;

        threads[i] = std::thread(func, startIdx, endIdx, i);

        startIdx += nPointEvalPerThread;
        endIdx   += nPointEvalPerThread;
    }

    startIdx = (nThreads - 1) * nPointEvalPerThread;
    endIdx   = m_evalGrid.size();

    partialResults[nThreads - 1] = 0.0;

    threads[nThreads - 1] = std::thread(func, startIdx, endIdx, nThreads - 1);



    for(std::size_t i = 0; i < nThreads; i++)
    {
        threads[i].join();
    }


    for(std::size_t i = 0; i < nThreads; i++)
    {
        if(m_error < partialResults[i])
        {
            m_error = partialResults[i];
        }
    }

    return m_error;
}



double RandomPointsError::computeRelativeError_parallel(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
														SmolyakMLS &smolyakMLSApprox)
{
    m_error = 0.0;

    std::size_t nThreads = thread::hardware_concurrency();
    std::vector<std::thread> threads;
    threads.resize(nThreads);

    std::vector<double> partialResults;
    partialResults.resize(nThreads);

    std::vector<double> maxValueF;
    maxValueF.resize(nThreads);

    std::size_t nPointEvalPerThread = m_nRandomPoints / nThreads;

    auto func = [&](std::size_t const &startIndex, std::size_t const &endIndex,
    		std::size_t const &threadIndex)
    {
         double tmpValue;
         double tmpValueF;
         std::vector<std::vector<double>> evalPoint(m_nDirections);

         for(std::size_t i = startIndex; i < endIndex; i++)
         {
             for(std::size_t direction = 0; direction < m_nDirections; direction++)
             {
            	 evalPoint[direction].resize(m_dimension[direction]);
            	 for(std::size_t dimension = 0; dimension < m_dimension[direction]; dimension++)
            	 {
            		 evalPoint[direction][dimension] = m_evalPoints[direction][i][dimension];
            	 }
             }
            //tmpValue = std::abs(f(m_evalPoints[i]) - smolyakMLSApprox.at(m_evalPoints[i],f));
             tmpValue = std::abs(f(evalPoint) - smolyakMLSApprox.at(evalPoint,f));
             tmpValueF = std::abs(f(evalPoint));

            if(partialResults[threadIndex] < tmpValue)
            {
                partialResults[threadIndex] = tmpValue;
            }
            if(maxValueF[threadIndex] < tmpValueF)
            {
            	maxValueF[threadIndex] = tmpValueF;
            }
        }
    };


    std::size_t startIdx = 0;
    std::size_t endIdx   = nPointEvalPerThread;
    // Threads work
    for(std::size_t i = 0; i < nThreads - 1; i++)
    {
        partialResults[i] = 0.0;
        maxValueF[i] = 0.0;

        threads[i] = std::thread(func, startIdx, endIdx, i);

        startIdx += nPointEvalPerThread;
        endIdx   += nPointEvalPerThread;
    }

    startIdx = (nThreads - 1) * nPointEvalPerThread;
    endIdx   = m_evalGrid.size();

    partialResults[nThreads - 1] = 0.0;
    maxValueF[nThreads - 1] = 0.0;

    threads[nThreads - 1] = std::thread(func, startIdx, endIdx, nThreads - 1);



    for(std::size_t i = 0; i < nThreads; i++)
    {
        threads[i].join();
    }

    double maxValueFTogether = 0.0;
    for(std::size_t i = 0; i < nThreads; i++)
    {
        if(m_error < partialResults[i])
        {
            m_error = partialResults[i];
        }
        if(maxValueFTogether < maxValueF[i])
        {
        	maxValueFTogether = maxValueF[i];
        }
    }

    return m_error / maxValueFTogether;
}



double RandomPointsError::computeError_sequential(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
												  SmolyakMLS &smolyakMLSApprox)
{
    m_error = 0.0;

    double tmp = 0.0;

    std::vector<std::vector<double>> evalPoint (m_nDirections);


    for(std::size_t i = 0; i < m_nRandomPoints; i++)
    {
    	for (std::size_t direction = 0; direction < m_nDirections; direction++)
    	{
    		evalPoint[direction].resize (m_dimension[direction]);
    		for (std::size_t dimension = 0; dimension < m_dimension[direction]; dimension++)
    		{
    			evalPoint[direction][dimension] =
    					m_evalPoints[direction][i][dimension];
    		}
    	}
        //tmp = std::abs(f(m_evalPoints[i]) - smolyakMLSApprox.at(m_evalPoints[i],f));
        tmp = std::abs(f(evalPoint) - smolyakMLSApprox.at(evalPoint,f));

        if(m_error < tmp)
        {
            m_error = tmp;
        }
    }

    return m_error;
}


double RandomPointsError::computeRelativeError_sequential(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
														  SmolyakMLS &smolyakMLSApprox)
{
    m_error = 0.0;

    double tmp = 0.0;
    double maxValueF;
    double tmpFValue = 0.0;
    std::vector<std::vector<double>> evalPoint (m_nDirections);


    for(std::size_t i = 0; i < m_nRandomPoints; i++)
    {
    	for (std::size_t direction = 0; direction < m_nDirections; direction++)
    	{
    		evalPoint[direction].resize (m_dimension[direction]);
    		for (std::size_t dimension = 0; dimension < m_dimension[direction]; dimension++)
    		{
    			evalPoint[direction][dimension] =
    					m_evalPoints[direction][i][dimension];
    		}
    	}
        //tmp = std::abs(f(m_evalPoints[i]) - smolyakMLSApprox.at(m_evalPoints[i],f));
        tmp = std::abs(f(evalPoint) - smolyakMLSApprox.at(evalPoint,f));
        tmpFValue = std::abs(f(evalPoint));

        if(m_error < tmp)
        {
            m_error = tmp;
        }
        if(maxValueF < tmpFValue)
        {
        	maxValueF = tmpFValue;
        }
    }

    return m_error / maxValueF;
}


//----------------------------------------------------------------------------------
/*
double RandomPointsError::computeRelativeError_parallel(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
														Approximation& approx)
{
    m_error = 0.0;

    std::size_t nThreads = thread::hardware_concurrency();
    std::vector<std::thread> threads;
    threads.resize(nThreads);

    std::vector<double> partialResults;
    partialResults.resize(nThreads);

    std::vector<double> maxValueF;
    maxValueF.resize(nThreads);

    std::size_t nPointEvalPerThread = m_nRandomPoints / nThreads;

    auto func = [&](std::size_t const &startIndex, std::size_t const &endIndex,
    		std::size_t const &threadIndex)
    {
         double tmpValue;
         double tmpValueF;
         std::vector<std::vector<double>> evalPoint(m_nDirections);

         for(std::size_t i = startIndex; i < endIndex; i++)
         {
             for(std::size_t direction = 0; direction < m_nDirections; direction++)
             {
            	 evalPoint[direction].resize(m_dimension[direction]);
            	 for(std::size_t dimension = 0; dimension < m_dimension[direction]; dimension++)
            	 {
            		 evalPoint[direction][dimension] = m_evalPoints[direction][i][dimension];
            	 }
             }
            //tmpValue = std::abs(f(m_evalPoints[i]) - smolyakMLSApprox.at(m_evalPoints[i],f));
             tmpValue = std::abs(f(evalPoint) - approx.at(evalPoint,f));
             tmpValueF = std::abs(f(evalPoint));

            if(partialResults[threadIndex] < tmpValue)
            {
                partialResults[threadIndex] = tmpValue;
            }
            if(maxValueF[threadIndex] < tmpValueF)
            {
            	maxValueF[threadIndex] = tmpValueF;
            }
        }
    };


    std::size_t startIdx = 0;
    std::size_t endIdx   = nPointEvalPerThread;
    // Threads work
    for(std::size_t i = 0; i < nThreads - 1; i++)
    {
        partialResults[i] = 0.0;
        maxValueF[i] = 0.0;

        threads[i] = std::thread(func, startIdx, endIdx, i);

        startIdx += nPointEvalPerThread;
        endIdx   += nPointEvalPerThread;
    }

    startIdx = (nThreads - 1) * nPointEvalPerThread;
    endIdx   = m_evalGrid.size();

    partialResults[nThreads - 1] = 0.0;
    maxValueF[nThreads - 1] = 0.0;

    threads[nThreads - 1] = std::thread(func, startIdx, endIdx, nThreads - 1);



    for(std::size_t i = 0; i < nThreads; i++)
    {
        threads[i].join();
    }

    double maxValueFTogether = 0.0;
    for(std::size_t i = 0; i < nThreads; i++)
    {
        if(m_error < partialResults[i])
        {
            m_error = partialResults[i];
        }
        if(maxValueFTogether < maxValueF[i])
        {
        	maxValueFTogether = maxValueF[i];
        }
    }

    return m_error / maxValueFTogether;
}



double RandomPointsError::computeError_sequential(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
												  Approximation& approx)
{
    m_error = 0.0;

    double tmp = 0.0;

    std::vector<std::vector<double>> evalPoint (m_nDirections);


    for(std::size_t i = 0; i < m_nRandomPoints; i++)
    {
    	for (std::size_t direction = 0; direction < m_nDirections; direction++)
    	{
    		evalPoint[direction].resize (m_dimension[direction]);
    		for (std::size_t dimension = 0; dimension < m_dimension[direction]; dimension++)
    		{
    			evalPoint[direction][dimension] =
    					m_evalPoints[direction][i][dimension];
    		}
    	}
        //tmp = std::abs(f(m_evalPoints[i]) - smolyakMLSApprox.at(m_evalPoints[i],f));
        tmp = std::abs(f(evalPoint) - approx.at(evalPoint,f));

        if(m_error < tmp)
        {
            m_error = tmp;
        }
    }

    return m_error;
}


double RandomPointsError::computeRelativeError_sequential(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
														  Approximation& approx)
{
    m_error = 0.0;

    double tmp = 0.0;
    double maxValueF;
    double tmpFValue = 0.0;
    std::vector<std::vector<double>> evalPoint (m_nDirections);


    for(std::size_t i = 0; i < m_nRandomPoints; i++)
    {
    	for (std::size_t direction = 0; direction < m_nDirections; direction++)
    	{
    		evalPoint[direction].resize (m_dimension[direction]);
    		for (std::size_t dimension = 0; dimension < m_dimension[direction]; dimension++)
    		{
    			evalPoint[direction][dimension] =
    					m_evalPoints[direction][i][dimension];
    		}
    	}
        //tmp = std::abs(f(m_evalPoints[i]) - smolyakMLSApprox.at(m_evalPoints[i],f));
        tmp = std::abs(f(evalPoint) - approx.at(evalPoint,f));
        tmpFValue = std::abs(f(evalPoint));

        if(m_error < tmp)
        {
            m_error = tmp;
        }
        if(maxValueF < tmpFValue)
        {
        	maxValueF = tmpFValue;
        }
    }

    return m_error / maxValueF;
}
*/
