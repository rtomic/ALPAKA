/*
 * LooError.cpp
 *
 *  Created on:   4.09.2019
 *  Last changed: 19.08.2025
 *
 *      Author: kempf
 */



#include "LooError.h"


LooError::LooError(std::vector<std::vector<double>> errorSites)
{
	m_grid = errorSites;

    //this->create();
}

LooError::LooError(std::shared_ptr<std::vector<std::vector<std::vector<double>>>> tensorErrorSites)
{
	m_tensorGrid_ptr = tensorErrorSites;
}



double LooError::getError()
{
	return m_error;
}

double LooError::getRelativeError()
{
	return m_relativeError;
}

void LooError::computeError(std::function<double(std::vector<double> const &evalPoint)> const &f,
		std::unique_ptr<Approximation> const &approximation)
{
	std::size_t nPoints = m_grid.size();
	std::size_t nThreads = std::thread::hardware_concurrency();
	std::vector<std::thread> threads;
	std::size_t chunkSize = (nPoints + nThreads - 1) / nThreads;

	double globalMaxTargetValue;

    std::mutex mutex;

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
	}

	m_relativeError = m_error / globalMaxTargetValue;
}

void LooError::computeError(std::function<double(std::vector<double> const &evalPoint)> &f,
		std::unique_ptr<Approximation> const &approximation)
{    
	std::size_t nPoints = m_grid.size();
	std::size_t nThreads = std::thread::hardware_concurrency();
	std::vector<std::thread> threads;
	std::size_t chunkSize = (nPoints + nThreads - 1) / nThreads;

	double globalMaxTargetValue;
    
    std::mutex mutex;

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
	}

	m_relativeError = m_error / globalMaxTargetValue;
}


void LooError::computeError(std::function<double(std::vector<std::vector<double>> const &evalPoint)> const &f,
		CombinationTechniqueApproximation const &approximation)
{
	std::size_t nDirections = m_tensorGrid_ptr->size();
	std::size_t nPoints = 1;
	for(std::size_t i = 0; i < nDirections; i++)
	{
		nPoints *= (*m_tensorGrid_ptr)[i].size();
	}
	std::size_t nThreads = std::thread::hardware_concurrency();
	std::vector<std::thread> threads;
	std::size_t chunkSize = (nPoints + nThreads - 1) / nThreads;

	double globalMaxTargetValue;

    std::mutex mutex;

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
				std::vector<std::vector<double>> pt(nDirections);
				for(std::size_t j = 0; j < nDirections; j++)
				{
					pt[i] = (*m_tensorGrid_ptr)[j][i];
				}
				double approxVal = approximation.at(pt,f);
				double targetVal = f(pt);
				double absError = std::abs(approxVal - targetVal);
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
	}

	m_relativeError = m_error / globalMaxTargetValue;
}

/*
    std::size_t nThreads = thread::hardware_concurrency();
    std::vector<std::thread> threads;
    threads.resize(nThreads);  
    
    std::vector<double> partialResults;
    partialResults.resize(nThreads);
    
    std::size_t nPointEvalPerThread = m_grid.size() / nThreads;
    
    auto func = [&](std::size_t const &startIndex, std::size_t const &endIndex, std::size_t const &threadIndex)
    {
        partialResults[threadIndex] = 0.0;
        
        double tmp = 0.0;
        
        for(std::size_t i = startIndex; i < endIndex; i++)
        {

            tmp = std::fabs(f(m_grid[i]) - approximation->at(m_grid[i]));
            
            if(partialResults[threadIndex] < tmp)
            {
                partialResults[threadIndex] = tmp;
            }
        }
        
    };
    
    std::size_t startIdx = 0;
    std::size_t endIdx = nPointEvalPerThread;
    
    for(std::size_t i = 0; i < nThreads - 1; i++)
    {
        threads[i] = std::thread(func, startIdx, endIdx, i);
        
        startIdx += nPointEvalPerThread;
        endIdx   += nPointEvalPerThread;
    }
    startIdx = (nThreads - 1) * nPointEvalPerThread;
    endIdx   = m_grid.size();
    threads[nThreads - 1] = std::thread (func, startIdx, endIdx, nThreads - 1);
    
    for(std::size_t i = 0; i < nThreads; i++)
    {
        threads[i].join();
    }
    
    for(std::size_t i = 0; i < partialResults.size(); i++)
    {
        if(m_error < partialResults[i])
        {
            m_error = partialResults[i];
        }
    }
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "Evaluating Loo-Error on " << m_grid.size() << " points" << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
    
    return m_error;

}
*/

/**************************************************************************************
 *
 *************************************************************************************/
/*
void LooError::create()
{
	std::vector<std::vector<double>> allOneDPoints;
	allOneDPoints.resize(m_dimension);

	for(std::size_t i = 0; i < m_dimension; i++)
	{
		allOneDPoints[i].resize(m_n1DPoints);
		for(std::size_t j = 0; j < m_n1DPoints; j++)
		{
			allOneDPoints[i][j] = m_domain->lo[i] + j * m_stepSize;
		}

		//allOneDPoints.push_back(oneDPoints);
	}
	m_grid = cartesian_product(allOneDPoints);
}
*/
/**************************************************************************************
 *
 *************************************************************************************/
/*
std::vector<std::vector<double>> LooError::cartesian_product(const std::vector<std::vector<double>>& vector) const
{
  std::vector<std::vector<double>> s = {{}};
  for(auto& u : vector)
  {
    std::vector<std::vector<double>> r;
    for(auto& x : s)
    {
      for(auto y : u)
      {
        r.push_back(x);
        r.back().push_back(y);
      }
    }
    s.swap(r);
  }

  return s;
}
*/


