/*
 * Multilevel.cpp
 *
 *  Created on: 16.06.2025
 *  Last changed: 20.08.2025
 *      Author: kempf
 *
 *      Implementation of class Multilevel
 */


//TO USE:
/*
int main() {
    std::unique_ptr<WendlandKernel> kernel_proto = std::make_unique<WendlandC0>();

    MultilevelSolver solver(std::move(kernel_proto))
*/


#include "DynamicDiscardingMultilevel.h"


#include <mutex>

DynamicDiscardingMultilevel::DynamicDiscardingMultilevel(MultilevelParams &params)
{
	m_motherKernel = std::move(params.motherKernel);
	m_centers = std::move(params.centers);
	m_nLevels = m_centers->size();
	m_indices.resize(m_nLevels);
	m_points.resize(m_nLevels);
	m_dimension = (*m_centers)[0][0].size();
	m_nPoints.resize(m_nLevels);
	m_radius.resize(m_nLevels);
	m_trees.resize(m_nLevels);
	m_coefficients.resize(m_nLevels);
	m_trees.resize(m_nLevels);
	m_residuals.resize(m_nLevels);

	double epsilon = 1e-6;

	for(std::size_t level = 0; level < m_nLevels; level++)
	{
		m_points[level] = (*m_centers)[level];
		m_nPoints[level] = m_points[level].size();
		m_indices[level].resize(m_nPoints[level]);
		std::iota(m_indices[level].begin(), m_indices[level].end(), 0);

		m_trees[level] = std::make_unique<KDtree>(m_points[level],
												  m_indices[level],
												  m_nPoints[level],
												  m_dimension);
		m_radius[level] = params.overlap * sqrt(m_dimension) * std::pow(static_cast<double>(m_nPoints[level]),
											 -1.0/static_cast<double>(m_dimension));
		double radius = m_radius[level];

		auto kernel = m_motherKernel->clone();
		kernel->setSupport(radius);
		kernel->setDimension(m_dimension);
		m_kernels.emplace_back(std::move(kernel));

		m_threshold[level] = epsilon * std::pow(radius, 0.5 * m_dimension);
	}

	for(std::size_t level = 0; level < m_nLevels-1; level++)
	{
		m_residuals[level].resize(m_nPoints[level+1]);
	}
	setUpLevels();
}

DynamicDiscardingMultilevel::DynamicDiscardingMultilevel(double const &overlap,
					   std::unique_ptr<std::vector<std::vector<std::vector<double>>>> centers,
					   std::unique_ptr<Kernel> motherKernel)
{
	m_motherKernel = std::move(motherKernel);
	m_centers = std::move(centers);
	m_nLevels = m_centers->size();
	m_indices.resize(m_nLevels);
	m_points.resize(m_nLevels);
	m_dimension = (*m_centers)[0][0].size();
	m_nPoints.resize(m_nLevels);
	m_radius.resize(m_nLevels);
	m_trees.resize(m_nLevels);
	m_coefficients.resize(m_nLevels);
	m_trees.resize(m_nLevels);
	m_residuals.resize(m_nLevels);

	double epsilon = 1e-6;

	for(std::size_t level = 0; level < m_nLevels; level++)
	{
		m_points[level] = (*m_centers)[level];
		m_nPoints[level] = m_points[level].size();
		m_indices[level].resize(m_nPoints[level]);
		std::iota(m_indices[level].begin(), m_indices[level].end(), 0);

		m_trees[level] = std::make_unique<KDtree>(m_points[level],
												  m_indices[level],
												  m_nPoints[level],
												  m_dimension);
		m_radius[level] = overlap * sqrt(m_dimension) * std::pow(static_cast<double>(m_nPoints[level]),
											 -1.0/static_cast<double>(m_dimension));
		double radius = m_radius[level];

		auto kernel = m_motherKernel->clone();
		kernel->setSupport(radius);
		kernel->setDimension(m_dimension);
		m_kernels.emplace_back(std::move(kernel));

		m_threshold[level] = epsilon * std::pow(radius, 0.5 * m_dimension);
	}

	for(std::size_t level = 0; level < m_nLevels-1; level++)
	{
		m_residuals[level].resize(m_nPoints[level+1]);
	}
	setUpLevels();
}

DynamicDiscardingMultilevel::DynamicDiscardingMultilevel(double const &overlap,
					   std::shared_ptr<std::vector<std::vector<std::vector<double>>>> centers,
					   std::unique_ptr<Kernel> motherKernel)
{
	m_motherKernel = std::move(motherKernel);
	m_centers_shared = centers;
	m_nLevels = m_centers_shared->size();
	m_indices.resize(m_nLevels);
	m_points.resize(m_nLevels);
	m_dimension = (*m_centers_shared)[0][0].size();
	m_nPoints.resize(m_nLevels);
	m_radius.resize(m_nLevels);
	m_trees.resize(m_nLevels);
	m_coefficients.resize(m_nLevels);
	m_trees.resize(m_nLevels);
	m_residuals.resize(m_nLevels);
	m_threshold.resize(m_nLevels);
	m_actuallyUsedIndices.resize(m_nLevels);

	double epsilon = 1e-6;

	for(std::size_t level = 0; level < m_nLevels; level++)
	{
		m_points[level] = (*m_centers_shared)[level];
		m_nPoints[level] = m_points[level].size();
		m_indices[level].resize(m_nPoints[level]);
		std::iota(m_indices[level].begin(), m_indices[level].end(), 0);

		m_trees[level] = std::make_unique<KDtree>(m_points[level],
												  m_indices[level],
												  m_nPoints[level],
												  m_dimension);
		m_radius[level] = overlap * sqrt(m_dimension) * std::pow(static_cast<double>(m_nPoints[level]),
											 -1.0/static_cast<double>(m_dimension));
		double radius = m_radius[level];

		auto kernel = m_motherKernel->clone();
		kernel->setSupport(radius);
		kernel->setDimension(m_dimension);

		m_kernels.emplace_back(std::move(kernel));

		m_threshold[level] = epsilon * std::pow(radius, 0.5 * m_dimension);
	}

	for(std::size_t level = 0; level < m_nLevels-1; level++)
	{
		m_residuals[level].resize(m_nPoints[level+1]);
	}
	setUpLevels();
}


DynamicDiscardingMultilevel::DynamicDiscardingMultilevel(double const &overlap,
					   std::vector<std::vector<std::vector<double>>> const &centers,
					   std::unique_ptr<Kernel> motherKernel)
{
	m_motherKernel = std::move(motherKernel);
	m_centers_vec = centers;
	m_nLevels = m_centers_vec.size();
	m_indices.resize(m_nLevels);
	m_points.resize(m_nLevels);
	m_dimension = m_centers_vec[0][0].size();
	m_nPoints.resize(m_nLevels);
	m_radius.resize(m_nLevels);
	m_trees.resize(m_nLevels);
	m_coefficients.resize(m_nLevels);
	m_trees.resize(m_nLevels);
	m_residuals.resize(m_nLevels);
	m_threshold.resize(m_nLevels);
	m_actuallyUsedIndices.resize(m_nLevels);

	double epsilon = 1e-6;

	for(std::size_t level = 0; level < m_nLevels; level++)
	{
		m_points[level] = m_centers_vec[level];
		m_nPoints[level] = m_points[level].size();
		m_indices[level].resize(m_nPoints[level]);
		std::iota(m_indices[level].begin(), m_indices[level].end(), 0);

		m_trees[level] = std::make_unique<KDtree>(m_points[level],
												  m_indices[level],
												  m_nPoints[level],
												  m_dimension);
		m_radius[level] = overlap * sqrt(m_dimension) * std::pow(static_cast<double>(m_nPoints[level]),
											 -1.0/static_cast<double>(m_dimension));
		double radius = m_radius[level];

		auto kernel = m_motherKernel->clone();
		kernel->setSupport(radius);
		kernel->setDimension(m_dimension);

		m_kernels.emplace_back(std::move(kernel));

		m_threshold[level] = epsilon * std::pow(radius, 0.5 * m_dimension);
	}

	for(std::size_t level = 0; level < m_nLevels-1; level++)
	{
		m_residuals[level].resize(m_nPoints[level+1]);
	}
	setUpLevels();
}

double DynamicDiscardingMultilevel::at(std::vector<double> const &x) const
{
	double value = 0.0;

	for(std::size_t level = 0; level < m_nLevels; level++)
	{
		value += atLevel(level, x);
	}

	return value;
}


double DynamicDiscardingMultilevel::atUpToLevel
		(std::vector<double> const &x, std::size_t const &givenLevel) const
{
	double value = 0.0;

	if(givenLevel >= m_nLevels)
	{
		std::cerr << "Given level is too high" << std::endl;
		return -100;
	}

	for(std::size_t level = 0; level <= givenLevel; level++)
	{
		value += atLevel(level,x);
	}
	return value;
}

void DynamicDiscardingMultilevel::setUpLevels()
{
	double kernelValue;
	for(std::size_t level = 0; level < m_nLevels; level++)
	{
		double radius = m_radius[level];
		std::vector<std::vector<std::size_t>> rows(m_nPoints[level]);
		rows = m_trees[level]->multipleRangeSearchs(m_points[level],
													m_indices[level],
													radius);

		SparseMatrix A(m_nPoints[level], m_nPoints[level]);

		std::size_t numberNonEmptyEntries = 0;
	    for(std::size_t i = 0; i < m_nPoints[level]; i++)
	    {
	        numberNonEmptyEntries += rows[i].size();
	    }

	    std::vector<std::tuple<size_t, size_t, double>> triplets(numberNonEmptyEntries);
	    std::size_t counter = 0;

	    for(std::size_t i = 0; i < m_nPoints[level]; i++)
	    {
	    	for(std::size_t j = 0; j < rows[i].size(); j++)
	    	{
	    		kernelValue = m_kernels[level]->at(m_points[level][i], m_points[level][rows[i][j]]);
	    		triplets[counter] = {i,rows[i][j],kernelValue};
	    		counter++;
	    	}
	    }
	    A.setFromTriplets(triplets);
	    m_kernelMatrices.emplace_back(std::move(A));

		if(CParam::getFirstValue<std::size_t>("MultilevelStatistics") == 0)
		{
			std::cout << "On " << level;
			std::cout << "we have " << m_nPoints[level] << " points" << std::endl;
			std::size_t nmax = 0;
			std::size_t naverag = 0;
			for (std::size_t i = 0; i < m_nPoints[level]; i++)
			{
				if(nmax > rows[i].size())
				{
					nmax = rows[i].size();
				}
				naverag += rows[i].size();
			}
			std::cout << "Level :" << level << " has neighbors (max/av)" << nmax << "/"
		    << naverag / m_nPoints[level] << std::endl;

		}
	}
}


void DynamicDiscardingMultilevel::solve(std::function<double(const std::vector<double> &)> const &targetFunc)
{
	{
	std::vector<double> RHS(m_nPoints[0]);
	for(std::size_t i = 0; i < m_nPoints[0]; i++)
	{
		RHS[i] = targetFunc(m_points[0][i]);
	}

	std::vector<double> tmp_coefficients = CG(m_kernelMatrices[0],RHS);

	for(std::size_t i = 0; i < m_nPoints[0]; i++)
	{
		if(std::fabs(tmp_coefficients[i]) > m_threshold[0])
		{
			m_coefficients[0].push_back(tmp_coefficients[i]);
			m_actuallyUsedIndices[0].push_back(i);
		}
	}

	std::cout << "Level: " << 0 << "actually using " << m_actuallyUsedIndices[0].size() << " points"
			<< std::endl;

	updateResiduals(0,targetFunc);
	}

	for(std::size_t level = 1; level < m_nLevels; level++)
	{
		std::vector<double> RHS;
		RHS = m_residuals[level-1];

		std::vector<double> tmp_coefficients = CG(m_kernelMatrices[level], RHS);

		for(std::size_t i = 0; i < m_nPoints[level]; i++)
		{
			if (std::fabs(tmp_coefficients[i]) > m_threshold[level])
			{
				m_coefficients[level].push_back(tmp_coefficients[i]);
				m_actuallyUsedIndices[level].push_back(i);
			}
		}

		std::cout << "Level: " << level << "actually using " << m_actuallyUsedIndices[level].size()
				<< " points" << std::endl;

		if(level < m_nLevels-1)
		{
			updateResiduals(level,targetFunc);
		}
	}
}


void DynamicDiscardingMultilevel::solve(std::vector<double> const &targetFunc)
{
	std::cout << "Multilevel Method with non-nested data does not work if given only one data vector"
			<< std::endl;
}



double DynamicDiscardingMultilevel::atLevel(std::size_t level, std::vector<double> const &x) const
{
	double value = 0.0;

	double radius = m_radius[level];
	std::vector<std::size_t> neighbourIndices = m_trees[level]->rangeSearch(x, radius);

	//Check if neighbourIndices appear in m_actuallyUsedIndices[level]
	std::vector<size_t> indices;
	indices.reserve(std::min(neighbourIndices.size(),m_actuallyUsedIndices[level].size()));

	for(std::size_t val : neighbourIndices)
	{
		if(std::binary_search(m_actuallyUsedIndices[level].begin(), m_actuallyUsedIndices[level].end(), val))
		{
			indices.push_back(val);
		}
	}

	for(std::size_t i = 0; i < indices.size(); i++)
	{
		value += m_coefficients[level][indices[i]]
				* m_kernels[level]->at(x,m_points[level][indices[i]]);
	}

	return value;
}


void DynamicDiscardingMultilevel::updateResiduals(std::size_t const &level,
								 std::function<double(std::vector<double> &)> const &targetFunc)
{
	auto func = [&](std::size_t i)
		{
			double approxValue = 0.0;
			for(std::size_t tmpLevel = 0; tmpLevel <= level; tmpLevel++)
			{
				approxValue += atLevel(tmpLevel,m_points[level+1][i]);
			}

			m_residuals[level][i] = targetFunc(m_points[level+1][i]) - approxValue;
		};

	parallelFor(func, 0, m_nPoints[level+1]);
}
