/*
 * nestedMultilevel.cpp
 *
 *  Created on: 23.05.2025
 *  Last changed: 19.11.2025
 *      Author: kempf
 *
 *      Implementation of class NestedMultilevel
 */


//TO USE:
/*
int main() {
    std::unique_ptr<WendlandKernel> kernel_proto = std::make_unique<WendlandC0>();

    MultilevelSolver solver(std::move(kernel_proto))
*/




#include "NestedMultilevel.h"

#include <unordered_set>
#include <mutex>



NestedMultilevel::NestedMultilevel(MultilevelParams &params)
{
		m_motherKernel = std::move(params.motherKernel);
		m_centers = std::move(params.nestedSites);
		m_indices = m_centers->getIndices();
		m_reverse = m_centers->getReverse();
		m_nLevels = m_centers->getNLevels();
		m_points = m_centers->getPoints();
		m_nPoints.resize(m_nLevels);
		m_dimension = m_points[0].size();
		m_radius.resize(m_nLevels);
		m_trees.resize(m_nLevels);
		m_coefficients.resize(m_nLevels);
		m_updateIndices.resize(m_nLevels);
		m_trees.resize(m_nLevels);


		for(std::size_t level = 0; level < m_nLevels; level++)
		{
			m_nPoints[level] = m_centers->getNPoints(level);

			m_trees[level] = std::make_unique<KDtree>(m_points, m_indices[level],m_nPoints[level],m_dimension);

			m_radius[level] = params.overlap * m_centers->getFillDistance()[level];

			auto kernel = m_motherKernel->clone();
			double radius = m_radius[level];
			kernel->setSupport(radius);
			kernel->setDimension(m_dimension);
			m_kernels.emplace_back(std::move(kernel));
		}
		m_funcValues.resize(m_nPoints[m_nLevels - 1]);

		setUpLevels();
}

NestedMultilevel::NestedMultilevel(double const &overlap, std::unique_ptr<NestedSitesGenerator> centers,
								   std::unique_ptr<Kernel> motherKernel)
{
	m_motherKernel = std::move(motherKernel);
	m_centers = std::move(centers);
	m_indices = m_centers->getIndices();
	m_reverse = m_centers->getReverse();
	m_nLevels = m_centers->getNLevels();
	m_points = m_centers->getPoints();
	m_nPoints.resize(m_nLevels);
	m_dimension = m_points[0].size();
	m_radius.resize(m_nLevels);
	m_trees.resize(m_nLevels);
	m_coefficients.resize(m_nLevels);
	m_updateIndices.resize(m_nLevels);
	m_trees.resize(m_nLevels);


	for(std::size_t level = 0; level < m_nLevels; level++)
	{
		m_nPoints[level] = m_centers->getNPoints(level);

		m_trees[level] = std::make_unique<KDtree>(m_points, m_indices[level],m_nPoints[level],m_dimension);

		m_radius[level] = overlap * m_centers->getFillDistance()[level];

		auto kernel = m_motherKernel->clone();
		double radius = m_radius[level];
		kernel->setSupport(radius);
		kernel->setDimension(m_dimension);
		m_kernels.emplace_back(std::move(kernel));
	}
	m_funcValues.resize(m_nPoints[m_nLevels - 1]);

	setUpLevels();
}


double NestedMultilevel::at(std::vector<double> const &x) const
{
	double value = 0.0;

	for(std::size_t level = 0; level < m_nLevels; level++)
	{
		value += atLevel(level, x);
	}

	return value;
}


double NestedMultilevel::atUpToLevel
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

void NestedMultilevel::setUpLevels()
{
	double radius;
	for(std::size_t level = 0; level < m_nLevels; level++)
	{

		m_updateIndices[level].resize(m_nPoints[m_nLevels-1]);

		radius = m_radius[level];
		m_updateIndices[level] = m_trees[level]->multipleRangeSearchs(m_points,
																	  m_indices[m_nLevels-1],
																	  //m_indices[level],
																	  radius);

		std::vector<std::vector<std::size_t>> rows;
		rows = m_trees[level]->multipleRangeSearchs(m_points,
													m_indices[level],
													radius);

		// Build global-to-local map for this level
		std::unordered_map<std::size_t, std::size_t> globalToLocal;
		globalToLocal.reserve(m_nPoints[level]);
		for (std::size_t local = 0; local < m_indices[level].size(); ++local) {
		    globalToLocal[m_indices[level][local]] = local;
		}

		std::vector<std::size_t> rowOffsets(m_nPoints[level] + 1, 0);

		// Step 1: Count entries per row (still fine to use global IDs here)
		for (std::size_t i = 0; i < m_nPoints[level]; ++i)
		    rowOffsets[i + 1] = rows[i].size();

		// Step 2: Exclusive prefix sum
		for (std::size_t i = 1; i <= m_nPoints[level]; ++i)
		    rowOffsets[i] += rowOffsets[i - 1];

		std::size_t nnz = rowOffsets[m_nPoints[level]];
		std::vector<std::size_t> colIndices(nnz);
		std::vector<double> values(nnz);

		// Step 3: Fill local colIndices and kernel values
		pFor([&](std::size_t i) {
		    std::size_t offset = rowOffsets[i];
		    for (std::size_t j = 0; j < rows[i].size(); ++j)
		    {
		        std::size_t globalCol = rows[i][j];
		        std::size_t localCol = globalToLocal.at(globalCol); // map global → local

		        colIndices[offset + j] = localCol;
		        values[offset + j] = m_kernels[level]->at(
		            m_points[m_indices[level][i]], // global point for row i
		            m_points[globalCol]            // global point for neighbor
		        );
		    }
		}, 0, m_nPoints[level]);

		// Step 4: Build matrix from CSR data
		SparseMatrix A(m_nPoints[level], m_nPoints[level]);
		A.buildFromCSRComponents(std::move(rowOffsets), std::move(colIndices), std::move(values));
		m_kernelMatrices.emplace_back(std::move(A));


	    std::cout << "On level " << level;
	    std::cout << " we have " << m_nPoints[level] << " points" << std::endl;
	    std::size_t nmax = 0;
	    std::size_t naverag = 0;
	    for (std::size_t i = 0; i < m_nPoints[level]; i++)
	    {
	    	if(nmax < rows[i].size())
	    	{
	    		nmax = rows[i].size();
	    	}
	    	naverag += rows[i].size();
	    }
	    std::cout << "Level: " << level << " has neighbors (max/av)" << nmax << "/"
		    << naverag / m_nPoints[level] << std::endl;
	}
}


void NestedMultilevel::solve(std::function<double(const std::vector<double> &)> const &targetFunc)
{
	for(std::size_t i = 0; i < m_nPoints[m_nLevels-1]; i++)
	{
		m_funcValues[i] = targetFunc(m_points[i]);
	}

	for(std::size_t level = 0; level < m_nLevels; level++)
	{
		std::cout << "Solving level " << level << std::endl;
		std::vector<double> RHS(m_nPoints[level]);
		for(std::size_t i = 0; i < m_nPoints[level]; i++)
		{
			RHS[i] = m_funcValues[m_indices[level][i]];
		}

		m_coefficients[level] = CG(m_kernelMatrices[level], RHS);

		if(level != m_nLevels-1)
		{
			updateResiduals(level);
		}
	}
}


void NestedMultilevel::solve(std::vector<double> const &targetFunc)
{
	m_funcValues = targetFunc;

	for(std::size_t level = 0; level < m_nLevels; level++)
	{
		std::cout << "Solving level " << level << std::endl;
		std::vector<double> RHS(m_nPoints[level]);
		for(std::size_t i = 0; i < m_nPoints[level]; i++)
		{
			RHS[i] = m_funcValues[m_indices[level][i]];
		}

		m_coefficients[level] = CG(m_kernelMatrices[level], RHS);

		if(level != m_nLevels-1)
		{
			updateResiduals(level);
		}

	}
}


double NestedMultilevel::atLevel(std::size_t level, std::vector<double> const &x) const
{
	double value = 0.0;

	double radius = m_radius[level];
	std::vector<std::size_t> neighbourIndices = m_trees[level]->rangeSearch(x, radius);

	for(std::size_t i = 0; i < neighbourIndices.size(); i++)
	{
		//value += m_coefficients[level](m_reverse[level][neighbourIndices[i]]) *
		value += m_coefficients[level][m_reverse[level][neighbourIndices[i]]] *
				m_kernels[level]->at(x,m_points[neighbourIndices[i]]);
	}
	return value;
}


void NestedMultilevel::updateResiduals(std::size_t const &level)
{
	auto func = [&](std::size_t i)
		{
			for(std::size_t j = 0; j < m_updateIndices[level][i].size(); j++)
			{
				m_funcValues[i] -= m_kernels[level]->at(m_points[i],m_points[m_updateIndices[level][i][j]])
						//* m_coefficients[level](m_reverse[level][m_updateIndices[level][i][j]]);
								* m_coefficients[level][m_reverse[level][m_updateIndices[level][i][j]]];

			}
		};

	parallelFor(func, 0, m_nPoints[m_nLevels-1]);
}

