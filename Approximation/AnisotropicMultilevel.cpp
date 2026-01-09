/*
 * AnisotropicMultilevel.cpp
 *
 *  Created on: 17.06.2025
 *  Last changed: 05.11.2025
 *      Author: kempf
 *
 *      Implementation of class AnisotropicMultilevel
 */


//TO USE:
/*
int main() {
    std::unique_ptr<WendlandKernel> kernel_proto = std::make_unique<WendlandC0>();

    MultilevelSolver solver(std::move(kernel_proto))
*/


#include "AnisotropicMultilevel.h"


AnisotropicMultilevel::AnisotropicMultilevel(MultilevelParams &params)
{
	m_motherKernel = std::move(params.motherKernel);
	m_centers = std::move(params.centers);
	m_nLevels = m_centers->size();
	m_indices.resize(m_nLevels);
	m_points.resize(m_nLevels);
	m_dimension = (*m_centers)[0][0].size();
	m_nPoints.resize(m_nLevels);
	m_scalingMatricesPairs.resize(m_nLevels);
	m_trees.resize(m_nLevels);
	m_coefficients.resize(m_nLevels);
	m_trees.resize(m_nLevels);
	m_residuals.resize(m_nLevels);

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

		auto kernel = m_motherKernel->clone();

		PCA principleComponents(m_points[level]);
		std::vector<std::vector<double>> Q = principleComponents.getPrincipalComponentsSTD();
		//std::vector<double> Dtmp = principleComponents.getEigenvaluesSTD();
		std::vector<double> Dtmp = principleComponents.getNormalizedEigenvaluesMin();

		std::cout << "Eigenvalues of level " << level << " from PCA:" << std::endl;
		std::vector<double> D(m_dimension);
		for(std::size_t i = 0; i < Dtmp.size(); i++)
		{
			std::cout << Dtmp[i] << std::endl;
			D[i] = Dtmp[i] * params.overlap;
		}

		std::cout << "Q:" << std::endl;
		for(std::size_t i = 0; i < Dtmp.size(); i++)
		{
			for(std::size_t j = 0; j < Dtmp.size(); j++)
			{
				std::cout << Q[i][j] << " ";
			}
			std::cout << std::endl;
		}

		std::cout << std::endl;

		m_scalingMatricesPairs[level] = std::make_pair(Q,D);

		kernel->setSupport(m_scalingMatricesPairs[level]);
		kernel->setDimension(m_dimension);
		m_kernels.emplace_back(std::move(kernel));
	}

	for(std::size_t level = 0; level < m_nLevels-1; level++)
	{
		m_residuals[level].resize(m_nPoints[level+1]);
	}

	setUpLevels();
}


AnisotropicMultilevel::AnisotropicMultilevel(double const &overlap,
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
	m_scalingMatricesPairs.resize(m_nLevels);
	m_trees.resize(m_nLevels);
	m_coefficients.resize(m_nLevels);
	m_trees.resize(m_nLevels);
	m_residuals.resize(m_nLevels);

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

		auto kernel = m_motherKernel->clone();

		PCA principleComponents(m_points[level]);
		std::vector<std::vector<double>> Q = principleComponents.getPrincipalComponentsSTD();
		//std::vector<double> Dtmp = principleComponents.getEigenvaluesSTD();
		std::vector<double> Dtmp = principleComponents.getNormalizedEigenvaluesMin();

		std::cout << "Eigenvalues of level " << level << " from PCA:" << std::endl;
		std::vector<double> D(m_dimension);
		for(std::size_t i = 0; i < Dtmp.size(); i++)
		{
			std::cout << Dtmp[i] << std::endl;
			D[i] = Dtmp[i] * overlap;
		}

		std::cout << "Q:" << std::endl;
		for(std::size_t i = 0; i < Dtmp.size(); i++)
		{
			for(std::size_t j = 0; j < Dtmp.size(); j++)
			{
				std::cout << Q[i][j] << " ";
			}
			std::cout << std::endl;
		}

		std::cout << std::endl;

		m_scalingMatricesPairs[level] = std::make_pair(Q,D);

		kernel->setSupport(m_scalingMatricesPairs[level]);
		kernel->setDimension(m_dimension);
		m_kernels.emplace_back(std::move(kernel));
	}

	for(std::size_t level = 0; level < m_nLevels-1; level++)
	{
		m_residuals[level].resize(m_nPoints[level+1]);
	}

	setUpLevels();
}


AnisotropicMultilevel::AnisotropicMultilevel(double const &overlap,
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
	m_scalingMatricesPairs.resize(m_nLevels);
	m_trees.resize(m_nLevels);
	m_coefficients.resize(m_nLevels);
	m_trees.resize(m_nLevels);
	m_residuals.resize(m_nLevels);

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

		auto kernel = m_motherKernel->clone();

		PCA principleComponents(m_points[level]);
		std::vector<std::vector<double>> Q = principleComponents.getPrincipalComponentsSTD();
		//std::vector<double> Dtmp = principleComponents.getEigenvaluesSTD();
		std::vector<double> Dtmp = principleComponents.getNormalizedEigenvaluesMin();

		std::cout << "Eigenvalues of level " << level << " from PCA:" << std::endl;
		std::vector<double> D(m_dimension);
		for(std::size_t i = 0; i < Dtmp.size(); i++)
		{
			std::cout << Dtmp[i] << std::endl;
			D[i] = Dtmp[i] * overlap;
		}

		std::cout << "Q:" << std::endl;
		for(std::size_t i = 0; i < Dtmp.size(); i++)
		{
			for(std::size_t j = 0; j < Dtmp.size(); j++)
			{
				std::cout << Q[i][j] << " ";
			}
			std::cout << std::endl;
		}

		std::cout << std::endl;

		m_scalingMatricesPairs[level] = std::make_pair(Q,D);

		kernel->setSupport(m_scalingMatricesPairs[level]);
		kernel->setDimension(m_dimension);
		m_kernels.emplace_back(std::move(kernel));
	}

	for(std::size_t level = 0; level < m_nLevels-1; level++)
	{
		m_residuals[level].resize(m_nPoints[level+1]);
	}

	setUpLevels();
}

AnisotropicMultilevel::AnisotropicMultilevel(double const &overlap,
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
	m_scalingMatricesPairs.resize(m_nLevels);
	m_trees.resize(m_nLevels);
	m_coefficients.resize(m_nLevels);
	m_trees.resize(m_nLevels);
	m_residuals.resize(m_nLevels);

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

		auto kernel = m_motherKernel->clone();

		PCA principleComponents(m_points[level]);
		std::vector<std::vector<double>> Q = principleComponents.getPrincipalComponentsSTD();
		//std::vector<double> Dtmp = principleComponents.getEigenvaluesSTD();
		std::vector<double> Dtmp = principleComponents.getNormalizedEigenvaluesMin();

		std::cout << "Eigenvalues of level " << level << " from PCA:" << std::endl;
		std::vector<double> D(m_dimension);
		for(std::size_t i = 0; i < Dtmp.size(); i++)
		{
			std::cout << Dtmp[i] << std::endl;
			D[i] = Dtmp[i] * overlap;
		}

		std::cout << "Q:" << std::endl;
		for(std::size_t i = 0; i < Dtmp.size(); i++)
		{
			for(std::size_t j = 0; j < Dtmp.size(); j++)
			{
				std::cout << Q[i][j] << " ";
			}
			std::cout << std::endl;
		}

		std::cout << std::endl;

		m_scalingMatricesPairs[level] = std::make_pair(Q,D);

		kernel->setSupport(m_scalingMatricesPairs[level]);
		kernel->setDimension(m_dimension);
		m_kernels.emplace_back(std::move(kernel));
	}

	for(std::size_t level = 0; level < m_nLevels-1; level++)
	{
		m_residuals[level].resize(m_nPoints[level+1]);
	}

	setUpLevels();
}

AnisotropicMultilevel::AnisotropicMultilevel(std::vector<double> const &overlap,
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
	m_scalingMatricesPairs.resize(m_nLevels);
	m_trees.resize(m_nLevels);
	m_coefficients.resize(m_nLevels);
	m_trees.resize(m_nLevels);
	m_residuals.resize(m_nLevels);

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

		auto kernel = m_motherKernel->clone();

		PCA principleComponents(m_points[level]);
		std::vector<std::vector<double>> Q = principleComponents.getPrincipalComponentsSTD();
		//std::vector<double> Dtmp = principleComponents.getEigenvaluesSTD();
		std::vector<double> Dtmp = principleComponents.getNormalizedEigenvaluesMin();

		std::cout << "Eigenvalues of level " << level << " from PCA:" << std::endl;
		std::vector<std::vector<double>> D(m_nLevels, std::vector<double>(m_dimension));
		for(std::size_t levels = 0; levels < m_nLevels; levels++)
		{
			for(std::size_t i = 0; i < Dtmp.size(); i++)
			{
			//std::cout << Dtmp[i] << std::endl;
				D[levels][i] = Dtmp[i] * overlap[levels];

			}
		}
		//std::cout << std::endl;

		m_scalingMatricesPairs[level] = std::make_pair(Q,D[level]);

		kernel->setSupport(m_scalingMatricesPairs[level]);
		kernel->setDimension(m_dimension);
		m_kernels.emplace_back(std::move(kernel));
	}

	for(std::size_t level = 0; level < m_nLevels-1; level++)
	{
		m_residuals[level].resize(m_nPoints[level+1]);
	}

	setUpLevels();
}


AnisotropicMultilevel::AnisotropicMultilevel(
		std::vector<std::pair<std::vector<std::vector<double>>,std::vector<double>>> const &scalingPairs,
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
	m_scalingMatricesPairs.resize(m_nLevels);
	m_trees.resize(m_nLevels);
	m_coefficients.resize(m_nLevels);
	m_trees.resize(m_nLevels);
	m_residuals.resize(m_nLevels);

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
		auto kernel = m_motherKernel->clone();

		std::vector<std::vector<double>> Q = scalingPairs[level].first;
		std::vector<double> D = scalingPairs[level].second;

		m_scalingMatricesPairs[level] = std::make_pair(Q,D);

		kernel->setSupport(m_scalingMatricesPairs[level]);
		kernel->setDimension(m_dimension);
		m_kernels.emplace_back(std::move(kernel));
	}

	for(std::size_t level = 0; level < m_nLevels-1; level++)
	{
		m_residuals[level].resize(m_nPoints[level+1]);
	}

	setUpLevels();
}



double AnisotropicMultilevel::at(std::vector<double> const &x) const
{
	double value = 0.0;

	for(std::size_t level = 0; level < m_nLevels; level++)
	{
		value += atLevel(level, x);
	}
	return value;
}


double AnisotropicMultilevel::atUpToLevel
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



void AnisotropicMultilevel::setUpLevels()
{
	for(std::size_t level = 0; level < m_nLevels; level++)
	{
		std::vector<std::vector<std::size_t>> rows(m_nPoints[level]);
		rows = m_trees[level]->multipleRangeSearchs(m_points[level],
													m_indices[level],
													m_scalingMatricesPairs[level]);

		//SparseMatrix A(m_nPoints[level], m_nPoints[level]);

		std::size_t numberNonEmptyEntries = 0;
	    for(std::size_t i = 0; i < m_nPoints[level]; i++)
	    {
	        numberNonEmptyEntries += rows[i].size();
	    }

	    std::vector<std::size_t> rowOffsets(m_nPoints[level] + 1, 0);

	    // Step 1: Count entries per row
	    for (std::size_t i = 0; i < m_nPoints[level]; ++i)
	        rowOffsets[i + 1] = rows[i].size();

	    // Step 2: Exclusive prefix sum to compute row offsets
	    for (std::size_t i = 1; i <= m_nPoints[level]; ++i)
	        rowOffsets[i] += rowOffsets[i - 1];

	    std::size_t nnz = rowOffsets[m_nPoints[level]];
	    std::vector<std::size_t> colIndices(nnz);
	    std::vector<double> values(nnz);

	    // Step 3: Parallel fill kernel values and col indices
	    pFor([&](std::size_t i) {
	        std::size_t offset = rowOffsets[i];
	        for (std::size_t j = 0; j < rows[i].size(); ++j) {
	            std::size_t col = rows[i][j];
	            colIndices[offset + j] = col;
	            values[offset + j] = m_kernels[level]->at(m_points[level][i], m_points[level][col]);
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


void AnisotropicMultilevel::solve(std::function<double(const std::vector<double>&)> const &targetFunc)
{
	{
	std::vector<double> RHS(m_nPoints[0]);
	for(std::size_t i = 0; i < m_nPoints[0]; i++)
	{
		RHS[i] = targetFunc(m_points[0][i]);
	}

	m_coefficients[0] = CG(m_kernelMatrices[0],RHS);

	if(m_nLevels != 1)
	{
		updateResiduals(0,targetFunc);
	}
	}

	for(std::size_t level = 1; level < m_nLevels; level++)
	{
		std::vector<double> RHS;
		RHS = m_residuals[level-1];

		m_coefficients[level] = CG(m_kernelMatrices[level], RHS);

		if(level < m_nLevels-1)
		{
			updateResiduals(level,targetFunc);
		}
	}
}


void AnisotropicMultilevel::solve(std::vector<double> const &targetFunc)
{
	std::cout << "Multilevel Method with non-nested data does not work if given only one data vector"
			<< std::endl;
}



double AnisotropicMultilevel::atLevel(std::size_t level, std::vector<double> const &x) const
{
	double value = 0.0;

	std::vector<std::size_t> neighbourIndices =
			m_trees[level]->rangeSearch(x,m_scalingMatricesPairs[level]);

	for(std::size_t i = 0; i < neighbourIndices.size(); i++)
	{
		value += m_coefficients[level][neighbourIndices[i]]
						* m_kernels[level]->at(x,m_points[level][neighbourIndices[i]]);
	}

	//std::cout << "atLevel value: " << value << std::endl;
	return value;
}


void AnisotropicMultilevel::updateResiduals(std::size_t const &level,
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
