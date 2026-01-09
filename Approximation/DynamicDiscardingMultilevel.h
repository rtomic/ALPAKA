/*
 * DynamicDiscardingMultilevel.h
 *
 *  Created on: 24.06.2025
 *  Last changed: 05.11.2025
 *      Author: kempf
 *
 *      Class for multilevel interpolation/approximation with general sites, using the dynamic discarding
 *      idea, i.e., ignoring coefficients smaller than a threshold
 */



#pragma once

#include <vector>
#include <functional>
#include <memory>
#include <atomic> // for safe counting of removed coefficients
#include <algorithm>

//#include <Eigen/Sparse>
//#include <Eigen/IterativeLinearSolvers>

#include "Approximation.h"
#include "WendlandKernel.h"
#include "Kernel.h"
#include "KDtree.h"
#include "LinearAlgebra.h"
#include "types.h"



class DynamicDiscardingMultilevel : public Approximation
{
	public:
		DynamicDiscardingMultilevel(MultilevelParams &params);
		DynamicDiscardingMultilevel(double const &overlap,
				   std::unique_ptr<std::vector<std::vector<std::vector<double>>>> centers,
				   std::unique_ptr<Kernel> motherKernel);

		DynamicDiscardingMultilevel(double const &overlap,
						   std::shared_ptr<std::vector<std::vector<std::vector<double>>>> centers,
						   std::unique_ptr<Kernel> motherKernel);

		DynamicDiscardingMultilevel(double const &overlap,
						   std::vector<std::vector<std::vector<double>>> const &centers,
						   std::unique_ptr<Kernel> motherKernel);

		~DynamicDiscardingMultilevel() {};

		double at(std::vector<double> const &x) const override;
		double atUpToLevel(std::vector<double> const &x, std::size_t const &givenLevel) const override;
		void solve(std::function<double(const std::vector<double> &)> const &targetFunc) override;
		void solve(std::vector<double> const &targetFunc) override;


	private:
		std::unique_ptr<Kernel> m_motherKernel;
		std::unique_ptr<std::vector<std::vector<std::vector<double>>>> m_centers;
		std::shared_ptr<std::vector<std::vector<std::vector<double>>>> m_centers_shared;
		std::vector<std::vector<std::vector<double>>> m_centers_vec;

		std::vector<std::unique_ptr<Kernel>> m_kernels;
		std::vector<std::unique_ptr<KDtree>> m_trees;

		//std::vector<Eigen::SparseMatrix<double>> m_kernelMatrices;
		std::vector<SparseMatrix> m_kernelMatrices;

		std::vector<std::vector<std::vector<double>>> m_points;

		//std::vector<Eigen::VectorXd> m_coefficients;
		std::vector<std::vector<double>> m_coefficients;

		std::vector<std::vector<double>> m_residuals;
		std::vector<std::vector<std::size_t>> m_indices;
		std::vector<std::vector<std::size_t>> m_actuallyUsedIndices;


		std::vector<double> m_radius;
		std::vector<double> m_threshold;
		std::vector<std::size_t> m_nPoints;

		std::size_t m_dimension;
		std::size_t m_nLevels;

		void setUpLevels();
		double atLevel(std::size_t level, std::vector<double> const &x) const;
		void updateResiduals(std::size_t const &level,
							 std::function<double(std::vector<double> &)> const &targetFunc);
};
