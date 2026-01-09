/*
 *  ApproximationFactory.h
 *
 *  Created on: 24.02.2025
 *  Last changed: 23.10.2025
 *      Author: kempf
 *
 *      Class ApproximationFactory: Factory design pattern
 *
 *
 *
 */


#pragma once

#include <memory>
#include <functional>
//#include "smolyakMovingLeastSquares.hpp"
#include "Approximation.h"
//#include "interpolation_old.h"
#include "Multilevel.h"
#include "AnisotropicMultilevel.h"
#include "DynamicDiscardingMultilevel.h"
#include "NestedMultilevel.h"
#include "LagrangeMultilevel.h"
#include "types.h"


class ApproximationFactory
{
public:
		static std::unique_ptr<Approximation> create(MultilevelParams &params)
		{
			switch (params.typeApproximation) {
				case ApproximationType::NestedMultilevel:
					return std::make_unique<NestedMultilevel>(params);
				case ApproximationType::Multilevel:
					return std::make_unique<Multilevel>(params);
				case ApproximationType::AnisotropicMultilevel:
					return std::make_unique<AnisotropicMultilevel>(params);
				case ApproximationType::DynamicDiscardingMultilevel:
					return std::make_unique<DynamicDiscardingMultilevel>(params);
				case ApproximationType::LagrangeMultilevel:
					return std::make_unique<LagrangeMultilevel>(params);
				default:
					throw std::invalid_argument("Unknown ApproximationType");
			}
		}

    static std::unique_ptr<Approximation> create(ApproximationType type,
							double const &overlap,
							std::unique_ptr<NestedSitesGenerator> centers,
							std::unique_ptr<Kernel> motherKernel)
    {
        switch (type) {
            case ApproximationType::NestedMultilevel:
            	return std::make_unique<NestedMultilevel>(overlap,
														  std::move(centers),
														  std::move(motherKernel));
            case ApproximationType::Multilevel:
            	throw std::invalid_argument("ApproximationType::Multilevel not implemented for these "
            			"parameters");
            case ApproximationType::AnisotropicMultilevel:
            	throw std::invalid_argument("ApproximationType::AnisotropicMultilevel not implemented "
            			"for these parameters");
            case ApproximationType::DynamicDiscardingMultilevel:
                        	throw std::invalid_argument("ApproximationType::DynamicDiscardingMultilevel"
                        			" not implemented for these parameters");
            default:
                throw std::invalid_argument("Unknown ApproximationType");
        }
    }

    static std::unique_ptr<Approximation> create(ApproximationType type,
								double const &overlap,
								std::unique_ptr<std::vector<std::vector<std::vector<double>>>> centers,
								std::unique_ptr<Kernel> motherKernel)
    {
    	switch(type)
    	{
    		case ApproximationType::NestedMultilevel:
    			throw std::invalid_argument("ApproximationType::NestedMultilevel not implemented for these "
    			            			"parameters");
            case ApproximationType::Multilevel:
            	return std::make_unique<Multilevel>(overlap,
													std::move(centers),
													std::move(motherKernel));
            case ApproximationType::AnisotropicMultilevel:
            	return std::make_unique<AnisotropicMultilevel>(overlap,
															   std::move(centers),
															   std::move(motherKernel));
            case ApproximationType::DynamicDiscardingMultilevel:
                        	return std::make_unique<DynamicDiscardingMultilevel>(overlap,
            															   std::move(centers),
            															   std::move(motherKernel));
            default:
            	throw std::invalid_argument("Unknown ApproximationType");
    	}
    }


    static std::unique_ptr<Approximation> create(ApproximationType type,
								double const &overlap,
								std::shared_ptr<std::vector<std::vector<std::vector<double>>>> centers,
								std::unique_ptr<Kernel> motherKernel)
    {
    	switch(type)
    	{
    		case ApproximationType::NestedMultilevel:
    			throw std::invalid_argument("ApproximationType::NestedMultilevel not implemented for these "
    			            			"parameters");
            case ApproximationType::Multilevel:
            	return std::make_unique<Multilevel>(overlap,
													centers,
													std::move(motherKernel));
            case ApproximationType::AnisotropicMultilevel:
            	return std::make_unique<AnisotropicMultilevel>(overlap,
															   centers,
															   std::move(motherKernel));
            case ApproximationType::DynamicDiscardingMultilevel:
                        	return std::make_unique<DynamicDiscardingMultilevel>(overlap,
            															   centers,
            															   std::move(motherKernel));
            default:
            	throw std::invalid_argument("Unknown ApproximationType");
    	}
    }


    static std::unique_ptr<Approximation> create(ApproximationType type,
								double const &overlap,
								std::vector<std::vector<std::vector<double>>> const &centers,
								std::unique_ptr<Kernel> motherKernel)
    {
    	switch(type)
    	{
    		case ApproximationType::NestedMultilevel:
    			throw std::invalid_argument("ApproximationType::NestedMultilevel not implemented for these "
    			            			"parameters");
            case ApproximationType::Multilevel:
            	return std::make_unique<Multilevel>(overlap,
													centers,
													std::move(motherKernel));
            case ApproximationType::AnisotropicMultilevel:
            	return std::make_unique<AnisotropicMultilevel>(overlap,
															   centers,
															   std::move(motherKernel));
            case ApproximationType::DynamicDiscardingMultilevel:
                        	return std::make_unique<DynamicDiscardingMultilevel>(overlap,
            															   centers,
            															   std::move(motherKernel));
            default:
            	throw std::invalid_argument("Unknown ApproximationType");
    	}
    }


    static std::unique_ptr<Approximation> create(ApproximationType type,
								std::vector<double> const &overlap,
								std::unique_ptr<std::vector<std::vector<std::vector<double>>>> centers,
								std::unique_ptr<Kernel> motherKernel)
    {
    	switch(type)
    	{
    		case ApproximationType::NestedMultilevel:
    			throw std::invalid_argument("ApproximationType::NestedMultilevel not implemented for these "
    			            			"parameters");
            case ApproximationType::Multilevel:
            	throw std::invalid_argument("ApproximationType::Multilevel not implemented for these "
            	    			            			"parameters");
            case ApproximationType::AnisotropicMultilevel:
            	return std::make_unique<AnisotropicMultilevel>(overlap,
															   std::move(centers),
															   std::move(motherKernel));
            case ApproximationType::DynamicDiscardingMultilevel:
            	throw std::invalid_argument("ApproximationType::DynamicDiscardringMultilevel not implemented for these "
            	    			            			"parameters");
            default:
            	throw std::invalid_argument("Unknown ApproximationType");
    	}
    }


    static std::unique_ptr<Approximation> create(ApproximationType type,
								std::vector<std::pair<std::vector<std::vector<double>>,std::vector<double>>> const &scalingPairs,
								std::unique_ptr<std::vector<std::vector<std::vector<double>>>> centers,
								std::unique_ptr<Kernel> motherKernel)
    {
    	switch(type)
    	{
    		case ApproximationType::NestedMultilevel:
    			throw std::invalid_argument("ApproximationType::NestedMultilevel not implemented for these "
    			            			"parameters");
            case ApproximationType::Multilevel:
            	throw std::invalid_argument("ApproximationType::Multilevel not implemented for these "
            	    			            			"parameters");
            case ApproximationType::AnisotropicMultilevel:
            	return std::make_unique<AnisotropicMultilevel>(scalingPairs,
															   std::move(centers),
															   std::move(motherKernel));
            case ApproximationType::DynamicDiscardingMultilevel:
            	throw std::invalid_argument("ApproximationType::DynamicDiscardringMultilevel not implemented for these "
            	    			            			"parameters");
            default:
            	throw std::invalid_argument("Unknown ApproximationType");
    	}
    }
};

