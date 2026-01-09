/*
 * ErrorFactory.h
 *
 *  Created on:  28.08.2019
 *  Last changed: 02.07.2025
 *
 *      Author: kempf
 */


#pragma once


#include <memory>

//#include "CError.h"
#include "LooError.h"
#include "RandomPointsError.h"


enum class ErrorType
{
		MaximumError,
		RandomPointsError
};


class ErrorFactory
{
    public:
        static std::unique_ptr<CError> create(ErrorType type, std::vector<std::vector<double>> errorSites)
        {
        	switch(type)
        	{
        		case ErrorType::MaximumError:
        			return std::make_unique<LooError>(errorSites);
        		case ErrorType::RandomPointsError:
        			return std::make_unique<RandomPointsError>(errorSites);

        		default:
        			throw std::invalid_argument("Unknown ErrorType");
        	};
        };

        static std::unique_ptr<CError> create(ErrorType type,
					std::shared_ptr<std::vector<std::vector<std::vector<double>>>> errorSites)
                {
                	switch(type)
                	{
                		case ErrorType::MaximumError:
                			return std::make_unique<LooError>(errorSites);
                		case ErrorType::RandomPointsError:
                			return std::make_unique<RandomPointsError>(errorSites);

                		default:
                			throw std::invalid_argument("Unknown ErrorType");
                	};
                };
};
