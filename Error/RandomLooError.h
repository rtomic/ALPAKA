/*
 * RandomLooError.h
 *
 *  Created on:   4.09.2019
 *  Last changed: 4.09.2019
 *
 *      Author: kempf
 */


#pragma once


#include "CError.h"



class RandomLooError: public CError
{
    public:
        double computeError(std::vector<HWfunction*> const &f, TensorInterpolant &interpolant);
        
    private:
        std::size_t                          m_dimension;
        std::size_t                          m_nRandomPoints;
        std::vector<HWorthRect>              m_domain;
        std::vector<std::vector<HWpoint>>    m_evalGrid;
        
        double                               m_error       = 0.0;
        
        void fillGrid();
};
