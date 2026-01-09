/*
 * RandomLooError.cpp
 *
 *  Created on:  28.08.2019
 *  Last changed: 3.09.2019
 *
 *      Author: kempf
 */

#include "RandomLooError.h"


RandomLooError::RandomLooError()
{
    m_dimension       =  CParam::getFirstValue<std::size_t>("dimension");
    //m_nDirections       =  CParam::getFirstValue<std::size_t>("nDirections");
    m_nRandomPoints   =  CParam::getFirstValue<std::size_t>("nRandomPoints");

    
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
    
    this->fillGrid();
}



void RandomLooError::fillGrid()
{
    m_evalGrid.resize(m_nRandomPoints);
    
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

