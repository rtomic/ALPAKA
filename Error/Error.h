/*
 * Error.h
 *
 *  Created on:  07.06.2019
 *  Last changed: 9.07.2019
 *
 *      Author: kempf
 */

#pragma once 

//General includes
#include <vector>
#include <list>
#include <iostream>
#include <cmath>
#include <thread>
#include <iterator>

//Tensor Problem includes
#include "IndexSet.hpp"
//#include "TensorInterpolant.h"
//#include "CSparseGridParallel.h"

//Multilevel includes
//#include "multilevel/hwmultilevel.h"
//#include "multilevel/hwnestedmultilevel.h"
//#include "general/hwfunction.h"
//#include "general/hwcenters.h"
//#include "general/hwutility.h"
//#include "general/hwbasics.h"
//#include "general/CParam.h"




class Error
{
    public:
        Error();
        Error(std::string const &filename);
        ~Error();//{std::cout << "In error destructor" << std::endl;};
        
        double getl2Error();
        double getlooError();
        void computel2Error(std::vector<HWfunction*> const &f, TensorInterpolant &interpolant);
        void computel2Error_parallel(std::vector<HWfunction*> const &f, TensorInterpolant &interpolant);
        
        void computelooError_parallel(std::vector<HWfunction*> const &f, TensorInterpolant &interpolant);
        

            
    private:
        std::size_t                          m_dimension;
        std::vector<HWorthRect>              m_domain;
        std::vector<std::vector<HWpoint>>    m_evalGrid;
        
        double                               m_l2Error       = 0.0;
        double                               m_looError      = 0.0;
        
        
        void fillL2EvalGrid();
        void fillLooEvalGrid();

};
