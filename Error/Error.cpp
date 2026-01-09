/*
 * Error.cpp
 *
 *  Created on:  07.06.2019
 *  Last changed:30.08.2019
 *
 *      Author: kempf
 */



#include "Error.h"




Error::Error()
{
    m_dimension       =  CParam::getFirstValue<std::size_t>("dimension");

    
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
            
//     this->fillLooEvalGrid();
    this->fillL2EvalGrid();
}


Error::Error(std::string const &filename)
{
    m_dimension       =  CParam::getFirstValue<std::size_t>("dimension");

    
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
        
    std::fstream file;
    file.open(filename.c_str(), std::ios::in);             // filename.c_str() funktioniert für alle OS

    if(!file)
    {
        std::cout<< "Error in Error-class! Please insert file to read eval-grid data!"<< std::endl;
        exit(1);
    }

    for(std::size_t i = 0; i < 4 + m_dimension; i++)
    {
        file.ignore(256);
    }

    std::size_t i = 0;
//     for(std::size_t i = 0; i < ????????; i++)                   // <----------überlegen, ob auch einlesen, oder irgendwie mit while(! file.end) und push_back...
    while(!file.eof())
    {
        m_evalGrid.resize(m_dimension);
        for (std::size_t j = 0; j < m_dimension; j++)
        {
            m_evalGrid[i][j] = new double [1];
            
            file >> m_evalGrid[i][j][0];
        }
        i++;
    }
    file.close();
            
//     std::cout << "Domain: " << std::endl;
//     for(std::size_t i = 0; i < m_dimension; i++)
//     {
//         std::cout << m_domain[i].lo[0] << " to " << m_domain[i].hi[0] << std::endl;
//     }    CSparseGridParallel sparseGrid(m_domain);


//     // Build evaluation Grid as d-dim sparse grid 
//     CSparseGridParallel sparseGrid(m_domain);
//     std::vector<std::vector<double>> tmpVec = sparseGrid.getSparseGrid();
// 
//     m_evalGrid.resize(tmpVec.size());
// 
//     for(std::size_t i = 0; i < tmpVec.size(); i++)
//     {
//         m_evalGrid[i].resize(m_dimension);
//         for(std::size_t j = 0; j < m_dimension; j++)
//         {
//             m_evalGrid[i][j] = new double [1];
//             m_evalGrid[i][j][0] = tmpVec[i][j];
//         }
//     }
//     
//     sparseGrid.writeSparseGrid();
    
    std::cout << "Building evaluation sparse grid with " << m_evalGrid.size() << " points" << std::endl;     
    
}

Error::~Error()
{
//     std::cout << "In error destructor" << std::endl;
    for(std::size_t i = 0; i < m_evalGrid.size(); i++)
    {
        for(std::size_t j = 0; j < m_dimension; j++)
        {
            delete[] m_evalGrid[i][j];
            m_evalGrid[i][j] = nullptr;
        }
    }
}


double Error::getl2Error()
{
        return m_l2Error;
}



double Error::getlooError()
{
    return m_looError;
}



void Error::computel2Error(std::vector<HWfunction*> const &f, TensorInterpolant &interpolant)
{
    
    this->fillL2EvalGrid();

    m_l2Error = 0.0;    


    double tmpProduct;
    double tmpValue;
    
    for(std::size_t i = 0; i < m_evalGrid.size(); i++)
    {
        tmpProduct = 1.0;
        for(std::size_t j = 0; j < m_dimension; j++)
        {
//             std::cout << m_evalGrid[i][j][0] << " ";
            tmpProduct *= f[j]->at(m_evalGrid[i][j]);
        }
//         std::cout << endl;
        tmpValue = interpolant.at(m_evalGrid[i]);

        m_l2Error += (tmpProduct - tmpValue) * (tmpProduct - tmpValue);
    }
//     delete[] tmpX;
}


void Error::computel2Error_parallel(std::vector<HWfunction*> const &f, TensorInterpolant &interpolant)
{
    this->fillL2EvalGrid();
    
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "Evaluating L2-Error on" << m_evalGrid.size() << "points" << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
    
    m_l2Error = 0.0;
    
    std::size_t nThreads = thread::hardware_concurrency();
    std::vector<std::thread> threads;
    threads.resize(nThreads);
//     std::vector<double> tmpProduct;
//     tmpProduct.resize(nThreads);
//     std::vector<double> tmpValue;
//     tmpValue.resize(nThreads);
    
    std::vector<double> partialResults;
    partialResults.resize(nThreads);
    
    std::size_t nPointEvalPerThread = m_evalGrid.size() / nThreads;
    
    auto func = [&](std::size_t const &startIndex, std::size_t const &endIndex, std::size_t const &threadIndex)
    {
         double tmpProduct;
         double tmpValue;
    
         for(std::size_t i = startIndex; i < endIndex; i++)
         {
            tmpProduct = 1.0;
            for(std::size_t j = 0; j < m_dimension; j++)
            {
//             std::cout << m_evalGrid[i][j][0] << " ";
                tmpProduct *= f[j]->at(m_evalGrid[i][j]);
            }
//         std::cout << endl;
            tmpValue = interpolant.at(m_evalGrid[i]);

            partialResults[threadIndex] += (tmpProduct - tmpValue) * (tmpProduct - tmpValue);
        }
    };
    
    
    std::size_t startIdx = 0;
    std::size_t endIdx   = nPointEvalPerThread;
    // Threads work
    for(std::size_t i = 0; i < nThreads - 1; i++)
    {
//         startIdx = i * nPointEvalPerThread;
//         endIdx   = (i + 1) * nPointEvalPerThread - 1;
    
        threads[i] = std::thread(func, startIdx, endIdx, i);
    
        startIdx += nPointEvalPerThread;
        endIdx   += nPointEvalPerThread;
    }
    
    startIdx = (nThreads - 1) * nPointEvalPerThread;
    endIdx   = m_evalGrid.size();
    threads[nThreads - 1] = std::thread(func, startIdx, endIdx, nThreads - 1);

    for(std::size_t i = 0; i < nThreads; i++)
    {
        threads[i].join();
    }
    
    
    for(std::size_t i = 0; i < nThreads; i++)
    {
        m_l2Error += partialResults[i];
    }
}


void Error::computelooError_parallel(std::vector<HWfunction*> const &f, TensorInterpolant &interpolant)
{
    //this->fillLooEvalGrid();
    
    m_looError = 0.0;
    
    std::size_t nThreads = thread::hardware_concurrency();
    std::vector<std::thread> threads;
    threads.resize(nThreads);  
    
    std::vector<double> partialResults;
    partialResults.resize(nThreads);
    
    std::size_t nPointEvalPerThread = m_evalGrid.size() / nThreads;
    
    auto func = [&](std::size_t const &startIndex, std::size_t const &endIndex, std::size_t const &threadIndex)
    {
        partialResults[threadIndex] = 0.0;
        
        double tmp = 0.0;
        double tmpProduct;
        
        for(std::size_t i = startIndex; i < endIndex; i++)
        {
            tmpProduct = 1.0;
            for(std::size_t j = 0; j < m_dimension; j++)
            {
                tmpProduct *= f[j]->at(m_evalGrid[i][j]);
            }
            tmp = std::abs(tmpProduct - interpolant.at(m_evalGrid[i]));
            
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
    endIdx   = m_evalGrid.size();
    threads[nThreads - 1] = std::thread (func, startIdx, endIdx, nThreads - 1);
    
    for(std::size_t i = 0; i < nThreads; i++)
    {
        threads[i].join();
    }
    
    for(std::size_t i = 0; i < partialResults.size(); i++)
    {
        if(m_looError < partialResults[i])
        {
            m_looError = partialResults[i];
        }
    }
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "Evaluating Loo-Error on " << m_evalGrid.size() << " points" << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;

}



void Error::fillL2EvalGrid()
{
    //******************************************************************//
    // Needed for l2-error, thinning in CSparseGridParallel is only     //
    // sequential, since ANN only allows sequential search...           //
    //******************************************************************//
    // Build evaluation Grid as d-dim sparse grid 
    CSparseGridParallel sparseGrid(m_domain);
    sparseGrid.parallel_create_sparse_grid();
    std::vector<std::vector<double>> tmpVec = sparseGrid.getSparseGrid();

    m_evalGrid.resize(tmpVec.size());

    for(std::size_t i = 0; i < tmpVec.size(); i++)
    {
        m_evalGrid[i].resize(m_dimension);
        for(std::size_t j = 0; j < m_dimension; j++)
        {
            m_evalGrid[i][j] = new double [1];
            m_evalGrid[i][j][0] = tmpVec[i][j];
        }
    }
    std::cout << "Writing l2 Eval Grid" << std::endl;
    sparseGrid.writeSparseGrid();
}    
    
    
void Error::fillLooEvalGrid()
{
    CSparseGridParallel sparseGrid(m_domain);
    std::vector<std::vector<double>> tmpVec = sparseGrid.get_unthinned_sparse_grid();

    m_evalGrid.resize(tmpVec.size());

    for(std::size_t i = 0; i < tmpVec.size(); i++)
    {
        m_evalGrid[i].resize(m_dimension);
        for(std::size_t j = 0; j < m_dimension; j++)
        {
            m_evalGrid[i][j] = new double [1];
            m_evalGrid[i][j][0] = tmpVec[i][j];
        }
    }
}