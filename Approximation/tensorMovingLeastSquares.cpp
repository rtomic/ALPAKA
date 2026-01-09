/*
 *  tensorMovingLeastSquares.cpp
 *
 *  Created on: 03.06.2024
 *  Last changed: 20.08.2024
 *      Author: kempf
 *
 *      Implementation of class TensorMovingLeastSquares
 */

#include "tensorMovingLeastSquares.hpp"


//TPcenters = TPcenters[direction][Punkte][Koordinaten]

TensorMovingLeastSquares::TensorMovingLeastSquares(std::vector<std::vector<std::vector<double>>> const &TPcenters,
												   std::vector<std::size_t> const &kernelTypes,
												   std::vector<double> const &deltas,
												   std::vector<std::size_t> const &polDegrees,
												   std::vector<std::vector<double>> const &shifts,
												   std::vector<double> const &scalings)
{
    m_nDirections = TPcenters.size();
    m_dimensions.resize(m_nDirections);
    
    m_tpPoints = TPcenters;
    m_nPts.resize(m_nDirections);
    
    m_deltas = deltas;
    m_polDegrees = polDegrees;
    m_shifts = shifts;
    m_scalings = scalings;
    m_kernelTypes = kernelTypes;
    
    m_trees.resize(m_nDirections);
    m_indices.resize(m_nDirections);
    m_dimPolSpaces.resize(m_nDirections);
    //m_pointsHW.resize(m_nDirections);

    for(std::size_t i = 0; i < m_nDirections; i++)
    {
        m_dimensions[i] = TPcenters[i][0].size();
        m_nPts[i] = TPcenters[i].size();
        
        m_dimPolSpaces[i] = this->binom(m_dimensions[i] + m_polDegrees[i], m_dimensions[i]);
        m_indices[i] = new int[m_nPts[i]];
        for(std::size_t j = 0; j < m_nPts[i]; j++)
        {
            m_indices[i][j] = j;
        }
        //m_pointsHW[i] = allocPts(m_nPts[i], m_dimensions[i]);
        //stlToHw(m_tpPoints[i], m_pointsHW[i]);
        
        //m_trees[i] = new HWtree(m_pointsHW[i],m_indices[i], m_nPts[i], m_dimensions[i]);
        m_trees[i] = new HWtree(m_tpPoints[i],m_indices[i], m_nPts[i], m_dimensions[i]);

        //deallocPts(pointsHW);
    }

}


TensorMovingLeastSquares::~TensorMovingLeastSquares()
{
    for(std::size_t i = 0; i < m_nDirections; i++)
    {
        delete[] m_indices[i];
        delete m_trees[i];
        //deallocPts(m_pointsHW[i]);
    }
}



double TensorMovingLeastSquares::at(std::vector<std::vector<double>> const &x,
									std::function<double(std::vector<std::vector<double>> const &x)> func)
{

	std::vector<HWidxArray> localIndices(m_nDirections);
    std::vector<int> nLocals(m_nDirections);
    for(std::size_t i = 0; i < m_nDirections; i++)
    {
        //hwX[i] = new double[m_dimensions[i]];
        //stlToHw(x[i],hwX[i]);
        //nLocals[i] = m_trees[i]->rangeSearch(hwX[i], m_deltas[i], localIndices[i]);
        nLocals[i] = m_trees[i]->rangeSearch(x[i], m_deltas[i], localIndices[i]);

    }

    std::vector<MLSApproximation> directionWiseMLS;
    directionWiseMLS.reserve(m_nDirections);
    
    std::vector<std::vector<double>> shapeFuncValues(m_nDirections);

    for(std::size_t i = 0; i < m_nDirections; i++)
    {
        directionWiseMLS.emplace_back(m_tpPoints[i],m_kernelTypes[i],m_deltas[i],m_polDegrees[i],
									  m_shifts[i],m_scalings[i]);
        shapeFuncValues[i].reserve(nLocals[i]);
        shapeFuncValues[i] = directionWiseMLS[i].shapeFuncValues(x[i]);
/*
        std::cout << std::endl;
        std::cout << "shape Func Values: " << std::endl;
        for(std::size_t j = 0; j < shapeFuncValues[i].size(); j++)
        {
        	std::cout << shapeFuncValues[i][j] << std::endl;
        }
        std::cout << std::endl;
  */
    }

    std::vector<int> k(m_nDirections);
    for(std::size_t i = 0; i < m_nDirections; i++)
    {
        k[i] = 0;
    }
    
    double prodOfShapeValues;
    double sum = 0.0;
    std::vector<std::vector<double>> evalPoint(m_nDirections);
    while(true)
    {
        prodOfShapeValues = 1.0;

        for(std::size_t i = 0; i < m_nDirections; i++)
        {
            prodOfShapeValues *= shapeFuncValues[i][k[i]];
            evalPoint[i].reserve(m_dimensions[i]);

            evalPoint[i] = m_tpPoints[i][localIndices[i][k[i]]];
        }

        sum += func(evalPoint) * prodOfShapeValues;
        
        for(std::size_t j = 0; j < m_nDirections; j++)
        {
            k[j] += 1;
            if(k[j] > nLocals[j]-1)
            {
                if(j == m_nDirections-1)
                {
                    return sum;
                }
                k[j] = 0;
            }
            else
            {
                break;
            }
        }
    }
}



std::size_t TensorMovingLeastSquares::binom(std::size_t const &n, std::size_t const &k)
{
        return
          (        k> n  )? 0 :          // out of range
          (k==0 || k==n  )? 1 :          // edge
          (k==1 || k==n-1)? n :          // first
          binom(n - 1, k - 1) * n / k;   // recursive
}

/*
double TensorMovingLeastSquares::evalPolynomial(std::vector<double> const &evalPoint, std::vector<double> const &shift, double const &scaling, MultiIndex const &degree, std::size_t dimension)
{
    std::vector<double> diff (dimension);
    
    for(std::size_t i = 0; i < dimension; i++)
    {
        diff[i] = evalPoint[i] - shift[i];
        diff[i] /= scaling;
    }
    
    double result = 1.0;
    double d1;
    double d2;
    
    for(std::size_t i = 0; i < dimension; i++)
    {
        d1 = diff[i];
        d2 = d1 * d1;
        
        switch(degree.getMultiIndex()[i]-1)
        {
            case 0:
                result *= 1.0;
                break;
            case 1:
                result *= d1;
                break;
            case 2:
                result *= d2;
                break;
            case 3:
                result *= d1 * d2;
                break;
            case 4:
                result *= d2 * d2;
                break;
            case 5:
                result *= d2 * d2 * d1;
                break;
            case 6:
                result *= d2 * d2 * d2;
                break;
            case 7:
                result *= d2 * d2 * d2 * d1;
                break;
            default:
                result *= pow(d1,degree.getMultiIndex()[i]-1);
        }
    }
    
    return result;
}



Eigen::VectorXd TensorMovingLeastSquares::forwardSubstitution(Eigen::MatrixXd const &L, Eigen::VectorXd const &RHS)
{
    std::size_t dim = L.rows();
    Eigen::VectorXd result(dim);
    
    double tmp;
    for(std::size_t i = 0; i < dim; i++)
    {
        tmp = RHS(i);
        for(std::size_t j = 0; j < i; j++)
        {
            tmp -= L(i,j) * result(j);
        }
        result(i) = tmp / L(i,i);
    }
    
    return result;
}


Eigen::VectorXd TensorMovingLeastSquares::backwardSubstitution(Eigen::MatrixXd const &R, Eigen::VectorXd const &RHS)
{
    std::size_t dim = R.rows();
    Eigen::VectorXd result(dim);
    
    double tmp;
    for(int i = dim - 1; i >= 0; i--)
    {
        tmp = 0;
        for(std::size_t j = i+1; j < dim; j++)
        {
            tmp += R(i,j) * result(j);
        }
        result(i) = (RHS(i) - tmp) / R(i,i);
    }
    
    return result;
}


void TensorMovingLeastSquares::nestedForLoops(std::size_t const &threshold)
{
    std::vector<std::size_t> k(m_nDirections);
    for(std::size_t i = 0; i < m_nDirections; i++)
    {
        k[i] = 1;
    }
    while(true)
    {
        // DO SOMETHING
        for(std::size_t j = 1; j <= m_nDirections; j++)
        {
            k[j] += 1;
            if(k[j] > threshold)
            {
                if(j==m_nDirections)
                {
                    return;
                }
                k[j] = 1;
            }
            else
            {
                break;
            }
        }
    }
}
*/
