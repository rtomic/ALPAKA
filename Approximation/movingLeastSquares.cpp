/*
 * movingLeastSquares.cpp
 *
 *  Created on: 29.04.2024
 *  Last changed: 17.02.2025
 *      Author: kempf
 *
 *      Implementation of the class MLSApproximation
 */

#include "movingLeastSquares.hpp"



//MLSApproximation::MLSApproximation(std::vector<std::vector<double>> const &centers, std::size_t const &kernelType, double const &delta, HWfunction const *f, std::size_t const &polDegree, std::vector<double> const &shift, double const &scaling)
MLSApproximation::MLSApproximation(std::vector<std::vector<double>> const &centers,
								   std::size_t const &kernelType, double const &delta,
								   std::size_t const &polDegree,
								   std::vector<double> const &shift,
								   double const &scaling)

{
    m_dimension = CParam::getFirstValue<std::size_t>("dimension");
    m_points = centers;
    m_nPts = m_points.size();
    m_delta = delta;
    
    m_polDegree = polDegree;
    m_dimPolSpace = this->binom(m_dimension + m_polDegree, m_dimension);
    m_shift = shift;
    m_scaling = scaling;
    
    //m_pointsHW = allocPts(m_nPts, m_dimension);
    //stlToHw(m_points, m_pointsHW);

    m_kernelType = kernelType;
    m_kernel = KernelFactory::create(m_kernelType, m_dimension, m_delta);
    /*
    m_rightHandSide.resize(m_nPts);
    for(std::size_t i = 0; i < m_nPts; i++)
    {
        m_rightHandSide(i) = f->at(m_pointsHW[i]);
    }
   */
    m_index = new int[m_nPts];
    for(std::size_t i = 0; i < m_nPts; i++)
    {
        m_index[i] = i;
    }

    m_tree = new HWtree(m_points, m_index, m_nPts, m_dimension);
}


/*
MLSApproximation::MLSApproximation(std::vector<std::vector<double>> const &centers, std::size_t const &kernelType,
                                 double const &delta, std::vector<double> const f, std::size_t const &polDegree, std::vector<double> const &shift, double const &scaling)
{
    m_dimension = centers[0].size();
    m_points = centers;
    m_nPts = m_points.size();
    m_delta = delta;
    
    m_polDegree = polDegree;
    m_dimPolSpace = this->binom(m_dimension + m_polDegree, m_dimension);
    m_shift = shift;
    m_scaling = scaling;

    m_pointsHW = allocPts(m_nPts, m_dimension);
    stlToHw(m_points, m_pointsHW);

    m_kernelType = kernelType;
    m_kernel = KernelFactory::create(m_dimension,m_kernelType, m_delta);
  
    m_rightHandSide.resize(m_nPts);
    for(std::size_t i = 0; i < m_nPts; i++)
    {
        m_rightHandSide(i) = f[i];
    }
 
    m_index = new int[m_nPts];
    for(std::size_t i = 0; i < m_nPts; i++)
    {
        m_index[i] = i;
    }
    m_tree = new HWtree(m_pointsHW, m_index, m_nPts, m_dimension);
    
    //deallocPts(m_pointsHW);
}
*/
MLSApproximation::~MLSApproximation()
{
    delete[] m_index;
    delete m_tree;
}

double MLSApproximation::at(std::vector<double> const &x, HWfunction const *func)
{
	m_shift = x;

    HWpoint hwX = new double[m_dimension];
    stlToHw(x, hwX);
    
    HWidxArray indices;
    //int nLocal = m_tree->rangeSearch(hwX, m_delta, indices);
    int nLocal = m_tree->rangeSearch(x, m_delta, indices);
    std::cout << "nLocal: " << nLocal << std::endl;
    std::cout << "m_polDegree: " << m_polDegree << std::endl;
    std::cout << "dimPolSpace: " << m_dimPolSpace << std::endl;
    if(nLocal < m_dimPolSpace)
    {
        std::cout << "Too few points, select lower polynomial degree!" << std::endl;
        return -100000000;
    }
    
    Eigen::VectorXd shapeFunctionValues = this->computeShapeFunctionsQR(x, indices, nLocal);
    //Eigen::VectorXd shapeFunctionValues = this->computeShapeFunctions(x, indices, nLocal);

    double result = 0.0;
    for(int i = 0; i < nLocal; i++)
    {
        //result += m_rightHandSide(indices[i]) * shapeFunctionValues(i);
        result += func->at(m_pointsHW[indices[i]]) * shapeFunctionValues(i);
    }

    delete[] hwX;
    
    return result;
}

double MLSApproximation::at(std::vector<double> const &x, std::function<double(std::vector<double> const &x)> f)
{
	m_shift = x;
    HWidxArray indices;
    int nLocal = m_tree->rangeSearch(x, m_delta, indices);
    std::cout << "nLocal: " << nLocal << std::endl;
    std::cout << "m_polDegree: " << m_polDegree << std::endl;
    std::cout << "dimPolSpace: " << m_dimPolSpace << std::endl;
    if(nLocal < m_dimPolSpace)
    {
        std::cout << "Too few points, select lower polynomial degree!" << std::endl;
        return -100000000;
    }
    
    Eigen::VectorXd shapeFunctionValues = this->computeShapeFunctionsQR(x, indices, nLocal);
    //Eigen::VectorXd shapeFunctionValues = this->computeShapeFunctions(x, indices, nLocal);

    		double result = 0.0;
    for(int i = 0; i < nLocal; i++)
    {
        //result += m_rightHandSide(indices[i]) * shapeFunctionValues(i);
        result += f(m_points[indices[i]]) * shapeFunctionValues(i);
    }

    return result;
}

std::vector<double> MLSApproximation::shapeFuncValues(std::vector<double> const &x)
{
/*
	int nLocal = 0;
	HWidxArray indices;

	if(nLocal < m_polDegree)
	{
		std::cout << "here" << std::endl;

		//HWidxArray indices;
		nLocal = m_tree->rangeSearch(x, 1.2 * m_delta, indices);
	}
*/
	HWidxArray indices;
    int nLocal = m_tree->rangeSearch(x, m_delta, indices);
    std::vector<double> returnVec(nLocal);

    // For poldegree = 0, we have easier representation
    if(m_polDegree == 0)
    {
    	double tmp = 0.0;
		for(int k = 0; k < nLocal; k++)
		{
			tmp += m_kernel->at(x,m_points[indices[k]]);
		}

    	for(int i = 0; i < nLocal; i++)
    	{

    		returnVec[i] = m_kernel->at(x,m_points[indices[i]]) / tmp;
    	}
    }
    // For poldegree > 0, we have to work
    else
    {
    	for(int i = 0; i < nLocal; i++)
		{
    		//returnVec[i] = this->computeShapeFunctions(x, indices, nLocal)(i);
		 	returnVec[i] = this->computeShapeFunctionsQR(x, indices, nLocal)(i);
		}
    }

    return returnVec;
}



std::vector<double> MLSApproximation::at(std::vector<std::vector<double>> const &multipleX, HWfunction const *func)
{
    std::size_t nEvalPoints = multipleX.size();
    
    std::vector<double> returnValues(nEvalPoints);
    
    for(std::size_t i = 0; i < nEvalPoints; i++)
    {
        returnValues[i] = this->at(multipleX[i],func);
    }
    
    return returnValues;
}


Eigen::VectorXd MLSApproximation::computeShapeFunctions(std::vector<double> const &x,
														HWidxArray const &indices,
                                                        int const &nLocal)
{
    std::vector<double> weights(m_dimension);
    for(std::size_t i = 0; i < m_dimension; i++)
    {
        weights[i] = 1.0;
    }
    RKIndexSet setOfPolynomialDegrees(weights, m_dimension, m_polDegree);
    std::list<MultiIndex> polynomialDegreeList = setOfPolynomialDegrees.getIndexSet();
        
    Eigen::MatrixXd P(nLocal, m_dimPolSpace);
    
    std::size_t counter;
    Eigen::VectorXd RHS(m_dimPolSpace);
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(nLocal,nLocal);
    for(int i = 0; i < nLocal; i++)
    {
        D(i,i) = m_kernel->at(x,m_points[indices[i]]);
        counter = 0;
        for(std::list<MultiIndex>::iterator iter = polynomialDegreeList.begin();
        		iter != polynomialDegreeList.end();
        		++iter)
        {
            P(i,counter) = this->evalPolynomial(m_points[indices[i]], m_shift, m_scaling, *iter);
            RHS(counter) =this->evalPolynomial(x, m_shift, m_scaling, *iter);
            counter++;
        }
    }
        
    Eigen::MatrixXd PTDP(m_dimPolSpace,m_dimPolSpace);
    PTDP = P.transpose() * D * P;
    
    Eigen::MatrixXd solution = PTDP.ldlt().solve(RHS);

    return D*P*solution;
}



Eigen::VectorXd MLSApproximation::computeShapeFunctionsQR(std::vector<double> const &x, HWidxArray const &indices,
                                                        int const &nLocal)
{
	m_shift = x;

    std::vector<double> weights(m_dimension);
    for(std::size_t i = 0; i < m_dimension; i++)
    {
        weights[i] = 1.0;
    }
    RKIndexSet setOfPolynomialDegrees(weights, m_dimension, m_polDegree);
    std::list<MultiIndex> polynomialDegreeList = setOfPolynomialDegrees.getIndexSet();
    
    Eigen::MatrixXd sqrtDP(nLocal, m_dimPolSpace);
    Eigen::MatrixXd DP(nLocal, m_dimPolSpace);
    Eigen::VectorXd RHS(m_dimPolSpace);

    std::size_t counter;
    double sqrtDi;
    double Di;
    double tmp;
    for(int i = 0; i < nLocal; i++)
    {
        counter = 0;
        Di = m_kernel->at(x,m_points[indices[i]]);
        sqrtDi = sqrt(Di);
        for(std::list<MultiIndex>::iterator iter = polynomialDegreeList.begin(); iter != polynomialDegreeList.end(); ++iter)
        {
            tmp = this->evalPolynomial(m_points[indices[i]], m_shift, m_scaling, *iter);
            sqrtDP(i,counter) = sqrtDi * tmp;
            DP(i,counter) = Di * tmp;
            if(i == 0)
            {
                RHS(counter) = this->evalPolynomial(x, m_shift, m_scaling, *iter);
            }
            
            counter++;
        }
    }
    
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(sqrtDP);
    qr.compute(sqrtDP);
    Eigen::MatrixXd temp = qr.matrixQR().triangularView<Eigen::Upper>();
    
    Eigen::MatrixXd R(m_dimPolSpace,m_dimPolSpace);
    for(int i = 0; i < m_dimPolSpace; i++)
    {
        for(int j = 0; j < m_dimPolSpace; j++)
        {
            R(i,j) = temp(i,j);
        }
    }
    
    Eigen::VectorXd b = this->forwardSubstitution(R.transpose(), RHS);
    Eigen::VectorXd a = this->backwardSubstitution(R, b);
    
    Eigen::VectorXd returnVec = DP * a;
    return returnVec;
}


std::size_t MLSApproximation::binom(std::size_t const &n, std::size_t const &k)
{
        return
          (        k> n  )? 0 :          // out of range
          (k==0 || k==n  )? 1 :          // edge
          (k==1 || k==n-1)? n :          // first
          binom(n - 1, k - 1) * n / k;   // recursive
}


double MLSApproximation::evalPolynomial(std::vector<double> const &evalPoint, std::vector<double> const &shift, double const &scaling, MultiIndex const &degree)
{
    std::vector<double> diff (m_dimension);
    
    for(std::size_t i = 0; i < m_dimension; i++)
    {
        diff[i] = evalPoint[i] - shift[i];
        diff[i] /= scaling;
    }
    
    double result = 1.0;
    double d1;
    double d2;
    
    for(std::size_t i = 0; i < m_dimension; i++)
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


Eigen::VectorXd MLSApproximation::forwardSubstitution(Eigen::MatrixXd const &L, Eigen::VectorXd const &RHS)
{
    std::size_t dim = L.rows();
    Eigen::VectorXd result(dim);
    
    double tmp;
    for(std::size_t i = 0; i < dim; i++)
    {
        //tmp = RHS(i);
        tmp = 0.0;
        for(std::size_t j = 0; j < i; j++)
        {
            //tmp -= L(i,j) * result(j);
            tmp = tmp + L(i,j) * result(j);
        }
        //result(i) = tmp / L(i,i);
        result(i) = (RHS(i)-tmp) / L(i,i);
    }
    
    return result;
}


Eigen::VectorXd MLSApproximation::backwardSubstitution(Eigen::MatrixXd const &R, Eigen::VectorXd const &RHS)
{
    std::size_t dim = R.rows();
    Eigen::VectorXd result(dim);
    
    result(dim-1) = RHS(dim-1) / R(dim-1,dim-1);
    for(int i = dim-2; i >= 0; i--)
    {
        result(i)= RHS(i);
        for(std::size_t j = i+1; j <= dim-1; j++)
        {
            result(i) = result(i) - R(i,j) * result(j);
        }
        result(i) = result(i) / R(i,i);
    }

    return result;
}
