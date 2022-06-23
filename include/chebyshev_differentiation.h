/**
 * @file chebyshev_differentiation.h
 * @author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * @brief This file contains the functions related to Chebyshev points and differentiation matrix
 * @version 0.1
 * @date 2022-06-23
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef CHEBYSHEV_DIFFERENTIATION_H
#define CHEBYSHEV_DIFFERENTIATION_H

#include <eigen3/Eigen/Dense>
#include <iostream>

/*!
 * \brief ComputeChebyshevPoints Computes the Chebyshev points in the given interval
 * \tparam t_N The number of Chebyshev points.
 * \tparam t_L The length of the interval. Default 1 for the interval [0, 1]
 * \return An std::array containing the Chebyshev points
 */
static std::vector<double> ComputeChebyshevPoints(const unsigned int t_number_of_chebyshev_nodes,
                                                                                const unsigned int t_L=1)
{
    std::vector<double> x(t_number_of_chebyshev_nodes);

    unsigned int j = 0;
    std::generate(x.begin(), x.end(), [&](){
        return (static_cast<double>(t_L)/2)*(1 +cos( M_PI * static_cast<double>(j++) / static_cast<double>(t_number_of_chebyshev_nodes-1) ));
    });

    return x;
}

/*!
 * \brief ComputeChebyshevPoints Computes the c coefficients used in the definition of the Chebyshev differentiation matrix
 * \tparam t_N The number of Chebyshev points.
 * \return An std::array containing the coeffieints
 */
static std::vector<double> GetCoefficients_c(const unsigned int t_number_of_chebyshev_nodes)
{
    std::vector<double> c(t_number_of_chebyshev_nodes);

    unsigned int i = 0;
    std::generate(c.begin(), c.end(), [&](){
        //  gain is 2 in the edges and 1 elsewhere
        const unsigned int gain = (i==0 or i==t_number_of_chebyshev_nodes-1) ? 2 : 1;

        //  Follows the formula
        return pow(-1, i++)*gain;
    });

    return c;
}

/*!
 * \brief getDn Computes the Chebyshev differentiation matrix
 * \tparam t_N The number of Chebyshev points.
 * \return The Chebyshev differentiation matrix
 */
static Eigen::MatrixXd getDn(const unsigned int t_number_of_chebyshev_nodes)
{

    //  Define the Chebyshev points on the unit circle
    const auto x = ComputeChebyshevPoints(t_number_of_chebyshev_nodes);


    //  Create a matrix every row filled with a point value
    Eigen::MatrixXd X(t_number_of_chebyshev_nodes, t_number_of_chebyshev_nodes);
    for(unsigned int i=0; i<X.rows(); i++)
        X(i, Eigen::all) = Eigen::RowVectorXd::Constant(1, X.cols(), x[i]);


    //  Now compute the array containing the coefficients used in the definition of Dn
    const auto c = GetCoefficients_c(t_number_of_chebyshev_nodes);


    //  Create the appropriate matrix of coefficients
    Eigen::MatrixXd C(t_number_of_chebyshev_nodes, t_number_of_chebyshev_nodes);
    for(unsigned int i=0; i<t_number_of_chebyshev_nodes;i++) {
        for(unsigned int j=0; j<t_number_of_chebyshev_nodes;j++) {
            C(i,j) = c[i]/c[j];
        }
    }

    //  Definition of the temporary matrix Y
    const auto dX = X - X.transpose() + Eigen::MatrixXd::Identity(t_number_of_chebyshev_nodes, t_number_of_chebyshev_nodes);

    //  Declare the differentiation matrix
    Eigen::MatrixXd  Dn(t_number_of_chebyshev_nodes, t_number_of_chebyshev_nodes);


    //  Obtain off diagonal element for the differentiation matrix
    for(unsigned int i=0; i<t_number_of_chebyshev_nodes;i++) {
        for(unsigned int j=0; j<t_number_of_chebyshev_nodes;j++) {
            Dn(i,j) = C(i, j) / dX(i, j);
        }
    }


    //  Remove row sum from the diagonal of Dn
    Dn.diagonal() -= Dn.rowwise().sum();

    //  Finally return the matrix
    return Dn;
}



#endif // CHEBYSHEV_DIFFERENTIATION_H
