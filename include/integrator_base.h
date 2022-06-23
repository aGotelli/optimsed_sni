/**
 * @file integrator_base.h
 * @author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * @brief This files contains the base class for the spectral numerical integration
 * @version 0.1
 * @date 2022-06-23
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef INTEGRATOR_BASE_H   
#define INTEGRATOR_BASE_H  


#include <eigen3/Eigen/Dense>

#include "chebyshev_differentiation.h"

struct IntegratorBase {

    IntegratorBase(const unsigned int t_number_of_chebychev_points) :   m_number_of_chebychev_points( t_number_of_chebychev_points ),
                                                                        m_chebyshev_points( ComputeChebyshevPoints(t_number_of_chebychev_points) ),
                                                                        m_Dn( getDn(t_number_of_chebychev_points) ),
                                                                        m_Dn_NN( m_Dn.block(0, 0, t_number_of_chebychev_points-1, t_number_of_chebychev_points-1) ),
                                                                        m_Dn_IN( m_Dn.block(0, t_number_of_chebychev_points-1, t_number_of_chebychev_points-1, 1) ){}

    virtual ~IntegratorBase()=default;


    virtual void updateA()=0;
    virtual void updateb(const Eigen::MatrixXd &t_parameters)=0;

    virtual void Integrate(const Eigen::VectorXd &t_init)=0;

    inline virtual Eigen::MatrixXd getStack()const=0;

    const unsigned int m_number_of_chebychev_points;
    const std::vector<double> m_chebyshev_points;

    const Eigen::MatrixXd m_Dn;
    const Eigen::MatrixXd m_Dn_NN;
    const Eigen::MatrixXd m_Dn_IN;

};



#endif  //  INTEGRATOR_BASE_H
