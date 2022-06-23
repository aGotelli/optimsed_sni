/**
 * @file quaternion_integrator.cpp
 * @author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * @brief This file contains the implementation of the spectral integration for the quaternions
 * @version 0.1
 * @date 2022-06-23
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "quaternion_integrator.h"

#include <eigen3/unsupported/Eigen/KroneckerProduct>


QuaternionIntegrator::QuaternionIntegrator(const unsigned int t_number_of_chebychev_points, const std::vector<Eigen::Vector3d>& t_K_stack) :
                            IntegratorBase(t_number_of_chebychev_points),
                            m_K_stack( t_K_stack ),
                            m_D( Eigen::KroneckerProduct(Eigen::MatrixXd::Identity(m_state_dimension, m_state_dimension), m_Dn) ),
                            m_D_NN( Eigen::KroneckerProduct(Eigen::MatrixXd::Identity(m_state_dimension, m_state_dimension), m_Dn_NN) ),
                            m_D_IN( Eigen::KroneckerProduct(Eigen::MatrixXd::Identity(m_state_dimension, m_state_dimension), m_Dn_IN) )
                        {}


void QuaternionIntegrator::updateA()
{

    for(unsigned int i=0; i<m_number_of_chebychev_points-1; i++){



        //  Compute the A matrix of Q' = 1/2 A(K) Q
        m_A_at_chebychev_point  <<          0       ,   -m_K_stack[i](0),   -m_K_stack[i](1),   -m_K_stack[i](2),
                                    m_K_stack[i](0) ,           0       ,    m_K_stack[i](2),   -m_K_stack[i](1),
                                    m_K_stack[i](1) ,   -m_K_stack[i](2),           0       ,    m_K_stack[i](0),
                                    m_K_stack[i](2) ,    m_K_stack[i](1),   -m_K_stack[i](0),           0       ;


        for (unsigned int row = 0; row < m_state_dimension; ++row) {
            for (unsigned int col = 0; col < m_state_dimension; ++col) {
                int row_index = row*(m_number_of_chebychev_points-1)+i;
                int col_index = col*(m_number_of_chebychev_points-1)+i;
                m_A_NN(row_index, col_index) = m_D_NN(row_index, col_index) - 0.5*m_A_at_chebychev_point(row, col);
            }
        }

    }

}


void QuaternionIntegrator::Integrate(const Eigen::VectorXd &t_Q_init)
{

    updateA();

    m_ivp = m_D_IN*t_Q_init;

    m_Q_stack = - m_A_NN.inverse() * m_ivp;

}
