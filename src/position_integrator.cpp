/**
 * @file position_integrator.cpp
 * @author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * @brief This file contains the source code for the implementation of the spectral integration for the positions
 * @version 0.1
 * @date 2022-06-23
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "position_integrator.h"

Positionintegrator::Positionintegrator(const unsigned int t_number_of_chebychev_points,
                                       const std::shared_ptr<const IntegratorBase> t_Q_integrator,
                                       const std::vector<Eigen::Vector3d>& t_Lambda_stack) :
                                            IntegratorBase(t_number_of_chebychev_points),
                                            m_Q_integrator( t_Q_integrator ),
                                            m_Lambda_stack( t_Lambda_stack )
                                        {}

void Positionintegrator::updateb(const Eigen::MatrixXd &t_Q_stack)
{
    Eigen::Quaterniond q;

    for (unsigned int i = 0; i < m_number_of_chebychev_points-1; ++i) {


        q = { t_Q_stack(i),
              t_Q_stack(i  +  (m_number_of_chebychev_points-1)),
              t_Q_stack(i + 2*(m_number_of_chebychev_points-1)),
              t_Q_stack(i + 3*(m_number_of_chebychev_points-1)) };


        m_b_NN.block(i, 0, 1, 3) = (q.toRotationMatrix()*Eigen::Vector3d(1, 0, 0)).transpose();

    }
}





void Positionintegrator::Integrate(const Eigen::VectorXd &t_r_init)
{

    for(unsigned int i=0; i<m_number_of_chebychev_points-1; i++)
        m_ivp.row(i) = m_Dn_IN(i, 0) * t_r_init.transpose();


    updateb( m_Q_integrator->getStack() );

    m_r_stack = m_Dn_NN_inv*(m_b_NN - m_ivp);


}
