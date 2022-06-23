/**
 * @file quaternion_integrator.h
 * @author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * @brief This file contains the implementation of the spectral integration for the quaternions
 * @version 0.1
 * @date 2022-06-23
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#ifndef QUATERNION_INTEGRATOR_H
#define QUATERNION_INTEGRATOR_H

#include "integrator_base.h"

class QuaternionIntegrator : public IntegratorBase {

public:
    QuaternionIntegrator(const unsigned int t_number_of_chebychev_points, const std::vector<Eigen::Vector3d>& t_K_stack);

    virtual void Integrate(const Eigen::VectorXd &t_Q_init) override;

    virtual void updateA() override;

    virtual void updateb(const Eigen::MatrixXd &t_parameters) override {}

    inline virtual Eigen::MatrixXd getStack()const override {return m_Q_stack;};

private:
    const unsigned int m_state_dimension { 4 };


    const std::vector<Eigen::Vector3d> m_K_stack;


    Eigen::MatrixXd m_A_at_chebychev_point { Eigen::MatrixXd::Zero(m_state_dimension, m_state_dimension) };



    const Eigen::MatrixXd m_D;
    const Eigen::MatrixXd m_D_NN;
    const Eigen::MatrixXd m_D_IN;

    Eigen::MatrixXd m_A_NN { m_D_NN };

    Eigen::VectorXd m_ivp { Eigen::VectorXd::Zero((m_number_of_chebychev_points-1)*m_state_dimension) };

    Eigen::VectorXd m_Q_stack { Eigen::VectorXd::Zero((m_number_of_chebychev_points-1)*m_state_dimension) };

};




#endif  //  QUATERNION_INTEGRATOR_H
