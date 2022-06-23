/**
 * @file position_integrator.h
 * @author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * @brief This file contains the implementation of the spectral integration for the positions
 * @version 0.1
 * @date 2022-06-23
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#ifndef POSITIONINTEGRATOR_H
#define POSITIONINTEGRATOR_H
#include <memory>

#include "integrator_base.h"
#include "quaternion_integrator.h"

class Positionintegrator : public IntegratorBase
{
public:
    Positionintegrator(const unsigned int t_number_of_chebychev_points, const std::shared_ptr<const IntegratorBase> t_Q_integrator,  const std::vector<Eigen::Vector3d>& t_Lambda_stack);

    virtual void Integrate(const Eigen::VectorXd &t_r_init) override;

    virtual void updateA() override {}
    virtual void updateb(const Eigen::MatrixXd &t_Q_stack) override;

    inline virtual Eigen::MatrixXd getStack()const override {return m_r_stack;};


private:

    const unsigned int m_state_dimension { 3 };


    const std::vector<Eigen::Vector3d> m_Lambda_stack;

    Eigen::Quaterniond m_q_at_chebushev_point;

    const Eigen::MatrixXd m_Dn_NN_inv { m_Dn_NN.inverse() };

    Eigen::MatrixXd m_ivp { Eigen::MatrixXd::Zero(m_number_of_chebychev_points-1, m_state_dimension) };


    Eigen::MatrixXd m_b_NN {  Eigen::MatrixXd::Zero(m_number_of_chebychev_points-1, 3) };
    Eigen::MatrixXd m_r_stack { Eigen::MatrixXd(m_number_of_chebychev_points-1, m_state_dimension) };

    const std::shared_ptr<const IntegratorBase> m_Q_integrator;

};

#endif // POSITIONINTEGRATOR_H
