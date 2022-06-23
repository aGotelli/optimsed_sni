/**
 * @file main.cpp
 * @author Andrea Gotelli (Andrea.Gotelli@ls2n.fr)
 * @brief This file contains benchmark for the spectral numerical integration
 * @version 0.1
 * @date 2022-06-23
 *
 * @copyright Copyright (c) 2022
 *
 */
#include <iostream>
#include <memory>


#include "integrator_base.h"
#include "quaternion_integrator.h"
#include "position_integrator.h"

#include <boost/math/special_functions/legendre.hpp>
#include <unsupported/Eigen/KroneckerProduct>

#include <benchmark/benchmark.h>

static constexpr unsigned int number_of_chebyshev_points = 16;

//  Define the Chebyshev points on the unit circle
static const auto x = ComputeChebyshevPoints(number_of_chebyshev_points);

static constexpr unsigned int ne = 3;
static constexpr unsigned int na = 3;

static constexpr unsigned int repetions = 5;



/*!
 * \brief Phi Compute the basis matrix Phi for a given X
 * \param t_X The coordinate in the rod normalized domain. Must be in [0, 1]
 * \param t_na The number of allowed strain coordinates
 * \param t_ne The number of modes per strain coordinate
 * \return The basis matrix Phi for a given X
 */
static const Eigen::MatrixXd Phi(const double t_X,
                                 const unsigned int t_na,
                                 const unsigned int t_ne,
                                 const double &t_begin=0,
                                 const double &t_end=1)
{

    //  The coordinate must be transposed into the Chebyshev domain [-1, 1];
    double x = ( 2 * t_X - ( t_end + t_begin) ) / ( t_end - t_begin );

    //  Compute the values of the polynomial for every element of the strain field
    Eigen::VectorXd Phi_i(t_ne, 1);
    for(unsigned int i=0; i<t_ne; i++)
        Phi_i[i] = boost::math::legendre_p(i, x);


    //  Define the matrix of bases
    Eigen::MatrixXd Phi = Eigen::KroneckerProduct(Eigen::MatrixXd::Identity(t_na, t_na), Phi_i.transpose());


    return Phi;
}




static void IntegrateQuaternion(benchmark::State& t_state)
{
    Eigen::Matrix<double, ne*na, 1> qe;
    //  Here we give some value for the strain
    qe <<   0,
            0,
            0,
            1.2877691307032,
           -1.63807499160786,
            0.437406679142598,
            0,
            0,
            0;


    std::vector<Eigen::Vector3d> K_stack(number_of_chebyshev_points);
    for(unsigned int i=0; i<number_of_chebyshev_points; i++)
        K_stack[i] = Phi(x[i], na, ne)*qe;


    std::shared_ptr<IntegratorBase> Q_integrator = std::make_shared<QuaternionIntegrator>( number_of_chebyshev_points, K_stack );


    Eigen::Vector4d Q_init(1, 0, 0, 0);


    while (t_state.KeepRunning())
        Q_integrator->Integrate(Q_init);


}
BENCHMARK(IntegrateQuaternion)->Repetitions(repetions);


static void IntegratePosition(benchmark::State& t_state)
{
    Eigen::Matrix<double, ne*na, 1> qe;
    //  Here we give some value for the strain
    qe <<   0,
            0,
            0,
            1.2877691307032,
           -1.63807499160786,
            0.437406679142598,
            0,
            0,
            0;


    std::vector<Eigen::Vector3d> K_stack(number_of_chebyshev_points);
    for(unsigned int i=0; i<number_of_chebyshev_points; i++)
        K_stack[i] = Phi(x[i], na, ne)*qe;

    std::vector<Eigen::Vector3d> Lambda_stack(number_of_chebyshev_points);
    for(unsigned int i=0; i<number_of_chebyshev_points; i++)
        Lambda_stack[i] = Eigen::Vector3d(1, 0, 0);


    std::shared_ptr<IntegratorBase> Q_integrator = std::make_shared<QuaternionIntegrator>( number_of_chebyshev_points, K_stack );
    std::shared_ptr<IntegratorBase> r_integrator = std::make_shared<Positionintegrator>( number_of_chebyshev_points, Q_integrator, Lambda_stack  );


    Eigen::Vector4d Q_init(1, 0, 0, 0);
    Q_integrator->Integrate(Q_init);


    Eigen::Vector3d r_init(0, 0, 0);


    while (t_state.KeepRunning())
        r_integrator->Integrate(r_init);



}
BENCHMARK(IntegratePosition)->Repetitions(repetions);




BENCHMARK_MAIN();
