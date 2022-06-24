#include "rod.h"

Eigen::Matrix4d getA(const Eigen::Vector3d t_k){
    Eigen::Matrix4d A;
    A   <<     0  , -t_k(0),  -t_k(1),  -t_k(2),
            t_k(0),     0  ,   t_k(2),  -t_k(1),
            t_k(1), -t_k(2),      0  ,   t_k(0),
            t_k(2),  t_k(1),  -t_k(0),      0  ;

    return A;
}

Eigen::Matrix3d skew(const Eigen::Vector3d &t_v) {

    Eigen::Matrix3d v_hat;
    v_hat <<  0   ,  -t_v(2),   t_v(1),
            t_v(2),     0   ,  -t_v(0),
           -t_v(1),   t_v(0),     0   ;

    return v_hat;
}



Eigen::MatrixXd Rod::GetRodShape(const Eigen::VectorXd &t_qe)
{
    //  Copy the strains generalized coordinates (and derivatives)
    m_qe = t_qe;

    //  Now the intial velocities ans accelerations are nu

    /*  The forward state has the form
     *  | Q |   w, x, y, z                  0-3
     *  | r |   x, y, z                     4-6
     */

    state_type y;
    y << 1, 0, 0, 0,
         Eigen::Vector3d::Zero();

    const double step = 0.005;
    const unsigned int steps = m_L/step;

    Eigen::MatrixXd stack_of_points(steps+1, 3);

    unsigned int index = 0;

//  Forward Integration with observer
    boost::numeric::odeint::integrate_const(stepper(),
                                            [this](const state_type& t_y, state_type& t_dyds, const double t_s){
                                                this->forwardODEs(t_y, t_dyds, t_s);}, y, 0.0, m_L, step,
    [&](const state_type& t_y, const double){
        stack_of_points(index, Eigen::all) = t_y.block<3,1>(4,0).transpose();
        index++;
    });


    return stack_of_points;
}






[[nodiscard]] Eigen::VectorXd Rod::getStrains(const double t_s)const
{
    const Eigen::VectorXd strain_a = m_strain_a0 + Phi(t_s)*m_qe;
    const Eigen::VectorXd strain = m_B*strain_a + m_B_bar*m_strain_c;

    return strain;
}

void Rod::forwardODEs(const state_type &t_y, state_type &t_dyds, const double t_s) const
{
/*  Preprocessing    */


//    Get the strain from the rod
    const Eigen::VectorXd strain = getStrains(t_s);

//  Decompose the strain
    const Eigen::Vector3d k = strain.block<3,1>(0,0);
    const Eigen::Vector3d gamma = strain.block<3,1>(3,0);


//  The state has the form
//  | Q |   w, x, y, z
//  | r |   x, y, z

    const Eigen::Quaterniond Q(t_y[0], t_y[1],t_y[2], t_y[3]);
    [[maybe_unused]] const Eigen::Vector3d r = t_y.block<3,1>(4,0);

    const Eigen::Matrix3d R = Q.toRotationMatrix();


//  Actual ODE
    //      ATTENTION, A(K) SHOULD BE CHQNGED TO BE DIRECTLY USED WITH Eigen
    //const Eigen::Vector4d Q_prime = 0.5*getA(k)*Q.coeffs();
    const Eigen::Vector4d Q_prime = 0.5*getA(k)*Eigen::Vector4d(Q.w(),Q.x(),Q.y(),Q.z());
    const Eigen::Vector3d r_prime = R*gamma;

//  Postprocessing
    t_dyds << Q_prime,
              r_prime;

}


