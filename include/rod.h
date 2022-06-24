#ifndef ROD_H
#define ROD_H

#include <eigen3/Eigen/Dense>
#include <boost/numeric/odeint.hpp>

#include <iostream>


/// \brief  Declare the state type which will be handled by the odeint integrator
typedef Eigen::Matrix<double, 7, 1 > state_type;



class Rod
{
public:
    Rod()
    {
        m_B = Eigen::MatrixXd(6, m_na);
        m_B << Eigen::Matrix3d::Identity(),
                Eigen::Matrix3d::Zero();
        m_B_bar = Eigen::MatrixXd(6, m_na);
        m_B_bar << Eigen::Matrix3d::Zero(),
                    Eigen::Matrix3d::Identity();

       


    }


    Eigen::MatrixXd GetRodShape(const Eigen::VectorXd &t_qe);

private:

    
    mutable double m_L { 1 };

    const unsigned int m_ne { 3 };
    const unsigned int m_na { 3 };

    Eigen::VectorXd m_strain_a0{ Eigen::Vector3d(0, 0, 0)};
    Eigen::VectorXd m_strain_c { Eigen::Vector3d(1, 0, 0)};


    Eigen::MatrixXd m_B;// {(Eigen::MatrixXd(6, m_na) << Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Zero()).finished()};
    Eigen::MatrixXd m_B_bar;// {(Eigen::MatrixXd(6, m_na) << Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Identity()).finished()};
    
    Eigen::VectorXd m_qe;
  
    //  The used integrator
    typedef boost::numeric::odeint::runge_kutta_dopri5< state_type, double,
                                                        state_type, double,
                                                        boost::numeric::odeint::vector_space_algebra> stepper;




    [[nodiscard]] Eigen::VectorXd getStrains(const double t_s)const;


    void forwardODEs(const state_type &t_y, state_type &t_dyds, const double t_s) const;




public:
    [[nodiscard]]const Eigen::MatrixXd Phi(const double t_s)const{


    const double X = t_s/m_L;

    // ---- Chebyshev ------
    Eigen::MatrixXd M_cheb(11,11);
    M_cheb  <<  1,    0,     0,      0,       0,        0,        0,        0,        0,        0,      0,
               -1,    2,     0,      0,       0,        0,        0,        0,        0,        0,      0,
                1,   -6,     6,      0,       0,        0,        0,        0,        0,        0,      0,
               -1,   12,   -30,     20,       0,        0,        0,        0,        0,        0,      0,
                1,  -20,    90,   -140,      70,        0,        0,        0,        0,        0,      0,
               -1,   30,  -210,    560,    -630,      252,        0,        0,        0,        0,      0,
                1,  -42,   420,  -1680,    3150,    -2772,      924,        0,        0,        0,      0,
               -1,   56,  -756,   4200,  -11550,    16632,   -12012,     3432,        0,        0,      0,
                1,  -72,  1260,  -9240,   34650,   -72072,    84084,   -51480,    12870,        0,      0,
               -1,   90, -1980,  18480,  -90090,   252252,  -420420,   411840,  -218790,    48620,      0,
                1, -110,  2970, -34320,  210210,  -756756,  1681680, -2333760,  1969110,  -923780, 184756;

    Eigen::VectorXd vec(11);
    vec << 1, X, pow(X, 2), pow(X, 3), pow(X, 4), pow(X, 5), pow(X, 6), pow(X, 7), pow(X, 8), pow(X, 9), pow(X, 10);

    Eigen::VectorXd res = (M_cheb*vec);
    Eigen::VectorXd Base = res.block(0,0,m_ne,1);



    Eigen::MatrixXd Phi = Eigen::MatrixXd::Zero(m_na, m_na*m_ne);

    for (unsigned int row=0;row<m_na;row++) {
        auto first_col = (row)*m_ne;
        //auto last_col = first_col + ne;
        Phi.block(row, first_col, 1, m_ne) << Base.transpose();
    }


    return Phi;
    }

    Eigen::MatrixXd Ke;


};

#endif // ROD_H
