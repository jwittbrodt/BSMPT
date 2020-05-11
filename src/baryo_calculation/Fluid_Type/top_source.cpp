/*
 * top_source.cpp
 *
 *  Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller

		This program is free software: you can redistribute it and/or modify
		it under the terms of the GNU General Public License as published by
		the Free Software Foundation, either version 3 of the License, or
		(at your option) any later version.

		This program is distributed in the hope that it will be useful,
		but WITHOUT ANY WARRANTY; without even the implied warranty of
		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
		GNU General Public License for more details.

		You should have received a copy of the GNU General Public License
		along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <BSMPT/baryo_calculation/Fluid_Type/top_source.h>

/**
 * @file
 */

namespace BSMPT{
namespace Baryo{

typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;


double top_source::Calc_nL(double z_start, double z_end) const {
    bool debug = false;
    if(debug) std::cout<<"start of debug in "<<__func__<<std::endl;
    if(debug) std::cout<<"top class called"<<std::endl;

        /*
        omega[0] -> q 
        omega[1] -> t
        omega[2] -> h1 
        omega[3] -> h2  
        omega[4] -> q_prime
        omega[5] -> t_prime
        omega[6] -> h1_prime
        omega[7] -> h2_prime
        */
    state_type mu(8);
    mu = {0,0,0,0,0,0,0,0};
    if(debug) std::cout<<"Before ODE Calc:"<<std::endl;
    if(debug) for(size_t i=0;i<mu.size();i++) std::cout<<"\tmu["<<i<<"] = "<<mu[i]<<std::endl;
    const double C_AbsErr = 1e-9;
    const double C_RelErr = 1e-5;
    double stepsize_initial;
    if(z_start<z_end)  stepsize_initial= 1e-8;
    if(z_start>z_end)  stepsize_initial = -1e-8;
    double abs_err = C_AbsErr;
    double rel_err =C_RelErr;
    integrate_adaptive(make_controlled( abs_err , rel_err , error_stepper_type() ) , *this , mu , z_start , z_end , stepsize_initial );
    if(debug) std::cout<<"After ODE Calc:"<<std::endl;
    if(debug) for(size_t i=0;i<mu.size();i++) std::cout<<"\tmu["<<i<<"] = "<<mu[i]<<std::endl;

    return mu[0];

}



void top_source::operator()(const state_type &omega , state_type &domega , const double z)
{
    bool debug=false;
    if(debug)std::cout<<"Begin of debug in "<<__func__<<std::endl;

    /*
        Definition of all transport coefficients
    */
        std::vector<double> quark_mass;
        std::vector<double> quark_mass_prime;
        top_func(z , quark_mass, quark_mass_prime);
        double mt = quark_mass[0];
        double sym_phase = gen_fluid::symmetric_CP_violating_phase;
        double brk_phase = gen_fluid::broken_CP_violating_phase; 
        auto theta_vec =   Calc_theta(z , sym_phase , brk_phase);

        double theta    =   theta_vec[0];
        double theta_prime = theta_vec[1];
        //?????
        //Calculation of kappa_t
        Calc_kappa_obj.set_class(Temp, mt);
        double num_int  =   NIntegrate_kappa(Calc_kappa_obj);
        double kappa_q  =   kappa_QL_0*num_int;
        double kappa_tR =   kappa_QR_0*num_int;

        //Rescaled chemical potential like in 1811.11104
        double mu_M     =   (omega[1]/kappa_tR - omega[0]/kappa_q);
        double bR       =   -(omega[1]+omega[0]);
        double mu_SS    =   2*omega[0]/kappa_q - omega[1]/kappa_tR - bR/kappa_QR_0;
        double mu_Y     =   omega[1]/kappa_tR - omega[0]/kappa_q - omega[2]/kappa_H_0-omega[3]/kappa_H_0;

        Calc_Gam_obj.set_class(Temp,vw,mt, msqrt_thermal_top,dmsqrt_thermal_top);
        double Gam_M    =   Nintegrate_GamM(Calc_Gam_obj);
        Calc_Scp_obj.set_class(Temp,vw,mt,theta_prime,msqrt_thermal_top,dmsqrt_thermal_top);
        double Scp      =   Nintegrate_Scp(Calc_Scp_obj);
        if(debug)
        {
            std::cout<<"\tyuk_q = "<<yuk_q<<std::endl;
            std::cout<<"\tTemp = " <<Temp<<std::endl;
            std::cout<<"\tGam_Y = "<<Gam_Y<<std::endl;
            std::cout<<"\tGam_SS = "<<Gam_SS<<std::endl;
            std::cout<<"\tGam_M = "<<Gam_M<<std::endl;
            std::cout<<"\tScp = " <<Scp<<std::endl;
            std::cout<<"\ttheta = " << theta<<std::endl;
            std::cout<<"\ttheta_prime = " << theta_prime<<std::endl;
        }
    /*
        omega[0] -> q 
        omega[1] -> t
        omega[2] -> h1 
        omega[3] -> h2  
        omega[4] -> q_prime
        omega[5] -> t_prime
        omega[6] -> h1_prime
        omega[7] -> h2_prime
        */
        domega[0] = omega[4];
        domega[1] = omega[5];
        domega[2] = omega[6];
        domega[3] = omega[7];
        domega[4] = (Scp - Gam_M*mu_M + 2*Gam_SS*mu_SS - Gam_Y*mu_Y + vw* omega[4])/Dq;
        domega[5] = (-Scp + Gam_M* mu_M + Gam_Y*mu_Y + vw* omega[5])/Dt;
        domega[6] = (-Gam_Y*mu_Y + vw* omega[6])/Dh;
        domega[7] = (- Gam_Y*mu_Y + vw*omega[7])/Dh; 

}

}
}
