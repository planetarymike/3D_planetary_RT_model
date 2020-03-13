//constants.h -- physical constants

#ifndef __CONSTANTS_H
#define __CONSTANTS_H

#include <cmath>

//planetary constants
const double rMars = 3395e5;// cm, radius of Mars
const double mMars = .1076 * 5.98e27;// gm, Mass of Mars (from fraction of Earth mass)
const double aMars_typical = 1.41; // AU, typical Mars-Sun distance

//physical  constants
const double G = 6.67e-8;// dyn cm^2 gm^-2, Cavendish constant
const double kB = 1.38e-16;// erg K^-1, Boltzmann constant
const double clight = 3e10;// cm s^-1, speed of light
const double mH = 1.673e-24;// gm, mass of hydrogen atom
const double mCO2 = 44 * mH;// gm, mass of CO2 molecule

//mathematical constants
const double pi = M_PI;

//radiative transfer parameters
const double lyman_alpha_lambda = 121.6e-7;// cm
const double lyman_alpha_cross_section_total = 2.647e-2 * 0.416;// cm^2 Hz
// ^^^ total cross section of Ly alpha pi*e^2/(m c) * f
const double lyman_alpha_line_center_cross_section_coef = lyman_alpha_cross_section_total/std::sqrt(2.0*pi*kB/mH)*lyman_alpha_lambda; // cm^2
// ^^^line center cross section coefficient, needs to be divided by sqrt(temp) to yield line center cross section
const double CO2_lyman_alpha_absorption_cross_section = 6.3e-20; //cm^2 CO2 cross section at Ly alpha



const double lyman_alpha_flux_Earth_typical = 4.5e15;//photons/s/m2/Angstrom line center flux, Earth (Solar Min)
const double lyman_alpha_flux_Mars_typical = (lyman_alpha_flux_Earth_typical/1e4  // photons/s/cm2/A
					      *1e8*lyman_alpha_lambda*lyman_alpha_lambda/clight //photons/s/cm2/Hz
					      /aMars_typical/aMars_typical); // photons/s/cm2/Hz at Mars-Sun distance
const double lyman_alpha_typical_g_factor = lyman_alpha_flux_Mars_typical*lyman_alpha_cross_section_total; // s-1




const double lyman_beta_lambda = 102.6e-7;// cm
const double lyman_beta_cross_section_total = 2.647e-2 * 0.079142;// cm^2 Hz
// ^^^ total cross section of Ly beta pi*e^2/(m c) * f
const double lyman_beta_branching_ratio = 0.8819;
const double lyman_beta_line_center_cross_section_coef = lyman_beta_cross_section_total/std::sqrt(2.0*pi*kB/mH)*lyman_beta_lambda; // cm^2
// ^^^line center cross section coefficient, needs to be divided by sqrt(temp) to yield line center cross section
const double CO2_lyman_beta_absorption_cross_section = 3.53e-17; //cm^2 CO2 cross section at Ly beta


const double lyman_beta_flux_Earth_typical = lyman_alpha_flux_Earth_typical/66.; //photons/s/m2/Angstrom solar flux line center Earth
const double lyman_beta_flux_Mars_typical = (lyman_beta_flux_Earth_typical/1e4  // photons/s/cm2/A
					      *1e8*lyman_beta_lambda*lyman_beta_lambda/clight //photons/s/cm2/Hz
					      /aMars_typical/aMars_typical); // photons/s/cm2/Hz at Mars-Sun distance
const double lyman_beta_typical_g_factor = lyman_beta_flux_Mars_typical*lyman_beta_cross_section_total; // s-1





#endif
