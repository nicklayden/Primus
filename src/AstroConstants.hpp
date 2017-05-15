#ifndef ASTROCONSTANTS_HPP
#define ASTROCONSTANTS_HPP

namespace constants{
  /// constants
  constexpr double mks_universal_G = 6.6726e-11; /// si units: Nm^2/kg^2
  constexpr double G = 6.6726e-11; /// si units: Nm^2/kg^2
  constexpr double au = 1.4960e11; /// meters
  constexpr double ly = 9.460730472580800e15; /// meters
  constexpr double parsec = 3.0857e16; /// meters
  constexpr double au_in_parsec = 206265; /// au per pc
  /// earth values
  constexpr double earth_tangent_vel = 29783.0; /// m/s
  constexpr double earth_mass = 5.972e24; /// kg
  constexpr double earth_radius = 6.371e6; /// meters
  constexpr double earth_eccentricity = 0.0167; // unitless
  constexpr double earth_perihelion = 0.9832687*au;
  constexpr double earth_aphelion = 1.01673*au;
  constexpr double earth_v_ap = 29291.35; // m/s
  constexpr double earth_density = 5514; // kg/m^3
  /// moon
  constexpr double moon_tangent_vel = 1023.0 + earth_tangent_vel; /// m/s
  constexpr double moon_mass = 7.347e22; /// kg
  constexpr double moon_init_pos = 384467000.0 + au; /// meters
  constexpr double moon_density = 3340; //kg/m^3
  /// sun values
  constexpr double solar_mass = 1.989e30; /// kg
  constexpr double teff_sun = 5840; /// Kelvin
  constexpr double solar_luminosity = 3.851e26; /// W
  constexpr double solar_flux_at_earth = 1369; /// W/m^2
  constexpr double solar_radius = 6.9598e8; /// meters
  constexpr double solar_density = 1410; // kg/m^3
  // Jupiter values
  constexpr double jupiter_mass = 1.8986e27; // kg
  constexpr double jupiter_aphelion = 5.45492*au; // meters
  constexpr double jupiter_semi_major = 778309574267; // meters
  constexpr double jupiter_v_ap = 12439.7; // m/s
  constexpr double jupiter_eccentricity = 0.048498; // unitless
  constexpr double jupiter_density = 1330; //kg/m^3
  //Mars
  constexpr double mars_mass= 6.4171e23; //kg
  constexpr double mars_aphelion = 1.6660*au;
  constexpr double mars_perihelion= 1.3814*au;
  constexpr double mars_avg_orb_vel = 24077; //m/s
  constexpr double mars_eccentricity = 0.0934;

  /// misc
  constexpr double pi = 3.14159265358979323846;
  constexpr double electron_rest_mass = 9.1094e-31; /// kg
  constexpr double electron_rest_energy_si = 8.1870e-14; /// J
  constexpr double electron_rest_energy_ev = 0.511e6; /// eV
  constexpr double proton_rest_mass = 1.6726e-27; /// kg
  constexpr double proton_rest_energy_si = 1.5033e-10; /// J
  constexpr double proton_rest_energy_ev = 0.93827e9; /// eV
  constexpr double amu = 1.6605e-27; /// kg
  ///more thermo/quantum constants
  constexpr double c = 2.9979e8; /// m/s
  constexpr double avogadro = 6.0221e23; /// 1/mol
  constexpr double boltzmann = 1.3807e-23; /// J/K
  constexpr double gas_constant = 8.3145; /// J/K/mol (R= k*Na)
  constexpr double planck_constant = 6.6261e-34; /// Js
  constexpr double stefan_boltzmann = 5.6705e-8; /// W/m^2K^4
  constexpr double fine_structure = 137.04; /// unitless
  constexpr double bohr_radius = 5.29818e-11; /// meters
  constexpr double thomson_cross = 6.6524e-29; ///m^2
  constexpr double hydrogen_ionization_energy_si = 2.1795e-18; /// J
  constexpr double hydrogen_ionization_energy_ev = 13.603; /// eV
  /// plancks quantum constants
  constexpr double planck_mass = 2.1767e-8; /// kg
  constexpr double planck_length = 1.6160e-35; /// meters
  constexpr double planck_time = 5.3907e-44; /// seconds

} /// namespace constants

#endif /* end of include guard: ASTROCONSTANTS_HPP */
