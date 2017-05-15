#ifndef PLANETS_HPP
#define PLANETS_HPP

/*
  Solar system planet properties
  including the sun, mercury, venus, earth, mars, jupiter, saturn, uranus,
  neptune
*/
namespace Planets {
namespace sun {
  // Assuming sun's position is (0,0,0), and velocity is (0,0,0)
  constexpr double mass = 1.989e30; /// kg
  constexpr double teff = 5840; /// Kelvin
  constexpr double radius = 6.9598e8; /// meters
  constexpr double density = 1410; // kg/m^3
  constexpr double luminosity = 3.851e26; /// W
  constexpr double flux_at_earth = 1369; /// W/m^2
} // sun

namespace earth {
  constexpr double orbitalvelocity = 29783.0; /// m/s
  constexpr double mass = 5.972e24; /// kg
  constexpr double radius = 6.371e6; /// au
  constexpr double eccentricity = 0.0167; // unitless
  constexpr double perihelion = 0.9832687; // au
  constexpr double aphelion = 1.01673; // au
  constexpr double v_aphelion = 29291.35; // m/s
  constexpr double density = 5514; // kg/m^3
} // earth

namespace jupiter {
  constexpr double mass = 1.8986e27; // kg
  constexpr double aphelion = 5.45492; // au
  constexpr double semi_major = 778309574267; // meters
  constexpr double v_aphelion = 12439.7; // m/s
  constexpr double eccentricity = 0.048498; // unitless
  constexpr double density = 1330; //kg/m^3
} // jupiter

namespace mars {
  constexpr double mass = 6.4171e23; // kg
  constexpr double aphelion = 1.6660; // au
  constexpr double perihelion = 1.3814; // au
  constexpr double orbitalvelocity = 24077; // m/s
  constexpr double eccentricity = 0.0934;
  constexpr double v_aphelion = 21972; // m/s
} // mars

namespace saturn {

} // saturn

namespace uranus {

} // uranus

namespace neptune {

} // neptune

} // planets





#endif /* end of include guard: PLANETS_HPP */
