#ifndef UNITS_HEADER
#define UNITS_HEADER

#include <boost/units/quantity.hpp>
#include <boost/units/systems/si.hpp>
#include <boost/units/cmath.hpp>

namespace boost
{
  namespace mpl
  {
    template<>
    struct divides_impl<boost::units::scale_dim_tag,boost::units::detail::static_rational_tag>
    {
      template<class T0, class T1>
      struct apply
      {
        typedef boost::units::scale_list_dim<
          boost::units::scale<
            (T0::base),
            typename mpl::divides<typename T0::exponent, T1>::type
          >
        > type;
      };
    };
  }
}

namespace photonflow {

namespace bu = boost::units;

using MicrometerBaseUnit = bu::scaled_base_unit<bu::si::meter_base_unit, bu::scale<10, bu::static_rational<-6>>>;
using MicrosecondBaseUnit = bu::scaled_base_unit<bu::si::second_base_unit, bu::scale<10, bu::static_rational<-6>>>;

using PhotonflowUnitSystem = bu::make_system<MicrometerBaseUnit, MicrosecondBaseUnit>::type;

using LengthUnit = bu::unit<bu::length_dimension, PhotonflowUnitSystem>;
using AreaUnit = bu::unit<bu::area_dimension, PhotonflowUnitSystem>;
using VolumeUnit = bu::unit<bu::volume_dimension, PhotonflowUnitSystem>;
using TimeUnit = bu::unit<bu::time_dimension, PhotonflowUnitSystem>;
using DimensionlessUnit = bu::unit<bu::dimensionless_type, PhotonflowUnitSystem>;

BOOST_UNITS_STATIC_CONSTANT(micrometer, LengthUnit);
BOOST_UNITS_STATIC_CONSTANT(micrometers, LengthUnit);
BOOST_UNITS_STATIC_CONSTANT(micrometre, LengthUnit);
BOOST_UNITS_STATIC_CONSTANT(micrometres, LengthUnit);

BOOST_UNITS_STATIC_CONSTANT(microsecond, TimeUnit);
BOOST_UNITS_STATIC_CONSTANT(microseconds, TimeUnit);

using Length = bu::quantity<LengthUnit>;
using Dimensionless = bu::quantity<DimensionlessUnit>;
using Area = bu::quantity<AreaUnit>;
using Volume = bu::quantity<VolumeUnit>;
using Time = bu::quantity<TimeUnit>;

// Literals

namespace literals {

#define BOOST_UNITS_LITERAL(suffix, unit, val, prefix, multiplier) \
inline bu::quantity<unit, double> operator "" _##prefix##suffix(long double x) \
{ \
    return bu::quantity<unit, double>(x * multiplier * val); \
} \
inline bu::quantity<unit, unsigned long long> operator "" _##prefix##suffix(unsigned long long x) \
{ \
    return bu::quantity<unit, unsigned long long>(x * multiplier * val); \
}


#define BOOST_UNITS_LITERAL_SET(suffix, unit, val) \
BOOST_UNITS_LITERAL(suffix, unit, val, Y, 1000000000000000000000000.0) \
BOOST_UNITS_LITERAL(suffix, unit, val, Z, 1000000000000000000000.0) \
BOOST_UNITS_LITERAL(suffix, unit, val, E, 1000000000000000000.0) \
BOOST_UNITS_LITERAL(suffix, unit, val, P, 1000000000000000.0) \
BOOST_UNITS_LITERAL(suffix, unit, val, T, 1000000000000.0) \
BOOST_UNITS_LITERAL(suffix, unit, val, G, 1000000000.0) \
BOOST_UNITS_LITERAL(suffix, unit, val, M, 1000000.0) \
BOOST_UNITS_LITERAL(suffix, unit, val, k, 1000.0) \
BOOST_UNITS_LITERAL(suffix, unit, val, h, 100.0) \
BOOST_UNITS_LITERAL(suffix, unit, val, da, 10.0) \
BOOST_UNITS_LITERAL(suffix, unit, val, , 1.0) \
BOOST_UNITS_LITERAL(suffix, unit, val, d, 0.1) \
BOOST_UNITS_LITERAL(suffix, unit, val, c, 0.01) \
BOOST_UNITS_LITERAL(suffix, unit, val, m, 0.001) \
BOOST_UNITS_LITERAL(suffix, unit, val, u, 0.000001) \
BOOST_UNITS_LITERAL(suffix, unit, val, n, 0.00000001) \
BOOST_UNITS_LITERAL(suffix, unit, val, p, 0.00000000001) \
BOOST_UNITS_LITERAL(suffix, unit, val, f, 0.00000000000001) \
BOOST_UNITS_LITERAL(suffix, unit, val, a, 0.00000000000000001) \
BOOST_UNITS_LITERAL(suffix, unit, val, z, 0.00000000000000000001) \
BOOST_UNITS_LITERAL(suffix, unit, val, y, 0.00000000000000000000001)

BOOST_UNITS_LITERAL_SET(m, photonflow::LengthUnit, 1000*1000*photonflow::micrometer)
//BOOST_UNITS_LITERAL_SET(g, si::mass, 0.001 * si::kilogram)
BOOST_UNITS_LITERAL_SET(s, photonflow::TimeUnit, 1000*1000*photonflow::microsecond)
//BOOST_UNITS_LITERAL_SET(A, si::current, si::ampere)
//BOOST_UNITS_LITERAL_SET(K, si::temperature, si::kelvin)
//BOOST_UNITS_LITERAL_SET(mol, si::amount, si::mole)
//BOOST_UNITS_LITERAL_SET(cd, si::luminous_intensity, si::candela)
//BOOST_UNITS_LITERAL_SET(Hz, si::frequency, si::hertz)
//BOOST_UNITS_LITERAL_SET(rad, si::plane_angle, si::radian)
//BOOST_UNITS_LITERAL_SET(sr, si::solid_angle, si::steradian)
//BOOST_UNITS_LITERAL_SET(N, si::force, si::newton)
//BOOST_UNITS_LITERAL_SET(Pa, si::pressure, si::pascal)
//BOOST_UNITS_LITERAL_SET(J, si::energy, si::joule)
//BOOST_UNITS_LITERAL_SET(W, si::power, si::watt)
//BOOST_UNITS_LITERAL_SET(C, si::electric_charge, si::coulomb)
//BOOST_UNITS_LITERAL_SET(V, si::electric_potential, si::volt)
//BOOST_UNITS_LITERAL_SET(F, si::capacitance, si::farad)
//BOOST_UNITS_LITERAL_SET(ohm, si::resistance, si::ohm)
//BOOST_UNITS_LITERAL_SET(S, si::conductance, si::siemens)
//BOOST_UNITS_LITERAL_SET(Wb, si::magnetic_flux, si::weber)
//BOOST_UNITS_LITERAL_SET(T, si::magnetic_flux_density, si::tesla)
//BOOST_UNITS_LITERAL_SET(H, si::inductance, si::henry)
//BOOST_UNITS_LITERAL_SET(degC, si::temperature, si::kelvin + 273.15 * si::kelvin)
//BOOST_UNITS_LITERAL_SET(lm, si::luminous_flux, si::lumen)
//BOOST_UNITS_LITERAL_SET(lx, si::illuminance, si::lux)
//BOOST_UNITS_LITERAL_SET(Bq, si::activity, si::becquerel)
//BOOST_UNITS_LITERAL_SET(Gy, si::absorbed_dose, si::gray)
//BOOST_UNITS_LITERAL_SET(Sv, si::dose_equivalent, si::sievert)
//BOOST_UNITS_LITERAL_SET(kat, si::catalytic_activity, si::katal)
//BOOST_UNITS_LITERAL_SET(min, si::time, 60.0 * si::second)
//BOOST_UNITS_LITERAL_SET(h, si::time, 60.0 * 60.0 * si::second)
//BOOST_UNITS_LITERAL_SET(day, si::time, 60.0 * 60.0 * 24.0 * si::second)
//BOOST_UNITS_LITERAL_SET(deg, si::plane_angle, M_PI / 180.0 * si::radian)
//BOOST_UNITS_LITERAL_SET(l, si::volume, 0.001 * si::cubic_meter)
//BOOST_UNITS_LITERAL_SET(L, si::volume, 0.001 * si::cubic_meter)
//BOOST_UNITS_LITERAL_SET(t, si::mass, 1000.0 * si::kilogram)

} // literals

std::ostream& operator<< (std::ostream& out, Length l);

} // namespace


#endif // UNITS

