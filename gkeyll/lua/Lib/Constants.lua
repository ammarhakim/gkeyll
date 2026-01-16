-- Gkyl ------------------------------------------------------------------------
--
-- Mathematical and physical constants in SI units. Taken from NIST
-- co-data website.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local _M = {}

_M.PI = 3.141592653589793238462643383279502884
_M.E = 2.718281828459045235360287471352662497
_M.SPEED_OF_LIGHT = 299792458.0 -- m/s
_M.PLANCKS_CONSTANT_H = 6.62606896e-34 -- joule*seconds
_M.ELECTRON_MASS = 9.10938215e-31 -- Kg
_M.PROTON_MASS = 1.672621637e-27 -- Kg
_M.MASS_UNIT = 1.66053907e-27 -- Kg
_M.ELEMENTARY_CHARGE = 1.602176487e-19 -- Coulomb
_M.BOLTZMANN_CONSTANT = 1.3806488e-23
_M.EPSILON0 = 8.854187817620389850536563031710750260608e-12 -- farad/meter
_M.MU0 = 12.56637061435917295385057353311801153679e-7 -- newtons/ampere/ampere
_M.EV2KELVIN = _M.ELEMENTARY_CHARGE/_M.BOLTZMANN_CONSTANT

return _M
