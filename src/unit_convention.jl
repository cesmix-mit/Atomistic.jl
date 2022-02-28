# The unit convention throughout the package and other packages in the CESMIX Julia Suite
# is to assume all unspecified units to be atomic units as defined in UnitfulAtomic.jl.

# Here we provide constants for the atomic units we use to make the code more readable.

MASS_UNIT = u"me_au"
LENGTH_UNIT = u"bohr"
ENERGY_UNIT = u"hartree"
TIME_UNIT = u"ħ_au / hartree"
VELOCITY_UNIT = u"bohr * hartree / ħ_au"
TEMPERATURE_UNIT = u"hartree / k_au"
FORCE_UNIT = u"hartree / bohr"

MASS_TYPE = typeof(zeros() * MASS_UNIT)
LENGTH_TYPE = typeof(zeros() * LENGTH_UNIT)
ENERGY_TYPE = typeof(zeros() * ENERGY_UNIT)
TIME_TYPE = typeof(zeros() * TIME_UNIT)
VELOCITY_TYPE = typeof(zeros() * VELOCITY_UNIT)
TEMPERATURE_TYPE = typeof(zeros() * TEMPERATURE_UNIT)
FORCE_TYPE = typeof(zeros() * FORCE_UNIT)
