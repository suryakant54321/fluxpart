from collections import namedtuple

MolecularWeight = namedtuple('MolecularWeight', 'dryair vapor co2')
SpecificGasConstant = namedtuple('SpecificGasConstant', 'dryair vapor co2')
SpecificHeatCapacity = namedtuple('SpecificHeatCapacity', 'dryair')

GRAVITY = 9.81             # m/s^2
VON_KARMAN = 0.40          # unitless
UNIVERSAL_GAS = 8.3144598  # J/K/mol

# kg/mol
MOLECULAR_WEIGHT = MolecularWeight(
    dryair=0.0289645,
    vapor=0.018016,
    co2=0.044010)

# J/K/mol == m^3 Pa/kg/K
SPECIFIC_GAS_CONSTANT = SpecificGasConstant(
    dryair=UNIVERSAL_GAS / MOLECULAR_WEIGHT.dryair,
    vapor=UNIVERSAL_GAS / MOLECULAR_WEIGHT.vapor,
    co2=UNIVERSAL_GAS / MOLECULAR_WEIGHT.co2)

# J/K/kg
SPECIFIC_HEAT_CAPACITY = SpecificHeatCapacity(
    dryair=1004.67)
