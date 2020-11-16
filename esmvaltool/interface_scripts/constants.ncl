; #############################################################################
; GENERAL COLLECTION OF PHYSICAL CONSTANTS
; #############################################################################

begin

  source = "http://physics.nist.gov/cuu/Constants/Table/allascii.txt"

  ; Avogadro number [mol^-1]
  Navo = 6.0221418e23
  Navo@source = source
  Navo@units = "mol^-1"
  Navo@name = "Avogadro number"

  ; Boltzmann constant [J K^-1]
  Kbol = 1.3806504e-23
  Kbol@source = source
  Kbol@units = "J K^-1"
  Kbol@name = "Boltzmann constant"

  ; Gas constant [J mol^-1 K^-1]
  Rgas = 8.314472
  Rgas@source = source
  Rgas@units = "J mol^-1 K^-1"
  Rgas@name = "Gas constant"

  ; Newtonian constant of gravitation [m^3 kg^-1 s^-2]
  Gnewt = 6.67428e-11
  Gnewt@source = source
  Gnewt@units = "m^3 kg^-1 s^-2"
  Gnewt@name = "Newtonian constant of gravitation"

  ; Speed of light in vacuum [m s^-1]
  Clig = 2.9979246e08
  Clig@source = source
  Clig@units = "m s^-1"
  Clig@name = "Speed of light in vacuum"

  ; Standard acceleration of gravity [m s^-2]
  Grav = 9.80665
  Grav@source = source
  Grav@units = "m s^-2"
  Grav@name = "Standard acceleration of gravity"

  ; Standard atmosphere [Pa]
  p0stp = 101325.
  p0stp@units = "Pa"
  p0stp@name = "Standard atmosphere"

  ; Standard temperature [K]
  T0stp = 273.15
  T0stp@units = "K"
  T0stp@name = "Standard temperature"

  ; Stefan-Boltzmann constant [W m^-2 K^-4]
  Sstebol = 5.670400e-08
  Sstebol@source = source
  Sstebol@units = "W m^-2 K^-4"
  Sstebol@name = "Stefan-Boltzmann constant"

  ; Pi [1]
  PI = 3.1415927
  PI@units = "1"
  PI@name = "Greek Pi"

  ; Earth radius [m]
  EARTHradius = 6371000.
  EARTHradius@units = "m"
  EARTHradius@name = "Earth radius"

  ; Earth surface area [m^2]
  EARTHarea = 4. * PI * EARTHradius ^ 2
  EARTHarea@units = "m2"
  EARTHarea@name = "Earth surface area"

  ; Molecular weight of the atmosphere [g mol^-1]
  MWair = 28.97
  MWair@source = "Seinfeld & Pandis - Atmospheric Chemistry " + \
    "and Physics - Wiley 2006"
  MWair@units = "g mol^-1"
  MWair@name = "Molecular weight of the atmosphere"

  ; Molecular weight of ozone [g mol^-1]
  MWh2o = 18.01528
  MWh2o@units = "g mol^-1"
  MWh2o@name = "Molecular mass of h2o"

  ; Molecular weight of ozone [g mol^-1]
  MWo3 = 47.9982
  MWo3@source = "http://pubchem.ncbi.nlm.nih.gov"
  MWo3@units = "g mol^-1"
  MWo3@name = "Molecular weight of ozone"

  ; Dobson unit
  DU = 2.1415e-05
  DU@units = "kg m^-2"
  DU@name = "conversion factor 1 Dobson unit = 2.1415E-5 kg m^-2"

end
