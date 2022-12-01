library(mpoly)

# radial component of the wave function
# this might need to be moved to its own file one day like "sph_harm.R
# according to https://keisan.casio.com/exec/system/1224054805, these values are correct
# THIS is the MAIN function
radial_comp = function(n, l, r) {
  # z is atomic number, might be able to change later
  Z = 1
  p = (2*Z) / n
  pr = p*r
  
  radical = (p**3) * (factorial(n-l-1))/(2*n*factorial(n+l))
  
  # associated laguerre function from mpoly
  # EVENTUALLY i would like to do these calculations in house
  assoc_lag = as.function(laguerre((n-l-1), alpha=(2*l+1)))
  
  Rr = sqrt(radical) * exp(-(pr/2)) * pr**l * assoc_lag(pr)
}

# radial component of the wave function taking into consideration the bohr radius
# in other programs the bohr radius seems to be left out:
# https://github.com/al2me6/evanescence/blob/f8d448de503f431e5a16548f7621a75932b4d64c/evanescence_core/src/orbital/atomic.rs
# https://dpotoyan.github.io/Chem324/H-atom-wavef.html
# https://keisan.casio.com/exec/system/1224054805
bohr_radial_comp = function(n, l, r) {
  # z is atomic number, might be able to change later
  Z = 1
  a0 = 0.0529 # in nm, the first bohr radius
  p = (2*Z) / (n*a0)
  pr = p*r
  
  radical = (p**3) * (factorial(n-l-1))/(2*n*factorial(n+l))
  
  # associated laguerre function from mpoly
  assoc_lag = as.function(laguerre((n-l-1), alpha=(2*l+1)))
  
  Rr = sqrt(radical) * exp(-(pr/2)) * pr**l * assoc_lag(pr)
}

