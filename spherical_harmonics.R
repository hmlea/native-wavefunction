source("associated_legendre.R")

# formulas for spherical harmonics were found in this paper: https://www.sciencedirect.com/science/article/pii/S0166128097001851
# all this was verified with wolfram and evenescene: https://github.com/al2me6/evanescence/blob/dev/evanescence_core/mathematica/real_spherical_harmonics.json

# inspiration taken from https://github.com/nishanmudalige/s.harmonic/blob/master/R/s.harmonic.R
# according to wolfram and mathematica, (what I have tried) are correct (as well as evenescense)
# these are for the complex spherical harmonics
sph_harm = function(l, m, theta, phi) {
  if (m < 0) {
    norm_factor = sqrt(((2*l+1)/(4*pi)) * (factorial(l-m)/factorial(l+m)))
    output = norm_factor * assoc_legendre(l, m, cos(theta)) * exp(sqrt(as.complex(-1))*m*phi)
  } else if (m == 0) {
    norm_factor = sqrt(((2*l+1)/(4*pi)))
    output = norm_factor * assoc_legendre(l, m, cos(theta))
  } else {
    norm_factor = sqrt(((2*l+1)/(4*pi)) * (factorial(l-m)/factorial(l+m)))
    output = norm_factor * assoc_legendre(l, m, cos(theta)) * exp(sqrt(as.complex(-1))*m*phi)
  }
}

# these are for the Y*lm(theta, phi) in the paper
# this is the complex conjugate of the spherical harmonics needed to get the real values
sph_harm_conj = function(l, m, theta, phi) {
  phase = (-1)**m
  m_dash = -m
  output = phase * sph_harm(l, m_dash, theta, phi)
}

# these calculate the normalized real spherical harmonics
# this is calculated the way done in the paper with the normalization factor and conjugates
  # see the bottom of page 2 in the paper
real_sph_harm = function(l, m, theta, phi) {
  if (m > 0) {
    norm_factor = ((-1)**m)/(sqrt(2))
    complex_harm = sph_harm(l, m, theta, phi)
    conj_harm = sph_harm_conj(l, m, theta, phi)
    real_harm = norm_factor*(complex_harm+conj_harm)
  } else if (m < 0) {
    norm_factor = ((-1)**m)/(sqrt(2)*sqrt(as.complex(-1)))
    m_abs = abs(m)
    complex_harm = sph_harm(l, m_abs, theta, phi)
    conj_harm = sph_harm_conj(l, m_abs, theta, phi)
    real_harm = norm_factor*(complex_harm-conj_harm)
  } else {
    real_harm = sph_harm(l, m, theta, phi)
  }
  
  # after the caluclations above, the imaginary part is always 0
  # so you can just take the real value here
  output = Re(real_harm)
}

# NOTE that:
  # -l <= m <= l
  # theta [0, pi] radians OR [0, 180] degrees
  # phi [0, 2pi) radians OR [0, 360) degrees

