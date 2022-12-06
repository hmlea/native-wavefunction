source("radial_component.R")
source("spherical_harmonics.R")

# ALL FUNCTIONS MUST BE CALLED WITH theta AND phi IN RADIANS
# for all (that I have tried) values of n and l with m=0, these have been confirmed
  # by wolfram (mathematica) and https://keisan.casio.com/exec/system/1224055293
# however, for non-zero values of m, none of these calculators or mine match up
  # i think this is because the complex nature of the sph_harm at non-zero m values
  # mine is calculating the real spher_harm based on paper by Blanco, Florez, and Bermejo
# NOTE that:
  # 0 < l < n-1
  # -l <= m <= l
  # theta [0, pi] radians OR [0, 180] degrees
  # phi [0, 2pi) radians OR [0, 360) degrees
psi = function(n, l, m, r, theta, phi) {
  radial = radial_comp(n, l, r)
  spherical = real_sph_harm(l, m, theta, phi)
  wave_func = radial * spherical
}

# converts degrees to radians
toRad = function(degrees) {
  degrees * (pi/180)
}

# converts radians to degrees
toDeg = function(radians) {
  radians * (180/pi)
}

