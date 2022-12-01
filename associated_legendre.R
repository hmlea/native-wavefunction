# formulas for legendres were found in this paper: https://faculty.washington.edu/seattle/physics227/reading/reading-26-27.pdf
# all this was verified with wolfram and evenescene: https://github.com/al2me6/evanescence/blob/dev/evanescence_core/mathematica/real_spherical_harmonics.json

# takes the mth derivative of a normal legendre polynomial as needed for the associated
# then to be used in the function below
derive_legendre = function(l, m, x_input) {
  coeff = ((-1)**m) / (2**l * factorial(l)) * (1-(x_input**2))**(m/2)
  leg_deriv = expression((x**2 - 1)**l)
  
  for (i in 1:(l+m)) {
    leg_deriv = D(leg_deriv, "x")
  }
  
  solved_deriv = eval(leg_deriv, envir=list(x=x_input))
  evaluated = coeff * solved_deriv
}

# values match with https://keisan.casio.com/exec/system/1180573406 and wolfram
# for wolfram alpha with l=3, m=-1, t=cos(1): LegendreP[3,-1,cos(1)] - same value
  # mathematica is the same as well
assoc_legendre = function(l, m, x) {
  if (m == 0 & l == 0) {
    output = 1
  } else if (m > 0) {
    output = derive_legendre(l, m, x)
  } else {
    m_pos = abs(m)
    Plm_pos = derive_legendre(l, m_pos, x)
    conversion_coeff = ((-1)**m_pos) * (factorial(l-m_pos)/factorial(l+m_pos))
    output = conversion_coeff * Plm_pos
  }
}

