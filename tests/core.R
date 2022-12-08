library(uuid)
source("monte_carlo.R")

# custom random number generator with weighted probabilities for larger nums:
  # 1/10 for 1
  # 2/10 for 2
  # 3/10 for 3
  # 4/10 for 4
rand_n = function(num, max_n) {
  cuts = c(0)
  for(i in 1:max_n) {
    cuts = append(cuts, tail(cuts, 1) + i)
  }
  
  # cuts and outputs the numbers
  nums = floor(runif(num, 1, max(cuts)))
  output = as.numeric(cut(nums, breaks=cuts))
}

# generate a seed based on the system time (UNIX time stamp)
gen_seed = function(date_time=NULL) {
  if(is.null(date_time)) date_time = Sys.time()
  as.numeric(as.POSIXct((date_time)))
}

# convert a seed to the time it was generated
seed_to_time = function(seed) {
  as.POSIXct(seed, origin="1970-01-01")
}

# generate test values: n, l, m, r, theta, phi, psi
gen_tests = function(num_tests=1000, seed=NULL, max_n=4, max_r=35) {
  # assign seed if need be
  if(is.null(seed)) {
    seed = gen_seed()
  }
  
  # set seed and unique object identifier (uoid)
  set.seed(seed)
  # uoids are meant to tell if objects that are equal are the same
  uoid = as.character(UUIDgenerate())
  seed_col = list(seed_uoid=c(as.character(seed), uoid,
                              rep(NA_character_, num_tests-2)))
  
  # init progress bar
  message("generating tests")
  counter = 0
  prog = txtProgressBar(0, num_tests+(num_tests/10), style=3, char="=")
  
  # generate random quantum numbers
  quant_nums = list(n=floor(runif(num_tests, 1, max_n+1)),
                    l=c(),
                    m=c())
  for(n_val in quant_nums$n) {
    l_val = floor(runif(1, 0, n_val))
    m_val = floor(runif(1, -l_val, l_val+1))
    quant_nums$l = append(quant_nums$l, l_val)
    quant_nums$m = append(quant_nums$m, m_val)
  }
  
  # increment progress bar after generating quantum numbers
  counter = counter + (num_tests / 10)
  setTxtProgressBar(prog, counter)
  
  # generate random spherical coordinates
  sph_coords = list(r=runif(num_tests, 0, 35),
                    theta=runif(num_tests, 0, pi),
                    phi=runif(num_tests, 0, 2*pi))
  
  # combine lists random quantum numbers and spherical coordinates
  vars = c(quant_nums, sph_coords)
  
  # generate the results list
  results = list(radial_component=c(),
                 associated_legendre=c(),
                 spherical_harmonics=c(),
                 wavefunction=c())
  for(i in seq_along(vars$r)) {
    # get current variables
    cur_n = vars$n[[i]]
    cur_l = vars$l[[i]]
    cur_m = vars$m[[i]]
    cur_r = vars$r[[i]]
    cur_theta = vars$theta[[i]]
    cur_phi = vars$phi[[i]]
    
    # calculate and append the orbital values
    results$radial_component = append(results$radial_component,
                                      radial_comp(cur_n, cur_l, cur_r))
    results$associated_legendre = append(results$associated_legendre,
                                         assoc_legendre(cur_l, cur_m, cos(cur_theta)))
    results$spherical_harmonics = append(results$spherical_harmonics,
                                         real_sph_harm(cur_l, cur_m, cur_theta, cur_phi))
    results$wavefunction = append(results$wavefunction,
                                  psi(cur_n, cur_l, cur_m, cur_r, cur_theta, cur_phi))
    
    # increment progress bar after calculating results
    counter = counter + 1
    setTxtProgressBar(prog, counter)
  }
  
  # close progress bar
  close(prog)
  rm(prog)
  
  # output a dataframe of the variables and the results
  output_df = as.data.frame(c(vars, results, seed_col))
  # change type of seed_uoid column
  output_df$seed_uoid = as.character(output_df$seed_uoid)
  output_df
}

# view the tests without the long obnoxious seed
view_tests = function(tests, print=TRUE, show_seed=FALSE) {
  if(!show_seed) display_vals = tests[1:10] else display_vals = tests
  if(print) print(display_vals) else View(display_vals)
}

# get the seed from tests
get_seed = function(tests) {
  seed = as.numeric(tests$seed_uoid[[1]])
}

# get the uoid from the tests
get_uoid = function(tests) {
  uoid = as.character(tests$seed_uoid[[2]])
}

# print the seed with all the 22 precision decimals as a string
  # by default all this precision is included in the int but R hides it
print_seed = function(tests=NULL, seed=NULL, precision=NULL) {
  # check for null params and respond
  if(is.null(tests) && is.null(seed)) stop("please supply tests or seed")
  else if(is.null(seed)) {
    seed = get_seed(tests)
  }
  
  # print the seed
  if(is.null(precision)) precision = nchar(seed) - 10
  sprintf(paste0("%.", precision, "f"), seed)
}

# check equality of the psi with the wavefunction package
check_equality = function(test_vals, tolerance=(10**-7)) {
  # calculate wavefunction using the wavefunction package (with Boost)
  wf_pkg_vals = mapply(wavefunction::psi,
                       test_vals$n,
                       test_vals$l,
                       test_vals$m,
                       test_vals$r,
                       test_vals$theta,
                       test_vals$phi)
  
  # check to see if wavefunction and native-wavefunction values are equal
  equal_apply = mapply(all.equal,
                       test_vals$wavefunction,
                       wf_pkg_vals,
                       tolerance=tolerance)
  
  # get the mean relative differences of each value as a number
  mean_dif = as.numeric(gsub("^\\D+", "", equal_apply))
  
  if(!all(is.na(mean_dif))) {
    # omit NA from mean_dif vector for statistical calculations
    omitted_dif = na.omit(mean_dif)
    
    # get stats on how different the numbers are
    num_dif_vals = length(omitted_dif)
    dif_val_index = which(!is.na(mean_dif))
    max_dif = max(omitted_dif)
    max_dif_index = which(mean_dif==max_dif)
    min_dif = min(omitted_dif)
    min_dif_index = which(mean_dif==min_dif)
    range_dif = max_dif - min_dif
    median_dif = median(omitted_dif)
    avg_dif = mean(omitted_dif)
    equal_val_index = which(is.na(mean_dif))
    num_equal_vals = length(equal_val_index)
    
  } else {
    # get stats on how different the numbers are
    num_dif_vals = 0
    dif_val_index = NA
    max_dif = NA
    max_dif_index = NA
    min_dif = NA
    min_dif_index = NA
    range_dif = NA
    median_dif = NA
    avg_dif = NA
    equal_val_index = which(is.na(mean_dif))
    num_equal_vals = length(equal_val_index)
  }
  
  # return a list on the equality stats
  list(total_tests=length(mean_dif),
       mapply_raw_return=equal_apply,
       equals_tolerance=tolerance,
       num_different=num_dif_vals,
       different_index=dif_val_index,
       max_difference=max_dif,
       max_index=max_dif_index,
       min_difference=min_dif,
       min_index=min_dif_index,
       difference_range=range_dif,
       median_difference=median_dif,
       average_difference=avg_dif,
       num_equal=num_equal_vals,
       equal_index=(if(num_equal_vals) equal_val_index else NA))
}

# summarize, print, and return the results of the tests
summarize_tests = function(tests, tolerance=10**-7, omit_index=TRUE) {
  # get a summary using by checking the equality with the above function
  summary = check_equality(tests, tolerance=tolerance)
  
  # get the tolerance needed for all values to be equal
  equal_tol = tolerance + summary$max_difference
  
  # message overview generic data (name, seed, uoid)
  message(paste0("TEST SUMMARY FOR \"", deparse(substitute(tests))), "\"\n")
  message(paste0("Seed: ", print_seed(tests=tests), "\n",
                 "UOID: ", get_uoid(tests)), "\n")
  
  # print the summary
  if(summary$num_different) {
    # send message if the tests are different in value
    message(paste0(summary$num_equal, "/", summary$total_tests,
                   " tests are exactly equal with a tolerance of ",
                   summary$equals_tolerance), "\n")
    message(paste0("The ", summary$num_different,
                   " tests that aren't exactly equal have an average relative difference of ",
                   summary$average_difference), "\n")
    message(paste0("The minimum average difference is ", summary$min_difference,
                   " and the maximum average difference is ",
                   summary$max_difference, " resulting in a range of ",
                   summary$difference_range, "\n"))
    message(paste0("The tolerance needed for all values to be considered equal is ",
                   equal_tol))
    
    # print indexes where values were unequal if given the parameter
    if(!omit_index) {
      message(paste0("The values that aren't exactly equal are at indexes ",
                     paste(summary$different_index, collapse=", ")))
    }
  } else {
    # send message if the tests are the same in value
    message(paste0("All ", length(summary$equal_index),
                  " values are exactly equal with a tolerance of ",
                  summary$equals_tolerance))
  }
  
  # return the full equality summary
  output = equal_tol
}

# returns if two tests are equal to each other
tests_are_equal = function(test1, test2) {
  # get number of rows
  test1_rows = nrow(test1)
  test2_rows = nrow(test2)
  
  # get uoids
  test1_uoid = get_uoid(test1)
  test2_uoid = get_uoid(test2)
  
  # remove uoids for equality check
  test1$seed_uoid = NA
  test2$seed_uoid = NA
  
  # if rows are inequal
  if(nrow(test1) != nrow(test2)) {
    if(test1_rows > test2_rows) {
      # check to see which test is larger
      greater_name = "test1"
      lesser_name = "test2"
    } else {
      greater_name = "test1"
      lesser_name = "test2"
    }
    
    # the the row number of the larger and smaller test
    greater_name_rows = get(paste0(greater_name, "_rows"))
    lesser_name_rows = get(paste0(lesser_name, "_rows"))
    
    # return a message saying that one was larger and by how much
    row_diff = greater_name_rows - lesser_name_rows
    msg = paste0(greater_name, " is larger than ", lesser_name, " by ",
                 row_diff, " rows")
    message(msg)
    
    # set return value to false
    return_bool = FALSE
  } else {
    # check if uoids are equal (same object)
    if(test1_uoid == test2_uoid) {
      message("test1 and test2 are the same object")
    }
    
    # otherwise check to see if they are all equal
    return_bool = all(na.omit(test1 == test2))
  }

  # return true or false
  return_bool
}

