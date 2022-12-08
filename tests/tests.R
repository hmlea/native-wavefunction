source("tests/core.R")

# this is the same as the default seed - if no param, it takes UNIX system time
# this seed was generated with gen_seed() at 2022-12-08 13:31:15 EST
# it can be reproduced with gen_seed("2022-12-08 13:31:15 EST")
seed_val = 1670524275.162657

# generate tests
t1 = gen_tests(num_tests=1000, seed=seed_val) # equal to t2
t2 = gen_tests(num_tests=1000, seed=seed_val) # equal to t1
t3 = gen_tests(num_tests=1000)                # unique

# summarize test results and get the tolerance for all values to be equal
tol1 = summarize_tests(t1)
tol2 = summarize_tests(t2)
tol3 = summarize_tests(t3)

# format and print the tests and tolerances
message("in respect to the wavefunction package:\n",
        "t1 is equal with a tolerance of ", tol1, "\n",
        "t2 is equal with a tolerance of ", tol2, "\n",
        "t3 is equal with a tolerance of ", tol3, "\n")

