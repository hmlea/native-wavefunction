# sources
  # [1] https://doi.org/10.1590/1806-9126-RBEF-2019-0073
  # [2] https://doi.org/10.1021/ed300393s

source("wave_function.R")
library(plotly)

# monte carlo simulator, checks and takes points less than max_psi
# from Henrique Fernandes Lobo et. al [1]
  # psi_lim taken from Tully et al. 2013
monte_carlo = function(n, l, m, psi_lim=0.0026, r_lim=35, num_points=2048) {
  # init data vectors
  r_points = c()
  theta_points = c()
  phi_points = c()
  psi = c()
  
  # init a progress bar
  message("calculating psi")
  prog = txtProgressBar(0, num_points, style=3, char="=")
  # loop until a sufficient number of points is collected
  counter = 0
  while (counter < num_points) {
    # generate a random radius, theta, and phi (spherical coords)
    r = runif(1, 0, r_lim)
    theta = runif(1, 0, pi)
    phi = runif(1, 0, 2*pi)
    
    # get the psi and |psi|^2, and a w
    psi_val = psi(n, l, m, r, theta, phi)
    psi2 = abs(psi_val)**2
    w = runif(1, 0, psi_lim)
    
    # if w is less than |psi|^2...
    if (w <= psi2) {
      # keep the point coordinates and psi value
      r_points = append(r_points, r)
      theta_points = append(theta_points, theta)
      phi_points = append(phi_points, phi)
      psi = append(psi, psi_val)
      
      # increment the counter and progress bar
      counter = counter + 1
      setTxtProgressBar(prog, counter)
      # print the counter (just to see how close is is)
      # print(counter)
    }
  }
  # close progress bar
  close(prog)
  rm(prog)
  
  # output the data points
  output = data.frame("r"=r_points, "theta"=theta_points, "phi"=phi_points, "psi"=psi)
}

# converts spherical coordinates to cartesian
convert_sph_coords = function(data) {
  # get spherical coordinates
  r = data[[1]]
  theta = data[[2]]
  phi = data[[3]]
  
  # convert to cartesian
  x = (r * sin(theta) * cos(phi))
  y = (r * sin(theta) * sin(phi))
  z = (r * cos(theta))
  
  # output the cartesian coordinates
  output = data.frame("x"=x, "y"=y, "z"=z, psi=data[[4]])
}

# returns a plot of the orbital data with color based on the psi value
plot_data = function(coords, pos_col="#d13010",
                     neg_col="#2d709f", bg_col="gray") {
  # get limit for the colorbar
  lim = max(max(coords$psi), abs(min(coords$psi)))
  # create the base figure
  fig = plot_ly(type="scatter3d", mode="markers")
  # add the scatter orbital to the figure
  fig = add_trace(fig, data=coords, x=~x, y=~y, z=~z, type="scatter3d",
                  mode="markers", size=1,
                  marker=list(size=1,
                              color=~psi,
                              colorscale=list(c(0, neg_col),
                                              c(0.5, "#ffffff"),
                                              c(1, pos_col)),
                              cauto=FALSE,
                              cmin=-lim,
                              cmaz=lim,
                              showscale=TRUE))
  # set the background color and hide the trace legend
  fig = layout(fig, paper_bgcolor=bg_col, showlegend=FALSE)
  # return the figure
  fig
}

sim = function(n, l, m, max_psi=0.0026, num_points=2048, r_lim=35,
               pos_col="#d13010", neg_col="#2d709f", bg_col="gray",
               benchmark=TRUE, plot_orb=TRUE) {
  # start the time benchmark if param is true
  if(benchmark) {
    start_time = proc.time()
  }
  
  # get the spherical coords from teh wave fucntion
  spherical_coords = monte_carlo(n, l, m, max_psi, num_points=2048, r_lim=35)
  # convert spherical coords to cartesian
  cart_coords = convert_sph_coords(spherical_coords)
  # plot the orbital data if param is true
  if(plot_orb) {
    fig = plot_data(cart_coords, pos_col, neg_col, bg_col)
  }
  
  # end the benchmark if param is true
  if(benchmark){
    print(proc.time() - start_time)
  }
  
  # return the figure or the cartesian coords depending on params
  if(plot_orb) {
    fig
  } else {
    cart_coords
  }
}

