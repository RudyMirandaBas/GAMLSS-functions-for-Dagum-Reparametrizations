library(VGAM)

find_delta <- function(model, theta, censorship, search_interval = c(0, 100)) {
  # Extract parameters
  mu    <- theta[1]
  sigma <- theta[2]
  nu    <- theta[3]
  a     <- censorship

  # Define the objective function whose root we want to find.
  # The goal is to find the 'delta' such that objective_function(delta) = 0.
  objective_function <- function(delta) {
    # Compute the integral
    integral_value <- integrate(
      function(x) { exp(-delta * x) * model(x, mu, sigma, nu) },
      lower = 0,
      upper = Inf # Use Inf for integration to infinity
    )$value
    
    # The equation to solve is: 1 - integral - a = 0
    # Note: The subtraction of 'a' is now OUTSIDE the integral, which is the most common formulation.
    return(1 - integral_value - a)
  }
  
  # Use uniroot() to find the root (the value of delta) within an interval.
  # uniroot searches for the point where the objective function crosses zero.
  solution <- tryCatch({
    uniroot(objective_function, interval = search_interval)
  }, error = function(e) {
    # Handle the case where no root is found within the interval
    cat("Error in uniroot: Could not find a root in the specified interval.\n")
    cat("Original error message:", e$message, "\n")
    return(NULL)
  })
  
  # Return only the root that was found
  if (!is.null(solution)) {
    return(solution$root)
  } else {
    return(NA) # Return NA if no solution was found
  }
}
