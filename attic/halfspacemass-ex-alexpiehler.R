# Function tain_depth takes numeric data and returns a halfspaces-list, which
# can be used in the plot_depth() and evaluate_depth() function
train_depth <- function(data,
                        # data to train halfspace depth
                        n_halfspace,
                        # number of halfspaces, which shall be sampled
                        subsample,
                        # fraction of of data which shall be used in each
                        scope = 1,
                        # the sampling region for the splitpoint within the
                        # convex set
                        seed = 111191) {
  ### Required Packages ###
  suppressMessages(require(tidyverse))
  suppressMessages(require(checkmate))


  ### Input Checks ###

  if (checkmate::test_list(data)) {
    data <- try(as.data.frame(data))
  }

  if (checkmate::test_data_frame(data)) {
    checkmate::assert_data_frame(data,
      types = c(
        "logical",
        "integer",
        "integerish",
        "double",
        "numeric"
      )
    )
  }

  data <- as.matrix(data)

  checkmate::assert_matrix(
    data,
    any.missing = FALSE,
    all.missing = FALSE,
    min.rows = 2,
    min.cols = 2,
  )

  checkmate::assert_count(n_halfspace)
  checkmate::assert_number(subsample, lower = 0, upper = 1)
  checkmate::assert_number(scope, lower = 1)
  checkmate::assert_count(seed)

  ### Computation of halfspaces ###

  set.seed(seed)

  # Create a list with the length number of halfspaces that should be sampled.
  # Each entry of the list contains the same numeric value: the dimensions of
  # the data (# features). The function 'get_unit_vector' is applied onto each
  # entry of the list. This is done because the 'get_unit_vector' function takes
  # as input a dimension and outputs a unit vector of the same dimension.
  unit_vector_list <- purrr::map(
    rep(list(ncol(data)), times = n_halfspace),
    get_unit_vector
  )

  # Create a list with the length number of halfspaces that should be sampled.
  # Each list entry initially contains the full dataset. Then to each list entry
  # the function 'sample_without_replacement' is applied, to obtain a subsample
  # of the data with a prespecified fraction.
  subset_list <- purrr::map(rep(list(data), times = n_halfspace),
    sample_without_replacement,
    fraction = subsample
  )

  # Each entry of the subset_list, containing the sampled subsets, is multiplied
  # with each entry of the unit_vector_list. As the the subsampled data of one
  # list entry has the dimensions (n x p) and the the unit vector, which is
  # multiplied from the right has the dimensions (p x 1) for every list entry
  # (n_halfspaces), we obtain a matrix (n x n_halfspaces). A matrix is chosen as
  # output, as this facilities easier subsequent computation.
  projections <- mapply(`%*%`, subset_list, unit_vector_list)

  # For every sampled halfspace the maximum, the minimum and the midpoint is
  # computed
  max_per_iter <- apply(projections, 2, max)
  min_per_iter <- apply(projections, 2, min)
  mid_per_iter <- (max_per_iter + min_per_iter) / 2

  # Compute the cutoff points in conxex set per sampled halfspace
  cutoff_points <- get_cutoff(
    max = max_per_iter,
    min = min_per_iter,
    mid = mid_per_iter,
    scope = scope
  )

  # The projections matrix (dimensions : (n x n_halfspace)) is split columnwise
  # into a list. We obtain a list with a length accoring to the number of
  # sampled halfspaces. Every entry of that list contains the projections of one
  # subsampled dataset (a vector with length accoring to the sampled
  # observations). This representation is required for the subsequent 'get_mass'
  projections_list <-
    split(projections, rep(1:ncol(projections), each = nrow(projections)))

  # The relative frequency of the points lying on the left hand side per sampled
  # halfspace is obtained with the function 'mass_left'
  mass_left <- get_mass_left(
    cutoffs = cutoff_points,
    projections = projections_list
  )

  # Since the total_mass = mass_right + mass_left = 1, we can subtract mass_left
  # from 1 to obtain the right mass for every sampled hyperplane
  mass_right <- 1 - mass_left

  halfspaces <- list(
    "cutoff_points" = cutoff_points,
    "hyperplanes" = unit_vector_list,
    "mass_left" = mass_left,
    "mass_right" = mass_right
  )

  # Return the list 'halfspaces'
  halfspaces
}

# evaluate_depth() takes as input the computed list in train_depth() and some
# data with the same columns that train_depth() was used on and returns a vector
# of either the halfspacemass or the halfspacedepth for each observation.
# The function is used in the plot_depth() function
evaluate_depth <- function(data,
                           halfspaces,
                           type = c("halfspacemass", "halfspacedepth")) {
  ### Input checks and data conversions ###
  if (checkmate::test_matrix(data)) {
    data <- try(as.data.frame(data))
  }

  checkmate::assert_data_frame(
    data,
    types = c("logical", "integer", "integerish", "double", "numeric"),
    min.rows = 1,
    min.cols = 1
  )

  checkmate::assert_list(halfspaces,
    len = 4
  )

  checkmate::assert_subset(
    names(halfspaces),
    c(
      "hyperplanes",
      "cutoff_points",
      "mass_left",
      "mass_right"
    )
  )

  checkmate::assert_character(type)

  type <- match.arg(type, c("halfspacemass", "halfspacedepth"))

  ### Execution of evaluation of depth ###

  projections_per_datum <-
    make_projections_per_datum(
      data = data, hyperplanes = halfspaces[["hyperplanes"]]
    )

  halfspace_mass <- get_halfspace_mass(
    projections_per_datum = projections_per_datum,
    cutoff_points = halfspaces[["cutoff_points"]],
    mass_left = halfspaces[["mass_left"]],
    mass_right = halfspaces[["mass_right"]],
    type = type
  )

  # Return an object containing the halfspacemass (1 x n)-vector
  halfspace_mass
}

# get_unit_vector() is used in the train_depth() function and returns a unit
# vector for a given dimension
get_unit_vector <- function(dimensions) {
  # Method according to Hicks, J. S. ad Wheeling, R. F. "An Efficient Method for
  # Generating Uniformly Distributed Points on the Surface of an n-Dimensional
  # Sphere." Comm. Assoc. Comput. Mach. 2, 13-15, 1959.

  randn <- rnorm(dimensions)
  normalized_vector <-
    matrix(randn, ncol = 1) * (1 / norm(randn, type = "2"))
  normalized_vector
}

# sample_without_replacement() takes data and samples a fraction of that data
# without replacement. The function is used in train_depth()
sample_without_replacement <- function(data,
                                       fraction) {
  rows_to_sample <-
    sample(nrow(data), size = floor(nrow(data) * fraction))

  data[rows_to_sample, , drop = FALSE]
}

# Get cutoff takes the minimum, the maximum, and the mid point of a convex set
# and samples a random cutoff point within this convex set. The function is used
# in the train_depth() function
get_cutoff <- function(max,
                       min,
                       mid,
                       scope) {
  # The lower and upper sampling bound is computed. The dimensions mid/max/min
  # are (1 x n_halfspace)
  lower_sampling_bound <- mid - (0.5 * scope * (max - min))
  upper_sampling_bound <- mid + (0.5 * scope * (max - min))

  # The cutoff points are sampled in the lower and upper sampling bound.
  cutoff_points <- runif(
    length(lower_sampling_bound),
    lower_sampling_bound,
    upper_sampling_bound
  )

  # Return a vector with the cutoff_points for every sampled halfspace
  cutoff_points
}

# get_mass_left() computes the relative frequency of points, which lie left
# w.r.t. a given cutoff point in a 1-d convex set. This function is used in
# the train_depth() function
get_mass_left <- function(cutoffs,
                          projections) {
  points_in_halfspace <-
    purrr::map2(.x = projections, .y = cutoffs, `<`)

  mass <- unlist(lapply(points_in_halfspace, mean))

  mass
}

# get_halfspace_mass() is used to compute the halfspacemass or halfspacedepth
# for points of a dataset. It is used in the evaluate_depth() function.
get_halfspace_mass <- function(projections_per_datum,
                               cutoff_points,
                               mass_left,
                               mass_right,
                               type = c("halfspacemass", "halfspacedepth")) {
  type <- match.arg(type, c("halfspacemass", "halfspacedepth"))

  # Initiate a matrix that has (n_halfspaces x n)-dimensions
  mass_matrix <- matrix(
    NA,
    nrow = length(projections_per_datum[[1]]),
    ncol = length(projections_per_datum)
  )

  # For every observation a column in the matrix is filled up with values of
  # mass right or mass left. A matrix with values between 0 and 1 is obtained
  for (datum in seq_len(length(projections_per_datum))) {
    mass_matrix[, datum] <-
      ifelse(projections_per_datum[[datum]] < cutoff_points,
        mass_left,
        mass_right
      )
  }

  if (type == "halfspacemass") {
    halfspace_mass <- apply(mass_matrix, 2, mean)
  } else {
    halfspace_mass <- apply(mass_matrix, 2, min)
  }

  halfspace_mass
}

# make_projections_per_datum() is used to project a given dataset onto a set
# of sampled unit vectors. It is used in the evaluate_depth() function
make_projections_per_datum <- function(data, hyperplanes) {
  # Convert data matrix into a list with length accoring to the amount of
  # observations. Every entry will be a d-dimensional vector.
  data_list <- purrr::modify(asplit(data, MARGIN = 1), as.matrix)

  # Transpose the sampled hyperplanes for subsequent computation of scalar
  # product
  hyperplanes_transposed <- purrr::modify(hyperplanes, t)

  # A empty list is initiated with the length according to the observations
  projections_per_datum <-
    vector(mode = "list", length = length(data_list))

  # For each observation (a [1 x d]-vector), the scalar product is computed with
  # every sampled hyperplane (a [d x 1]-vector). The output per interation is a
  # vector containing all scaler products with current observation and all the
  # scalar products (a [1 x n_halfspaces]-vector).
  for (datum in seq_len(length(data_list))) {
    projections_per_datum[[datum]] <-
      purrr::map_dbl(hyperplanes_transposed, `%*%`, data_list[[datum]])
  }

  projections_per_datum
}

# plot_depth() uses the evaluate_depth() function to output a plot the
# halfspacemass or halfspacedepth of a certain 2-d dataset. The 'data' argument
# takes the original data, which was used to train the 'halfspaces' object.
# This is done to obtain a range for the x- and y-axis and to plot the
# datapoints. One can also choose to already specify a certain grid over which
# the data should be plotted.
plot_depth <- function(data,
                       halfspaces,
                       type = c("halfspacemass", "halfspacedepth"),
                       prespecified_grid = NULL) {
  ### Input checks ###

  checkmate::assert_list(halfspaces,
    len = 4
  )

  checkmate::assert_subset(
    names(halfspaces),
    c(
      "hyperplanes",
      "cutoff_points",
      "mass_left",
      "mass_right"
    )
  )

  if (!is.null(prespecified_grid)) {
    prespecified_grid <- try(as.data.frame(prespecified_grid))
    checkmate::assert_data_frame(
      prespecified_grid,
      types = c(
        "logical",
        "integer",
        "integerish",
        "double",
        "numeric"
      ),
      ncols = 2,
      min.rows = 1,
      any.missing = FALSE
    )
  }

  data <- try(as.data.frame(data))
  checkmate::assert_data_frame(
    data,
    types = c(
      "logical",
      "integer",
      "integerish",
      "double",
      "numeric"
    ),
    ncols = 2,
    min.rows = 1,
    any.missing = FALSE
  )

  checkmate::assert_character(type)
  type <- match.arg(type, c("halfspacemass", "halfspacedepth"))

  ### Computations for plot ###

  # Get the endpoints of the data to expand a grid of values
  minima <- apply(data, 2, min) - apply(data, 2, sd)
  maxima <- apply(data, 2, max) + apply(data, 2, sd)

  # If the data of a variable is a constant, get a range of 10 for this axis
  if (any(minima == maxima)) {
    loc <- which(minima == maxima)
    minima[loc] <- -5
    maxima[loc] <- 5
  }

  if (is.null(prespecified_grid)) {
    # Expand the grid of values
    grid_for_plot <- expand.grid(
      x = seq(minima[[1]], maxima[[1]], length.out = 50),
      y = seq(minima[[2]], maxima[[2]], length.out = 50)
    )
  } else {
    grid_for_plot <- prespecified_grid
  }

  halfspace_mass_grid <-
    evaluate_depth(
      data = grid_for_plot,
      halfspaces = halfspaces,
      type = type
    )

  data_for_plot <-
    as.data.frame(cbind(grid_for_plot, halfspace_mass_grid))

  # Obtain column names or assign column names if the were not specified
  if (is.null(colnames(data))) {
    colnames(data_for_plot) <- c("v1", "v2", type)
    colnames(data) <- c("v1", "v2")
  } else {
    colnames(data_for_plot) <- c(colnames(data), type)
  }

  ### Plot ###

  spectralcolors <- c(
    "darkblue",
    "blue",
    "cyan",
    "lightgreen",
    "yellow",
    "orange",
    "red",
    "darkred"
  )

  theme_plot <-
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )


  plot_val <-
    ggplot(data_for_plot, aes_string(
      x = colnames(data_for_plot)[[1]],
      y = colnames(data_for_plot)[[2]]
    )) +
    geom_tile(aes_string(
      fill = type,
      colour = type
    )) +
    scale_fill_gradientn(type, colors = spectralcolors) +
    scale_colour_gradientn(type, colors = spectralcolors) +
    geom_point(
      data = as.data.frame(data),
      aes_string(
        x = colnames(data)[[1]],
        y = colnames(data)[[2]]
      ),
      color = "white"
    ) +
    theme_minimal() +
    theme_plot

  plot_val
}