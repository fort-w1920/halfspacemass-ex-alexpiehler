train_depth <- function(data, # data to train halfspace depth
                        n_halfspace, # number of halfspaces, which shall be sampled
                        subsample, # fraction of of data which shall be used in each iteration
                        scope = 1,
                        plot_yes = FALSE,
                        seed = 111191) {
  ### Input Checks ###

  # if (checkmate::test_list(data)) {
  #   data <- try(as.data.frame(data))
  # }

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

  column_names <- colnames(data)
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

  # Create a list with the length number of halfspaces that should be sampled. Each entry of the list contains the same numeric value: the dimensions of the data (# features). The function 'get_unit_vector' is applied onto each entry of the list. This is done because the 'get_unit_vector' function takes as input a dimension and outputs a unit vector of the same dimension.
  # This is done because the
  unit_vector_list <- purrr::map(
    rep(list(ncol(data)), times = n_halfspace),
    get_unit_vector
  )

  # Create a list with the length number of halfspaces that should be sampled. Eacth list entry initially contains the full dataset. Then to each list entry the function 'sample_without_replacement' is applied, to obtain a subsample of the data with a prespecified fraction.
  subset_list <- purrr::map(rep(list(data), times = n_halfspace),
    sample_without_replacement,
    fraction = subsample
  )

  # Each entry of the subset_list, containing the sampled subsets, is multiplied with each entry of the unit_vector_list. As the the subsampled data of one list entry has the dimensions (n x p) and the the unit vector, which is multiplied from the right has the dimensions (p x 1) for every list entry (n_halfspaces), we obtain a matrix (n x n_halfspaces). A matrix is chosen as output, as this facilities easier computation in the following.
  projections <- mapply(`%*%`, subset_list, unit_vector_list)

  # For every sampled halfspace the maximum, the minimum and the midpoint is computed
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

  # The projections matrix (dimensions : (n x n_halfspace)) is split columnwise into a list. We obtain a list with a length accoring to the number of sampled halfspaces. Every entry of that list contains the projections of one subsampled dataset (a vector with length accoring to the sampled observations). This representation is required for the subsequent 'get_mass' function.
  projections_list <-
    split(projections, rep(1:ncol(projections), each = nrow(projections)))


  mass_left <- get_mass(
    cutoffs = cutoff_points,
    projections = projections_list,
    side = "left",
    subsample_size = floor(nrow(data) * subsample)
  )

  mass_right <- get_mass(
    cutoffs = cutoff_points,
    projections = projections_list,
    side = "right",
    subsample_size = floor(nrow(data) * subsample)
  )

  if (plot_yes) {
    if (ncol(data) != 2) {
      warning(
        "Results can only be plotted for 2-dimensional data. Returning trained model without plot"
      )
      return(
        list(
          "cutoff_points" = unlist(cutoff_points),
          "hyperplanes" = unit_vector_list,
          "mass_left" = unlist(mass_left),
          "mass_right" = unlist(mass_right)
        )
      )
    }

    minima <- apply(data, 2, min) - apply(data, 2, sd)
    maxima <- apply(data, 2, max) + apply(data, 2, sd)

    grid_for_plot <- expand.grid(
      x = seq(minima[1], maxima[1], length.out = 50),
      y = seq(minima[2], maxima[2], length.out = 50)
    )

    projections_per_datum <- make_projections_per_datum(
      data = grid_for_plot,
      hyperplanes = unit_vector_list
    )

    halfspace_mass_grid <- get_halfspace_mass(
      projections_per_datum = projections_per_datum,
      cutoff_points = unlist(cutoff_points),
      mass_left = unlist(mass_left),
      mass_right = unlist(mass_right)
    )

    data_for_plot <- as.data.frame(cbind(grid_for_plot, halfspace_mass_grid))
    colnames(data_for_plot) <- c("v1", "v2", "halfspacemass")

    spectralcolors <- c(
      "darkblue", "blue", "cyan", "lightgreen",
      "yellow", "orange", "red", "darkred"
    )

    halfspace_mass_plot <-
      ggplot(data_for_plot, aes(x = v1, y = v2)) +
      geom_tile(aes(fill = halfspacemass, colour = halfspacemass)) +
      scale_fill_gradientn("mass", colors = spectralcolors) +
      scale_colour_gradientn("mass", colors = spectralcolors) +
      # geom_point(data = as.data.frame(data), aes(as.name(column_names[1])), as.name[column_names[2]])
      theme_minimal()

    print(halfspace_mass_plot)

    list(
      "cutoff_points" = unlist(cutoff_points),
      "hyperplanes" = unit_vector_list,
      "mass_left" = unlist(mass_left),
      "mass_right" = unlist(mass_right),
      "plot" = halfspace_mass_plot
    )
  } else {
    list(
      "cutoff_points" = unlist(cutoff_points),
      "hyperplanes" = unit_vector_list,
      "mass_left" = unlist(mass_left),
      "mass_right" = unlist(mass_right)
    )
  }
}

evaluate_depth <- function(data, halfspaces) {
  if (checkmate::test_matrix(data)) {
    data <- try(as.data.frame(data))
  }


  checkmate::assert_data_frame(data,
    types = c("logical", "integer", "integerish", "double", "numeric"),
    min.rows = 1,
    min.cols = 1
  )

  projections_per_datum <- make_projections_per_datum(
    data = data,
    hyperplanes = halfspaces[["hyperplanes"]]
  )

  halfspace_mass <- get_halfspace_mass(
    projections_per_datum = projections_per_datum,
    cutoff_points = halfspaces[["cutoff_points"]],
    mass_left = halfspaces[["mass_left"]],
    mass_right = halfspaces[["mass_right"]]
  )

  halfspace_mass
}

get_unit_vector <- function(dimensions) {
  # Method according to Hicks, J. S. ad Wheeling, R. F. "An Efficient Method for
  # Generating Uniformly Distributed Points on the Surface of an n-Dimensional
  # Sphere." Comm. Assoc. Comput. Mach. 2, 13-15, 1959.

  randn <- rnorm(dimensions)
  normalized_vector <- matrix(randn, ncol = 1) * (1 / norm(randn, type = "2"))
  normalized_vector
}

sample_without_replacement <- function(data,
                                       fraction) {
  rows_to_sample <- sample(nrow(data), size = floor(nrow(data) * fraction))

  data[rows_to_sample, , drop = FALSE]
}

get_cutoff <- function(max,
                       min,
                       mid,
                       scope,
                       output_list = TRUE) { # Computes the cutoff per sampled halfspace.

  lower_sampling_bound <- mid - (0.5 * scope * (max - min))
  upper_sampling_bound <- mid + (0.5 * scope * (max - min))

  if (output_list) {
    cutoff_points <- list(
      rep(1, length(lower_sampling_bound)),
      lower_sampling_bound,
      upper_sampling_bound
    ) %>% pmap(runif)
  } else {
    cutoff_points <- matrix(runif(
      length(lower_sampling_bound),
      lower_sampling_bound,
      upper_sampling_bound
    ),
    nrow = 1
    )
  }

  cutoff_points
}

get_mass <- function(cutoffs,
                     projections,
                     side = c("left", "right"),
                     subsample_size) {
  points_in_halfspace <- switch(
    side,
    "left" = purrr::map2(.x = projections, .y = cutoffs, `<`),
    "right" = purrr::map2(.x = projections, .y = cutoffs, `>=`)
  )

  mass <- purrr::map2(
    lapply(points_in_halfspace, sum),
    subsample_size,
    `/`
  )

  mass
}

get_halfspace_mass <- function(projections_per_datum,
                               cutoff_points,
                               mass_left,
                               mass_right) {
  mass_matrix <- matrix(NA,
    nrow = length(projections_per_datum[[1]]),
    ncol = length(projections_per_datum)
  )

  for (datum in seq_len(length(projections_per_datum))) {
    mass_matrix[, datum] <-
      ifelse(projections_per_datum[[datum]] < cutoff_points,
        mass_left,
        mass_right
      )
  }

  halfspace_sums <- apply(mass_matrix, 2, sum)
  halfspace_mass <- halfspace_sums / nrow(mass_matrix)
  halfspace_mass
}

make_projections_per_datum <- function(data, hyperplanes) {
  data_list <- purrr::modify(asplit(data, MARGIN = 1), as.matrix) # Convert data matrix into a list with length accoring to the amount of observations. Every entry will be a d-dimensional vector.
  hyperplanes_transposed <- purrr::modify(hyperplanes, t) # Transpose the sampled hyperplanes for subsequent computation of scalar product

  projections_per_datum <- vector(mode = "list", length = length(data_list)) # A empty list is initiated with the length according to the observations

  # For each observation (a [1 x d]-vector), the scalar product is computed with every sampled hyperplane (a [d x 1]-vector). The output per interation is a vector containing all scaler products with current observation and all the scalar products (a [1 x n_halfspaces]-vector).
  for (datum in seq_len(length(data_list))) {
    projections_per_datum[[datum]] <-
      purrr::map_dbl(hyperplanes_transposed, `%*%`, data_list[[datum]])
  }

  projections_per_datum
}
