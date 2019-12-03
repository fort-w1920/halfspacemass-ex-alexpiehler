train_depth <- function(data, # data to train halfspace depth
                        n_halfspace, # number of halfspaces, which shall be sampled
                        subsample, # fraction of of data which shall be used in each iteration
                        scope,
                        seed = 111191) {
  ### Input Checks ###

  if (checkmate::test_list(data)) {
    data <- try(as.data.frame(data))
  }

  if (checkmate::test_data_frame(data)) {
    checkmate::assert_data_frame(data,
      types = c("logical", "integer", "integerish", "double", "numeric")
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

  unit_vector_list <- purrr::map(
    rep(list(ncol(data)), times = n_halfspace),
    get_unit_vector
  )

  subset_list <- purrr::map(rep(list(data), times = n_halfspace),
    sample_without_replacement,
    fraction = subsample
  )


  projections <- mapply(`%*%`, subset_list, unit_vector_list)

  max_per_iter <- apply(projections, 2, max)
  min_per_iter <- apply(projections, 2, min)
  mid_per_iter <- (max_per_iter + min_per_iter) / 2

  cutoff_points <- get_cutoff(
    max = max_per_iter,
    min = min_per_iter,
    mid = mid_per_iter,
    scope = scope
  )

  projections_list <-
    split(projections, rep(1:ncol(projections), each = nrow(projections)))

  mass_left <- get_mass(
    cutoff_points = cutoff_points,
    projections = projections_list,
    side = "left",
    subsample_size = floor(nrow(data) * subsample)
  )

  mass_right <- get_mass(
    cutoff_points = cutoff_points,
    projections = projections_list,
    side = "right",
    subsample_size = floor(nrow(data) * subsample)
  )

  list(
    "cutoff_points" = cutoff_points,
    "halfspaces" = unit_vector_list,
    "mass_left" = mass_left,
    "mass_right" = mass_right
  )
}

get_unit_vector <- function(dimensions){
    randn_unif <- runif(dimensions, -1, 1)
    normalized_vector <- matrix(randn_unif, ncol = 1) * (1/norm(randn_unif, type = "2"))
    normalized_vector
  
}


sample_without_replacement <- function(data, fraction){
  rows_to_sample <- sample(nrow(data), size = floor(nrow(data) * fraction))
  
  data[rows_to_sample, , drop = FALSE]
}

get_cutoff <- function(max, min, mid, scope, output_list = TRUE){
  
  lower_sampling_bound <- mid_per_iter - (0.5 * scope * (max_per_iter - min_per_iter))
  upper_sampling_bound <- mid_per_iter + (0.5 * scope * (max_per_iter - min_per_iter))
  
  if (output_list) {
    cutoff_points <- purrr::pmap(
      list(rep(1, length(lower_sampling_bound)),
        lower_sampling_bound,
        upper_sampling_bound
      ),
      runif
    )
  } else {
    cutoff_points <- matrix(runif(
      length(lower_sampling_bound),
      lower_sampling_bound,
      upper_sampling_bound
    ),
    nrow = 1)
  }

  cutoff_points

}

get_mass <- function(cutoff_points,
                     projections,
                     side = c("left", "right"),
                     subsample_size) {
  
  points_in_halfspace <- switch(
    side,
    "left" = purrr::map2(.x = projections_list, .y = cutoff_points, `<`),
    "right" = purrr::map2(.x = projections_list, .y = cutoff_points, `>=`)
  )
  
  mass <- purrr::map2(
    lapply(points_in_halfspace, sum),
    subsample_size,
    `/`
  )
  
  mass
  
}
