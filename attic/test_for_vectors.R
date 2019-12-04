# Test to check the projection

get_normalized_vector(2)

n <- 500

data <- matrix(
  c(rnorm(n/2, mean = sqrt(1:n), sd = log(1:n)),rexp(n/2,seq(from = 1, to = 2, length.out = 1000))),
  ncol = 2,
  nrow = n/2,
  byrow = F
)

plot(data)
debugonce(train_depth)
halfspaces <- train_depth(data_fig3, n_halfspace = 1e4, subsample = 1, scope = 1, plot_yes = TRUE, seed = 4163)

plot(data)
points(data[which.max(evaluate_depth(data = data, halfspaces = halfspaces)), ], col = "red")
points(data[which.min(evaluate_depth(data = data, halfspaces = halfspaces)), ], col = "red")


# 
# normal <- get_normalized_vector(2)
# intercept <- -1.2
# normal_1 <- normal + intercept
# test %*% normal
# (test-intercept) %*% normal
# 
# plot(test)
# 
# arrows(0,0,normal[1,],normal[2,], col = "red")
# abline(coef = c(0, normal[1]/(-normal[2])), col = "red")
# 
# arrows(normal[1,],normal[2,], normal_1[1,], normal_1[2,], col = "blue")
# abline(coef = c(intercept/(-normal[2]), normal[1]/(-normal[2])), col = "blue")
# 
# 
# arrows(0,0,normal_1[1,],normal_1[2,], col = "blue")
# 
# arrows(-0.5,-1,normal_1[1,],normal_1[2,], col = "blue")
# 
# dev.off()
# 
# halfspaces <- train_depth(test, n_halfspace = 10, subsample = 1, scope = 1)
