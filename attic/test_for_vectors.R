# Test to check the projection

get_normalized_vector(2)

test <- matrix(
  c(2,2,-2,2,2,-2,-2,-2),
  ncol = 2,
  nrow = 4,
  byrow = T
)
normal <- get_normalized_vector(2)
intercept <- -1.2
normal_1 <- normal + intercept
test %*% normal
(test-intercept) %*% normal

plot(test)

arrows(0,0,normal[1,],normal[2,], col = "red")
abline(coef = c(0, normal[1]/(-normal[2])), col = "red")

arrows(normal[1,],normal[2,], normal_1[1,], normal_1[2,], col = "blue")
abline(coef = c(intercept/(-normal[2]), normal[1]/(-normal[2])), col = "blue")


arrows(0,0,normal_1[1,],normal_1[2,], col = "blue")

arrows(-0.5,-1,normal_1[1,],normal_1[2,], col = "blue")

dev.off()
