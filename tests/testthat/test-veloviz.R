test_that("veloviz makes a graph", {
  data(vel)
  expect_length(buildVeloviz(
    curr = vel$current, proj = vel$projected,
    normalize.depth = TRUE,
    use.ods.genes = FALSE,
    alpha = 0.05,
    pca = TRUE,
    nPCs = 3,
    center = TRUE,
    scale = TRUE,
    k = 10,
    similarity.threshold = -1,
    distance.weight = 1,
    distance.threshold = 1,
    weighted = TRUE,
    verbose = FALSE
  ),
  3)

  expect_equal(class(buildVeloviz(
    curr = vel$current, proj = vel$projected,
    normalize.depth = TRUE,
    use.ods.genes = FALSE,
    alpha = 0.05,
    pca = TRUE,
    nPCs = 3,
    center = TRUE,
    scale = TRUE,
    k = 10,
    similarity.threshold = -1,
    distance.weight = 1,
    distance.threshold = 1,
    weighted = TRUE,
    verbose = FALSE
  )[[1]]),
  "igraph")

  expect_equal(ncol(buildVeloviz(
    curr = vel$current, proj = vel$projected,
    normalize.depth = TRUE,
    use.ods.genes = FALSE,
    alpha = 0.05,
    pca = TRUE,
    nPCs = 3,
    center = TRUE,
    scale = TRUE,
    k = 10,
    similarity.threshold = -1,
    distance.weight = 1,
    distance.threshold = 1,
    weighted = TRUE,
    verbose = FALSE
  )[[2]]),
  2)

})
