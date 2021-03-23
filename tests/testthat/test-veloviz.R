test_that("veloviz makes a graph", {
  expect_length(buildVeloviz(
    curr = pancreas$vel$current, proj = pancreas$vel$projected,
    normalize.depth = TRUE,
    use.ods.genes = TRUE,
    alpha = 0.05,
    pca = TRUE,
    nPCs = 20,
    center = TRUE,
    scale = TRUE,
    k = 5,
    similarity.threshold = 0.25,
    distance.weight = 1,
    distance.threshold = 0.5,
    weighted = TRUE,
    seed = 0,
    verbose = FALSE
  ),
  3)
  
  expect_equal(class(buildVeloviz(
    curr = pancreas$vel$current, proj = pancreas$vel$projected,
    normalize.depth = TRUE,
    use.ods.genes = TRUE,
    alpha = 0.05,
    pca = TRUE,
    nPCs = 20,
    center = TRUE,
    scale = TRUE,
    k = 5,
    similarity.threshold = 0.25,
    distance.weight = 1,
    distance.threshold = 0.5,
    weighted = TRUE,
    seed = 0,
    verbose = FALSE
  )[[1]]),
  "igraph")
  
  expect_equal(ncol(buildVeloviz(
    curr = pancreas$vel$current, proj = pancreas$vel$projected,
    normalize.depth = TRUE,
    use.ods.genes = TRUE,
    alpha = 0.05,
    pca = TRUE,
    nPCs = 20,
    center = TRUE,
    scale = TRUE,
    k = 5,
    similarity.threshold = 0.25,
    distance.weight = 1,
    distance.threshold = 0.5,
    weighted = TRUE,
    seed = 0,
    verbose = FALSE
  )[[2]]),
  2)
  
})