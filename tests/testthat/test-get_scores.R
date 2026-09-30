test_that("waist z-scores are calculated correctly", {

  result <- get_scores(
    variable = "waist",
    sex = "f",
    age = c(5, 6),
    values = c(50, 52),
    return_values = "z.score"
  )

  expect_named(result, "z.score")
  expect_length(result$z.score, 2)

  expect_equal(
    result$z.score,
    c(-0.360472866517135, 0.345562615961339),
    tolerance = 1e-8
  )
})


test_that("CRP scores are calculated correctly", {

  result <- get_scores(
    variable = "crp",
    sex = c("f", "f"),
    age = c(3, 12),
    values = c(4, 0.5)
  )

  expect_named(result, c("percentile", "z.score"))
  expect_length(result$percentile, 2)
  expect_length(result$z.score, 2)

  expect_equal(
    result$percentile,
    c(0.956611555583334, 0.700291111444043),
    tolerance = 1e-8
  )

  expect_equal(
    result$z.score,
    c(1.71265032592801, 0.525237963355671),
    tolerance = 1e-8
  )
})


test_that("DBP scores are calculated correctly with height dependence", {

  result <- get_scores(
    variable = "dbp",
    sex = c("f", "m"),
    age = c(5, 5),
    height = c(120, 110),
    values = c(70, 60)
  )

  expect_named(result, c("percentile", "z.score"))
  expect_length(result$percentile, 2)
  expect_length(result$z.score, 2)

  expect_equal(
    result$percentile,
    c(0.840000948105481, 0.401475695007478),
    tolerance = 1e-8
  )

  expect_equal(
    result$z.score,
    c(0.994461779886976, -0.249529286280741),
    tolerance = 1e-8
  )
})


