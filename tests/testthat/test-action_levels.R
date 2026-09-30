test_that("action_levels classifies waist percentiles correctly", {

  df <- data.frame(
    waist_percentile = c(0.85, 0.96)
  )

  result <- action_levels(df)

  expect_named(result, "adiposity.level")

  expect_equal(
    as.character(result$adiposity.level),
    c("none", "action")
  )
})


test_that("action_levels calculates multiple domains correctly", {

  df <- data.frame(
    hdl_percentile   = c(0.1, 0.5),
    homa_percentile  = c(0.4, 0.9),
    trg_percentile   = c(0.6, 0.5),
    waist_percentile = c(0.9, 0.99),
    sbp_percentile   = c(0.8, 0.01)
  )

  result <- action_levels(df)

  expect_equal(
    as.character(result$adiposity.level),
    c("none", "action")
  )

  expect_equal(
    as.character(result$blood_pressure.level),
    c("none", "none")
  )

  expect_equal(
    as.character(result$blood_lipids.level),
    c("none", "none")
  )

  expect_equal(
    as.character(result$blood_glu_insu.level),
    c("none", "none")
  )

  expect_equal(
    as.character(result$overall.level),
    c("none", "none")
  )
})


test_that("action_levels calculates sex-specific CRP levels correctly", {

  df <- data.frame(
    hdl_percentile   = c(0.1, 0.5),
    homa_percentile  = c(0.4, 0.9),
    trg_percentile   = c(0.6, 0.5),
    crp_percentile   = c(0.95, 0.9),
    waist_percentile = c(0.9, 0.99),
    sbp_percentile   = c(0.8, 0.01)
  )

  result <- action_levels(
    df,
    sex = c("m", "m")
  )

  expect_equal(
    as.character(result$crp.level),
    c("Elevated", "Not elevated")
  )

  expect_equal(
    as.character(result$adiposity.level),
    c("none", "action")
  )

  expect_equal(
    as.character(result$blood_pressure.level),
    c("none", "none")
  )

  expect_equal(
    as.character(result$blood_lipids.level),
    c("none", "none")
  )

  expect_equal(
    as.character(result$blood_glu_insu.level),
    c("none", "none")
  )

  expect_equal(
    as.character(result$overall.level),
    c("none", "none")
  )
})
