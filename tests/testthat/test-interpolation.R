test_that("DBP boys sigma interpolation is stable", {

  # approx_param_functions <- get(
  #   "approx_param_functions",
  #   envir = asNamespace("IDEFICS.scalc")
  # )

  result <- approx_param_functions$dbp_boys_sigma(
    age = 5,
    height = 120
  )

  expect_equal(
    result,
    0.0975872601856281,
    tolerance = 1e-10
  )
})


test_that("DBP mu interpolation is stable across sex and height", {

  result <- approx_param_functions$dbp_mu(
    sex = c("f", "m", "f", "m"),
    age = c(5, 5, 5, 5),
    height = c(120, 120, 121, 121)
  )

  expect_equal(
    result,
    c(63.5974441138056, 62.2938403110434, 63.6883607297003, 62.3727561518278
    ),
    tolerance = 1e-10
  )
})


test_that("BMI sigma interpolation is stable across sex", {

  result <- approx_param_functions$bmi_sigma(
    sex = c("f", "m"),
    age = c(5, 5)
  )

  expect_equal(
    result,
    c(0.0901899999999977, 0.0835299999999977),
    tolerance = 1e-10
  )
})
