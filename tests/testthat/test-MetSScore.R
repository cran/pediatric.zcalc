test_that("MetSScore calculates the metabolic syndrome score correctly", {

  df <- data.frame(
    waist_z.score = c(1.2, 0.5),
    homa_z.score  = c(0.8, 0.6),
    sbp_z.score   = c(0.7, 0.3),
    dbp_z.score   = c(0.6, 0.4),
    trg_z.score   = c(1.0, 0.9),
    hdl_z.score   = c(-0.5, -0.2)
  )

  result <- MetSScore(df)

  expect_equal(
    result,
    c(3.4, # 1.2 + 0.8 + 0.5 * (0.7 + 0.6 + 1.0 - (-0.5)) = 3.4
      2.0  # 0.5 + 0.6 + 0.5 * (0.3 + 0.4 + 0.9 - (-0.2)) = 2.0
      ),
    tolerance = 1e-10
  )
})
