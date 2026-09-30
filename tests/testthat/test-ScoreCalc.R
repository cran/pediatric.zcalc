test_that("ScoreCalc calculates percentiles and action levels correctly", {

  df <- data.frame(
    sex = c("f", "m"),
    age = c(8, 15),
    height = c(120, 125),
    waist = c(55, 60),
    homa = c(1.2, 1.4),
    sbp = c(100, 105),
    dbp = c(65, 70),
    crp = c(6, 2),

    # IDEFICS lipid reference values are expressed in mg/dL.
    # Use a high triglyceride value for the first child so that
    # the blood-lipid action level is exercised deliberately.
    trg = c(150, 50),
    hdl = c(52, 52)
  )

  result <- ScoreCalc(
    df,
    return_values = c("percentile", "cutoff.levels")
  )

  expect_equal(nrow(result), 2)

  expect_true(all(c(
    "waist_percentile",
    "height_percentile",
    "homa_percentile",
    "sbp_percentile",
    "dbp_percentile",
    "crp_percentile",
    "trg_percentile",
    "hdl_percentile"
  ) %in% names(result)))

  # CRP supports the older age.
  expect_false(is.na(result$crp_percentile[2]))

  # Age 15 is outside the IDEFICS reference range for these variables.
  expect_true(is.na(result$waist_percentile[2]))
  expect_true(is.na(result$height_percentile[2]))
  expect_true(is.na(result$homa_percentile[2]))
  expect_true(is.na(result$sbp_percentile[2]))
  expect_true(is.na(result$dbp_percentile[2]))
  expect_true(is.na(result$trg_percentile[2]))
  expect_true(is.na(result$hdl_percentile[2]))

  # All available percentiles must be valid probabilities.
  percentile_cols <- c(
    "waist_percentile",
    "height_percentile",
    "homa_percentile",
    "sbp_percentile",
    "dbp_percentile",
    "crp_percentile",
    "trg_percentile",
    "hdl_percentile"
  )

  for (column in percentile_cols) {
    values <- result[[column]]
    values <- values[!is.na(values)]

    expect_true(all(is.finite(values)))
    expect_true(all(values >= 0 & values <= 1))
  }

  # Keep exact regression checks for variables whose inputs have not changed.
  expect_equal(
    result$height_percentile,
    c(0.0451617720789443, NA_real_),
    tolerance = 1e-8
  )

  expect_equal(
    result$homa_percentile,
    c(0.67496379032534, NA_real_),
    tolerance = 1e-8
  )

  expect_equal(
    result$waist_percentile,
    c(0.493334606610262, NA_real_),
    tolerance = 1e-8
  )

  expect_equal(
    result$sbp_percentile,
    c(0.523532114847159, NA_real_),
    tolerance = 1e-8
  )

  expect_equal(
    result$dbp_percentile,
    c(0.616152122915867, NA_real_),
    tolerance = 1e-8
  )

  expect_equal(
    result$crp_percentile,
    c(0.988553314401349, 0.906628776100773),
    tolerance = 1e-8
  )

  expect_equal(
    as.character(result$adiposity.level),
    c("none", NA)
  )

  expect_equal(
    as.character(result$crp.level),
    c("Elevated", "Not elevated")
  )

  expect_equal(
    as.character(result$blood_pressure.level),
    c("none", NA)
  )

  # The deliberately high triglyceride value should trigger action.
  expect_equal(
    as.character(result$blood_lipids.level),
    c("action", NA)
  )

  expect_equal(
    as.character(result$blood_glu_insu.level),
    c("none", NA)
  )

  expect_equal(
    as.character(result$overall.level),
    c("none", NA)
  )

  # Action-level outputs should remain ordered factors.
  expect_true(is.ordered(result$adiposity.level))
  expect_true(is.ordered(result$crp.level))
  expect_true(is.ordered(result$blood_pressure.level))
  expect_true(is.ordered(result$blood_lipids.level))
  expect_true(is.ordered(result$blood_glu_insu.level))
  expect_true(is.ordered(result$overall.level))
})


test_that("ScoreCalc reproduces published IDEFICS MetS examples", {

  # Validation data are the four example children reported in Appendix A of:
  #
  # Ahrens W, Moreno LA, Marild S, et al. (2014).
  # "Metabolic syndrome in young children: definitions and results
  # of the IDEFICS study."
  # International Journal of Obesity 38(Suppl 2): S4-S14.
  # doi:10.1038/ijo.2014.130
  #
  # Appendix A reports the measured values, component z-scores and
  # resulting MetS scores. Heights for the same four example children
  # are provided in Example_Data.csv of the original BIPS
  # IDEFICS Score Calculator.
  #
  # TRG and HDL are supplied in mg/dL, matching the IDEFICS
  # reference parameterisation.

  df <- data.frame(
    sex = c("m", "f", "f", "m"),
    age = c(6.5, 5.8, 5.2, 5.5),
    height = c(126, 118, 112, 119),
    waist = c(59.5, 60.4, 53.0, 57.5),
    homa = c(0.65, 1.59, 1.00, 0.99),
    sbp = c(109.0, 92.0, 94.0, 99.0),
    dbp = c(67.5, 52.5, 58.0, 66.0),
    trg = c(48, 50, 45, 78),
    hdl = c(38, 52, 61, 45)
  )

  result <- ScoreCalc(
    df,
    return_values = c("z.score", "MetS")
  )

  expect_equal(nrow(result), 4)

  # Published z-scores from Appendix A.
  #
  # The publication reports z-scores rounded to two decimal places,
  # so a tolerance appropriate for the published precision is used.

  expect_equal(
    result$waist_z.score,
    c(1.66, 2.20, 0.58, 1.64),
    tolerance = 0.02
  )

  expect_equal(
    result$homa_z.score,
    c(-0.15, 1.39, 0.68, 0.71),
    tolerance = 0.02
  )

  expect_equal(
    result$sbp_z.score,
    c(0.85, -0.96, -0.51, -0.13),
    tolerance = 0.02
  )

  expect_equal(
    result$dbp_z.score,
    c(0.72, -1.77, -0.81, 0.59),
    tolerance = 0.02
  )

  expect_equal(
    result$trg_z.score,
    c(0.29, 0.17, -0.73, 1.15),
    tolerance = 0.02
  )

  # Subject 3 has TRG = 45 mg/dL, which corresponds to the
  # triglyceride detection limit reported for the IDEFICS study.
  # Triglyceride measurements below this limit were handled specially
  # when deriving the reference distribution.
  expect_true(is.finite(result$trg_z.score[3]))

  expect_equal(
    result$hdl_z.score,
    c(-1.22, 0.10, 0.88, -0.54),
    tolerance = 0.02
  )

  # Published MetS scores from Appendix A.
  #
  # For subjects 1, 2 and 4, the reported component z-scores reproduce
  # the published MetS score using the stated formula.
  #
  # Subject 3 is not included in this comparison because Appendix A is
  # internally inconsistent for this row: the reported component z-scores
  # yield a MetS score of approximately -0.205, whereas the table reports
  # -0.07.
  expect_equal(
    result$MetS[c(1, 2, 4)],
    c(3.05, 2.26, 3.43),
    tolerance = 0.02
  )

  # The standardized MetS z-score is not reported in Appendix A,
  # but it should be finite for all four published examples.
  expect_true(all(is.finite(result$MetS_z.score)))

  # Internal consistency: MetS must follow the published formula:
  #
  # z_waist + z_HOMA +
  #   (z_SBP + z_DBP) / 2 +
  #   (z_TRG - z_HDL) / 2

  expected_MetS <-
    result$waist_z.score +
    result$homa_z.score +
    0.5 * (
      result$sbp_z.score +
        result$dbp_z.score +
        result$trg_z.score -
        result$hdl_z.score
    )

  expect_equal(
    result$MetS,
    expected_MetS,
    tolerance = 1e-10
  )
})
