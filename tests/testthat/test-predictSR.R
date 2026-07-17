context("test-predictSR")

fitSpider <- function() {
    data(eSpider, package = "gllvm")
    y <- as.matrix(eSpider$abund[eSpider$nonNA,])
    X <- scale(eSpider$X[eSpider$nonNA,])
    fit <- gllvm(y = y, X = X, formula = ~ConWate,
                           family = poisson(), seed = 999)
    fit
  }

test_that("predictSR returns a valid predictSR.gllvm object", {
  fit <- fitSpider()
  predSR <- predictSR(fit, se.fit = FALSE)

  expect_s3_class(predSR, "predictSR.gllvm")
  expect_true(all(c("predicted", "expected", "spp") %in% names(predSR)))

  # pmf: one column per possible richness value 0, 1, ..., p
  pmf <- predSR$predicted$fit
  expect_equal(dim(pmf), c(nrow(fit$y), ncol(fit$y)+1))
  expect_equal(colnames(pmf), paste0("SR_", 0:ncol(fit$y)))
  # each row is a proper probability mass function
  expect_true(all(pmf >= 0))
  expect_true(all(abs(rowSums(pmf) - 1) < 1e-6))

  # expected richness: one value per site, within [0, p]
  eSR <- predSR$expected$fit
  expect_equal(length(eSR), nrow(fit$y))
  expect_true(all(is.finite(eSR)))
  expect_true(all(eSR >= 0 & eSR <= ncol(fit$y)))
})

test_that("predictSR returns a simulated CI on expected richness", {
  fit <- fitSpider()

  predSR <- predictSR(fit, se.fit = 50, seed = 1)

  expect_true(all(c("lower", "upper") %in% names(predSR$expected)))
  expect_equal(length(predSR$expected$lower), nrow(fit$y))
  expect_equal(length(predSR$expected$upper), nrow(fit$y))
  expect_true(all(is.finite(predSR$expected$lower)))
  expect_true(all(is.finite(predSR$expected$upper)))
  expect_true(all(predSR$expected$lower <= predSR$expected$upper))
})

test_that("predictSR respects the spp argument", {
  fit <- fitSpider()

  predSR <- predictSR(fit, spp = 1:4, se.fit = FALSE)

  expect_equal(predSR$spp, 1:4)
  # richness of 4 species can only be 0, 1, ..., 4
  expect_equal(ncol(predSR$predicted$fit), 5)
  expect_true(all(abs(rowSums(predSR$predicted$fit) - 1) < 1e-6))
})

test_that("predictSR predicts for new sites with level = 0", {
  fit <- fitSpider()

  newX <- data.frame(ConWate = seq(min(fit$X[, "ConWate"]),
                                   max(fit$X[, "ConWate"]), length.out = 20))
  predSR <- predictSR(fit, newX = newX, level = 0, se.fit = FALSE)

  expect_equal(dim(predSR$predicted$fit), c(20, ncol(fit$y) + 1))
  expect_true(all(abs(rowSums(predSR$predicted$fit) - 1) < 1e-6))
})

test_that("residuals.predictSR.gllvm returns fitted values and residuals", {
  fit <- fitSpider()

  predSR <- predictSR(fit, se.fit = FALSE)
  res <- residuals(predSR, model = fit)

  expect_true(all(c("fitted", "residuals") %in% names(res)))
  expect_equal(length(res$fitted), nrow(fit$y))
  expect_equal(length(res$residuals), nrow(fit$y))
  expect_true(all(is.finite(res$residuals)))
})

test_that("plot.predictSR.gllvm produces residual diagnostics", {
  fit <- fitSpider()

  predSR <- predictSR(fit, se.fit = FALSE)

  # send graphics to a temporary device so no window/file pollutes the run
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp)
  on.exit({ dev.off(); unlink(tmp) }, add = TRUE)

  expect_error(plot(predSR, object = fit), NA)              # both panels
  expect_error(plot(predSR, object = fit, which = 1), NA)   # residuals vs fitted
  expect_error(plot(predSR, object = fit, which = 2), NA)   # QQ-plot
})
