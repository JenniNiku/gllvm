context("test-glmmVA")

# Shared test data: spider abundance in long format
make_spider_long <- function() {
  data(eSpider)
  abund <- eSpider$abund[eSpider$nonNA, ]
  X     <- eSpider$X[eSpider$nonNA, ]
  p <- ncol(abund); n <- nrow(abund)
  data.frame(
    abund    = c(abund),
    Species  = factor(rep(colnames(abund), each = n)),
    id       = factor(rep(seq_len(n), times = p)),
    BareSand = rep(X[, "BareSand"], times = p),
    ConWate  = rep(X[, "ConWate"],  times = p)
  )
}

test_that("glmmVA: fixed effects only returns glmmVA object", {
  dat <- make_spider_long()
  m <- glmmVA(abund ~ BareSand + Species, data = dat, family = "poisson",
              sd.errors = FALSE)
  expect_s3_class(m, "glmmVA")
  expect_s3_class(m, "gllvm")
  expect_null(m$params$row.params.random)
})

test_that("glmmVA: random intercept per group returns random effects", {
  dat <- make_spider_long()
  m <- glmmVA(abund ~ BareSand + Species + (1|id), data = dat,
              family = "poisson", sd.errors = FALSE)
  expect_s3_class(m, "glmmVA")
  expect_false(is.null(m$params$row.params.random))
  # one random intercept per level (minus reference); may vary by implementation
  expect_gt(length(m$params$row.params.random), 0L)
})

test_that("glmmVA: continuous random slope gives 1x1 sigmaijr", {
  dat <- make_spider_long()
  m <- glmmVA(abund ~ Species + (0 + BareSand|id), data = dat,
              family = "poisson", sd.errors = FALSE)
  trmsize <- m$TMBfn$env$data$trmsize
  expect_equal(unname(trmsize[1, 1]), 1L)
  expect_equal(m$corP$cstruc, "diag")
})

test_that("glmmVA: factor LHS gives nc x nc sigmaijr (bug fix)", {
  dat   <- make_spider_long()
  nspec <- nlevels(dat$Species)
  m <- glmmVA(abund ~ BareSand + Species + (0 + Species|id), data = dat,
              family = "poisson", sd.errors = FALSE)
  trmsize <- m$TMBfn$env$data$trmsize
  expect_equal(unname(trmsize[1, 1]), nspec,
               info = "trmsize should equal nlevels(Species), not 1")
  expect_equal(m$corP$cstruc, "ustruc")
  expect_equal(dim(m$params$sigmaijr), c(nspec, nspec))
})

test_that("glmmVA: correlated slopes give correct sigmaijr dimensions", {
  dat <- make_spider_long()
  m <- glmmVA(abund ~ Species + (0 + BareSand + ConWate|id), data = dat,
              family = "poisson", sd.errors = FALSE)
  trmsize <- m$TMBfn$env$data$trmsize
  expect_equal(unname(trmsize[1, 1]), 2L)
  expect_equal(m$corP$cstruc, "ustruc")
  expect_equal(dim(m$params$sigmaijr), c(2L, 2L))
})

test_that("row.eff formula: continuous LHS stays diag", {
  data(eSpider)
  y  <- eSpider$abund[eSpider$nonNA, ]
  X  <- as.data.frame(scale(eSpider$X[eSpider$nonNA, ]))
  sd <- cbind(data.frame(id = factor(seq_len(nrow(y)))), X)
  m <- gllvm(y, studyDesign = sd,
             row.eff = ~(0 + BareSand|id),
             family = "negative.binomial", num.lv = 0, sd.errors = FALSE)
  expect_equal(m$corP$cstruc, "diag")
  expect_equal(unname(m$TMBfn$env$data$trmsize[1, 1]), 1L)
})

test_that("formula: random slope gives correct sigmaB", {
  data(eSpider)
  y <- eSpider$abund[eSpider$nonNA, ]
  X <- as.data.frame(scale(eSpider$X[eSpider$nonNA, ]))
  m <- gllvm(y, X = X,
             formula = ~(0 + BareSand|1),
             family = "negative.binomial", num.lv = 0,
             beta0com = TRUE, sd.errors = FALSE)
  expect_false(is.null(m$params$sigmaB))
  expect_equal(dim(m$params$sigmaB), c(1L, 1L))
})

test_that("ranef.glmmVA: returns object without condVar by default", {
  dat <- make_spider_long()
  m <- glmmVA(abund ~ BareSand + Species + (1|id), data = dat,
              family = "poisson", sd.errors = FALSE)
  re <- ranef(m)
  expect_null(attr(re, "condVar"))
  expect_true(length(re) > 0L)
})

test_that("ranef.glmmVA: condVar = TRUE attaches condVar attribute", {
  dat <- make_spider_long()
  m <- glmmVA(abund ~ BareSand + Species + (1|id), data = dat,
              family = "poisson", sd.errors = TRUE)
  re <- ranef(m, condVar = TRUE)
  expect_false(is.null(attr(re, "condVar")))
})

# ---- formula random effects (species-specific slopes) ----------------------

test_that("formula: uncorrelated random slopes give diagonal sigmaB", {
  data(eSpider)
  y <- eSpider$abund[eSpider$nonNA, ]
  X <- as.data.frame(scale(eSpider$X[eSpider$nonNA, ]))
  # 0+ excludes the intercept so sigmaB is 2x2 (one row per covariate)
  m <- gllvm(y, X = X,
             formula = ~diag(0 + BareSand + ConWate|1),
             family = "negative.binomial", num.lv = 0,
             beta0com = TRUE, sd.errors = FALSE)
  expect_false(is.null(m$params$sigmaB))
  expect_equal(dim(m$params$sigmaB), c(2L, 2L))
  expect_equal(unname(m$params$sigmaB[1, 2]), 0)
})

test_that("formula: correlated random slopes give non-diagonal sigmaB", {
  data(eSpider)
  y <- eSpider$abund[eSpider$nonNA, ]
  X <- as.data.frame(scale(eSpider$X[eSpider$nonNA, ]))
  # 0+ excludes the intercept so sigmaB is 2x2
  m <- gllvm(y, X = X,
             formula = ~(0 + BareSand + ConWate|1),
             family = "negative.binomial", num.lv = 0,
             beta0com = TRUE, sd.errors = FALSE)
  expect_false(is.null(m$params$sigmaB))
  expect_equal(dim(m$params$sigmaB), c(2L, 2L))
})

test_that("formula: random intercept (1|1) gives scalar sigmaB", {
  data(eSpider)
  y <- eSpider$abund[eSpider$nonNA, ]
  m <- gllvm(y, formula = ~(1|1),
             family = "negative.binomial", num.lv = 0,
             beta0com = TRUE, sd.errors = FALSE)
  expect_false(is.null(m$params$sigmaB))
  expect_equal(dim(m$params$sigmaB), c(1L, 1L))
})

# ---- lv.formula random effects (effects in the ordination) -----------------

test_that("lv.formula: random slope in ordination gives sigmaLvXcoef", {
  data(eSpider)
  y <- eSpider$abund[eSpider$nonNA, ]
  X <- as.data.frame(scale(eSpider$X[eSpider$nonNA, ]))
  # lv.formula with random slopes requires num.lv.c > 0 and randomB
  m <- gllvm(y, X = X, lv.formula = ~(0 + BareSand|1),
             family = "negative.binomial", num.lv.c = 2,
             randomB = "LV", sd.errors = FALSE,
             control.start = list(starting.val = "zero"))
  expect_false(is.null(m$params$sigmaLvXcoef))
  expect_equal(m$randomB, "LV")
})

test_that("lv.formula: fixed effects in ordination give LvXcoef", {
  data(eSpider)
  y <- eSpider$abund[eSpider$nonNA, ]
  X <- as.data.frame(scale(eSpider$X[eSpider$nonNA, ]))
  # num.RR cannot exceed predictors in lv.formula; use 1
  m <- gllvm(y, X = X, lv.formula = ~BareSand,
             family = "negative.binomial", num.RR = 1,
             sd.errors = FALSE)
  expect_false(is.null(m$params$LvXcoef))
  expect_equal(m$num.RR, 1L)
})

# ---- issue #241: diag(nc>1|grp) fix -------------------------------------------

make_spider_long_dummy <- function() {
  data(eSpider)
  abund <- eSpider$abund[eSpider$nonNA, ]
  X     <- as.data.frame(scale(eSpider$X[eSpider$nonNA, ]))
  p <- ncol(abund); n <- nrow(abund)
  data.frame(
    abund    = c(abund),
    Species  = factor(rep(colnames(abund), each = n)),
    dummy    = factor(1),
    BareSand = rep(X[, "BareSand"], times = p),
    ConWate  = rep(X[, "ConWate"],  times = p)
  )
}

test_that("glmmVA: diag(a+b|grp) gives same logLik as (0+a|grp)+(0+b|grp) (issue #241)", {
  dat <- make_spider_long_dummy()
  ctrl <- list(optimizer = "nlminb", max.iter = 2000, maxit = 2000)
  m1 <- glmmVA(abund ~ Species + diag(0 + BareSand + ConWate|dummy),
               family = "negative.binomial", data = dat,
               sd.errors = FALSE, control = ctrl)
  m2 <- glmmVA(abund ~ Species + (0 + BareSand|dummy) + (0 + ConWate|dummy),
               family = "negative.binomial", data = dat,
               sd.errors = FALSE, control = ctrl)
  expect_equal(logLik(m1), logLik(m2), tolerance = 1e-3,
               info = "diag(a+b|grp) and (0+a|grp)+(0+b|grp) must give identical log-likelihoods")
  expect_equal(unname(m1$params$sigma), unname(m2$params$sigma), tolerance = 1e-3)
})

test_that("glmmVA: diag(nc>1|grp) uses nc lg_Ar entries (not nl)", {
  dat <- make_spider_long_dummy()
  m <- glmmVA(abund ~ Species + diag(0 + BareSand + ConWate|dummy),
              family = "negative.binomial", data = dat, sd.errors = FALSE)
  trmsize <- m$TMBfn$env$data$trmsize
  nc <- unname(trmsize[1, 1]); nl <- unname(trmsize[2, 1])
  expect_equal(length(m$TMBfn$env$parameters$lg_Ar), nc * nl,
               info = "lg_Ar must have nc*nl entries for diag(nc>1|grp)")
  # Two sigma values (one per covariate), not one
  expect_equal(length(m$params$sigma), 2L)
})

# ---- issue #250: 0+ intercept suppression fix ---------------------------------

test_that("glmmVA: 0+ suppresses global intercept (beta0 == 0)", {
  dat <- make_spider_long()
  m <- glmmVA(abund ~ 0 + Species + (0 + BareSand|id), data = dat,
              family = "poisson", sd.errors = FALSE)
  expect_equal(unname(m$params$beta0), rep(0, length(m$params$beta0)),
               tolerance = 1e-10, info = "beta0 must be zero when formula starts with 0+")
})

test_that("glmmVA: 0+ model random effects are non-trivial", {
  dat <- make_spider_long()
  m0 <- glmmVA(abund ~ 0 + Species + (0 + BareSand|id), data = dat,
               family = "poisson", sd.errors = FALSE)
  expect_gt(length(m0$params$row.params.random), 0L)
  expect_false(all(m0$params$row.params.random == 0))
})

# ---- structured row effects: corExp and propto --------------------------------

test_that("row.eff corExp fits and produces Scale/range sigma parameters", {
  data(eSpider)
  y  <- eSpider$abund[eSpider$nonNA, ]
  n  <- nrow(y)
  set.seed(42)
  coords <- matrix(sort(runif(n)), ncol = 1)
  sd <- data.frame(id = factor(seq_len(n)))
  m <- gllvm(y, studyDesign = sd,
             row.eff = ~corExp(1|id),
             dist = list(coords),
             family = "negative.binomial", num.lv = 0, sd.errors = FALSE)
  expect_equal(m$corP$cstruc, "corExp")
  expect_true(any(grepl("\\.Scale$", names(m$params$sigma))),
              info = "corExp must produce a sigma entry named '<term>.Scale'")
  expect_equal(length(m$params$sigma), 2L,
               info = "corExp produces exactly 2 sigma values (scale + range)")
})

# ---- mixed-response (per-observation family vector) ---------------------------

# Two response groups from the spider data: site-level random intercept,
# first column Poisson, second column negative binomial.
make_spider_mixed <- function() {
  data(eSpider)
  abund <- eSpider$abund[eSpider$nonNA, ]
  X     <- eSpider$X[eSpider$nonNA, ]
  n <- nrow(abund)
  # Collapse species into 2 groups to keep the model small
  grp1 <- rowSums(abund[, 1:6])
  grp2 <- rowSums(abund[, 7:12])
  data.frame(
    y        = c(grp1, grp2),
    site     = factor(c(seq_len(n), seq_len(n))),
    family   = c(rep("poisson", n), rep("negative.binomial", n)),
    BareSand = c(X[, "BareSand"], X[, "BareSand"])
  )
}

test_that("glmmVA mixed: auto-infers response.group from per-obs family vector", {
  dat <- make_spider_mixed()
  m <- glmmVA(y ~ (1|site), family = dat$family, data = dat, sd.errors = FALSE)
  expect_s3_class(m, "glmmVA")
  # y stored as wide matrix with NAs
  expect_true(is.matrix(m$y))
  expect_equal(ncol(m$y), 2L)
  expect_true(any(is.na(m$y)))
})

test_that("glmmVA mixed: nobs counts non-NA cells, not all matrix cells", {
  dat <- make_spider_mixed()
  n   <- nlevels(dat$site)
  m   <- glmmVA(y ~ (1|site), family = dat$family, data = dat, sd.errors = FALSE)
  # Total observations = 2*n (one per long-format row), not 4*n (prod of matrix dims)
  expect_equal(nobs(m), nrow(dat))
  expect_equal(attributes(logLik(m))$nobs, nrow(dat))
})

test_that("glmmVA mixed: beta0com collapses intercepts to 1 shared value", {
  dat <- make_spider_mixed()
  m   <- glmmVA(y ~ (1|site), family = dat$family, data = dat, sd.errors = FALSE)
  expect_true(m$beta0com)
  # Both columns share the same intercept
  expect_equal(length(unique(m$params$beta0)), 1L)
})

test_that("glmmVA mixed: 0+ suppresses intercept without error", {
  dat <- make_spider_mixed()
  m   <- glmmVA(y ~ 0 + (1|site), family = dat$family, data = dat, sd.errors = FALSE)
  expect_equal(unname(m$params$beta0), rep(0, length(m$params$beta0)),
               tolerance = 1e-10,
               info = "beta0 must be zero when 0+ is used in mixed mode")
})

test_that("predict.glmmVA mixed: returns long-format vector, not wide matrix", {
  dat <- make_spider_mixed()
  m   <- glmmVA(y ~ (1|site), family = dat$family, data = dat, sd.errors = FALSE)
  preds <- predict(m, type = "response")
  expect_false(is.matrix(preds), info = "predict() must collapse wide matrix to vector")
  expect_equal(length(preds), nrow(dat))
  expect_true(all(is.finite(preds)))
})

test_that("row.eff propto fits and produces sigma parameter", {
  data(eSpider)
  y  <- eSpider$abund[eSpider$nonNA, ]
  n  <- nrow(y)
  set.seed(42)
  M_propto_test <<- crossprod(matrix(rnorm(n * n), n, n))
  on.exit(suppressWarnings(rm("M_propto_test", envir = globalenv())), add = TRUE)
  sd <- data.frame(id = factor(seq_len(n)))
  m <- gllvm(y, studyDesign = sd,
             row.eff = ~propto(1|id, M_propto_test),
             family = "negative.binomial", num.lv = 0, sd.errors = FALSE)
  expect_equal(m$corP$cstruc, "propto")
  expect_true(length(m$params$sigma) > 0L,
              info = "propto must produce at least one sigma value")
})
