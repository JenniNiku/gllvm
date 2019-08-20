context("test-fitgllvm")

test_that("basic data fitting works", {
  data(spider)
  f0<-gllvm(spider$abund, family = poisson(), seed = 999)
  f1<-gllvm(spider$abund, family = "negative.binomial", seed = 999)
  f2<-gllvm(spider$abund, spider$x, formula = ~soil.dry + bare.sand, family = "negative.binomial", seed = 999)
  expect_equal(round(f0$logL, digits = 1),-873.8)
  expect_equal(round(f1$logL, digits = 1),-761.7)
  expect_equal(round(f2$logL, digits = 1),-712.8)
})


test_that("fourth corner models works", {
  data(antTraits)
  ff0<-gllvm(antTraits$abund,antTraits$env, antTraits$traits, family = "negative.binomial", seed = 999)
  ff1<-gllvm(antTraits$abund,antTraits$env, antTraits$traits, formula = ~Bare.ground + Bare.ground:Pilosity, family = "negative.binomial", seed = 999)
  expect_equal(round(ff0$logL, digits = 1),-1851.3)
  expect_equal(round(ff1$logL, digits = 1),-1888.1)
})

test_that("row effects works", {
  data("spider")
  fr0<-gllvm(spider$abund, family = "negative.binomial", seed = 999, row.eff = "fixed", num.lv = 1)
  fr1<-gllvm(spider$abund, family = "negative.binomial", seed = 999, row.eff = "random", num.lv = 0)
  result<-c(0.21, 1.60)
  names(result)<-c("Row2", "sigma")
  expect_equal(round(fr0$params$row.params[2], digits = 2), result[1])
  expect_equal(round(fr1$params$sigma, digits = 2), result[2])
})

test_that("binomial works", {
  data(spider)
  y01<-(spider$abund>0)*1
  fb0<-gllvm(y01, family = binomial(link = "logit"), seed = 999, method = "LA")
  fb1<-gllvm(y01, family = binomial(link = "probit"), seed = 999, method = "LA")
  fb2<-gllvm(y01, family = binomial(link = "probit"), seed = 999)
  expect_equal(round(fb0$logL, digits = 0),-94)
  expect_equal(round(fb1$logL, digits = 1),-94.2)
  expect_equal(round(fb2$logL, digits = 1),-184.3)
})

test_that("ZIP works", {
  data(spider)
  fz0<-gllvm(spider$abund[1:10,c(1:4,6)], family = "ZIP", seed = 999, method = "LA")
  expect_equal( length(fz0$params$beta0), 5 )
  expect_true( is.finite(fz0$logL))
})
