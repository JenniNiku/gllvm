context("test-fitgllvm")

test_that("basic data fitting works", {
  data(microbialdata)
  X <- microbialdata$Xenv[1:30,]
  y <- microbialdata$Y[1:30, order(colMeans(microbialdata$Y > 0), decreasing = TRUE)[21:35]]
  f0<-gllvm(y, family = poisson(), seed = 999)
  f1<-gllvm(y, family = "negative.binomial", seed = 999)
  f2<-gllvm(y, X, formula = ~pH + Phosp, family = "negative.binomial", seed = 999)
  expect_true(round(mean(f0$params$beta0), digits = 1)-2<0.01)
  expect_true(round(mean(f1$params$beta0), digits = 1)-2<0.01)
  expect_true(round(mean(f2$params$Xcoef[,1]), digits = 1)-0.1<0.01)
})


test_that("fourth corner models works", {
  data(microbialdata)
  X <- microbialdata$Xenv[1:30,2:3]
  y <- microbialdata$Y[1:30, order(colMeans(microbialdata$Y > 0), decreasing = TRUE)[21:35]]
  TR <- matrix(rnorm(15)); colnames(TR) <- "t1"
  ff0<-gllvm(y, X, TR=TR, family = "negative.binomial", seed = 999)
  ff1<-gllvm(y, X, TR=TR, formula = ~pH + pH:t1, family = "negative.binomial", seed = 999)
  expect_true(is.finite(ff0$logL))
  expect_true(is.finite(ff1$logL))
})

test_that("row effects works", {
  data(microbialdata)
  y <- microbialdata$Y[1:30, order(colMeans(microbialdata$Y > 0), decreasing = TRUE)[21:35]]
  fr0<-gllvm(y, family = "negative.binomial", seed = 999, row.eff = "fixed", num.lv = 1)
  fr1<-gllvm(y, family = "negative.binomial", seed = 999, row.eff = "random", num.lv = 0)
  result<-c(0.34, 0.29)
  names(result)<-c("AB3", "sigma")
  expect_true(round(fr0$params$row.params[2], digits = 2)- result[1]<0.1)
  expect_true(round(fr1$params$sigma, digits = 2)- result[2]<0.1)
})

test_that("binomial works", {
  data(microbialdata)
  y <- microbialdata$Y[1:30, order(colMeans(microbialdata$Y > 0), decreasing = TRUE)[21:35]]
  y01<-(y>0)*1
  fb0<-gllvm(y01, family = binomial(link = "logit"), seed = 999, method = "LA", num.lv = 1)
  fb2<-gllvm(y01, family = binomial(link = "probit"), seed = 999)
  expect_true(is.finite(fb0$logL))
  expect_true(is.finite(fb2$logL))
})

test_that("ZIP works", {
  data(microbialdata)
  y <- microbialdata$Y[1:10, order(colMeans(microbialdata$Y > 0), decreasing = TRUE)[301:306]]
  fz0<-gllvm(y, family = "ZIP", seed = 999, method = "LA")
  expect_equal( length(fz0$params$beta0), 6 )
  expect_true( is.finite(fz0$logL))
})
