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

  # Structured
  X <- microbialdata$Xenv[1:30,1:3]
  StudyDesign = data.frame(Site = factor(microbialdata$Xenv$Site[1:30]), Soiltype = microbialdata$Xenv$Soiltype[1:30])
  fr2<-gllvm(y, family = "negative.binomial", seed = 999, studyDesign = StudyDesign, row.eff = ~(1|Site), num.lv = 1)
  fr3<-gllvm(y, family = "negative.binomial", seed = 999, studyDesign = StudyDesign, row.eff = ~(1|Site)+(1|Soiltype), num.lv = 1)
  StudyDesign = data.frame(Site = factor(microbialdata$Xenv$Site[1:30]), pH = microbialdata$Xenv$pH[1:30])
  fr4<-gllvm(y,X, family = "negative.binomial", seed = 999, studyDesign = StudyDesign,  row.eff = ~pH, num.lv = 1)
  fr5<-gllvm(y,X, family = "negative.binomial", seed = 999, studyDesign = StudyDesign,  row.eff = ~(0+pH|Site), num.lv = 1)
  
  result<-c(-0.34, 0.29, 0.41, -0.02, 0)
  names(result)<-c("AB3", "1|site", "1|Site", "pH", "pH|Site")
  resultf3<-c(0.30, 0.29)
  names(resultf3)<-c("1|Site", "1|Soiltype")
  
  expect_true(round(fr0$params$row.params[1], digits = 2)- result[1]<0.1)
  expect_true(round(fr1$params$sigma, digits = 2)- result[2]<0.1)
  expect_true(round(fr2$params$sigma, digits = 2)- result[3]<0.1)
  expect_true(all(round(fr3$params$sigma, digits = 2)- resultf3<0.1))
  expect_true(round(fr4$params$row.params.fixed, digits = 2)- result[4]<0.1)
  expect_true(round(fr5$params$sigma, digits = 2)- result[5]<0.1)
})

test_that("binomial works", {
  data(microbialdata)
  y <- microbialdata$Y[, order(colMeans(microbialdata$Y > 0), decreasing = TRUE)[160:175]]
  y01<-(y>0)*1
  fb0<-gllvm(y01, family = binomial(link = "logit"), seed = 999, method = "LA", num.lv = 1)
  fb2<-gllvm(y01, family = binomial(link = "probit"), seed = 999)
  expect_true(is.finite(fb0$logL))
  expect_true(is.finite(fb2$logL))
})

test_that("ZIP works", {
  data(eSpider)
  spider <- eSpider
  spider$abund <- eSpider$abund[eSpider$nonNA,]
  y <- spider$abund[order(rowSums(spider$abund>0), decreasing = TRUE)[1:20],1:8]
  fz0<-gllvm(y, family = "ZIP", seed = 999, method = "LA", num.lv = 1)
  expect_equal( length(fz0$params$beta0), 8 )
  expect_true( is.finite(fz0$logL))
})

test_that("quadratic models work", {
  data(eSpider)
  spider <- eSpider
  spider$abund <- spider$abund[spider$nonNA,]
  spider$x <- spider$X[spider$nonNA,]
  X <- scale(spider$x)
  y <- spider$abund
  fq0<-gllvm(y, num.lv = 2, family = "poisson", seed = 999)
  fq1<-gllvm(y, X, num.lv = 2, family = "poisson", seed = 999)
  fq2<-gllvm(y, X, num.lv = 2, family = "poisson", row.eff="random", seed = 999)
  expect_true(is.finite(fq0$logL))
  expect_true(is.finite(fq1$logL))
  expect_true(is.finite(fq2$logL))
})

test_that("constrained ordination models work", {
  data(eSpider)
  spider <- eSpider
  spider$abund <- spider$abund[spider$nonNA,]
  spider$x <- spider$X[spider$nonNA,]
  X <- scale(spider$x)
  y <- spider$abund
  suppressWarnings({fc0<-gllvm(y, X, num.RR = 2, family = "poisson", seed = 999)})
  fc1<-gllvm(y, X, num.RR = 2, family = "poisson", seed = 999, randomB="LV")
  fc2<-gllvm(y, X, num.RR = 2, family = "poisson", seed = 999, randomB="LV", row.eff="random")
  fc3<-gllvm(y, X, num.RR = 2, quadratic=T, family = "poisson", seed = 9226, randomB="LV", row.eff="random")
  expect_true(is.finite(fc0$logL))
  expect_true(is.finite(fc1$logL))
  expect_true(is.finite(fc2$logL))
  expect_true(is.finite(fc3$logL))
})

test_that("concurrent ordination models work", {
  data(eSpider)
  spider <- eSpider
  spider$abund <- spider$abund[spider$nonNA,]
  spider$x <- spider$X[spider$nonNA,]
  X <- scale(spider$x)
  y <- spider$abund
  fc0<-gllvm(y, X, num.lv.c = 2, family = "poisson", seed = 999)
  fc1<-gllvm(y, X, num.lv.c = 2, family = "poisson", seed = 999, randomB="LV")
  #this has a warning for overfitting that can be ignored
  suppressWarnings(fc2<-gllvm(y, X, num.lv.c = 2, family = "poisson", seed = 999, randomB="LV", row.eff="random"))
  fc3<-gllvm(y, X, num.lv.c = 2, quadratic=T, family = "poisson", seed = 999, randomB="LV")
  expect_true(is.finite(fc0$logL))
  expect_true(is.finite(fc1$logL))
  expect_true(is.finite(fc2$logL))
  expect_true(is.finite(fc3$logL))
})

test_that("phylogenetic models work", {
  data(eSpider)
  spider <- eSpider
  spider$abund <- spider$abund[spider$nonNA,]
  spider$x <- spider$X[spider$nonNA,]
  X <- scale(spider$x)
  colMat=matrix(rnorm(12*12),12,12)
  colMat<-colMat%*%t(colMat)
  dist=abs(colMat)
  colnames(dist)<-colnames(colMat)<-1:12
  
  expect_error({model<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),num.lv=0, sd.errors = FALSE)})
  
  colnames(colMat)<-row.names(colMat)<-colnames(dist)<-row.names(dist)<-colnames(spider$abund)
  expect_error({model<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=colMat,colMat.rho.struct="term",num.lv=0, sd.errors=FALSE)})
  
  expect_warning({model<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=colMat,Ab.struct="diagonal", num.lv=0, sd.errors = FALSE)})

  #trait models are not yet tested in the following, just basic infrastructure for phylogenetic rando effects
  suppressWarnings({suppressMessages(invisible(capture.output({
  model11<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="single",Ab.struct="diagonal",Ab.struct.rank=12,beta0com=TRUE)
  model12<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="single",Ab.struct="blockdiagonal",Ab.struct.rank=12,beta0com=TRUE)
  model13<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="single",Ab.struct="MNdiagonal",Ab.struct.rank=12,beta0com=TRUE)
  model14<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="single",Ab.struct="MNunstructured",Ab.struct.rank=12,beta0com=TRUE)
  model15<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="single",Ab.struct="diagonalCL2",Ab.struct.rank=12,beta0com=TRUE)
  model16<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="single",Ab.struct="diagonalCL1",Ab.struct.rank=12,beta0com=TRUE)
  model17<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="single",Ab.struct="CL1",Ab.struct.rank=12,beta0com=TRUE)
  model18<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="single",Ab.struct="CL2",Ab.struct.rank=12,beta0com=TRUE)
  model19<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="single",Ab.struct="unstructured",Ab.struct.rank=1e3,beta0com=TRUE)

  model21<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="term",Ab.struct="diagonal",Ab.struct.rank=12,beta0com=TRUE)
  model22<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="term",Ab.struct="blockdiagonal",Ab.struct.rank=12,beta0com=TRUE)
  model23<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="term",Ab.struct="MNdiagonal",Ab.struct.rank=12,beta0com=TRUE)
  model24<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="term",Ab.struct="MNunstructured",Ab.struct.rank=12,beta0com=TRUE)
  model25<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="term",Ab.struct="diagonalCL2",Ab.struct.rank=12,beta0com=TRUE)
  model26<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="term",Ab.struct="diagonalCL1",Ab.struct.rank=12,beta0com=TRUE)
  model27<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="term",Ab.struct="CL1",Ab.struct.rank=12,beta0com=TRUE)
  model28<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="term",Ab.struct="CL2",Ab.struct.rank=12,beta0com=TRUE)
  model29<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="term",Ab.struct="unstructured",Ab.struct.rank=1e3,beta0com=TRUE)

  model31<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="single",Ab.struct="diagonal",Ab.struct.rank=1,beta0com=TRUE)
  model32<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="single",Ab.struct="blockdiagonal",Ab.struct.rank=1,beta0com=TRUE)
  model33<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="single",Ab.struct="MNdiagonal",Ab.struct.rank=1,beta0com=TRUE)
  model34<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="single",Ab.struct="MNunstructured",Ab.struct.rank=1,beta0com=TRUE)
  model35<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="single",Ab.struct="diagonalCL2",Ab.struct.rank=1,beta0com=TRUE)
  model36<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="single",Ab.struct="diagonalCL1",Ab.struct.rank=1,beta0com=TRUE)
  model37<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="single",Ab.struct="CL1",Ab.struct.rank=1,beta0com=TRUE)
  model38<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="single",Ab.struct="CL2",Ab.struct.rank=1,beta0com=TRUE)
  model39<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="single",Ab.struct="unstructured",Ab.struct.rank=1,beta0com=TRUE)

  model41<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="term",Ab.struct="diagonal",Ab.struct.rank=1,beta0com=TRUE)
  model42<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="term",Ab.struct="blockdiagonal",Ab.struct.rank=1,beta0com=TRUE)
  model43<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="term",Ab.struct="MNdiagonal",Ab.struct.rank=1,beta0com=TRUE)
  model44<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="term",Ab.struct="MNunstructured",Ab.struct.rank=1,beta0com=TRUE)
  model45<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="term",Ab.struct="diagonalCL2",Ab.struct.rank=1,beta0com=TRUE)
  model46<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="term",Ab.struct="diagonalCL1",Ab.struct.rank=1,beta0com=TRUE)
  model47<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="term",Ab.struct="CL1",Ab.struct.rank=1,beta0com=TRUE)
  model48<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="term",Ab.struct="CL2",Ab.struct.rank=1,beta0com=TRUE)
  model49<-gllvm(spider$abund,X=X,formula=~nocorr(ConWate+CovMoss|1),family="poisson",colMat=list(colMat,dist=dist),nn.colMat=3,num.lv=0,colMat.rho.struct="term",Ab.struct="unstructured",Ab.struct.rank=1,beta0com=TRUE)
  })))})
  expect_true(is.finite(model11$logL))
  expect_true(is.finite(model12$logL))
  expect_true(is.finite(model13$logL))
  expect_true(is.finite(model14$logL))
  expect_true(is.finite(model15$logL))
  expect_true(is.finite(model16$logL))
  expect_true(is.finite(model17$logL))
  expect_true(is.finite(model18$logL))
  expect_true(is.finite(model19$logL))
  expect_true(is.finite(model21$logL))
  expect_true(is.finite(model22$logL))
  expect_true(is.finite(model23$logL))
  expect_true(is.finite(model24$logL))
  expect_true(is.finite(model25$logL))
  expect_true(is.finite(model26$logL))
  expect_true(is.finite(model27$logL))
  expect_true(is.finite(model28$logL))
  expect_true(is.finite(model29$logL))
  expect_true(is.finite(model31$logL))
  expect_true(is.finite(model32$logL))
  expect_true(is.finite(model33$logL))
  expect_true(is.finite(model34$logL))
  expect_true(is.finite(model35$logL))
  expect_true(is.finite(model36$logL))
  expect_true(is.finite(model37$logL))
  expect_true(is.finite(model38$logL))
  expect_true(is.finite(model39$logL))
  expect_true(is.finite(model41$logL))
  expect_true(is.finite(model42$logL))
  expect_true(is.finite(model43$logL))
  expect_true(is.finite(model44$logL))
  expect_true(is.finite(model45$logL))
  expect_true(is.finite(model46$logL))
  expect_true(is.finite(model47$logL))
  expect_true(is.finite(model48$logL))
  expect_true(is.finite(model49$logL))
  
})
