treatment_validation_methods <- c("drdid_panel", "drdid_imp_panel", "reg_did_panel",
  "ipw_did_panel", "std_ipw_did_panel", "twfe_did_panel", "drdid_rc",
  "drdid_rc1", "drdid_imp_rc", "drdid_imp_rc1", "reg_did_rc", "ipw_did_rc",
  "std_ipw_did_rc", "twfe_did_rc")
treatment_validation_fit <- function(method,D) {
  post <- rep(rep(0:1,each=20),2)
  y <- sin(seq_len(80))+rep(0:1,each=40)*post
  args <- list(D=D,covariates=NULL,inffunc=TRUE)
  if(grepl("panel",method)) args <- c(args,list(y1=y,y0=rep(0,80))) else
    args <- c(args,list(y=y,post=post))
  do.call(get(method,mode="function"),args)
}
test_that("low-level estimators reject nonbinary treatment without recoding", {
  D <- rep(0:1,each=40)
  bad <- list(D+1,2*D,D-.5,replace(D,1,Inf),replace(D,1,NA_real_),factor(D))
  for(method in treatment_validation_methods) for(d in bad)
    expect_error(treatment_validation_fit(method,d),"D must contain only 0|Panel inputs must be complete")
})
test_that("supported logical treatment in OR estimators retains its behavior", {
  D <- rep(0:1,each=40)
  for(method in c("reg_did_panel", "reg_did_rc")) {
    a <- treatment_validation_fit(method,D); b <- treatment_validation_fit(method,as.logical(D))
    expect_equal(a$ATT,b$ATT,tolerance=1e-12)
    expect_equal(a$att.inf.func,b$att.inf.func,tolerance=1e-12)
  }
})
test_that("public wrappers reject alternative numeric treatment codes", {
  d <- expand.grid(id=1:40,time=0:1)
  d$group <- as.numeric(d$id>20)+1; d$outcome <- sin(d$id)+d$time*(d$group-1)
  for(panel in c(TRUE,FALSE)) for(method in c("drdid","ordid","ipwdid")) {
    args <- list(yname="outcome",tname="time",idname="id",dname="group",data=d,panel=panel)
    expect_error(do.call(get(method,mode="function"),args),"D must contain only 0")
  }
})

test_that("numeric column-matrix treatment retains vector behavior", {
  D <- rep(0:1,each=40)
  for(method in treatment_validation_methods) {
    a <- treatment_validation_fit(method,D)
    b <- treatment_validation_fit(method,matrix(D,ncol=1))
    expect_equal(a$ATT,b$ATT,tolerance=1e-12)
    expect_equal(a$att.inf.func,b$att.inf.func,tolerance=1e-12)
  }
})
