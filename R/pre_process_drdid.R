#' @importFrom Rcpp sourceCpp
NULL

# Function that pre-process the data to use with drdid, ipwdid and ordid

pre_process_drdid <- function(yname,
                              tname,
                              idname,
                              dname,
                              xformla = NULL,
                              data,
                              panel = TRUE,
                              normalized = TRUE,
                              estMethod = c("imp", "trad"),
                              weightsname = NULL,
                              boot = FALSE,
                              boot.type =  c("weighted", "multiplier"),
                              nboot = NULL,
                              inffunc = FALSE) {
  #-----------------------------------------------------------------------------
  # Data pre-processing and error checking
  #-----------------------------------------------------------------------------
  # ####
  # set estimation method
  estMethod <- estMethod[1]
  # set bootstrap type
  boot.type <- boot.type[1]


  # Flag for estMethod
  if ( (estMethod!="imp") && (estMethod!="trad")) {
    warning("estMethod = ",estMethod,  " is not supported. Using 'imp'.")
    estMethod <- "imp"
  }
  # Flag for boot.type
  if (boot){
    if ( (boot.type!="weighted") && (boot.type!="multiplier")) {
      warning("boot.type = ",boot.type,  " is not supported. Using 'weighted'.")
      boot.type <- "weighted"
    }
  }

  # Flag for normalized
  if ( (normalized != TRUE) && (normalized != FALSE)) {
    warning("normalized = ",normalized,  " is not supported. Using 'TRUE'.")
    normalized <- TRUE
  }

  # make sure dataset is a data.frame
  dta <- data
  # this gets around RStudio's default of reading data as tibble
  if (!all( class(dta) == "data.frame")) {
    #warning("class of data object was not data.frame; converting...")
    dta <- as.data.frame(dta)
  }


  # Flag for yname
  if ( !is.element(yname, base::colnames(dta))) {
    stop("yname = ",yname,  " could not be found in the data provided.")
  }
  # Flag for tname
  if ( !is.element(tname, base::colnames(dta))) {
    stop("tname = ",tname,  " could not be found in the data provided.")
  }
  # Flag for dname
  if ( !is.element(dname, base::colnames(dta))) {
    stop("dname = ",dname,  " could not be found in the data provided.")
  }
  # Flag for idname
  if ( !is.element(idname, base::colnames(dta))) {
    stop("idname = ",idname,  " could not be found in the data provided.")
  }
  #
  # if(panel){
  #   if ( !is.element(idname, base::colnames(dta))) {
  #     stop("idname = ",idname,  " could not be found in the data provided.")
  #   }
  # }
  # Flag for weightsname
  if(!is.null(weightsname)){
    if ( !is.element(weightsname, base::colnames(dta))) {
      stop("weightsname = ",weightsname,  " could not be found in the data provided.")
    }
  }
  # Flag for covariates
  #if (!is.null(xformla)) {
  #  if(!all(is.element(BMisc::rhs.vars(xformla), base::colnames(dta)))){
  #    stop("Some of the covariate names could not be found in the data provided. Double check the covariate names in xformla.")
  #  }
  #}




  # make sure time periods are numeric
  if (! (is.numeric(dta[, tname])) ) {
    stop("data[, tname] must be numeric. Please convert it.")

  }

  #  make sure dname is numeric
  if (! (is.numeric(dta[, dname])) ) {
    stop("data[, dname] must be numeric. Please convert it.")

  }

  validate_treatment_indicator(dta[, dname], allow_missing = TRUE)

  #  make sure id is numeric
  if (! is.null(idname)){
    #  make sure id is numeric
    if (! (is.numeric(dta[, idname])) ) {
      stop("data[, idname] must be numeric. Please convert it.")
    }

    # Check if idname is unique by tname
    n_id_year = base::all( base::table(dta[, idname], dta[, tname]) <= 1)
    if (! n_id_year) stop("The value of idname must be the unique (by tname)")

    # Check if dname is unique by idname
    checkTreatmentUniqueness(dta, idname, dname)
  }

  # figure out the time periods
  # list of dates from smallest to largest
  tlist <- unique(dta[,tname])[base::order(unique(dta[,tname]))]
  if ( length(tlist)!=2) {
    stop("drdid package only covers the cases with two time periods (pre and post) and two treatment groups (d=1 if treated at post, and d=0 if not treated in both pre and post).
           See package `did' for the cases with multiple groups and/or multiple time periods.")
  }

  # list of groups from smallest to largest
  glist <- unique(dta[,dname])[base::order(unique(dta[,dname]))]
  if ( length(glist)!=2) {
    stop("drdid only work for case with two time periods and two treatment groups (d=1 if treated at post, and d=0 if not treated in both pre and post).
         See package `did' for the cases with multiple groups and/or multiple time periods.")
  }


  # put in blank xformla if no covariates
  if (is.null(xformla)) {
    xformla <- ~1
  }

  # # If repeated cross section, allow for null idname
  # if(is.null(idname) & (panel == FALSE)){
  #   dta$id <- seq(1:nrow(dta))
  #   idname <- "id"
  # }

  # sort data with respect to id and time
  dta <- dta[base::order(dta[,idname], dta[,tname]),]

  # If panel, make it a balanced data set
  if (panel) {
    dta <- BMisc::make_balanced_panel(dta, idname, tname)
  }


  # Keep estimator inputs separate from user columns. Formula evaluation must
  # see the original data, including columns named y, D, post, or w.
  y <- dta[[yname]]
  D <- dta[[dname]]
  post <- as.numeric(dta[[tname]] == tlist[2])
  w <- if (is.null(weightsname)) rep(1, nrow(dta)) else dta[[weightsname]]
  # matrix of covariates
  covariates <- stats::model.matrix(xformla,
                                    stats::model.frame(xformla,
                                                       data = dta,
                                                       na.action=na.pass) )

  # check against very small groups
  gsize <- stats::aggregate(D,
                            by=list(D),
                            function(x) length(x)/length(tlist))

  # how many in each group before give warning
  # 5 is just a buffer, could pick something else, but seems to work fine
  reqsize <- length(BMisc::rhs_vars(xformla)) + 5

  # which groups to warn about
  gsize <- subset(gsize,
                  gsize$x < reqsize) # x is name of column from gsize

  # warn if some groups are small
  if (nrow(gsize) > 0) {
    stop(paste0(" Either treatment or the comparison group in your dataset is very small. Inference is not feasible "))
  }
  #-----------------------------------------------------------------------------
  # setup data in panel case
  if (panel) {
    # make it a balanced data set
    dta <- BMisc::make_balanced_panel(dta, idname, tname)

    # Keep IDs separate so they cannot enter the covariate matrix.
    panel_ids <- dta[[idname]]

    # Only use this smaller dataset
    dta <- as.data.frame(cbind(y = y,
                               D = D,
                               post = post,
                               w = w,
                               covariates))

    # Check if we have missing values
    missing_X_flag <- base::anyNA(dta[,-c(1:4)])
    missing_D_flag <- base::anyNA(dta$D)
    missing_post_flag <- base::anyNA(dta$post)
    missing_y_flag <- base::anyNA(dta$y)
    missing_w_flag <- base::anyNA(dta$w)

    # Missingness in either period excludes the whole panel unit. Removing
    # individual rows would destroy the alignment of the pre/post outcomes.
    incomplete_ids <- unique(panel_ids[!stats::complete.cases(dta)])
    keep_rows <- !(panel_ids %in% incomplete_ids)
    missing_message <- paste0(" We are dropping all observations for ",
                              length(incomplete_ids), " affected panel unit(s).")

    # Warning messages
    if(missing_X_flag) warning("Some covariates are missing (NA).", missing_message)
    if(missing_D_flag) warning("Some treatment indicators are missing (NA).", missing_message)
    if(missing_post_flag) warning("Some post-treatment indicators are missing (NA).", missing_message)
    if(missing_w_flag) warning("Some weights are missing (NA).", missing_message)
    if(missing_y_flag) warning("Some outcomes are missing (NA).", missing_message)

    dta <- dta[keep_rows, , drop = FALSE]
    panel_ids <- panel_ids[keep_rows]

    # Match by ID rather than relying on equal pre/post vector lengths.
    pre_rows <- which(dta$post == 0)
    post_rows <- which(dta$post == 1)
    pre_match <- match(panel_ids[post_rows], panel_ids[pre_rows])
    if (length(pre_rows) != length(post_rows) || anyNA(pre_match)) {
      stop("Error: panel observations could not be matched by id across time periods.")
    }
    pre_rows <- pre_rows[pre_match]

    # Count retained units, including groups that disappeared after filtering.
    retained_gsize <- table(factor(dta$D[post_rows], levels = glist))
    if (any(retained_gsize < reqsize)) {
      stop("Either treatment or the comparison group in your dataset is very small after removing incomplete panel units. Inference is not feasible.")
    }

    covariates <- base::as.matrix(dta[, -c(1:4), drop = FALSE])
    # Check invariance only after removing incomplete units, and compare each
    # retained unit with its matched observation in the other period. Building
    # separate model matrices would omit missing rows before matching IDs.
    if (!all(covariates[pre_rows, , drop = FALSE] ==
             covariates[post_rows, , drop = FALSE])) {
      stop("Error: covariates should be time invariant, and there should be no missing data.")
    }
    if (!all(dta$D[pre_rows] == dta$D[post_rows])) {
      stop("Error: group indicator must be time invariant.")
    }
    if (!all(dta$w[pre_rows] == dta$w[post_rows])) {
      stop("Error: weights must be time invariant.")
    }
    # Drop collinear covariates
    qr.covariates <-  base::qr(covariates,
                               tol=1e-6,
                               LAPACK = FALSE)
    rnk_covariates <- qr.covariates$rank
    if (rnk_covariates == 0L) {
      stop("The covariate matrix has no linearly independent columns.")
    }
    keep_x <- qr.covariates$pivot[seq_len(rnk_covariates)]
    # Check if any covariates were dropped due to collinearity
    dropped_covariates <- setdiff(colnames(covariates), colnames(covariates)[keep_x])
    if (length(dropped_covariates) > 0) {
      warning("The following covariates were dropped due to collinearity: ", paste(dropped_covariates, collapse = ", "))
    }

    covariates <- covariates[, keep_x, drop = FALSE]




    validate_sampling_weights(dta$w[post_rows], length(post_rows))

    out <- list(y1 = dta$y[post_rows],
                y0 = dta$y[pre_rows],
                D = dta$D[post_rows],
                covariates = covariates[post_rows, , drop = FALSE],
                # Normalize weights
                i.weights = dta$w[post_rows]/mean(dta$w[post_rows]),
                panel = panel,
                estMethod = estMethod,
                normalized = normalized,
                boot = boot,
                boot.type = boot.type,
                nboot = nboot,
                inffunc = inffunc

    )
  }
  else {
    # Only use this smaller dataset
    dta <- as.data.frame(cbind(y = y,
                               D = D,
                               post = post,
                               w = w,
                               covariates))

    # Check if we have missing values
    missing_X_flag <- base::anyNA(dta[,-c(1:4)])
    missing_D_flag <- base::anyNA(dta$D)
    missing_post_flag <- base::anyNA(dta$post)
    missing_y_flag <- base::anyNA(dta$y)
    missing_w_flag <- base::anyNA(dta$w)

    # Warning messages
    if(missing_X_flag) warning("Some covariates are missing (NA). We are dropping those observations")
    if(missing_D_flag) warning("Some treatment indicators are missing (NA). We are dropping those observations")
    if(missing_post_flag) warning("Some post-treatment indicators are missing (NA). We are dropping those observations")
    if(missing_w_flag) warning("Some weights are missing (NA). We are dropping those observations")
    if(missing_y_flag) warning("Some outcomes are missing (NA). We are dropping those observations")

    # Remove NAs
    dta <- dta[stats::complete.cases(dta), ]
    covariates <- base::as.matrix(dta[,-c(1:4)])

    # Drop collinear covariates
    qr.covariates <-  base::qr(covariates,
                               tol=1e-6,
                               LAPACK = FALSE)
    rnk_covariates <- qr.covariates$rank
    keep_x <- qr.covariates$pivot[seq_len(rnk_covariates)]
    # Check if any covariates were dropped due to collinearity (match the panel path,
    # which warns; the RC path previously dropped them silently)
    dropped_covariates <- setdiff(colnames(covariates), colnames(covariates[, keep_x]))
    if (length(dropped_covariates) > 0) {
      warning("The following covariates were dropped due to collinearity: ", paste(dropped_covariates, collapse = ", "))
    }
    covariates <- covariates[,keep_x]


    validate_sampling_weights(dta$w, nrow(dta))

    out <- list(y = dta$y,
                D = dta$D,
                post = dta$post,
                covariates = covariates,
                # Normalize weights
                i.weights = dta$w/mean(dta$w),
                panel = panel,
                estMethod = estMethod,
                normalized = normalized,
                boot = boot,
                boot.type = boot.type,
                nboot = nboot,
                inffunc = inffunc
    )
  }

  return(out)
}




