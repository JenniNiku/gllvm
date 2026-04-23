# Data frame of allowed families, methods and links

raw <- data.frame(
  family = c(
    "binomial","binomial",
    "poisson","negative.binomial","negative.binomial1","ZIP","ZINB",
    "ZIB","ZIB",
    "ZNIB","ZNIB",
    "tweedie","ordinal","ordinal",
    "gaussian",
    "gamma","exponential",
    "beta",
    "orderedBeta","orderedBeta",
    "betaH", "beta.binomial"
  ),
  method = c(
    "VA/LA","EVA",
    "VA/LA","VA/EVA/LA","VA/LA","VA/LA","VA/LA",
    "VA","LA",
    "VA","LA",
    "EVA/VA/LA","VA/LA", "EVA",
    "VA/LA",
    "VA/LA","VA/LA",
    "LA/EVA",
    "VA","EVA",
    "EVA/VA/LA",
    "LA"
  ),
  link = c(
    "probit/logit/cloglog","probit/logit", 
    "log","log","log","log","log",
    "probit/logit/cloglog","probit/logit",
    "probit/logit/cloglog","probit/logit",
    "log","probit/logit","logit",
    "identity",
    "log","log",
    "probit/logit",
    "probit/logit","logit",
    "probit/logit",
    "probit/logit/cloglog"
  ),
  stringsAsFactors = FALSE
)

split_rows <- function(df) {
  res <- list()
  k <- 1
  
  for (i in seq_len(nrow(df))) {
    methods <- unlist(strsplit(df$method[i], "/"))
    links   <- unlist(strsplit(df$link[i], "/"))
    
    comb <- expand.grid(
      family = df$family[i],
      method = methods,
      link   = links,
      stringsAsFactors = FALSE
    )
    
    res[[k]] <- comb
    k <- k + 1
  }
  
  do.call(rbind, res)
}

gllvmFML_allowed <- split_rows(raw)

gllvmFML_allowed$key <- paste(gllvmFML_allowed$family, gllvmFML_allowed$method, gllvmFML_allowed$link, sep = "|")

usethis::use_data(gllvmFML_allowed, internal = TRUE, overwrite = TRUE)

