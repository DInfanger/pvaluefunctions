# Internal helpers shared by conf_dist() and the cdist_*() functions.

# TRUE if all supplied objects have the same length.
same_length <- function(...) {
  length(unique(lengths(list(...)))) == 1L
}

# TRUE if `trans` is the identity transformation, given either as a string
# (case-insensitive) or as a function.
is_identity_trans <- function(trans) {
  if (is.function(trans)) {
    return(identical(trans, base::identity))
  }
  is.character(trans) && length(trans) == 1L && identical(tolower(trans), "identity")
}

# Resolve the `trans` argument of conf_dist() to a list with
#   name: "identity", "exp" or the (lower-case) name of a user function
#   fun:  the function that is applied to the values
# `trans` can be a function or the name of a function. Names are looked up in
# `env` first (usually the caller of conf_dist()), then in the package
# namespace and its parents. For backward compatibility, the lower-case version
# of the name is tried if the name as given cannot be found.
resolve_trans <- function(trans, env = parent.frame(), call = rlang::caller_env()) {
  if (is.function(trans)) {
    name <- if (identical(trans, base::exp)) {
      "exp"
    } else if (identical(trans, base::identity)) {
      "identity"
    } else {
      "custom"
    }

    return(list(name = name, fun = trans))
  }

  if (!is.character(trans) || length(trans) != 1L || is.na(trans)) {
    cli::cli_abort(
      c(
        "{.arg trans} must be a function or the name of a function.",
        x = "You supplied {.cls {class(trans)}}."
      ),
      call = call
    )
  }

  key <- tolower(trans)

  if (key %in% c("identity", "exp")) {
    return(list(name = key, fun = get(key, envir = baseenv(), mode = "function")))
  }

  for (candidate in unique(c(trans, key))) {
    fun <- get0(candidate, envir = env, mode = "function", inherits = TRUE)

    if (is.null(fun)) {
      fun <- get0(
        candidate,
        envir = environment(resolve_trans),
        mode = "function",
        inherits = TRUE
      )
    }

    if (!is.null(fun)) {
      return(list(name = key, fun = fun))
    }
  }

  cli::cli_abort("Function {.fn {trans}} was not found.", call = call)
}
