#This document is shared across cobalt, WeightIt, and optweight as a symbolic link.
#Any edits will be automatically synced across all folders. Make sure functions work
#in all packages!
#The original file is in cobalt/R/.

#Strings
word.list <- function(word.list = NULL, and.or = c("and", "or"), is.are = FALSE, quotes = FALSE) {
    #When given a vector of strings, creates a string of the form "a and b"
    #or "a, b, and c"
    #If is.are, adds "is" or "are" appropriately
    L <- length(word.list)
    if (quotes) word.list <- vapply(word.list, function(x) paste0("\"", x, "\""), character(1L))
    if (L == 0) {
        out <- ""
        attr(out, "plural") = FALSE
    }
    else {
        word.list <- word.list[!word.list %in% c(NA_character_, "")]
        L <- length(word.list)
        if (L == 0) {
            out <- ""
            attr(out, "plural") = FALSE
        }
        else if (L == 1) {
            out <- word.list
            if (is.are) out <- paste(out, "is")
            attr(out, "plural") = FALSE
        }
        else {
            and.or <- match_arg(and.or)
            if (L == 2) {
                out <- paste(word.list, collapse = paste0(" ", and.or," "))
            }
            else {
                out <- paste(paste(word.list[seq_len(L-1)], collapse = ", "),
                             word.list[L], sep = paste0(", ", and.or," "))
                
            }
            if (is.are) out <- paste(out, "are")
            attr(out, "plural") = TRUE
        }
        
        
    }
    return(out)
}
expand.grid_string <- function(..., collapse = "") {
    return(apply(expand.grid(...), 1, paste, collapse = collapse))
}
num_to_superscript <- function(x) {
    nums <- setNames(c("\u2070",
                       "\u00B9",
                       "\u00B2",
                       "\u00B3",
                       "\u2074",
                       "\u2075",
                       "\u2076",
                       "\u2077",
                       "\u2078",
                       "\u2079"),
                     as.character(0:9))
    x <- as.character(x)
    splitx <- strsplit(x, "", fixed = TRUE)
    supx <- sapply(splitx, function(y) paste0(nums[y], collapse = ""))
    return(supx)
}
ordinal <- function(x) {
    if (!is.numeric(x) || !is.vector(x) || is_null(x)) stop("x must be a numeric vector.")
    if (length(x) > 1) return(vapply(x, ordinal, character(1L)))
    else {
        x0 <- abs(x)
        out <- paste0(x0, switch(substring(x0, nchar(x0), nchar(x0)),
                                 "1" = "st",
                                 "2" = "nd",
                                 "3" = "rd",
                                 "th"))
        if (sign(x) == -1) out <- paste0("-", out)
        
        return(out)
    }
}
round_df_char <- function(df, digits, pad = "0", na_vals = "") {
    if (NROW(df) == 0) return(df)
    nas <- is.na(df)
    if (!is.data.frame(df)) df <- as.data.frame.matrix(df, stringsAsFactors = FALSE)
    infs <- sapply(df, function(x) !is.na(x) & (x == Inf | x == -Inf), simplify = "array")
    rn <- rownames(df)
    cn <- colnames(df)
    df <- as.data.frame(lapply(df, function(col) {
        if (suppressWarnings(all(!is.na(as.numeric(as.character(col)))))) {
            as.numeric(as.character(col))
        } else {
            col
        }
    }), stringsAsFactors = FALSE)
    nums <- vapply(df, is.numeric, logical(1))
    o.negs <- sapply(1:NCOL(df), function(x) if (nums[x]) df[[x]] < 0 else rep(FALSE, length(df[[x]])))
    df[nums] <- round(df[nums], digits = digits)
    
    df[nas | infs] <- ""
    
    df <- as.data.frame(lapply(df, format, scientific = FALSE, justify = "none"), stringsAsFactors = FALSE)
    
    for (i in which(nums)) {
        if (any(grepl(".", df[[i]], fixed = TRUE))) {
            s <- strsplit(df[[i]], ".", fixed = TRUE)
            lengths <- lengths(s)
            digits.r.of.. <- sapply(seq_along(s), function(x) {
                if (lengths[x] > 1) nchar(s[[x]][lengths[x]])
                else 0 })
            df[[i]] <- sapply(seq_along(df[[i]]), function(x) {
                if (df[[i]][x] == "") ""
                else if (lengths[x] <= 1) {
                    paste0(c(df[[i]][x], rep(".", pad == 0), rep(pad, max(digits.r.of..) - digits.r.of..[x] + as.numeric(pad != 0))),
                           collapse = "")
                }
                else paste0(c(df[[i]][x], rep(pad, max(digits.r.of..) - digits.r.of..[x])),
                            collapse = "")
            })
        }
    }
    
    df[o.negs & df == 0] <- paste0("-", df[o.negs & df == 0])
    
    # Insert NA placeholders
    df[nas] <- na_vals
    df[infs] <- "N/A"
    
    if (length(rn) > 0) rownames(df) <- rn
    if (length(cn) > 0) names(df) <- cn
    
    return(df)
}
text.box.plot <- function(range.list, width = 12) {
    full.range <- range(unlist(range.list))
    ratio = diff(full.range)/(width+1)
    rescaled.range.list <- lapply(range.list, function(x) round(x/ratio))
    rescaled.full.range <- round(full.range/ratio)
    d <- as.data.frame(matrix(NA_character_, ncol = 3, nrow = length(range.list),
                              dimnames = list(names(range.list), c("Min", paste(rep(" ", width + 1), collapse = ""), "Max"))),
                       stringsAsFactors = FALSE)
    d[,"Min"] <- vapply(range.list, function(x) x[1], numeric(1L))
    d[,"Max"] <- vapply(range.list, function(x) x[2], numeric(1L))
    for (i in seq_len(nrow(d))) {
        spaces1 <- rescaled.range.list[[i]][1] - rescaled.full.range[1]
        #|
        dashes <- max(0, diff(rescaled.range.list[[i]]) - 2)
        #|
        spaces2 <- max(0, diff(rescaled.full.range) - (spaces1 + 1 + dashes + 1))
        
        d[i, 2] <- paste0(paste(rep(" ", spaces1), collapse = ""), "|", paste(rep("-", dashes), collapse = ""), "|", paste(rep(" ", spaces2), collapse = ""))
    }
    return(d)
}
equivalent.factors <- function(f1, f2) {
    return(nunique(f1) == nunique(interaction(f1, f2)))
}
paste. <- function(..., collapse = NULL) {
    #Like paste0 but with sep = ".'
    paste(..., sep = ".", collapse = collapse)
}

#Numbers
check_if_zero <- function(x) {
    # this is the default tolerance used in all.equal
    tolerance <- .Machine$double.eps^0.5
    # If the absolute deviation between the number and zero is less than
    # the tolerance of the floating point arithmetic, then return TRUE.
    # This means, to me, that I can treat the number as 0 rather than
    # -3.20469e-16 or some such.
    abs(x - 0) < tolerance
}
between <- function(x, range, inclusive = TRUE, na.action = FALSE) {
    if (!all(is.numeric(x))) stop("x must be a numeric vector.", call. = FALSE)
    if (length(range) != 2) stop("range must be of length 2.", call. = FALSE)
    if (anyNA(range) || !is.numeric(range)) stop("range must contain numeric entries only.", call. = FALSE)
    range <- sort(range)
    
    if (anyNA(x)) {
        if (length(na.action) != 1 || !is.atomic(na.action)) stop("na.action must be an atomic vector of length 1.", call. = FALSE)
    }
    if (inclusive) out <- ifelse(is.na(x), na.action, x >= range[1] & x <= range[2])
    else out <- ifelse(is.na(x), na.action, x > range[1] & x < range[2])
    
    return(out)
}

#Statistics
binarize <- function(variable) {
    nas <- is.na(variable)
    if (!is_binary(variable[!nas])) stop(paste0("Cannot binarize ", deparse(substitute(variable)), ": more than two levels."))
    if (is.character(variable)) variable <- factor(variable)
    variable.numeric <- as.numeric(variable)
    if (0 %in% variable.numeric) zero <- 0
    else zero <- min(variable.numeric, na.rm = TRUE)
    newvar <- setNames(ifelse(!nas & variable.numeric==zero, 0L, 1L), names(variable))
    newvar[nas] <- NA_integer_
    return(newvar)
}
ESS <- function(w) {
    sum(w)^2/sum(w^2)
}
center <- function(x, na.rm = TRUE, at = NULL) {
    if (!is.numeric(x)) warning("x is not numeric and will not be centered.")
    else if (is.matrix(x)) x <- apply(x, 2, center, na.rm = na.rm, at = at)
    else {
        if (is_null(at)) at <- mean(x, na.rm = na.rm)
        else if (!is.numeric(at)) stop("at must be numeric.")
        x <- x - at
    }
    return(x)
}
w.m <- function(x, w = NULL, na.rm = TRUE) {
    if (is_null(w)) w <- as.numeric(!is.na(x))
    return(sum(x*w, na.rm=na.rm)/sum(w, na.rm=na.rm))
}
w.v <- function(x, w = NULL) {
    # return(sum(w*(x-w.m(x, w))^2, na.rm=TRUE)/(sum(w, na.rm=TRUE)-1))
    # return(sum(w*(x-w.m(x, w))^2, na.rm=TRUE) / w.cov.scale(w))
    w.cov(x, x, w = w)
}
w.cov <- function(x, y, w = NULL) {
    if (is_null(w)) w <- as.numeric(!is.na(x))
    wmx <- w.m(x, w)
    wmy <- w.m(y, w)
    w.cov.scale <- (sum(w, na.rm = TRUE)^2 - sum(w^2, na.rm = TRUE)) / sum(w, na.rm = TRUE)
    #wcov <- sum(w*(x - wmx)*(y - wmy), na.rm = TRUE)/sum(w, na.rm = TRUE)
    wcov <- sum(w*(x - wmx)*(y - wmy), na.rm = TRUE) / w.cov.scale
    return(wcov)
}
col.w.m <- function(mat, w = NULL, na.rm = TRUE) {
    if (is_null(w)) {
        w <- 1
        w.sum <- apply(mat, 2, function(x) sum(!is.na(x)))
    }
    else {
        w.sum <- rep(sum(w, na.rm = na.rm), ncol(mat))
    }
    return(colSums(mat*w, na.rm = na.rm)/w.sum)
}
col.w.v <- function(mat, w = NULL, na.rm = TRUE) {
    if (is_null(w)) {
        w <- rep(1, nrow(mat))
    }
    w.cov.scale <- (sum(w, na.rm = TRUE)^2 - sum(w^2, na.rm = TRUE)) / sum(w, na.rm = TRUE)
    return(colSums(t((t(mat) - col.w.m(mat, w, na.rm = na.rm))^2) * w, na.rm = na.rm) / w.cov.scale)
}
coef.of.var <- function(x, pop = TRUE, na.rm = TRUE) {
    if (na.rm) x <- x[!is.na(x)]
    if (pop) sqrt(mean((x-mean(x))^2))/mean(x)
    else sd(x)/mean(x)
}
mean.abs.dev <- function(x, na.rm = TRUE) {
    if (na.rm) x <- x[!is.na(x)]
    mean(abs(x - mean(x)))
}
geom.mean <- function(y, na.rm = TRUE) {
    exp(mean(log(y[is.finite(log(y))]), na.rm = na.rm))
}


#Formulas
is.formula <- function(f, sides = NULL) {
    res <- is.name(f[[1]])  && deparse(f[[1]]) %in% c( '~', '!') &&
        length(f) >= 2
    if (is_not_null(sides) && is.numeric(sides) && sides %in% c(1,2)) {
        res <- res && length(f) == sides + 1
    }
    return(res)
}

#treat/covs
get.covs.and.treat.from.formula <- function(f, data = NULL, terms = FALSE, sep = "", ...) {
    A <- list(...)
    
    #Check if data exists
    if (is_not_null(data) && is.data.frame(data)) {
        data.specified <- TRUE
    }
    else {
        data <- environment(f)
        data.specified <- FALSE
    }
    
    env <- environment(f)
    
    tt <- terms(f)
    
    #Check if response exists
    if (is.formula(tt, 2)) {
        resp.vars.mentioned <- as.character(tt)[2]
        resp.vars.failed <- vapply(resp.vars.mentioned, function(v) {
            null.or.error(try(eval(parse(text=v)[[1]], data, env), silent = TRUE))
        }, logical(1L))
        
        if (any(resp.vars.failed)) {
            if (is_null(A[["treat"]])) stop(paste0("The given response variable, \"", as.character(tt)[2], "\", is not a variable in ", word.list(c("data", "the global environment")[c(data.specified, TRUE)], "or"), "."), call. = FALSE)
            tt <- delete.response(tt)
        }
    }
    else resp.vars.failed <- TRUE
    
    if (any(!resp.vars.failed)) {
        treat.name <- resp.vars.mentioned[!resp.vars.failed][1]
        treat <- eval(parse(text=treat.name)[[1]], data, env)
    }
    else {
        treat <- A[["treat"]]
        treat.name <- NULL
    }
    
    #Check if RHS variables exist
    tt.covs <- delete.response(tt)
    rhs.vars.mentioned.lang <- attr(tt.covs, "variables")[-1]
    rhs.vars.mentioned <- vapply(rhs.vars.mentioned.lang, deparse, character(1L))
    rhs.vars.failed <- vapply(rhs.vars.mentioned, function(v) {
        null.or.error(try(eval(parse(text=v)[[1]], data, env), silent = TRUE))
    }, logical(1L))
    
    if (any(rhs.vars.failed)) {
        stop(paste0(c("All variables in formula must be variables in data or objects in the global environment.\nMissing variables: ",
                      paste(rhs.vars.mentioned[rhs.vars.failed], collapse=", "))), call. = FALSE)
        
    }
    
    rhs.term.labels <- attr(tt.covs, "term.labels")
    rhs.term.orders <- attr(tt.covs, "order")
    
    rhs.df <- vapply(rhs.vars.mentioned, function(v) {
        d <- try(eval(parse(text=v)[[1]], data, env), silent = TRUE)
        is.data.frame(d) || is.matrix(d)
    }, logical(1L))
    
    if (any(rhs.df)) {
        if (any(rhs.vars.mentioned[rhs.df] %in% unlist(lapply(rhs.term.labels[rhs.term.orders > 1], function(x) strsplit(x, ":", fixed = TRUE))))) {
            stop("Interactions with data.frames are not allowed in the input formula.", call. = FALSE)
        }
        addl.dfs <- setNames(lapply(rhs.vars.mentioned[rhs.df], function(x) {as.data.frame(eval(parse(text=x)[[1]], data, env))}),
                             rhs.vars.mentioned[rhs.df])
        
        for (i in rhs.term.labels[rhs.term.labels %in% rhs.vars.mentioned[rhs.df]]) {
            ind <- which(rhs.term.labels == i)
            rhs.term.labels <- append(rhs.term.labels[-ind],
                                      values = names(addl.dfs[[i]]),
                                      after = ind - 1)
        }
        new.form <- as.formula(paste("~", paste(rhs.term.labels, collapse = " + ")))
        
        tt.covs <- terms(new.form)
        if (data.specified) data <- do.call("cbind", unname(c(addl.dfs, list(data))))
        else data <- do.call("cbind", unname(addl.dfs))
    }
    
    #Get model.frame, report error
    mf.covs <- quote(stats::model.frame(tt.covs, data,
                                        drop.unused.levels = TRUE,
                                        na.action = "na.pass"))
    
    tryCatch({covs <- eval(mf.covs)},
             error = function(e) {stop(conditionMessage(e), call. = FALSE)})
    
    if (is_not_null(treat.name) && treat.name %in% names(covs)) stop("The variable on the left side of the formula appears on the right side too.", call. = FALSE)
    
    if (is_null(rhs.vars.mentioned)) {
        covs <- data.frame(Intercept = rep(1, if (is_null(treat)) 1 else length(treat)))
    }
    else attr(tt.covs, "intercept") <- 0
    
    if (s <- !identical(sep, "")) {
        if (!is.character(sep) || length(sep) > 1) stop("sep must be a string of length 1.", call. = FALSE)
        original.covs.levels <- setNames(vector("list", ncol(covs)), names(covs))
        for (i in names(covs)) {
            if (is.character(covs[[i]])) covs[[i]] <- factor(covs[[i]])
            if (is.factor(covs[[i]])) {
                original.covs.levels[[i]] <- levels(covs[[i]])
                levels(covs[[i]]) <- paste0(sep, original.covs.levels[[i]])
            }
        }
    }
    
    #Get full model matrix with interactions too
    covs.matrix <- model.matrix(tt.covs, data = covs,
                                contrasts.arg = lapply(Filter(is.factor, covs),
                                                       contrasts, contrasts=FALSE))
    
    if (s) {
        for (i in names(covs)) {
            if (is.factor(covs[[i]])) {
                levels(covs[[i]]) <- original.covs.levels[[i]]
            }
        }
    }
    
    if (!terms) attr(covs, "terms") <- NULL
    
    return(list(reported.covs = covs,
                model.covs = covs.matrix,
                treat = treat,
                treat.name = treat.name))
}
get.treat.type <- function(treat) {
    #Returns treat with treat.type attribute
    nunique.treat <- nunique(treat)
    if (nunique.treat == 2) {
        treat.type <- "binary"
    }
    else if (nunique.treat < 2) {
        stop("The treatment must have at least two unique values.", call. = FALSE)
    }
    else if (is.factor(treat) || is.character(treat)) {
        treat.type <- "multinomial"
        treat <- factor(treat)
    }
    else {
        treat.type <- "continuous"
    }
    attr(treat, "treat.type") <- treat.type
    return(treat)
}
process.s.weights <- function(s.weights, data = NULL) {
    #Process s.weights
    if (is_not_null(s.weights)) {
        if (!(is.character(s.weights) && length(s.weights) == 1) && !is.numeric(s.weights)) {
            stop("The argument to s.weights must be a vector or data frame of sampling weights or the (quoted) names of variables in data that contain sampling weights.", call. = FALSE)
        }
        if (is.character(s.weights) && length(s.weights)==1) {
            if (is_null(data)) {
                stop("s.weights was specified as a string but there was no argument to data.", call. = FALSE)
            }
            else if (s.weights %in% names(data)) {
                s.weights <- data[[s.weights]]
            }
            else stop("The name supplied to s.weights is not the name of a variable in data.", call. = FALSE)
        }
    }
    return(s.weights)
}

#Uniqueness
nunique <- function(x, nmax = NA, na.rm = TRUE) {
    if (is_null(x)) return(0)
    else {
        if (na.rm) x <- x[!is.na(x)]
        if (is.factor(x)) return(nlevels(x))
        else return(length(unique(x, nmax = nmax)))
    }
    
}
nunique.gt <- function(x, n, na.rm = TRUE) {
    if (missing(n)) stop("n must be supplied.")
    if (n < 0) stop("n must be non-negative.")
    if (is_null(x)) FALSE
    else {
        if (na.rm) x <- x[!is.na(x)]
        if (n == 1 && is.numeric(x)) !check_if_zero(max(x) - min(x))
        else if (length(x) < 2000) nunique(x) > n
        else tryCatch(nunique(x, nmax = n) > n, error = function(e) TRUE)
    }
}
all_the_same <- function(x) !nunique.gt(x, 1)
is_binary <- function(x) !nunique.gt(x, 2)

#R Processing
is_ <- function(x, types, stop = FALSE, arg.to = FALSE) {
    s1 <- deparse(substitute(x))
    if (is_not_null(x)) {
        for (i in types) {
            if (is_not_null(get0(paste.("is", i)))) {
                it.is <- get0(paste.("is", i))(x)
            }
            else it.is <- inherits(x, i)
            if (it.is) break
        }
    }
    else it.is <- FALSE
    
    if (stop) {
        if (!it.is) {
            s0 <- ifelse(arg.to, "The argument to ", "")
            s2 <- ifelse(any(types %in% c("factor", "character", "numeric", "logical")),
                         "vector", "")
            stop(paste0(s0, s1, " must be a ", word.list(types, and.or = "or"), " ", s2, "."), call. = FALSE)
        }
    }
    else {
        return(it.is)
    }
}
is_null <- function(x) length(x) == 0L
is_not_null <- function(x) !is_null(x)
clear_null <- function(x) {
    x[vapply(x, is_null, logical(1L))] <- NULL
    return(x)
}
probably.a.bug <- function() {
    fun <- paste(deparse(sys.call(-1)), collapse = "\n")
    stop(paste0("An error was produced and is likely a bug. Please let the maintainer know a bug was produced by the function\n",
                fun), call. = FALSE)
}
`%nin%` <- function(x, table) is.na(match(x, table, nomatch = NA_integer_))
null.or.error <- function(x) {is_null(x) || class(x) == "try-error"}
match_arg <- function(arg, choices, several.ok = FALSE) {
    #Replaces match.arg() but gives cleaner error message and processing
    #of arg.
    if (missing(arg))
        stop("No argument was supplied to match_arg.", call. = FALSE)
    arg.name <- deparse(substitute(arg))
    
    if (missing(choices)) {
        formal.args <- formals(sys.function(sysP <- sys.parent()))
        choices <- eval(formal.args[[as.character(substitute(arg))]],
                        envir = sys.frame(sysP))
    }
    
    if (is.null(arg))
        return(choices[1L])
    else if (!is.character(arg))
        stop(paste0("'", arg.name, "' must be NULL or a character vector"), call. = FALSE)
    if (!several.ok) {
        if (identical(arg, choices))
            return(arg[1L])
        if (length(arg) > 1L)
            stop(paste0("'", arg.name, "' must be of length 1"), call. = FALSE)
    }
    else if (is_null(arg))
        stop(paste0("'", arg.name, "' must be of length >= 1"), call. = FALSE)
    
    i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
    if (all(i == 0L))
        stop(paste0("'", arg.name, "' should be one of ", word.list(choices, and.or = "or", quotes = TRUE), "."),
             call. = FALSE)
    i <- i[i > 0L]
    if (!several.ok && length(i) > 1)
        stop("there is more than one match in 'match_arg'")
    choices[i]
}
