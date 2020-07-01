#This document is shared across cobalt, WeightIt, and optweight

#Strings
word_list <- function(word.list = NULL, and.or = c("and", "or"), is.are = FALSE, quotes = FALSE) {
    #When given a vector of strings, creates a string of the form "a and b"
    #or "a, b, and c"
    #If is.are, adds "is" or "are" appropriately
    L <- length(word.list)
    if (quotes) {
        if (as.integer(quotes) == 2) word.list <- vapply(word.list, function(x) paste0("\"", x, "\""), character(1L))
        else if (as.integer(quotes) == 1) word.list <- vapply(word.list, function(x) paste0("\'", x, "\'"), character(1L))
        else stop("'quotes' must be boolean, 1, or 2.")
    }
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
firstup <- function(x) {
    #Capitalize first letter
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
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
    if (!is.numeric(x) || !is.vector(x) || is_null(x)) stop("'x' must be a numeric vector.")
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
text_box_plot <- function(range.list, width = 12) {
    full.range <- range(unlist(range.list))
    ratio = diff(full.range)/(width+1)
    rescaled.range.list <- lapply(range.list, function(x) round(x/ratio))
    rescaled.full.range <- round(full.range/ratio)
    d <- make_df(c("Min", paste(rep(" ", width + 1), collapse = ""), "Max"),
                 names(range.list),
                 "character")
    d[["Min"]] <- vapply(range.list, function(x) x[1], numeric(1L))
    d[["Max"]] <- vapply(range.list, function(x) x[2], numeric(1L))
    for (i in seq_len(nrow(d))) {
        spaces1 <- rescaled.range.list[[i]][1] - rescaled.full.range[1]
        #|
        dashes <- max(c(0, diff(rescaled.range.list[[i]]) - 2))
        #|
        spaces2 <- max(c(0, diff(rescaled.full.range) - (spaces1 + 1 + dashes + 1)))

        d[i, 2] <- paste0(paste(rep(" ", spaces1), collapse = ""), "|", paste(rep("-", dashes), collapse = ""), "|", paste(rep(" ", spaces2), collapse = ""))
    }
    return(d)
}
equivalent.factors <- function(f1, f2) {
    nu1 <- nunique(f1)
    nu2 <- nunique(f2)
    if (nu1 == nu2) {
        return(nu1 == nunique(paste.(f1, f2)))
    }
    else {
        return(FALSE)
    }
}
equivalent.factors2 <- function(f1, f2) {
    return(qr(cbind(1, as.numeric(f1), as.numeric(f2)))$rank == 2)
}
paste. <- function(..., collapse = NULL) {
    #Like paste0 but with sep = ".'
    paste(..., sep = ".", collapse = collapse)
}
wrap <- function(s, nchar, ...) {
    vapply(s, function(s_) {
        x <- strwrap(s_, width = nchar, ...)
        paste(x, collapse = "\n")
    }, character(1L))
}
strsplits <- function(x, splits, fixed = TRUE, ...) {
    #Link strsplit but takes multiple split values.
    #Only works for one string at a time (in x).
    for (split in splits) x <- unlist(strsplit(x, split, fixed = TRUE, ...))
    return(x[x != ""]) # Remove empty values
}
c.factor <- function(..., recursive=TRUE) {
    #c() for factors
    unlist(list(...), recursive=recursive)
}
can_str2num <- function(x) {
    nas <- is.na(x)
    suppressWarnings(x_num <- as.numeric(as.character(x[!nas])))
    return(!anyNA(x_num))
}
str2num <- function(x) {
    nas <- is.na(x)
    suppressWarnings(x_num <- as.numeric(as.character(x)))
    x_num[nas] <- NA
    return(x_num)
}

#Numbers
check_if_zero <- function(x) {
    # this is the default tolerance used in all.equal
    tolerance <- .Machine$double.eps^0.5
    abs(x) < tolerance
}
between <- function(x, range, inclusive = TRUE, na.action = FALSE) {
    if (!all(is.numeric(x))) stop("'x' must be a numeric vector.", call. = FALSE)
    if (length(range) != 2) stop("'range' must be of length 2.", call. = FALSE)
    if (anyNA(range) || !is.numeric(range)) stop("'range' must contain numeric entries only.", call. = FALSE)

    if (range[2] < range[1]) range <- c(range[2], range[1])

    if (anyNA(x)) {
        if (length(na.action) != 1 || !is.atomic(na.action)) stop("'na.action' must be an atomic vector of length 1.", call. = FALSE)
    }
    if (inclusive) out <- ifelse(is.na(x), na.action, x >= range[1] & x <= range[2])
    else out <- ifelse(is.na(x), na.action, x > range[1] & x < range[2])

    return(out)
}
max_ <- function(..., na.rm = TRUE) {
    if (!any(is.finite(unlist(list(...))))) NA_real_
    else max(..., na.rm = na.rm)
}
min_ <- function(..., na.rm = TRUE) {
    if (!any(is.finite(unlist(list(...))))) NA_real_
    else min(..., na.rm = na.rm)
}
check_if_int <- function(x) {
    #Checks if integer-like
    if (is.integer(x)) rep(TRUE, length(x))
    else if (is.numeric(x)) check_if_zero(x - round(x))
    else rep(FALSE, length(x))
}

#Statistics
binarize <- function(variable, zero = NULL, one = NULL) {
    if (!is_binary(variable)) stop(paste0("Cannot binarize ", deparse1(substitute(variable)), ": more than two levels."))
    if (is.character(variable)) {
        variable <- factor(variable, nmax = 2)
        unique.vals <- levels(variable)
    }
    else unique.vals <- unique(variable, nmax = 2)

    variable.numeric <- as.numeric(variable)

    if (is_null(zero)) {
        if (is_null(one)) {
            if (0 %in% variable.numeric) zero <- 0
            else zero <- min(variable.numeric, na.rm = TRUE)
            return(setNames(as.integer(variable.numeric != zero), names(variable)))
        }
        else {
            if (one %in% unique.vals) return(setNames(as.integer(variable == one), names(variable)))
            else stop("The argument to 'one' is not the name of a level of variable.", call. = FALSE)
        }
    }
    else {
        if (zero %nin% unique.vals) stop("The argument to 'zero' is not the name of a level of variable.", call. = FALSE)
        return(setNames(as.integer(variable != zero), names(variable)))
    }
}
ESS <- function(w) {
    sum(w)^2/sum(w^2)
}
center <- function(x, at = NULL, na.rm = TRUE) {
    if (is.data.frame(x)) {
        x <- as.matrix.data.frame(x)
        type <- "df"
    }
    if (!is.numeric(x)) stop("'x' must be numeric.")
    else if (is.array(x) && length(dim(x)) > 2) stop("'x' must be a numeric or matrix-like (not array).")
    else if (!is.matrix(x)) {
        x <- matrix(x, ncol = 1)
        type <- "vec"
    }
    else type <- "matrix"
    if (is_null(at)) at <- colMeans(x, na.rm = na.rm)
    else if (length(at) %nin% c(1, ncol(x))) stop("'at' is not the right length.")
    out <- x - matrix(at, byrow = TRUE, ncol = ncol(x), nrow = nrow(x))
    if (type == "df") out <- as.data.frame.matrix(out)
    else if (type == "vec") out <- drop(out)
    return(out)
}
w.m <- function(x, w = NULL, na.rm = TRUE) {
    if (is_null(w)) w <- rep(1, length(x))
    if (anyNA(x)) w[is.na(x)] <- NA
    return(sum(x*w, na.rm=na.rm)/sum(w, na.rm=na.rm))
}
col.w.m <- function(mat, w = NULL, na.rm = TRUE) {
    if (is_null(w)) w <- 1
    w.sum <- colSums(w*!is.na(mat))
    return(colSums(mat*w, na.rm = na.rm)/w.sum)
}
col.w.v <- function(mat, w = NULL, bin.vars = NULL, na.rm = TRUE) {
    if (!is.matrix(mat)) {
        if (is.data.frame(mat)) {
            if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
                stop("'mat' must be a numeric matrix.")
            }
            else mat <- data.matrix(mat)
        }
        else if (is.numeric(mat)) {
            mat <- matrix(mat, ncol = 1)
        }
        else stop("'mat' must be a numeric matrix.")
    }

    if (is_null(bin.vars)) bin.vars <- rep(FALSE, ncol(mat))
    else if (length(bin.vars) != ncol(mat) || anyNA(as.logical(bin.vars))) {
        stop("'bin.vars' must be a logical vector with length equal to the number of columns of 'mat'.", call. = FALSE)
    }
    bin.var.present <- any(bin.vars)
    non.bin.vars.present <- any(!bin.vars)

    var <- setNames(numeric(ncol(mat)), colnames(mat))
    if (is_null(w)) {
        if (non.bin.vars.present) {
            den <- colSums(!is.na(mat[, !bin.vars, drop = FALSE])) - 1
            var[!bin.vars] <- colSums(center(mat[, !bin.vars, drop = FALSE])^2, na.rm = na.rm)/den
        }
        if (bin.var.present) {
            means <- colMeans(mat[, bin.vars, drop = FALSE], na.rm = na.rm)
            var[bin.vars] <- means * (1 - means)
        }
    }
    else if (na.rm && anyNA(mat)) {
        # n <- nrow(mat)
        w <- array(w, dim = dim(mat))
        w[is.na(mat)] <- NA
        s <- colSums(w, na.rm = na.rm)
        w <- mat_div(w, s)
        if (non.bin.vars.present) {
            x <- sqrt(w[, !bin.vars, drop = FALSE]) * center(mat[, !bin.vars, drop = FALSE],
                                                             at = colSums(w[, !bin.vars, drop = FALSE] * mat[, !bin.vars, drop = FALSE], na.rm = na.rm))
            var[!bin.vars] <- colSums(x*x, na.rm = na.rm)/(1 - colSums(w[, !bin.vars, drop = FALSE]^2, na.rm = na.rm))
        }
        if (bin.var.present) {
            means <- colSums(w[, bin.vars, drop = FALSE] * mat[, bin.vars, drop = FALSE], na.rm = na.rm)
            var[bin.vars] <- means * (1 - means)
        }
    }
    else {
        if (is_null(w)) w <- rep(1, nrow(mat))
        w <- w/sum(w)
        if (non.bin.vars.present) {
            x <- sqrt(w) * center(mat[, !bin.vars, drop = FALSE],
                                  at = colSums(w * mat[, !bin.vars, drop = FALSE], na.rm = na.rm))
            var[!bin.vars] <- colSums(x*x, na.rm = na.rm)/(1 - sum(w^2))
        }
        if (bin.var.present) {
            means <- colSums(w * mat[, bin.vars, drop = FALSE], na.rm = na.rm)
            var[bin.vars] <- means * (1 - means)
        }
    }
    return(var)
}
col.w.cov <- function(mat, y, w = NULL, na.rm = TRUE) {
    if (!is.matrix(mat)) {
        if (is_null(w)) return(cov(mat, y, use = if (na.rm) "pair" else "everything"))
        else mat <- matrix(mat, ncol = 1)
    }
    if (is_null(w)) {
        y <- array(y, dim = dim(mat))
        if (anyNA(mat)) y[is.na(mat)] <- NA
        if (anyNA(y)) mat[is.na(y)] <- NA
        den <- colSums(!is.na(mat*y)) - 1
        cov <- colSums(center(mat, na.rm = na.rm)*center(y, na.rm = na.rm), na.rm = na.rm)/den
    }
    else if (na.rm && anyNA(mat)) {
        n <- nrow(mat)
        w <- array(w, dim = dim(mat))
        w[is.na(mat)] <- NA_real_
        s <- colSums(w, na.rm = na.rm)
        w <- mat_div(w, s)
        x <- w * center(mat, at = colSums(w * mat, na.rm = na.rm))
        cov <- colSums(x*y, na.rm = na.rm)/(1 - colSums(w^2, na.rm = na.rm))
    }
    else {
        n <- nrow(mat)
        w <- w/sum(w)
        x <- w * center(mat, at = colSums(w * mat, na.rm = na.rm))
        cov <- colSums(x*y, na.rm = na.rm)/(1 - sum(w^2))
    }
    return(cov)
}
col.w.r <- function(mat, y, w = NULL, s.weights = NULL, bin.vars = NULL, na.rm = TRUE) {
    if (is_null(w) && is_null(s.weights)) return(cor(mat, y, w, use = if (na.rm) "pair" else "everything"))
    else {
        cov <- col.w.cov(mat, y = y, w = w, na.rm = na.rm)
        den <- sqrt(col.w.v(mat, w = s.weights, bin.vars = bin.vars, na.rm = na.rm)) *
            sqrt(col.w.v(y, w = s.weights, na.rm = na.rm))
        return(cov/den)
    }
}
coef.of.var <- function(x, pop = TRUE) {
    if (pop) sqrt(mean_fast((x-mean_fast(x, TRUE))^2, TRUE))/mean_fast(x, TRUE)
    else sd(x)/mean_fast(x, TRUE)
}
mean.abs.dev <- function(x) {
    mean_fast(abs(x - mean_fast(x, TRUE)), TRUE)
}
rms <- function(x) {
    sqrt(mean_fast(x^2))
}
geom.mean <- function(y) {
    exp(mean_fast(log(y[is.finite(log(y))]), TRUE))
}
mat_div <- function(mat, vec) {
    mat/vec[col(mat)]
}
abs_ <- function(x, ratio = FALSE) {
    if (ratio) pmax(x, 1/x)
    else (abs(x))
}
mean_fast <- function(x, nas.possible = FALSE) {
    #Equal to mean(x, na.rm = TRUE) but faster
    #Set no.nas = FALSE if it's possible there are NAs
    if (nas.possible && anyNA(x)) {
        s <- sum(x, na.rm = TRUE)
        n <- sum(!is.na(x))
        return(s/n)
    }
    s <- sum(x)
    n <- length(x)
    return(s/n)
}
bw.nrd <- function(x) {
    #R's bw.nrd doesn't always work, but bw.nrd0 does
    bw.nrd0(x)*1.06/.9
}

#Formulas
subbars <- function(term) {
    if (is.name(term) || !is.language(term))
        return(term)
    if (length(term) == 2) {
        term[[2]] <- subbars(term[[2]])
        return(term)
    }

    if (is.call(term) && (term[[1]] == as.name("|") || term[[1]] == as.name("||"))) {
        term[[1]] <- as.name("+")
    }
    for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
    return(term)
}

#treat/covs
get.covs.and.treat.from.formula <- function(f, data = NULL, terms = FALSE, sep = "", ...) {
    A <- list(...)

    #Check if data exists
    if (is_not_null(data)) {
        if (is.data.frame(data)) {
            data.specified <- TRUE
        }
        else {
            warning("The argument supplied to data is not a data.frame object. This may causes errors or unexpected results.", call. = FALSE)
            data <- environment(f)
            data.specified <- FALSE
        }
    }
    else {
        data <- environment(f)
        data.specified <- FALSE
    }

    env <- environment(f)

    if (!is.formula(f)) stop("'f' must be a formula.")

    eval.model.matrx <- identical(f, f <- subbars(f))

    tryCatch(tt <- terms(f, data = data),
             error = function(e) {
                 if (conditionMessage(e) == "'.' in formula and no 'data' argument") {
                     stop("'.' is not allowed in formulas.", call. = FALSE)
                 }
                 else stop(conditionMessage(e), call. = FALSE)
             })

    #Check if response exists
    if (is.formula(tt, 2)) {
        resp.vars.mentioned <- as.character(tt)[2]
        resp.vars.failed <- vapply(resp.vars.mentioned, function(v) {
            test <- tryCatch(eval(parse(text=v), data, env), error = function(e) e)
            if (inherits(test, "simpleError")) {
                if (conditionMessage(test) == paste0("object '", v, "' not found")) return(TRUE)
                else stop(test)
            }
            else if (is_null(test)) return(TRUE)
            else return(FALSE)
        }, logical(1L))

        if (any(resp.vars.failed)) {
            if (is_null(A[["treat"]])) stop(paste0("The given response variable, \"", as.character(tt)[2], "\", is not a variable in ", word_list(c("data", "the global environment")[c(data.specified, TRUE)], "or"), "."), call. = FALSE)
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
    rhs.vars.mentioned <- vapply(rhs.vars.mentioned.lang, deparse1, character(1L))
    rhs.vars.failed <- vapply(rhs.vars.mentioned, function(v) {
        test <- tryCatch(eval(parse(text=v), data, env), error = function(e) e)
        if (inherits(test, "simpleError")) {
            if (conditionMessage(test) == paste0("object '", v, "' not found")) return(TRUE)
            else stop(test)
        }
        else if (is_null(test)) return(TRUE)
        else return(FALSE)
    }, logical(1L))

    if (any(rhs.vars.failed)) {
        stop(paste0(c("All variables in 'formula' must be variables in 'data' or objects in the global environment.\nMissing variables: ",
                      paste(rhs.vars.mentioned[rhs.vars.failed], collapse=", "))), call. = FALSE)

    }

    rhs.term.labels <- attr(tt.covs, "term.labels")
    rhs.term.orders <- attr(tt.covs, "order")

    rhs.df <- setNames(vapply(rhs.vars.mentioned, function(v) {
        is_(try(eval(parse(text=v)[[1]], data, env), silent = TRUE),
            c("data.frame", "matrix", "rms"))
    }, logical(1L)), rhs.vars.mentioned)

    rhs.term.labels.list <- setNames(as.list(rhs.term.labels), rhs.term.labels)
    if (any(rhs.df)) {
        if (any(rhs.vars.mentioned[rhs.df] %in% unlist(lapply(rhs.term.labels[rhs.term.orders > 1], function(x) strsplit(x, ":", fixed = TRUE))))) {
            stop("Interactions with data.frames are not allowed in the input formula.", call. = FALSE)
        }
        addl.dfs <- setNames(lapply(rhs.vars.mentioned[rhs.df], function(x) {
            df <- eval(parse(text=x)[[1]], data, env)
            if (is_(df, "rms")) {
                if (length(dim(df)) == 2L) class(df) <- "matrix"
                df <- setNames(as.data.frame(as.matrix(df)), attr(df, "colnames"))
            }
            else if (can_str2num(colnames(df))) colnames(df) <- paste(x, colnames(df), sep = sep)
            return(as.data.frame(df))
        }),
        rhs.vars.mentioned[rhs.df])

        for (i in rhs.term.labels[rhs.term.labels %in% rhs.vars.mentioned[rhs.df]]) {
            ind <- which(rhs.term.labels == i)
            rhs.term.labels <- append(rhs.term.labels[-ind],
                                      values = names(addl.dfs[[i]]),
                                      after = ind - 1)
            rhs.term.labels.list[[i]] <- names(addl.dfs[[i]])
        }

        if (data.specified) data <- do.call("cbind", unname(c(addl.dfs, list(data))))
        else data <- do.call("cbind", unname(addl.dfs))
    }

    if (is_null(rhs.term.labels)) {
        new.form <- as.formula("~ 1")
        tt.covs <- terms(new.form)
        covs <- data.frame(Intercept = rep(1, if (is_null(treat)) 1 else length(treat)))
        if (is_not_null(treat.name) && treat.name == "Intercept") {
            names(covs) <- "Intercept_"
        }
    }
    else {
        new.form.char <- paste("~", paste(vapply(names(rhs.term.labels.list), function(x) {
            if (x %in% rhs.vars.mentioned[rhs.df]) paste0("`", rhs.term.labels.list[[x]], "`", collapse = " + ")
            else rhs.term.labels.list[[x]]
            # try.form <- try(as.formula(paste("~", x)), silent = TRUE)
            # if (null_or_error(try.form) || (grepl("^", x, fixed = TRUE) && !startsWith(x, "I("))) {
            #     paste0("`", x, "`")
            # }
            # else x
        } , character(1L)), collapse = " + "))

        new.form <- as.formula(new.form.char)
        tt.covs <- terms(new.form)
        attr(tt.covs, "intercept") <- 0

        #Get model.frame, report error
        mf.covs <- quote(stats::model.frame(tt.covs, data,
                                            drop.unused.levels = TRUE,
                                            na.action = "na.pass"))

        tryCatch({covs <- eval(mf.covs)},
                 error = function(e) {stop(conditionMessage(e), call. = FALSE)})

        if (is_not_null(treat.name) && treat.name %in% names(covs)) stop("The variable on the left side of the formula appears on the right side too.", call. = FALSE)
    }

    if (eval.model.matrx) {
        if (s <- !identical(sep, "")) {
            if (!is.character(sep) || length(sep) > 1) stop("'sep' must be a string of length 1.", call. = FALSE)
            original.covs.levels <- make_list(names(covs))
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
    }
    else {
        covs.matrix <- NULL
    }

    if (!terms) attr(covs, "terms") <- NULL

    return(list(reported.covs = covs,
                model.covs = covs.matrix,
                treat = treat,
                treat.name = treat.name))
}
assign.treat.type <- function(treat, use.multi = FALSE) {
    #Returns treat with treat.type attribute
    nunique.treat <- nunique(treat)

    if (nunique.treat < 2) {
        stop("The treatment must have at least two unique values.", call. = FALSE)
    }
    else if (!use.multi && nunique.treat == 2) {
        treat.type <- "binary"
    }
    else if (use.multi || is_(treat, c("factor", "character"))) {
        treat.type <- "multinomial"
        if (!is_(treat, "processed.treat")) treat <- factor(treat)
    }
    else {
        treat.type <- "continuous"
    }
    attr(treat, "treat.type") <- treat.type
    return(treat)
}
get.treat.type <- function(treat) {
    return(attr(treat, "treat.type"))
}
has.treat.type <- function(treat) {
    is_not_null(get.treat.type(treat))
}

#Input processing
process.bin.vars <- function(bin.vars, mat) {
    if (missing(bin.vars)) bin.vars <- is_binary_col(mat)
    else if (is_null(bin.vars)) bin.vars <- rep(FALSE, ncol(mat))
    else {
        if (is.logical(bin.vars)) {
            bin.vars[is.na(bin.vars)] <- FALSE
            if (length(bin.vars) != ncol(mat)) stop("If 'bin.vars' is logical, it must have length equal to the number of columns of 'mat'.")
        }
        else if (is.numeric(bin.vars)) {
            bin.vars <- bin.vars[!is.na(bin.vars) & bin.vars != 0]
            if (any(bin.vars < 0) && any(bin.vars > 0)) stop("Positive and negative indices cannot be mixed with 'bin.vars'.")
            if (any(abs(bin.vars) > ncol(mat))) stop("If 'bin.vars' is numeric, none of its values can exceed the number of columns of 'mat'.")
            logical.bin.vars <- rep(any(bin.vars < 0), ncol(mat))
            logical.bin.vars[abs(bin.vars)] <- !logical.bin.vars[abs(bin.vars)]
            bin.vars <- logical.bin.vars
        }
        else if (is.character(bin.vars)) {
            bin.vars <- bin.vars[!is.na(bin.vars) & bin.vars != ""]
            if (is_null(colnames(mat))) stop("If 'bin.vars' is character, 'mat' must have column names.")
            if (any(bin.vars %nin% colnames(mat))) stop("If 'bin.vars' is character, all its values must be column names of 'mat'.")
            bin.vars <- colnames(mat) %in% bin.vars
        }
        else stop("'bin.vars' must be a logical, numeric, or character vector.")
    }
    return(bin.vars)
}
process.s.weights <- function(s.weights, data = NULL) {
    #Process s.weights
    if (is_not_null(s.weights)) {
        if (!(is.character(s.weights) && length(s.weights) == 1) && !is.numeric(s.weights)) {
            stop("The argument to 's.weights' must be a vector or data frame of sampling weights or the (quoted) names of the variable in 'data' that contains sampling weights.", call. = FALSE)
        }
        if (is.character(s.weights) && length(s.weights)==1) {
            if (is_null(data)) {
                stop("'s.weights' was specified as a string but there was no argument to 'data'.", call. = FALSE)
            }
            else if (s.weights %in% names(data)) {
                s.weights <- data[[s.weights]]
            }
            else stop("The name supplied to 's.weights' is not the name of a variable in 'data'.", call. = FALSE)
        }
    }
    else s.weights <- NULL
    return(s.weights)
}

#Uniqueness
nunique <- function(x, nmax = NA, na.rm = TRUE) {
    if (is_null(x)) return(0)
    else {
        if (na.rm && anyNA(x)) x <- na.rem(x)
        if (is.factor(x)) return(nlevels(x))
        else return(length(unique(x, nmax = nmax)))
    }

}
nunique.gt <- function(x, n, na.rm = TRUE) {
    if (missing(n)) stop("'n' must be supplied.")
    if (n < 0) stop("'n' must be non-negative.")
    if (is_null(x)) FALSE
    else {
        if (n == 1) !all_the_same(x, na.rm)
        else if (length(x) < 2000) nunique(x, na.rm = na.rm) > n
        else tryCatch(nunique(x, nmax = n, na.rm = na.rm) > n, error = function(e) TRUE)
    }
}
all_the_same <- function(x, na.rm = TRUE) {
    if (anyNA(x)) {
        x <- na.rem(x)
        if (!na.rm) return(is_null(x))
    }
    if (is.double(x)) check_if_zero(max(x) - min(x))
    else all(x == x[1])
}
is_binary <- function(x, na.rm = TRUE) {
    if (na.rm && anyNA(x)) x <- na.rem(x)
    !all_the_same(x) && all_the_same(x[x != x[1]])
}
is_binary_col <- function(dat, na.rm = TRUE) {
    if (length(dim(dat)) != 2) stop("is_binary_col cannot be used with objects that don't have 2 dimensions.")
    apply(dat, 2, is_binary)
}

#R Processing
make_list <- function(n) {
    if (length(n) == 1L && is.numeric(n)) {
        vector("list", as.integer(n))
    }
    else if (is_(n, "atomic")) {
        setNames(vector("list", length(n)), as.character(n))
    }
    else stop("'n' must be an integer(ish) scalar or an atomic variable.")
}
make_df <- function(ncol, nrow = 0, types = "numeric") {
    if (length(ncol) == 1L && is.numeric(ncol)) {
        col_names <- NULL
        ncol <- as.integer(ncol)
    }
    else if (is_(ncol, "atomic")) {
        col_names <- as.character(ncol)
        ncol <- length(ncol)
    }
    if (length(nrow) == 1L && is.numeric(nrow)) {
        row_names <- NULL
        nrow <- as.integer(nrow)
    }
    else if (is_(nrow, "atomic")) {
        row_names <- as.character(nrow)
        nrow <- length(nrow)
    }
    df <- as.data.frame.matrix(matrix(NA_real_, nrow = nrow, ncol = ncol))
    colnames(df) <- col_names
    rownames(df) <- row_names
    if (is_not_null(types)) {
        if (length(types) %nin% c(1, ncol)) stop("'types' must be equal to the number of columns.")
        if (any(types %nin% c("numeric", "integer", "logical", "character", NA))) {
            stop("'types' must be an acceptable type. For factors, use NA.")
        }
        if (length(types) == 1) types <- rep(types, ncol)
        for (i in seq_len(ncol)) if (!is.na(types)[i] && types[i] != "numeric") df[[i]] <- get(types[i])(nrow)
    }
    return(df)
}
ifelse_ <- function(...) {
    dotlen <- ...length()
    if (dotlen %% 2 == 0) stop("ifelse_ must have an odd number of arguments: pairs of test/yes, and one no.")
    out <- ...elt(dotlen)
    if (dotlen > 1) {
        if (!is_(out, "atomic")) stop("The last entry to ifelse_ must be atomic.")
        if (length(out) == 1) out <- rep(out, length(..1))
        n <- length(out)
        for (i in seq_len((dotlen - 1)/2)) {
            test <- ...elt(2*i - 1)
            yes <- ...elt(2*i)
            if (length(yes) == 1) yes <- rep(yes, n)
            if (length(yes) != n || length(test) != n) stop("All entries must have the same length.")
            if (!is.logical(test)) stop(paste("The", ordinal(2*i - 1), "entry to ifelse_ must be logical."))
            if (!is_(yes, "atomic")) stop(paste("The", ordinal(2*i), "entry to ifelse_ must be atomic."))
            pos <- which(test)
            out[pos] <- yes[pos]
        }
    }
    else {
        if (!is_(out, "atomic")) stop("The first entry to ifelse_ must be atomic.")
    }
    return(out)
}
is_ <- function(x, types, stop = FALSE, arg.to = FALSE) {
    s1 <- deparse1(substitute(x))
    if (is_not_null(x)) {
        for (i in types) {
            if (i == "list") it.is <- is.list(x) && !is.data.frame(x)
            else if (is_not_null(get0(paste0("is_", i)))) {
                it.is <- get0(paste0("is_", i))(x)
            }
            else if (is_not_null(get0(paste.("is", i)))) {
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
            stop(paste0(s0, "'", s1, "' must be a ", word_list(types, and.or = "or"), " ", s2, "."), call. = FALSE)
        }
    }
    return(it.is)
}
is_null <- function(x) length(x) == 0L
is_not_null <- function(x) !is_null(x)
if_null_then <- function(x1 = NULL, x2 = NULL, ...) {
    if (is_not_null(x1)) x1
    else if (is_not_null(x2)) x2
    else if (...length() > 0) {
        for (k in seq_len(...length())) {
            if (is_not_null(...elt(k))) return(...elt(k))
        }
        return(..1)
    }
    else return(x1)
}
clear_null <- function(x) {
    x[vapply(x, is_null, logical(1L))] <- NULL
    return(x)
}
clear_attr <- function(x, all = FALSE) {
    if (all) {
        attributes(x) <- NULL
    }
    else {
        dont_clear <- c("names", "class", "dim", "dimnames", "row.names")
        attributes(x)[names(attributes(x)) %nin% dont_clear] <- NULL
    }
    return(x)
}
probably.a.bug <- function() {
    fun <- paste(deparse1(sys.call(-1)), collapse = "\n")
    stop(paste0("An error was produced and is likely a bug. Please let the maintainer know a bug was produced by the function\n",
                fun), call. = FALSE)
}
`%nin%` <- function(x, table) is.na(match(x, table, nomatch = NA_integer_))
`%pin%` <- function(x, table) {
    #Partial in. TRUE if x uniquely identifies values in table.
    !is.na(pmatch(x, table))
}
`%cin%` <- function(x, table) {
    #Partial in w/ charmatch. TRUE if x at all in table.
    !is.na(charmatch(x, table))
}
null_or_error <- function(x) {is_null(x) || any(class(x) == "try-error")}
match_arg <- function(arg, choices, several.ok = FALSE) {
    #Replaces match.arg() but gives cleaner error message and processing
    #of arg.
    if (missing(arg))
        stop("No argument was supplied to match_arg.", call. = FALSE)
    arg.name <- deparse1(substitute(arg))

    if (missing(choices)) {
        formal.args <- formals(sys.function(sysP <- sys.parent()))
        choices <- eval(formal.args[[as.character(substitute(arg))]],
                        envir = sys.frame(sysP))
    }

    if (is.null(arg))
        return(choices[1L])
    else if (!is.character(arg))
        stop(paste0("The argument to '", arg.name, "' must be NULL or a character vector"), call. = FALSE)
    if (!several.ok) {
        if (identical(arg, choices))
            return(arg[1L])
        if (length(arg) > 1L)
            stop(paste0("The argument to '", arg.name, "' must be of length 1"), call. = FALSE)
    }
    else if (is_null(arg))
        stop(paste0("The argument to '", arg.name, "' must be of length >= 1"), call. = FALSE)

    i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
    if (all(i == 0L))
        stop(paste0("The argument to '", arg.name, "' should be ", if (length(choices) > 1) {if (several.ok) "at least one of " else "one of "} else "",
                    word_list(choices, and.or = "or", quotes = 2), "."),
             call. = FALSE)
    i <- i[i > 0L]
    if (!several.ok && length(i) > 1)
        stop("There is more than one match in 'match_arg'")
    choices[i]
}
last <- function(x) {
    x[[length(x)]]
}
`last<-` <- function(x, value) {
    `[[<-`(x, length(x), value)
}
len <- function(x, recursive = TRUE) {
    if (is_null(x)) 0L
    else if (length(dim(x)) > 1) NROW(x)
    else if (is.list(x) && recursive) vapply(x, len, numeric(1L), recursive = FALSE)
    else length(x)
}
na.rem <- function(x) {
    #A faster na.omit for vectors
    x[!is.na(x)]
}
check.package <- function(package.name, alternative = FALSE) {
    packages.not.installed <- package.name[!vapply(package.name, requireNamespace, logical(1L),
                                                   quietly = TRUE)]
    if (is_not_null(packages.not.installed)) {
        if (alternative) return(FALSE)
        else {
            plural <- length(packages.not.installed) > 1
            stop(paste0("Package", if (plural) "s " else " ",
                        word_list(packages.not.installed, quotes = 1, is.are = TRUE),
                        " needed for this function to work. Please install ",
                        if (plural) "them" else "it","."),
                 call. = FALSE)
        }
    }
    else return(invisible(TRUE))
}
check_if_call_from_fun <- function(fun) {
    # Check if called from within function f
    if (missing(fun) || !exists(deparse1(substitute(fun)), mode = "function")) return(FALSE)
    sp <- sys.parents()
    sys.funs <- lapply(sp, sys.function)
    for (x in sys.funs) {
        if (identical(fun, x)) return(TRUE)
    }
    FALSE
}

if (getRversion() < 3.6) str2expression <- function(text) parse(text=text, keep.source=FALSE)

#Not used cobalt; replaced with rlang
is.formula <- function(f, sides = NULL) {
    #Replaced by rlang::is_formula
    res <- inherits(f, "formula") && is.name(f[[1]]) && deparse1(f[[1]]) %in% c( '~', '!') &&
        length(f) >= 2
    if (is_not_null(sides) && is.numeric(sides) && sides %in% c(1,2)) {
        res <- res && length(f) == sides + 1
    }
    return(res)
}
