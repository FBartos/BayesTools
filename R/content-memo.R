# Session memo of results that depend on the content of their inputs only.
#
# A validator or a metadata builder is a pure function of the objects it is
# given, and one workflow hands the same unmodified catalog, hypothesis AST,
# label parts, or scaling metadata to it again and again. The memo keeps a
# reference to each recent input and recognises a later input only when it is
# identical() to a kept one, so a modified copy never matches: it is computed
# (or validated) again, and an invalid one is rejected again. Nothing is
# recorded for a call that fails. Holding the reference also keeps the input
# from being modified in place. Like the parameter-map runtime cache the memo
# lives in this session only, never in the objects, so saved objects carry no
# marker that could outlive a change of the code.

# Entries kept per kind; a lookup compares against each of them.
.bt_content_memo_limit <- function(){
  16L
}

.bt_content_memo_same <- function(x, y){
  identical(x, y, num.eq = FALSE, single.NA = FALSE)
}

# The result of 'compute()' for the inputs 'key' (the object, or the list of
# objects, whose content decides the result): the kept result of an entry
# identical() to 'key', otherwise a fresh one that is kept. A recognised entry
# moves to the front so that frequently reused inputs stay. A 'keep' function
# decides from a fresh result whether it is kept at all, so that results that
# are cheap to recompute do not displace expensive ones.
.bt_content_memo <- function(kind, key, compute, limit = .bt_content_memo_limit(),
                             keep = NULL){

  entries <- .BayesTools_private$content_memo[[kind]]
  for(i in seq_along(entries)){
    if(.bt_content_memo_same(entries[[i]]$key, key)){
      if(i > 1L){
        .bt_content_memo_write(kind, c(entries[i], entries[-i]))
      }
      return(entries[[i]]$value)
    }
  }

  value <- compute()
  if(!is.null(keep) && !isTRUE(keep(value))){
    return(value)
  }
  entries <- c(list(list(key = key, value = value)), entries)
  if(length(entries) > limit){
    entries <- entries[seq_len(limit)]
  }
  .bt_content_memo_write(kind, entries)

  value
}

.bt_content_memo_write <- function(kind, entries){

  memo <- .BayesTools_private$content_memo
  if(is.null(memo)){
    memo <- list()
  }
  memo[[kind]] <- entries
  .BayesTools_private$content_memo <- memo

  invisible(NULL)
}

# Runs 'validate()' unless an object identical to 'key' already passed the
# validator of 'kind'. A validator that fails stops here and leaves the memo
# unchanged.
.bt_validate_once <- function(kind, key, validate){

  .bt_content_memo(kind, key, function(){
    validate()
    TRUE
  })

  invisible(TRUE)
}
