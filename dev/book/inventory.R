args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
rows <- list()
anonymous <- 0L
for (file in list.files('R', pattern = '\\.R$', full.names = TRUE)) {
  p <- getParseData(parse(file, keep.source = TRUE))
  if (is.null(p)) next
  lines <- readLines(file, warn = FALSE)
  node <- function(id) p[match(id, p$id), , drop = FALSE]
  source_text <- function(id) {
    x <- node(id)
    s <- lines[x$line1:x$line2]
    s[length(s)] <- substr(s[length(s)], 1, x$col2)
    s[1] <- substring(s[1], x$col1)
    paste(s, collapse = '\n')
  }
  fs <- p[p$token == 'FUNCTION', , drop = FALSE]
  fn_ids <- fs$parent
  named <- list()
  for (i in seq_len(nrow(fs))) {
    fn <- node(fs$parent[i])
    assignment <- node(fn$parent)
    children <- p[p$parent == assignment$id, , drop = FALSE]
    ops <- children[children$token %in% c('LEFT_ASSIGN', 'EQ_ASSIGN', 'RIGHT_ASSIGN'), , drop = FALSE]
    name <- NULL
    if (nrow(ops) == 1L) {
      others <- children[children$token == 'expr' & children$id != fn$id, , drop = FALSE]
      if (nrow(others) == 1L) name <- source_text(others$id)
    }
    if (is.null(name)) { anonymous <- anonymous + 1L; next }
    name <- gsub('^`|`$', '', name)
    named[[as.character(fn$id)]] <- name
    parents <- character()
    ancestor <- fn$parent
    while (!is.na(ancestor) && ancestor > 0L) {
      if (ancestor %in% fn_ids) parents <- c(as.character(ancestor), parents)
      ancestor <- node(ancestor)$parent
    }
    rows[[length(rows) + 1L]] <- list(name = name, file = file,
      line = fn$line1, function_id = fn$id, parent_ids = parents,
      nested = length(parents) > 0L)
  }
  for (j in seq_along(rows)) if (rows[[j]]$file == file) {
    parents <- vapply(rows[[j]]$parent_ids, function(id) {
      if (is.null(named[[id]])) '<callback>' else named[[id]]
    }, character(1))
    rows[[j]]$scope <- paste(parents, collapse = ' / ')
  }
}
jsonlite::write_json(list(functions = rows, anonymous_callbacks = anonymous),
 args[[1]], auto_unbox = TRUE, pretty = TRUE)
