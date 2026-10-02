# test-adopt-aux.R
#
# Round-9 adoption of the audit's "aux-scripts" group (group-reports/
# aux-scripts.md, findings AUX-F1..F20 and the test-harness review, section 3):
# offline tests for the build / development scripts that were fixed in round 1
# and had no regression test, apart from dev/test-attribution-guard.sh itself.
#
# Everything that reads repository files (dev/, evals/, .githooks/, .github/,
# tools/, benchmarks/, data-raw/, docs/) calls skip_if_no_source() first: none of
# it exists in an installed package. Shell tests skip on Windows and when the
# interpreter (bash, python3, git, jq) is missing. All writes go to tempdir();
# the scripts under test are run either against a throw-away git repository
# built in tempdir() (so no branch/worktree/log of the real checkout is touched)
# or with their temp/audit locations redirected into tempdir().

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
.aux_path <- function(...) {
  normalizePath(testthat::test_path("..", "..", ...), mustWork = FALSE)
}

.aux_need_tools <- function(...) {
  testthat::skip_on_os("windows")
  for (tl in c(...)) {
    testthat::skip_if(!nzchar(Sys.which(tl)), paste(tl, "is not available"))
  }
}

# run a command; returns list(out = <stdout + stderr lines>, status = <int>).
# Output goes through a file (not a pipe) so that non-ASCII text from the scripts
# (arrows, check marks) cannot trip R's re-encoding in a C locale.
.aux_run <- function(cmd, args = character(), env = character(), dir = NULL,
                     input = NULL) {
  out_f <- tempfile("auxout_")
  on.exit(unlink(out_f), add = TRUE)
  wrapped <- c("-c", "\"$@\" > \"$AUX_OUT\" 2>&1", "aux", cmd, args)
  run <- function() {
    suppressWarnings(system2("bash", shQuote(wrapped), input = input,
                             env = c(env, paste0("AUX_OUT=", shQuote(out_f)))))
  }
  st <- if (is.null(dir)) run() else withr::with_dir(dir, run())
  out <- if (file.exists(out_f)) .aux_lines(out_f) else character()
  list(out = out, status = as.integer(st))
}

.aux_lines <- function(f) suppressWarnings(readLines(f, warn = FALSE, encoding = "UTF-8"))

.aux_env <- function(...) {
  kv <- c(...)
  paste0(names(kv), "=", shQuote(unname(kv)))
}

.aux_git <- function(dir, ...) {
  .aux_run("git", c("-C", dir, "-c", "user.name=aux", "-c", "user.email=aux@example.org",
                    "-c", "core.hooksPath=/dev/null", "-c", "commit.gpgsign=false", ...))
}

# a throw-away git repository holding copies of the scripts under test
.aux_cache <- new.env()
.aux_fixture <- function() {
  if (!is.null(.aux_cache$dir)) return(.aux_cache$dir)
  .aux_need_tools("bash", "git", "jq", "python3")
  skip_if_no_source("evals", "run.sh")
  skip_if_no_source("dev", "dual.sh")
  root <- .aux_path()
  d <- tempfile("auxrepo_")
  dir.create(d)
  files <- c("evals/run.sh", "evals/mutations.json", "dev/dual.sh", "dev/debate.sh",
             "dev/audit-all.sh", "dev/lib/audit.sh", "dev/lib/verdict.sh",
             "dev/lib/worktree.sh", "dev/lib/parse_verdict.py",
             ".githooks/commit-msg", ".githooks/ai-patterns", "R/select_ocs.R",
             "R/select_ind.R", "R/select_usefulness.R", "docs/THEORY_REVIEW.md")
  for (f in files) {
    dest <- file.path(d, f)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(root, f), dest, copy.mode = TRUE)
  }
  testthat::skip_if(.aux_git(d, "init", "-q")$status != 0, "git init failed")
  .aux_git(d, "add", "-A")
  r <- .aux_git(d, "commit", "-q", "-m", "fixture")
  testthat::skip_if(r$status != 0, "git commit failed in the fixture repository")
  .aux_cache$dir <- d
  d
}

.aux_fake <- function(dir, name, body) {
  f <- file.path(dir, name)
  writeLines(c("#!/bin/bash", body), f)
  Sys.chmod(f, "0755")
  f
}

.aux_agree <- c("echo 'looks fine'", "echo '```json'",
                "echo '{\"verdict\":\"AGREE\",\"open\":[],\"confidence\":0.5,\"summary\":\"ok\"}'",
                "echo '```'")
.aux_block <- c("echo 'objection at O1 and S1'", "echo '```json'",
                "echo '{\"verdict\":\"BLOCK\",\"open\":[\"O1\"],\"confidence\":0.9,\"summary\":\"bad\"}'",
                "echo '```'")

# ---------------------------------------------------------------------------
# AUX-F4: attribution guard (hook, shared pattern file, CI)
# ---------------------------------------------------------------------------
test_that("AUX-F4: the guard's own regression script passes (>= 40 cases)", {
  .aux_need_tools("bash")
  skip_if_no_source("dev", "test-attribution-guard.sh")
  skip_if_no_source(".githooks", "commit-msg")
  tmp <- tempfile("guard_")
  dir.create(tmp)
  r <- .aux_run("bash", .aux_path("dev", "test-attribution-guard.sh"),
                env = .aux_env(TMPDIR = tmp))
  expect_identical(r$status, 0L)
  last <- tail(r$out, 1)
  expect_match(last, "^attribution guard: [0-9]+/[0-9]+ cases as expected$")
  nums <- as.integer(strsplit(sub("^attribution guard: ([0-9]+/[0-9]+).*$", "\\1", last), "/")[[1]])
  expect_identical(nums[1], nums[2])
  expect_gte(nums[2], 40L)
})

test_that("AUX-F4: ai_attribution_scan blocks further assistant names and signatures and passes human text", {
  .aux_need_tools("bash")
  skip_if_no_source(".githooks", "ai-patterns")
  pat <- .aux_path(".githooks", "ai-patterns")
  scan <- function(msg) {
    .aux_run("bash", c("-c", ". \"$1\"; printf '%b' \"$2\" | ai_attribution_scan", "scan", pat, msg))
  }
  blocked <- c(
    "fix: x\\n\\nCo-authored-by: Windsurf <w@example.org>\\n",
    "fix: x\\n\\nGenerated with Gemini CLI\\n",
    "fix: x\\n\\nCreated by Cursor\\n",
    "fix: x\\n\\nCo-authored-by: aider (gpt-4o) <noreply@aider.chat>\\n",
    "fix: x\\n\\nAssisted-by: GitHub Copilot\\n",
    "fix: x\\n\\nbody text\\n\\nAuthored with ChatGPT\\n",
    "fix: x\\n\\nCo-authored-by: Claude Code <claude@anthropic.com>\\n",
    "fix: x\\n\\nPowered by Codex\\n",
    "fix: x\\n\\nco-authored-by: Codeium <c@example.org>\\n"
  )
  for (m in blocked) {
    r <- scan(m)
    expect_identical(r$status, 1L, info = m)
    expect_true(any(grepl("offending|attribution|trailer", r$out, ignore.case = TRUE)), info = m)
  }
  clean <- c(
    "fix: x\\n\\nCo-authored-by: Samuel Fernandes <samuelbf@uark.edu>\\n",
    "docs: reviewed by the maintainers\\n",
    "feat: add maintain() helper\\n",
    "fix: tables made from the real data by hand\\n",
    "fix: x\\n\\nReviewed-by: Alexander Lipka <alipka@illinois.edu>\\n",
    "chore: bump\\n\\nSigned-off-by: Samuel Fernandes <samuelbf@uark.edu>\\n"
  )
  for (m in clean) {
    expect_identical(scan(m)$status, 0L, info = m)
  }
})

test_that("AUX-F4: the commit-msg hook rejects with the reason, ignores '#' lines, and also works on /dev/stdin (dev/setup.sh self-test)", {
  .aux_need_tools("bash")
  skip_if_no_source(".githooks", "commit-msg")
  hook <- .aux_path(".githooks", "commit-msg")
  expect_true(file.access(hook, 1L) == 0L)                     # executable, or git cannot run it
  f <- tempfile()
  withr::defer(unlink(f))
  writeLines(c("feat: x", "", "Co-Authored-By: Claude <noreply@anthropic.com>"), f)
  r <- .aux_run("bash", c(hook, f))
  expect_identical(r$status, 1L)
  expect_true(any(grepl("commit blocked", r$out)))
  expect_true(any(grepl("offending", r$out)))
  expect_true(any(grepl("rule #1", r$out)))
  writeLines(c("feat: x", "# Co-Authored-By: Claude <noreply@anthropic.com>", "# \U0001F916"), f)
  expect_identical(.aux_run("bash", c(hook, f))$status, 0L)
  r <- .aux_run("bash", c("-c", "printf 'test\\n\\nCo-Authored-By: Claude <noreply@anthropic.com>\\n' | \"$1\" /dev/stdin", "x", hook))
  expect_identical(r$status, 1L)
})

test_that("AUX-F4: local hook and CI workflow share ONE pattern file (no drift), and CI runs the guard self-test", {
  skip_if_no_source(".githooks", "commit-msg")
  skip_if_no_source(".github", "workflows", "attribution-guard.yml")
  hook <- paste(.aux_lines(.aux_path(".githooks", "commit-msg")), collapse = "\n")
  ci <- paste(.aux_lines(.aux_path(".github", "workflows", "attribution-guard.yml")),
              collapse = "\n")
  expect_match(hook, "ai-patterns", fixed = TRUE)
  expect_match(ci, ". .githooks/ai-patterns", fixed = TRUE)
  expect_match(ci, "dev/test-attribution-guard.sh", fixed = TRUE)
  # neither defines its own pattern list any more
  expect_false(grepl("AI_NAMES=", hook, fixed = TRUE))
  expect_false(grepl("AI_NAMES=", ci, fixed = TRUE))
  # the over-broad term of the old hook (any GitHub noreply address) is gone
  pats <- paste(.aux_lines(.aux_path(".githooks", "ai-patterns")), collapse = "\n")
  expect_false(grepl("AI_NAMES='[^']*users\\.noreply", pats))
})

test_that("AUX-F4: the CI scan step, run against commit ranges in a scratch repository, passes clean history and fails on an AI trailer", {
  .aux_need_tools("bash", "git")
  skip_if_no_source(".github", "workflows", "attribution-guard.yml")
  skip_if_no_source(".githooks", "ai-patterns")
  lines <- .aux_lines(.aux_path(".github", "workflows", "attribution-guard.yml"))
  i0 <- grep("- name: Scan commit messages", lines, fixed = TRUE)
  expect_length(i0, 1L)
  i1 <- i0 + which(grepl("^\\s+run: \\|\\s*$", lines[(i0 + 1):length(lines)]))[1]
  body <- character()
  for (k in (i1 + 1):length(lines)) {
    if (nzchar(trimws(lines[k])) && !grepl("^ {10}", lines[k])) break
    body <- c(body, sub("^ {10}", "", lines[k]))
  }
  script <- paste(body, collapse = "\n")
  expect_match(script, "ai_attribution_scan", fixed = TRUE)

  repo <- tempfile("cirepo_")
  dir.create(repo)
  dir.create(file.path(repo, ".githooks"))
  file.copy(.aux_path(".githooks", "ai-patterns"), file.path(repo, ".githooks", "ai-patterns"))
  skip_if(.aux_git(repo, "init", "-q")$status != 0)
  sha <- character()
  commit <- function(msg) {
    writeLines(msg, file.path(repo, "msg.txt"))
    writeLines(as.character(length(sha)), file.path(repo, "f.txt"))
    .aux_git(repo, "add", "-A")
    r <- .aux_git(repo, "commit", "-q", "-F", "msg.txt")
    skip_if(r$status != 0, "commit failed in the scratch repository")
    sha <<- c(sha, .aux_run("git", c("-C", repo, "rev-parse", "HEAD"))$out[1])
  }
  commit("feat: first")
  commit("fix: second")
  commit("fix: third\n\nCo-Authored-By: Claude <noreply@anthropic.com>")
  scan <- function(event, before = "", base = "", head = "", sha_now) {
    s <- script
    s <- gsub("${{ github.event_name }}", event, s, fixed = TRUE)
    s <- gsub("${{ github.event.before }}", before, s, fixed = TRUE)
    s <- gsub("${{ github.event.pull_request.base.sha }}", base, s, fixed = TRUE)
    s <- gsub("${{ github.event.pull_request.head.sha }}", head, s, fixed = TRUE)
    s <- gsub("${{ github.sha }}", sha_now, s, fixed = TRUE)
    sf <- tempfile("ciscan_", fileext = ".sh")
    writeLines(s, sf, useBytes = TRUE)
    .aux_run("bash", sf, dir = repo)
  }
  ok <- scan("push", before = sha[1], sha_now = sha[2])
  expect_identical(ok$status, 0L)
  expect_true(any(grepl("no AI attribution found", ok$out)))
  bad <- scan("push", before = sha[1], sha_now = sha[3])
  expect_identical(bad$status, 1L)
  expect_true(any(grepl("::error::", bad$out, fixed = TRUE)))
  # pull request: base..head
  expect_identical(scan("pull_request", base = sha[1], head = sha[2], sha_now = sha[2])$status, 0L)
  expect_identical(scan("pull_request", base = sha[1], head = sha[3], sha_now = sha[3])$status, 1L)
  # first push of a branch (all-zero "before"): only the head commit is scanned
  zero <- strrep("0", 40)
  expect_identical(scan("push", before = zero, sha_now = sha[2])$status, 0L)
  expect_identical(scan("push", before = zero, sha_now = sha[3])$status, 1L)
})

test_that("AUX-F4: git itself refuses a commit through the hook (end to end), and accepts a human one", {
  .aux_need_tools("bash", "git")
  skip_if_no_source(".githooks", "commit-msg")
  skip_if_no_source(".githooks", "ai-patterns")
  repo <- tempfile("hookrepo_")
  dir.create(file.path(repo, ".githooks"), recursive = TRUE)
  for (f in c("commit-msg", "ai-patterns")) {
    file.copy(.aux_path(".githooks", f), file.path(repo, ".githooks", f), copy.mode = TRUE)
  }
  skip_if(.aux_git(repo, "init", "-q")$status != 0)
  commit <- function(msg) {
    writeLines(msg, file.path(repo, "msg.txt"))
    .run <- function(...) {
      .aux_run("git", c("-C", repo, "-c", "user.name=aux", "-c", "user.email=aux@example.org",
                        "-c", "core.hooksPath=.githooks", "-c", "commit.gpgsign=false", ...))
    }
    .run("commit", "-q", "--allow-empty", "-F", "msg.txt")
  }
  expect_gt(commit("feat: x\n\nCo-Authored-By: Claude <noreply@anthropic.com>")$status, 0L)
  expect_identical(commit("feat: x\n\nCo-authored-by: Jane Doe <12345+jd@users.noreply.github.com>")$status, 0L)
})

# ---------------------------------------------------------------------------
# AUX-F1 / F5 / F6: evals (golden set + runner)
# ---------------------------------------------------------------------------
test_that("evals/mutations.json: every seeded bug applies exactly once to the current source, changes it, and still parses", {
  skip_if_no_source("evals", "mutations.json")
  skip_if_no_source("docs", "THEORY_REVIEW.md")
  fromJSON <- optional_fun("jsonlite", "fromJSON")
  m <- fromJSON(.aux_path("evals", "mutations.json"), simplifyDataFrame = FALSE)
  expect_gte(length(m), 5L)
  ids <- vapply(m, function(x) x$id, "")
  expect_identical(anyDuplicated(ids), 0L)
  doc <- paste(.aux_lines(.aux_path("docs", "THEORY_REVIEW.md")), collapse = "\n")
  for (x in m) {
    expect_true(all(c("id", "file", "rubric", "find", "replace", "bug", "expect") %in% names(x)),
                info = x$id)
    expect_identical(x$expect, "BLOCK")
    expect_false(identical(x$find, x$replace), info = x$id)
    f <- .aux_path(x$file)
    skip_if_not(file.exists(f), paste(x$file, "not present"))
    src <- paste(.aux_lines(f), collapse = "\n")
    hits <- gregexpr(x$find, src, fixed = TRUE)[[1]]
    expect_identical(sum(hits > 0), 1L, info = paste(x$id, "find string must occur exactly once"))
    mutated <- sub(x$find, x$replace, src, fixed = TRUE)
    expect_false(identical(src, mutated), info = x$id)
    expect_no_error(parse(text = mutated))                       # a logic bug, not a syntax error
    expect_true(grepl(paste0("\\*\\*", x$rubric, "\\*\\*"), doc), info = paste("rubric", x$rubric))
  }
})

test_that("AUX-F1/F5: evals/run.sh --check validates the golden set without writing to any source file", {
  d <- .aux_fixture()
  tmp <- tempfile("evtmp_")
  dir.create(tmp)
  targets <- c("R/select_ocs.R", "R/select_ind.R", "R/select_usefulness.R")
  before <- tools::md5sum(file.path(d, targets))
  r <- .aux_run("bash", "evals/run.sh", dir = d, env = .aux_env(PIPELINE_TMPDIR = tmp))
  expect_identical(r$status, 0L)
  expect_true(any(grepl("check: 5/5 present mutations valid", r$out)))
  expect_identical(unname(tools::md5sum(file.path(d, targets))), unname(before))
  expect_identical(length(list.files(tmp)), 0L)                  # its own mktemp file is removed
  # --check reports a drifted target (find string no longer present) with a non-zero exit
  f <- file.path(d, "R/select_ocs.R")
  orig <- readBin(f, "raw", file.size(f))
  withr::defer(writeBin(orig, f))
  txt <- .aux_lines(f)
  writeLines(sub("2 * sum(p * (1 - p))", "sum(p * (1 - p)) * 2", txt, fixed = TRUE), f)
  r2 <- .aux_run("bash", "evals/run.sh", dir = d, env = .aux_env(PIPELINE_TMPDIR = tmp))
  expect_gt(r2$status, 0L)
  expect_true(any(grepl("find string not present", r2$out)))
})

test_that("AUX-F1: evals/run.sh --eval never leaves a seeded bug behind, whatever the reviewer does (crash / AGREE / BLOCK)", {
  d <- .aux_fixture()
  bin <- tempfile("evbin_")
  dir.create(bin)
  tmp <- tempfile("evtmp_")
  dir.create(tmp)
  fake_fail <- .aux_fake(bin, "fail.sh", c("echo 'reviewer crashed' >&2", "exit 1"))
  fake_agree <- .aux_fake(bin, "agree.sh", .aux_agree)
  fake_block <- .aux_fake(bin, "block.sh", .aux_block)
  targets <- file.path(d, c("R/select_ocs.R", "R/select_ind.R", "R/select_usefulness.R"))
  before <- tools::md5sum(targets)
  eval_run <- function(fake) {
    .aux_run("bash", c("evals/run.sh", "--eval"), dir = d,
             env = .aux_env(PIPELINE_TMPDIR = tmp, REVIEWER_CMD = fake))
  }
  # 1. the reviewer crashes on every call
  r <- eval_run(fake_fail)
  expect_gt(r$status, 0L)
  expect_true(any(grepl("reviewer errors (non-zero exit): 6", r$out, fixed = TRUE)))
  # 2. the reviewer approves everything: recall 0/5, but the clean control passes
  r <- eval_run(fake_agree)
  expect_gt(r$status, 0L)
  expect_true(any(grepl("recall   (bugs blocked):     0/5 present", r$out, fixed = TRUE)))
  expect_true(any(grepl("precision(clean not blocked): PASS (AGREE)", r$out, fixed = TRUE)))
  # 3. the reviewer blocks everything: recall 5/5, but the clean control is blocked
  r <- eval_run(fake_block)
  expect_gt(r$status, 0L)
  expect_true(any(grepl("recall   (bugs blocked):     5/5 present", r$out, fixed = TRUE)))
  expect_true(any(grepl("precision(clean not blocked): SUSPECT (BLOCK)", r$out, fixed = TRUE)))
  # in every case: working-tree sources untouched, worktree and branch removed
  expect_identical(unname(tools::md5sum(targets)), unname(before))
  wt <- .aux_run("git", c("-C", d, "worktree", "list"))
  expect_length(wt$out, 1L)
  expect_identical(length(grep("^agent/", .aux_run("git", c("-C", d, "branch", "--format=%(refname:short)"))$out)), 0L)
  expect_identical(.aux_git(d, "status", "--porcelain", "--", "R", "evals", "docs", ".githooks",
                            "dev/lib", "dev/dual.sh")$out, character(0))     # no source changed
})

# ---------------------------------------------------------------------------
# dev/lib/parse_verdict.py: the contract dual.sh / debate.sh / audit-all.sh rely on
# ---------------------------------------------------------------------------
test_that("AUX-INFO-1: parse_verdict.py contract (verdict extraction, last JSON wins, template echo, --json, --field)", {
  .aux_need_tools("python3")
  skip_if_no_source("dev", "lib", "parse_verdict.py")
  py <- .aux_path("dev", "lib", "parse_verdict.py")
  pv <- function(text, ...) .aux_run("python3", c(py, ...), input = text)$out
  fence <- function(j) c("analysis", "```json", j, "```")
  expect_identical(pv(fence('{"verdict":"AGREE","open":[],"confidence":0.9,"summary":"ok"}')), "AGREE")
  expect_identical(pv(fence('{"verdict":"agree","open":[]}')), "AGREE")                  # case-insensitive
  expect_identical(pv(c(fence('{"verdict":"AGREE","open":[]}'), fence('{"verdict":"BLOCK","open":["O1"]}'))),
                   "BLOCK")                                                              # last one wins
  expect_identical(pv('{"verdict":"AGREE|BLOCK","open":["O1","O2"],"confidence":0.0}'), "UNKNOWN")  # echoed template
  expect_identical(pv("no structured verdict at all"), "UNKNOWN")
  expect_identical(pv('{"verdict":"MAYBE"}'), "UNKNOWN")
  expect_identical(pv("Overall, Verdict: AGREE"), "AGREE")                                 # documented prose fallback
  expect_identical(pv('{"verdict":"BLOCK","meta":{},"summary":"x"}'), "BLOCK")
  j <- jsonlite_or_skip <- optional_fun("jsonlite", "fromJSON")
  obj <- j(paste(pv(fence('{"verdict":"BLOCK","open":["O1"],"summary":"a|b"}'), "--json"), collapse = ""))
  expect_identical(obj$summary, "a|b")                                                    # '|' preserved (audit-all strips it)
  expect_identical(pv(fence('{"verdict":"BLOCK","open":["O1","O2"]}'), "--field", "open"), "O1,O2")
  expect_identical(pv("nothing", "--json"), "{}")
})

# ---------------------------------------------------------------------------
# AUX-F14 / F15 / F16: dual.sh, debate.sh, audit-all.sh
# ---------------------------------------------------------------------------
test_that("AUX-F14/F16: dual.sh review survives a missing path, sends the file contents, writes ONE transcript and ONE audit record", {
  d <- .aux_fixture()
  bin <- tempfile("dualbin_")
  dir.create(bin)
  prompt_file <- file.path(bin, "prompt.txt")
  # the audit record also runs `codex --version`: only capture the review call
  .aux_fake(bin, "codex", c("if [ \"$1\" = exec ]; then",
                           sprintf("  printf '%%s' \"$5\" > '%s'", prompt_file),
                           .aux_agree, "else echo 'fake-codex 1.0'; fi"))
  audit <- tempfile("audit_")
  tmp <- tempfile("dualtmp_")
  dir.create(tmp)
  transcript <- file.path(tmp, "one-transcript.md")
  r <- .aux_run("bash", c("dev/dual.sh", "review", "R/does_not_exist.R", "R/select_ind.R"), dir = d,
                env = .aux_env(PATH = paste(bin, Sys.getenv("PATH"), sep = ":"),
                               PIPELINE_TMPDIR = tmp, AUDIT_DIR = audit,
                               REVIEW_TRANSCRIPT = transcript, REVIEW_LABEL = "audit:test"))
  expect_identical(r$status, 0L)                                 # was a silent exit 1 under set -e
  expect_true(any(grepl("path not found (skipped in contents): R/does_not_exist.R", r$out, fixed = TRUE)))
  expect_true(any(grepl("\"verdict\":\"AGREE\"", r$out, fixed = TRUE)))
  prompt <- paste(.aux_lines(prompt_file), collapse = "\n")
  expect_match(prompt, "Scope: R/does_not_exist.R R/select_ind.R", fixed = TRUE)
  expect_match(prompt, "current contents", fixed = TRUE)
  expect_match(prompt, "intensity", fixed = TRUE)                # the reviewed file's text reached the prompt
  expect_true(file.exists(transcript))
  log <- .aux_lines(file.path(audit, "log.jsonl"))
  expect_length(log, 1L)
  expect_match(log, "\"verdict\":\"AGREE\"", fixed = TRUE)
  expect_match(log, "audit:test: R/does_not_exist.R", fixed = TRUE)
  expect_length(list.files(file.path(audit, "transcripts")), 0L)  # the single transcript is the one named
  # nothing in the scratch checkout changed
  expect_identical(.aux_git(d, "status", "--porcelain", "--", "R", "evals", "docs", ".githooks",
                            "dev/lib", "dev/dual.sh")$out, character(0))
})

test_that("AUX-F15: debate.sh ISOLATE=1 keeps its provenance in the main checkout, not in the throw-away worktree", {
  d <- .aux_fixture()
  bin <- tempfile("debbin_")
  dir.create(bin)
  .aux_fake(bin, "codex", .aux_agree)
  .aux_fake(bin, "claude", .aux_agree)
  .aux_fake(bin, "Rscript", "exit 0")                              # the objective gate (devtools::test) "passes"
  tmp <- tempfile("debtmp_")
  dir.create(tmp)
  r <- .aux_run("bash", c("dev/debate.sh", "R/select_ind.R"), dir = d,
                env = .aux_env(PATH = paste(bin, Sys.getenv("PATH"), sep = ":"),
                               PIPELINE_TMPDIR = tmp, ISOLATE = "1"))
  # clean the worktree/branch the script deliberately keeps for review
  wt <- grep("^worktree kept for review: ", r$out, value = TRUE)
  if (length(wt)) {
    path <- sub("^worktree kept for review: (.*) \\(branch (.*)\\); remove with.*$", "\\1", wt[1])
    br <- sub("^worktree kept for review: (.*) \\(branch (.*)\\); remove with.*$", "\\2", wt[1])
    .aux_run("git", c("-C", d, "worktree", "remove", "--force", path))
    .aux_run("git", c("-C", d, "branch", "-D", br))
  }
  expect_identical(r$status, 0L)
  expect_true(any(grepl("CONSENSUS", r$out)))
  log <- file.path(d, "dev", ".audit", "log.jsonl")
  expect_true(file.exists(log))
  expect_match(paste(.aux_lines(log), collapse = ""), "\"kind\":\"debate\"", fixed = TRUE)
  expect_match(paste(.aux_lines(log), collapse = ""), "\"verdict\":\"AGREE\"", fixed = TRUE)
})

test_that("AUX-F10/F16: dev/audit-all.sh --list resolves every group to existing files (hash.rs audited) without model calls", {
  .aux_need_tools("bash", "git", "jq")
  skip_if_no_source("dev", "audit-all.sh")
  skip_if_no_source("R", "legacy_vQTL.R")
  root <- .aux_path()
  skip_if(.aux_run("git", c("-C", root, "rev-parse", "--show-toplevel"))$status != 0,
          "not a git checkout")
  tmp <- tempfile("audtmp_")
  dir.create(tmp)
  r <- .aux_run("bash", c(.aux_path("dev", "audit-all.sh"), "--list"), dir = root,
                env = .aux_env(PIPELINE_TMPDIR = tmp, AUDIT_DIR = file.path(tmp, "audit")))
  expect_identical(r$status, 0L)
  grp <- grep("^  [a-z-]+ +[0-9]+ files: ", r$out, value = TRUE)
  expect_gte(length(grp), 8L)
  expect_false(any(grepl("no files present", r$out)))
  for (g in grp) {
    n <- as.integer(sub("^  [a-z-]+ +([0-9]+) files: .*$", "\\1", g))
    files <- strsplit(sub("^  [a-z-]+ +[0-9]+ files: ", "", g), " ")[[1]]
    expect_length(files, n)
    expect_true(all(file.exists(file.path(root, files))), info = g)
  }
  rust <- grep("^  rust-core", grp, value = TRUE)
  expect_length(rust, 1L)
  expect_match(rust, "src/rust/src/hash.rs", fixed = TRUE)
  # nothing was written to the checkout's own audit folder by this call
  expect_true(dir.exists(file.path(tmp, "audit")))
})

test_that("shell scripts and hooks parse (bash -n / python ast), and the entry points are executable", {
  .aux_need_tools("bash", "python3")
  skip_if_no_source("dev", "dual.sh")
  root <- .aux_path()
  sh <- c(list.files(file.path(root, "dev"), "\\.sh$", full.names = TRUE),
          list.files(file.path(root, "dev", "lib"), "\\.sh$", full.names = TRUE),
          file.path(root, "evals", "run.sh"),
          file.path(root, ".githooks", c("commit-msg", "ai-patterns")))
  sh <- sh[file.exists(sh)]
  expect_gte(length(sh), 8L)
  for (f in sh) {
    expect_identical(.aux_run("bash", c("-n", f))$status, 0L, info = f)
  }
  for (f in file.path(root, c("dev/dual.sh", "dev/debate.sh", "dev/audit-all.sh", "dev/setup.sh",
                              "dev/test-attribution-guard.sh", "evals/run.sh", ".githooks/commit-msg"))) {
    if (file.exists(f)) expect_true(file.access(f, 1L) == 0L, info = paste("not executable:", f))
  }
  r <- .aux_run("python3", c("-c", "import ast,sys; ast.parse(open(sys.argv[1]).read())",
                             file.path(root, "dev", "lib", "parse_verdict.py")))
  expect_identical(r$status, 0L)
})

# ---------------------------------------------------------------------------
# tools/msrv.R (the Rust version gate run by configure)
# ---------------------------------------------------------------------------
test_that("tools/msrv.R: error branches, MSRV comparison and the Cargo.toml rust-version cross-check", {
  .aux_need_tools("bash")
  skip_if_no_source("tools", "msrv.R")
  rscript <- file.path(R.home("bin"), "Rscript")
  msrv <- .aux_path("tools", "msrv.R")
  bin <- tempfile("rustbin_")
  dir.create(bin)
  mk_rust <- function(rustc = "1.80.1", cargo = "1.80.0") {
    .aux_fake(bin, "rustc", sprintf("echo 'rustc %s (aaaaaaaaa 2024-08-01)'", rustc))
    .aux_fake(bin, "cargo", sprintf("echo 'cargo %s (bbbbbbbbb 2024-07-01)'", cargo))
  }
  run <- function(sysreq, cargo_toml = NULL, no_sysreq = FALSE) {
    w <- tempfile("msrvwd_")
    dir.create(w)
    desc <- c("Package: x", "Version: 1.0")
    if (!no_sysreq) desc <- c(desc, paste0("SystemRequirements: ", sysreq))
    writeLines(desc, file.path(w, "DESCRIPTION"))
    if (!is.null(cargo_toml)) {
      dir.create(file.path(w, "src", "rust"), recursive = TRUE)
      writeLines(cargo_toml, file.path(w, "src", "rust", "Cargo.toml"))
    }
    .aux_run(rscript, msrv, dir = w, env = .aux_env(PATH = paste(bin, "/usr/bin:/bin", sep = ":")))
  }
  mk_rust()
  # missing / incomplete SystemRequirements
  r <- run(NULL, no_sysreq = TRUE)
  expect_gt(r$status, 0L)
  expect_true(any(grepl("`SystemRequirements` not found", r$out, fixed = TRUE)))
  r <- run("rustc >= 1.65.0, xz")
  expect_gt(r$status, 0L)
  expect_true(any(grepl("Cargo", r$out)))
  r <- run("Cargo (Rust's package manager), xz")
  expect_gt(r$status, 0L)
  expect_true(any(grepl("rustc", r$out)))
  # MSRV satisfied / not satisfied by the installed toolchain
  r <- run("Cargo (Rust's package manager), rustc >= 1.65.0, xz")
  expect_identical(r$status, 0L)
  expect_true(any(grepl("Using cargo 1.80.0", r$out)))
  expect_true(any(grepl("Using rustc 1.80.1", r$out)))
  r <- run("Cargo (Rust's package manager), rustc >= 1.99.0, xz")
  expect_gt(r$status, 0L)
  expect_true(any(grepl("UNSUPPORTED RUST VERSION", r$out)))
  expect_true(any(grepl("Minimum supported Rust version is 1.99.0", r$out, fixed = TRUE)))
  # the crate's own rust-version is the effective minimum when it is larger
  toml <- c("[package]", "name = 'x'", "rust-version = '1.85'")
  r <- run("Cargo (Rust's package manager), rustc >= 1.65.0, xz", cargo_toml = toml)
  expect_gt(r$status, 0L)                                        # installed 1.80.1 < 1.85
  expect_true(any(grepl("Minimum supported Rust version is 1.85", r$out, fixed = TRUE)))
  mk_rust(rustc = "1.90.0")
  r <- run("Cargo (Rust's package manager), rustc >= 1.65.0, xz", cargo_toml = toml)
  expect_identical(r$status, 0L)
  expect_true(any(grepl("lists rustc >= 1.65.0 but src/rust/Cargo.toml requires rust-version 1.85", r$out, fixed = TRUE)))
  # a missing toolchain gives the install hint, not a stack trace
  unlink(file.path(bin, c("rustc", "cargo")))
  r <- run("Cargo (Rust's package manager), rustc >= 1.65.0, xz")
  expect_gt(r$status, 0L)
  expect_true(any(grepl("RUST NOT FOUND|CARGO NOT FOUND|rustc|cargo", r$out)))
})

test_that("tools/msrv.R: the real DESCRIPTION and Cargo.toml state the same minimum Rust version", {
  skip_if_no_source("DESCRIPTION")
  skip_if_no_source("src", "rust", "Cargo.toml")
  d <- read.dcf(.aux_path("DESCRIPTION"), fields = "SystemRequirements")[1, 1]
  ver <- sub(".*rustc >= ([0-9.]+).*", "\\1", d)
  expect_match(ver, "^[0-9]+\\.[0-9]+(\\.[0-9]+)?$")
  ct <- .aux_lines(.aux_path("src", "rust", "Cargo.toml"))
  rv <- sub("^[^=]*=\\s*['\"]?([0-9.]+).*$", "\\1", grep("^\\s*rust-version\\s*=", ct, value = TRUE)[1])
  expect_match(rv, "^[0-9]+\\.[0-9]+(\\.[0-9]+)?$")
  # DESCRIPTION must not understate what the crate needs (msrv.R would only note it at build time)
  expect_gte(utils::compareVersion(ver, rv), 0L)
})

# ---------------------------------------------------------------------------
# AUX-F9 / F17 / F16: package metadata drift
# ---------------------------------------------------------------------------
test_that("AUX-F9: CITATION.cff mirrors DESCRIPTION (release number and date)", {
  skip_if_no_source("CITATION.cff")
  skip_if_no_source("DESCRIPTION")
  desc <- read.dcf(.aux_path("DESCRIPTION"), fields = c("Version", "Date"))
  cff <- .aux_lines(.aux_path("CITATION.cff"))
  cff_ver <- sub("^version:\\s*['\"]?([^'\"]+)['\"]?\\s*$", "\\1", grep("^version:", cff, value = TRUE)[1])
  cff_date <- sub("^date-released:\\s*['\"]?([^'\"]+)['\"]?\\s*$", "\\1", grep("^date-released:", cff, value = TRUE)[1])
  # development suffixes (.9000+) advance between releases; the release number must agree
  base <- function(v) paste(strsplit(v, ".", fixed = TRUE)[[1]][1:3], collapse = ".")
  expect_identical(base(cff_ver), base(unname(desc[1, "Version"])))
  expect_identical(cff_date, unname(desc[1, "Date"]))
  expect_match(desc[1, "Date"], "^[0-9]{4}-[0-9]{2}-[0-9]{2}$")
})

test_that("AUX-F17: .Rbuildignore patterns are valid regular expressions", {
  skip_if_no_source(".Rbuildignore")
  pats <- .aux_lines(.aux_path(".Rbuildignore"))
  pats <- pats[nzchar(trimws(pats)) & !grepl("^\\s*#", pats)]
  expect_gt(length(pats), 20L)
  for (p in pats) {
    expect_no_error(grepl(p, "x", perl = TRUE))
  }
})

test_that("AUX-F17 (AUX-N1): .Rbuildignore lists every pattern once", {
  skip_if_no_source(".Rbuildignore")
  pats <- .aux_lines(.aux_path(".Rbuildignore"))
  pats <- pats[nzchar(trimws(pats)) & !grepl("^\\s*#", pats)]
  expect_identical(anyDuplicated(pats), 0L)
})

test_that("AUX-F17 (AUX-N2): every relative image embedded in README.md ships in the package", {
  skip_if_no_source("README.md")
  skip_if_no_source(".Rbuildignore")
  md <- paste(.aux_lines(.aux_path("README.md")), collapse = "\n")
  imgs <- unique(c(regmatches(md, gregexpr("<img src=\"[^\"]+\"", md))[[1]],
                   regmatches(md, gregexpr("!\\[[^]]*\\]\\([^)]+\\)", md))[[1]]))
  imgs <- sub("^<img src=\"([^\"]+)\"$", "\\1", imgs)
  imgs <- sub("^!\\[[^]]*\\]\\(([^)]+)\\)$", "\\1", imgs)
  imgs <- imgs[!grepl("^https?://", imgs)]
  expect_gt(length(imgs), 0L)
  pats <- .aux_lines(.aux_path(".Rbuildignore"))
  pats <- pats[nzchar(trimws(pats)) & !grepl("^\\s*#", pats)]
  excluded <- imgs[vapply(imgs, function(i) any(vapply(pats, function(p) grepl(p, i, perl = TRUE), NA)), NA)]
  for (i in imgs) expect_true(file.exists(.aux_path(i)), info = paste("missing:", i))
  expect_length(excluded, 0L)
})

test_that("AUX-F16: the CI workflows are well formed, the gates named in the docs are wired, and no branch filter names a stale branch", {
  skip_if_no_source(".github", "workflows", "R-CMD-check.yml")
  wf <- list.files(.aux_path(".github", "workflows"), "\\.ya?ml$", full.names = TRUE)
  expect_gte(length(wf), 3L)
  txt <- lapply(wf, function(f) paste(.aux_lines(f), collapse = "\n"))
  names(txt) <- basename(wf)
  for (nm in names(txt)) {
    expect_match(txt[[nm]], "^name: ", info = nm)
    expect_match(txt[[nm]], "\njobs:\n", fixed = TRUE, info = nm)
    expect_match(txt[[nm]], "runs-on: ", fixed = TRUE, info = nm)
    expect_false(grepl("\t", txt[[nm]], fixed = TRUE), info = paste(nm, "tab in YAML"))
  }
  expect_false(grepl("restructure", txt[["R-CMD-check.yml"]], fixed = TRUE))
  expect_match(txt[["eval-check.yml"]], "bash evals/run.sh --check", fixed = TRUE)
  expect_match(txt[["R-CMD-check.yml"]], "check-r-package", fixed = TRUE)
})

test_that("AUX-F16: every workflow parses as YAML with jobs that have steps", {
  skip_if_no_source(".github", "workflows", "R-CMD-check.yml")
  yaml_load <- optional_fun("yaml", "yaml.load")
  for (f in list.files(.aux_path(".github", "workflows"), "\\.ya?ml$", full.names = TRUE)) {
    y <- suppressWarnings(yaml_load(paste(.aux_lines(f), collapse = "\n")))
    expect_true(is.list(y$jobs) && length(y$jobs) >= 1L, info = basename(f))
    for (job in y$jobs) {
      expect_true(is.list(job$steps) && length(job$steps) >= 1L, info = basename(f))
      for (st in job$steps) expect_true(!is.null(st$uses) || !is.null(st$run), info = basename(f))
    }
  }
})

# ---------------------------------------------------------------------------
# AUX-F12 / F13 / benchmark + data-raw helpers
# ---------------------------------------------------------------------------
test_that("dev scripts, benchmarks, data-raw and fixture-capture scripts all parse", {
  skip_if_no_source("benchmarks")
  skip_if_no_source("data-raw")
  root <- .aux_path()
  scripts <- c(list.files(file.path(root, "benchmarks"), "\\.R$", full.names = TRUE),
               list.files(file.path(root, "data-raw"), "\\.R$", full.names = TRUE),
               list.files(file.path(root, "dev"), "\\.R$", full.names = TRUE),
               list.files(file.path(root, "tools"), "\\.R$", full.names = TRUE),
               list.files(file.path(root, "inst", "extdata"), "^capture.*\\.R$",
                          full.names = TRUE, recursive = TRUE))
  expect_gte(length(scripts), 10L)
  for (f in scripts) expect_no_error(parse(f, keep.source = FALSE))
})

test_that("AUX-F12/F13: benchmark scripts - no top-level on.exit(), no T/F, scratch script binds exactly formals(create_phenotypes)", {
  skip_if_no_source("benchmarks", "benchmark_as_numeric.R")
  skip_if_no_source("benchmarks", "scratch_create_phenotypes_args.R")
  toplevel_calls <- function(f) {
    ex <- parse(f, keep.source = FALSE)
    vapply(ex, function(e) if (is.call(e)) as.character(e[[1]])[1] else "", "")
  }
  expect_false("on.exit" %in% toplevel_calls(.aux_path("benchmarks", "benchmark_as_numeric.R")))
  for (f in list.files(.aux_path("benchmarks"), "\\.R$", full.names = TRUE)) {
    pd <- utils::getParseData(parse(f, keep.source = TRUE))
    tf <- pd[pd$token == "SYMBOL" & pd$text %in% c("T", "F"), ]
    expect_identical(nrow(tf), 0L, info = paste(basename(f), "uses T/F"))
  }
  f <- .aux_path("benchmarks", "scratch_create_phenotypes_args.R")
  ex <- parse(f, keep.source = FALSE)
  expect_false(any(vapply(ex, function(e) is.call(e) && identical(as.character(e[[1]])[1], "library") &&
                            identical(as.character(e[[2]]), "here"), NA)))
  assigned <- unlist(lapply(ex, function(e) {
    if (is.call(e) && as.character(e[[1]])[1] %in% c("=", "<-") && is.name(e[[2]])) as.character(e[[2]])
  }))
  expect_setequal(assigned, names(formals(create_phenotypes)))
})

test_that("data-raw/render_docs.R renders only sources that exist and writes the shipped copies", {
  skip_if_no_source("data-raw", "render_docs.R")
  src <- .aux_lines(.aux_path("data-raw", "render_docs.R"))
  code <- src[!grepl("^\\s*#", src)]
  inputs <- unlist(regmatches(code, gregexpr("\"[^\"]+\\.Rmd\"", code)))
  inputs <- gsub("\"", "", inputs)
  expect_gte(length(inputs), 3L)
  for (i in inputs) expect_true(file.exists(.aux_path(i)), info = i)
  outs <- gsub("\"", "", unlist(regmatches(code, gregexpr("output_file\\s*=\\s*\"[^\"]+\"", code))))
  outs <- sub("^output_file\\s*=\\s*", "", outs)
  dirs <- gsub("\"", "", sub("^output_dir\\s*=\\s*", "",
                             unlist(regmatches(code, gregexpr("output_dir\\s*=\\s*\"[^\"]+\"", code)))))
  expect_gte(length(outs), 3L)
  for (o in outs) {
    expect_true(file.exists(.aux_path(o)) || any(file.exists(.aux_path(dirs, o))), info = o)
  }
})

test_that("dev/accept-nonadditive-cor.R calls only internals that still exist and arguments simulate_phenotype() still has", {
  skip_if_no_source("dev", "accept-nonadditive-cor.R")
  code <- .aux_lines(.aux_path("dev", "accept-nonadditive-cor.R"))
  code <- code[!grepl("^\\s*#", code)]
  ns <- asNamespace("simplePHENOTYPES")
  internals <- unique(unlist(regmatches(code, gregexpr("\\.[a-z][A-Za-z0-9_]+(?=\\()", code, perl = TRUE))))
  internals <- grep("^\\.(genetic|component|[a-z]+_[a-z_]+)", internals, value = TRUE)
  for (nm in internals) expect_true(exists(nm, envir = ns, inherits = FALSE), info = nm)
  expect_true(all(c("cor") %in% names(formals(simulate_phenotype))) ||
                "..." %in% names(formals(simulate_phenotype)))
})

# ---------------------------------------------------------------------------
# AUX-F2: documented calls use the real API (including eval = FALSE chunks that
# nothing executes)
# ---------------------------------------------------------------------------
test_that("AUX-F2: every call to an exported function in the vignettes and README uses only arguments that function has", {
  skip_if_no_source("vignettes")
  skip_if_no_source("README.Rmd")
  rmd <- c(list.files(.aux_path("vignettes"), "\\.Rmd$", full.names = TRUE), .aux_path("README.Rmd"))
  expect_gte(length(rmd), 6L)
  exports <- getNamespaceExports("simplePHENOTYPES")
  ns <- asNamespace("simplePHENOTYPES")
  bad <- character()
  n_checked <- 0L
  for (f in rmd) {
    lines <- .aux_lines(f)
    open <- grep("^```+\\s*\\{r", lines)
    close <- grep("^```+\\s*$", lines)
    for (o in open) {
      cl <- close[close > o][1]
      if (is.na(cl) || cl <= o + 1L) next
      code <- lines[(o + 1L):(cl - 1L)]
      ex <- tryCatch(parse(text = code, keep.source = FALSE), error = function(e) NULL)
      if (is.null(ex)) next
      walk <- function(e) {
        if (!is.call(e)) return(invisible())
        h <- e[[1]]
        fn <- if (is.name(h)) as.character(h) else
          if (is.call(h) && length(h) == 3L && as.character(h[[1]]) %in% c("::", ":::") &&
              identical(as.character(h[[2]]), "simplePHENOTYPES")) as.character(h[[3]]) else NA_character_
        if (!is.na(fn) && fn %in% exports) {
          fo <- get(fn, envir = ns)
          if (is.function(fo)) {
            fm <- names(formals(fo))
            nm <- names(as.list(e))[-1]
            nm <- nm[!is.na(nm) & nzchar(nm)]
            n_checked <<- n_checked + 1L
            if (!("..." %in% fm)) {
              extra <- setdiff(nm, fm)
              # partial matching of argument names is legal R
              extra <- extra[!vapply(extra, function(a) sum(startsWith(fm, a)) == 1L, NA)]
              if (length(extra)) bad <<- c(bad, paste0(basename(f), ": ", fn, "(", paste(extra, collapse = ", "), ")"))
            }
          }
        }
        for (k in seq_along(e)) {
          if (identical(e[[k]], quote(expr = ))) next
          walk(e[[k]])
        }
      }
      for (e in ex) walk(e)
    }
  }
  expect_gt(n_checked, 50L)
  expect_identical(bad, character(0))
})

# ---------------------------------------------------------------------------
# AUX-F10: architecture document vs the Rust sources it lists
# ---------------------------------------------------------------------------
test_that("AUX-F10: every Rust source file named in docs/ARCHITECTURE.md exists, and every file in src/rust/src is named", {
  skip_if_no_source("docs", "ARCHITECTURE.md")
  skip_if_no_source("src", "rust", "src")
  doc <- .aux_lines(.aux_path("docs", "ARCHITECTURE.md"))
  named <- unique(unlist(regmatches(doc, gregexpr("\\b[a-z_]+\\.rs\\b", doc))))
  expect_gte(length(named), 4L)
  actual <- list.files(.aux_path("src", "rust", "src"), "\\.rs$")
  expect_setequal(named, actual)
})

# ---------------------------------------------------------------------------
# Documentation objects: roxygen output is complete and consistent
# ---------------------------------------------------------------------------
test_that("man/: every exported object is documented, usage matches the formals, no Rd error-level problems", {
  skip_if_no_source("man")
  skip_if_no_source("DESCRIPTION")
  root <- .aux_path()
  expect_length(unlist(tools::undoc(dir = root)), 0L)
  expect_length(unlist(tools::checkDocFiles(dir = root)), 0L)
  rd <- list.files(file.path(root, "man"), "\\.Rd$", full.names = TRUE)
  expect_gte(length(rd), 50L)
  for (f in rd) {
    msgs <- tools::checkRd(f)
    lev <- as.integer(sub("^checkRd: \\((-?[0-9]+)\\).*$", "\\1", as.character(msgs)))
    expect_true(all(is.na(lev) | lev >= -5L), info = paste(basename(f), paste(msgs, collapse = "; ")))
  }
})

# ---------------------------------------------------------------------------
# AUX-F8 + audit section 3: the harness runs under edition 3; the contract list
# in test-backend-contract.R (copied by hand) stays in step with its document
# ---------------------------------------------------------------------------
test_that("AUX-F8: DESCRIPTION declares testthat edition 3 and setup.R does not claim a different one", {
  skip_if_no_source("DESCRIPTION")
  d <- read.dcf(.aux_path("DESCRIPTION"), fields = "Config/testthat/edition")[1, 1]
  expect_identical(unname(d), "3")
  expect_identical(testthat::edition_get(), 3L)
  su <- .aux_lines(testthat::test_path("setup.R"))
  expect_false(any(grepl("edition 2|2e\\b", su)))
})

test_that("audit section 3: every function in test-backend-contract.R's hand-copied list is named in docs/BACKEND_CONTRACT.md", {
  skip_if_no_source("docs", "BACKEND_CONTRACT.md")
  tf <- testthat::test_path("test-backend-contract.R")
  skip_if_not(file.exists(tf))
  ex <- parse(tf, keep.source = FALSE)
  assign_call <- Filter(function(e) {
    is.call(e) && as.character(e[[1]])[1] == "<-" && identical(as.character(e[[2]]), "backend_contract")
  }, as.list(ex))
  expect_length(assign_call, 1L)
  listed <- eval(assign_call[[1]][[3]], baseenv())
  expect_gt(length(listed), 40L)
  doc <- paste(.aux_lines(.aux_path("docs", "BACKEND_CONTRACT.md")), collapse = "\n")
  absent <- listed[!vapply(listed, function(nm) grepl(paste0("`", nm, "`"), doc, fixed = TRUE) ||
                             grepl(paste0("`", nm, "("), doc, fixed = TRUE), NA)]
  expect_identical(absent, character(0))
  # ... and the list itself is exported (so the document cannot name a removed function)
  expect_identical(setdiff(listed, getNamespaceExports("simplePHENOTYPES")), character(0))
})

# ---------------------------------------------------------------------------
# Test-harness quality (audit section 3): no test block without an expectation
# ---------------------------------------------------------------------------
test_that("harness: every test_that() block asserts something (expect_*, skip, or a helper that does)", {
  dir <- testthat::test_path()
  files <- list.files(dir, "^(test|helper)-.*\\.R$", full.names = TRUE)
  expect_gte(length(files), 20L)
  parsed <- lapply(files, function(f) tryCatch(parse(f, keep.source = FALSE), error = function(e) NULL))
  names(parsed) <- basename(files)
  calls_in <- function(x) {
    out <- character()
    walk <- function(e) {
      if (is.call(e)) {
        h <- e[[1]]
        if (is.name(h)) out <<- c(out, as.character(h))
        else if (is.call(h) && length(h) == 3L && as.character(h[[1]]) %in% c("::", ":::")) {
          out <<- c(out, as.character(h[[3]]))
        }
        for (k in seq_along(e)) {
          if (identical(e[[k]], quote(expr = ))) next                 # empty argument, e.g. x[, 1]
          walk(e[[k]])
        }
      } else if (is.pairlist(e) || is.list(e)) {
        for (k in seq_along(e)) walk(e[[k]])
      }
    }
    walk(x)
    out
  }
  assertive <- function(nms) any(grepl("^(expect|skip|fail$|succeed$|stopifnot$|local_snapshot|verify_)", nms))
  # top-level helper functions that contain assertions (one level deep is enough here)
  helper_with_expect <- character()
  for (ex in parsed) {
    for (e in ex) {
      if (is.call(e) && as.character(e[[1]])[1] %in% c("<-", "=") && is.name(e[[2]]) &&
          is.call(e[[3]]) && as.character(e[[3]][[1]])[1] == "function") {
        if (assertive(calls_in(e[[3]]))) helper_with_expect <- c(helper_with_expect, as.character(e[[2]]))
      }
    }
  }
  empty <- character()
  for (nm in names(parsed)) {
    if (!startsWith(nm, "test-") || is.null(parsed[[nm]])) next
    for (e in parsed[[nm]]) {
      if (is.call(e) && identical(as.character(e[[1]])[1], "test_that")) {
        nms <- calls_in(e)
        if (!assertive(nms) && !any(nms %in% helper_with_expect)) {
          empty <- c(empty, paste0(nm, ": ", e[[2]]))
        }
      }
    }
  }
  expect_identical(empty, character(0), info = paste(empty, collapse = "\n"))
})
