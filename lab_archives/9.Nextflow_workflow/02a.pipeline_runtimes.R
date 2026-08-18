#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# Nextflow workflow - step 02a
#
# Wall-clock time, CPU time and peak memory per process per run.
#
# No trace file is configured in nextflow.config - only `report` and `dag` are
# enabled - so there is no tab-separated trace to read. Nextflow's HTML
# execution report, however, embeds the whole per-task trace as JSON in a
# `window.data.trace = [...]` assignment, which is what this step parses. Every
# number here is Nextflow's own accounting, not a re-measurement.
#
# Fields used:
#   realtime  - wall clock the task itself ran, in ms
#   duration  - realtime plus scheduling/staging overhead, in ms
#   peak_rss  - high-water memory, in bytes
#   cpus      - allocated, from the process directive
#
# CPU time is realtime * cpus. It is allocated CPU time, not consumed - a
# process given 16 cores that uses one still bills 16.
#
# Run from the scripts/ directory - paths are relative to it.
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(jsonlite)
})

options(scipen = 999)

RES_NF  <- "../results_nf"
out_dir <- "../results/benchmarking/nextflow_workflow"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

reports <- list.files(RES_NF, pattern = "^execution_report_.*\\.html$", full.names = TRUE)
reports <- reports[!grepl("(^|/)\\._", reports)]

if (!length(reports)) stop("No execution_report_<run>.html found in ", RES_NF)

message("Execution reports: ", length(reports))

# ----------------------------
# Pull window.data.trace out of the HTML
# ----------------------------
# The assignment is a single very long line. Take everything between the first
# '[' after the assignment and its matching ']' by bracket counting, so a ']'
# inside a string value cannot truncate it.

extract_trace <- function(path) {

  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")

  # Anchor on the assignment `window.data = { "trace":[ ... ] }`. Note that
  # `window.data.trace` also appears further down as a DataTable reference -
  # matching that instead finds no array and silently yields nothing.
  m <- str_locate(txt, "window\\.data\\s*=\\s*\\{")
  if (any(is.na(m))) {
    warning("No window.data assignment in ", basename(path), call. = FALSE)
    return(NULL)
  }

  key <- str_locate(str_sub(txt, m[2]), '"trace"\\s*:\\s*\\[')
  if (any(is.na(key))) {
    warning("No trace array in ", basename(path), call. = FALSE)
    return(NULL)
  }

  start <- m[2] + key[2] - 1
  chars <- str_split(str_sub(txt, start), "")[[1]]

  depth <- 0; in_str <- FALSE; esc <- FALSE; end <- NA_integer_
  for (i in seq_along(chars)) {
    ch <- chars[i]
    if (esc) { esc <- FALSE; next }
    if (ch == "\\") { esc <- TRUE; next }
    if (ch == '"') { in_str <- !in_str; next }
    if (in_str) next
    if (ch == "[") depth <- depth + 1
    if (ch == "]") {
      depth <- depth - 1
      if (depth == 0) { end <- i; break }
    }
  }

  if (is.na(end)) {
    warning("Could not find the end of the trace array in ", basename(path), call. = FALSE)
    return(NULL)
  }

  json <- str_sub(txt, start, start + end - 1)

  # Nextflow embeds each task's shell script in the `script` field using
  # JavaScript string escaping, which allows \' - that is not legal JSON and
  # makes the whole array unparseable. Unescape it before handing to fromJSON.
  json <- str_replace_all(json, "\\\\'", "'")

  out <- try(fromJSON(json, flatten = TRUE), silent = TRUE)
  if (inherits(out, "try-error")) {
    warning("Could not parse the trace JSON in ", basename(path), ": ",
            conditionMessage(attr(out, "condition")), call. = FALSE)
    return(NULL)
  }

  run <- str_remove(str_remove(basename(path), "^execution_report_"), "\\.html$")
  as_tibble(out) %>% mutate(run = run)
}

trace <- map_dfr(reports, extract_trace)

if (!nrow(trace)) stop("No trace rows recovered from any execution report.")

message("Tasks recovered: ", nrow(trace))

num <- function(x) suppressWarnings(as.numeric(x))

# NOTE ON status: Nextflow reports CACHED for any task it resumed from a
# previous run's work directory rather than re-executing. CACHED is a success,
# not a failure - treating anything that is not COMPLETED as failed marks a
# resumed run as almost entirely broken. Failure is FAILED/ABORTED, or a
# non-zero exit code.
#
# The timestamps on a CACHED task come from the ORIGINAL execution, not from
# the resumed one, so the timeline below is the schedule of whichever run
# actually did the work.

tasks <- trace %>%
  transmute(
    run,
    process   = as.character(process),
    name      = as.character(name),
    status    = as.character(status),
    exit      = as.character(exit),
    cached    = status == "CACHED",
    failed    = status %in% c("FAILED", "ABORTED") |
                (!is.na(suppressWarnings(as.numeric(exit))) &
                   suppressWarnings(as.numeric(exit)) != 0),
    cpus      = num(cpus),
    realtime_s = num(realtime) / 1000,
    duration_s = num(duration) / 1000,
    peak_rss_gb = num(peak_rss) / 1024^3,
    submit_ms   = num(submit),
    start_ms    = num(start),
    complete_ms = num(complete)
  ) %>%
  filter(!is.na(realtime_s))

write_csv(tasks, file.path(out_dir, "15_task_runtimes.csv"))

# ----------------------------
# Per process per run
# ----------------------------

fmt_hms <- function(s) {
  s <- round(s)
  sprintf("%02d:%02d:%02d", s %/% 3600, (s %% 3600) %/% 60, s %% 60)
}

by_process <- tasks %>%
  group_by(run, process) %>%
  summarise(
    n_tasks        = n(),
    cached         = sum(cached),
    failed         = sum(failed),
    wall_total_s   = sum(realtime_s),
    wall_longest_s = max(realtime_s),
    cpu_hours      = sum(realtime_s * cpus, na.rm = TRUE) / 3600,
    peak_rss_gb    = max(peak_rss_gb, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    wall_total   = fmt_hms(wall_total_s),
    wall_longest = fmt_hms(wall_longest_s),
    cpu_hours    = round(cpu_hours, 2),
    peak_rss_gb  = round(peak_rss_gb, 2)
  ) %>%
  arrange(run, desc(wall_total_s))

write_csv(by_process, file.path(out_dir, "16_runtime_by_process_by_run.csv"))

# Wide view: one row per process, one column per run
wide <- by_process %>%
  dplyr::select(process, run, wall_total) %>%
  pivot_wider(names_from = run, values_from = wall_total)

write_csv(wide, file.path(out_dir, "17_runtime_by_process_wide.csv"))

# ----------------------------
# Per run totals
# ----------------------------

by_run <- tasks %>%
  group_by(run) %>%
  summarise(
    n_tasks      = n(),
    cached       = sum(cached),
    failed       = sum(failed),
    wall_total_s = sum(realtime_s),
    cpu_hours    = sum(realtime_s * cpus, na.rm = TRUE) / 3600,
    .groups = "drop"
  ) %>%
  mutate(
    wall_total_if_serial = fmt_hms(wall_total_s),
    cpu_hours = round(cpu_hours, 2)
  ) %>%
  dplyr::select(run, n_tasks, cached, failed, wall_total_if_serial, cpu_hours)

write_csv(by_run, file.path(out_dir, "18_runtime_by_run.csv"))

# ----------------------------
# Figures
# ----------------------------

proc_order <- by_process %>%
  group_by(process) %>%
  summarise(t = sum(wall_total_s), .groups = "drop") %>%
  arrange(t) %>%
  pull(process)

p_time <- by_process %>%
  mutate(process = factor(process, levels = proc_order),
         wall_total_h = wall_total_s / 3600) %>%
  ggplot(aes(process, wall_total_h, fill = run)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.75) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 11) +
  labs(
    title = "Wall-clock time per process, summed over that process's tasks",
    subtitle = "tasks within a process run in parallel, so this is total work, not elapsed time",
    x = NULL, y = "Hours", fill = NULL
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "03_runtime_by_process.pdf"), p_time, width = 10, height = 7)

p_cpu <- by_process %>%
  mutate(process = factor(process, levels = proc_order)) %>%
  ggplot(aes(process, cpu_hours, fill = run)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.75) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 11) +
  labs(
    title = "Allocated CPU hours per process",
    subtitle = "realtime x allocated cpus - what the process costs the cluster, not what it uses",
    x = NULL, y = "CPU hours", fill = NULL
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "04_cpu_hours_by_process.pdf"), p_cpu, width = 10, height = 7)

p_mem <- tasks %>%
  ggplot(aes(reorder(process, peak_rss_gb, FUN = function(x) max(x, na.rm = TRUE)),
             peak_rss_gb, colour = run)) +
  geom_point(size = 2, alpha = 0.8) +
  coord_flip() +
  scale_colour_brewer(palette = "Set2") +
  theme_bw(base_size = 11) +
  labs(title = "Peak memory per task", x = NULL, y = "Peak RSS (GB)", colour = NULL) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "05_peak_memory_by_process.pdf"), p_mem, width = 10, height = 7)

# ----------------------------
# Timeline
# ----------------------------
# Gantt of the actual schedule: one bar per task, from start to complete,
# positioned relative to the first submit of its run. Overlapping bars are what
# ran concurrently, so the width of a panel is elapsed wall clock while the sum
# of the bars is the work done.

timeline <- tasks %>%
  filter(!is.na(start_ms), !is.na(complete_ms)) %>%
  group_by(run) %>%
  mutate(
    t0      = min(submit_ms, na.rm = TRUE),
    start_h = (start_ms - t0) / 3600000,
    end_h   = (complete_ms - t0) / 3600000
  ) %>%
  ungroup()

if (nrow(timeline)) {

  proc_by_start <- timeline %>%
    group_by(process) %>%
    summarise(first_start = min(start_h), .groups = "drop") %>%
    arrange(desc(first_start)) %>%
    pull(process)

  timeline <- timeline %>% mutate(process = factor(process, levels = proc_by_start))

  write_csv(
    timeline %>%
      dplyr::select(run, process, name, cached, failed,
                    start_h, end_h, realtime_s, cpus, peak_rss_gb) %>%
      arrange(run, start_h),
    file.path(out_dir, "19_task_timeline.csv")
  )

  p_timeline <- ggplot(timeline) +
    geom_segment(
      aes(x = start_h, xend = end_h, y = process, yend = process, colour = cached),
      linewidth = 3.2, lineend = "butt"
    ) +
    facet_wrap(~ run, ncol = 1, scales = "free_x") +
    scale_colour_manual(
      values = c(`FALSE` = "#1B9E77", `TRUE` = "grey65"),
      labels = c(`FALSE` = "executed", `TRUE` = "resumed from cache"),
      name = NULL
    ) +
    theme_bw(base_size = 11) +
    labs(
      title = "Pipeline timeline: when each task ran",
      subtitle = paste("hours from the first submit of that run, one bar per task.",
                       "Cached tasks carry the timestamps of the run that did the work."),
      x = "Hours from start of run", y = NULL
    ) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())

  ggsave(file.path(out_dir, "06_pipeline_timeline.pdf"), p_timeline,
         width = 11, height = 11)

  elapsed <- timeline %>%
    group_by(run) %>%
    summarise(
      elapsed_h = max(end_h),
      work_h    = sum(realtime_s) / 3600,
      .groups = "drop"
    ) %>%
    mutate(
      concurrency = round(work_h / elapsed_h, 2),
      elapsed_h   = round(elapsed_h, 2),
      work_h      = round(work_h, 2)
    )

  write_csv(elapsed, file.path(out_dir, "20_elapsed_vs_work.csv"))

  cat("\n=== Elapsed against work done ===\n")
  print(as.data.frame(elapsed))
}

# ----------------------------
# Report
# ----------------------------

cat("\n=== Per run ===\n")
print(as.data.frame(by_run))

cat("\n=== Per process per run (wall clock, longest first) ===\n")
print(as.data.frame(
  by_process %>% dplyr::select(run, process, n_tasks, cached, failed, wall_total,
                               wall_longest, cpu_hours, peak_rss_gb)
))

cat("\nDone. Output directory:", out_dir, "\n")
