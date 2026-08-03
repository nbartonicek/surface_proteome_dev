#!/bin/bash

# ------------------------------------------------------------------
# LK1 pilot run - render helper (not a numbered analysis step)
#
# Renders the report Rmd to HTML, then prints that HTML to PDF.
#
# The PDF exists because LabArchives previews PDFs inline in the entry, so the
# report is readable in the notebook itself, whereas an attached .html is only
# ever a download. Keep both: the HTML is the reproducible artifact, the PDF is
# what people actually read.
#
# PDF conversion goes through headless Chrome rather than pdf_document, because
# that renders the real page - same flatly theme, same base64-embedded figures -
# and needs no LaTeX. tinytex is installed on this laptop but has no LaTeX
# distribution, so pdf_document would fail.
#
# pandoc is not on PATH here; it lives inside the RStudio bundle.
#
# Usage:
#   ./render_report.sh                      # render HTML then PDF
#   ./render_report.sh --pdf-only           # skip the render, just redo the PDF
# ------------------------------------------------------------------

set -uo pipefail

cd "$(dirname "$0")" || exit 1

RMD="19.LK1_pilot_run_report.Rmd"
OUT="LK1_pilot_run_report"

export RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64"
CHROME="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"

# ----------------------------
# 1. Rmd -> HTML
# ----------------------------

if [ "${1:-}" != "--pdf-only" ]; then
  echo "rendering $RMD -> $OUT.html"
  Rscript -e "rmarkdown::render('$RMD', output_file='$OUT.html', quiet=TRUE)" || {
    echo "render failed"; exit 1; }
fi

[ -f "$OUT.html" ] || { echo "no $OUT.html to convert"; exit 1; }

# ----------------------------
# 2. HTML -> PDF
# ----------------------------
# --virtual-time-budget lets the floating TOC's javascript settle before the
# page is printed; without it the TOC can come out empty.

echo "printing $OUT.html -> $OUT.pdf"
"$CHROME" --headless=new --disable-gpu --no-sandbox --no-pdf-header-footer \
  --virtual-time-budget=30000 --run-all-compositor-stages-before-draw \
  --print-to-pdf="$PWD/$OUT.pdf" "file://$PWD/$OUT.html" 2>/dev/null

# ----------------------------
# 3. Report what came out
# ----------------------------

if [ -f "$OUT.pdf" ]; then
  Rscript -e "
    i <- pdftools::pdf_info('$OUT.pdf')
    t <- pdftools::pdf_text('$OUT.pdf')
    cat(sprintf('  %s  %.1f MB, %d pages, %d characters of text\n',
                '$OUT.pdf', file.size('$OUT.pdf')/1048576, i\$pages, sum(nchar(t))))
  " 2>/dev/null
else
  echo "PDF conversion produced nothing"
  exit 1
fi

echo "done."
