#!/bin/bash

# ------------------------------------------------------------------
# Nextflow workflow - render helper (not a numbered analysis step)
#
# Renders the report Rmd to HTML.
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
# There is also a .docx / .rtf route, for pasting the report straight into the
# LabArchives rich text editor rather than attaching a file. Both carry the
# figures (all 37 land in the file), which the HTML route does not survive -
# rich text editors generally strip the base64 data: URIs the self-contained
# HTML uses. Of the two, .docx is the one to try first: copying from Word puts
# both RTF and HTML on the clipboard, which is what web editors handle best.
#
# These are converted from the rendered HTML, not from a second knit, so the
# content is identical to the HTML and the Rmd does not need a word_document
# entry. Word styling is plainer than the flatly HTML - pandoc maps to Word's
# own styles - but the text, tables and figures all come through.
#
# Usage:
#   ./render_report.sh                      # render the Rmd to HTML
#   ./render_report.sh --docx               # ... and write .docx for pasting
#   ./render_report.sh --rtf                # ... and write .rtf for pasting
#   ./render_report.sh --pdf-only --docx    # no re-render, refresh PDF + docx
# ------------------------------------------------------------------

set -uo pipefail

cd "$(dirname "$0")" || exit 1

PDF_ONLY=0; WANT_DOCX=0; WANT_RTF=0
for a in "$@"; do
  case "$a" in
    --pdf-only) PDF_ONLY=1 ;;
    --docx)     WANT_DOCX=1 ;;
    --rtf)      WANT_RTF=1 ;;
    *) echo "unknown option: $a"; exit 1 ;;
  esac
done

RMD="03.nextflow_workflow_report.Rmd"
OUT="nextflow_workflow_report"

export RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64"
CHROME="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"

# ----------------------------
# 1. Rmd -> HTML
# ----------------------------

if [ "$PDF_ONLY" -eq 0 ]; then
  echo "rendering $RMD -> $OUT.html"
  Rscript -e "rmarkdown::render('$RMD', output_file='$OUT.html', quiet=TRUE)" || {
    echo "render failed"; exit 1; }
fi

[ -f "$OUT.html" ] || { echo "no $OUT.html to convert"; exit 1; }

# ----------------------------
# 2. Report what came out
# ----------------------------
# HTML only. No PDF step - the report is read in a browser, and the headless
# Chrome print stage was both slow and prone to hanging.

if [ -f "$OUT.html" ]; then
  echo "  $OUT.html  $(du -h "$OUT.html" | cut -f1)"
else
  echo "no $OUT.html produced"
  exit 1
fi

echo "done."
