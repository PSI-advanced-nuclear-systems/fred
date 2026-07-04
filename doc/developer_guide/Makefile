# developer_guide/Makefile
#
# Usage:
#   make          — scan sources, then compile LaTeX (two passes)
#   make scan     — re-run coderef scanner only
#   make pdf      — compile LaTeX only (assumes coderef_defs.tex exists)
#   make clean    — remove generated and build artefacts

PYTHON   ?= python3
SCANNER  := latex-code-ref/scan_coderefs.py
SRC_DIR  := ../../src
EXTS     := "*.hpp" "*.cpp"

TEX_FILE  := main.tex
PDF_FILE  := main.pdf
DB_JSON   := coderef_db.json
DB_TEX    := coderef_defs.tex

# Let pdflatex find coderef.sty in the latex-code-ref/ subdirectory.
TEXINPUTS := .:latex-code-ref:

.PHONY: all scan pdf clean

all: scan pdf

scan: $(DB_TEX)

$(DB_TEX):
	$(PYTHON) $(SCANNER) $(SRC_DIR) --extensions $(EXTS) \
	    --json-output $(DB_JSON) --latex-output $(DB_TEX)

pdf: $(DB_TEX)
	TEXINPUTS=$(TEXINPUTS) pdflatex -halt-on-error $(TEX_FILE)
	TEXINPUTS=$(TEXINPUTS) pdflatex -halt-on-error $(TEX_FILE)

clean:
	rm -f $(DB_JSON) $(DB_TEX) \
	      *.aux *.log *.out *.toc *.nav *.snm $(PDF_FILE)
