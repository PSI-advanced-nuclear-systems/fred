# Makefile for the FRED 2.0 Theory Manual
#
# Requires: pdflatex, biber (TeXLive 2020+), Python 3 with matplotlib
#
# Usage:
#   make          → compile PDF (default target)
#   make figures  → regenerate figure PDFs from Python scripts
#   make clean    → remove auxiliary files
#   make distclean → remove auxiliary files and output PDF

MAIN     = main
TEXSRC   = $(MAIN).tex introduction.tex thermal_conductivity_stress_strain.tex \
		   fred_ox.tex fred_m_na.tex references.bib

FIGURES  = images/rz_discretization.pdf images/gap_schematic.pdf images/mech_ring_compatibility.pdf \
	       images/waltar_reynolds_zones.pdf

PDFLATEX = pdflatex -interaction=nonstopmode -halt-on-error
BIBER    = biber

.PHONY: all figures clean distclean

all: $(MAIN).pdf

# Full LaTeX + Biber compile cycle
$(MAIN).pdf: $(TEXSRC) $(FIGURES)
	$(PDFLATEX) $(MAIN)
	$(BIBER)    $(MAIN)
	$(PDFLATEX) $(MAIN)
	$(PDFLATEX) $(MAIN)

# Generate figures from Python scripts
figures: $(FIGURES)

# Remove auxiliary files produced by pdflatex / biber
clean:
	rm -f $(MAIN).aux $(MAIN).bbl $(MAIN).bcf $(MAIN).blg \
	      $(MAIN).lof $(MAIN).lot $(MAIN).out $(MAIN).run.xml \
	      $(MAIN).toc $(MAIN).log

distclean: clean
	rm -f $(MAIN).pdf
