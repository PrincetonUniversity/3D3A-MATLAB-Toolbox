TEX_FILES := $(wildcard *.tex)
PDF_FILES := $(patsubst %.tex, %.pdf, $(TEX_FILES))
INTERMEDIATES := *.out *.aux *.log *.synctex.gz *.bbl *.blg *.nlo
FIG_INTERMEDIATES := *-eps-converted-to.pdf **/*-eps-converted-to.pdf

.PHONY: mostlyclean

# By default, make all files and then clean up
default: all mostlyclean

# Generate PDFs, leaving any intermediate files
all: $(PDF_FILES)

# Generate PDF for a given TeX file
%.pdf: %.tex
	if ls *.bib 1> /dev/null 2>&1; then pdflatex $<; fi;
	if ls *.bib 1> /dev/null 2>&1; then bibtex $*; fi;
	pdflatex $<
	pdflatex $<
	
# Remove all TeX intermediates
mostlyclean:
	$(RM) $(INTERMEDIATES)

# Remove TeX intermediates and figure intermediates
clean: mostlyclean
	$(RM) $(FIG_INTERMEDIATES)

# Remove all intermediates and final PDFs
distclean: clean
	$(RM) $(PDF_FILES)