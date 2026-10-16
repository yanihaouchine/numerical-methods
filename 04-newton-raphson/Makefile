LATEX = pdflatex

all: rapport.pdf

rapport: rapport.pdf
	evince rapport.pdf &

rapport.pdf: rapport.tex img/*
	$(LATEX) rapport.tex
	$(LATEX) rapport.tex

clean:
	rm -f *.aux *.log *.out *.pdf *.toc

.PHONY: clean rapport