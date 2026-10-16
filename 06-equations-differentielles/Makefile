LATEX = pdflatex

all: rapport

rapport: rapport.pdf
	evince rapport.pdf &

rapport.pdf: rapport/rapport.tex rapport/img/*.png
	$(LATEX) rapport/rapport.tex
	$(LATEX) rapport/rapport.tex

clean:
	rm -f *.aux *.log *.out *.pdf *.toc

.PHONY: clean rapport
