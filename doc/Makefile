CC = pdflatex


all: SRS_CodingConventions.pdf SRS_PythonAPI.pdf
	echo 'Party!'



%.pdf : %.tex
	$(CC) $< 





clean:
	rm -r -f *.aux *.dvi *.log *.out *.pdf *.ps *.blg *.bbl
