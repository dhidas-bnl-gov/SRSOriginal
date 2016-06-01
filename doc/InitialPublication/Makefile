CC = latex


all:
	echo 'What do you want to make?'



%.dvi : %.tex
	$(CC) $< && $(CC) $< 

%.ps : %.dvi
	dvips $< -o $@
#	dvips $< -o $@ && gv $@

#%.pdf : %.ps
#	ps2pdf $< && open $@

%.pdf : %.tex
	pdflatex $< && open $@




clean:
	rm -r -f *.aux *.dvi *.log *.out *.pdf *.ps *.blg *.bbl *.cls
