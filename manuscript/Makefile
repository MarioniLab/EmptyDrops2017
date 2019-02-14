FIGURES=Figure1.pdf Figure2.pdf Figure3.pdf
TEXT=LunEtAl_Manuscript.pdf LunEtAl_Supplementary.pdf

all: $(TEXT) $(FIGURES)

Figure1.pdf: pics/figure1.tex
Figure2.pdf: pics/figure2.tex
Figure3.pdf: pics/figure3.tex

$(FIGURES):
	cd pics; pdflatex `basename $<`
	pdfcrop `echo $< | sed "s/tex/pdf/"` tmp.pdf
	mv tmp.pdf $@

LunEtAl_Manuscript.pdf: debrief.tex
LunEtAl_Supplementary.pdf: additional.tex

$(TEXT):
	pdflatex $<
	bibtex `echo $< | sed "s/\\.tex//"`
	pdflatex $<
	pdflatex $<
	mv `echo $< | sed "s/tex/pdf/"` $@
