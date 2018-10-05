all:
	make tagm.pdf
	make si.pdf

tagm.tex: tagm.Rnw
	Rscript -e "knitr::knit('tagm.Rnw')"
tagm.pdf: tagm.tex
	pdflatex tagm.tex
	pdflatex tagm.tex


si.tex: si.Rnw
	Rscript -e "knitr::knit('si.Rnw')"
si.pdf: si.tex
	pdflatex si.tex
	pdflatex si.tex
