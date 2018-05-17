tagm.tex: tagm.Rnw
	Rscript -e "knitr::knit('tagm.Rnw')"
tagm.pdf: tagm.tex
	pdflatex tagm.tex
