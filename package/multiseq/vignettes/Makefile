# Makefile to use knitr for package vignettes
# put all PDF targets here, separated by spaces
#PDFS= multiseq-intro.pdf
HTMLS=multiseq.html

#all: $(PDFS)

all: $(HTMLS)

clean:	
	rm -rf *.R *.tex *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ cache/

%.html:	%.Rmd
	$(R_HOME)/bin/Rscript -e "if (getRversion() < '3.0.0') knitr::knit2html('$*.Rmd')"

#%.pdf:  %.Rmd
#	$(R_HOME)/bin/Rscript -e "rmarkdown::render('$*.Rmd', '$*')"

