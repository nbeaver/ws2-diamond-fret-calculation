IPYNB :=$(wildcard *.ipynb)
HTML  :=$(patsubst %.ipynb,%.html, $(IPYNB))
PY    :=$(patsubst %.ipynb,%.py, $(IPYNB))
PDF   :=$(patsubst %.ipynb,%.pdf, $(IPYNB))

ZIP :=known_dyes_fret_calculation.zip

.PHONY: all
all : $(HTML) $(PY) $(PDF)

.PHONY: ipython-notebook
ipython-notebook:
	jupyter notebook

.PHONY: jupyter-notebook
jupyter-notebook:
	jupyter notebook

%.html : %.ipynb
	jupyter nbconvert --to html $<

%.pdf : %.ipynb
	jupyter nbconvert --to pdf $<

%.py : %.ipynb
	jupyter nbconvert --to script $<

.PHONY: clean
clean :
	rm -f -- $(HTML) $(PY) $(PDF) $(ZIP)

.PHONY: zip
zip : $(ZIP)

.PHONY: $(ZIP)
$(ZIP): Makefile .gitignore $(IPYNB) ATTO_390_emission.csv ATTO_425_absorption.csv fpbase/FPbase_FRET_Calculator.png ae/ae_UV-Vis_Spectroscopy_Application.png
	zip --filesync --quiet $@ $^

# man /usr/share/man/man1/jupyter-nbconvert.1.gz
# https://ipython.org/ipython-doc/dev/notebook/nbconvert.html
