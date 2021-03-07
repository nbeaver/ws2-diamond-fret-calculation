IPYNB :=$(wildcard *.ipynb)
HTML  :=$(patsubst %.ipynb,%.html, $(IPYNB))
PY    :=$(patsubst %.ipynb,%.py, $(IPYNB))
PDF   :=$(patsubst %.ipynb,%.pdf, $(IPYNB))

ZIP := fret_calculation.zip

.PHONY: all
all : $(HTML) $(PY) $(PDF)

%.html : %.ipynb
	jupyter nbconvert --to html $<

%.pdf : %.ipynb
	jupyter nbconvert --to pdf $<

%.py : %.ipynb
	jupyter nbconvert --to script $<


.PHONY: clean
clean :
	rm -f -- $(HTML) $(PY) $(PDF)

.PHONY: zip
zip : $(ZIP)

.PHONY: $(ZIP)
$(ZIP): Makefile .gitignore WS2_NV-_center_FRET_estimation.ipynb Makefile NV_center_absorption.txt WS2_emission.txt
	zip --filesync --quiet $@ $^

# man /usr/share/man/man1/jupyter-nbconvert.1.gz
# https://ipython.org/ipython-doc/dev/notebook/nbconvert.html
