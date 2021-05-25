IPYNB :=quantum_dots_NV-_center_FRET_estimation.ipynb
HTML  :=quantum_dots_NV-_center_FRET_estimation.html
PY    :=$(patsubst %.ipynb,%.html, $(IPYNB))
PY    :=$(patsubst %.ipynb,%.py, $(IPYNB))
PDF   :=$(patsubst %.ipynb,%.pdf, $(IPYNB))

.PHONY: all
all : $(HTML) $(PY) $(PDF)

.PHONY: ipython-notebook
ipython-notebook:
	jupyter notebook

.PHONY: jupyter-notebook
jupyter-notebook:
	jupyter notebook

$(HTML) : $(IPYNB)
	jupyter nbconvert --to html $< --stdout > $@

$(PDF) : $(IPYNB)
	jupyter nbconvert --to pdf $<

$(PY) : $(IPYNB)
	jupyter nbconvert --to script $<


.PHONY: clean
clean :
	rm -f -- $(HTML) $(PY) $(PDF)

# man /usr/share/man/man1/jupyter-nbconvert.1.gz
# https://ipython.org/ipython-doc/dev/notebook/nbconvert.html
