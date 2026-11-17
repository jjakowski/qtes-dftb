DOC = qtes-dftb-manual

$(DOC).pdf: $(DOC).tex
	pdflatex -interaction=nonstopmode -halt-on-error $(DOC).tex
	pdflatex -interaction=nonstopmode -halt-on-error $(DOC).tex
	pdflatex -interaction=nonstopmode -halt-on-error $(DOC).tex

clean:
	rm -f $(DOC).aux $(DOC).log $(DOC).out $(DOC).toc

.PHONY: clean
