
SPHINXOPTS    ?=
SPHINXBUILD   ?= jupyter-book build
BUILDDIR      = _build

.PHONY: html clean

html:
	$(SPHINXBUILD) . $(SPHINXOPTS)

clean:
	-rm -rf $(BUILDDIR)
