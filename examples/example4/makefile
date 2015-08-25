EXAMPLEDIR=$(PWD)
BIN=$(EXAMPLEDIR)/../../bin
SOLVER=$(BIN)/stellarstructure
PLOTTING=$(EXAMPLEDIR)/../../plotting
PLOTFILE=$(PLOTTING)/stellarprofile_savefig.py

all: example4_figure.png

example4_figure.png: stellarstructure_inner.txt	stellarstructure_outer.txt
	python $(PLOTFILE) inout example4_figure.png

stellarstructure_inner.txt stellarstructure_outer.txt :star.inp $(SOLVER)
	$(SOLVER)

clean:
	rm -f stellarstructure_inner.txt stellarstructure_outer.txt
	rm -f example4_figure.png
