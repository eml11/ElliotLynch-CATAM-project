EXAMPLEDIR=$(PWD)
BIN=$(EXAMPLEDIR)/../../bin
SOLVER=$(BIN)/stellar_adiabat
PLOTTING=$(EXAMPLEDIR)/../../plotting
PLOTFILE=$(PLOTTING)/adiabatprofile_savefig.py

all: example1_figure.png

example1_figure.png: stellarstructure_inner.txt	stellarstructure_outer.txt
	python $(PLOTFILE) example1_figure.png

stellarstructure_inner.txt stellarstructure_outer.txt :star.inp $(SOLVER)
	$(SOLVER)

clean:
	rm -f stellarstructure_inner.txt stellarstructure_outer.txt
	rm -f example1_figure.png
