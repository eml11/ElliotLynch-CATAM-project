EXAMPLEDIR=$(PWD)
BIN=$(EXAMPLEDIR)/../../bin
SOLVER=$(BIN)/shooting_solution
PLOTTING=$(EXAMPLEDIR)/../../plotting
PLOTFILE=$(PLOTTING)/stellarprofile_savefig.py

all: example2_figure.png

example2_figure.png: stellarstructure_inner.txt	stellarstructure_outer.txt
	python $(PLOTFILE) in example2_figure.png

stellarstructure_inner.txt stellarstructure_outer.txt :star.inp $(SOLVER)
	$(SOLVER)

clean:
	rm -f stellarstructure_outer.txt
	rm -f example2_figure.png
