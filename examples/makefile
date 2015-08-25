EXAMPLEDIR=$(PWD)
EXAMPLES=example1 example2 example3 example4

all:
	for EXAMPLE in $(EXAMPLES); do \
	  cd $(PWD)/$$EXAMPLE ; $(MAKE); \
	done

.PHONY: clean
clean:

	for EXAMPLE in $(EXAMPLES); do \
	  cd $(PWD)/$$EXAMPLE ; $(MAKE) clean ; \
	done 
