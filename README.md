# ElliotLynch-CATAM-project
CATAM Stellar Structure project

To compile type:
make

To run examples type:
make examples

binaries will be located in : bin/

binaries:

  stellar\_adiabat : solves for adiabatic equations of
                     stellar structure using shooting method

  shooting\_solution : solves full equations of stellar stucture
                       using shooting method, can preform inward,
                       and outward integration or to the midpoint
 
  stellarstructure : refines parameters using jacobian iteration
                     before solving full equations using shooting
                     method

sample input files are found in: examples/

examples:
  example1: sample input for question4 with adiabatic eos

  example2: sample input for question5 full equations
            inward integration

  example3 sample input for question5 full equations
           shooting method, inward and outward integration

  example4 sample input for quation6 full equations
           using jacobian iteration method

running make in an example directory will run the appropriate
binary and plot the results:
  stellarstructure\_inner.txt : results of inner integral (if
  present)
  stellarstructure\_outer.txt : results of outer integral

  example\*\_figure.png : plot of stellar properties 
