# BAP_Int

This project is the source code for the article <Buffer Allocation Problem in production flow lines: a new Benders-decomposition-based exact solution approach>, accepted in IISE Transactions.

A simple example of solving the buffer allocation problem (BAP) of production flow line can be found in src/main/java/test_mainfunc/TestAlter6.java. Specifically, 
  1. an object of serial line is constructed by reading the input file /INPUT/SerialLine_test_6stage.txt
  2. parameters of BAP, such as buffer lowerbound, buffer upperbound, target throughput and simulation length, are defined
  3. the samples of processing time tij are generated. 
  4. an object myReversedAlter6 of class BendersIntModelAlter6ReversedCut is constructed by providing an object of SerialLine, target throughput, buffer lowerbound, buffer upperbound and simulation length.
  5. the method solveBAPWithIntModel(tij,false) is called to solve the BAP.
  6. the solution is then printed in the output file /OUTPUT/Alter6_6stage.txt.
