# Fast-Fourier-Transform-for-Polynomial-Multiplication
Fast Fourier Transform for Polynomial Multiplication

Written in Java.
Uncomment the code between the horizontal dotted lines for corresponding solutions in Java file.
 
Change the length variable to change the number of input objects 
All solutions are written in seperate functions. To execute each uncomment the correponding section.

                 ------------------------------------------------------------------------

Folders -

Raw Data - Inputs used for testing all the algorithms

Java Code - Code used for the algos

Graphs - Graphs with different parameters to measure the efficiency of the algorithm

Report - Detailed information about graphs and observations.

                 ------------------------------------------------------------------------

Implemented recursive and dynamic programming version of the FFT and apply it to perform polynomial multiplication.
Compared the performance of the FFT method with the previous algebraic techniques both w.r.t speed and accuracy.

1)  Implemented the recursive FFT with the following optimizations: 
  
  (a) pre-compute the omega table (do not count this as part of the run time) and use the indexing technique to compute the x2 values as needed.

2) Applied your FFT method to polynomial multiplication of two polynomials P Q size n by 
      
      (a) padding the polynomials with high-order 0’s to make them size 2n, 
      
      (b) evaluate P and Q at 2n values (powers of omega), 
      
      (c) multiple P(xi)*Q(xi) to obtain samples of PQ(xi), 
      
      (d) use the inverse FFT to interpolate the coefficients of PQ.

3)  Checked that the values computed by the FFT sampling technique are the same as those obtained by the previous approaches. There will be some slight rounding errors due limited resolution of the doubles (or longs).

4)  Implemented FFT using the optimized DP method: 
      
      (a) No dynamic allocation, 
      
      (b) pre-computed bit shuffle for coefficient lookup (don’t count in the timing), 
      
      (c) 2 row space optimization.

5)  Executed timing studies comparing the performance of recursive with DP FFT for the largest problems that can be feasibly solved.
      Plot the results on a log-log graph, measure the vertical offset between the two methods and compute the corresponding factor.

6)  Compared the run-time performance of FFT with the previous methods. Ploted the time performance for the four algorithms on one graph. The maximum problem size solvable will increase with algorithm performance. Given a table of crossing points.

7)  Compared the accuracy of FFT with the algebraic methods. For larger problem sizes computed the mean absolute error between the coefficient values obtained by FFT and the n1.59 method as a function of problem size. The x axis us problem size (increasing geometrically) vs. the mean absolute error (on a linear scale). 

For the three empirical studies above (5, 6, 7) did randomized repeats and ploted the average value on the graphs.

8)  Implemented two overloads of the complex number class. One where the complex multiply uses 3 real multiplies and one where it uses 4. Compared the performance using these two implementations running the DP FFT for the largest problem size.

