# Least-Squares-Regression
Least-squares regression polynomial fit calculator. Done in 2017.

### Purpose
This program allows you to input a data set of (x, y) pairs and returns a set of coefficients for the polynomial of best fit with order 
specified by the user.

### Method
The program accepts x and y values and turns them into column matrices **x** and **y** and calculates the corresponding Vandermonde matrix
**V**. Then we calculate the coefficient column matrix **c** using the formula: 

**c** = (**V**<sup>T</sup>**V**)<sup>-1</sup>**V**<sup>T</sup>**y**.

The first entry of **c** will be the coefficient of the highest order term; coefficients of terms are listed by order in descending 
sequence.
