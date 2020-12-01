CLS
INPUT "Desired mean"; MEAN
INPUT "Desired standard deviation"; STDEV
SUM = 0
REM we are getting the sum of 12 random numbers
FOR i = 1 TO 12
SUM = SUM + RND
NEXT i
REM now we will generate a Gaussian random number
GAUSSRANDOM = (SUM - 6) * STDEV + MEAN
PRINT GAUSSRANDOM
END

