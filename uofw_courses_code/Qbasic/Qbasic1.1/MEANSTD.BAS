CLS
10 REM find the mean and standard deviation
20 REM define the variable names
30 REM MEAN= mean of variable x
40 REM STD = standard deviation
50 REM SUMX = the sum of x
60 REM SUMSQ = the sum of X squared
100 DIM x(1000): REM create an array with slots for 100 numbers to be stored.
110 PRINT "Calculation of the mean and standard deviation"
120 GOSUB 1000: REM input the data
130 GOSUB 5000: REM find the mean and standard deviation
140 PRINT
150 PRINT "Arithmetic mean ="; mean
160 PRINT "Standard deviation ="; STD
170 PRINT "Sample size, n ="; n
180 END

1000 REM input the data
1010 INPUT "How many data points"; n
1020 IF n > 1000 THEN PRINT "sample size too large try a smaller number"
1030 IF n > 1000 THEN GOTO 1010
1100 FOR i = 1 TO n
1110 PRINT i
INPUT x(i)
1130 NEXT i
1140 RETURN

5000 REM find the mean and standard deviation
5010 sumx = 0: sumsq = 0
5030 FOR i = 1 TO n
5040 sumx = sumx + x(i)
5050 sumsq = sumsq + (x(i) ^ 2)
5060 NEXT i
5100 mean = sumx / n
5110 STD = SQR((sumsq - (sumx ^ 2) / n) / (n - 1))
5200 RETURN



