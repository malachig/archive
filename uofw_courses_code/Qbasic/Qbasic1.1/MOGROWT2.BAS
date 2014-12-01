0 CLS
1 PRINT "A MODEL OF POPULATION GROWTH"
10 REM N0 = The initial population size
20 REM Nt = The final population size
30 REM r = The population growth rate
40 REM t = The length of time in years
50 REM EXP = 2.7128 to the X
60 PRINT ""

100 INPUT "The initial population size"; N0
110 INPUT "The population growth rate"; r
120 INPUT "The length of time in years"; t
130 Nt = N0 * EXP(r * t)
140 PRINT "Population size at time t ="; Nt
150 PRINT ""

200 INPUT "The initial population size"; N0
210 INPUT "The population growth rate"; r
220 INPUT "The length of time in years"; t
230 Nt = N0 * EXP(r * t)
240 PRINT "Population size at time t ="; Nt
250 PRINT ""

300 INPUT "The initial population size"; N0
310 INPUT "The population growth rate"; r
320 INPUT "The length of time in years"; t
330 Nt = N0 * EXP(r * t)
340 PRINT "Population size at time t ="; Nt
350 PRINT ""

1000 END


