0 CLS
1 PRINT "A MODEL OF POPULATION GROWTH"
10 REM N0 = The initial population size
20 REM Nt = The final population size
30 REM r = The population growth rate
40 REM t = The length of time in years
50 REM EXP = 2.7128 to the X

100 N0 = 100
110 r = 2
120 t = 1
200 PRINT "Initial population size ="; N0
210 PRINT "r="; r
220 PRINT "Number of years ="; t
300 Nt = N0 * EXP(r * t)
310 PRINT "Population size at time t ="; Nt
1000 END


