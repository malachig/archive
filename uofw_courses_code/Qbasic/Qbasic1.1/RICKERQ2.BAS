REM THE MODEL VARIABLE ARE:
REM R = recruits (or the population size in generation t+1)
REM S = the number of parent spawners in the current population

REM THE MODEL PARAMETERS ARE: (determine the shape of the curve)
REM a = the rate of growth of a very small population
REM B = the carrying capacity of spwners (also a scale parameter)

REM e = the base of the natural logarithm 2.7128

10 CLS
20 B = 1000000
30 PRINT "The Ricker Curve for Predicting Number of Recruits"
40 PRINT "The Carrying Capacity ="; B
50 S = 10
60 PRINT "Initial Spawners ="; S
70 INPUT "Rate of Growth, a ="; a
80 INPUT "The number of Generation, t ="; t
90 FOR t = 1 TO t
100 R = S * EXP(a * (1 - (S / B)))
110 PRINT "The Number of Recruits (Population at t + 1) ="; R
120 S = R
130 NEXT t
140 PRINT "Rate of growth ="; a; "Spawners ="; S; "Recruits ="; R
150 END




