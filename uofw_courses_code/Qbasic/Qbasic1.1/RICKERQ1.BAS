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
50 FOR a = 1 TO 4
60 FOR S = 100000 TO 2000000 STEP 100000
70 R = S * EXP(a * (1 - (S / B)))
80 PRINT "Rate of growth ="; a; "Spawners ="; S; "Recruits ="; R
90 NEXT S
100 INPUT "Begin next value of a"; Answer$
110 IF Answer$ = "yes" GOTO 200 ELSE GOTO 7000
200 GOSUB 1000
210 NEXT a
7000 END

1000 PRINT "The Ricker Curve for Predicting Number of Recruits"
1010 PRINT "The Carrying Capacity ="; B
1020 RETURN

