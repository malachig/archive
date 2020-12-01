CLS
REM When mother cease menstruation, this is called amenorrhea.
REM When breast feeding, mothers experience periods of amenorrhea.
REM The following statistical relationship was developed by work
REM done in Taiwan, Senegal, the Philippines and the United States.

REM A = the mean duration of post-partum amenorrhea in months.
REM B = the mean duration of breast feeding in months.
REM e = the base of the natural logarithm, 2.71...
REM Calculate A for B = 0 to 60

PRINT "This program calculates the mean duration of amenorrhea for a given"
PRINT "duration of breast feeding"
OPEN "c:\qbasic\assign1.txt" FOR OUTPUT AS #1
FOR B = 0 TO 60 STEP 1
A = 1.73 * EXP((.14 * B) - (.0019 * (B ^ 2)))
PRINT A
PRINT #1, A
NEXT B
END

