REM THIS PROGRAM IS A MUCH FASTER WAY TO SOLVE THE PROBLEM
REM IT ALSO CALCULATES THE # OF ATTEMPTS REQUIRED TO MAKE THE SENTENCE
REM "To boldly go where no man has gone before"
REM There are 41 characters in this statement.
REM 26 letter + space are the possibilities.
REM Therefore there is a 1/27 probability for each character

REM generate random numbers between 1 and 27
REM Assign each a letter or space
REM Use a selection process to generate the sentence from random numbers.

CLS
RANDOMIZE TIMER
OPEN "C:\Monkey2.txt" FOR OUTPUT AS #1

DIM array1(41)
DATA 20,15,27,2,15,12,4,12,25,27,7,15,27,23,8,5,18,5,27,14,15,27,15,14,5
DATA 27,8,1,19,27,7,15,14,5,27,2,5,6,15,18,5
FOR i = 1 TO 41
READ array1(i)
NEXT i

FOR i = 1 TO 41
ATTEMPTS = 1
GOSUB 1000
WHILE RL <> array1(i)
ATTEMPTS = ATTEMPTS + 1
GOSUB 1000
WEND
PRINT #1, RL, RL$, "ATTEMPTS =", ATTEMPTS
TOTAL = TOTAL + ATTEMPTS
NEXT i
PRINT #1, "TOTAL =", TOTAL
CLOSE #1
END

1000
U = 27
L = 1
RL = INT((U - L + 1) * RND + L)
IF RL = 1 THEN RL$ = "A"
IF RL = 2 THEN RL$ = "B"
IF RL = 3 THEN RL$ = "C"
IF RL = 4 THEN RL$ = "D"
IF RL = 5 THEN RL$ = "E"
IF RL = 6 THEN RL$ = "F"
IF RL = 7 THEN RL$ = "G"
IF RL = 8 THEN RL$ = "H"
IF RL = 9 THEN RL$ = "I"
IF RL = 10 THEN RL$ = "J"
IF RL = 11 THEN RL$ = "K"
IF RL = 12 THEN RL$ = "L"
IF RL = 13 THEN RL$ = "M"
IF RL = 14 THEN RL$ = "N"
IF RL = 15 THEN RL$ = "O"
IF RL = 16 THEN RL$ = "P"
IF RL = 17 THEN RL$ = "Q"
IF RL = 18 THEN RL$ = "R"
IF RL = 19 THEN RL$ = "S"
IF RL = 20 THEN RL$ = "T"
IF RL = 21 THEN RL$ = "U"
IF RL = 22 THEN RL$ = "V"
IF RL = 23 THEN RL$ = "W"
IF RL = 24 THEN RL$ = "X"
IF RL = 25 THEN RL$ = "Y"
IF RL = 26 THEN RL$ = "Z"
IF RL = 27 THEN RL$ = " "
RETURN

