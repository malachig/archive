REM THIS PROGRAM IS THE MOST INELEGANT WAY TO SOLVE THE PROBLEM
REM "To boldly go where no man has gone before"
REM There are 41 characters in this statement.
REM 26 letter + space are the possibilities.
REM Therefore there is a 1/27 probability for each character

REM generate random numbers between 1 and 27
REM Assign each a letter or space
REM Use a selection process to generate the sentance from random numbers.

CLS
RANDOMIZE TIMER
OPEN "C:\Monkey1.txt" FOR OUTPUT AS #1
GOSUB 1000

WHILE RL$ <> "T"
GOSUB 1000
WEND
RL1$ = RL$

WHILE RL$ <> "O"
GOSUB 1000
WEND
RL2$ = RL$

WHILE RL$ <> " "
GOSUB 1000
WEND
RL3$ = RL$
PRINT #1, RL1$ + RL2$ + RL3$

WHILE RL$ <> "B"
GOSUB 1000
WEND
RL4$ = RL$
WHILE RL$ <> "O"
GOSUB 1000
WEND
RL5$ = RL$
WHILE RL$ <> "L"
GOSUB 1000
WEND
RL6$ = RL$
WHILE RL$ <> "D"
GOSUB 1000
WEND
RL7$ = RL$
WHILE RL$ <> "L"
GOSUB 1000
WEND
RL8$ = RL$
WHILE RL$ <> "Y"
GOSUB 1000
WEND
RL9$ = RL$
WHILE RL$ <> " "
GOSUB 1000
WEND
RL10$ = RL$
PRINT #1, RL1$ + RL2$ + RL3$ + RL4$ + RL5$ + RL6$ + RL7$ + RL8$ + RL9$ + RL10$

WHILE RL$ <> "G"
GOSUB 1000
WEND
RL11$ = RL$
WHILE RL$ <> "O"
GOSUB 1000
WEND
RL12$ = RL$
WHILE RL$ <> " "
GOSUB 1000
WEND
RL13$ = RL$
PRINT #1, RL1$ + RL2$ + RL3$ + RL4$ + RL5$ + RL6$ + RL7$ + RL8$ + RL9$ + RL10$ + RL11$ + RL12$ + RL13$

WHILE RL$ <> "W"
GOSUB 1000
WEND
RL14$ = RL$
WHILE RL$ <> "H"
GOSUB 1000
WEND
RL15$ = RL$
WHILE RL$ <> "E"
GOSUB 1000
WEND
RL16$ = RL$
WHILE RL$ <> "R"
GOSUB 1000
WEND
RL17$ = RL$
WHILE RL$ <> "E"
GOSUB 1000
WEND
RL18$ = RL$
WHILE RL$ <> " "
GOSUB 1000
WEND
RL19$ = RL$
PRINT #1, RL1$ + RL2$ + RL3$ + RL4$ + RL5$ + RL6$ + RL7$ + RL8$ + RL9$ + RL10$ + RL11$ + RL12$ + RL13$ + RL14$ + RL15$ + RL16$ + RL17$ + RL18$ + RL19$

WHILE RL$ <> "N"
GOSUB 1000
WEND
RL20$ = RL$
WHILE RL$ <> "O"
GOSUB 1000
WEND
RL21$ = RL$
WHILE RL$ <> " "
GOSUB 1000
WEND
RL22$ = RL$
PRINT #1, RL1$ + RL2$ + RL3$ + RL4$ + RL5$ + RL6$ + RL7$ + RL8$ + RL9$ + RL10$ + RL11$ + RL12$ + RL13$ + RL14$ + RL15$ + RL16$ + RL17$ + RL18$ + RL19$ + RL20$ + RL21$ + RL22$

WHILE RL$ <> "O"
GOSUB 1000
WEND
RL23$ = RL$
WHILE RL$ <> "N"
GOSUB 1000
WEND
RL24$ = RL$
WHILE RL$ <> "E"
GOSUB 1000
WEND
RL25$ = RL$
WHILE RL$ <> " "
GOSUB 1000
WEND
RL26$ = RL$
PRINT #1, RL1$ + RL2$ + RL3$ + RL4$ + RL5$ + RL6$ + RL7$ + RL8$ + RL9$ + RL10$ + RL11$ + RL12$ + RL13$ + RL14$ + RL15$ + RL16$ + RL17$ + RL18$ + RL19$ + RL20$ + RL21$ + RL22$ + RL23$ + RL24$ + RL25$ + RL26$

WHILE RL$ <> "H"
GOSUB 1000
WEND
RL27$ = RL$
WHILE RL$ <> "A"
GOSUB 1000
WEND
RL28$ = RL$
WHILE RL$ <> "S"
GOSUB 1000
WEND
RL29$ = RL$
WHILE RL$ <> " "
GOSUB 1000
WEND
RL30$ = RL$
PRINT #1, RL1$ + RL2$ + RL3$ + RL4$ + RL5$ + RL6$ + RL7$ + RL8$ + RL9$ + RL10$ + RL11$ + RL12$ + RL13$ + RL14$ + RL15$ + RL16$ + RL17$ + RL18$ + RL19$ + RL20$ + RL21$ + RL22$ + RL23$ + RL24$ + RL25$ + RL26$ + RL27$ + RL28$ + RL29$ + RL30$

WHILE RL$ <> "G"
GOSUB 1000
WEND
RL31$ = RL$
WHILE RL$ <> "O"
GOSUB 1000
WEND
RL32$ = RL$
WHILE RL$ <> "N"
GOSUB 1000
WEND
RL33$ = RL$
WHILE RL$ <> "E"
GOSUB 1000
WEND
RL34$ = RL$
WHILE RL$ <> " "
GOSUB 1000
WEND
RL35$ = RL$
PRINT #1, RL1$ + RL2$ + RL3$ + RL4$ + RL5$ + RL6$ + RL7$ + RL8$ + RL9$ + RL10$ + RL11$ + RL12$ + RL13$ + RL14$ + RL15$ + RL16$ + RL17$ + RL18$ + RL19$ + RL20$ + RL21$ + RL22$ + RL23$ + RL24$ + RL25$ + RL26$ + RL27$ + RL28$ + RL29$ + RL30$ + RL31$ + RL32$ + RL33$ + RL34$ + RL35$

WHILE RL$ <> "B"
GOSUB 1000
WEND
RL36$ = RL$
WHILE RL$ <> "E"
GOSUB 1000
WEND
RL37$ = RL$
WHILE RL$ <> "F"
GOSUB 1000
WEND
RL38$ = RL$
WHILE RL$ <> "O"
GOSUB 1000
WEND
RL39$ = RL$
WHILE RL$ <> "R"
GOSUB 1000
WEND
RL40$ = RL$
WHILE RL$ <> "E"
GOSUB 1000
WEND
RL41$ = RL$
PRINT #1, RL1$ + RL2$ + RL3$ + RL4$ + RL5$ + RL6$ + RL7$ + RL8$ + RL9$ + RL10$ + RL11$ + RL12$ + RL13$ + RL14$ + RL15$ + RL16$ + RL17$ + RL18$ + RL19$ + RL20$ + RL21$ + RL22$ + RL23$ + RL24$ + RL25$ + RL26$ + RL27$ + RL28$ + RL29$ + RL30$ + RL31$ + RL32$ + RL33$ + RL34$ + RL35$ + RL36$ + RL37$ + RL38$ + RL39$ + RL40$ + RL41$

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

