REM This program retrieves data from a permanent file and prints
REM it to the screen
CLS
OPEN "c:\qbasic\input2.txt" FOR INPUT AS #1
FOR i = 1 TO 10
INPUT #1, x
PRINT x
NEXT i
CLOSE #1
END


