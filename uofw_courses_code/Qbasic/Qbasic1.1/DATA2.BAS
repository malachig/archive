REM This program enters data into an array and saves the data
REM to a permanent text file
CLS
OPEN "c:\qbasic\input2.txt" FOR OUTPUT AS #1
DATA 2,4,8,16,32,64,128,264,528,1056
DIM x(10)
FOR i = 1 TO 10
READ x(i)
PRINT x(i)
PRINT #1, x(i)
NEXT i
CLOSE #1
END


