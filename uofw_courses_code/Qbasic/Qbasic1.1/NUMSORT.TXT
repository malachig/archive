CLS
RANDOMIZE TIMER
n = 10
DIM x(n)
FOR i = 1 TO n
x(i) = INT(RND * 50) + 1: REM this generates a random number and places it in an array variable
NEXT i
GOSUB 1000: REM Prints out unsorted numbers
FOR i = 1 TO n
FOR k = i + 1 TO n
IF x(i) > x(k) THEN
SWAP x(i), x(k)
END IF
NEXT k
NEXT i
GOSUB 1000: REM Prints out sorted numbers
END

1000 REM This is the printing routine
FOR i = 1 TO n
PRINT x(i)
NEXT i
PRINT
RETURN

