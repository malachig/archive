REM This program simple reads data and prints it
CLS
PRINT "The data number", "The sum so far"
DATA 2,4,8,16,32,64,128,264,528,1056
FOR i = 1 TO 10
READ x
z = z + x
PRINT x, , z
NEXT i
END

