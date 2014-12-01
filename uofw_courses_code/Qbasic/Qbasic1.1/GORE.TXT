CLS
REM 1 = Gore
RANDOMIZE TIMER
Winner = 0
FOR j = 1 TO 10
Gore = 0
FOR i = 1 TO 10000
x = RND
IF x > .4989 THEN y = 1 ELSE y = 0
IF y = 1 THEN Gore = Gore + 1
NEXT i
IF Gore > 50000 THEN Winner = Winner + 1
NEXT j
PRINT Winner
END

