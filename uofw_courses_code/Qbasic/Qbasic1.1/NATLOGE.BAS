CLS
PRINT "The value of the natural log, e ="
PRINT
INPUT "Number of iterations"; n
FOR z# = 1 TO n
GOSUB 2000

x# = (1 / i#)
e# = e# + x#
NEXT z#
e# = e# + 1
PRINT "e="; e#
END


2000
i# = 1
FOR f# = 1 TO z#
i# = i# * f#
NEXT f#
RETURN

