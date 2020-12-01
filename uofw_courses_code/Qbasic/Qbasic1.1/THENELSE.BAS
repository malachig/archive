CLS
INPUT "x="; x
INPUT "y="; y
IF x > 100 AND y < 10 THEN z = 0
IF x > 100 AND y >= 10 THEN z = 1
IF x < 100 THEN z = 2 ELSE z = 99
IF z = 0 THEN PRINT "Bart"
IF z = 1 THEN PRINT "lisa"
IF z = 2 THEN PRINT "Homer"
IF z = 99 THEN PRINT "Marge"
END


