CLS
b = 4
100 GOSUB 1000
110 GOSUB 2000
120 a = (b * c) / d
PRINT a
END

1000 c = SQR(b)
1010 RETURN

2000 d = b + b
2010 RETURN

