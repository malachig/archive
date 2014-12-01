CLS
PRINT "The Molecular Weight of a Molecule"
PRINT ""
PRINT "Note: Use Capital Letters"
PRINT "To Finish Simply Type DONE"
PRINT ""

FOR i = 1 TO 10
INPUT "Atom"; ATOM$
IF ATOM$ = "DONE" GOTO 2000
GOSUB 1000
INPUT "Number"; Number
Total = MW * Number
sumTotal = sumTotal + Total
NEXT i
2000 PRINT ""
PRINT "The Molecular Weight of the Molecule ="; sumTotal
END

1000 IF ATOM$ = "H" THEN MW = 1.0079
IF ATOM$ = "C" THEN MW = 12.011
IF ATOM$ = "O" THEN MW = 15.9994
IF ATOM$ = "N" THEN MW = 14.0067
IF ATOM$ = "Na" THEN MW = 22.98977
IF ATOM$ = "K" THEN MW = 39.098
IF ATOM$ = "P" THEN MW = 30.97376
IF ATOM$ = "S" THEN MW = 32.06
IF ATOM$ = "Cl" THEN MW = 35.453
IF ATOM$ = "F" THEN MW = 18.9984
IF ATOM$ = "I" THEN MW = 126.9045
IF ATOM$ = "Si" THEN MW = 28.086
IF ATOM$ = "He" THEN MW = 4.0026
IF ATOM$ = "Mg" THEN MW = 24.305
IF ATOM$ = "Ca" THEN MW = 40.08
IF ATOM$ = "Li" THEN MW = 6.941
IF ATOM$ = "Al" THEN MW = 26.98154
IF ATOM$ = "Fe" THEN MW = 55.847
IF ATOM$ = "Cu" THEN MW = 63.546
IF ATOM$ = "Zn" THEN MW = 65.38
IF ATOM$ = "Ni" THEN MW = 58.69
RETURN







