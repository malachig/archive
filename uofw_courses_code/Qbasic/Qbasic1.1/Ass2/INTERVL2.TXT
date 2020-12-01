CLS
OPEN "C:\QBASIC\ASS2\MEXICO.TXT" FOR OUTPUT AS #1
PRINT #1, "A birth interval program for a Mexican population"
PRINT #1,
'The birth intervals 0-19, 19-37, 37+
DATA 0,19,37

'The infant mortality per 1000
DATA 93,69,30

'Enter data into array in columns
DIM Array1(3, 3)
FOR j = 1 TO 2
FOR i = 1 TO 3
READ Array1(i, j)
NEXT i
NEXT j

'Calculate the population mean
FOR i = 1 TO 3
sum1 = sum1 + Array1(i, 2)
NEXT i
mean = sum1 / 3
PRINT #1, "The population mean mortality="; mean
PRINT #1,

'Calculate the relative risk
FOR i = 1 TO 3
risk = Array1(i, 2) / 30
Array1(i, 3) = risk
NEXT i

PRINT #1, "The birth interval", "Mortality of Interval", "Relative Risk"
FOR i = 1 TO 3
PRINT #1, Array1(i, 1), , Array1(i, 2), , Array1(i, 3)
NEXT i
PRINT #1,
'Relative risk for each interval = mortality for that interval
'divided by the mortality for the optimal, 30.
'Now we want to compare the mortality for each interval to the
'mortality for the optimal interval of 37+ months.

'Since the average mortality is 64 deaths/1000 then if every one
'used the optimal interval with a mortality of 30 then 30/64
a = (30 / 64) * 100
PRINT #1, "Proportion of deaths that could be avoided ="; a; "%"
'accounts for only 46.9% of the deaths, so if everyone used this interval
'53.1% of the deaths could be avoided.




