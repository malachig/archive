CLS
OPEN "C:\QBASIC\ASS2\PAKISTAN.TXT" FOR OUTPUT AS #1
PRINT #1, "A birth interval program using data from a Pakistani"
PRINT #1, "Population."
PRINT #1,
'The birth intervals 0-11, 12-17, 18-23, etc.
DATA 0,12,18,24,30,36,42,48,54,64

'The infant mortality per 1000
DATA 258,200,151,120,101,75,95,81,85,64

'The frequency of each interval
DATA 712,2804,4262,5084,3590,2516,1234,788,390,1032

'Enter data into array in 3 columns
DIM array1(10, 5)
FOR j = 1 TO 3
FOR i = 1 TO 10
READ array1(i, j)
NEXT i
NEXT j

REM Calculate the relative frequency of each interval.
sum1 = 0
FOR i = 1 TO 10
sum1 = sum1 + array1(i, 3)
NEXT i

REM Enter RF values into array column four.
FOR i = 1 TO 10
RF = array1(i, 3) / sum1
array1(i, 4) = RF
NEXT i

'Calculate infant mortality * the RF of each interval.
'This is like a weighted mean for the population by each interval.
'An extreme result that only has an occurance of one should not mean much.
'If there is only one occurence of a particular interval we should account
'for this when considering the population average.

'Weighted Mean = Column 2 * column 4 = column 5

FOR i = 1 TO 10
WMean = array1(i, 2) * array1(i, 4)
array1(i, 5) = WMean
sum2 = sum2 + WMean
NEXT i

PRINT #1, "Interval", "Interval Mortality Rate", "Relative Risk"
FOR i = 1 TO 10
risk = array1(i, 2) / 64
PRINT #1, array1(i, 1), array1(i, 2); , , risk
NEXT i
PRINT #1,
PRINT #1, "The weighted population mean ="; sum2

'Once we have the weighted mean for each segment of the population
'to get the total population mean we simply added these up.

'Relative risk for each interval = mortality for that interval
'divided by the mortality for the optimal, 64.
'Now we want to compare the mortality for each interval to the
'mortality for the optimal interval of 24-26 months.

'Since the average mortality is 126.3 deaths/1000 then if every one
'used the optimal interval with a mortality of 64 then 64/126.3
'accounts for only 50.7% of the deaths, so if everyone used this interval
'49.3% of the deaths could be avoided.
PRINT #1, "The proportion of deaths that could be avoided = 49.3%"



