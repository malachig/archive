CLS
OPEN "C:\QBASIC\ASS2\HUNGARY.TXT" FOR OUTPUT AS #1
PRINT #1, "A birth interval program using data from a Hungarian"
PRINT #1, "Population"
PRINT #1,
'The birth intervals 0-11, 12-17, 18-23, etc.
DATA 0,12,18,24,36,48,60

'The infant mortality per 1000
DATA 63.8,25.2,18.8,20.4,21.9,24.0,31.3

'The relative frequency of each interval
DATA .023,.13,.125,.201,.129,.104,.288

'Enter data into array in 3 columns
DIM array1(7, 4)
FOR j = 1 TO 3
FOR i = 1 TO 7
READ array1(i, j)
NEXT i
NEXT j

'Calculate infant mortality * the RF of each interval.
'This is like a weighted mean for the population by each interval.
'An extreme result that only has an occurance of one should not mean much.
'If there is only one occurence of a particular interval we should account
'for this when considering the population average.

'Column 2 * column 3 = column 4
FOR i = 1 TO 7
WMean = array1(i, 2) * array1(i, 3)
array1(i, 4) = WMean
sum2 = sum2 + WMean
NEXT i

PRINT #1, "Interval", "Interval Mortality Rate", "Relative Risk"
FOR i = 1 TO 7
risk = array1(i, 2) / 18.8
PRINT #1, array1(i, 1), array1(i, 2); , , risk
NEXT i
PRINT #1,
PRINT #1, "The weighted population mean ="; sum2

'Once we have the weighted mean for each segment of the population
'to get the total population mean we simply added these up.

'Relative risk for each interval = mortality for that interval
'divided by the mortality for the optimal, 13.5.
'Now we want to compare the mortality for each interval to the
'mortality for the optimal interval of 24-26 months.

'Since the average mortality is 25.5 deaths/1000 then if every one
'used the optimal interval with a mortality of 18.8 then 18.8/25.5
'accounts for only 74% of the deaths, so if everyone used this interval
'26% of the deaths could be avoided.
PRINT #1, "The proportion of deaths that could be avoided= 26%"



