CLS
OPEN "C:\QBASIC\ASS2\BANGDESH.TXT" FOR OUTPUT AS #1
PRINT #1, "Birth Interval Program for a Bangladesh Population"
PRINT #1,
'The birth intervals 0-19, 19-37, 38+.
DATA 0,20,38

'The infant mortality per 1000
DATA 151,101,64

'The frequency of each interval
DATA 765,2309,2109

'Enter data into array in 3 columns
DIM array1(3, 5)
FOR j = 1 TO 3
FOR i = 1 TO 3
READ array1(i, j)
NEXT i
NEXT j

'Calculate the relative frequency of each interval
sum1 = 0
FOR i = 1 TO 3
sum1 = sum1 + array1(i, 3)
NEXT i

'Enter RF values into array column 4
FOR i = 1 TO 3
RF = array1(i, 3) / sum1
array1(i, 4) = RF
NEXT i

'Calculate infant mortality * the RF of each interval.
'This is like a weighted mean for the population by each interval.
'An extreme result that only has an occurance of one should not mean much.
'If there is only one occurence of a particular interval we should account
'for this when considering the population average.

'Column 2 * column 4 = column 5
FOR i = 1 TO 3
WMean = array1(i, 2) * array1(i, 4)
array1(i, 5) = WMean
sum2 = sum2 + WMean
NEXT i

PRINT #1, "Interval", "Interval Mortality Rate", "Relative Risk"
FOR i = 1 TO 3
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
'mortality for the optimal interval of 38+ months.

'Since the average mortality is 93.3 deaths/1000 then if every one
'used the optimal interval with a mortality of 64 then 64/93.3
P = (64 / 93.3) * 100
PRINT #1, "The proportion of deaths which could be avoided ="; P; "%"
'accounts for only 17% of the deaths, so if everyone used this interval
'83% of the deaths could be avoided.




