DECLARE SUB Newfeat ()
'
'                         Q B a s i c   N i b b l e s
'
'                   Copyright (C) Microsoft Corporation 1990
'
' Nibbles is a game for one or two players.  Navigate your snakes
' around the game board trying to eat up numbers while avoiding
' running into walls or other snakes.  The more numbers you eat up,
' the more points you gain and the longer your snake becomes.
'
' To run this game, press Shift+F5.
'
' To exit QBasic, press Alt, F, X.
'
' To get help on a BASIC keyword, move the cursor to the keyword and press
' F1 or click the right mouse button.
'

'Set default data type to integer for faster game play
DEFINT A-Z

'User-defined TYPEs
TYPE snakeBody
    row AS INTEGER
    col AS INTEGER
END TYPE

'This type defines the player's snake
TYPE snaketype
    head      AS INTEGER
    length    AS INTEGER
    row       AS INTEGER
    col       AS INTEGER
    direction AS INTEGER
    lives     AS INTEGER
    score     AS INTEGER
    scolor    AS INTEGER
    alive     AS INTEGER
END TYPE

'This type is used to represent the playing screen in memory
'It is used to simulate graphics in text mode, and has some interesting,
'and slightly advanced methods to increasing the speed of operation.
'Instead of the normal 80x25 text graphics using chr$(219) "€", we will be
'using chr$(220)"‹" and chr$(223) "ﬂ" and chr$(219) "€" to mimic an 80x50
'pixel screen.
'Check out sub-programs SET and POINTISTHERE to see how this is implemented
'feel free to copy these (as well as arenaType and the DIM ARENA stmt and the
'initialization code in the DrawScreen subprogram) and use them in your own
'programs
TYPE arenaType
    realRow     AS INTEGER        'Maps the 80x50 point into the real 80x25
    acolor      AS INTEGER        'Stores the current color of the point
    sister      AS INTEGER        'Each char has 2 points in it.  .SISTER is
END TYPE                          '-1 if sister point is above, +1 if below

'Sub Declarations
DECLARE SUB SpacePause (text$)
DECLARE SUB PrintScore (NumPlayers%, score1%, score2%, score3%, lives1%, lives2%, lives3%)
DECLARE SUB Intro ()
DECLARE SUB GetInputs (NumPlayers, speed, diff$, monitor$)
DECLARE SUB DrawScreen ()
DECLARE SUB PlayNibbles (NumPlayers, speed, diff$)
DECLARE SUB Set (row, col, acolor)
DECLARE SUB center (row, text$)
DECLARE SUB DoIntro ()
DECLARE SUB Initialize ()
DECLARE SUB SparklePause ()
DECLARE SUB Level (WhatToDO, sammy() AS snaketype)
DECLARE SUB InitColors ()
DECLARE SUB EraseSnake (snake() AS ANY, snakeBod() AS ANY, snakeNum%)
DECLARE FUNCTION StillWantsToPlay ()
DECLARE FUNCTION PointIsThere (row, col, backColor)

'Constants
CONST TRUE = -1
CONST FALSE = NOT TRUE
CONST MAXSNAKELENGTH = 1000
CONST STARTOVER = 1             ' Parameters to 'Level' SUB
CONST SAMELEVEL = 2
CONST NEXTLEVEL = 3

'Global Variables
DIM SHARED arena(1 TO 50, 1 TO 80) AS arenaType
DIM SHARED curLevel, colorTable(11)

    RANDOMIZE TIMER
'    GOSUB ClearKeyLocks
    Intro
    GetInputs NumPlayers, speed, diff$, monitor$
    GOSUB SetColors
    DrawScreen

    DO
      PlayNibbles NumPlayers, speed, diff$
    LOOP WHILE StillWantsToPlay

    GOSUB RestoreKeyLocks
    COLOR 15, 0
    CLS
END

ClearKeyLocks:
    DEF SEG = 0                     ' Turn off CapLock, NumLock and ScrollLock
    KeyFlags = PEEK(1047)
    POKE 1047, &H0
    DEF SEG
    RETURN

RestoreKeyLocks:
    DEF SEG = 0                     ' Restore CapLock, NumLock and ScrollLock states
    POKE 1047, KeyFlags
    DEF SEG
    RETURN

SetColors:
    IF monitor$ = "M" THEN
        RESTORE mono
    ELSE
        RESTORE normal
    END IF

    FOR a = 1 TO 6
        READ colorTable(a)
    NEXT a
    RETURN

           'snake1     snake2   snake3   Walls  Background  Dialogs-Fore  Back
mono:   DATA 15,         7,       2,       7,     0,          15,            0
normal: DATA 14,         13,      10,      12,    1,          15,            4
END

'Center:
'  Centers text on given row
SUB center (row, text$)
    LOCATE row, 41 - LEN(text$) / 2
    PRINT text$;
END SUB

'DrawScreen:
'  Draws playing field
SUB DrawScreen

    'initialize screen
    VIEW PRINT
    COLOR colorTable(1), colorTable(5)
    CLS

    'Print title & message
    center 2, "Nibbles!"
    center 11, "Initializing arena..."
    
    'Initialize arena array
    FOR row = 1 TO 50
        FOR col = 1 TO 80
            arena(row, col).realRow = INT((row + 1) / 2)
            arena(row, col).sister = (row MOD 2) * 2 - 1
        NEXT col
    NEXT row
END SUB

'EraseSnake:
'  Erases snake to facilitate moving through playing field
SUB EraseSnake (snake() AS snaketype, snakeBod() AS snakeBody, snakeNum)

    FOR c = 0 TO 9
        FOR b = snake(snakeNum).length - c TO 0 STEP -10
            tail = (snake(snakeNum).head + MAXSNAKELENGTH - b) MOD MAXSNAKELENGTH
            Set snakeBod(tail, snakeNum).row, snakeBod(tail, snakeNum).col, colorTable(5)
        NEXT b
    NEXT c
    
END SUB

'GetInputs:
'  Gets player inputs
SUB GetInputs (NumPlayers, speed, diff$, monitor$)

    COLOR 7, 0
    CLS

    DO
        LOCATE 4, 47: PRINT SPACE$(34);
        LOCATE 4, 20
        INPUT "How many players (1, 2 or 3)"; num$
    LOOP UNTIL VAL(num$) = 1 OR VAL(num$) = 2 OR VAL(num$) = 3 OR VAL(num$) = q
    NumPlayers = VAL(num$)
        IF VAL(num$) = q THEN
        END
        END IF
        IF VAL(num$) = 3 THEN
    DO
        LOCATE 6, 2
        INPUT "Game 1 (2 survivors) or 2 (1 Survivor)"; gameno$
        LOOP UNTIL VAL(gameno$) = 1 OR VAL(gameno$) = 2
        END IF
    LOCATE 8, 21: PRINT "Skill level (1 to 100)"
    LOCATE 9, 22: PRINT "1   = Laughable, yet quietly pathetic"
    LOCATE 10, 22: PRINT "90  = Zippier than a very zippy thing"
    LOCATE 11, 22: PRINT "100 = Whoa! Twiddle fingers!"
    LOCATE 12, 15: PRINT "(Computer speed may affect your skill level, especially with 3 players)"
    DO
        LOCATE 8, 44: PRINT SPACE$(35);
        LOCATE 8, 43
        INPUT gamespeed$
    LOOP UNTIL VAL(gamespeed$) >= 1 AND VAL(gamespeed$) <= 100
    speed = VAL(gamespeed$)
  
    speed = (100 - speed) * 2 + 1

    DO
        LOCATE 15, 56: PRINT SPACE$(25);
        LOCATE 15, 15
        INPUT "Increase game speed during play (Y or N)"; diff$
        diff$ = UCASE$(diff$)
    LOOP UNTIL diff$ = "Y" OR diff$ = "N"

    DO
        LOCATE 17, 46: PRINT SPACE$(34);
        LOCATE 17, 17
        INPUT "Monochrome or color monitor (M or C)"; monitor$
        monitor$ = UCASE$(monitor$)
    LOOP UNTIL monitor$ = "M" OR monitor$ = "C" OR monitor$ = ""
    IF monitor$ = "" THEN monitor$ = "C"

    startTime# = TIMER                          ' Calculate speed of system
    FOR I# = 1 TO 1000: NEXT I#                 ' and do some compensation
    stopTime# = TIMER
    speed = speed * .5 / (stopTime# - startTime#)

END SUB

'InitColors:
'Initializes playing field colors
SUB InitColors
    
    FOR row = 1 TO 50
        FOR col = 1 TO 80
            arena(row, col).acolor = colorTable(5)
        NEXT col
    NEXT row

    CLS
   
    'Set (turn on) pixels for screen border
    FOR col = 1 TO 80
        Set 3, col, colorTable(4)
        Set 50, col, colorTable(4)
    NEXT col

    FOR row = 4 TO 49
        Set row, 1, colorTable(4)
        Set row, 80, colorTable(4)
    NEXT row

END SUB

'Intro:
'  Displays game introduction
SUB Intro
    SCREEN 0
    WIDTH 80, 25
    COLOR 15, 0
    CLS

    center 3, "Q B a s i c   N i b b l e s   M k . I I"
    COLOR 7
    center 5, "Copyright (C) Microsoft Corporation 1990"
    center 6, ""
    PRINT "    Nibbles is one groovy game for up to THREE players, so get all your"
    PRINT "  friends round and enough beer for the month. Whizz your snakes around"
    PRINT "  the arena trying to eat up numbers while avoiding running into walls"
    PRINT "  or other snakes.  The more numbers you eat up, the more points you"
    PRINT "  gain and the longer your snake becomes. Master this game."
    center 12, "All thanks for the 3 player addition must go to Sammy. &B-)"
    center 14, " Game Controls: "
    PRINT
    COLOR 14
    LOCATE 16, 13: PRINT "Player 1"
    COLOR 13
    LOCATE 16, 36: PRINT "Player 2":
    COLOR 10
    LOCATE 16, 59: PRINT "Player 3"
    COLOR 7
    PRINT "               " + CHR$(24) + "                      W                      8"
    PRINT "             " + CHR$(27) + "   " + CHR$(26) + "                  A   D                  4   6"
    PRINT "               " + CHR$(25) + "                      S                      2   (Num Pad)"
    center 21, "In game: P for Pause, Q for quit"
    center 25, "Let battle commence..."

    PLAY "MBT160O1L8CDEDCDL4ECC"
    SparklePause

END SUB

'Level:
'Sets game level
SUB Level (WhatToDO, sammy() AS snaketype) STATIC
    
    SELECT CASE (WhatToDO)

    CASE STARTOVER
        curLevel = 1
    CASE NEXTLEVEL
        curLevel = curLevel + 1
    END SELECT

    sammy(1).head = 1                       'Initialize Snakes
    sammy(1).length = 2
    sammy(1).alive = TRUE
    sammy(2).head = 1
    sammy(2).length = 2
    sammy(2).alive = TRUE
    sammy(3).head = 1
    sammy(3).length = 2
    sammy(3).alive = TRUE
    
    InitColors
    
    SELECT CASE curLevel
    CASE 1
        sammy(1).row = 25: sammy(2).row = 25: sammy(3).row = 25
        sammy(1).col = 50: sammy(2).col = 30: sammy(3).col = 40
        sammy(1).direction = 4: sammy(2).direction = 3: sammy(3).direction = 2


    CASE 2
        FOR I = 20 TO 60
            Set 25, I, colorTable(4)
        NEXT I
        sammy(1).row = 7: sammy(2).row = 43: sammy(3).row = 43
        sammy(1).col = 60: sammy(2).col = 20: sammy(3).col = 60
        sammy(1).direction = 3: sammy(2).direction = 4: sammy(3).direction = 3

    CASE 3
        FOR I = 10 TO 40
            Set I, 20, colorTable(4)
            Set I, 60, colorTable(4)
        NEXT I
        sammy(1).row = 25: sammy(2).row = 25: sammy(3).row = 27
        sammy(1).col = 50: sammy(2).col = 30: sammy(3).col = 40
        sammy(1).direction = 1: sammy(2).direction = 2: sammy(3).direction = 4

    CASE 4
        FOR I = 4 TO 30
            Set I, 20, colorTable(4)
            Set 53 - I, 60, colorTable(4)
        NEXT I
        FOR I = 2 TO 40
            Set 38, I, colorTable(4)
            Set 15, 81 - I, colorTable(4)
        NEXT I
        sammy(1).row = 7: sammy(2).row = 43: sammy(3).row = 20
        sammy(1).col = 60: sammy(2).col = 20: sammy(3).col = 4
        sammy(1).direction = 3: sammy(2).direction = 4: sammy(3).direction = 2
   
    CASE 5
        FOR I = 13 TO 39
            Set I, 21, colorTable(4)
            Set I, 59, colorTable(4)
        NEXT I
        FOR I = 23 TO 57
            Set 11, I, colorTable(4)
            Set 41, I, colorTable(4)
        NEXT I
        sammy(1).row = 25: sammy(2).row = 25: sammy(3).row = 27
        sammy(1).col = 50: sammy(2).col = 30: sammy(3).col = 40
        sammy(1).direction = 1: sammy(2).direction = 2: sammy(3).direction = 4

    CASE 6
        FOR I = 4 TO 49
            IF I > 30 OR I < 23 THEN
                Set I, 10, colorTable(4)
                Set I, 20, colorTable(4)
                Set I, 30, colorTable(4)
                Set I, 40, colorTable(4)
                Set I, 50, colorTable(4)
                Set I, 60, colorTable(4)
                Set I, 70, colorTable(4)
            END IF
        NEXT I
        sammy(1).row = 7: sammy(2).row = 43: sammy(3).row = 7
        sammy(1).col = 65: sammy(2).col = 15: sammy(3).col = 40
        sammy(1).direction = 2: sammy(2).direction = 1: sammy(3).direction = 2
   
    CASE 7
        FOR I = 6 TO 47
            Set I, I, colorTable(4)
            Set I, I + 28, colorTable(4)
        NEXT I
        sammy(1).row = 40: sammy(2).row = 15: sammy(3).row = 28
        sammy(1).col = 75: sammy(2).col = 5: sammy(3).col = 15
        sammy(1).direction = 1: sammy(2).direction = 2: sammy(3).direction = 4
      
    CASE 8
        FOR I = 4 TO 49 STEP 2
            Set I, 40, colorTable(4)
        NEXT I
        sammy(1).row = 7: sammy(2).row = 43: sammy(3).row = 6
        sammy(1).col = 65: sammy(2).col = 15: sammy(3).col = 40
        sammy(1).direction = 2: sammy(2).direction = 1: sammy(3).direction = 4

    
    CASE 9
        FOR I = 4 TO 40
            Set I, 10, colorTable(4)
            Set 53 - I, 20, colorTable(4)
            Set I, 30, colorTable(4)
            Set 53 - I, 40, colorTable(4)
            Set I, 50, colorTable(4)
            Set 53 - I, 60, colorTable(4)
            Set I, 70, colorTable(4)
        NEXT I
        sammy(1).row = 7: sammy(2).row = 43: sammy(3).row = 7
        sammy(1).col = 65: sammy(2).col = 15: sammy(3).col = 38
        sammy(1).direction = 2: sammy(2).direction = 1: sammy(3).direction = 2

    CASE ELSE
        FOR I = 4 TO 49 STEP 2
            Set I, 10, colorTable(4)
            Set I + 1, 20, colorTable(4)
            Set I, 30, colorTable(4)
            Set I + 1, 40, colorTable(4)
            Set I, 50, colorTable(4)
            Set I + 1, 60, colorTable(4)
            Set I, 70, colorTable(4)
        NEXT I
        sammy(1).row = 7: sammy(2).row = 43: sammy(3).row = 7
        sammy(1).col = 65: sammy(2).col = 15: sammy(3).col = 40
        sammy(1).direction = 2: sammy(2).direction = 1: sammy(3).direction = 2

    END SELECT
END SUB

'Newfeatures:
'Displays new features added by Sammy the man
    SUB Newfeat
    SCREEN 0
    WIDTH 80, 25
    COLOR 15, 0
    CLS

    center 3, "Nibbles Mk. II - New Features"
    COLOR 7
    center 5, "Added 1995, by Sammy the Man " + CHR$(2)
    center 6, ""
    PRINT "    The original Nibbles was good, but (hopefully) this is better."
    PRINT "  All the original bits are in here, I've just added some extra"
    PRINT "  stuff in. Some are just to make life easier, others are there to"
    PRINT "  make it fun. 3 players especially. So just give it a go..."
    PRINT ""
    PRINT "  So, the new features then:"
    PRINT "       3rd player option, quit when press q, Improved intro"
    PRINT "  screen, this bit, change of snake names."
    PRINT ""
    PRINT "  Hopefully features to be added:"
    PRINT "       Improved 3 player option - 1 survivor/match, more levels."
    PRINT " "
    SparklePause

END SUB

'PlayNibbles:
'  Main routine that controls game play
SUB PlayNibbles (NumPlayers, speed, diff$)

    'Initialize Snakes
    DIM sammyBody(MAXSNAKELENGTH - 1, 1 TO 3) AS snakeBody
    DIM sammy(1 TO 3) AS snaketype
    sammy(1).lives = 5
    sammy(1).score = 0
    sammy(1).scolor = colorTable(1)
    sammy(2).lives = 5
    sammy(2).score = 0
    sammy(2).scolor = colorTable(2)
    sammy(3).lives = 5
    sammy(3).score = 0
    sammy(3).scolor = colorTable(3)

    Level STARTOVER, sammy()
    startRow1 = sammy(1).row: startCol1 = sammy(1).col
    startRow2 = sammy(2).row: startCol2 = sammy(2).col
    startRow3 = sammy(3).row: startCol3 = sammy(3).col

    curSpeed = speed

    'play Nibbles until finished

    SpacePause "  Enter arena" + STR$(curLevel) + ", Push Space  "
    gameOver = FALSE
    DO
        IF NumPlayers = 1 THEN
            sammy(2).row = 0
            sammy(3).row = 0
        END IF

        IF NumPlayers = 2 THEN
            sammy(3).row = 0
        END IF

        number = 1          'Current number that snakes are trying to run into
        nonum = TRUE        'nonum = TRUE if a number is not on the screen
                
        playerDied = FALSE
        PrintScore NumPlayers, sammy(1).score, sammy(2).score, sammy(3).score, sammy(1).lives, sammy(2).lives, sammy(3).lives
        PLAY "T160O1>L20CDEDCDL10ECC"

        DO
            'Print number if no number exists
            IF nonum = TRUE THEN
                DO
                    numberRow = INT(RND(1) * 47 + 3)
                    NumberCol = INT(RND(1) * 78 + 2)
                    sisterRow = numberRow + arena(numberRow, NumberCol).sister
                LOOP UNTIL NOT PointIsThere(numberRow, NumberCol, colorTable(5)) AND NOT PointIsThere(sisterRow, NumberCol, colorTable(5))
                numberRow = arena(numberRow, NumberCol).realRow
                nonum = FALSE
                COLOR colorTable(1), colorTable(5)
                LOCATE numberRow, NumberCol
                PRINT RIGHT$(STR$(number), 1);
                count = 0
            END IF

            'Delay game
            FOR a# = 1 TO curSpeed:  NEXT a#

            'Get keyboard input & Change direction accordingly
            kbd$ = INKEY$
            SELECT CASE kbd$
                CASE "w", "W": IF sammy(2).direction <> 2 THEN sammy(2).direction = 1
                CASE "s", "S": IF sammy(2).direction <> 1 THEN sammy(2).direction = 2
                CASE "a", "A": IF sammy(2).direction <> 4 THEN sammy(2).direction = 3
                CASE "d", "D": IF sammy(2).direction <> 3 THEN sammy(2).direction = 4
                CASE CHR$(0) + "H": IF sammy(1).direction <> 2 THEN sammy(1).direction = 1
                CASE CHR$(0) + "P": IF sammy(1).direction <> 1 THEN sammy(1).direction = 2
                CASE CHR$(0) + "K": IF sammy(1).direction <> 4 THEN sammy(1).direction = 3
                CASE CHR$(0) + "M": IF sammy(1).direction <> 3 THEN sammy(1).direction = 4
                CASE "8": IF sammy(3).direction <> 2 THEN sammy(3).direction = 1
                CASE "5": IF sammy(3).direction <> 1 THEN sammy(3).direction = 2
                CASE "4": IF sammy(3).direction <> 4 THEN sammy(3).direction = 3
                CASE "6": IF sammy(3).direction <> 3 THEN sammy(3).direction = 4
                CASE "p", "P": SpacePause " Game Paused ... Push Space  "
                CASE "q", "Q": LOCATE 11, 27: INPUT "  Really quit? (y/n)   ", quit$
                    IF quit$ = "y" THEN
                    CLS
                    END
                    ELSE
                    FOR I = 10 TO 12
                        FOR j = 26 TO 52
                            Set I, j, arena(I, j).acolor
                        NEXT j
                    NEXT I
                    SpacePause " Push Space to continue  "
                    END IF
                CASE ELSE
            END SELECT

            FOR a = 1 TO NumPlayers
                'Move Snake
                SELECT CASE sammy(a).direction
                    CASE 1: sammy(a).row = sammy(a).row - 1
                    CASE 2: sammy(a).row = sammy(a).row + 1
                    CASE 3: sammy(a).col = sammy(a).col - 1
                    CASE 4: sammy(a).col = sammy(a).col + 1
                END SELECT

                'If snake hits number, respond accordingly
                IF numberRow = INT((sammy(a).row + 1) / 2) AND NumberCol = sammy(a).col THEN
                    PLAY "MBO0L16>CCCE"
                    IF sammy(a).length < (MAXSNAKELENGTH - 30) THEN
                        sammy(a).length = sammy(a).length + number * 4
                    END IF
                    sammy(a).score = sammy(a).score + number
                    
                    PrintScore NumPlayers, sammy(1).score, sammy(2).score, sammy(3).score, sammy(1).lives, sammy(2).lives, sammy(3).lives
                    number = number + 1
                    IF number = 10 THEN
                        EraseSnake sammy(), sammyBody(), 1
                        EraseSnake sammy(), sammyBody(), 2
                        EraseSnake sammy(), sammyBody(), 3
                        LOCATE numberRow, NumberCol: PRINT " "
                        Level NEXTLEVEL, sammy()
                        PrintScore NumPlayers, sammy(1).score, sammy(2).score, sammy(3).score, sammy(1).lives, sammy(2).lives, sammy(3).lives
                        SpacePause "  Enter Arena" + STR$(curLevel) + ",  Push Space"
                        IF NumPlayers = 1 THEN sammy(2).row = 0 AND sammy(3).row = 0
                        IF NumPlayers = 2 THEN sammy(3).row = 0
                        number = 1
                        IF diff$ = "P" THEN speed = speed - 10: curSpeed = speed
                    END IF
                    nonum = TRUE
                    IF curSpeed < 1 THEN curSpeed = 1
                END IF
            NEXT a

            FOR a = 1 TO NumPlayers
                'If player runs into any point, or the head of the other snake, it dies
IF PointIsThere(sammy(a).row, sammy(a).col, colorTable(5)) OR (sammy(1).row = sammy(2).row AND sammy(1).col = sammy(2).col) OR (sammy(1).row = sammy(3).row AND sammy(1).col = sammy(3).col) OR (sammy(2).row = sammy(3).row AND sammy(2).col = sammy(3).col) THEN
                    PLAY "MBO0L32EFGEFDC"
                    COLOR , colorTable(5)
                    LOCATE numberRow, NumberCol
                    PRINT " "
                   
                    playerDied = TRUE
                    sammy(a).alive = FALSE
                    sammy(a).lives = sammy(a).lives - 1
'----            
'            FOR a = 1 TO NumPlayers
'                'If player runs into any point, or the head of the other snake, it dies
'IF PointIsThere(sammy(a).row, sammy(a).col, colorTable(5)) OR (sammy(1).row = sammy(2).row AND sammy(1).col = sammy(2).col) OR (sammy(1).row = sammy(3).row AND sammy(1).col = sammy(3).col) OR (sammy(2).row = sammy(3).row AND sammy(2).col = sammy(3).col) THEN
'                    PLAY "MBO0L32EFGEFDC"
'                    COLOR , colorTable(5)
'                    LOCATE numberRow, NumberCol
'                    PRINT " "
'             
'                    playerDied = TRUE
'----
                'Otherwise, move the snake, and erase the tail
                ELSE
                    sammy(a).head = (sammy(a).head + 1) MOD MAXSNAKELENGTH
                    sammyBody(sammy(a).head, a).row = sammy(a).row
                    sammyBody(sammy(a).head, a).col = sammy(a).col
                    tail = (sammy(a).head + MAXSNAKELENGTH - sammy(a).length) MOD MAXSNAKELENGTH
                    Set sammyBody(tail, a).row, sammyBody(tail, a).col, colorTable(5)
                    sammyBody(tail, a).row = 0
                    Set sammy(a).row, sammy(a).col, sammy(a).scolor
                END IF
            NEXT a

        LOOP UNTIL playerDied

        curSpeed = speed                ' reset speed to initial value
       
        FOR a = 1 TO NumPlayers
            EraseSnake sammy(), sammyBody(), a

            'If dead, then erase snake in really cool way
            IF sammy(a).alive = FALSE THEN
                'Update score
                sammy(a).score = sammy(a).score - 10
                PrintScore NumPlayers, sammy(1).score, sammy(2).score, sammy(3).score, sammy(1).lives, sammy(2).lives, sammy(3).lives
                
                IF a = 1 THEN
                    SpacePause "     Sammy dies! Happy now? "
                ELSE
                    IF a = 2 THEN
                    SpacePause "     Jake Dies! Satisfied? "
                ELSE
                    SpacePause "      Dez dies! Pleased?  "
                    END IF
                END IF
            END IF
        NEXT a

        Level SAMELEVEL, sammy()
        PrintScore NumPlayers, sammy(1).score, sammy(2).score, sammy(3).score, sammy(1).lives, sammy(2).lives, sammy(3).lives
     
    'Play next round, until either of snake's lives have run out.
    LOOP UNTIL sammy(1).lives = 0 OR sammy(2).lives = 0 OR sammy(3).lives = 0

END SUB

'PointIsThere:
'  Checks the global  arena array to see if the boolean flag is set
FUNCTION PointIsThere (row, col, acolor)
    IF row <> 0 THEN
        IF arena(row, col).acolor <> acolor THEN
            PointIsThere = TRUE
        ELSE
            PointIsThere = FALSE
        END IF
    END IF
END FUNCTION

'PrintScore:
'  Prints players scores and number of lives remaining
SUB PrintScore (NumPlayers, score1, score2, score3, lives1, lives2, lives3)
    COLOR 15, colorTable(5)
   
    IF NumPlayers = 3 THEN
        LOCATE 1, 33
        PRINT USING "Dez:  Lives: #     #,###,#00"; lives3; score3
       
        LOCATE 1, 50
        PRINT USING "Jake-->  Lives: #     #,###,#00"; lives2; score2
    END IF

    IF NumPlayers = 2 THEN
        LOCATE 1, 50
        PRINT USING "Jake-->  Lives: #     #,###,#00"; lives2; score2
    END IF

    LOCATE 1, 1
    PRINT USING "#,###,#00  Lives: #  <--Sammy"; score1; lives1

END SUB

'Set:
'  Sets row and column on playing field to given color to facilitate moving
'  of snakes around the field.
SUB Set (row, col, acolor)
    IF row <> 0 THEN
        arena(row, col).acolor = acolor             'assign color to arena
        realRow = arena(row, col).realRow           'Get real row of pixel
        topFlag = arena(row, col).sister + 1 / 2    'Deduce whether pixel
                                                    'is on topﬂ, or bottom‹
        sisterRow = row + arena(row, col).sister    'Get arena row of sister
        sisterColor = arena(sisterRow, col).acolor  'Determine sister's color

        LOCATE realRow, col

        IF acolor = sisterColor THEN                'If both points are same
            COLOR acolor, acolor                           'Print chr$(219) "€"
            PRINT CHR$(219);
        ELSE
            IF topFlag THEN                         'Since you cannot have
                IF acolor > 7 THEN                  'bright backgrounds
                    COLOR acolor, sisterColor       'determine best combo
                    PRINT CHR$(223);                'to use.
                ELSE
                    COLOR sisterColor, acolor
                    PRINT CHR$(220);
                END IF
            ELSE
                IF acolor > 7 THEN
                    COLOR acolor, sisterColor
                    PRINT CHR$(220);
                ELSE
                    COLOR sisterColor, acolor
                    PRINT CHR$(223);
                END IF
            END IF
        END IF
    END IF
END SUB

'SpacePause:
'  Pauses game play and waits for space bar to be pressed before continuing
SUB SpacePause (text$)

    COLOR colorTable(6), colorTable(7)
    center 11, "€ﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂ€"
    center 12, "€ " + LEFT$(text$ + SPACE$(29), 29) + " €"
    center 13, "€‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹€"
    WHILE INKEY$ <> "": WEND
    WHILE INKEY$ <> " ": WEND
    COLOR 15, colorTable(5)

    FOR I = 21 TO 26            ' Restore the screen background
        FOR j = 24 TO 56
            Set I, j, arena(I, j).acolor
        NEXT j
    NEXT I

END SUB

'SparklePause:
'  Creates flashing border for intro screen
SUB SparklePause
    
    COLOR 4, 0
    a$ = "*    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    "
    IF INKEY$ = "n" OR INKEY$ = "N" THEN
    Newfeat
    END IF
    WHILE INKEY$ <> "": WEND   'Clear keyboard buffer

    WHILE INKEY$ = ""
        FOR a = 1 TO 5
            LOCATE 1, 1                             'print horizontal sparkles
            PRINT MID$(a$, a, 80);
            LOCATE 22, 1
            PRINT MID$(a$, 6 - a, 80);

            FOR b = 2 TO 21                         'Print Vertical sparkles
                c = (a + b) MOD 5
                IF c = 1 THEN
                    LOCATE b, 80
                    PRINT "*";
                    LOCATE 23 - b, 1
                    PRINT "*";
                ELSE
                    LOCATE b, 80
                    PRINT " ";
                    LOCATE 23 - b, 1
                    PRINT " ";
                END IF
            NEXT b
        NEXT a
    WEND

END SUB

'StillWantsToPlay:
'  Determines if users want to play game again.
FUNCTION StillWantsToPlay

    COLOR colorTable(6), colorTable(7)
    center 10, "€ﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂﬂ€"
    center 11, "€       G A M E   O V E R       €"
    center 12, "€                               €"
    center 13, "€        Rematch?  (Y/N)        €"
    center 14, "€‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹€"

    WHILE INKEY$ <> "": WEND
    DO
        kbd$ = UCASE$(INKEY$)
    LOOP UNTIL kbd$ = "Y" OR kbd$ = "N"

    COLOR 15, colorTable(5)
    center 10, "                                 "
    center 11, "                                 "
    center 12, "                                 "
    center 13, "                                 "
    center 14, "                                 "

    IF kbd$ = "Y" THEN
        StillWantsToPlay = TRUE
    ELSE
        StillWantsToPlay = FALSE
        COLOR 7, 0
        CLS
    END IF

END FUNCTION

