#!MC 1410
$!VarSet |MFBD| = '/home/pritam/Desktop/Github/Cpp/AdaptiveMethods/Output/Refinement'
$!VarSet |T| = 0
$!LOOP 100
$!PICK SETMOUSEMODE
  MOUSEMODE = SELECT
$!PAGE NAME = 'Untitled'
$!PAGECONTROL CREATE
$!PICK SETMOUSEMODE
  MOUSEMODE = SELECT
$!OPENLAYOUT  "|MFBD|/Template-Refinement.lpk"
$!READDATASET  '"|MFBD|/Field-|T|.tec" '
  READDATAOPTION = REPLACE
  RESETSTYLE = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  VARNAMELIST = '"X" "Y"'
$!PRINTSETUP PALETTE = COLOR
$!EXPORTSETUP EXPORTFORMAT = JPEG
$!EXPORTSETUP IMAGEWIDTH = 186
$!IF |T| <= 9
  $!EXPORTSETUP EXPORTFNAME = '|MFBD|/Frame-000|T|.jpeg'
$!ELSE
  $!IF |T| <= 99
    $!EXPORTSETUP EXPORTFNAME = '|MFBD|/Frame-00|T|.jpeg'
  $!ELSE
    $!EXPORTSETUP EXPORTFNAME = '|MFBD|/Frame-0|T|.jpeg'
  $!ENDIF
$!ENDIF
$!EXPORT 
  EXPORTREGION = CURRENTFRAME
$!VarSet |T|+=1
$!ENDLOOP
$!RemoveVar |MFBD|
