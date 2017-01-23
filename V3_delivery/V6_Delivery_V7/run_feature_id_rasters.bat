SET DSMPATH="V:\Demo\DSM_subset"
SET DTMPATH="V:\Demo\DTM_subset"
SET ZUNITS="Meter"
SET FEATS2EXTRACT="Buildings and Trees"
SET RESULTS="V:\Demo\utah_demo_pro_debug1"
SET PROCLEN=5000
SET POST=true
SET DEL=false
SET REG=false

SET pathtoProPy="C:\Python34\python.exe"
SET pyscript="V:\Demo\V3_Delivery\ExtractBuildingsTreesAutomation_V5.py"

%pathToProPy% %pyscript% %DSMPATH% %DTMPATH% %ZUNITS% %FEATS2EXTRACT% %RESULTS% %PROCLEN% %POST% %DEL% %REG%



