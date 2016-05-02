The batch file should be edited as follows



If you wish to operate on rasters, use run_feature_id_rasters.bat and edit the variables accordingly

SET DSMPATH="path to DSM folder"
SET DTMPATH="path to DTM folder"
SET ZUNITS="Meter" or "Feet"
SET FEATS2EXTRACT="Buildings and Trees" or "Buildings only"
SET RESULTS="output folder" THIS MUST EXIST
SET PROCLEN=5000 
SET POST=true or false
SET DEL=false or true
SET REG=false or true

SET pathtoProPy="path to python.exe" Tend to use 64 bit version if running from CMD
SET pyscript="path to ExtractBuildingsTreesAutomation_v5.py"

%pathToProPy% %pyscript% %DSMPATH% %DTMPATH% %ZUNITS% %FEATS2EXTRACT% %RESULTS% %PROCLEN% %POST% %DEL% %REG%



If you wish to operate on a folder of .las files, use run_feature_id_las.bat and edit the variables accordingly

SET LASPATH="V:\Demo\las"
SET RETURNS="Last Returns" or "ALL Returns
SET ZUNITS="Meter" or "Feet"
SET FEATS2EXTRACT="Buildings and Trees" or "Buildings only"
SET RESULTS="output folder" THIS MUST EXIST
SET PROCLEN=5000 
SET POST=true or false
SET DEL=false or true
SET REG=false or true

SET pathtoProPy="path to python.exe" Tend to use 64 bit version if running from CMD
SET pyscript="path to ExtractBuildingsTreesAutomation_v5.py"

%pathToProPy% %pyscript% %LASPATH% %RETURNS% %ZUNITS% %FEATS2EXTRACT% %RESULTS% %PROCLEN% %POST% %DEL% %REG%