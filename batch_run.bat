@ECHO OFF
FOR /L %%n IN (1,1,20) DO (
    ECHO
    ECHO ------------------
    ECHO   Run number: %%n
    ECHO ------------------
    C:\Python27\python.exe runmodel.py
)
