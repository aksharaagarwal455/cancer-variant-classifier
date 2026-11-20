@echo off
echo ========================================================
echo   Cancer Variant Classifier - Web GUI
echo   URL: http://127.0.0.1:8000
echo ========================================================
python -m uvicorn app.api.main:app --app-dir backend --port 8000 --reload
pause
