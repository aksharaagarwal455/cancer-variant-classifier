@echo off
echo Starting Cancer Variant Classifier GUI...
echo Access at: http://127.0.0.1:8000
python -m uvicorn app.api.main:app --app-dir backend --port 8000 --reload
pause
