# PowerShell script to create virtual environment and install requirements
# Run from folder: c:\proga\REPEATABILITY\python

# 1) Check python
python --version
if ($LASTEXITCODE -ne 0) {
    Write-Host "Python not found. Install Python from https://www.python.org/downloads/ or use 'winget install --id Python.Python.3'." -ForegroundColor Yellow
    exit 1
}

# 2) Create virtual environment
python -m venv .venv

# 3) Ensure pip is up-to-date using the venv python directly
& .\.venv\Scripts\python.exe -m pip install --upgrade pip

# 4) Install requirements
& .\.venv\Scripts\python.exe -m pip install -r requirements.txt

# 5) Quick import test
& .\.venv\Scripts\python.exe -c "import pandas, numpy; import Bio; print('IMPORT_OK')"

Write-Host "Done. If you saw IMPORT_OK then environment is ready." -ForegroundColor Green
