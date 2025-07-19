@echo off
REM Build UDT manuscript PDF using Pandoc
REM Implements reviewer requirement for automated PDF generation

echo === UDT Manuscript PDF Builder ===
echo.

REM Check if pandoc is available
pandoc --version >nul 2>&1
if errorlevel 1 (
    echo Error: Pandoc not found. 
    echo Please install Pandoc from https://pandoc.org/installing.html
    echo.
    echo Alternative: View manuscript as Markdown in Universal_Distance_Dilation_Theory_Manuscript.md
    pause
    exit /b 1
)

REM Check if LaTeX is available (for PDF output)
pdflatex --version >nul 2>&1
if errorlevel 1 (
    echo Warning: pdflatex not found. Attempting HTML output instead...
    echo.
    echo Generating HTML manuscript...
    pandoc Universal_Distance_Dilation_Theory_Manuscript.md -o manuscript.html --defaults=pandoc-config.yaml --to=html --self-contained
    if errorlevel 1 (
        echo Error: Failed to generate HTML manuscript
        pause
        exit /b 1
    ) else (
        echo Success: manuscript.html generated
        echo You can open this file in any web browser
    )
) else (
    echo Generating PDF manuscript with full citations...
    echo.
    
    REM Generate PDF with bibliography
    pandoc Universal_Distance_Dilation_Theory_Manuscript.md -o manuscript.pdf --defaults=pandoc-config.yaml
    
    if errorlevel 1 (
        echo Error: Failed to generate PDF manuscript
        echo Trying simplified PDF generation...
        pandoc Universal_Distance_Dilation_Theory_Manuscript.md -o manuscript-simple.pdf --pdf-engine=pdflatex --number-sections
        if errorlevel 1 (
            echo Error: PDF generation failed completely
            echo Generating HTML as fallback...
            pandoc Universal_Distance_Dilation_Theory_Manuscript.md -o manuscript.html --number-sections --toc --self-contained
            echo Fallback HTML generated: manuscript.html
        ) else (
            echo Simplified PDF generated: manuscript-simple.pdf
        )
    ) else (
        echo Success: manuscript.pdf generated with full bibliography
    )
)

echo.
echo === Build Complete ===
echo.
echo Available outputs:
if exist "manuscript.pdf" echo   - manuscript.pdf (full version with citations)
if exist "manuscript-simple.pdf" echo   - manuscript-simple.pdf (simplified version)
if exist "manuscript.html" echo   - manuscript.html (web version)
echo   - Universal_Distance_Dilation_Theory_Manuscript.md (source markdown)
echo.

pause