rmdir /s /q xmath_unit_test.vs2022
rem rmdir /s /q ..\dependencies
cmake ../ -G "Visual Studio 17 2022" -A x64 -B xmath_unit_test.vs2022
pause