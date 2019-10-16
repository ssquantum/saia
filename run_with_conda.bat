rem -- set conda and python as environment variables. Check if python is installed for all users first.
if EXIST %ALLUSERSPROFILE%\Anaconda3 ( 
	call %ALLUSERSPROFILE%\Anaconda3\Scripts\activate
) else ( 
	call %USERPROFILE%\Anaconda3\Scripts\activate
)


rem -- activate the Python environment to load modules
call conda activate saiaenv

rem -- now start SAIA
python main.py

call deactivate