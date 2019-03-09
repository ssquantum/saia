for /l %%n in (0,1,49) do (
	echo %%n
	echo %%n > "C:\Users\qgtx64\DocumentsCDrive\Python Scripts\SAIA\currentfile.txt"
	copy "C:\Users\qgtx64\DocumentsCDrive\Python Scripts\TestROIImages\TestIm%%n.asc" "C:\Users\qgtx64\DocumentsCDrive\Python Scripts\SAIA\CameraImages\CameraIm.asc"
	ping 127.0.0.1 -n 1 -w 500> nul

)