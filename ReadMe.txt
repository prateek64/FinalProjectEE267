
Project : Estimating Lighthouse Intrinsic Parameters using Levenberg Marquardt Optimization

Team Members : Prateek Murgai 

Contact : pmurgai@ccrma.stanford.edu



Section-1 : ** Data Files :** 

1. photoDiodeData_4thrun.txt ( Contains the most data points )
2. photoDiodeData_2ndrun.txt
3. photoDiodeData.txt


Each row in the data file is 8 points [x1 , y1 , x2 , y2 , x3 , y3 , x4 , y4] where (xi , yi) are the photodiode locations of the ith photodiode

As we have only 4 photodiodes thus 8 values in a row.

You can record data in any format you want but remember to modify saveData.py in that case



Section-2 : **To directly run the system (Linux / MacOS) : **

calibrateLightHouse.cpp - Run this file to estimate the intrinsic parameters of the lighthouse , this file uses the above data files which were recorded offline.

calibrateLightHouseSimulated.cpp - Run this file to estimate intrinsic parameters of the lighthouse using simulated Lighthouse data given in simulatedLighthouseData.h


You need a C++ compiler to run the files . I used g++ to run them but you can use something else too

To build an output file run the following on terminal : g++ calibrateLightHouse.cpp -o calibrateLightHouse.out

If you are using Windows please use a IDE/Visual compiler to do this

To run the compiled file write the following on terminal : ./calibrateLightHouse.out

Once you run this command the optimization should happen automatically and once the code finishes note down the final values of the intrinsic parameters ( F , XP , YP ) . An example output on terminal looks like this : 


#################################

ThX : 12.8336
ThY : 25.184
ThZ : -1.64781
TrX : -207.066
TrY : 382.715
TrZ : -1007.45
F : 0.885537
XP : -0.17977
YP : 0.256887

Residual : 2.53754e-06
#################################

							  					


Section-3 : **If recording your own data :** 

1. First run saveData.py
2. Then run preProcessData.py on the file saved in the previous step 
3. Then use the output file from preProcessData.py to run calibrateLightHouse files

** Note that the serial port should be outputting only photodiode data when you run saveData.py else your system will save garbage data. Modify these files according to your own use, these files cannot be made generic enough. 


