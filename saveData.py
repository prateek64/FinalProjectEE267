

## establish connection to the serial port that your arduino 
## is connected to.

# Python script that reads the serial port of your teensy and saves the serial port outputted data into a text file.
# For my experiments the only data that was written to the serial port was the photodiode location data
# 

import serial 

connected = False
locations=['/dev/ttyUSB0','/dev/ttyUSB1','/dev/ttyUSB2','/dev/ttyUSB3' , '/dev/tty.usbmodem2815011']

for device in locations:
    try:
        print "Trying...",device
        ser = serial.Serial(device, 115200)
        break
    except:
        print "Failed to connect on",device

## loop until the arduino tells us it is ready
while not connected:
    serin = ser.read()
    connected = True

## open text file to store the current 
##gps co-ordinates received from the rover    
text_file = open("positionPhotoDiodes_4thrun.txt", 'w')
## read serial data from arduino and 
## write it to the text file 'position.txt'
while 1:
    if ser.inWaiting():
        x=ser.read()
        print(x) 
        text_file.write(x)
        # if x=="\n":
        #      text_file.seek(0)
        #      text_file.truncate()
        # text_file.flush()

## close the serial connection and text file
text_file.close()
ser.close()