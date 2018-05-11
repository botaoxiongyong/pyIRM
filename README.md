# pyIRM
this is for rock magnetic irm acquisition curves decompose

it was wrote in python3, 

with interface wrote in pyQt5

New features:
Now you can manually adjust all the parameters and see the results immediately, and afterwards, you could try 'refit' button to better refine your fittings.



For mac users

an APP can runing on Mac OS (not up to date)

https://drive.google.com/open?id=0B9ERsMj7djRxQWhLZlktRlJCQ0E

the package list:

        matplotlib
        pandas
        PyQt5
        lea
        scipy
        numpy
        sklearn
        lmfit

the most easy way to insall all of this above is using pip

on Apple Mac

     brew install python3
     brew install pip
     pip3 install matplotlib, PyQt5........
     
on Linux 

     sudo apt-get install python3
     sudo apt-get install pip
     pip3 install matplotlib, PyQt5........
    
open the terminal

     cd /home/Documents/ (if you put the pyIRM.py file in /home/Documents/)
     python3 pyIRM.py
     
ps: 

if you save the fig and data, they will be save in the path where you open your measured IRM data.

the initial value play a key role in fitting, so if the fit result is not promising, just try fit more times or change the component numbers.
