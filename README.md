# pyIRM
![alt text](https://github.com/botaoxiongyong/pyIRM/blob/master/example/Screenshot%20from%202018-09-04%2013-38-17.png)
this is for rock magnetic irm acquisition curves decompose

it was wrote in python3,

with interface wrote in pyQt5

New features:
Now you can manually adjust all the parameters and see the results immediately, and afterwards, you could try 'refit' button to better refine your fittings.



//For mac users

//an APP can run on Mac OS (not up to date)

//        https://drive.google.com/open?id=0B9ERsMj7djRxQWhLZlktRlJCQ0E

# Installation
1. Clone this repository with Git:

        $ git clone https://github.com/botaoxiongyong/pyIRM
2. install dependences:

        $ cd pyIRM
        $ pip3 install -r requirements.txt
3. install pyIRM

        # in the same derectory
        $ python3 setup.py install

# Runing

        $ python3 -m pyIRM

ps:

if you save the fig and data, they will be save in the path where you open your measured IRM data.

the initial value play a key role in fitting, so if the fit result is not promising, just try fit more times or change the component numbers.
