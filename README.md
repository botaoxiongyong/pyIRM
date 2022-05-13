# pyIRM
![alt text](https://github.com/botaoxiongyong/pyIRM/blob/master/example/py_irm_gui.png)
this is for rock magnetic irm acquisition curves decompose


with interface wrote in pyQt5

New features:

1. You can manually adjust all the parameters and see the results immediately, and afterwards, you could try 'refit' button to better refine your fittings.

2. All components have Gauss uncertainties.

# License
MIT License

# Data Format
AGM2900/3900

VSM8600

txt file with two columns of Field and Moment, with delimeter of comma

# Installation
Python3 version: 3.9
        $ pip install git+https://github.com/botaoxiongyong/pyIRM.git

# Runing
        $ python3 -m pyIRM

# Jupyter notebook
        $ check the example: .../pyIRM/pyIRM.ipynb

ps:

if you save the fig and data, they will be save in the path where you open your measured IRM data.

the initial value play a key role in fitting, so if the fit result is not promising, just try fit more times or change the component numbers.
