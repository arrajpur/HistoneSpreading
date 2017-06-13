## MainAnim.py
## Author: Aparna Rajpurkar
import Input
import Chromatin
from MyEnum import Divisions
from Animate import animate_from_file

def main():
    inputs = Input.get_input()

    # initialize Chromatin object
    chromatin = Chromatin.Chromatin(inputs)
    # run simulation
    chromatin.timesim(inputs['n'], 0)
    
    # get number of divisions
    div = 0
    if inputs['d'] != Divisions.NONE:
        div = inputs['data']['divisions']

    # animate plot using output file from simulation
    animate_from_file(inputs['o'] + '_0.txt', inputs['n'], inputs['f'], inputs['t'], inputs['o'], div)

# run main
main()
