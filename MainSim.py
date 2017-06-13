## MainSim.py
## Author: Aparna Rajpurkar

import Input
import Chromatin

def main():
    inputs = Input.get_input()

    for i in range(100):
        print("On sim", i)
        chromatin = Chromatin.Chromatin(inputs)
        chromatin.timesim(inputs['n'], i)

main()
