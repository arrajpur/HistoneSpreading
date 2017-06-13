## MainRunOnce.py
## Author: Aparna Rajpurkar

import Input
import Chromatin

def main():
    inputs = Input.get_input()

    chromatin = Chromatin.Chromatin(inputs)
    chromatin.timesim(inputs['n'], 0)

main()
