#!/usr/bin/env python

import sys

# input comes from STDIN (standard input)
for line in sys.stdin:
    line = line.strip()
    pairs = line.split('\n')
    # extract x and y coordinates
    for pair in pairs:
        pair = pair.strip()
        pair = pair.split()
        x = float(pair[0])
        y = float(pair[1])
        
        # find the lower bounds of the bin
        if round(x, 1) < x:
            x_low = round(x, 1)
        elif round(x, 1) > x:
            x_low = round(x, 1) - 0.1
        else:
            x_low = x - 0.1
            
        if round(y, 1) < y:
            y_low = round(y, 1)
        elif round(y, 1) > y:
            y_low = round(y, 1) - 0.1
        else:
            y_low = y - 0.1
        
        # combine the lower bounds of both x and y
        xy_low = str(x_low) + str('_') + str(y_low)
            
        print '%s\t%s' % (xy_low, 1)

