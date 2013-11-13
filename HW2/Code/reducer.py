#!/usr/bin/env python

from operator import itemgetter
import sys

current_xy_low = None
current_count = 0
xy_low = None

# input comes from STDIN
for line in sys.stdin:
    # remove leading and trailing whitespace
    line = line.strip()
    # parse the input we got from mapper.py
    xy_low, count = line.split('\t', 1)

    # convert count (currently a string) to int
    try:
        count = int(count)
    except ValueError:
        # count was not a number, so silently
        # ignore/discard this line
        continue

    # if the key is the same, sum the counts
    if current_xy_low == xy_low:
        current_count += count
    else:
        if current_xy_low:
            # write result to STDOUT in csv format
            x_low, y_low = current_xy_low.split('_', 1)
            out_str = x_low + ',' + str(float(x_low)+0.1) + ',' + y_low + ',' + str(float(y_low)+0.1)
            print '%s,%s' % (out_str, current_count)
        current_count = count
        current_xy_low = xy_low

# output the last word if needed
if current_xy_low == xy_low:
    x_low, y_low = current_xy_low.split('_', 1)
    out_str = x_low + ',' + str(float(x_low)+0.1) + ',' + y_low + ',' + str(float(y_low)+0.1)
    print '%s,%s' % (out_str, current_count)