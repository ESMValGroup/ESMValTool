#! /usr/local/sci/bin/python

# test_rms
# Simple code to test the rms class

import pdb

import rms

rms_list = rms.start('antie', 'aohyf')

rms.add_example(rms_list, 'a', 'First field')
rms.add_example(rms_list, 'b', 'First field')
rms.add_example(rms_list, 'a', 'Second field')
rms.add_example(rms_list, 'b', 'Second field')

control_values = rms_list(region='global')(letter='a')
first_values = rms_list(region='global')(description='First field')
second_exper = rms_list(region='africa_land')(description='Second field', letter='b')

print rms_list

pdb.set_trace()
