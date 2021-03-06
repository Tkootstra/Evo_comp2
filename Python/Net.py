#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 17:08:11 2020
Author: Alex Hoogerbrugge (@higher-bridge)
"""

class Net(object):
    def __init__(self, pin1, pin2):
        self.nodes = (pin1, pin2)
        self.indices = tuple(sorted([pin1.index, pin2.index]))

#        self.is_cut = pin1.belongs_to_partition != pin2.belongs_to_partition
        
    def is_cut(self):               # If 0 and 1 are in different partitions
        return self.nodes[0].belongs_to_partition != self.nodes[1].belongs_to_partition
        
    def __eq__(self, other):        # Allows it to be a set
        return (self.indices[0] == other.indices[0] and self.indices[1] == other.indices[1]) or \
               (self.indices[1] == other.indices[0] and self.indices[0] == other.indices[1])

    def __hash__(self):             # Allows it to be a set
        return hash(self.indices)
    