#!/usr/bin/env pypy
# -*- coding: utf-8 -*-

class Node(object):
    
    def __init__(self, index_number, n_connections, connection_locs):
        self.index = index_number
        self.number_of_connections = n_connections
        self.connection_locations = connection_locs
        self.is_fixed = False
        self.gain = 0;
        self.belongs_to_partition = None
    
    def flip_partition(self):
        self.belongs_to_partition = 1 if self.belongs_to_partition == 0 else 0
    
    def __eq__(self, other):
        return self.index == other.index
        
    
   
    