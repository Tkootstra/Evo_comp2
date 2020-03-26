#!/usr/bin/env pypy
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 09:17:12 2020

@author: timo
"""


class Graph(object):
    
    def __init__(self, node_list):
        self.node_list = node_list
    
    def define_partitions(self, solution):
        for i in range(len(self.node_list)):
            self.node_list[i].belongs_to_partition = solution[i]
    
    def count_connections(self):
        total = 0
        for node in self.node_list:
            if node.belongs_to_partition == 0:
                for loc in node.connection_locations:
                    if self.node_list[loc].belongs_to_partition == 1:
                        total += 1
        return total
    
    def compute_gains(self, node_index):
        tsize = 0
        fsize = 0
        
        checknode = self.node_list[node_index]
        t = checknode.belongs_to_partition
        f = not t
        
#        tsize += 1 # (because checknode is on t)
        
        for net in checknode.connection_locations:
            pin = self.node_list[net]
            fixed = pin.is_fixed
            
            if pin.belongs_to_partition == t:
                tsize += 1
            else:
                fsize += 1
        
            if tsize == 0:
                if not fixed:
                    pin.gain += 1
                    checknode.gain += 1
            elif tsize == 1:
                if not fixed:
                    pin.gain -= 1
                    checknode.gain -= 1
                    
            if fsize == 0:
                if not fixed:
                    pin.gain -= 1
                    checknode.gain -= 1  
            elif fsize == 1:
                if not fixed and pin.belongs_to_partition == f:
                    pin.gain += 1
                    checknode.gain += 1
                    
    
    def compute_initial_gains(self):
        for node in self.node_list:
            p = node.belongs_to_partition
            
            for idx in node.connection_locations:
                pin = self.node_list[idx]
                if p != pin.belongs_to_partition:
                    node.gain += 1
                    pin.gain += 1
                else:
                    node.gain -= 1
                    pin.gain -= 1
        
    
#    def count_single_cell_connections(self, node_index, total):
##        total = total_value
#        checknode = self.node_list[node_index]
#        
#        for loc in checknode.connection_locations:
#            
#            if not self.node_list[loc].is_fixed:
#                if self.node_list[loc].belongs_to_partition != checknode.belongs_to_partition:
#                    total += 1
#                else:
#                    total -= 1
#                        
#        return total
#                    
            
            
            
    