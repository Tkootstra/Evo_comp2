#!/usr/bin/env pypy
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 09:17:12 2020

@author: timo
"""
from Net import Net

class Graph(object):
    
    def __init__(self, node_list):
        self.node_list = node_list
        self.net_list = set()
        
    def create_nets(self):
        nets = []
    
        for n in self.node_list:
            for p in n.connection_locations:
                new_net = Net(n, self.node_list[p])
                
                n.localnets.append(new_net)
                nets.append(new_net)
        
        self.net_list = set(nets)
    
    
    def define_partitions(self, solution):
        for i in range(len(self.node_list)):
            self.node_list[i].belongs_to_partition = solution[i]
            
        self.create_nets()
    
    def count_connections(self):
        total = 0
        for node in self.node_list:
            if node.belongs_to_partition == 0:
                for loc in node.connection_locations:
                    if self.node_list[loc].belongs_to_partition == 1:
                        total += 1
        return total
    
    def compute_gains(self, node_index):
        checknode = self.node_list[node_index]
        t = checknode.belongs_to_partition
        f = not t
        
        
        for net in checknode.localnets:
            tsize = 0
            fsize = 0
            
            for n in net.nodes:
                if not n.belongs_to_partition == t:
                    tsize += 1
                if n.belongs_to_partition == f:
                    fsize += 1
            
            if tsize == 0:
                for n in net.nodes:
                    if not n.is_fixed:
                        n.gain += 1
            elif tsize == 1:
                for n in net.nodes:
                    if not n.is_fixed and not n.belongs_to_partition == t:
                        n.gain -= 1
                    
            if fsize == 0:
                for n in net.nodes:
                    if not n.is_fixed:
                        n.gain -= 1
            elif fsize == 1:
                for n in net.nodes:
                    if not n.is_fixed and n.belongs_to_partition == f:
                        n.gain += 1
                    
    
    def compute_initial_gains(self):
        self.reset_gains()
        all_gains = []
        
        for net in self.net_list:
            if net.is_cut():
                net.nodes[0].gain += 1
                net.nodes[1].gain += 1
            else:
                net.nodes[0].gain -= 1
                net.nodes[1].gain -= 1
                
            all_gains.append(net.nodes[0].gain)
            all_gains.append(net.nodes[1].gain)
            
        print(f'Max init: {max(all_gains)}')
        
    def reset_gains(self):
        new_nodes = []
        
        for node in self.node_list:
            node.gain = 0
            node.is_fixed = False
            new_nodes.append(node)
            
        self.node_list = new_nodes


#        
#        # Make sure we don't count both ways
#        to_check = [n for n in self.node_list if n.belongs_to_partition == 0]
#        
#        for node in to_check:
#            p = node.belongs_to_partition
#            
#            for idx in node.connection_locations:
#                pin = self.node_list[idx]
#                if p != pin.belongs_to_partition:
#                    node.gain += 1
#                    pin.gain += 1 
#                else:
#                    node.gain -= 1
#                    pin.gain -= 1
        
    
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
            
            
            
    