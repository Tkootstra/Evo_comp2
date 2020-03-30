#!/usr/bin/env pypy
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 09:17:12 2020

@author: timo
"""
from Net import Net
from copy import deepcopy
from Node import Node
import time

class Graph(object):
    
    def __init__(self, node_list, solution):
        self.node_list = node_list
        self.net_list = set()
        self.free_nodes = len(node_list)
        
        self.define_partitions(solution)
        
    # def __deepcopy__(self, memo):
    #     cls = self.__class__
    #     result = cls.__new__(cls)
    #     memo[id(self)] = result
    #     for k, v in self.__dict__.items():
    #         setattr(result, k, deepcopy(v, memo))
    #     return result
        
    def create_nets(self):
        # Creates nets and puts them in the node
        nets = []
    
        for n in self.node_list:
            for p in n.connection_locations:
                new_net = Net(n, self.node_list[p])
        
                n.localnets.append(new_net)
                nets.append(new_net)
                        
        self.net_list = set(nets)
    
    
    def define_partitions(self, solution):
        self.reset_gains() 
        
        for i in range(len(self.node_list)):
            self.node_list[i].belongs_to_partition = solution[i]
        
        self.create_nets()
    
    def count_connections(self):
        self.create_nets()
        
        total = 0
        for net in self.net_list:
            if net.is_cut():
                total += 1
                
        return total
        
#        total = 0
#        for node in self.node_list:
#            if node.belongs_to_partition == 0:
#                for loc in node.connection_locations:
#                    if self.node_list[loc].belongs_to_partition == 1:
#                        total += 1
#        return total
    
    def compute_gains(self, checknode):
#        checknode = self.node_list[node_index]
        t = checknode.belongs_to_partition
        f = 0 if t == 1 else 1
        
        
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
        all_gains = []      # For testing purposes
        
        for net in self.net_list:
            if net.is_cut():
                net.nodes[0].gain += 1
                net.nodes[1].gain += 1
            else:
                net.nodes[0].gain -= 1
                net.nodes[1].gain -= 1
                
            all_gains.append(net.nodes[0].gain)
            all_gains.append(net.nodes[1].gain)
            
#        print(f'Max init: {max(all_gains)}')
        
    def reset_gains(self):
        new_nodes = []
        
        for node in self.node_list:
            node.gain = 0
            node.is_fixed = False
            new_nodes.append(node)
            
        self.node_list = new_nodes
        self.free_nodes = 250
    
    
    def get_best_node(self, partition):
        bucket = [n for n in self.node_list if n.belongs_to_partition == partition]
            
        gains = [(n.gain, n) for n in bucket if not n.is_fixed]
        best_node = max(gains, key=lambda item:item[0])[1]
        
        self.free_nodes -= 1
        best_node.is_fixed = True
        self.node_list[best_node.index].is_fixed = True
        
        return best_node

    def calc_gain_sum(self):
        return sum([n.gain for n in self.node_list if not n.is_fixed])
    
    def flip_partition(self, node_index):
        self.node_list[node_index].flip_partition()
        self.node_list[node_index].is_fixed = True
#        self.create_nets()
        

    

            
            
    