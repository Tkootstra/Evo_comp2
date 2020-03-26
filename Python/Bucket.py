#!/usr/bin/env pypy
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 10:27:25 2020

@author: timo
"""


class Bucket(object):
    
    def __init__(self, dict_a, dict_b, maxcon):
        self.bucketA = dict_a
        self.bucketB = dict_b
        self.bucketA_max_pointer = -19
        self.bucketB_max_pointer = -19
        self.bucketA_size = 250
        self.bucketB_size = 250
        
    def add_to_bucket(self, partition, key, node):
        if partition == 0:
            #  dit heb ik toegevoegd, zodat die 
            # automatische een nieuwe key toevoegt als die nog niet in de dict zit
#            if key not in self.bucketA.keys():
#                self.bucketA[key] = []
            
            self.bucketA[key].append(node)
            
            if key > self.bucketA_max_pointer:
                self.bucketA_max_pointer = key
        
        else:
            #  same here
#            if key not in self.bucketB.keys():
#                self.bucketB[key] = []
                
            self.bucketB[key].append(node)
            
            if key > self.bucketB_max_pointer:
                self.bucketB_max_pointer = key
                
    def update_bucket(self, partition, new_key, node):
        # yankadoodle
        if partition == 0:
            for key in self.bucketA.keys():
                for item in self.bucketA[key]:
                    if item == node:
                        self.bucketA[key].remove(item)
                        self.add_to_bucket(partition, new_key, node)
                        return
                        
        else:
            for key in self.bucketB.keys():
                for item in self.bucketB[key]:
                    if item == node:
                        self.bucketB[key].remove(item)
                        self.add_to_bucket(partition, new_key, node)
                        return
    
    def pop_from_bucket_key(self, partition):
        if partition == 0:
            if len(self.bucketA[self.bucketA_max_pointer]) > 0:
                item = self.bucketA[self.bucketA_max_pointer].pop(0)
                self.bucketA_size -= 1
                item.is_fixed = True
                return item
            else:
#                print(self.bucketA_max_pointer, ' new loop A')
                self.bucketA_max_pointer -= 1
#                print(self.bucketA_max_pointer, ' updated A')
                return self.pop_from_bucket_key(partition)
        else:
            if len(self.bucketB[self.bucketB_max_pointer]) > 0:
                item = self.bucketB[self.bucketB_max_pointer].pop(0)
                self.bucketB_size -= 1
                item.is_fixed = True
                return item
            else:
#                print(self.bucketB_max_pointer, ' new loop B')
                self.bucketB_max_pointer -= 1
#                print(self.bucketB_max_pointer, ' updated B')
                return self.pop_from_bucket_key(partition)
            
    def gain_sum(self):
        gain = 0
        for key in self.bucketA.keys():
            gain = gain + (key * len(self.bucketA[key]))
            
        for key in self.bucketB.keys():
            gain = gain + (key * len(self.bucketB[key])) 
        
        return gain
    
#    def init_gain(self, graph):
#        # This works fine
#        for i in range(len(graph.node_list)):
#            current = graph.node_list[i]
#            
#            if not current.is_fixed:
#                graph.node_list[current.index].flip_partition()
#                gain = graph.count_single_cell_connections(current.index, current.gain)
#                graph.node_list[current.index].flip_partition()
#                
#                if gain < -16 or gain > 16:
#                    print(gain)
#                
#                self.add_to_bucket(current.belongs_to_partition, gain, current)
#                graph.node_list[i].gain = gain
#        
#    def update_gain(self, graph, neighbors):
#        # Gains get out of hand, usually towards -20 or -30. -16 should be minimum
#        for idx in neighbors:
#            current = graph.node_list[idx]
#            
#            if not current.is_fixed:
#                graph.node_list[current.index].flip_partition()                         # Temp flip
#                gain = graph.compute_gain(current.index, current.gain)
#                graph.node_list[current.index].flip_partition()                         # Flip back
#                
#                if gain < -16 or gain > 16:
#                    print(f'{idx}: {current.gain} -> {gain} = {gain - current.gain}')
#                
#                self.update_bucket(current.belongs_to_partition, gain, current)
#                graph.node_list[idx].gain = gain   
                
    
    def update_gain_test(self, graph, node_index):
        graph.compute_gains(node_index)
        
        for idx in graph.node_list[node_index].connection_locations:
            current = graph.node_list[idx]
            
            if current.gain < -16 or current.gain > 16:
                print(f'{idx}: {current.gain}')
            
            self.update_bucket(current.belongs_to_partition, current.gain, current)
            
    def init_gain_test(self, graph):
        graph.compute_initial_gains()
        
        for node in graph.node_list:
            if node.gain < -16 or node.gain > 16:
                print(f'{node.index}: {node.gain} initial')
            
            self.add_to_bucket(node.belongs_to_partition, node.gain, node)
            
        
        
        
        
            
                    
            