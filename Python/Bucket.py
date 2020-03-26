#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 10:27:25 2020

@author: timo
"""


class Bucket(object):
    
    def __init__(self, dict_a, dict_b):
        self.bucketA = dict_a
        self.bucketB = dict_b
        self.bucketA_max_pointer = -999
        self.bucketB_max_pointer = -999
        self.bucketA_size = 250
        self.bucketB_size = 250
        
    def add_to_bucket(self, partition, key, node):
        if partition == 0:
            #  dit heb ik toegevoegd, zodat die 
            # automatische een nieuwe key toevoegt als die nog niet in de dict zit
            if key not in self.bucketA.keys():
                self.bucketA[key] = []
            self.bucket[key].append(node)
            if key > self.bucketA_max_pointer:
                self.bucketA_max_pointer = key
        else:
            #  same here
            if key not in self.bucketB.keys():
                self.bucketB[key] = []
                
            self.bucketB[key].append(node)
            if key > self.bucketB_max_pointer:
                self.bucketB_max_pointer = key
                
    def update_bucket(self, partition, node):
        # yankadoodle
        if partition == 0:
            for key in self.bucketA.keys():
                for item in self.bucketA[key]:
                    if item == node:
                        self.bucketA[key].remove(item)
                        self.add_to_bucket(partition, key, node)
                        
        else:
            for key in self.bucketB.keys():
                for item in self.bucketB[key]:
                    if item == node:
                        self.bucketB[key].remove(item)
                        self.add_to_bucket(partition, key, node)
    
    def pop_from_bucket_key(self, partition):
        if partition == 0:
            if len(self.bucketA[self.bucketA_max_pointer]) > 0:
                item = self.bucketA[self.bucketA_max_pointer].pop()
                self.bucketA_size -= 1
                item.is_fixed = True
                return item
            else:
                self.bucketA_max_pointer -= 1
                return self.pop_from_bucket_key(partition)
        else:
            if len(self.bucketB[self.bucketB_max_pointer]) > 0:
                item = self.bucketB[self.bucketB_max_pointer].pop()
                self.bucketB_size -= 1
                item.is_fixed = True
                return item
            else:
                self.bucketB_max_pointer -= 1
                return self.pop_from_bucket_key(partition)
            
    def gain_sum(self):
        gain = 0
        for key in self.bucketA.keys():
            gain = gain + (key * len(self.bucketA[key]))
            
        for key in self.bucketB.keys():
            gain = gain + (key * len(self.bucketB[key])) 
        
        return gain
        
            
            
            
        
        
        
        
            
                    
            