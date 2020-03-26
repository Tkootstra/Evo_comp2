#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
                        total+=1
        return total
    
    def count_single_cell_connections(self, node_index, total_value):
        total = total_value
        checknode = self.node_list[node_index]
        for loc in checknode.connection_locations:
            if not self.node_list[loc].is_fixed:
                if self.node_list[loc].belongs_to_partition != checknode.belongs_to_partition:
                    total+= 1
                else:
                    total+= 1
                    
        return total
                    
            
            
            
    