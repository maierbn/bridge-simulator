#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import sys
import copy

def draw(node_positions_, elements, fixed_dofs, loads, id_matrix, u, color, linewidth, zorder, with_loads, label):
  
  node_positions = copy.deepcopy(node_positions_)
  
  # draw elements
  for (element_no, element) in enumerate(elements):
    displacement = 4*[0.0]
    for i in range(4):
      displacement[i] = u[id_matrix[i,element_no]]
      
    # element is between nodes node0 and node1
    node0_pos = node_positions[element[0]]
    node1_pos = node_positions[element[1]]

    node0_x = node0_pos[0] + displacement[0]
    node0_y = node0_pos[1] + displacement[1]
    node1_x = node1_pos[0] + displacement[2]
    node1_y = node1_pos[1] + displacement[3]
    
    if element_no == 0:   # plot with label for the first element
      plt.plot([node0_x, node1_x], [node0_y, node1_y], "-", lw=linewidth, color=color, zorder=zorder, label=label)
    else:
      plt.plot([node0_x, node1_x], [node0_y, node1_y], "-", lw=linewidth, color=color, zorder=zorder)

  # draw nodes
  for (node_index, node_position) in enumerate(node_positions):
    dof0 = node_index*2 + 0
    dof1 = node_index*2 + 1
    node_position[0] += u[dof0]
    node_position[1] += u[dof1]
    
    #plt.plot(node_position[0], node_position[1], 'o', color=color)
    circle = plt.Circle([node_position[0], node_position[1]], 0.05, lw=linewidth, color=color, fill=True, fc='w', zorder=zorder)
    ax = plt.gca()
    ax.add_artist(circle)
    
  # draw loads
  if with_loads:
    for load in loads:
      dof = load[0]
      value = load[1]
      node_no = int(dof/2)
      node_position = node_positions[node_no]
      direction = np.array([1,0]) if dof%2==0 else np.array([0,1])
      
      ax.arrow(node_position[0], node_position[1], direction[0]*value, direction[1]*value, head_width=0.05, head_length=0.1, fc='r', ec='r')
    
  # draw fixed dofs
  for fixed_dof in fixed_dofs:
    node_no = int(fixed_dof/2)
    node_position = node_positions[node_no]
    
    size = 0.1
    offset = 0.05
    points = []
    if fixed_dof%2 == 0:    # horizontal
      points = np.array([[node_position[0]+offset, node_position[1]],\
                        [node_position[0]+offset+1.5*size, node_position[1]+size],\
                        [node_position[0]+offset+1.5*size, node_position[1]-size]])
    else:     # vertical
      points = np.array([[node_position[0], node_position[1]-offset],\
                        [node_position[0]+size, node_position[1]-offset-1.5*size],\
                        [node_position[0]-size, node_position[1]-offset-1.5*size]])
    
    polygon = plt.Polygon(points, color=color, fill=False, lw=4, joinstyle='round', zorder=zorder)
    ax.add_artist(polygon)



def output_to_file(scenario, number_nonfixed_dofs, number_elements, number_nodes, s=""):
  ax = plt.gca()
  plt.title("Truss with "+str(number_nonfixed_dofs)+" nonfixed degrees of freedom")
  plt.axis('equal')
  plt.xlim([ax.get_xlim()[0]-2, ax.get_xlim()[1]+2])    # set x limits to have a margin of 2
  plt.ylim([ax.get_ylim()[0]-2, ax.get_ylim()[1]+2])    # set y limits to have a margin of 2

  plt.legend(loc='best')
  filename = 'scenario_'+str(scenario)+s+'_dof'+str(number_nonfixed_dofs)+"_el"+str(number_elements)+'_n'+str(number_nodes)+'.png'
  plt.savefig(filename)
  print "Created file \""+filename+"\"."
