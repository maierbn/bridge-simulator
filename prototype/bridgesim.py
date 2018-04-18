#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import sys
import draw_truss
import solve_truss
from matplotlib import animation

ea = 10         # material parameter (E modulus * cross area)
scenario = 11    # select predefined scenario
rho = 1.0

if len(sys.argv) > 1:
  scenario = int(sys.argv[1])
  print "Scenario",scenario
  print "=========="

if scenario == 0:
  #------ 3-Gelenkträger------
  node_positions = [[0,0], [5,0], [2,2]]    # for each node [posx, posy]
  elements = [[0, 2], [1, 2]]               # the node numbers [node0, node1]
  loads = [[4, 1], [5, 0.5]]                # [dofno, value]
  fixed_dofs = [0, 1, 2, 3]

elif scenario == 1:
  #----- 1 Stab ------
  node_positions = [[0,0], [1,0]]    # for each node [posx, posy]
  elements = [[0, 1]]               # the node numbers [node0, node1]
  loads = [[2, 1]]                # [dofno, value]
  fixed_dofs = [0, 1, 3]

elif scenario == 2:
  #------ more complex -----
  node_positions = [[0,0], [2,5], [3,2], [4,7], [5,2]]    # for each node [posx, posy]
  elements = [[0,1], [0,2], [1,2], [1,3], [2,3], [2,4], [3,4]]               # the node numbers [node0, node1]
  loads = [[4,0.2], [5,-0.5]]                # [dofno, value]
  fixed_dofs = [0, 1, 9]

elif scenario == 3:
  #------ 3-Gelenkträger einfach statisch überbestimmt------
  node_positions = [[0,0], [5,0], [2,2]]    # for each node [posx, posy]
  elements = [[0, 2], [1, 2], [0, 1]]        # the node numbers [node0, node1]
  loads = [[4, 1], [5, 0.5]]                # [dofno, value]
  fixed_dofs = [0, 1, 2, 3]
  
elif scenario == 4:
  #------ aufgesetzter 3-Gelenkträger einfach statisch überbestimmt------
  node_positions = [[0,0], [2,0], [1,1], [3,0], [5,0], [4,2], [2,4]]    # for each node [posx, posy]
  elements = [[0, 2], [1, 2], [3, 5], [4,5], [2,5], [2,6], [5,6]]        # the node numbers [node0, node1]
  loads = [[4,-0.5], [12, 1], [13, 0.5]]                # [dofno, value]
  fixed_dofs = [0, 1, 2, 3, 6, 7, 8, 9]

elif scenario == 5:
  #------ same as scenario 4, different order ------
  node_positions = [[0,0], [2,0], [1,1], [3,0], [5,0], [4,2], [2,4]]    # for each node [posx, posy]
  elements = [[2,5], [0, 2], [1, 2], [3, 5], [4,5], [2,6], [5,6]]        # the node numbers [node0, node1]
  loads = [[4,-0.5], [12, 1], [13, 0.5]]                # [dofno, value]
  fixed_dofs = [0, 1, 2, 3, 6, 7, 8, 9]

elif scenario == 6:
  #------ same as scenario 4, different order ------
  node_positions = [[0,0], [2,0], [1,1], [3,0], [5,0], [4,2], [2,4]]    # for each node [posx, posy]
  elements = [[2,5], [0, 2], [1, 2], [3, 5], [4,5], [2,6], [5,6], [1,6]]        # the node numbers [node0, node1]
  loads = [[4,-0.5], [12, 1], [13, 0.5]]                # [dofno, value]
  fixed_dofs = [0, 1, 2, 3, 6, 7, 8, 9]

elif scenario == 7:
  #------ static example ------
  node_positions = [[0,0], [2,0], [1,1]]    # for each node [posx, posy]
  elements = [[0,2], [1, 2]]        # the node numbers [node0, node1]
  loads = [[4,0.5], [5, 1]]                # [dofno, value]
  fixed_dofs = [0, 1, 2, 3]

elif scenario == 8:
  #------ dynamic example, free member diagonal up right ------
  node_positions = [[0,0], [2,0], [1,1], [3,2]]    # for each node [posx, posy]
  elements = [[0,2], [1, 2], [2, 3]]        # the node numbers [node0, node1]
#  loads = [[6,0.5], [7, 1]]                # [dofno, value]
  loads = [[6,0.5]]                # [dofno, value]
  fixed_dofs = [0, 1, 2, 3]

elif scenario == 9:
  #------ dynamic example, free member horizontal to the right ------
  node_positions = [[0,0], [2,0], [1,1], [3,1]]    # for each node [posx, posy]
  elements = [[0,2], [1, 2], [2, 3]]        # the node numbers [node0, node1]
#  loads = [[6,0.5], [7, 1]]                # [dofno, value]
  loads = [[6,0.5]]                # [dofno, value]
  fixed_dofs = [0, 1, 2, 3]

elif scenario == 10:
  #------ dynamic example, single member ------
  node_positions = [[0.,0.], [1.,-0.5]]    # for each node [posx, posy]
  elements = [[0,1]]        # the node numbers [node0, node1]
  loads = [[2,1.0]]                # [dofno, value]
  fixed_dofs = [0, 1]

elif scenario == 11:
  #------ dynamic example, more complex ------
  node_positions = [[0,0], [2,0], [1,1], [3,0], [5,0], [4,2], [2,4]]    # for each node [posx, posy]
  elements = [[2,5], [0, 2], [1, 2], [3, 5], [4,5], [2,6], [5,6], [1,6]]        # the node numbers [node0, node1]
  loads = [[4,-0.5], [12, 1], [13, 0.5], [7,-0.2], [9,-0.3]]                # [dofno, value]
  fixed_dofs = [0, 1, 3]

elif scenario == 12:
  #------ dynamic example, single horizontal member, 1 DOF ------
  node_positions = [[0.,0.], [1.,0.]]    # for each node [posx, posy]
  elements = [[0,1]]        # the node numbers [node0, node1]
  loads = [[2,1.0]]                # [dofno, value]
  fixed_dofs = [0, 1,3]
  
number_nodes = len(node_positions)
number_elements = len(elements)
number_dofs = number_nodes * 2    # every nodes has two displacements (degrees of freedom = dof)
number_nonfixed_dofs = number_dofs - len(fixed_dofs)

def animate(frame_no):
  global node_positions, elements, loads, fixed_dofs, ea, rho, scenario, ax1, t_start, initial_velocity_reduced
  
  global_time = frame_no/10.
  current_time = global_time - t_start
  
  if frame_no == 0:
    return
  
  #print "frame_no=",frame_no,", t_start=",t_start,", current_time=",current_time
  
  options = {
    "ea": ea, "rho": rho, "scenario": scenario, 
    "n_ode_iterations": 500, 
    "damping_factor": 0.2,
    "t_end": current_time, 
    "initial_velocity_reduced": initial_velocity_reduced
  }
  result = solve_truss.solve_truss_dynamic(node_positions, elements, loads, fixed_dofs, options)
  u = result["u"]
  velocity_reduced = result["velocity_reduced"]
  id_matrix = result["id_matrix"]
  u_reduced = result["u_reduced"]
  ku = result["ku"]

  ax1.clear()
  plt.sca(ax1)

  # draw truss in reference configuration                                 u,               color, linewidth, zorder, with_loads
  draw_truss.draw(node_positions, elements, fixed_dofs, loads, id_matrix, np.zeros([number_dofs]), 'k', 3, 0, False, "reference configuration")

  # draw deformed configuration
  draw_truss.draw(node_positions, elements, fixed_dofs, loads, id_matrix, u, 'g', 3, 1, True, "deformed configuration")

  if frame_no%1 == 0: 
    print "u:",u,", velocity: ", velocity_reduced
    t_start = global_time
    initial_velocity_reduced = velocity_reduced
    
    # update node positions
    for i in range(number_nodes):
      dof0 = 2*i
      dof1 = 2*i + 1
      node_positions[i][0] += u[dof0]
      node_positions[i][1] += u[dof1]

    print "create new reference configuration at time ",t_start

if scenario < 7:    # static examples
  force_method = False
  (u, id_matrix) = solve_truss.solve_truss(node_positions, elements, loads, fixed_dofs, ea, scenario, force_method)
  
  print "----- plot -----"
  fig = plt.figure("main", figsize=(14,10))

  # draw truss in reference configuration                                 u,               color, linewidth, zorder, with_loads
  draw_truss.draw(node_positions, elements, fixed_dofs, loads, id_matrix, np.zeros([number_dofs]), 'k', 3, 0, False, "reference configuration")

  # draw deformed configuration
  draw_truss.draw(node_positions, elements, fixed_dofs, loads, id_matrix, u, 'g', 3, 1, True, "deformed configuration")

  # create file
  draw_truss.output_to_file(scenario, number_nonfixed_dofs, number_elements, number_nodes)
  #plt.show()

else:   # dynamic scenarios
  fig = plt.figure("main", figsize=(14,10))
  ax1 = fig.gca()
  #ax2 = fig.add_subplot(122)
  #ax1 = fig.add_subplot(121)
  
  initial_velocity_reduced = None
  t_start = 0

  ax1.set_aspect('equal', 'datalim')
  ax1.set_xlim([ax1.get_xlim()[0]-2, ax1.get_xlim()[1]+2])    # set x limits to have a margin of 2
  ax1.set_ylim([ax1.get_ylim()[0]-2, ax1.get_ylim()[1]+2])    # set y limits to have a margin of 2

  #ax2.set_xlabel('u')
  #ax2.set_ylabel('F')

  anim = animation.FuncAnimation(fig, animate, frames=10000, interval=20)
  
  # save to animation file
  #anim.save(animation_name)
  #print "saved to {}".format(animation_name)
  
  plt.show()
