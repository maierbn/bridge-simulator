#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import sys
import draw_truss
import scipy
import scipy.integrate

figure_cntr = 0
def solve_truss(node_positions, elements, loads, fixed_dofs, ea, scenario, force_method):
  """
    solve the truss, if a solution exists, returns vector of nodal displacements (u,id_matrix)
    if no solution is found, returns (None,id_matrix)
    if force_method is True, then the force method is used (mainly useful for hand calculations)
    if force_method is False, the direct stiffness method = displacement method is used
  """
  global figure_cntr

  print ""
  print "======= function solve truss with {} elements =======".format(len(elements))

  number_nodes = len(node_positions)
  number_elements = len(elements)
  number_dofs = number_nodes * 2    # every nodes has two displacements (degrees of freedom = dof)
  number_nonfixed_dofs = number_dofs - len(fixed_dofs)

  print "---------- Input parameters ------------"
  print len(node_positions), "node_positions [x, y]:          ",node_positions
  print len(elements),       "elements [node start, node end]:",elements
  print len(loads),          "loads [dof, value]:             ",loads
  print len(fixed_dofs),     "fixed_dofs:                     ",fixed_dofs
  print "parameter EA =",ea

  # create ID-matrix
  id_matrix = np.zeros([4, number_elements], int)
  for (index, element) in enumerate(elements):
    id_matrix[0, index] = element[0]*2 + 0
    id_matrix[1, index] = element[0]*2 + 1
    id_matrix[2, index] = element[1]*2 + 0
    id_matrix[3, index] = element[1]*2 + 1

  print ""
  print "id_matrix:"
  print id_matrix
  
  # evaluate counting criterion
  value = len(fixed_dofs) + number_elements - 2*number_nodes
  print "counting criterion: ",len(fixed_dofs),"+",number_elements,"- 2 *",number_nodes,"=",value
  if value == 0:
    print " -> statically determined (statisch bestimmt)"
  elif value > 0 and force_method:
    print " -> hyperstatic (statisch überbestimmt)"
    
    # determine element to be replaced by X force
    for element_no in range(number_elements):
      subsystem_node_positions = list(node_positions)
      subsystem_elements = list(elements)
      subsystem_loads = list(loads)
      subsystem_fixed_dofs = list(fixed_dofs)
    
      # remove selected element
      subsystem_elements.pop(element_no)
      
      print "Try if subsystem without element {} is solvable.".format(element_no)
      
      # test if subsystem is solvable
      (subsystem_u, _) = solve_truss(subsystem_node_positions, subsystem_elements, subsystem_loads, subsystem_fixed_dofs, ea, scenario, force_method)
    
      if subsystem_u is None:
        print ""
        print "Subsystem is not solvable. Remove a different element instead."
        continue
        
      # subsystem had a solution
      # store solution of 0-system
      subsystem_0_u = subsystem_u.copy()
      
      # compute solution of X-system
      subsystem_x_loads = []        # no external loads
      # add X=1 force inside removed element
      node_start = elements[element_no][0]
      node_end = elements[element_no][1]
      initial_position_node_start = np.array(node_positions[elements[element_no][0]])
      initial_position_node_end = np.array(node_positions[elements[element_no][1]])
      angle = np.arctan2(initial_position_node_end[1]-initial_position_node_start[1], initial_position_node_end[0]-initial_position_node_start[0])
      
      subsystem_x_loads = [
        [2*node_start+0, np.cos(angle)],
        [2*node_start+1, np.sin(angle)],
        [2*node_end+0, -np.cos(angle)],
        [2*node_end+1, -np.sin(angle)]
      ]
      
      print "Subsystem without element {} was solved and is now \
        considered as O-system. Also compute X-system without that element.".format(element_no)
      print ""
      print "----- setup loads for X-system (force X=1 in x-element) -----"
      print "x-element is element {}, between nodes {} and {}".format(element_no, node_start, node_end)
      print "initial positions of nodes: ", initial_position_node_start, initial_position_node_end
      print "angle of element:", angle
      print "subsystem_x_loads =", subsystem_x_loads
      print "now solve X-system"
      
      (subsystem_x_u, subsystem_id_matrix) = solve_truss(subsystem_node_positions, subsystem_elements, subsystem_x_loads, subsystem_fixed_dofs, ea, scenario, force_method)
    
      # compute lengthening of removed element
      # length difference in 0-system
      displacement_x_node_start = subsystem_0_u[2*elements[element_no][0]+0] 
      displacement_y_node_start = subsystem_0_u[2*elements[element_no][0]+1] 
      displacement_x_node_end = subsystem_0_u[2*elements[element_no][1]+0] 
      displacement_y_node_end = subsystem_0_u[2*elements[element_no][1]+1] 
      actual_position_node_start = initial_position_node_start + np.array([displacement_x_node_start, displacement_y_node_start])
      actual_position_node_end = initial_position_node_end + np.array([displacement_x_node_end, displacement_y_node_end])
      
      initial_length = np.linalg.norm(initial_position_node_end - initial_position_node_start)
      actual_length = np.linalg.norm(actual_position_node_end - actual_position_node_start)
      d0 = actual_length - initial_length
      
      print ""
      print "----- compute X -----"
      print "0-system (system with all external loads and element {} removed):".format(element_no)
      print "   displacements:"
      print "     node {}, x: {}, y: {}".format(node_start, displacement_x_node_start, displacement_y_node_start)
      print "     node {}, x: {}, y: {}".format(node_end, displacement_x_node_end, displacement_y_node_end)
      print "   initial distance between nodes: initial_length =", initial_length
      print "   actual distance between nodes: actual_length =", actual_length
      print "   extension d0 = actual_length-initial_length = {}".format(d0)
      
      
      # length difference in X-system
      displacement_x_node_start = subsystem_x_u[2*elements[element_no][0]+0] 
      displacement_y_node_start = subsystem_x_u[2*elements[element_no][0]+1] 
      displacement_x_node_end = subsystem_x_u[2*elements[element_no][1]+0] 
      displacement_y_node_end = subsystem_x_u[2*elements[element_no][1]+1] 
      actual_position_node_start = initial_position_node_start + np.array([displacement_x_node_start, displacement_y_node_start])
      actual_position_node_end = initial_position_node_end + np.array([displacement_x_node_end, displacement_y_node_end])
      actual_length = np.linalg.norm(actual_position_node_end - actual_position_node_start)
      
      # sigma = E*epsilon <=>  F/A = E*u/l <=> EA/l*u = F <=> u = F*l/EA
      
      extension_x_element =  1*initial_length / ea
      d1 = actual_length - initial_length - extension_x_element
      
      print "X-system (system without external loads and with element {} cut and member force X=1):".format(element_no)
      print "   displacements:"
      print "     node {}, x: {}, y: {}".format(node_start, displacement_x_node_start, displacement_y_node_start)
      print "     node {}, x: {}, y: {}".format(node_end, displacement_x_node_end, displacement_y_node_end)
      print "   initial distance between nodes: initial_length =", initial_length
      print "   actual distance between nodes: actual_length =", actual_length
      print "   additional extension of element {}: extension_x_element = {}".format(element_no, extension_x_element)
      print "   extension d1 = actual_length-initial_length-extension_x_element = {}".format(d1)
      
      # d0 = -X*d1  =>  X = -d0/d1
      x_factor = -d0/d1
      
      # superpose displacements of 0-system and X-system
      u = subsystem_0_u + x_factor * subsystem_x_u

      # scale loads of X-system (for plot)
      for i in range(len(subsystem_x_loads)):
        subsystem_x_loads[i][1] *= x_factor
      
      print "X = -d0/d1 = {}".format(-d0/d1)
      print ""
      print "----- superpose displacements -----"
      print "u = subsystem_0_u + X * subsystem_x_u:"
      print u
            
      # output file with truss in subsystems                                u,               color, linewidth, zorder, with_loads
      plt.figure("subsystem-"+str(figure_cntr), figsize=(14,10))
      draw_truss.draw(subsystem_node_positions, subsystem_elements, subsystem_fixed_dofs, subsystem_loads, subsystem_id_matrix, np.zeros([number_dofs]), 'k', 3, 0, False, "reference configuration")
      draw_truss.draw(subsystem_node_positions, subsystem_elements, subsystem_fixed_dofs, subsystem_loads, subsystem_id_matrix, subsystem_0_u, 'b', 3, 1, True, "0-system deformed")
      draw_truss.draw(subsystem_node_positions, subsystem_elements, subsystem_fixed_dofs, subsystem_x_loads, subsystem_id_matrix, x_factor*subsystem_x_u, 'r', 3, 1, True, "X-system deformed (with actual X="+str(x_factor)+")")
      draw_truss.output_to_file(scenario, len(subsystem_fixed_dofs), len(subsystem_elements), number_nodes, str(unichr(97+figure_cntr)))
      figure_cntr += 1
      
      return (u, id_matrix)
    
  else:
    print " -> hypostatic (statisch unterbestimmt, verschieblich)"

  u = np.zeros([number_dofs])   # vector of displacements, this is seeked, it is the solution of k_global * u = f
  f = np.zeros([number_dofs])   # vector of loads, this is the right hand side
  k_global = np.zeros([number_dofs, number_dofs])   # the global stiffness matrix

  print ""
  print "----- assemble global stiffness matrix from local matrices -----"

  # local stiffness matrix
  k_local = np.array([[ 1, 0, -1, 0],\
                      [ 0, 0,  0, 0],\
                      [-1, 0,  1, 0],\
                      [ 0, 0,  0, 0]])

  # rotation matrix
  def rotation_matrix(alpha):
    return np.array([[ np.cos(alpha), np.sin(alpha),             0,             0],\
                     [-np.sin(alpha), np.cos(alpha),             0,             0],\
                     [             0,             0, np.cos(alpha), np.sin(alpha)],\
                     [             0,             0,-np.sin(alpha), np.cos(alpha)]])

  # assemble global stiffness matrix
  for (element_no, element) in enumerate(elements):
    
    # compute angle of member
    node0_pos = node_positions[element[0]]
    node1_pos = node_positions[element[1]]
    alpha = np.arctan2(node1_pos[1] - node0_pos[1], node1_pos[0] - node0_pos[0])
    
    # compute length of member
    l = np.linalg.norm(np.array(node1_pos)-np.array(node0_pos))
      
    # rotate local stiffness matrix
    k_local_rotated = rotation_matrix(alpha).transpose().dot(k_local*ea/l).dot(rotation_matrix(alpha))
    
    print "   - element ", element_no,", angle ",180./np.pi*alpha,", length ",l
    if False:
      print "k_local:"
      print k_local
      print "rotation_matrix:"
      print rotation_matrix(alpha)
      print "k_local_rotated:"
      print k_local_rotated
    
    
    # add contribution of rotated local stiffness matrix to global stiffness matrix
    for i in range(4):
      for j in range(4): 
        dof_i = int(id_matrix[i, element_no])
        dof_j = int(id_matrix[j, element_no])
        k_global[dof_i, dof_j] += k_local_rotated[i, j]

  # assemble load vector
  f = np.zeros([number_dofs, 1])
  for load in loads:
    f[load[0]] = load[1]
    
  print ""
  print "global stiffness matrix k_global:"
  print k_global
  print ""
  print "global right hand side load vector f:"
  print f
  print ""

  # remove entries that are fixed by boundary conditions
  f_reduced = np.delete(f, fixed_dofs, 0)
  k_global_reduced = np.delete(k_global, fixed_dofs, 0)
  k_global_reduced = np.delete(k_global_reduced, fixed_dofs, 1)

  print "----- reduce system by removing fixed dofs -----"
  print "k_global_reduced:"
  print k_global_reduced
  print ""
  print "f reduced:"
  print f_reduced

  print ""
  print "----- solve system k*u = f for u -----"
  print "number of nonfixed dofs: ", number_nonfixed_dofs, ", rank of stiffness matrix: ",np.linalg.matrix_rank(k_global_reduced)
  print "(these numbers should be equal, otherwise the system is not statically determined)"

  if number_nonfixed_dofs > np.linalg.matrix_rank(k_global_reduced):
    print "Matrix rank is too low. No solution could be found."
    return (None, id_matrix)

  # solve system
  try:
    u_reduced = np.linalg.solve(k_global_reduced, f_reduced)
  except:
    print "Linear system has no solution - system in not statically determined!"
    return (None, id_matrix)
    #draw_truss.draw(node_positions, elements, fixed_dofs, loads, id_matrix, np.zeros([number_dofs]), 'k', 3, 0, False)
    #plt.show()
    #sys.exit(0)

  # recover u from u_reduced
  u = np.zeros([number_dofs])
  for dof in range(number_dofs):
    if not (dof in fixed_dofs):
      reduced_index = dof - sum([fixed_dof < dof for fixed_dof in fixed_dofs])
      #print [fixed_dof < dof for fixed_dof in fixed_dofs]
      #print "reduced_index:", reduced_index
      u[dof] = u_reduced[reduced_index]

  print "solution (u reduced):"
  print u_reduced
  print ""
  print "----- recover full sized u by inserting 0's -----"
  print ", u recovered:"
  print u
  print ""
  
  return (u, id_matrix)
  
def solve_truss_dynamic(node_positions, elements, loads, fixed_dofs, options):
  """
    solve the truss as dynamic system if it is hypostatic, if a solution exists, returns vector of nodal displacements (u,id_matrix)
    if no solution is found, returns (None,id_matrix)
  """
  global figure_cntr
  
  output = False

  ea = options["ea"] if "ea" in options else 1
  rho = options["rho"] if "rho" in options else 1
  scenario = options["scenario"] if "scenario" in options else 1
  n_ode_iterations = options["n_ode_iterations"] if "n_ode_iterations" in options else 5
  t_end = options["t_end"] if "t_end" in options else 1.0
  initial_velocity_reduced = options["initial_velocity_reduced"] if "initial_velocity_reduced" in options else None
  damping_factor = options["damping_factor"] if "damping_factor" in options else 0.0

  if output:
    print ""
    print "======= function solve truss with {} elements =======".format(len(elements))

  number_nodes = len(node_positions)
  number_elements = len(elements)
  number_dofs = number_nodes * 2    # every nodes has two displacements (degrees of freedom = dof)
  number_nonfixed_dofs = number_dofs - len(fixed_dofs)

  if output:
    print "---------- Input parameters ------------"
    print len(node_positions), "node_positions [x, y]:          ",node_positions
    print len(elements),       "elements [node start, node end]:",elements
    print len(loads),          "loads [dof, value]:             ",loads
    print len(fixed_dofs),     "fixed_dofs:                     ",fixed_dofs
    print "parameter EA =",ea

  # create ID-matrix
  id_matrix = np.zeros([4, number_elements], int)
  for (index, element) in enumerate(elements):
    id_matrix[0, index] = element[0]*2 + 0
    id_matrix[1, index] = element[0]*2 + 1
    id_matrix[2, index] = element[1]*2 + 0
    id_matrix[3, index] = element[1]*2 + 1

  if output:
    print ""
    print "id_matrix:"
    print id_matrix
  
  # evaluate counting criterion
  value = 2*number_nodes - (len(fixed_dofs) + number_elements)
  
  if output:
    print "counting criterion: 2*{} - ({} + {}) = {}".format(number_nodes,len(fixed_dofs),number_elements,value)
    if value == 0:
      print " -> statically determined (statisch bestimmt)"
    elif value < 0:
      print " -> hyperstatic (statisch überbestimmt)"
    else:
      print " -> hypostatic (statisch unterbestimmt, verschieblich)"

  u = np.zeros([number_dofs])   # vector of displacements, this is seeked, it is the solution of k_global * u = f
  f = np.zeros([number_dofs])   # vector of loads, this is the right hand side
  m_global = np.zeros([number_dofs, number_dofs])   # the global mass matrix

  if initial_velocity_reduced is None:
    initial_velocity_reduced = np.zeros(number_nonfixed_dofs)   # the initial velocities
    initial_velocity_reduced.shape = (number_nonfixed_dofs,1)

  if output:
    print ""
    print "----- assemble global stiffness matrix from local matrices -----"

  # local stiffness matrix
  k_local = np.array([[ 1, 0, -1, 0],\
                      [ 0, 0,  0, 0],\
                      [-1, 0,  1, 0],\
                      [ 0, 0,  0, 0]])

  # rotation matrix
  def rotation_matrix(alpha):
    return np.array([[ np.cos(alpha), np.sin(alpha),             0,             0],\
                     [-np.sin(alpha), np.cos(alpha),             0,             0],\
                     [             0,             0, np.cos(alpha), np.sin(alpha)],\
                     [             0,             0,-np.sin(alpha), np.cos(alpha)]])
  def rotation_vector(alpha):
    return np.array([[ np.cos(alpha), np.sin(alpha)],\
                     [-np.sin(alpha), np.cos(alpha)]])

  # assemble global mass matrix
  for (element_no, element) in enumerate(elements):
  
    # compute length of member
    node0_pos = np.array(node_positions[element[0]])
    node1_pos = np.array(node_positions[element[1]])
      
    l = np.linalg.norm(node1_pos - node0_pos)
    
    for i in range(4):
      dof_i = int(id_matrix[i, element_no])
      m_global[dof_i, dof_i] += l/2.*rho

  # assemble global stiffness matrix
  def global_stiffness_matrix(u):
    k_global = np.zeros([number_dofs, number_dofs])   # the global stiffness matrix
    
    for (element_no, element) in enumerate(elements):
    
      # compute angle of member
      node0_pos = np.array(node_positions[element[0]], dtype='float')
      node1_pos = np.array(node_positions[element[1]], dtype='float')
      
      alpha = np.arctan2(node1_pos[1] - node0_pos[1], node1_pos[0] - node0_pos[0])
      
      # reference length
      if len(element) == 2:
        L = np.linalg.norm(node1_pos - node0_pos)
        elements[element_no] = element + [L]
      else:
        L = element[2]
      
      # rotate local stiffness matrix
      k_local_rotated = rotation_matrix(alpha).transpose().dot(k_local*ea/L).dot(rotation_matrix(alpha))
      
      if output or False:
        print "   - (lin) element ", element_no,", angle ",180./np.pi*alpha,", length l=",l,", L=",L
        if True:
          print "k_local:"
          print k_local
          print "rotation_matrix:"
          print rotation_matrix(alpha)
          print "k_local_rotated:"
          print k_local_rotated
      
      
      # add contribution of rotated local stiffness matrix to global stiffness matrix
      for i in range(4):
        dof_i = int(id_matrix[i, element_no])
        for j in range(4): 
          dof_j = int(id_matrix[j, element_no])
          k_global[dof_i, dof_j] += k_local_rotated[i, j]
        
    return k_global

  def compute_distance_to_initial_lengths(u):
	  
    d_initial_lengths = np.zeros((number_dofs,1))
    
    for (element_no, element) in enumerate(elements):
      
      # compute angle of member
      node0_pos = np.array(node_positions[element[0]], dtype='float')
      node1_pos = np.array(node_positions[element[1]], dtype='float')
      
      alpha = np.arctan2(node1_pos[1] - node0_pos[1], node1_pos[0] - node0_pos[0])
      
      # reference length
      if len(element) == 2:
        L = np.linalg.norm(node1_pos - node0_pos)
        elements[element_no] = element + [L]
      else:
        L = element[2]
      
      
      u_node0 = np.array([float(u[id_matrix[0, element_no]]), float(u[id_matrix[1, element_no]])])
      u_node1 = np.array([float(u[id_matrix[2, element_no]]), float(u[id_matrix[3, element_no]])])
      
      node0_pos += u_node0
      node1_pos += u_node1
      
      
      # current length
      l = np.linalg.norm(node1_pos - node0_pos)
      uu = l-L
      
      d_initial_lengths[id_matrix[0, element_no]] = -uu/2. * np.cos(alpha)
      d_initial_lengths[id_matrix[1, element_no]] = -uu/2. * np.sin(alpha)
      
      d_initial_lengths[id_matrix[2, element_no]] = uu/2. * np.cos(alpha)
      d_initial_lengths[id_matrix[3, element_no]] = uu/2. * np.sin(alpha)
      
    return d_initial_lengths

  def global_stiffness_matrix_times_u(u):
    
    debug = False
    
    ku_global = np.zeros((number_dofs,1))
    
    for (element_no, element) in enumerate(elements):
      
      # compute angle of member
      node0_pos = np.array(node_positions[element[0]], dtype='float')
      node1_pos = np.array(node_positions[element[1]], dtype='float')
      
      alpha = np.arctan2(node1_pos[1] - node0_pos[1], node1_pos[0] - node0_pos[0])
      
      
      # reference length
      if len(element) == 2:
        L = np.linalg.norm(node1_pos - node0_pos)
        elements[element_no] = element + [L]
        if debug:
          print "store L=",L
      else:
        L = element[2]
        if debug:
          print "retrieve L=",L
      
      
      u_node0 = np.array([float(u[id_matrix[0, element_no]]), float(u[id_matrix[1, element_no]])])
      u_node1 = np.array([float(u[id_matrix[2, element_no]]), float(u[id_matrix[3, element_no]])])
      
      node0_pos += u_node0
      node1_pos += u_node1
      
      
      # current length
      l = np.linalg.norm(node1_pos - node0_pos)
      uu = l-L
      
      # strain
      eps = uu/L
      
      # Euler-Lagrange strain E = 0.5*(l**2 - L**2) / (L**2) = 0.5*((L+u)**2 - L**2) / (L**2) 
      # E = 0.5*(((L+u)/L)**2 - 1) = 0.5*((1+eps)**2 - 1) = eps + 1./2*eps**2
      E = 0.5*(l**2 - L**2) / (L**2)
      
      # natural strain
      e = np.log(1 + eps)
      
      # force
      #f = ea*e
      f = ea*eps
      
      # add contribution to global ku
      ku_global[id_matrix[0, element_no]] += -f * np.cos(alpha)
      ku_global[id_matrix[1, element_no]] += -f * np.sin(alpha)
      
      ku_global[id_matrix[2, element_no]] += f * np.cos(alpha)
      ku_global[id_matrix[3, element_no]] += f * np.sin(alpha)
      
      if output or False:
        print "   - (nli) element ", element_no,", angle ",180./np.pi*alpha,", length l=",l,"L=",L
        print "     u=",uu,",u*ea/L=f=",f
        print "     add ku entry ",id_matrix[0, element_no]," value ",-f * np.cos(alpha)
        print "     add ku entry ",id_matrix[1, element_no]," value ",-f * np.sin(alpha)
        print "     add ku entry ",id_matrix[2, element_no]," value ",f * np.cos(alpha)
        print "     add ku entry ",id_matrix[3, element_no]," value ",f * np.sin(alpha)
        if True:
          print "strain:"
          print e
          
    return ku_global

  #k_global = global_stiffness_matrix(np.zeros(number_dofs))
  k_global = np.zeros([number_dofs, number_dofs])

  # assemble load vector
  f = np.zeros([number_dofs, 1])
  for load in loads:
    f[load[0]] = load[1]
    
  if output:
    print ""
    print "global stiffness matrix k_global:"
    print k_global
    print ""
    print "global mass matrix m_global:"
    print m_global
    print ""
    print "global right hand side load vector f:"
    print f
    print ""

  # remove entries that are fixed by boundary conditions
  def reduced(matrix):
    """ Remove rows and columns in matrix that correspond to fixed dofs """
    matrix_reduced = np.delete(matrix, fixed_dofs, 0)
    matrix_reduced = np.delete(matrix_reduced, fixed_dofs, 1)
    return matrix_reduced
  
  def reduced_vector(vector):
    """ Remove entries in vector that correspond to fixed dofs """
    vector_reduced = np.delete(vector, fixed_dofs, 0)
    return vector_reduced
    
  def recovered(vector_reduced):
    """ Extend vector with number_nonfixed_dofs entries to full vector with number_dofs entries
        New entries have the value 0
    """
    result = np.zeros([number_dofs])
    
    for dof in range(number_dofs):
      if not (dof in fixed_dofs):
        reduced_index = dof - sum([fixed_dof < dof for fixed_dof in fixed_dofs])
        #print [fixed_dof < dof for fixed_dof in fixed_dofs]
        #print "reduced_index:", reduced_index
        result[dof] = vector_reduced[reduced_index]
    return result
    
  f_reduced = reduced_vector(f)
  m_global_reduced = reduced(m_global)
  k_global_reduced = reduced(k_global)

  if output:
    print "----- reduce system by removing fixed dofs -----"
    #print "k_global_reduced:"
    #print k_global_reduced
    print ""
    print "m_global_reduced:"
    print m_global_reduced
    print ""
    print "f reduced:"
    print f_reduced


  # compute the inverse matrix of the reduced mass matrix, m_reduced_inv
  m_reduced_inv = np.linalg.inv(m_global_reduced)

  if output:
    print ""
    print "m_reduced_inv:"
    print m_reduced_inv

    print ""
    print "----- solve ODE system M*u'' = -k*u + -----"
    print "number of nonfixed dofs: ", number_nonfixed_dofs, ", rank of stiffness matrix: ",np.linalg.matrix_rank(k_global_reduced)

    if number_nonfixed_dofs > np.linalg.matrix_rank(k_global_reduced):
      print "Matrix rank {} is lower than number of nonfixed dofs ({}).".format(np.linalg.matrix_rank(k_global_reduced), number_nonfixed_dofs)

  damping = damping_factor

  # solve system until end time t_end
  def rhs(w, t, number_nonfixed_dofs):
    
    debug = False
    
    d = np.array(w[0:number_nonfixed_dofs])
    v = np.array(w[number_nonfixed_dofs:])
    
    # transform to column vectors
    d.shape = (number_nonfixed_dofs, 1)
    v.shape = (number_nonfixed_dofs, 1)
    
    mode = "nonlinear" # linear does not work like this
    
    # linear, i.e. with explicit stiffness_matrix*vector multiplication
    if mode == "linear":
      k_global = global_stiffness_matrix(recovered(d))
      k_global_reduced = reduced(k_global)
      
      # the displacement vector needs to be adjusted to have the displacement from a setting with initial lengths
      d_from_initial_lengths = compute_distance_to_initial_lengths(recovered(d))
      
      if debug:
        print ""
        print "lin ku:",k_global_reduced.dot(d),", k_global_reduced:",k_global_reduced,"u:",d, "from initial lengths:",d_from_initial_lengths
        
      ku = k_global_reduced.dot(reduced_vector(d_from_initial_lengths))
    else:
      # nonlinear, i.e. stiffness matrix is not explicitly assembled, but k*u is computed directly
      ku = reduced_vector(global_stiffness_matrix_times_u(recovered(d)))
      
      if debug:
        print "nli ku:", ku
    rhs = (-damping*v - ku + f_reduced)
    
    if debug:
      print "ku: ",ku
      print "f_reduced: ", f_reduced
      print "rhs:  ",rhs
    
    rhs = m_reduced_inv.dot(rhs)
    
    rhs.shape = (number_nonfixed_dofs)
    v.shape = (number_nonfixed_dofs)
    
    if debug:
      print "acceleration rhs: ", rhs
      print "velocity v:",v
    
    result = np.concatenate((v, rhs))
    return result
    
  # ODE solver parameters
  abserr = 1.0e-10
  relerr = 1.0e-8
  stoptime = t_end
  n_iterations = options["n_ode_iterations"]

  # Create the time samples for the output of the ODE solver.
  t = [stoptime * float(i) / (n_iterations - 1) for i in range(n_iterations)]

  initial_velocity_reduced.shape = (number_nonfixed_dofs)

  w0 = np.concatenate((np.zeros(number_nonfixed_dofs),initial_velocity_reduced))
  
  if output:
    print "stoptime=",stoptime,",initial values:",w0
    print "solve for timesteps ",t
  
  wsol = scipy.integrate.odeint(rhs, w0, t, args=(number_nonfixed_dofs,),atol=abserr, rtol=relerr)

  if output:  
    print "wsol:",wsol
  
  u_reduced = wsol[-1][0:number_nonfixed_dofs]
  u_reduced.shape = (number_nonfixed_dofs,1)      # make column vector
  velocity_reduced = wsol[-1][number_nonfixed_dofs:]


  # recover u from u_reduced
  u = recovered(u_reduced)

  if output:
    print "solution (u reduced):"
    print u_reduced
    print ""
    print "----- recover full sized u by inserting 0's -----"
    print ", u recovered:"
    print u
    print ""
    
  if output:
    print "k_00:", k_global_reduced[0,0],", k_11:", k_global_reduced[1,1],", k_01:",k_global_reduced[0,1],k_global_reduced[1,0]
    
    print "k*d: "
    print k_global_reduced.dot(u_reduced)
    print "f:"
    print f_reduced
    print "m*d'' + D*d + k*d = f    => f-k*d = m*d'' + D*d,  d'' = (f-k*d)/m:"
    print m_reduced_inv*(f_reduced - k_global_reduced.dot(u_reduced))
  
  #  (k_00 k_01)  (u0) =  (f0)         k_00*u0   +   k_01*u1   =   f0,  k_00*(u0p + u0) + k_01*(u1p + u1) = f0
  #  (k_10 k_11)  (u1)    (f1)
  
  # k_00*(u0p + u0) + k_01*(u1p + u1) = f0,  k_00*u0p + k_01*u1p    + k_00*u0 + k_01*u1 = f0
  
  result = dict()
  result["u"] = u
  result["velocity_reduced"] = velocity_reduced
  result["id_matrix"] = id_matrix
  result["u_reduced"] = u_reduced
  result["ku"] =  k_global_reduced.dot(u_reduced)
  
  return result
  
