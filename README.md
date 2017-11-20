Manipulator Performance Constraints for human-robot Cooperation

# Description:
Use any local performance index of a robotic manipulator to constrain the robot above a certain threshold.
This is particularly useful for singularity avoidance in physical human-robot interaction

A detailed description can be found in the following RCIM paper:
- Dimeas, Fotios, Vassilis C. Moulianitis, and Nikos Aspragathos. "Manipulator performance constraints in human-robot cooperation." Robotics and Computer-Integrated Manufacturing (2017).

and also in ICRA 2016 (early version):
- Dimeas, Fotios, et al. "Manipulator performance constraints in Cartesian admittance control for human-robot cooperation." Robotics and Automation (ICRA), 2016 IEEE International Conference on. IEEE, 2016.

A video demonstration: 
<a href="http://www.youtube.com/watch?feature=player_embedded&v=1zTDmiDjDOA
" target="_blank"><img src="http://img.youtube.com/vi/1zTDmiDjDOA/0.jpg" 
alt="IMAGE ALT TEXT HERE" width="240" height="180" border="10" /></a>


This repository includes the Matlab code that was used in the ICRA 2016 paper and the extended method in C++ that is presented in the RCIM journal. 

The C++ implementation includes performance constraints in both translational and rotational axes. It also includes 3 ways of calculating them:
* Serial: Each direction is calculated after the other. It might take a while
* Parallel: Each direction is calculated in its own thread worker. Function returns when all calculations are over. This usually works in less than 1ms. Better have a multi-thread CPU.
* Parallel non-blocking: Same as above but the function returns immediately. This can work in slow computers. WARNING! Be very careful when you use it, synchronization issues of the constraint forces may arise.

# Usage:
## C++
> g++ demo.cpp performanceConstraints.cpp Jacobian.cpp  -o demo -larmadillo -pthread -std=c++11

The Jacobian matrix that is provided in symbolic form is for the KUKA LWR 4+

### Requirements:
- Armadillo library

## Matlab
1. Set the desired parameters for the robot and the performance
constraints
2. Run the simulation
3. View the simulated motion of the robot and the plots

### Requirements:
- Robotics Toolbox for Matlab (http://www.petercorke.com/Robotics_Toolbox.html)

### Notes: 
- The simulation does not consider the joint limits so the robot
might behave weird
- This code has been tested with Robotics Toolbox 9.10 and Matlab R2014a


Authors: Fotios Dimeas, Charalambos Papakonstantinou

Copyright 2015-2017 Fotios Dimeas
