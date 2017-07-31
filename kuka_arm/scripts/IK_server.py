#!/usr/bin/env python

# Copyright (C) 2017 Electric Movement Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *

# Define DH param symbols for Joint angle
q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8') #theta_i
d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')

# Modified DH params
s = {alpha0:     0,  a0:      0, d1:  0.75,
     alpha1: -pi/2,  a1:   0.35, d2:     0,  q2: q2-pi/2,
     alpha2:     0,  a2:   1.25, d3:     0,
     alpha3: -pi/2,  a3: -0.054, d4:  1.50,
     alpha4:  pi/2,  a4:      0, d5:     0,
     alpha5: -pi/2,  a5:      0, d6:     0,
     alpha6:     0,  a6:      0, d7: 0.303,  q7: 0}

# Define Modified DH Transformation matrix
T0_1 = Matrix([[             cos(q1),            -sin(q1),            0,              a0],
               [ sin(q1)*cos(alpha0), cos(q1)*cos(alpha0), -sin(alpha0), -sin(alpha0)*d1],
               [ sin(q1)*sin(alpha0), cos(q1)*sin(alpha0),  cos(alpha0),  cos(alpha0)*d1],
               [                   0,                   0,            0,               1]])
T0_1 = T0_1.subs(s)

T1_2 = Matrix([[             cos(q2),            -sin(q2),            0,              a1],
               [ sin(q2)*cos(alpha1), cos(q2)*cos(alpha1), -sin(alpha1), -sin(alpha1)*d2],
               [ sin(q2)*sin(alpha1), cos(q2)*sin(alpha1),  cos(alpha1),  cos(alpha1)*d2],
               [                   0,                   0,            0,               1]])
T1_2 = T1_2.subs(s)

T2_3 = Matrix([[             cos(q3),            -sin(q3),            0,              a2],
               [ sin(q3)*cos(alpha2), cos(q3)*cos(alpha2), -sin(alpha2), -sin(alpha2)*d3],
               [ sin(q3)*sin(alpha2), cos(q3)*sin(alpha2),  cos(alpha2),  cos(alpha2)*d3],
               [                   0,                   0,            0,               1]])
T2_3 = T2_3.subs(s)

T3_4 = Matrix([[             cos(q4),            -sin(q4),            0,              a3],
               [ sin(q4)*cos(alpha3), cos(q4)*cos(alpha3), -sin(alpha3), -sin(alpha3)*d4],
               [ sin(q4)*sin(alpha3), cos(q4)*sin(alpha3),  cos(alpha3),  cos(alpha3)*d4],
               [                   0,                   0,            0,               1]])
T3_4 = T3_4.subs(s)

T4_5 = Matrix([[             cos(q5),            -sin(q5),            0,              a4],
               [ sin(q5)*cos(alpha4), cos(q5)*cos(alpha4), -sin(alpha4), -sin(alpha4)*d5],
               [ sin(q5)*sin(alpha4), cos(q5)*sin(alpha4),  cos(alpha4),  cos(alpha4)*d5],
               [                   0,                   0,            0,               1]])
T4_5 = T4_5.subs(s)

T5_6 = Matrix([[             cos(q6),            -sin(q6),            0,              a5],
               [ sin(q6)*cos(alpha5), cos(q6)*cos(alpha5), -sin(alpha5), -sin(alpha5)*d6],
               [ sin(q6)*sin(alpha5), cos(q6)*sin(alpha5),  cos(alpha5),  cos(alpha5)*d6],
               [                   0,                   0,            0,               1]])
T5_6 = T5_6.subs(s)

T6_G = Matrix([[             cos(q7),            -sin(q7),            0,              a6],
               [ sin(q7)*cos(alpha6), cos(q7)*cos(alpha6), -sin(alpha6), -sin(alpha6)*d7],
               [ sin(q7)*sin(alpha6), cos(q7)*sin(alpha6),  cos(alpha6),  cos(alpha6)*d7],
               [                   0,                   0,            0,               1]])
T6_G = T6_G.subs(s)

# Create individual transformation matrices
T0_2 = (T0_1 * T1_2)  # baselink to link2
T0_3 = (T0_2 * T2_3)  # baselink to link3
T0_4 = (T0_3 * T3_4)  # baselink to link4
T0_5 = (T0_4 * T4_5)  # baselink to link5
T0_6 = (T0_5 * T5_6)  # baselink to link6
T3_6 = (T3_4 * T4_5 * T5_6)
T0_G = (T0_6 * T6_G)  # baselink to linkG
#print("T0_G = ",T0_G)
#print('T0_3 = ', T0_3.evalf(subs={q1: 0, q2: 0, q3: 0}))
#print('T3_6 = ', T3_6.evalf(subs={q4: 0, q5: 0, q6: 0})) 
#print("\n")

### Correction needed to account for orientation difference between
#    definition of Gripper link in urdf vs. DH Convention
R_z = Matrix([[   cos(pi), -sin(pi),         0, 	0],
              [   sin(pi),  cos(pi),         0, 	0],
              [  	0,        0,         1, 	0],
              [ 	0,        0,         0, 	1]])

R_y = Matrix([[ cos(-pi/2),	      0, sin(-pi/2), 	0],
              [          0,	      1,          0, 	0],
              [-sin(-pi/2),	      0, cos(-pi/2), 	0],
              [	         0,	      0,          0, 	1]])

R_corr = (R_z * R_y)
#print('R_corr = ', R_corr)

### Numerically evaluate transforms (compare to tf_echo output)
#print("T0_1 = ",T0_1.evalf(subs={q1: 0}))
#print("T0_2 = ",T0_2.evalf(subs={q1: 0, q2: 0}))
#print("T0_3 = ",T0_3.evalf(subs={q1: 0, q2: 0, q3: 0}))
#print("T0_4 = ",T0_4.evalf(subs={q1: 0, q2: 0, q3: 0, q4: 0}))
#print("T0_5 = ",T0_5.evalf(subs={q1: 0, q2: 0, q3: 0, q4: 0, q5: 0}))
#print("T0_6 = ",T0_6.evalf(subs={q1: 0, q2: 0, q3: 0, q4: 0, q5: 0, q6: 0}))
#print("T0_G = ",T0_G.evalf(subs={q1: 0, q2: 0, q3: 0, q4: 0, q5: 0, q6: 0}))
#print("\n")

### Total Homogeneous Transform between baselink and Gripper link
#    with orientation correction applied
T_total = simplify(T0_G * R_corr)
#print('T_total = ', T_total.evalf(subs={q1: 0, q2: 0, q3: 0, q4: 0, q5: 0, q6: 0}))
#print('T_total = ', T_total)


def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        # Initialize service response
        joint_trajectory_list = []
        for x in xrange((len(req.poses)-3), len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

            # Extract end-effector position = px,py,pz
            px = (req.poses[x].position.x - 0.325) #simple correction for Wrist Center
            py = req.poses[x].position.y
            pz = req.poses[x].position.z
	    #print('px = ',px)	
	    #print('py = ',py)	
	    #print('pz = ',pz)	
	    #print("\n")

            # Calculate joint angles using Geometric IK method
	    rtd = 180/pi  #radians to degrees

            #Solve for end-effector POSITION - theta 1, 2, 3
	    theta1 = atan2(py, px) 
            #th1deg = simplify(theta1*rtd)
	    #print('theta1 = ', th1deg)

            pf = (px**2 + py**2)**0.5
	    theta2 = (atan2(pf,pz)-acos((1.5**2-1.196**2-(pf-0.35)**2-(pz-0.75)**2)
			/(-2*1.5*(((pf-0.35)**2+(pz-0.75)**2)**0.5))))+.22
            #th2deg = theta2*rtd
	    #print('theta2 = ', simplify(th2deg))

            D = (((pf-0.35)**2 + (pz-0.75)**2 - 1.196**2 - 1.5**2)/(2*1.196*1.5))
	    theta3 = -pi/2 + acos(D)+.2
            #th3deg = theta3*rtd
	    #print('theta3 = ', simplify(th3deg))
	    #print("\n")


            #Decouple the wrist to solve for end-effector ORIENTATION - theta 4, 5, 6
            R0_3 = (T0_3).evalf(subs={q1:theta1,q2:theta2,q3:theta3}) 
            #print'R0_3 = ', simplify(R0_3)
	    #print("\n")

            I0_3 = R0_3**-1
            #print'I0_3 = ', simplify(I0_3)
	    #print("\n")

            # Extract end-effector orientation = roll, pitch, yaw  
            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, 
		 req.poses[x].orientation.y,
                 req.poses[x].orientation.z, 
		 req.poses[x].orientation.w])
	    #print'roll = ',roll, 'pitch = ',pitch, 'yaw = ',yaw
	    #print("\n")

            R_roll =  Matrix([[ 1,         0,          0, 0],
                              [ 0, cos(roll), -sin(roll), 0],
                              [ 0, sin(roll),  cos(roll), 0],
                              [ 0,         0,          0, 1]])

            R_pitch = Matrix([[  cos(pitch),    0, sin(pitch), 0],
                              [           0,    1,          0, 0],
                              [ -sin(pitch),    0, cos(pitch), 0],
                              [           0,    0,          0, 1]])

            R_yaw  =  Matrix([[ cos(yaw), -sin(yaw),    0, 0],
                              [ sin(yaw),  cos(yaw),    0, 0],
                              [        0,         0,    1, 0],
                              [        0,         0,    0, 1]])

            Rrpy = R_yaw * R_pitch * R_roll * R_corr
            #print'Rrpy = ', simplify(Rrpy)
	    #print("\n")
           
            R3_6 = R0_3.inv() * Rrpy  #R0_3**-1
            #print'R3_6 = ', simplify(R3_6)
	    #print("\n")

            r11 = R3_6[0,0]
            r12 = R3_6[0,1]
            r13 = R3_6[0,2]
            r21 = R3_6[1,0]
            r22 = R3_6[1,1]
            r23 = R3_6[1,2]
            r31 = R3_6[2,0]
            r32 = R3_6[2,1]
            r33 = R3_6[2,2]

            theta4 = atan2(r33, -r13) 
            #th4deg = theta4*rtd
	    #print'theta4 = ', simplify(th4deg)

            theta5 = atan2(sqrt(r13**2 + r33**2),r23) 
            #th5deg = theta5*rtd
	    #print'theta5 = ', simplify(th5deg)

            theta6 = atan2(-r22, r21) 
            #th6deg = theta6*rtd
	    #print'theta6 = ', simplify(th6deg)

            print "FW Kinematics Result:"
            print(T_total.evalf(subs={q1:theta1, q2:theta2, q3:theta3, q4:theta4, q5:theta5, q6:theta6}))

            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
	    joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
	    joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
