"""
Python Quaternion Functions
Written By: Dean Keithly
Written On: 06/25/2019

References to (x.xx) are from Markley and Crassidis Fundamentals of Spacecraft Attitude Determination and Control

Variables of Form:
x = np.asarray([[0.,0.,1.]])
e = np.asarray([[1.,0.,0.]])
theta = np.pi/2.
"""
import numpy as np


def cross_product_matrix(Y):
    """  [Y x] (2.55) where x here is a cross product pg 29
    Y (array) - 3x1 matrix
    cpm (array) - 3x3 cross product matrix
    """
    cpm = np.asarray([  [0,    -Y[2],Y[1]],\
                        [Y[2], 0,    -Y[0]],\
                        [-Y[1],Y[0], 0]])
    return cpm

def q_otimes(q):
    """ [q otimes] (2.85)
    Args:
        q (array) - 1x4 array
    Return:
        qotimes (array) - 4x4 quaternion multiplication array
    """
    qotimes = np.hstack((q[0][3]*np.identity(3) - cross_product_matrix(q[0][0:3]),q.T[0:3]))
    bottom_row = np.hstack((-q[0][0:3],q[0][3]))
    qotimes = np.vstack((qotimes,bottom_row))
    return qotimes

def q_odot(q):
    """ [q odot] (2.86)
    Args:
        q (array) - 1x4 quaternion
    Return:
        qodot (array) - 4x4 quaternion multiplication array
    """
    qodot = np.hstack((q[0][3]*np.identity(3) + cross_product_matrix(q[0][0:3]),q.T[0:3]))
    bottom_row = np.hstack((-q[0][0:3],q[0][3]))
    qodot = np.vstack((qodot,bottom_row))
    return qodot

def q_conjugate(q):
    """ (2.91)
    Args:
        q (array) - 1x4 array
    Returns:
        q_conjugate (array) - 1x4 array
    """
    q_conjugate = q
    q_conjugate[0:3] = -q[0:3]
    return q_conjugate

def eulerRotationMatrix(e,theta):
    """ (2.110)
    Args:
        e (array) - 1x3 array to rotate about
        theta (float) - angle to rotate about in radians
    Return:
        rot_mat (array) - 3x3 rotation matrix
    """
    rot_mat = np.asarray([[np.cos(theta)+(1.-np.cos(theta))*e[0]**2.,        (1.-np.cos(theta))*e[0]*e[1]+np.sin(theta)*e[2],    (1.-np.cos(theta))*e[0]*e[2]-np.sin(theta)*e[1]],\
                          [(1.-np.cos(theta))*e[1]*e[0]-np.sin(theta)*e[2],  np.cos(theta)+(1.-np.cos(theta))*e[1]**2.,          (1.-np.cos(theta))*e[1]*e[2]+np.sin(theta)*e[0]],\
                          [(1.-np.cos(theta))*e[2]*e[0]+np.sin(theta)*e[1],  (1.-np.cos(theta))*e[2]*e[1]-np.sin(theta)*e[0],    np.cos(theta)+(1.-np.cos(theta))*e[2]**3.]])
    return rot_mat

def q_rep(e,theta):
    """ (2.124)
    Args:
        e (array) - 1x3 array of vector to rotate about
        theta (scalar) - angle to rotate about
    Returns:
        q (array) - 1x4 array
    """
    q = np.vstack((np.sin(theta/2.)*e.T,np.cos(theta/2.)))
    return q.T

def rot_quaternion(q,x):
    """ (2.128)
    Args:
    Returns:
        x (array) - new position after rotation
    """
    y = np.matmul(np.matmul(q_odot(q).T,q_otimes(q)),np.vstack((x.T,np.asarray([[0.]]))))
    return y[0:3]

