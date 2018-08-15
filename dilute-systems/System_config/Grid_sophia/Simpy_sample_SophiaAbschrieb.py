#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 14:27:36 2018

@author: sophia
"""

import simpy
def fy ():
    x,y = 1,4
    yield x, y, x+y
    print ("after yield")
    
    z = 6
    print("after z = 6")
    yield z/y
    print("after yield z/x 1", z/x)
    

def main():
    f = fy()
    print( f.next())
    print( f.next())
if __name__ == '__main__':
    main()

def gy ():
    print( "before gy")