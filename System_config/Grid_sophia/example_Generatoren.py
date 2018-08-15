#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 16:47:50 2018

@author: sophia
"""

#ranodm Generator
def abc_generator():
    yield("a")
    print('hello hello there')
    yield("b")
    yield("c")

x = abc_generator()
print(x.next())
print(x.next())
print(x.next())

x = abc_generator()
print(x.next()) 
#-----------------------------------------------
#fibonacci Generator 
def fibonacci(n):
    ''' ein fibonacci zahlengenerator '''
    a,b, counter = 0,1,0
    while True:
        if (counter > n): return
        yield a
        a,b = b, a+b
        counter += 1
f = fibonacci(5)
for k in f:
    print(k),
#-----------------------------------------------
# permutation Generator
def permutaions(items):
    n = len(items)
    if n == 0: 
        yield [] 
    else: 
        for i in range(len(items)): 
            print i
            for cc in permutaions(items[:i] + items[i+1:]):
                print cc
                yield[items[i]]+cc

for p in permutaions(['r', 'e', 'd']): print ' '.join(p)
for p in permutaions(list('game')): print ' '.join(p)

