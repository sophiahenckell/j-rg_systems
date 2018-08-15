# yield_exampe.py example of yield, return in generator function
import simpy

def fy():
	print("before fy")
	x = 1
	print("after x=1")
	y = 1
	print("after y= 1")
	yield x,y,x+y
	print("after yield x,y,x+y,1")
	z = 6
	print("after z = 12")
	yield z/x
	print("after yield z/x 1",z/x) 

def gy():
	print("before gy")
	x = 2
	print("after x = 2")
	y = 2
	print("after y = 2")
	yield x,y,x+y
	print("after yield x,y,x+y")
	z=12
	print("after z = 12")
	yield z/x
	print("after yield z/x 2",z/2)

def main():
	f = fy()
	g = gy()
	print(f.__next__())
	print(g.__next__())
	print(f.__next__())
	print(g.__next__())

if __name__ == '__main__':
	main()
