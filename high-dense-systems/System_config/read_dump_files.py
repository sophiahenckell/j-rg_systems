import numpy as np
f = open('output/elli2.dump','r')
print("this is for 512 particle system")
#read timestep zero first

# f.readline()
# time = f.readline()
# time = time.strip()
# time = float(time)
# headeratoms = f.readline()
# natoms = f.readline()
# natoms = natoms.strip()
# natoms = int(natoms)
# headerbox = f.readline()
# Lx = f.readline()
# Lx = Lx.strip()
# Lx = float(Lx.split()[1]) - float(Lx.split()[0])
# Ly = f.readline()
# Ly = Ly.strip()
# Ly= float(Ly.split()[1]) - float(Ly.split()[0])
# Lz = f.readline()
# Lz = Lz.strip()
# Lz = float(Lz.split()[1]) - float(Lz.split()[0])
# print(Lx,Ly,Lz,time,natoms)
# f.readline()
#  
# counter = 0
# for line in f:
# 	if counter < natoms:
# 		counter +=1
# 		line = line.strip()
# 		columns = line.split()
# 		id = int(columns[0])
# 		type = int(columns[1])
# 		posx = float(columns[2])
# 		posy = float(columns[3])
# 		posz = float(columns[4])
# 		qx = float(columns[5])
# 		qy = float(columns[6])
# 		qz = float(columns[7])
# 		qw = float(columns[8])
# 		aaxis = float(columns[9])
# 		baxis = float(columns[10])
# 		caxis = float(columns[11])
# 		print(id,type,posx,posy,posz,qx,qy,qz,aaxis,baxis,caxis)
# 
# print('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
#read the rest of the files
counter = 0
id = []
type = []
r = []
q = []
time = 0

for line in f:
	counter += 1
	if counter == 2:
		time = float(line.strip())
		print(time)
	if counter == 4:
		natoms = int(line.strip())
	if counter == 6:
		L = line.strip()
		L = float(L.split()[1])
	if counter > 9: #number of header lines
		line = line.strip()
		columns = line.split()
		id.append(int(columns[0]))
		type.append(int(columns[1]))
		r.append(np.array([float(columns[2]),float(columns[3]),float(columns[4])]))
		q.append(np.array([float(columns[5]),float(columns[6]),float(columns[7]),float(columns[8])]))
		aaxis = float(columns[9])
		baxis = float(columns[10])
		caxis = float(columns[11])
	if counter > 520 : #number of lines in one time step
		counter = 0
		
		
 		
k_space = np.linspace(start=0, stop=40, num=41, endpoint=True).tolist() # num+1 because endpoint counts as well...
k_data = [ [float(kx),float(ky),float(kz)] for kx in k_space for ky in k_space for kz in k_space]
k_data = np.delete(k_data, 0, 0)
sk_static = []
 
klist = np.zeros((int(np.sqrt(3.*len(k_space)**2)),2))
 
for k in k_data:
	normk = np.linalg.norm(k)
 	 
	if normk != 0:
 		 
		if normk % 1 == 0:
 			 
			bin = klist[int(normk)][0]
			rhoq2_per_time = klist[int(normk)][1]
 			
			bin += 1
 			
			cosq = 0
			sinq = 0
 			
			for rr in r:
				cosq  += np.cos( np.dot(k,rr)*2.*np.pi/L )
				sinq += np.cos( np.dot(k,rr)*2.*np.pi/L )
 			
# 			cosq =  np.sum( [np.cos( np.dot(k,rr)*2.*np.pi/L )  for rr in r ] ) 
# 			sinq =  np.sum( [np.sin( np.dot(k,rr)*2.*np.pi/L )  for rr in r ] ) 
 			 
			rhoq2_per_time += (cosq*cosq + sinq*sinq)/natoms 
 
# 			rhoq_per_time = (np.sum( [ np.exp(-1j * np.dot(k,rr) * 2*np.pi/L) for rr in r ] ) )
# 			rhoq2_per_time += np.real(rhoq_per_time*np.conj(rhoq_per_time)) / natoms
 			 
 			 
			klist[int(normk)][0] = bin
			klist[int(normk)][1] = rhoq2_per_time
 			 
 		 
sk_static = []
  
for (bins,rhos) in klist.tolist():
	if bins != 0:
		sk_static.append(rhos/bins)
	else:
		sk_static.append(rhos)
 
 
  
fdat = open('elli.dat','w')
 
for sks in sk_static:
 
	fdat.write( str(format(sks, '.6f')) + '\n')
		
	
	


# headertimestep = f.readline()
# time = f.readline()
# time = time.strip()
# headeratoms = f.readline()
# natoms = f.readline()
# natoms = int(natoms.strip())
# headerbox = f.readline()
# Lx = f.readline()
# Lx = Lx.strip()
# Lx = float(Lx.split()[1]) - float(Lx.split()[0])
# Ly = f.readline()
# Ly = Ly.strip()
# Ly= float(Ly.split()[1]) - float(Ly.split()[0])
# Lz = f.readline()
# Lz = Lz.strip()
# Lz = float(Lz.split()[1]) - float(Lz.split()[0])
# print(Lx,Ly,Lz,time,natoms)
# f.readline()
# 
# counter = 0
# if counter % natoms == 0:
# 	counter = 0
# for line in f:
# 	if counter < natoms:
# 		counter +=1
# 		line = line.strip()
# 		columns = line.split()
# 		id = int(columns[0])
# 		type = int(columns[1])
# 		posx = float(columns[2])
# 		posy = float(columns[3])
# 		posz = float(columns[4])
# 		qx = float(columns[5])
# 		qy = float(columns[6])
# 		qz = float(columns[7])
# 		qw = float(columns[8])
# 		aaxis = float(columns[9])
# 		baxis = float(columns[10])
# 		caxis = float(columns[11])
# 		print(id,type,posx,posy,posz,qx,qy,qz,aaxis,baxis,caxis)
