import numpy as np
from collections import deque

f = open('512part.dump','r')
print("this is for 512 particle system")

counter = 0
time = 0

natoms = 512

#1 Create k-space grid
k_space = np.linspace(start=0, stop=40, num=41, endpoint=True).tolist() # num+1 because endpoint counts as well...
k_data = [ [float(kx),float(ky),float(kz)] for kx in k_space for ky in k_space for kz in k_space]
k_data = np.delete(k_data, 0, 0)

#2 Sort the wave vectors according to increasing magnitude for averaging
seen = set() # use of set will make sure no element is repeated
seenlist = []
 
for (kx,ky,kz) in k_data:
    normk = np.linalg.norm([kx,ky,kz])
    if normk not in seen:
        seen.add(normk)
        seenlist.append(normk)
        
if len(seen) != len(seenlist):
    sys.exit("list error")
        
seenlist = sorted([normk for normk in seenlist]) 
klist = ([[normk] for normk in seenlist]) # convert each element to a list so it can be appended later

for (kx,ky,kz) in k_data:
    normk = np.linalg.norm([kx,ky,kz])
    ind = seenlist.index(normk)
    klist[ind].append([kx,ky,kz])
    

klist = deque(klist)     # list of (normk, [wavevector1],[wavevector2]...)
 
sk_ave = []
 
id = []
type = []
r = []
q = []
sk_static = []
sk_ave = []
mags =[]
 
for line in f:
    counter += 1
    if counter == 2:
        time = float(line.strip())
        print("time",time)
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
    if counter > natoms + 8 : #number of lines in one time step
        # start computing stuff
        sk_static = []     
        ## klist = np.zeros((int(np.sqrt(3.*len(k_space)**2)),2))
 
        sk = []
  
        for sets in klist:
            sets = deque(sets) # so you can popleft
            mag_k = sets.popleft() # normalized q, ie distance from the origin
            rhoq2_per_time = 0.
            knum = len(sets) # number of wave vectors with the same norm 
            while sets: # for all the wave vectors with the same norm
                k = sets.popleft()
#                 rhoq_per_time = (np.sum( [ exp(-1j * np.dot(q,r) * 2*np.pi/L) for r in rmin ] ) )
                cosq = np.sum( np.cos(np.dot(k,rr) * 2.*np.pi/L) for rr in r )
                sinq = np.sum( np.cos(np.dot(k,rr) * 2.*np.pi/L) for rr in r )
                rhoq_per_time = cosq*cosq + sinq*sinq
                rhoq2_per_time += rhoq_per_time/natoms
#                 rhoq2_per_time += np.real(rhoq_per_time*np.conj(rhoq_per_time))/len(d0)             
            sk_static.append(rhoq2_per_time/knum)
        
        #reset everything to zero for the next time frame calculation
        counter = 0
        id = []
        type = []
        r = []
        q = []
            
        sk_ave.append(sk_static) # list of sk per time frame

if sk_ave != []:
    sk_ave = np.array(sk_ave)
    sks = np.mean(sk_ave, axis=0)

    sks = list(sks)
        
fdat = open('512parttime.dat','w')
 
for (mag,sk) in zip(seenlist,sks):
        
    fdat.write(str(format(mag, '.6f')) + '\t'+ str(format(sk, '.6f')) + '\n')