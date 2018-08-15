import numpy as np
f = open('output/elli3.dump','r')
print("this is for 512 particle system")

counter = 0
time = 0

natoms = 512

k_space = np.linspace(start=0, stop=40, num=41, endpoint=True).tolist() # num+1 because endpoint counts as well...
k_data = [ [float(kx),float(ky),float(kz)] for kx in k_space for ky in k_space for kz in k_space]
k_data = np.delete(k_data, 0, 0)

sk_ave = []

id = []
type = []
r = []
q = []
sk_static = []
sk_ave = []

for line in f:
    counter += 1
    if counter == 2:
        time = float(line.strip())
        print("time",time,counter)
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
                         
                    rhoq2_per_time += (cosq*cosq + sinq*sinq)/natoms
            
                    klist[int(normk)][0] = bin
                    klist[int(normk)][1] = rhoq2_per_time
        
        count = 0
        for (bins,rhos) in klist.tolist():
            count += 1
            
            if bins != 0:
                sk_static.append(rhos/bins)
            else:
                sk_static.append(rhos)
     
        counter = 0
        id = []
        type = []
        r = []
        q = []
     
    if sk_static != []:
        sk_ave.append(sk_static)
       
if sk_ave != []:
    sk_ave = np.array(sk_ave)
   
    sks = np.mean(sk_ave, axis=0)
 
    fdat = open('elli3.dat','w')
    
    for skss in sks:
    
        fdat.write( str(format(skss, '.6f')) + '\n')
