import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import sys

inner = np.loadtxt(sys.argv[1] + '_inner.txt',skiprows=2)
outer = np.loadtxt(sys.argv[1] + '_outer.txt',skiprows=2)[::-1]

data = np.concatenate((inner,outer),axis=0)

inner = np.loadtxt(sys.argv[2] + '_inner.txt',skiprows=2)
outer = np.loadtxt(sys.argv[2] + '_outer.txt',skiprows=2)[::-1]

dataf = np.concatenate((inner,outer),axis=0)

inner = np.loadtxt(sys.argv[3] + '_inner.txt',skiprows=2)
outer = np.loadtxt(sys.argv[3] + '_outer.txt',skiprows=2)[::-1]

datab = np.concatenate((inner,outer),axis=0)

#data = outer

gamma = 5./3.

rsun = 6.9598e8
msun = 1.9891e30
gbar = 1.0e14

#put in temperature and luminosity later
mass_max = data[-1,0]
#radiusvar_max = data[0,1]
#pressurevar_max = data[0,2]
#tmax = data[0,3]
#lummax =  data[-1,4]

#Temperature = data[:,3]

print
print 'Stellar Parameters'
print
print 'Surface'
#print 'Temperature: ', Temperature[-1]
#print 'Radius: ', radiusvar_max
print
print 'Core'
#print 'Temperature: ', Temperature[0]
#print 'Pressure: ', pressurevar_max
print
print 'Error'
print 'Pressure: ', np.abs((inner[-1,1] - outer[0,1])/inner[-1,1])
print 'Radius: ',  np.abs((inner[-1,2] - outer[0,2])/inner[-1,2])

#data[:,1] = np.where(data[:,1]>1.0,1.0,data[:,1])
#data[:,2] = np.where(data[:,2]>1.0,1.0,data[:,2])

fig, ax1 = plt.subplots()

plot1 = ax1.plot(data[:,0],data[:,1],'k',linewidth=2.0,label='Radius [Solar Radii]')
ax1.plot(dataf[:,0],dataf[:,1],'b',linewidth=2.0)
ax1.plot(datab[:,0],datab[:,1],'r',linewidth=2.0)

ax1.set_ylabel('Radius',fontsize=14)
ax1.set_xlabel('Mass [Solar Masses]',fontsize=14)

ax2 = ax1.twinx()

plot2 = ax2.plot(data[:,0],data[:,2],'k:',linewidth=2.0,label='Pressure [Gbar]')
ax2.plot(dataf[:,0],dataf[:,2],'b:',linewidth=2.0)
ax2.plot(datab[:,0],datab[:,2],'r:',linewidth=2.0)
ax2.set_ylabel('Pressure',fontsize=14)

plots = plot1 + plot2
labls = [l.get_label() for l in plots]
ax1.legend(plots,labls,loc='best')

ax1.set_ylim(0,1.5)
ax2.set_ylim(0,15)
ax1.set_xlim(0,3.0)
ax2.set_xlim(0,3.0)

#pb.plot(data[:,0]/msun,data[:,1]/rsun,'k',linewidth=2.0,label='Radius [Solar Radii]')
#pb.plot(data[:,0]/msun,data[:,2]/gbar,'k:',linewidth=2.0,label='Pressure [Gbar]')
#pb.plot(data[:,0]/mass_max,data[:,3]/tmax,'k--',linewidth=2.0,label='Normalised T')
#pb.plot(data[:,0]/mass_max,data[:,4]/lummax,'k-.',linewidth=2.0,label='Normalised L')
#pb.plot(data[:,0]/mass_max,linewidth=2.0,label='Normalised T^4')

#pb.xlabel('Mass [Solar Masses]')
#pb.ylim((0.0,1.0))

#pb.legend(loc='best')
plt.show()

