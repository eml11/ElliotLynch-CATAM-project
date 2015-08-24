import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

inner = np.loadtxt('./stellarstructure_inner.txt',skiprows=2)
outer = np.loadtxt('./stellarstructure_outer.txt',skiprows=2)[::-1]

data = np.concatenate((inner,outer),axis=0)

gamma = 5./3.

rsun = 6.9598e8
msun = 1.9891e30
gbar = 1.0e14

mass_max = data[-1,0]

print
print 'Stellar Parameters'
print
print 'Surface'
print
print 'Core'
print
print 'Error'
print 'Pressure: ', np.abs((inner[-1,1] - outer[0,1])/inner[-1,1])
print 'Radius: ',  np.abs((inner[-1,2] - outer[0,2])/inner[-1,2])

fig, ax1 = plt.subplots()

plot1 = ax1.plot(data[:,0],data[:,1],'k',linewidth=2.0,label='Radius [Solar Radii]')
ax1.set_ylabel('Radius',fontsize=14)
ax1.set_xlabel('Mass [Solar Masses]',fontsize=14)

ax2 = ax1.twinx()

plot2 = ax2.plot(data[:,0],data[:,2],'k:',linewidth=2.0,label='Pressure [Gbar]')
ax2.set_ylabel('Pressure',fontsize=14)

plots = plot1 + plot2
labls = [l.get_label() for l in plots]
ax1.legend(plots,labls,loc='best')

ax1.set_ylim(0,1.5)
ax2.set_ylim(0,15)
ax1.set_xlim(0,3.0)
ax2.set_xlim(0,3.0)

plt.show()

