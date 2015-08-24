import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import sys

sysargs = sys.argv

if len(sysargs) > 1:
  if sysargs[1] == 'in':
    outer = np.loadtxt('./stellarstructure_outer.txt',skiprows=2)[::-1]
    data = outer
  elif sysargs[1] == 'out':
    inner = np.loadtxt('./stellarstructure_inner.txt',skiprows=2) 
    data = inner   
  else:
    inner = np.loadtxt('./stellarstructure_inner.txt',skiprows=2)
    outer = np.loadtxt('./stellarstructure_outer.txt',skiprows=2)[::-1]

    data = np.concatenate((inner,outer),axis=0)
else:
  inner = np.loadtxt('./stellarstructure_inner.txt',skiprows=2)
  outer = np.loadtxt('./stellarstructure_outer.txt',skiprows=2)[::-1]

  data = np.concatenate((inner,outer),axis=0)

gamma = 5./3.

rsun = 6.9598e8
msun = 1.9891e30
gbar = 1.0e14
mkelvin = 1.0e6
lsun=3.8515e26

print
print
print
if len(sysargs) == 1:
  print 'Error'
  print 'Radius: ', np.abs((inner[-1,1] - outer[0,1])/inner[-1,1])
  print 'Pressure: ',  np.abs((inner[-1,2] - outer[0,2])/inner[-1,2])
  print 'Temperature: ',  np.abs((inner[-1,3] - outer[0,3])/inner[-1,3])
  print 'Luminosity: ',  np.abs((inner[-1,4] - outer[0,4])/inner[-1,4])

fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax3 = ax1.twinx()

axes = [ax1,ax2,ax3]

fig.subplots_adjust(right=0.75)

axes[-1].spines['right'].set_position(('axes', 1.2))

axes[-1].set_frame_on(True)
axes[-1].patch.set_visible(False)

plot1 = ax1.plot(data[:,0],data[:,1],'k',linewidth=2.0,label='Radius [Solar Radii]')
plot2 = ax2.plot(data[:,0],data[:,4],'k--',linewidth=2.0,label='Luminosity [Solar Luminosity]')
ax1.set_ylabel('Radius',fontsize=14)
ax1.set_xlabel('Mass [Solar Masses]',fontsize=14)

ax2.set_ylabel('Luminosity',fontsize=14)

plot3 = ax3.plot(data[:,0],data[:,2],'k:',linewidth=2.0,label='Pressure [Gbar]')
plot4 = ax3.plot(data[:,0],data[:,3],'k-.',linewidth=2.0,label='Temperature [MK]')
ax3.set_ylabel('Pressure,Temperature',fontsize=14)

plots = plot1 + plot2 + plot3 + plot4
labls = [l.get_label() for l in plots]
ax3.legend(plots,labls,loc='best')

ax1.grid()
ax1.minorticks_on()
ax1.set_ylim(0,1.5)
ax2.set_ylim(0,220)
ax3.set_ylim(0,65)
ax1.set_xlim(0,3.0)

plt.show()
