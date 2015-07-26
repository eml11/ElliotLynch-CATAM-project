import numpy as np
import pylab as pb

inner = np.loadtxt('./stellarstructure_inner.txt',skiprows=2)
outer = np.loadtxt('./stellarstructure_outer.txt',skiprows=2)[::-1]

data = np.concatenate((inner,outer),axis=1)

gamma = 5./3.

mass_max = data[-1,0]
radiusvar_max = data[-1,1]
pressurevar_max = data[0,2]

Temperature = (data[:,2]/pressurevar_max)**(1 - 1/gamma)

print
print 'Stellar Parameters'
print
print 'Surface'
print 'Temperature: ', Temperature[-1]
print 'Radius: ', radiusvar_max
print
print 'Core'
print 'Temperature: ', Temperature[0]
print 'Pressure: ', pressurevar_max
print
print 'Error'
print 'Pressure: ', np.abs((inner[-1,1] - outer[0,1])/inner[-1,1])
print 'Radius: ',  np.abs((inner[-1,2] - outer[0,2])/inner[-1,2])



pb.plot(data[:,0]/mass_max,data[:,1]/radiusvar_max,linewidth=2.0,label='Normalised r^3')
pb.plot(data[:,0]/mass_max,data[:,2]/pressurevar_max,linewidth=2.0,label='Normalised lnP')
pb.plot(data[:,0]/mass_max,linewidth=2.0,label='Normalised T^4')

pb.xlabel('Normalised Mass')

pb.legend(loc='best')
pb.show()

