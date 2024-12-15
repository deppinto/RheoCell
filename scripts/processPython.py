import sys
import numpy as np
from math import *
import matplotlib.pyplot as plt


if len(sys.argv)!=3:
        print(sys.argv[0]," [scripts] [variable]")
        sys.exit(1)

scripts=int(float(sys.argv[1]))
variable=int(float(sys.argv[2]))

'''
#-------------------------------------optimization plot
plt.figure(figsize=(5.452423529,4.089317647))
fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/optimization.txt","r")
cores=[]
runtime=[]
for line in fileoutput:
    save=line.split()
    cores.append((float(save[0])))
    runtime.append((float(save[1])))
fileoutput.close()

plt.plot(cores, runtime, '--o')
#plt.yscale('log')
#plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.ylabel('Runtime(s)', fontsize=18)
plt.xlabel('# cores', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/runtime_cores.png", transparent=True)
plt.show()
exit(1)
#-------------------------------------optimization plot
'''

'''
filedata=open("/home/diogo/MEGA/cenas/SAT-assembly/PatchyParticles_fork/Results/finite_size_shells/scripts"+str(scripts)+"/dados.txt","r")
#filedata=open("/home/diogo/MEGA/cenas/SAT-assembly/PatchyParticles_fork/Results/scripts"+str(scripts)+"/dados_total.txt","r")
N=[]
rho=[]
angle=[]
temp=[]
cosmax=[]
jobs=[]
jobs_seq=[0]
for line in filedata:
	save=line.split()
	N.append(float(save[0]))
	rho.append(float(save[1]))
	angle.append(float(save[2]))
	temp.append(float(save[3]))
	cosmax.append(float(save[4]))
	jobs.append(float(save[5]))
	jobs_seq.append(float(save[5])+jobs_seq[len(jobs_seq)-1])
filedata.close()
'''

# read output file
plt.figure(figsize=(5.452423529,4.089317647))
avg_vel=[]
for job in range(1, 6, 2):
    fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/MSD.txt","r")
    time=[]
    MSD=[]

    for line in fileoutput:
        save=line.split()
        time.append((float(save[0])))
        MSD.append((float(save[1])+1))

    fileoutput.close()

    fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/v_width.txt","r")
    y=[]
    vel=[]

    for line in fileoutput:
        save=line.split()
        y.append((float(save[2])))
        vel.append((float(save[1])))

    fileoutput.close()

    fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/mean_velocity.txt","r")
    for line in fileoutput:
        save=line.split()
        avg_vel.append((float(save[0])))

    #plt.plot(time, MSD, 'o')
    #plt.plot(y, vel, 'o')

xvals=[0.01, 0.1, 0.4]
plt.plot(xvals, avg_vel, '--o')

'''
fitline2=[]
fitline1=[]
for i in range(len(time)):
    fitline2.append(time[i]*time[i]/3e4)
    fitline1.append(time[i]/4.5e1)
plt.plot(time, fitline2, '--', lw=2)
plt.plot(time, fitline1, '--', lw=2)
'''

#plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
#plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
plt.xscale('log')
#plt.yscale('log')
#plt.ylabel('MSD', fontsize=18)
#plt.xlabel('Time', fontsize=18)
#plt.xlim(1e3, 1e4)
#plt.ylabel('Velocity', fontsize=18)
#plt.xlabel('Width', fontsize=18)
plt.ylabel('Mean velocity', fontsize=18)
plt.xlabel(r'$\omega$', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#plt.legend([r'$\omega$=0.01',r'$\omega$=0.1',r'$\omega$=0.4'], fontsize=14, loc=(0.,0.65), frameon=False)
plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
#plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/MSD_2.png", transparent=True)
plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/mean_velocity_1.png", transparent=True)
plt.show()

exit(1)

'''
calc1=[i / j for i, j in zip(y1, jobs)]
calc2=[i / j for i, j in zip(y2, jobs)]
calc3=[i / j for i, j in zip(y3, jobs)]
calc4=[i / j for i, j in zip(y4, jobs)]
calc5=[i / j for i, j in zip(y5, jobs)]
#plt.plot(angle[0:11], calc[0:11], 'o-', color='darkred')
#plt.plot(angle[11:22], calc[11:22], 's-', color='darkgreen')
#plt.plot(angle[22:33], calc[22:33], '^-', color='darkblue')
#plt.legend(['cosmax=0.98', 'cosmax=0.95', 'cosmax=0.925'])
#plt.figure(figsize=(8,6))
plt.figure(figsize=(5.452423529,4.089317647))
if variable < 2:
	plt.plot(x1[:], y1[:], '^', color='darkred')
	plt.plot(x2[:], y2[:], 's', color='darkgreen')
	plt.plot(x3[:], y3[:], 'p', color='darkblue')
	if variable==0:
		plt.plot(x4[:], y4[:], 'x', color='black')
	plt.plot(x5[:], y5[:], '*', color='gray')
elif variable==2:
	#plt.plot(x1[0:30:3], calc1[0:30:3], '--^', color='darkred')
	#plt.plot(x2[1:30:3], calc2[1:30:3], '--s', color='darkgreen')
	#plt.plot(x3[2:30:3], calc3[2:30:3], '--p', color='darkblue')

	#first=3
	#second=6
	#print (calc1[first:second], calc2[first:second], calc3[first:second], jobs)
	#plt.plot(x1[first:second], calc1[first:second], '--^', color='darkred')
	#plt.plot(x2[first:second], calc2[first:second], '--s', color='darkgreen')
	#plt.plot(x3[first:second], calc3[first:second], '--p', color='darkblue')

        #plt.plot(x1[:], calc1[:], '--^', color='darkred')
        #plt.plot(x2[:], calc2[:], '--s', color='darkgreen')
        #plt.plot(x3[:], calc3[:], '--p', color='darkblue')

        plt.plot(x1[0:11], calc1[0:11], '--^', color='darkred')
        plt.plot(x1[22:33], calc1[22:33], '--s', color='darkgreen')
        plt.plot(x1[44:55], calc1[44:55], '--p', color='darkblue')
	
        #avg=[(calc1[0]+calc1[30])/2, calc1[3], (calc1[6]+calc1[32])/2, (calc1[9]), (calc1[12]+calc1[34])/2, (calc1[15]), (calc1[18]+calc1[36])/2, (calc1[21]), (calc1[24]+calc1[38])/2, (calc1[27]+calc1[40])/2]	
        #plt.plot(x1[0:30:3], avg[:], '--^', color='darkred')
        #plt.plot(x1[0:30:3], calc1[0:30:3], '--^', color='darkred')
        #plt.plot(x2[1:30:3], calc2[1:30:3], '--^', color='darkred')
        #plt.plot(x1[30:42:2], calc1[30:42:2], '--s', color='darkgreen')
        #plt.plot(x2[31:42:2], calc2[31:42:2], '--s', color='darkgreen')
        #plt.plot(x1[42:54:2], calc1[42:54:2], '--p', color='darkblue')
        #plt.plot(x2[43:54:2], calc2[43:54:2], '--p', color='darkblue')
        #plt.plot(x1[54:66:2], calc1[54:66:2], '--h', color='black')
        #plt.plot(x2[55:66:2], calc2[55:66:2], '--h', color='black')
        #plt.plot(x1[66:78:2], calc1[66:78:2], '--8', color='gray')
        #plt.plot(x2[67:78:2], calc2[67:78:2], '--8', color='gray')
        #plt.plot(x1[78:90:2], calc1[78:90:2], '--D', color='goldenrod')
        #plt.plot(x2[79:90:2], calc2[79:90:2], '--D', color='goldenrod')
        #plt.plot(x1[90:102:2], calc1[90:102:2], '--d', color='peru')
        #plt.plot(x2[91:102:2], calc2[91:102:2], '--d', color='peru')

        #plt.plot(x3[2:33:3], calc3[2:33:3], '--^', color='darkred')
        #plt.plot(x3[33:39], calc3[33:39], '--s', color='darkgreen')
        #plt.plot(x3[39:45], calc3[39:45], '--p', color='darkblue')
        #plt.plot(x3[45:51], calc3[45:51], '--h', color='black')
        #plt.plot(x3[51:57], calc3[51:57], '--8', color='gray')
        #plt.plot(x3[57:63], calc3[57:63], '--D', color='goldenrod')
        #plt.plot(x3[63:69], calc3[63:69], '--d', color='peru')

	#plt.plot(x4[0:11], calc4[0:11], '--x', color='black')
	#plt.plot(x5[0:11], calc5[0:11], '--*', color='gray')
	#plt.plot(x1[22:33], calc1[22:33], '--^', color='darkred')
	#plt.plot(x2[22:33], calc2[22:33], '--s', color='darkgreen')
	#plt.plot(x3[22:33], calc3[22:33], '--p', color='darkblue')
	#plt.plot(x4[22:33], calc4[22:33], '--x', color='black')
	#plt.plot(x5[22:33], calc5[22:33], '--*', color='gray')

#print (x2[1:30:3], calc2[1:30:3])

if variable==0:
	plt.legend(['icosahedron', 'snubcube', 'snubdodecahedron', 'incomplete', 'polymorphs'], fontsize=14, loc=(1.5,1.5))
elif variable==1:
	plt.legend(['icosahedron', 'snubcube', 'snubdodecahedron'], fontsize=14, loc='center')
elif variable==2:
	#plt.legend(['icosahedron', 'snubcube', 'snubdodecahedron', 'incomplete', 'polymorphs'], fontsize=14, loc=(0.65, 0.05), frameon=False)
	plt.legend(['0.98', '0.95', '0.925'], fontsize=14, loc=(0.05, 0.25), frameon=False)
        #plt.legend(['icosahedron, C2(1)', 'snubcube, C2(1)', 'icosahedron, C2(2)', 'snubcube, C2(2)', 'icosahedron, C3(1)', 'snubcube, C3(1)', 'icosahedron, C3(2)', 'snubcube, C3(2)', 'icosahedron, C4(1)', 'snubcube, C4(1)', 'icosahedron, C4(2)', 'snubcube, C4(2)', 'icosahedron, C5', 'snubcube, C5'], fontsize=14, loc=(0.7, 0.05))
	#plt.legend(['icosahedron, C2(1)', 'icosahedron, C2(2)', 'icosahedron, C3(1)', 'icosahedron, C3(2)', 'icosahedron, C4(1)', 'icosahedron, C4(2)', 'icosahedron, C5'], fontsize=14, loc=(0.7, 0.05))
	#plt.legend(['snubcube, C2(1)', 'snubcube, C2(2)', 'snubcube, C3(1)', 'snubcube, C3(2)', 'snubcube, C4(1)', 'snubcube, C4(2)', 'snubcube, C5'], fontsize=14, loc=(0.7, 0.05))
	#plt.legend(['C2(1)', 'C2(2)', 'C3(1)', 'C3(2)', 'C4(1)', 'C4(2)', 'C5'], fontsize=14, loc=(0.8, 0.5), frameon=False)
	#plt.legend(['C2(1), C2(2)', 'C3(1)', 'C3(2)', 'C4(1)', 'C4(2)', 'C5'], fontsize=14, loc=(0.705, 0.5), frameon=False)
if variable==1:
	plt.xlabel('T', fontsize=18, fontname='Times New Roman')
elif variable==0:
	plt.xlabel(r'$\gamma$', fontsize=18, fontname='Times New Roman')
elif variable==2:
	#plt.xlabel(r'$\rho$', fontsize=18, fontname='Times New Roman')
	plt.xlabel(r'$\gamma$', fontsize=18, fontname='Times New Roman')
	#plt.xlabel('In plane angle', fontsize=18, fontname='Times New Roman')
if variable==0:
	plt.ylabel(r'$cos\theta_{max}$', fontsize=18, fontname='Times New Roman')
elif variable==1:
	plt.ylabel(r'$\rho$', fontsize=18, fontname='Times New Roman')
elif variable==2:
	plt.ylabel('Yield', fontsize=18, fontname='Times New Roman')

#plt.ylim([-0.05,1.05])
#plt.xlim([0.00075,5])
#plt.xlim([0,0.65])
#plt.xscale('log')
#plt.title('T=0.1', fontsize=18, fontname='Times New Roman')
plt.xticks(fontsize=18, fontname='Times New Roman')
plt.yticks(fontsize=18, fontname='Times New Roman')
plt.subplots_adjust(left=0.175, bottom=0.175, right=0.95, top=0.95)

plt.show()
'''
