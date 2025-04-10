import sys
import numpy as np
from math import *
import matplotlib.pyplot as plt
import os.path


if len(sys.argv)!=5:
        print(sys.argv[0]," [scripts] [start] [end] [variable]")
        sys.exit(1)

scripts=int(float(sys.argv[1]))
start=int(float(sys.argv[2]))
end=int(float(sys.argv[3]))
variable=int(float(sys.argv[4]))

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

#filedata=open("/scratch/pinto/Phase_Field/RheoCell/Work/Results/scripts"+str(scripts)+"/dados.txt","r")
filedata=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/params","r")
#JOBS	N	LX	LY	EQTIME	TIME	DT	PRTCONF	JO	FRICELL	FRIC	NEMSELF	NEMINTR	GAMMA	KAPPA	LAMBDA	MU	OMEGA	CGTOL	CORES	WSLIP	SHEAR
jobs=[]
jobs_seq=[1]
N=[]
lx=[]
ly=[]
eqtime=[]
time=[]
dt=[]
deltat=[]
J0=[]
friccell=[]
fric=[]
nemself=[]
nemintr=[]
gamma=[]
kappa=[]
llambda=[]
mu=[]
omega=[]
tolerance=[]
cores=[]
wallslip=[]
shear_rate=[]

F = []

cont_line=0
for line in filedata:
    save=line.split()
    if cont_line == end:
        break

    jobs.append(int(float(save[0])))
    jobs_seq.append(int(float(save[0]))+jobs_seq[len(jobs_seq)-1])

    N.append(int(float(save[1])))
    lx.append(int(float(save[2])))
    ly.append(int(float(save[3])))
    eqtime.append(float(save[4]))
    time.append(float(save[5]))
    dt.append(float(save[6]))
    deltat.append(float(save[7]))
    J0.append(float(save[8]))
    friccell.append(float(save[9]))
    fric.append(float(save[10]))
    nemself.append(float(save[11]))
    nemintr.append(float(save[12]))
    gamma.append(float(save[13]))
    kappa.append(float(save[14]))
    llambda.append(float(save[15]))
    mu.append(float(save[16]))
    omega.append(float(save[17]))
    tolerance.append(float(save[18]))
    cores.append(int(float(save[19])))
    wallslip.append(float(save[20]))
    shear_rate.append(float(save[21]))

    F.append( sqrt( float(save[10])/( float(save[11]) )) ) 

    cont_line+=1

filedata.close()


# read output file
n_holes = np.zeros(end-start)
lifetime_holes = np.zeros(end-start)
max_area_holes = np.zeros(end-start)

#voids_histogram_area.txt  voids_histogram_circularity.txt  voids_histogram_radius_speed_time.txt  voids_histogram_radius_speed.txt  voids_histogram_radius.txt  voids_stats.txt

fig_histograms , ax = plt.subplots(2, 3, figsize=(10,4), constrained_layout=True)
#fig_histograms.tight_layout()

for traj in range(start, end):

    area_histogram = np.zeros(int(0.5 * time[traj] / deltat[traj]))
    circularity_histogram = np.zeros(int(0.5 * time[traj] / deltat[traj]))
    radius_speed_histogram = np.zeros(int(0.5 * time[traj] / deltat[traj]))
    ani_histogram = np.zeros(int(0.5 * time[traj] / deltat[traj]))
    rspeed_histogram = np.zeros(int(lx[traj]))

    avg_stress_histogram = np.zeros(10)
    max_stress_histogram = np.zeros(10)
    #with open('voids_stress_histogram_tau10.txt', 'w') as f:

    l_bins = 0.1
    x_l = [i * l_bins for i in range(int(1./l_bins))]
    lifetime_hist = np.zeros(int(1. / l_bins))
    a_bins = 0.02
    x_a = [i * a_bins for i in range(int(1./a_bins))]
    max_area_hist = np.zeros(int(1. / a_bins))
    total_hist_counts = 0

    count_jobs = 0
    for job in range(jobs_seq[traj], jobs_seq[traj+1]):

        if job==144 and scripts==74:
            continue
        if job==125 and scripts==75:
            continue

        fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_stats.txt","r")
        nn = 0
        for line in fileoutput:
            save=line.split()
            nn = int(float(save[2]))
            n_holes[traj - start] += nn / jobs[traj]
            hlife = (float(save[3]) / (time[traj] * dt[traj]))
            maxA = float(save[4]) * pi * 8 * 8 / (lx[traj] * ly[traj])
            lifetime_holes[traj - start] += hlife / jobs[traj]
            max_area_holes[traj - start] += maxA / jobs[traj]

            for i in range(nn):
                hlife = (float(save[6+i]) / (time[traj] * dt[traj])) / l_bins
                maxA = float(save[6+nn+i+1]) * pi * 8 * 8 / (lx[traj] * ly[traj] * a_bins)

                lifetime_hist[int(hlife)] += 1
                max_area_hist[int(maxA)] += 1
                total_hist_counts += 1

        fileoutput.close()
        if nn > 0:
            count_jobs += 1


        fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_histogram_area.txt","r")
        for line in fileoutput:
            save=line.split()
            area_histogram[int(float(save[0]) * 50.01)] += (float(save[1]) * pi * 8 * 8) / (lx[traj] * ly[traj])
        fileoutput.close()

        fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_histogram_circularity.txt","r")
        for line in fileoutput:
            save=line.split()
            circularity_histogram[int(float(save[0]) * 50.01)] += float(save[1])
        fileoutput.close()

        fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_histogram_radius_speed_time.txt","r")
        for line in fileoutput:
            save=line.split()
            radius_speed_histogram[int(float(save[0]) * 50.01)] += float(save[1])
        fileoutput.close()

        fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_histogram_ani.txt","r")
        for line in fileoutput:
            save=line.split()
            ani_histogram[int(float(save[0]) * 50.01)] += float(save[1])
        fileoutput.close()

        fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_histogram_radius_speed.txt","r")
        for line in fileoutput:
            save=line.split()
            rspeed_histogram[int(float(save[0]))] += float(save[1])
        fileoutput.close()


        fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_stress_histogram_tau10.txt","r")
        for line in fileoutput:
            save=line.split()
            avg_stress_histogram[int(float(save[0]))] += float(save[1])
            max_stress_histogram[int(float(save[0]))] += float(save[2])
        fileoutput.close()


    if count_jobs > 0:
        area_histogram[:] = area_histogram[:] / count_jobs
        circularity_histogram[:] = circularity_histogram[:] / count_jobs
        radius_speed_histogram[:] = radius_speed_histogram[:] / count_jobs
        rspeed_histogram[:] = rspeed_histogram[:] / count_jobs
        ani_histogram[:] = ani_histogram[:] / count_jobs

        avg_stress_histogram[:] = avg_stress_histogram[:] / count_jobs
        max_stress_histogram[:] = max_stress_histogram[:] / count_jobs

    if variable == 1:
        if total_hist_counts>0:
            lifetime_hist[:] = lifetime_hist[:] / total_hist_counts
            max_area_hist[:] = max_area_hist[:] / total_hist_counts
        ax[1,0].plot(x_l, lifetime_hist, '-o')
        ax[1,1].plot(x_a, max_area_hist, '-s', label=nemself[traj])

    if variable == 2:
        ax[0,0].plot(area_histogram, '-o')
        ax[0,1].plot(circularity_histogram, '-s')
        ax[1,0].plot(radius_speed_histogram, '-^')
        ax[1,1].plot(rspeed_histogram, '-v')
        ax[1,2].plot(ani_histogram, '-D', label=fric[traj])

    if variable == 3:
        ax[0,0].plot(avg_stress_histogram, '-o')
        ax[0,1].plot(max_stress_histogram, '-s', label=omega[traj])


if variable == 2:
    ax[0,0].set_ylabel('Area')
    ax[0,1].set_ylabel('Circularity')
    ax[1,0].set_ylabel('Radius speed')
    ax[1,1].set_ylabel('Radius speed')
    ax[1,2].set_ylabel('Anisotropy ratio')
    ax[0,0].set_xlabel('Time')
    ax[0,1].set_xlabel('Time')
    ax[1,0].set_xlabel('Time')
    ax[1,1].set_xlabel('Radius')
    ax[1,2].set_xlabel('Time')

    fig_histograms.delaxes(ax[0,2])
    fig_histograms.legend(loc=(0.7, 0.6), ncols=2, frameon='false')
    #plt.legend(loc="upper right", frameon='false')
    plt.show()
    exit (1)


if variable==1:
    ax[0,0].set_ylabel('Lifetime')
    ax[0,1].set_ylabel('Max area')
    ax[1,0].set_ylabel('P(Lifetime)')
    ax[1,1].set_ylabel('P(Max area)')
    ax[0,0].set_xlabel('X')
    ax[0,1].set_xlabel('X')
    ax[1,0].set_xlabel('Normalized lifetime')
    ax[1,1].set_xlabel('Normalized max area')
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox='false', ncol=5)
    #fig.legend(loc=7)
    #ax[0,0].plot(nemself[0:5], lifetime_holes[0:5], '--p')
    #ax[0,0].plot(mu[0:4], lifetime_holes[0:4], '--p')
    #ax[0,0].plot(gamma[4:8], lifetime_holes[0:4], '--p')
    #ax[0,0].plot(F[0:5], lifetime_holes[0:5], '--p')
    #ax[0,0].plot(F[5:10], lifetime_holes[5:10], '--p')
    #ax[0,0].plot(F[10:15], lifetime_holes[10:15], '--p')

    #ax[0,1].plot(nemself[0:5], max_area_holes[0:5], '--v')
    #ax[0,1].plot(mu[0:4], max_area_holes[0:4], '--v')
    #ax[0,1].plot(gamma[4:8], max_area_holes[0:4], '--v')
    #ax[0,1].plot(F[0:5], max_area_holes[0:5], '--v')
    #ax[0,1].plot(F[5:10], max_area_holes[5:10], '--v')
    #ax[0,1].plot(F[10:15], max_area_holes[10:15], '--v')
    ax[0,0].plot(fric, lifetime_holes, '--p')
    ax[0,1].plot(fric, max_area_holes, '--v')
    #ax[0,0].set_xscale('log')
    #ax[0,1].set_xscale('log')

    #ax[0,2].plot(nemself[0:5], n_holes[0:5], '--^')
    #ax[0,2].plot(mu[0:4], n_holes[0:4], '--^')
    #ax[0,2].plot(gamma[4:8], n_holes[0:4], '--^')
    #ax[0,2].plot(F[0:5], n_holes[0:5], '--^')
    #ax[0,2].plot(F[5:10], n_holes[5:10], '--^')
    #ax[0,2].plot(F[10:15], n_holes[10:15], '--^')
    ax[0,2].plot(fric, n_holes, '--^')
    #ax[0,2].set_xscale('log')
    #ax[0,2].set_ylabel('# voids')
    #ax[0,2].set_xlabel('X')

    fig_histograms.delaxes(ax[1,2])
    fig_histograms.legend(loc=(0.75,0.15), ncols=2, frameon='false')

    #fig1 = plt.figure(figsize=(5.452423529,4.089317647))
    #plt.plot(nemself[0:5], n_holes[0:5], '--o')
    #plt.plot(nemself[5:10], n_holes[0:5], '--o')
    #plt.plot(fric, n_holes, '--^')
    #plt.xscale('log')
    #plt.ylabel('# voids', fontsize=18)
    #plt.xlabel('P', fontsize=18)
    plt.show()
    exit(1)

if variable == 3:
    ax[0,0].set_ylabel('Avg stress')
    ax[0,1].set_ylabel('Max stress')
    ax[0,0].set_xlabel('Time')
    ax[0,1].set_xlabel('Time')
    fig_histograms.delaxes(ax[0,2])
    fig_histograms.legend(loc=(0.7, 0.6), ncols=2, frameon='false')
    plt.show()
    exit(1)


#first figure-----------------------------------------------------------------------------------------------
fig1 = plt.figure(figsize=(5.452423529,4.089317647))
plt.plot(nemself[0:5], n_holes[0:5], '--o')
#plt.plot(nemself[5:10], n_holes[5:10], '--s')
#plt.plot(nemself[10:15], n_holes[10:15], '--^')
#plt.legend([r'$\xi$=0.1',r'$\xi$=0.01',r'$\xi$=0.001'], fontsize=14, frameon=False)

#plt.plot(fric, n_holes, '--^')
#plt.xscale('log')
#plt.yscale('log')

#plt.plot(omega[0:7], n_holes[0:7], '--o')
#plt.plot(omega[7:14], n_holes[7:14], '--s')
#plt.legend([r'$\xi$=0.01',r'$\xi$=0.001'], fontsize=14, frameon=False)

plt.ylabel('# voids', fontsize=18)
plt.xlabel('P', fontsize=18)




#second figure----------------------------------------------------------------------------------------------
#fig2 = plt.figure(figsize=(5.452423529,4.089317647))
#plt.plot(nemself[0:5], lifetime_holes[0:5], '--o')
#plt.plot(nemself[5:10], lifetime_holes[5:10], '--s')
#plt.plot(nemself[10:15], lifetime_holes[10:15], '--^')
#plt.legend([r'$\xi$=0.1',r'$\xi$=0.01',r'$\xi$=0.001'], fontsize=14, frameon=False)

#plt.plot(fric, lifetime_holes, '--^')
#plt.xscale('log')
#plt.yscale('log')

#plt.plot(omega[0:7], lifetime_holes[0:7], '--o')
#plt.plot(omega[7:14], lifetime_holes[7:14], '--s')
#plt.legend([r'$\xi$=0.01',r'$\xi$=0.001'], fontsize=14, frameon=False)

#plt.ylabel('Lifetime', fontsize=18)
#plt.xlabel('P', fontsize=18)


#third figure-----------------------------------------------------------------------------------------------
#fig3 = plt.figure(figsize=(5.452423529,4.089317647))
#plt.plot(nemself[0:5], max_area_holes[0:5], '--o')
#plt.plot(nemself[5:10], max_area_holes[5:10], '--s')
#plt.plot(nemself[10:15], max_area_holes[10:15], '--^')
#plt.legend([r'$\xi$=0.1',r'$\xi$=0.01',r'$\xi$=0.001'], fontsize=14, frameon=False)

#plt.plot(fric, max_area_holes, '--^')
#plt.xscale('log')
#plt.yscale('log')

#plt.plot(omega[0:7], max_area_holes[0:7], '--o')
#plt.plot(omega[7:14], max_area_holes[7:14], '--s')
#plt.legend([r'$\xi$=0.01',r'$\xi$=0.001'], fontsize=14, frameon=False)

#plt.ylabel('Max Area', fontsize=18)
#plt.xlabel('P', fontsize=18)


#fourth figure----------------------------------------------------------------------------------------------
#fig4 = plt.figure(figsize=(5.452423529,4.089317647))
#plt.plot(F[0:5], max_area_holes[0:5], '--o')
#plt.plot(F[5:10], max_area_holes[5:10], '--s')
#plt.plot(F[10:15], max_area_holes[10:15], '--^')
#plt.legend([r'$\xi$=0.1',r'$\xi$=0.01',r'$\xi$=0.001'], fontsize=14, frameon=False)

#plt.plot(fric, max_area_holes, '--^')
#plt.xscale('log')
#plt.yscale('log')

#plt.ylabel('Max Area', fontsize=18)
#plt.xlabel('F', fontsize=18)

plt.show()


exit(1)


'''
#stress_time plot
plt.figure(figsize=(5.452423529,4.089317647))
avg_vel=[]
for job in range(1, 6, 1):
    #fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/stress_time.txt","r")
    if job==4:
        fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts50/Job_4/stress_time.txt","r")
    else:
        fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/stress_time.txt","r")
    time=[]
    stress=[]
    cont=0

    for line in fileoutput:
        if cont%3!=0:
            cont+=1
            continue
        save=line.split()
        time.append((float(save[0])) * 0.1 * 500)
        stress.append((float(save[1])))
        cont+=1

    plt.plot(time, stress, '-o')


plt.show()
'''

#correlations below
filedata=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/params","r")
jobs = []
jobs_seq = [0]
lx = []
N = []
for line in filedata:
    save=line.split()
    if len(save) == 0:
        break
    jobs.append(float(save[0]))
    jobs_seq.append(float(save[0])+jobs_seq[len(jobs_seq)-1])
    N.append(float(save[1]))
    lx.append(float(save[2]))
filedata.close()

plt.figure(figsize=(5.452423529,4.089317647))
total_params = len(jobs)

if start==0:
    i=start
    size_R = int(lx[i]/2)
    corr_R = [i for i in range(size_R)]
    corr_vel = [0. for i in range(size_R)]
    avg_value = 0

    for j in range(0 + 1, 3 + 1):
        #fname = "/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts72/Job_"+str(j)+"/velocity_correlations.txt"
        fname = "/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts72/Job_"+str(j)+"/nematic_correlations.txt"
        if os.path.isfile(fname) == False:
            continue

        fileoutput=open(fname,"r")
        #print(fname)

        avg_value += 1
        for line in fileoutput:
            save=line.split()
            index = int(float(save[0]))
            value = float(save[1])
            if index<size_R:
                corr_vel[index] += value
        fileoutput.close()

    newList = [x/avg_value for x in corr_vel]
    plt.plot(corr_R, newList, "-o")

if start==8:
    i=start
    size_R = int(lx[i]/2)
    corr_R = [i for i in range(size_R)]
    corr_vel = [0. for i in range(size_R)]
    avg_value = 0

    for j in range(3 + 1, 6 + 1):
        #fname = "/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts72/Job_"+str(j)+"/velocity_correlations.txt"
        fname = "/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts72/Job_"+str(j)+"/nematic_correlations.txt"
        if os.path.isfile(fname) == False:
            continue

        fileoutput=open(fname,"r")
        #print(fname)

        avg_value += 1
        for line in fileoutput:
            save=line.split()
            index = int(float(save[0]))
            value = float(save[1])
            if index<size_R:
                corr_vel[index] += value
        fileoutput.close()

    newList = [x/avg_value for x in corr_vel]
    plt.plot(corr_R, newList, "-o")



for i in range(start, end):
    size_R = int(lx[i]/2)
    corr_R = [i for i in range(size_R)]
    corr_vel = [0. for i in range(size_R)]
    avg_value = 0

    for j in range(int(jobs_seq[i]) + 1, int(jobs_seq[i+1]) + 1):
        #fname = "/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(j)+"/velocity_correlations.txt"
        fname = "/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(j)+"/nematic_correlations.txt"
        if os.path.isfile(fname) == False:
            continue

        fileoutput=open(fname,"r")
        #print(fname)

        avg_value += 1
        for line in fileoutput:
            save=line.split()
            index = int(float(save[0]))
            value = float(save[1])
            if index<size_R:
                corr_vel[index] += value
        fileoutput.close()

    newList = [x/avg_value for x in corr_vel]
    plt.plot(corr_R, newList, "-o")

    if i==1:
        size_R = int(lx[i]/2)
        corr_R = [i for i in range(size_R)]
        corr_vel = [0. for i in range(size_R)]
        avg_value = 0
        for j in range(0 + 1, 3 + 1):
            #fname = "/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts71/Job_"+str(j)+"/velocity_correlations.txt"
            fname = "/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts71/Job_"+str(j)+"/nematic_correlations.txt"
            if os.path.isfile(fname) == False:
                continue

            fileoutput=open(fname,"r")
            #print(fname)

            avg_value += 1
            for line in fileoutput:
                save=line.split()
                index = int(float(save[0]))
                value = float(save[1])
                if index<size_R:
                    corr_vel[index] += value
            fileoutput.close()

        newList = [x/avg_value for x in corr_vel]
        plt.plot(corr_R, newList, "-o")

    if i==9:
        size_R = int(lx[i]/2)
        corr_R = [i for i in range(size_R)]
        corr_vel = [0. for i in range(size_R)]
        avg_value = 0
        for j in range(3 + 1, 6 + 1):
            #fname = "/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts71/Job_"+str(j)+"/velocity_correlations.txt"
            fname = "/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts71/Job_"+str(j)+"/nematic_correlations.txt"
            if os.path.isfile(fname) == False:
                continue

            fileoutput=open(fname,"r")
            #print(fname)

            avg_value += 1
            for line in fileoutput:
                save=line.split()
                index = int(float(save[0]))
                value = float(save[1])
                if index<size_R:
                    corr_vel[index] += value
            fileoutput.close()

        newList = [x/avg_value for x in corr_vel]
        plt.plot(corr_R, newList, "-o")

#plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
#plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
#plt.xscale('log')
#plt.yscale('log')
#plt.ylabel('MSD', fontsize=18)
#plt.xlabel('Time', fontsize=18)
#plt.xlim(1e3, 1e4)
#plt.ylabel('Velocity', fontsize=18)
#plt.xlabel('Width', fontsize=18)
#plt.ylabel(r'$C_v$', fontsize=18)
plt.ylabel(r'$C_Q$', fontsize=18)
plt.xlabel('R', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#plt.legend([r'$\omega$=0.01',r'$\omega$=0.1',r'$\omega$=0.4'], fontsize=14, loc=(0.,0.65), frameon=False)
plt.legend([r'$\xi_{cell}$=0', r'$\xi$=1',r'$\xi$=0.1',r'$\xi$=0.01', r'$\xi$=0.001',r'$\xi$=0.0001'], fontsize=14, loc=(0.575,0.32), frameon=False)
plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
#plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/C_Q_soft.png", transparent=True)
#plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/C_Q_soft.svg", transparent=True)
#plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/mean_velocity_1.png", transparent=True)
plt.show()


'''
#other details for plotting
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
