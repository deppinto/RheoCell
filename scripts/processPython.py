import sys
import numpy as np
from math import *
import matplotlib.pyplot as plt
import os.path
import matplotlib.cm as cm
from scipy.stats import linregress


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
    if len(save) == 0:
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
avg_FTLE = np.zeros(end-start)

#voids_histogram_area.txt  voids_histogram_circularity.txt  voids_histogram_radius_speed_time.txt  voids_histogram_radius_speed.txt  voids_histogram_radius.txt  voids_stats.txt

fig_histograms , ax = plt.subplots(2, 3, figsize=(10,4), constrained_layout=True)
#fig_histograms.tight_layout()


if scripts==74:
    jobs[14] -= 1
if scripts==75:
    jobs[12] -= 1
if scripts==76:
    jobs[3] -= 1
    jobs[4] -= 3


for traj in range(start, end):

    area_histogram = np.zeros(int(0.25 * time[traj] / deltat[traj]))
    circularity_histogram = np.zeros(int(0.25 * time[traj] / deltat[traj]))
    radius_speed_histogram = np.zeros(int(0.25 * time[traj] / deltat[traj]))
    ani_histogram = np.zeros(int(0.25 * time[traj] / deltat[traj]))
    rspeed_histogram = np.zeros(int(lx[traj]))

    avg_stress_histogram = np.zeros(10)
    max_stress_histogram = np.zeros(10)

    avg_distance_velocity_defects = np.zeros(10)

    l_bins = 0.1
    x_l = [i * l_bins for i in range(int(1./l_bins))]
    lifetime_hist = np.zeros(int(1. / l_bins))
    a_bins = 0.02
    x_a = [i * a_bins for i in range(int(1./a_bins))]
    max_area_hist = np.zeros(int(1. / a_bins))
    total_hist_counts = 0

    sp_bins = 0.04
    x_sp = [i * sp_bins for i in range(int(1./sp_bins))]
    survival_probability_hist = np.zeros(int(1. / sp_bins))
    survival_probability_counts = np.zeros(int(1. / sp_bins))

    strain_bins = 10
    strain_min = 0.
    strain_max = 0.004
    strain_delta = strain_max - strain_min
    delta_bin = strain_delta / strain_bins
    x_st = [(i * delta_bin + strain_min + delta_bin/2) for i in range(strain_bins)]
    strain_hist = np.zeros(strain_bins)
    strain_count = np.zeros(strain_bins)


    strain_bins_20 = 20
    x_st_20 = np.zeros(strain_bins_20)
    strain_hist_20 = np.zeros(strain_bins_20)
    strain_count_20 = np.zeros(strain_bins_20)

    strain_bins_30 = 30
    x_st_30 = np.zeros(strain_bins_30)
    strain_hist_30 = np.zeros(strain_bins_30)
    strain_count_30 = np.zeros(strain_bins_30)

    strain_bins_40 = 40
    x_st_40 = np.zeros(strain_bins_40)
    strain_hist_40 = np.zeros(strain_bins_40)
    strain_count_40 = np.zeros(strain_bins_40)


    count_jobs = 0
    for job in range(jobs_seq[traj], jobs_seq[traj+1]):

        if job==144 and scripts==74:
            continue
        if job==125 and scripts==75:
            continue
        if job==39 and scripts==76:
            continue
        if job==43 and scripts==76:
            continue
        if job==49 and scripts==76:
            continue
        if job==50 and scripts==76:
            continue
        if job==247 and scripts==85:
            continue
        if job==250 and scripts==85:
            continue
        if job==349 and scripts==85:
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
                sp = (float(save[6+i]) / (time[traj] * dt[traj])) / sp_bins

                survival_probability_hist[0:int(sp)+1] += 1.
                survival_probability_counts[int(sp)] += 1.
                lifetime_hist[int(hlife)] += 1
                max_area_hist[int(maxA)] += 1
                total_hist_counts += 1

        fileoutput.close()
        if nn > 0:
            count_jobs += 1


        fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_histogram.txt","r")
        for line in fileoutput:
            save=line.split()
            hist_step=int(int(float(save[0]) * 50.01) * 0.5)
            area_histogram[hist_step] += (float(save[1]) * pi * 8 * 8) / (lx[traj] * ly[traj]) * 0.5
            circularity_histogram[hist_step] += float(save[2]) * 0.5
            radius_speed_histogram[hist_step] += float(save[3]) * 0.5
            ani_histogram[hist_step] += float(save[4]) * 0.5
        fileoutput.close()

        '''
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
        '''

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

        fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_velocity_histogram_tau10.txt","r")
        for line in fileoutput:
            save=line.split()
            avg_distance_velocity_defects[int(float(save[0]))] += float(save[1])
        fileoutput.close()

        fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_velocity_stats.txt","r")
        for line in fileoutput:
            save=line.split()
            avg_FTLE[traj - start] += float(save[2]) / jobs[traj]
        fileoutput.close()


        fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_strain_histogram_20bins.txt","r")
        cont_line = 0
        for line in fileoutput:
            save=line.split()
            x_st_20[cont_line] = float(save[0])
            strain_hist_20[cont_line] += float(save[1])
            strain_count_20[cont_line] += float(save[2])
            cont_line+=1
        fileoutput.close()

        fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_strain_histogram_30bins.txt","r")
        cont_line = 0
        for line in fileoutput:
            save=line.split()
            x_st_30[cont_line] = float(save[0])
            strain_hist_30[cont_line] += float(save[1])
            strain_count_30[cont_line] += float(save[2])
            cont_line+=1
        fileoutput.close()

        fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_strain_histogram_40bins.txt","r")
        cont_line = 0
        for line in fileoutput:
            save=line.split()
            x_st_40[cont_line] = float(save[0])
            strain_hist_40[cont_line] += float(save[1])
            strain_count_40[cont_line] += float(save[2])
            cont_line+=1
        fileoutput.close()


        fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_strain_defectsplusone.txt","r")
        old_t = 0
        stmax = 0.
        old_index = 0
        nh = 0
        for line in fileoutput:
            save=line.split()
            if int(float(save[2])) > old_t + 3 and old_t > 0:
                if old_index >= strain_bins:
                    old_index = strain_bins - 1
                #strain_hist[old_index] += 1.
                strain_hist[old_index:strain_bins] += 1
                nh += 1
            #elif int(float(save[2])) == old_t + 1 and float(save[1])>0:
                #if old_index >= strain_bins:
                    #old_index = strain_bins - 1
                #strain_hist[old_index] += 1.
                #strain_hist[old_index:strain_bins] += 1

            if float(save[0])>strain_min:
                index = (float(save[0]) - strain_min)/delta_bin
                if index >= strain_bins:
                    index = strain_bins - 1
                strain_count[int(index)] += 1.
                if int(float(save[2])) > old_t:
                    stmax = 0.
                if float(save[0]) > stmax:
                    old_index = int(index)
                    stmax = float(save[0])
            old_t = int(float(save[2]))
        fileoutput.close()
        if old_t < 100:
            #strain_hist[old_index] += 1
            strain_hist[old_index:strain_bins] += 1
            nh += 1
        #print(nn, nh)

    if count_jobs > 0:
        area_histogram[:] = area_histogram[:] / count_jobs
        circularity_histogram[:] = circularity_histogram[:] / count_jobs
        radius_speed_histogram[:] = radius_speed_histogram[:] / count_jobs
        rspeed_histogram[:] = rspeed_histogram[:] / count_jobs
        ani_histogram[:] = ani_histogram[:] / count_jobs

        avg_stress_histogram[:] = avg_stress_histogram[:] / count_jobs
        max_stress_histogram[:] = max_stress_histogram[:] / count_jobs

        avg_distance_velocity_defects[:] = avg_distance_velocity_defects[:] / count_jobs

        '''
        for q in range(len(strain_count_20)):
            if strain_count_20[q]>0:
                strain_hist_20[q] = strain_hist_20[q] / strain_count_20[q]

        for q in range(len(strain_count_30)):
            if strain_count_30[q]>0:
                strain_hist_30[q] = strain_hist_30[q] / strain_count_30[q]

        for q in range(len(strain_count_40)):
            if strain_count_40[q]>0:
                strain_hist_40[q] = strain_hist_40[q] / strain_count_40[q]
        '''


    if variable == 1:
        if total_hist_counts>0:
            lifetime_hist[:] = lifetime_hist[:] / total_hist_counts
            max_area_hist[:] = max_area_hist[:] / total_hist_counts
        ax[1,0].plot(x_l, lifetime_hist, '-o')
        ax[1,1].plot(x_a, max_area_hist, '-s', label=nemself[traj])

    if variable == 2:
        x_lt = []
        for rr in range(len(area_histogram)):
            x_lt.append(rr/len(area_histogram))

        ax[0,0].plot(x_lt, area_histogram, '-o')
        ax[0,1].plot(x_lt, circularity_histogram, '-s')
        ax[1,0].plot(x_lt, radius_speed_histogram, '-^')

        #ax[1,1].plot(rspeed_histogram, '-v')

        '''
        for q in range(len(strain_hist)):
            #if strain_count[q] > 0:
                #strain_hist[q] /= strain_count[q]
            if total_hist_counts>0:
                strain_hist[q] /= total_hist_counts
            if 1 - strain_hist[q] > 0:
                strain_hist[q] = 1 - strain_hist[q]
            else:
                strain_hist[q] = 0
            strain_hist[q] = 1 - strain_hist[q]
        ax[1,1].plot(x_st, strain_hist, '-v')
        '''

        xx = x_st_40[0:]
        yy = strain_hist_40[0:]
        #for q in range(len(yy)):
            #if yy[q] > 0:
                #yy[q] = -log(yy[q])
        ax[1,1].plot(xx, yy, '-v')

        ax[1,2].plot(x_lt, ani_histogram, '-o', label=nemself[traj], color='firebrick')
        ax2 = ax[1,2].twinx()
        ax2.plot(x_lt, area_histogram, '--o', color='forestgreen')
        ax2.set_ylabel('Normalized area')
        ax[1,2].tick_params(axis='y', colors='firebrick')
        ax2.tick_params(axis='y', colors='forestgreen')
        ax[1,2].yaxis.label.set_color('firebrick')
        ax2.yaxis.label.set_color('forestgreen')
        ax2.spines['right'].set_color('forestgreen')
        ax2.spines['left'].set_color('firebrick')

    if variable == 3:
        if total_hist_counts>0:
            survival_probability_hist[:] = survival_probability_hist[:] / total_hist_counts

        ax[0,0].plot(avg_stress_histogram, '-o')
        ax[0,1].plot(max_stress_histogram, '-s')
        ax[1,0].plot(avg_distance_velocity_defects, '-^', label=nemself[traj])

        xx = x_sp
        yy = np.zeros(len(x_sp))
        Ntotal = sum(survival_probability_counts[:])
        final_bin = -1
        if Ntotal > 0:
            for q in range(len(x_sp)):
                yy[q] = (Ntotal - sum(survival_probability_counts[0:q])) / Ntotal
                if yy[q] < 0.065 and final_bin == -1:
                    final_bin = q

        ax[1,2].plot(xx, yy, '-o')
        final_bin=-3
        if final_bin > -2:
            S_sigma = yy[1:final_bin]
            S_sigma[S_sigma <= 1e-10] = 1e-10
            lnlnS = np.log(-np.log(S_sigma))
            #lnlnS = np.log(S_sigma)
            lnsigma = xx[1:final_bin]
            lnsigma = np.log(lnsigma)
            mask = ~np.isnan(lnlnS)
            slope, intercept, _, _, _ = linregress(lnsigma[mask], lnlnS[mask])
            m = slope
            sigma_0 = np.exp(-intercept / slope) 
            print(f"First Weibull shape (m): {m:.2f}")
            print(f"First Weibull scale (σ₀): {sigma_0:.5f}")

            #for q in range(len(yy)):
                #if yy[q] > 0:
                    #yy[q] = -log(yy[q])

            ax[1,2].plot(lnsigma, lnlnS, '-v')
        #ax[1,2].plot(xx, yy, '-o')
        #ax[1,2].plot(x_sp, survival_probability_hist, '-o')


        xx = x_st_40[0:]
        #yy = 1 - strain_hist_20[0:]

        yy = np.zeros(strain_bins_40)
        Ntotal = sum(strain_hist_40[:])
        final_bin = -1
        if Ntotal > 0:
            for q in range(strain_bins_40):
                yy[q] = (Ntotal - sum(strain_hist_40[0:q])) / Ntotal
                if yy[q] < 0.065 and final_bin == -1:
                    final_bin = q

        
        #ax[1,1].plot(xx, yy, '-v')
        #final_bin = -3
        if final_bin > -2:
            S_sigma = yy[1:final_bin]
            S_sigma[S_sigma <= 1e-10] = 1e-10
            minuslnS = -np.log(S_sigma)
            lnlnS = np.log(-np.log(S_sigma))
            sigma = xx[1:final_bin]
            lnsigma = np.log(sigma)
            mask = ~np.isnan(lnlnS)
            #print(lnsigma[mask], lnlnS[mask])
            slope, intercept, _, _, _ = linregress(lnsigma[mask], lnlnS[mask])
            m = slope
            sigma_0 = np.exp(-intercept / slope) 
            print(f"Weibull shape (m): {m:.2f}")
            print(f"Weibull scale (σ₀): {sigma_0:.5f}")

            #for q in range(len(yy)):
                #if yy[q] > 0:
                    #yy[q] = -log(yy[q])

            #ax[1,1].plot(xx, yy, '-v')
            #ax[1,1].plot(lnsigma, lnlnS, '-v')
            ax[1,1].plot(sigma, minuslnS, '-v')

            line_plot = []
            if traj == end - 1:
                for q in range(len(lnsigma)):
                    line_plot.append(np.exp(2 * lnsigma[q] + intercept - 3))
                ax[1,1].plot(sigma, line_plot, '-')


'''
PD_Activity = []
PD_Friction = []
PD_Phase = []
for i in range(len(lifetime_holes)):
    if lifetime_holes[i]<0.01:
        PD_Activity.append(nemself[i])
        PD_Friction.append(fric[i])
        PD_Phase.append(0)
    elif lifetime_holes[i] < 0.5:
        PD_Activity.append(nemself[i])
        PD_Friction.append(fric[i])
        PD_Phase.append(1)
    else:
        PD_Activity.append(nemself[i])
        PD_Friction.append(fric[i])
        PD_Phase.append(2)

print(PD_Activity)
print(PD_Friction)
print(PD_Phase)
'''

if variable == 2:
    ax[0,0].set_ylabel('Area')
    ax[0,1].set_ylabel('Circularity')
    ax[1,0].set_ylabel('Radius speed')
    #ax[1,1].set_ylabel('Radius speed')
    ax[1,1].set_ylabel(r'F($\varepsilon$)')
    ax[1,2].set_ylabel('Anisotropy ratio')
    ax[0,0].set_xlabel('Lifetime')
    ax[0,1].set_xlabel('Lifetime')
    ax[1,0].set_xlabel('Lifetime')
    #ax[1,1].set_xlabel('Radius')
    ax[1,1].set_xlabel(r'$\varepsilon$')
    ax[1,2].set_xlabel('Lifetime')

    ax[1,1].set_yscale('log')
    #ax[1,1].set_xscale('log')

    fig_histograms.delaxes(ax[0,2])
    fig_histograms.legend(loc=(0.7, 0.6), ncols=2, frameon=False)
    #plt.legend(loc="upper right", frameon='false')
    plt.show()
    exit (1)


if variable==1:
    ax[0,0].set_ylabel('Lifetime')
    ax[0,1].set_ylabel('Max area')
    ax[1,0].set_ylabel('P(Lifetime)')
    ax[1,1].set_ylabel('P(Max area)')
    ax[0,0].set_xlabel(r'$\gamma$')
    ax[0,1].set_xlabel(r'$\gamma$')
    ax[1,0].set_xlabel('Normalized lifetime')
    ax[1,1].set_xlabel('Normalized max area')
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox='false', ncol=5)
    #fig.legend(loc=7)
    ax[0,0].plot(nemself[0:5], lifetime_holes[0:5], '--p')
    ax[0,0].plot(nemself[5:10], lifetime_holes[5:10], '--p')
    #ax[0,0].plot(nemself[10:15], lifetime_holes[10:15], '--p')
    #ax[0,0].plot(lx[5:8], lifetime_holes[0:3], '--p')
    #ax[0,0].plot(mu[0:4], lifetime_holes[0:4], '--p')
    #ax[0,0].plot(gamma[4:8], lifetime_holes[0:4], '--p')
    #ax[0,0].plot(omega[0:6], lifetime_holes[0:6], '--p')
    #ax[0,0].plot(omega[7:13], lifetime_holes[7:13], '--p')
    #ax[0,0].plot(fric, lifetime_holes, '--p')
    #ax[0,0].set_xscale('log')
    #ax[0,0].plot(F[0:5], lifetime_holes[0:5], '--p')
    #ax[0,0].plot(F[5:10], lifetime_holes[5:10], '--p')
    #ax[0,0].plot(F[10:15], lifetime_holes[10:15], '--p')

    ax[0,1].plot(nemself[0:5], max_area_holes[0:5], '--v')
    ax[0,1].plot(nemself[5:10], max_area_holes[5:10], '--v')
    #ax[0,1].plot(nemself[10:15], max_area_holes[10:15], '--v')
    #ax[0,1].plot(lx[5:8], max_area_holes[0:3], '--v')
    #ax[0,1].plot(omega[0:6], max_area_holes[0:6], '--v')
    #ax[0,1].plot(omega[7:13], max_area_holes[7:13], '--v')
    #ax[0,1].plot(mu[0:4], max_area_holes[0:4], '--v')
    #ax[0,1].plot(gamma[4:8], max_area_holes[0:4], '--v')
    #ax[0,1].plot(F[0:5], max_area_holes[0:5], '--v')
    #ax[0,1].plot(F[5:10], max_area_holes[5:10], '--v')
    #ax[0,1].plot(F[10:15], max_area_holes[10:15], '--v')
    #ax[0,1].plot(fric, max_area_holes, '--v')
    #ax[0,1].set_xscale('log')

    ax[0,2].plot(nemself[0:5], n_holes[0:5], '--^')
    ax[0,2].plot(nemself[5:10], n_holes[5:10], '--^')
    #ax[0,2].plot(nemself[10:15], n_holes[10:15], '--^')
    #ax[0,2].plot(lx[5:8], n_holes[0:3], '--^')
    #ax[0,2].plot(omega[0:6], n_holes[0:6], '--^')
    #ax[0,2].plot(omega[7:13], n_holes[7:13], '--^')
    #ax[0,2].plot(mu[0:4], n_holes[0:4], '--^')
    #ax[0,2].plot(gamma[4:8], n_holes[0:4], '--^')
    #ax[0,2].plot(F[0:5], n_holes[0:5], '--^')
    #ax[0,2].plot(F[5:10], n_holes[5:10], '--^')
    #ax[0,2].plot(F[10:15], n_holes[10:15], '--^')
    #ax[0,2].plot(fric, n_holes, '--^')
    #ax[0,2].set_xscale('log')
    ax[0,2].set_ylabel('# voids')
    ax[0,2].set_xlabel(r'$\gamma$')

    fig_histograms.delaxes(ax[1,2])
    fig_histograms.legend(loc=(0.75,0.15), ncols=3, frameon=False)

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
    ax[0,1].set_ylabel('Normalized max strain')
    ax[0,0].set_xlabel('-Time')
    ax[0,1].set_xlabel('-Time')
    ax[1,0].set_ylabel('Vortex distance')
    ax[1,0].set_xlabel('-Time')
    #ax[1,1].set_ylabel('Average FTLE')
    #ax[1,1].set_xlabel(r'$\zeta$')
    ax[1,2].set_ylabel(r'S(t)')
    ax[1,2].set_xlabel('t')
    #ax[1,1].set_ylabel(r'S($\varepsilon$)')
    #ax[1,1].set_xlabel(r'$\varepsilon$')
    ax[1,1].set_ylabel(r'$-log(S(\varepsilon))$')
    ax[1,1].set_xlabel(r'$\varepsilon$')

    #ax[1,2].set_yscale('log')
    #ax[1,2].set_xscale('log')

    ax[1,1].set_yscale('log')
    ax[1,1].set_xscale('log')
    ax[1,1].set_ylim([1e-3, 10])

    #ax[1,1].plot(nemself[0:5], avg_FTLE[0:5], '--p')
    #ax[1,1].plot(nemself[5:10], avg_FTLE[5:10], '--p')
    #ax[1,1].plot(nemself[10:15], avg_FTLE[10:15], '--p')
    #ax[1,1].plot(lx[5:8], avg_FTLE[0:3], '--p')
    #ax[1,1].plot(fric, avg_FTLE, '--p')
    #ax[1,1].set_xscale('log')
    #ax[1,1].plot(omega[0:6], avg_FTLE[0:6], '--p')
    #ax[1,1].plot(mu[0:4], avg_FTLE[0:4], '--p')
    #ax[1,1].plot(gamma[4:8], avg_FTLE[0:4], '--p')
    fig_histograms.delaxes(ax[0,2])
    fig_histograms.legend(loc=(0.7, 0.6), ncols=2, frameon=False)
    plt.show()
    exit(1)

plt.close()

if variable == 4:

    fig = plt.figure(figsize=(5.452423529,4.089317647), constrained_layout=True)
    PD_activity = [0.1, 0.2, 0.3, 0.4, 0.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.1, 0.2, 0.3, 0.4, 0.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.1, 0.2, 0.3, 0.4, 0.5]
    PD_friction = [0.1, 0.1, 0.1, 0.1, 0.1, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.75, 0.5, 0.25, 0.1, 0.075, 0.05, 0.025, 0.01, 0.0075, 0.005, 0.0025, 0.001, 0.05, 0.05, 0.05, 0.05, 0.05, 0.025, 0.025, 0.025, 0.025, 0.025, 0.0075, 0.0075, 0.0075, 0.0075, 0.0075, 0.005, 0.005, 0.005, 0.005, 0.005, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025]
    PD_phase = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 2, 0, 0, 1, 1, 2, 0, 1, 1, 2, 2]

    colors = cm.Set1(np.linspace(0, 1, 9))
    for i in range(len(PD_friction)):
        plt.scatter(PD_activity[i], PD_friction[i], color=colors[PD_phase[i]])

    plt.ylabel(r'$\chi$', fontsize=18)
    plt.xlabel(r'$\zeta$', fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xticks(np.arange(0.1,0.6,0.1))
    plt.yscale('log')
    plt.ylim([0.0007, 0.15])
    plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/Slides/Results7/Phase_diagram_1.png", transparent=True)
    plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/Slides/Results7/Phase_diagram_1.svg", transparent=True)
    plt.show()

    
#first figure-----------------------------------------------------------------------------------------------
#fig1 = plt.figure(figsize=(5.452423529,4.089317647))
#plt.plot(nemself[0:5], n_holes[0:5], '--o')
#plt.plot(nemself[5:10], n_holes[5:10], '--s')
#plt.plot(nemself[10:15], n_holes[10:15], '--^')
#plt.legend([r'$\xi$=0.1',r'$\xi$=0.01',r'$\xi$=0.001'], fontsize=14, frameon=False)

#plt.plot(fric, n_holes, '--^')
#plt.xscale('log')
#plt.yscale('log')

#plt.plot(omega[0:7], n_holes[0:7], '--o')
#plt.plot(omega[7:14], n_holes[7:14], '--s')
#plt.legend([r'$\xi$=0.01',r'$\xi$=0.001'], fontsize=14, frameon=False)

#plt.ylabel('# voids', fontsize=18)
#plt.xlabel('P', fontsize=18)




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

#plt.show()


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
