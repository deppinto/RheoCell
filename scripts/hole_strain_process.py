import sys
import numpy as np
from math import *
import matplotlib.pyplot as plt
import os.path
import matplotlib.cm as cm
from scipy.stats import linregress
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import font_manager, rcParams
# Load font from file
#font_path = "/home/p/pinto/Fonts/Times_New_Roman_Normal.ttf"
font_path = "/home/p/pinto/Fonts/times.ttf"
italic_font_path = "/home/p/pinto/Fonts/timesi.ttf"
bold_font_path = "/home/p/pinto/Fonts/timesbd.ttf"

custom_font = FontProperties(fname=font_path)
legend_font = FontProperties(fname=font_path, size=12)
font = font_manager.FontProperties(fname=font_path)
fonti = font_manager.FontProperties(fname=italic_font_path)
fontbd = font_manager.FontProperties(fname=bold_font_path)

# Register font with a name
font_manager.fontManager.addfont(font_path)
font_manager.fontManager.addfont(italic_font_path)
font_manager.fontManager.addfont(bold_font_path)

# Set custom mathtext font to match your font
rcParams['mathtext.fontset'] = 'custom'

# Set roman (upright), italic, and bold versions (all Times New Roman if needed)
#rcParams['mathtext.rm'] = font.get_name()  # e.g. "Times New Roman"
#rcParams['mathtext.it'] = fonti.get_name()
#rcParams['mathtext.bf'] = fontbd.get_name()
rcParams['mathtext.it'] = 'Times New Roman:italic'
rcParams['mathtext.rm'] = 'Times New Roman'
rcParams['mathtext.bf'] = 'Times New Roman:bold'


if len(sys.argv)!=5:
        print(sys.argv[0]," [scripts] [start] [end] [variable]")
        sys.exit(1)

scripts=int(float(sys.argv[1]))
start=int(float(sys.argv[2]))
end=int(float(sys.argv[3]))
variable=int(float(sys.argv[4]))

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
J_Q=[]
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
    #J_Q.append(float(save[22]))

    #F.append( sqrt( float(save[10])/( float(save[11]) )) ) 

    cont_line+=1
filedata.close()


#MPF holes paper----------------------------------------------------------------------------------------------
# read output file

st_all_x = []
st_all_y = []
log_st_all_x = []
log_st_all_y = []

fit_x = []
fit_y = []
slope_fit = []


Avg_rate = []
Area = 306.25
Time = 1000 * 0.1 

for traj in range(start, end):

    avg_distance_velocity_defects = np.zeros(10)

    strain_bins = 40
    strain_min = 0.
    strain_max = 0.005
    strain_delta = strain_max - strain_min
    delta_bin = strain_delta / strain_bins
    x_st = [(i * delta_bin + strain_min + delta_bin/2) for i in range(strain_bins)]
    N_strain_hole = np.zeros(strain_bins)
    N_strain = np.zeros(strain_bins)
    N_all = 0
    N_holes = 0

    count_jobs = 0
    n_hole = 0


    new_N_strain_hole = np.zeros(strain_bins)
    new_N_strain = np.zeros(strain_bins)
    new_N_all = 0
    new_N_holes = 0

    Rate_traj = []

    for job in range(jobs_seq[traj], jobs_seq[traj+1]):

        if job==144 and scripts==74:
            continue
        if scripts==85 and os.path.isfile("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_strain_histogram_"+str(strain_bins)+"bins.txt") == False:
            continue
        if scripts==83 and os.path.isfile("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_strain_histogram_"+str(strain_bins)+"bins.txt") == False:
            continue
        if scripts==84 and os.path.isfile("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_strain_histogram_"+str(strain_bins)+"bins.txt") == False:
            continue

        fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_test_plusonestrain.txt","r")
        cont_line = 0
        strain_values = []
        hole_strain = []
        hole_size = []
        time_step = []
        hole_posx = []
        hole_posy = []
        list_strain = []
        last_timeframe = 0
        hole_formed = 0
        for line in fileoutput:
            save=line.split()
            if save[0] == "hole":
                n_hole += 1

                index = len(strain_values) - 1
                timeframe = time_step[index]

                #print(timeframe, last_timeframe, traj, job)
                if(timeframe >  last_timeframe + 1):
                    Rate_traj.append(1/(Area * (timeframe - last_timeframe) * Time))
                hole_formed = 1

                max_strain = 0
                posx = 0.
                posy = 0.
                while time_step[index] == timeframe:
                    if strain_values[index] > max_strain:
                        max_strain = strain_values[index]
                        posx = hole_posx[index]
                        posy = hole_posy[index]
                    index -= 1
                    if index < 0:
                        break
                        #print("ERROR1 ", timeframe)
                        #exit(1)

                if index < 0:
                    continue

                list_strain.append(max_strain)
                N_holes += 1
                N_strain_hole[int( (max_strain - strain_min)/delta_bin )] += 1
                timeframe = time_step[index]

                while hole_size[index] > 0:
                    while time_step[index] == timeframe:
                        distx = posx - hole_posx[index]
                        if distx >= lx[traj] / 2:
                            distx -= lx[traj]
                        if distx <= -lx[traj] / 2:
                            distx += lx[traj]

                        disty = posy - hole_posy[index]
                        if disty >= ly[traj] / 2:
                            disty -= ly[traj]
                        if disty <= -ly[traj] / 2:
                            disty += ly[traj]

                        dist = distx * distx + disty * disty
                        if dist < 64:
                            list_strain.append(strain_values[index])
                            N_holes += 1
                            N_strain_hole[int( (strain_values[index] - strain_min)/delta_bin )] += 1
                            posx = hole_posx[index]
                            posy = hole_posy[index]

                        index -= 1
                        if index < 0:
                            print("ERROR2")
                            exit (1)

                    timeframe = time_step[index]


                timeframe = time_step[index]
                while time_step[index] == timeframe:
                    distx = posx - hole_posx[index]
                    if distx >= lx[traj] / 2:
                        distx -= lx[traj]
                    if distx <= -lx[traj] / 2:
                        distx += lx[traj]

                    disty = posy - hole_posy[index]
                    if disty >= ly[traj] / 2:
                        disty -= ly[traj]
                    if disty <= -ly[traj] / 2:
                        disty += ly[traj]

                    dist = distx * distx + disty * disty
                    if dist < 64:
                        list_strain.append(strain_values[index])
                        N_holes += 1
                        N_strain_hole[int( (strain_values[index] - strain_min)/delta_bin )] += 1
                        posx = hole_posx[index]
                        posy = hole_posy[index]

                    index -= 1
                    if index < 0:
                        break
                        #print("ERROR3")
                        #exit (1)


            elif save[0] == "final:":
                #print("end of file")
                break
            else:
                if hole_formed == 1:
                    hole_formed = 0
                    last_timeframe = int(float(save[2]))

                if float(save[0]) > strain_min:
                    strain_read = float(save[0])
                    strain_values.append(strain_read)
                    N_all += 1
                    #print(strain_read, job, traj)
                    N_strain[int( (strain_read - strain_min)/delta_bin )] += 1
                    hole_size.append(int(float(save[1])))
                    time_step.append(int(float(save[2])))
                    hole_posx.append(float(save[3]))
                    hole_posy.append(float(save[4]))
            cont_line+=1

        fileoutput.close()


        #fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_strain_histogram_40bins.txt","r")
        fileoutput=open("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_strain_histogram_"+str(strain_bins)+"bins.txt","r")
        cont_line = 0
        for line in fileoutput:
            save=line.split()
            new_N_strain_hole[cont_line] += float(save[1])
            new_N_holes += float(save[1])
            new_N_strain[cont_line] += float(save[2])
            new_N_all += float(save[2])
            #print(cont_line, float(save[1]), float(save[2]))
            cont_line+=1
        fileoutput.close()

    total_number_holes = len(Rate_traj)
    total_rate = 0.
    for i in Rate_traj:
        total_rate += i
    Avg_rate.append(total_rate / total_number_holes)


    P_final = []
    x_final = []
    PP_final = []
    xx_final = []
    for i in range(strain_bins):
        '''
        if N_strain[i] > 0:
            P_hole_E = N_strain_hole[i] / N_strain[i]
            #P_hole_E = N_strain_hole[i] / N_holes
            #P_E = N_strain[i] / N_all
            #P_final.append(P_hole_E / P_E)
            P_final.append(P_hole_E)
            x_final.append(x_st[i])
            if N_strain_hole[i] > 0:
                #PP_final.append(P_hole_E / P_E)
                PP_final.append(P_hole_E)
                #P_final.append(P_hole_E)
                xx_final.append(x_st[i])
        '''
        if new_N_strain[i] > 0 and x_st[i]<0.005:
            P_hole_E = new_N_strain_hole[i] / new_N_strain[i]
            #print(new_N_strain_hole[i], new_N_strain[i])
            #P_hole_E = new_N_strain_hole[i] / new_N_holes
            #P_E = new_N_strain[i] / new_N_all
            #P_final.append(P_hole_E / P_E)
            P_final.append(P_hole_E)
            x_final.append(x_st[i])
            if new_N_strain_hole[i] > 0:
                #PP_final.append(P_hole_E / P_E)
                PP_final.append(P_hole_E)
                #P_final.append(P_hole_E)
                xx_final.append(x_st[i])

    st_all_x.append(x_final)
    st_all_y.append(P_final)
    log_st_all_x.append(xx_final)
    log_st_all_y.append(PP_final)

    xx = np.array(xx_final)
    yy = np.array(PP_final)
    final_bin = len(PP_final)
    if final_bin > -2:
        S_sigma = yy[1:final_bin] #/ (nemself[traj]/fric[traj]**0.1)
        #S_sigma[S_sigma <= 1e-10] = 1e-10
        lnS = np.log10(S_sigma)
        #lnS = S_sigma
        sigma = xx[1:final_bin] #/ (nemself[traj]/fric[traj]**0.1)
        #lnsigma = np.log10(sigma)
        lnsigma = sigma
        mask = ~np.isnan(lnS)

        if final_bin > 5:
            #slope, intercept, _, _, _ = linregress(sigma[mask], lnS[mask])
            slope, intercept, _, _, _ = linregress(lnsigma[mask], lnS[mask])
            m = slope
            slope_fit.append(m)
            #sigma_0 = np.exp(-intercept / slope) 
            #sigma_0 = intercept 
            #print(f"Weibull shape (m): {m:.2f}")
            #print(f"Weibull scale (σ₀): {sigma_0:.5f}")
            print("slope:\t", m, "\tintercept:\t", intercept, "\tFriction:\t",  fric[traj], "\tActivity:\t", nemself[traj])

        fit_x.append(lnsigma[mask])
        fit_y.append(lnS[mask])

        #line_plot = []
        #if traj == end - 1:
            #for q in range(len(lnsigma)):
                #line_plot.append(np.exp(2 * lnsigma[q] + intercept - 3))


plt.figure(figsize=(5.452423529,4.089317647))
#plt.figure(figsize=(3*40.179/25.4, 3*25.142/25.4))
#plt.plot(fit_x[2], fit_y[2], ':o', color='firebrick', markerfacecolor='none')
#plt.plot(fit_x[3], fit_y[3], ':^', color='green', markerfacecolor='none')
#plt.plot(fit_x[4], fit_y[4], ':s', color='royalblue', markerfacecolor='none')
#plt.plot(fit_x[5], fit_y[5], '--p', color='goldenrod', label='0.1')
#plt.plot(fit_x[6], fit_y[6], '--h', color='peru', label='0.2')
#plt.plot(fit_x[7], fit_y[7], '--o', color='firebrick', label='0.3')
#plt.plot(fit_x[8], fit_y[8], '--^', color='green', label='0.4')
#plt.plot(fit_x[9], fit_y[9], '--s', color='royalblue', label='0.5')


#fit_x = nemself
#for i in range(len(nemself)):
#    fit_x[i] = nemself[i] / fric[i]**0.12
#fit_y = Avg_rate

#plt.plot(fit_x[8], fit_y[8], ':o', color='firebrick', markerfacecolor='none')
#plt.plot(fit_x[7], fit_y[7], ':^', color='green', markerfacecolor='none')
#plt.plot(fit_x[6], fit_y[6], ':s', color='royalblue', markerfacecolor='none')
#plt.plot(fit_x[5], fit_y[5], '--o', color='firebrick', alpha=0.5)
#plt.plot(fit_x[4], fit_y[4], '--^', color='green', alpha=0.5)
#plt.plot(fit_x[3], fit_y[3], '--s', color='royalblue', alpha=0.5)
#plt.plot(fit_x[2], fit_y[2], '-.o', color='firebrick', label='0.5')
#plt.plot(fit_x[1], fit_y[1], '-.^', color='green', label='0.4')
#plt.plot(fit_x[0], fit_y[0], '-.s', color='royalblue', label='0.3')

#plt.plot(st_all_x[8], st_all_y[8], ':o', color='firebrick', markerfacecolor='none')
#plt.plot(st_all_x[7], st_all_y[7], ':^', color='green', markerfacecolor='none')
#plt.plot(st_all_x[6], st_all_y[6], ':s', color='royalblue', markerfacecolor='none')
#plt.plot(st_all_x[5], st_all_y[5], '--o', color='firebrick', alpha=0.5)
#plt.plot(st_all_x[4], st_all_y[4], '--^', color='green', alpha=0.5)
#plt.plot(st_all_x[3], st_all_y[3], '--s', color='royalblue', alpha=0.5)
#plt.plot(st_all_x[2], st_all_y[2], '-.o', color='firebrick', label='0.5')
#plt.plot(st_all_x[1], st_all_y[1], '-.^', color='green', label='0.4')
#plt.plot(st_all_x[0], st_all_y[0], '-.s', color='royalblue', label='0.3')

'''
log_st_all_x[0] = [6.25e-05, 0.0001875, 0.0003125, 0.0004375, 0.0005625000000000001, 0.0006875, 0.0008125000000000001, 0.0009375, 0.0010625, 0.0011875000000000002, 0.0013125, 0.0014375, 0.0015625, 0.0016875000000000002, 0.0018125, 0.0019375, 0.0020625, 0.0021875, 0.0023125000000000003, 0.0024375, 0.0025625, 0.0026875000000000002, 0.0028125, 0.0029375, 0.0030625, 0.0031875000000000002, 0.0033125000000000003, 0.0034375, 0.0035625, 0.0036875000000000002, 0.0038125, 0.0039375, 0.0040625, 0.0041875, 0.0043125, 0.0044375000000000005, 0.004562500000000001, 0.0046875, 0.0048125, 0.0049375]
log_st_all_y[0] = [np.float64(0.001581275084586966), np.float64(0.001599664278250175), np.float64(0.0017300817160367723), np.float64(0.0017922866263147797), np.float64(0.0019957831759532907), np.float64(0.0021539636229950935), np.float64(0.0023431274900398407), np.float64(0.0027304551998172614), np.float64(0.0028905201222472936), np.float64(0.003180059629996184), np.float64(0.0034964211319789607), np.float64(0.003998583150445322), np.float64(0.004236712853891307), np.float64(0.004883534563293353), np.float64(0.00537892754873887), np.float64(0.005814842316764218), np.float64(0.005815086343376282), np.float64(0.006249957074470643), np.float64(0.0071737130202891315), np.float64(0.007947067175040852), np.float64(0.00853772279561737), np.float64(0.009163842992047136), np.float64(0.009793793686357112), np.float64(0.010849937733499378), np.float64(0.012101309916905547), np.float64(0.013286242852337706), np.float64(0.015519852995952422), np.float64(0.016795226313491486), np.float64(0.017966308925417743), np.float64(0.021493166795084957), np.float64(0.02293186771534161), np.float64(0.02519003495533485), np.float64(0.028346456692913385), np.float64(0.03033182862120852), np.float64(0.03345169127265932), np.float64(0.037388331311348846), np.float64(0.034677525285695524), np.float64(0.0412703530549734), np.float64(0.03764288529759559), np.float64(0.041251778093883355)]

log_st_all_x[1] = [8.75e-05, 0.0002625, 0.0004375, 0.0006125, 0.0007875, 0.0009625, 0.0011374999999999998, 0.0013124999999999999, 0.0014874999999999999, 0.0016625, 0.0018375, 0.0020125, 0.0021874999999999998, 0.0023625, 0.0025375, 0.0027125, 0.0028875, 0.0030624999999999997, 0.0032375, 0.0034124999999999997, 0.0035875, 0.0037625, 0.0039375, 0.0041125, 0.0042875, 0.0044624999999999995, 0.0046375, 0.0048125, 0.0049875]
log_st_all_y[1] = [np.float64(0.000655326900963207), np.float64(0.0007433493077300022), np.float64(0.0008661321313998125), np.float64(0.0010438137519158607), np.float64(0.0011099513038385423), np.float64(0.001158513038352879), np.float64(0.0015486525592331344), np.float64(0.0015814544299777708), np.float64(0.0019622010653972573), np.float64(0.0019847364317274294), np.float64(0.0023479151219856654), np.float64(0.002523111497012105), np.float64(0.0025328554360812424), np.float64(0.003403820788835468), np.float64(0.004153005464480874), np.float64(0.0038714194518234797), np.float64(0.003911869368800949), np.float64(0.005698540923038387), np.float64(0.0063503769255981645), np.float64(0.008047373215912541), np.float64(0.009196662403313235), np.float64(0.011383081615951194), np.float64(0.013151532733778164), np.float64(0.016552716804605974), np.float64(0.012138188608776844), np.float64(0.013675845499907595), np.float64(0.019464163040989236), np.float64(0.017265385686438317), np.float64(0.01698043558508675)]

log_st_all_x[2] = [8.75e-05, 0.0002625, 0.0004375, 0.0006125, 0.0007875, 0.0009625, 0.0011374999999999998, 0.0013124999999999999, 0.0014874999999999999, 0.0016625, 0.0018375, 0.0020125, 0.0021874999999999998, 0.0023625, 0.0025375, 0.0027125, 0.0028875, 0.0030624999999999997, 0.0032375, 0.0034124999999999997, 0.0035875, 0.0037625, 0.0039375, 0.0041125, 0.0042875, 0.0044624999999999995, 0.0046375, 0.0048125, 0.0049875]
log_st_all_y[2] = [np.float64(0.0011534836613259435), np.float64(0.001580005745475438), np.float64(0.0015948003490506425), np.float64(0.001999750031246094), np.float64(0.0021750506966703737), np.float64(0.00277360657185761), np.float64(0.0032396180660806306), np.float64(0.0032616997997564046), np.float64(0.003200980726009671), np.float64(0.0039036416867841185), np.float64(0.004996949536620087), np.float64(0.005549899839085744), np.float64(0.008581257324311042), np.float64(0.008486199211383508), np.float64(0.010504096094888676), np.float64(0.011499847885609978), np.float64(0.012232639906240845), np.float64(0.015401920913731937), np.float64(0.012687188019966721), np.float64(0.017498776309348995), np.float64(0.01985035883340968), np.float64(0.022761474793077503), np.float64(0.021869742850276375), np.float64(0.029929034248688677), np.float64(0.034278959810874705), np.float64(0.02766599597585513), np.float64(0.02792956891317547), np.float64(0.02888700084961767), np.float64(0.03333333333333333)]
'''

#print(log_st_all_x[1])
#print(log_st_all_y[1])

plt.plot(log_st_all_x[8], log_st_all_y[8], ':o', color='firebrick', markerfacecolor='none')
plt.plot(log_st_all_x[7], log_st_all_y[7], ':^', color='green', markerfacecolor='none')
plt.plot(log_st_all_x[6], log_st_all_y[6], ':s', color='royalblue', markerfacecolor='none')
plt.plot(log_st_all_x[5], log_st_all_y[5], '--o', color='firebrick', alpha=0.5)
plt.plot(log_st_all_x[4], log_st_all_y[4], '--^', color='green', alpha=0.5)
plt.plot(log_st_all_x[3], log_st_all_y[3], '--s', color='royalblue', alpha=0.5)
plt.plot(log_st_all_x[2], log_st_all_y[2], '-.o', color='firebrick', label='0.5')
plt.plot(log_st_all_x[1], log_st_all_y[1], '-.^', color='green', label='0.4')
plt.plot(log_st_all_x[0], log_st_all_y[0], '-.s', color='royalblue', label='0.3')

#plt.plot(log_st_all_x[2], log_st_all_y[2], '--o', color='firebrick', label='400')
#plt.plot(log_st_all_x[1], log_st_all_y[1], '--^', color='green', label='225')
#plt.plot(log_st_all_x[0], log_st_all_y[0], '--s', color='royalblue', label='100')

#plt.legend(prop={'family':'Times New Roman', 'size':'12'}, loc=(0.02, 0.65), ncols=2, frameon=False)
plt.ylabel(r'$P(hole \vert \dot{\varepsilon})$', fontname='Times New Roman', fontsize=18)
plt.xlabel(r'$\dot{\varepsilon}$', fontname='Times New Roman', fontsize=18)
plt.xticks(fontname='Times New Roman', fontsize=18)
plt.yticks(fontname='Times New Roman',  fontsize=18)
#plt.xscale('log')
plt.yscale('log')
#plt.ylim(-0.05, 100.05)
plt.xlim(-0.0001, 0.0055)
plt.subplots_adjust(left=0.18, bottom=0.2, right=0.955, top=0.97)

plt.text(0.0011, 0.07, r'$\zeta$', fontsize=18, fontname="Times New Roman")
#plt.text(0.0011, 0.025, 'N', fontsize=18, fontname="Times New Roman")
plt.legend(prop={'family':'Times New Roman', 'size':'12'}, loc=(0.02, 0.65), ncols=2, frameon=False)

fig = plt.gcf()
ax = plt.gca()
inset_ax = fig.add_axes([0.66, 0.29, 0.28, 0.24])
inset_ax.plot(fric[8]**0.1/nemself[8], slope_fit[8], 'o', color='firebrick', ms=5, markerfacecolor='none')
inset_ax.plot(fric[7]**0.1/nemself[7], slope_fit[7], '^', color='green', ms=5, markerfacecolor='none')
inset_ax.plot(fric[6]**0.1/nemself[6], slope_fit[6], 's', color='royalblue', ms=5, markerfacecolor='none')
inset_ax.plot(fric[5]**0.1/nemself[5], slope_fit[5], 'o', color='firebrick', ms=5, alpha=0.5)
inset_ax.plot(fric[4]**0.1/nemself[4], slope_fit[4], '^', color='green', ms=5, alpha=0.5)
inset_ax.plot(fric[3]**0.1/nemself[3], slope_fit[3], 's', color='royalblue', ms=5, alpha=0.5)
inset_ax.plot(fric[2]**0.1/nemself[2], slope_fit[2], 'o', color='firebrick', ms=5)
inset_ax.plot(fric[1]**0.1/nemself[1], slope_fit[1], '^', color='green', ms=5)
inset_ax.plot(fric[0]**0.1/nemself[0], slope_fit[0], 's', color='royalblue', ms=5)
inset_ax.set_ylabel('Slope', fontname='Times New Roman', fontsize=12)
inset_ax.set_xlabel(r'$\chi^{0.1}/\zeta$', labelpad=-13, fontname='Times New Roman', fontsize=12)
inset_ax.set_xticks([1, 2])
for label in inset_ax.get_xticklabels():
    label.set_fontproperties('Times New Roman')
    label.set_fontsize(12)
for label in inset_ax.get_yticklabels():
    label.set_fontproperties('Times New Roman')
    label.set_fontsize(12)

plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/Paper/fig2_new.png", transparent=True)
plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/Paper/fig2_new.svg", transparent=True)
#plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/Paper/figSM6_new.png", transparent=True)
#plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/Paper/figSM6_new.svg", transparent=True)
#plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig2.png", transparent=True)
#plt.savefig("/home/p/pinto/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig2.svg", transparent=True)
plt.show()
exit(1)

#MPF holes paper----------------------------------------------------------------------------------------------
