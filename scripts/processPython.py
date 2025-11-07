#!/home/diogo/PythonEnv/bin/python
import sys
import numpy as np
from math import *
import matplotlib.pyplot as plt
import os.path
import matplotlib.cm as cm
from scipy.stats import linregress
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats

from statsmodels.tsa.stattools import acf

from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import font_manager, rcParams


# Load font from file
#font_path = "/home/diogo/Fonts/Times_New_Roman_Normal.ttf"
font_path = "/home/diogo/Fonts/times.ttf"
italic_font_path = "/home/diogo/Fonts/timesi.ttf"
bold_font_path = "/home/diogo/Fonts/timesbd.ttf"

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

def estimate_period(t, y, acf_threshold=0.3):
    """
    Estimate the period of a time series using FFT and autocorrelation.
    
    Parameters
    ----------
    t : array-like
        Time values (must be evenly spaced).
    y : array-like
        Signal values.
    acf_threshold : float
        Minimum autocorrelation peak height to consider as significant.
        
    Returns
    -------
    period : float or None
        Estimated period, or None if no clear period is found.
    method : str
        Method used ("fft", "acf", or "none").
    """
    
    dt = t[1] - t[0]  # assume uniform sampling
    y = np.asarray(y)
    y_detrended = y - np.mean(y)
    
    # === Method 1: FFT ===
    freqs = np.fft.rfftfreq(len(t), d=dt)
    spectrum = np.abs(np.fft.rfft(y_detrended))
    
    if len(spectrum) > 1:
        freqs, spectrum = freqs[1:], spectrum[1:]  # remove DC component
        dominant_freq = freqs[np.argmax(spectrum)]
        if dominant_freq > 0:
            period_fft = 1.0 / dominant_freq
        else:
            period_fft = None
    else:
        period_fft = None
    
    # === Method 2: Autocorrelation ===
    acorr = acf(y_detrended, nlags=len(y)//2, fft=True)
    
    # Find peaks (ignoring lag=0)
    peaks = np.where((acorr[1:-1] > acorr[:-2]) & (acorr[1:-1] > acorr[2:]))[0] + 1
    
    period_acf = None
    if len(peaks) > 0:
        best_peak = peaks[np.argmax(acorr[peaks])]
        if acorr[best_peak] > acf_threshold:
            period_acf = best_peak * dt
    
    # === Decision logic ===
    if period_fft is not None:
        return period_fft, "fft"
    elif period_acf is not None:
        return period_acf, "acf"
    else:
        return None, "none"


'''
#-------------------------------------optimization plot
plt.figure(figsize=(5.452423529,4.089317647))
fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/optimization.txt","r")
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
plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/runtime_cores.png", transparent=True)
plt.show()
exit(1)
#-------------------------------------optimization plot
'''

#filedata=open("/scratch/pinto/Phase_Field/RheoCell/Work/Results/scripts"+str(scripts)+"/dados.txt","r")
filedata=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/params","r")
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
    J_Q.append(float(save[22]))

    #F.append( sqrt( float(save[10])/( float(save[11]) )) ) 

    cont_line+=1
filedata.close()


'''
# Example datasets
circle_minus   = [
0.9544999999999999,
0.862,
0.9560000000000002,
0.9375,
0.9180000000000003,
0.9624999999999999,
0.9250000000000002,
0.91,
0.9554999999999999,
0.9605
]

circle_plus    = [
0.8825000000000001,
0.9365,
0.9440000000000001,
0.931,
0.931,
0.9265000000000001,
0.9085000000000001,
0.8925000000000001,
0.954,
0.926
]

square_minus   = [
0.792,
0.8015000000000001,
0.514,
0.7150000000000001,
0.49199999999999994,
0.5599999999999999,
0.7020000000000001,
0.48399999999999993,
0.5235,
0.541
]

square_plus    = [
0.8145,
0.9,
0.881,
0.9365,
0.9279999999999998,
0.9225000000000001,
0.925,
0.9290000000000002,
0.864,
0.898
]

triangle_minus = [
0.3505,
0.34950000000000003,
0.3205,
0.41450000000000004,
0.3025,
0.264,
0.2245,
0.272,
0.23600000000000002,
0.356
]

triangle_plus  = [
0.5565,
0.329,
0.36349999999999993,
0.46599999999999997,
0.3435,
0.6030000000000001,
0.633,
0.35350000000000004,
0.6525,
0.5945
]


datasets = [
    [circle_minus, circle_plus],
    [square_minus, square_plus],
    [triangle_minus, triangle_plus]
]

# Colors: "-" darker, "+" brighter
colors = [
    ['darkred', 'red'],        # Circle
    ['darkgreen', 'limegreen'],      # Square
    ['brown', 'orange']            # Triangle
]

shapes = ['circle', 'square', 'triangle']
markers = ['o', 's', '^']
sub_labels = ['–', '+']
sub_offset = [-0.15, 0.15]

#fig, ax = plt.subplots(figsize=(7,5))
fig, ax = plt.subplots(figsize=(5.452423529/2,4.089317647/2))

for i, group in enumerate(datasets):
    for j, data in enumerate(group):
        xpos = i + sub_offset[j]
        
        # jittered scatter points
        x = np.random.normal(xpos, 0.04, size=len(data))
        ax.plot(x, data, markers[i], color=colors[i][j], alpha=0.7, markersize=5)

        # mean and SEM
        mean = np.mean(data)
        sem  = np.std(data, ddof=1) / np.sqrt(len(data))

        # horizontal line at mean
        ax.hlines(mean, xpos-0.1, xpos+0.1, color='k', linewidth=1.5)

        # error bar
        ax.errorbar(xpos, mean, yerr=sem, color='k', capsize=5, lw=1.5)

# Build x tick positions and labels
tick_positions = []
tick_labels = []
for i in range(len(datasets)):
    for j, lbl in enumerate(sub_labels):
        tick_positions.append(i + sub_offset[j])
        tick_labels.append(lbl)

ax.set_xticks(tick_positions)
ax.set_xticklabels(tick_labels, fontsize=18)
ax.set_yticks(np.arange(0, 1.01, 0.2))
ax.tick_params(axis='y', labelsize=18)

# Add shapes below the axis (outside the plotting area)
ymin, ymax = ax.get_ylim()
y_offset = ymin - (ymax - ymin) * 0.1  # 10% below axis

for i, shape in enumerate(shapes):
    xpos = i  # group center
    if shape == 'circle':
        ax.plot(xpos, y_offset, 'o', color='red', markersize=12, clip_on=False)
    elif shape == 'square':
        ax.plot(xpos, y_offset, 's', color='green', markersize=12, clip_on=False)
    elif shape == 'triangle':
        ax.plot(xpos, y_offset, '^', color='darkorange', markersize=12, clip_on=False)

# Keep the x-axis as standard
ax.set_ylim(ymin, ymax)  # restore normal axis limits
ax.set_ylabel('Junction persistence', fontsize=18)

plt.tight_layout()
plt.show()
exit(1)
'''


markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'x']
colors = plt.cm.tab10.colors  # 10 distinct colors


final_x = []
final_y = []
final_yy = []
omega = []
time_array = np.linspace(0, 1000000, 1000)
#plt.figure(figsize=(10, 6))
for traj in range(start, end):
#for traj in [start, end]:
    #if traj==start:
        #theta_5 = []
    #else:
        #theta_1 = []

    last_yy = []
    xx = []
    theta_5 = []
    theta_1 = []
    theta_2 = []
    theta_3 = []
    theta_6 = []
    theta_7 = []
    S = []
    phi_all = []
    for job in range(jobs_seq[traj], jobs_seq[traj+1]):

        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/order_parameters.txt","r")
        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/theta_shape.txt","r")
        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/elongation_shape.txt","r")
        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/elongation_minor_shape.txt","r")
        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/MSD.txt","r")
        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/time_rotation.txt","r")
        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/time_interface_rotation_tangent.txt","r")
        index_count = 0
        mean_orientation = 0.
        mean_phi = 0.
        sum_time = 0.
        for line in fileoutput:
            save=line.split()
            theta_5.append((2*(float(save[variable]))*pi/180))
            #sum_time += float(save[variable])
            #theta_5.append(float(save[variable]))

            #plt.plot(
                    #float(save[variable+1]), float(save[variable]),
                    #marker=markers[job % len(markers)],
                    #color=colors[job % len(colors)],
                    #linestyle=''
                    #)

            #if job == jobs_seq[traj]:
                #theta_5.append(float(save[variable]) / jobs[traj])
                #theta_1.append(float(save[variable+1]) / jobs[traj])
                #xx.append(float(save[0]))
            #else:
                #theta_5[index_count] += float(save[variable]) / jobs[traj]
                #theta_1[index_count] += float(save[variable+1]) / jobs[traj]
            #index_count+=1
            #theta_1.append(float(save[variable-1]))

        #for qq in range(len(theta_5)-50, len(theta_5)):
            #last_yy.append(theta_5[qq])


            '''
            num_rows = 0
            sin_sum = 0.
            cos_sum = 0.
            for q in range(3, len(save) - 3):
                mean_orientation += float(save[q])
                sin_sum += sin(2 * (float(save[q])*pi/180))
                cos_sum += cos(2 * (float(save[q])*pi/180))
                num_rows += 1

            mean_orientation = mean_orientation / num_rows
            mean_phi = 0.5 * atan2(sin_sum, cos_sum)
            phi_all.append(mean_phi * 180 / pi)
            avg_value_S = 0.
            for q in range(3, len(save) - 3):
                #avg_value_S += cos(2*(float(save[q]) - mean_orientation))
                avg_value_S += cos(2*( (float(save[q])*pi/180) - mean_phi))
            S.append(avg_value_S / num_rows)
            '''

            #if traj == start:
                #theta_5.append(float(save[variable]))
            #else:
                #theta_1.append(float(save[variable]))

        fileoutput.close()
        #last_yy.append(np.mean(theta_5[-50:]))


        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/elongation_shape.txt","r")
        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/aspect_ratio_shape.txt","r")
        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/time_jamm.txt","r")
        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/order_parameters.txt","r")
        for line in fileoutput:
            save=line.split()
            if job == jobs_seq[traj]:
                #theta_1.append(float(save[variable]) / jobs[traj])
                theta_1.append(float(save[4]) / jobs[traj])
                #theta_1.append(10*((float(save[4])-42)/20))
            else:
                theta_1[index_count] += float(save[variable]) / jobs[traj]
            index_count+=1
        fileoutput.close()


        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_21/elongation_shape.txt","r")
        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_21/order_parameters.txt","r")
        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_21/aspect_ratio_shape.txt","r")
        for line in fileoutput:
            save=line.split()
            theta_2.append(float(save[variable]))
        fileoutput.close()
        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_21/theta_shape.txt","r")
        for line in fileoutput:
            save=line.split()
            theta_3.append((2*(float(save[variable]))*pi/180))
        fileoutput.close()


        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_23/order_parameters.txt","r")
        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_23/elongation_shape.txt","r")
        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_23/aspect_ratio_shape.txt","r")
        for line in fileoutput:
            save=line.split()
            theta_6.append(float(save[variable]))
        fileoutput.close()
        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_23/theta_shape.txt","r")
        for line in fileoutput:
            save=line.split()
            theta_7.append((2*(float(save[variable]))*pi/180))
        fileoutput.close()


    #print("-------------------------")
    #omega.append(sum(theta_5)/(sum(theta_5)+sum(theta_1)))
    #omega.append(sum(theta_5)/jobs[traj])


    #last_yy = theta_5[-50:]
    #last_xx = xx[-50:]
    #slope, intercept, _, _, _ = linregress(last_xx, last_yy)
    #final_x.append((N[traj] * 8 * 8) / ((lx[traj] - 2 * 6)/2)**2)
    #final_y.append(slope)
    #final_x.append(shear_rate[traj])
    #final_y.append(np.mean(last_yy))
    #final_yy.append(stats.sem(last_yy))
    #final_yy.append(np.std(last_yy))

    #plt.plot(theta_5, '--o', label=shear_rate[traj])
    #plt.plot(theta_5, '--o', label=variable)
    #plt.plot(theta_1, '--s', label=variable-1)
    #plt.plot(S, '--o', label=traj)
    #plt.plot(phi_all, '--o', label=traj)
    #if traj==end:
    #theta_diff = [abs(theta_5[i] - theta_1[i]) for i in range(len(theta_5))]
    #plt.plot(theta_diff, '--o', label=traj)


    #period, method = estimate_period(time_array, np.array(theta_5))
    #print("Estimated period:", period, " (method:", method, ")")
    #print(period)
    #plt.clf()

    if end == start:
        break

#print(final_x)
#print(final_y)
#print(final_yy)

#final_x = [1.0, 0.75, 0.5, 0.25, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001]
#final_y = [0.49992082742607136, 0.48058329369018765, 0.46116523160021705, 0.4377948871250764, 0.4428061377749363, 0.5159620665712491, 0.523863398829233, 0.5828769623310098, 0.2893244788742896, 0.40499525754918736, 0.03450679770844628]
#final_yy = [np.float64(0.9509520336270171), np.float64(0.9029565600880146), np.float64(0.7558733701227552), np.float64(0.746363071735866), np.float64(0.4716647352837026), np.float64(0.6990953993090971), np.float64(0.8111007830280228), np.float64(0.8917666537386556), np.float64(0.5177900353375421), np.float64(0.24271438919906277), np.float64(0.040002335139040814)]

#final_x = np.array([0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0007, 0.0008, 0.0013, 0.0017, 0.0024, 0.0033, 0.0045, 0.0062, 0.0085, 0.0117, 0.0161, 0.0221, 0.0304, 0.0418, 0.0574, 0.0788, 0.1083, 0.1487, 0.2043, 0.2807, 0.3857, 0.5298, 0.7279, 1.0])
#final_y = np.array([np.float64(0.2896697218998192), np.float64(0.33630445246330626), np.float64(0.43400328779488384), np.float64(0.5311565138185511), np.float64(0.36390157063579237), np.float64(0.4825153593432849), np.float64(0.5456727494937862), np.float64(0.6290565442735762), np.float64(0.7141504126736063), np.float64(0.6562556393808907), np.float64(0.6472854206159414), np.float64(0.6003395495408252), np.float64(0.5664608817871121), np.float64(0.5448686791337419), np.float64(0.4660333152296762), np.float64(0.4173409962674875), np.float64(0.4488264016467213), np.float64(0.4981183325771869), np.float64(0.4900805944770621), np.float64(0.47759117672996), np.float64(0.4374513703054221), np.float64(0.460042828926506), np.float64(0.4700894656485616), np.float64(0.4471321150942393), np.float64(0.4672926325416104), np.float64(0.4742098447723692), np.float64(0.4724620589223277), np.float64(0.43883697137926064), np.float64(0.41935366245236233)])
#final_y_err_sem = np.array([np.float64(2.4779534323162614e-05), np.float64(3.408960752434948e-05), np.float64(0.0001535934146319442), np.float64(0.0010893234448635936), np.float64(0.00322176268479729), np.float64(0.0004006262118523729), np.float64(0.0001690090859186959), np.float64(0.001843101630475457), np.float64(0.0038740933990313157), np.float64(0.0038788802323136856), np.float64(0.01112146317191135), np.float64(0.005624149335539616), np.float64(0.006976612355473009), np.float64(0.01249706029438529), np.float64(0.013778587488239059), np.float64(0.014268167255011852), np.float64(0.01611260699501519), np.float64(0.01580071682927956), np.float64(0.011599906972495107), np.float64(0.009983951939427839), np.float64(0.010675663616198369), np.float64(0.01774433056979585), np.float64(0.007883718438772218), np.float64(0.011526042237994569), np.float64(0.011757790853711076), np.float64(0.009544023984836125), np.float64(0.012037905923731749), np.float64(0.012598008015292714), np.float64(0.007171705254389642)])
#final_y_err_std = np.array([np.float64(6.659588414731947e-05), np.float64(0.0021466228624031905), np.float64(0.0024526157987470928), np.float64(0.009710083120257984), np.float64(0.037068463729466955), np.float64(0.017306937546626322), np.float64(0.00628556824839815), np.float64(0.03333037203776919), np.float64(0.014079632541178892), np.float64(0.05281092285016849), np.float64(0.01290754822973717), np.float64(0.046402252195048815), np.float64(0.04431935018851303), np.float64(0.0440890124301868), np.float64(0.07066040716066056), np.float64(0.08275386080895042), np.float64(0.024882519261225266), np.float64(0.05708737575473875), np.float64(0.03525734632507903), np.float64(0.03343552763563954), np.float64(0.05082809682281099), np.float64(0.04285401805692968), np.float64(0.06610697814841611), np.float64(0.08009137715956625), np.float64(0.07903461795566741), np.float64(0.06348856133079687), np.float64(0.03258013474861658), np.float64(0.008039474890279357), np.float64(0.008740268641538517)])
#final_yy = [np.float64(0.1770427643629094), np.float64(0.19647260244224415), np.float64(0.2991794492672249), np.float64(0.4789462290131415), np.float64(0.4251314850451014), np.float64(0.5809837510369256), np.float64(0.6519316011323135), np.float64(0.7874967499823531), np.float64(0.8991038844949556), np.float64(0.8971982943468012), np.float64(0.9527027378747758), np.float64(0.8822926685862289), np.float64(0.8467586730150047), np.float64(0.8332334497422386), np.float64(0.769568120856856), np.float64(0.7706259080359582), np.float64(0.7932582112437027), np.float64(0.7610656248840395), np.float64(0.7450465472190528), np.float64(0.7097550125187485), np.float64(0.6940237228965223), np.float64(0.6189111241656535), np.float64(0.6157832941421842), np.float64(0.6747623316267313), np.float64(0.6362965311781632), np.float64(0.7547420485208346), np.float64(0.8219757103885442), np.float64(0.8704084163937043), np.float64(0.8979139352454599)]
#final_yy_err_sem = [np.float64(9.51369773533135e-06), np.float64(0.00030666040891474147), np.float64(0.000350373685535299), np.float64(0.0013871547314654262), np.float64(0.00529549481849528), np.float64(0.002472419649518046), np.float64(0.0008979383211997357), np.float64(0.0047614817196813125), np.float64(0.0020113760773112703), np.float64(0.007544417550024068), np.float64(0.001843935461391024), np.float64(0.0066288931707212595), np.float64(0.006331335741216147), np.float64(0.006298430347169542), np.float64(0.010094343880094366), np.float64(0.011821980115564344), np.float64(0.003554645608746467), np.float64(0.008155339393534107), np.float64(0.005036763760725575), np.float64(0.004776503947948506), np.float64(0.007261156688972999), np.float64(0.006122002579561383), np.float64(0.009443854021202302), np.float64(0.011441625308509463), np.float64(0.011290659707952486), np.float64(0.009069794475828126), np.float64(0.004654304964088084), np.float64(0.001148496412897051), np.float64(0.001248609805934074)]
#final_yy_err_std = np.array([np.float64(6.659588414731947e-05), np.float64(0.0021466228624031905), np.float64(0.0024526157987470928), np.float64(0.009710083120257984), np.float64(0.037068463729466955), np.float64(0.017306937546626322), np.float64(0.00628556824839815), np.float64(0.03333037203776919), np.float64(0.014079632541178892), np.float64(0.05281092285016849), np.float64(0.01290754822973717), np.float64(0.046402252195048815), np.float64(0.04431935018851303), np.float64(0.0440890124301868), np.float64(0.07066040716066056), np.float64(0.08275386080895042), np.float64(0.024882519261225266), np.float64(0.05708737575473875), np.float64(0.03525734632507903), np.float64(0.03343552763563954), np.float64(0.05082809682281099), np.float64(0.04285401805692968), np.float64(0.06610697814841611), np.float64(0.08009137715956625), np.float64(0.07903461795566741), np.float64(0.06348856133079687), np.float64(0.03258013474861658), np.float64(0.008039474890279357), np.float64(0.008740268641538517)])


#final_y = np.array([np.float64(0.2896697218998192), np.float64(0.3363044524633062), np.float64(0.43400328779488395), np.float64(0.5311565138185511), np.float64(0.36390157063579237), np.float64(0.4825153593432848), np.float64(0.5456727494937861), np.float64(0.629056544273576), np.float64(0.7141504126736063), np.float64(0.6562556393808908), np.float64(0.6472854206159414), np.float64(0.6003395495408252), np.float64(0.566460881787112), np.float64(0.5448686791337419), np.float64(0.46603331522967634), np.float64(0.41734099626748744), np.float64(0.44882640164672116), np.float64(0.4981183325771868), np.float64(0.4900805944770621), np.float64(0.47759117672996), np.float64(0.43745137030542214), np.float64(0.46004282892650616), np.float64(0.4700894656485616), np.float64(0.4471321150942392), np.float64(0.4672926325416104), np.float64(0.4742098447723692), np.float64(0.47246205892232773), np.float64(0.43883697137926064), np.float64(0.4193536624523622)])
#final_y_err_std = np.array([np.float64(0.14988819660637157), np.float64(0.12840999060439892), np.float64(0.16753510967171173), np.float64(0.22749771731515872), np.float64(0.18721095598679557), np.float64(0.22566023416979358), np.float64(0.13732808225565007), np.float64(0.19930413291475382), np.float64(0.11341902604432726), np.float64(0.0613008472313453), np.float64(0.04512467677189514), np.float64(0.038754018218931124), np.float64(0.045388992236570405), np.float64(0.05314927074891524), np.float64(0.03457109744525126), np.float64(0.06331653129309672), np.float64(0.07010703529098548), np.float64(0.04867610184204162), np.float64(0.0342553098757067), np.float64(0.032667163651612194), np.float64(0.019475273910254148), np.float64(0.030732971330311054), np.float64(0.03862256518604238), np.float64(0.04761534197680471), np.float64(0.020846434948343397), np.float64(0.044639250362093534), np.float64(0.021467925660586123), np.float64(0.031651481523000444), np.float64(0.0359232120502082)])

#final_yy = np.array([np.float64(0.1770427643629094), np.float64(0.19647260244224413), np.float64(0.29917944926722495), np.float64(0.4789462290131416), np.float64(0.42513148504510134), np.float64(0.5809837510369256), np.float64(0.6519316011323136), np.float64(0.7874967499823531), np.float64(0.8991038844949555), np.float64(0.8971982943468013), np.float64(0.9527027378747759), np.float64(0.882292668586229), np.float64(0.846758673015005), np.float64(0.8332334497422386), np.float64(0.7695681208568559), np.float64(0.770625908035958), np.float64(0.7932582112437029), np.float64(0.7610656248840396), np.float64(0.7450465472190528), np.float64(0.7097550125187485), np.float64(0.694023722896522), np.float64(0.6189111241656535), np.float64(0.6157832941421841), np.float64(0.6747623316267314), np.float64(0.6362965311781632), np.float64(0.7547420485208345), np.float64(0.821975710388544), np.float64(0.8704084163937044), np.float64(0.89791393524546)])
#final_yy_err_std = np.array([np.float64(0.10863810861394337), np.float64(0.09212537713171122), np.float64(0.11954129797131798), np.float64(0.17723006672867725), np.float64(0.12901985836815003), np.float64(0.13963142290636962), np.float64(0.09001271021915588), np.float64(0.05598211021176802), np.float64(0.08756261418090687), np.float64(0.04683217950636677), np.float64(0.03305875910738458), np.float64(0.06380806273706008), np.float64(0.07056540913058833), np.float64(0.035674676212209974), np.float64(0.07130928449810396), np.float64(0.08101935232552417), np.float64(0.0758855297325258), np.float64(0.11596248087355591), np.float64(0.08230799634161359), np.float64(0.0878132979496909), np.float64(0.07182992654035648), np.float64(0.1294069883158062), np.float64(0.08048978295632345), np.float64(0.06983955838366351), np.float64(0.0662632732697333), np.float64(0.046758489540577264), np.float64(0.06629134036082554), np.float64(0.030609695672286545), np.float64(0.028092017828046656)])


plt.plot(theta_5, '--o', color='firebrick')
plt.xlim([200,400])
plt.ylim([0,2])
#plt.plot(theta_1, '--s', color='forestgreen')
ax1 = plt.gca()
ax2 = ax1.twinx()
ax2.plot(theta_1, '--s', color='forestgreen')
#plt.plot(theta_6, '--^', color='royalblue')
ax2.set_xlim([200,400])
ax2.set_ylim([40,60])


'''
from scipy.optimize import curve_fit
x = np.linspace(0, 1000, 1000)
#y = 8*np.cos(0.075 * shear_rate[21] * x)**2 + 44.2
def f(x, k):
    return 8 * np.cos(k * shear_rate[start] * x)**2 + 44.2
popt, _ = curve_fit(f, x, theta_5, p0=0.1, maxfev=10000)
print("params:", popt)
plt.plot(x, f(x, *popt), color='red', label='Cosine fit')


def abs_cos_model(x, k, p):
    return 8 * (np.abs(np.cos(k * shear_rate[start] * x))**1) + 44.2
p0 = [0.1, 2]
popt, pcov = curve_fit(abs_cos_model, x, theta_5, p0=p0, maxfev=10000)
#print("params:", popt)
#plt.plot(x, abs_cos_model(x, *popt), 'r', lw=2)
'''


#plt.ylabel(r'$\theta$', fontsize=18)
#plt.xlabel('t', fontsize=18)

#plt.plot(omega, '--o')

'''
gamma = [0.0418, 0.0574, 0.0788, 0.1083, 0.1487, 0.2043, 0.2807, 0.3857, 0.5298, 0.7279, 1]
#T_period = [5561.116672, 5005.005005, 3925.494122, 2860.00286, 2002.002002, 1668.335002, 1112.223334, 834.1675008, 625.6256256, 455.000455, 333.6670003]
T_period = [7500, 5005.005005, 3925.494122, 2860.00286, 2002.002002, 1668.335002, 1112.223334, 834.1675008, 625.6256256, 455.000455, 333.6670003]
color_T = ["royalblue", "royalblue", "royalblue", "forestgreen", "royalblue", "royalblue", "royalblue", "royalblue", "royalblue", "royalblue", "firebrick"]
gamma = np.array(gamma)
T_period = np.array(T_period)

# Fit in log-log space
log_gamma = np.log10(gamma)
log_T = np.log10(T_period)

# Solve for intercept b (since slope = -1 is fixed)
b = np.mean(log_T + log_gamma)

# Best-fit line with fixed slope
x_fit = np.logspace(np.log10(gamma.min()), np.log10(gamma.max()), 400)
y_fit = 10**b / x_fit   # since slope = -1

fig, ax = plt.subplots(figsize=(7,5))

for i in range(len(gamma)):
    ax.plot(gamma[i], T_period[i], 'o', color=color_T[i], markersize=10)

ax.plot(x_fit, y_fit, '--', color='k', linewidth=1.5)

ax.set_ylabel('T', fontsize=18, fontname='Times New Roman')
ax.set_xlabel(r'$\dot{\gamma}$', fontsize=18, fontname='Times New Roman')
ax.tick_params(axis='both', labelsize=18)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontname('Times New Roman')
    label.set_fontsize(18)   # keep same size as tick_params
ax.set_xscale('log')
ax.set_yscale('log')

# Draw triangle (all black)
ax.plot([0.06, 0.075], [5570, 5570], 'k-', linewidth=1.2)
ax.plot([0.075, 0.075], [5570, 4420], 'k-', linewidth=1.2)
# Add slope label at the corner opposite to hypotenuse
ax.text(0.08, 5580, "-1", fontsize=16, ha='center', va='bottom', fontname='Times New Roman')

# ---- Add inset polar plots ----

#phase_space_plots
# Axes through center
size_r1 = max(theta_1)
size_r2 = max(theta_2)

# Compute velocity = sqrt((Δθ)^2 + (Δr)^2)
dtheta = np.diff(theta_5)
dr = np.diff(theta_1)
vel = np.sqrt(dtheta**2 + dr**2)
vel = vel / (dt[0]*deltat[0])

ax_inset1 = fig.add_axes([0.65, 0.605, 0.25, 0.25], polar=True)  # top-right
ax_inset2 = fig.add_axes([0.2, 0.275, 0.25, 0.25], polar=True)  # bottom-left

# Line in polar coordinates
ax_inset1.scatter(theta_5[0:500], theta_1[0:500], c=vel[0:500], cmap="RdBu_r", alpha=0.2, s=2)
#ax_inset1.scatter(theta_5[-500:], theta_1[-500:], color='firebrick', s=2)
sc1 = ax_inset1.scatter(theta_5[-500:], theta_1[-500:], c=vel[-500:], cmap="RdBu_r" , s=2)

# Colorbar for inset1
cbar1 = fig.colorbar(sc1, ax = ax_inset1, fraction=0.046, pad=0.2)
cbar1.set_label(r"$\omega$", labelpad=-25, y=1.2, rotation=0, fontsize=12, fontname="Times New Roman")
for ticklabel in cbar1.ax.get_yticklabels():   # or get_xticklabels() if horizontal
    ticklabel.set_fontsize(12)
    ticklabel.set_fontname("Times New Roman")


# Styling
ax_inset1.set_rmax(1.2 * size_r1)
ax_inset1.set_rticks([1.2*size_r1])  # Fewer radial ticks
for label in ax_inset1.get_yticklabels() + ax_inset1.get_xticklabels():
    label.set_fontname('Times New Roman')
    label.set_fontsize(10)
ax_inset1.set_thetagrids([0, 90, 180, 270], labels=[r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$"], fontsize=14, fontname='Times New Roman')
ax_inset1.grid(True)
ax_inset1.spines['polar'].set_color('firebrick')

# Compute velocity = sqrt((Δθ)^2 + (Δr)^2)
dtheta1 = np.diff(theta_3)
dr1 = np.diff(theta_2)
vel1 = np.sqrt(dtheta1**2 + dr1**2)
vel1 = vel1 / (dt[0]*deltat[0])


# Line in polar coordinates
ax_inset2.scatter(theta_3[0:500], theta_2[0:500], c=vel1[0:500], cmap="RdBu_r", alpha=0.2, s=2)
#ax_inset2.plot(theta_3[-500:], theta_2[-500:], color='forestgreen')
sc2 = ax_inset2.scatter(theta_3[-500:], theta_2[-500:], c=vel1[-500:], cmap="RdBu_r", s=2)

# Colorbar for inset2
cbar2 = fig.colorbar(sc2, ax=ax_inset2, fraction=0.046, pad=0.2)
cbar2.set_label(r"$\omega$", labelpad=-25, y=1.2, rotation=0, fontsize=12, fontname="Times New Roman")
for ticklabel in cbar2.ax.get_yticklabels():   # or get_xticklabels() if horizontal
    ticklabel.set_fontsize(12)
    ticklabel.set_fontname("Times New Roman")


# Styling
ax_inset2.set_rmax(1.2 * size_r2)
ax_inset2.set_rticks([1.2*size_r2])  # Fewer radial ticks
for label in ax_inset2.get_yticklabels() + ax_inset2.get_xticklabels():
    label.set_fontname('Times New Roman')
    label.set_fontsize(10)
ax_inset2.set_thetagrids([0, 90, 180, 270], labels=[r"$0$", r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$"], fontsize=14, fontname='Times New Roman')
#ax_inset2.set_yticklabels(['',''])
ax_inset2.grid(True)
ax_inset2.spines['polar'].set_color('forestgreen')

plt.tight_layout()
'''


'''
#psi6 and psi2
plt.plot(final_x, final_y, '--o', color='firebrick')
plt.fill_between(final_x, final_y - final_y_err_std, final_y + final_y_err_std, color="firebrick", alpha=0.1)
plt.ylabel(r'$\psi_6$', fontsize=18)
plt.xlabel(r'$\dot{\gamma}$', fontsize=18)
#fig = plt.gcf()
ax1 = plt.gca()
ax2 = ax1.twinx()
ax2.plot(final_x, final_yy, '--o', color='forestgreen')
ax2.fill_between(final_x, final_yy - final_yy_err_std, final_yy + final_yy_err_std, color="forestgreen", alpha=0.1)
ax2.set_ylabel(r'$\psi^L_2$', fontsize=18)
#ax2.plot(theta_1, '--o', color='forestgreen')
#ax2.set_ylabel('r', fontsize=18)
ax1.tick_params(axis='y', colors='firebrick')
ax2.tick_params(axis='y', colors='forestgreen')
ax1.yaxis.label.set_color('firebrick')
ax2.yaxis.label.set_color('forestgreen')
ax2.spines['right'].set_color('forestgreen')
ax2.spines['left'].set_color('firebrick')
ax2.yaxis.set_major_locator(MaxNLocator(nbins=5))
for tick in ax2.yaxis.get_ticklabels():
    tick.set_fontsize(18)
    tick.set_fontname('Times New Roman')
for tick in ax1.yaxis.get_ticklabels():
    tick.set_fontsize(18)
    tick.set_fontname('Times New Roman')
for tick in ax1.xaxis.get_ticklabels():
    tick.set_fontsize(18)
    tick.set_fontname('Times New Roman')
plt.xscale('log')
plt.tight_layout()
'''


#plt.ylabel(r'$\theta_i$', fontsize=18)
#plt.ylabel('MSD', fontsize=18)
#plt.ylabel('S', fontsize=18)
#plt.xlabel('Time', fontsize=18)
#plt.ylabel('D', fontsize=18)
#plt.xlabel(r'$\phi$', fontsize=18)
#plt.xlim([0,500])
#plt.ylim([0,3])
#plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
#plt.subplots_adjust(left=0.15, bottom=0.15, right=0.85, top=0.85)
#plt.legend(ncols=1, frameon=False, loc='upper left')
#plt.xscale('log')
#plt.yscale('log')
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Shear/new_psi6_psiN_shear.png", transparent=True)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Shear/new_psi6_psiN_shear.svg", transparent=True)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Shear/period_shear.png", transparent=True)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Shear/period_shear.svg", transparent=True)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Shear/lattice_phi_time_series_fit.png", transparent=True)
plt.show()
exit(1)



#Sumesh two tensor models----------------------------------------------------------------------------------------------
'''
valsY1 = [0. for i in range(end)]
valsY2 = [0. for i in range(end)]
valsY3 = [0. for i in range(end)]
valsX = [0. for i in range(end)]

jobs[27] -= 1

plt.figure(figsize=(5.452423529,4.089317647))
'''
#for traj in range(start, end, 2):
#for traj in [5,9]:
    #pdf_values = []

    #size_R = int(lx[traj]/2)
    #corr_R = [i for i in range(size_R)]

    #corr_Q = [0. for i in range(size_R)]

'''
    for job in range(jobs_seq[traj], jobs_seq[traj+1]):
        if job == 83:
            continue

        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/defects_velocity_nematic.txt","r")
        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/defects_velocity_shape.txt","r")
        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/histogram_QS.txt","r")
        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/histogram_SV.txt","r")

        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/histogram_QV.txt","r")
        #for line in fileoutput:
            #save=line.split()
            #pdf_values.append(float(save[0]))
        #fileoutput.close()

        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/histogram_QSV_stats.txt","r")
        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/defects_stats.txt","r")
        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/correlations_nematic.txt","r")
        #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/correlations_shape.txt","r")
        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/correlations_velocity.txt","r")
        for line in fileoutput:
            save=line.split()
            #valsY1[traj] += float(save[0]) / jobs[traj]
            #valsY2[traj] += float(save[1]) / jobs[traj]
            #valsY3[traj] += float(save[2]) / jobs[traj]
            #valsX[traj] = J_Q[traj]+0.0000001

            if int(float(save[0])) < size_R:
                corr_Q[int(float(save[0]))] += float(save[2]) / jobs[traj]
        fileoutput.close()


    plt.plot(corr_R, corr_Q, '--o', label=J_Q[traj])
'''
    

'''
    counts, bins = np.histogram(pdf_values)
    bin_width = abs(bins[1] - bins[0]) / 2
    bin_length = len(bins)
    total_counts = sum(counts)
    probability = []
    for i in range(len(counts)):
        probability.append(float(counts[i]) / float(total_counts))
        bins[i] = (bins[i] + bin_width) * 180 / pi
        #bins[i] = (bins[i] + bin_width)
    #plt.stairs(counts, bins)
    plt.plot(bins[0:bin_length-1], probability, '--o', label=J_Q[traj])
    #plt.plot(bins[0:bin_length-1], probability, '--o', label=gamma[traj])
    #plt.plot(bins[0:bin_length-1], probability, '--o', label=omega[traj])
'''


'''
plt.ylabel(r'$PDF$', fontsize=18)
#plt.xlabel(r'Velocity $+1/2$ Shape', fontsize=18)
plt.xlabel(r'$\theta_{QV}$', fontsize=18)
plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
plt.legend(loc=(0.4, 0.7), ncols=1, frameon=False)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/mean_velocity_1.png", transparent=True)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/ResultsSumesh/PDF_theta_QV_gamma006_act05_fric001_JQ.png", transparent=True)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/ResultsSumesh/PDF_theta_QV_gamma006_act05_fric001_JQ.svg", transparent=True)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/ResultsSumesh/PDF_velocity_defS_gamma006_act05_fric001_JQ.png", transparent=True)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/ResultsSumesh/PDF_velocity_defS_gamma006_act05_fric001_JQ.svg", transparent=True)
plt.show()
'''

'''
save = valsX[17]
valsX[17] = valsX[19]
valsX[19] = save
save = valsY1[17]
valsY1[17] = valsY1[19]
valsY1[19] = save
save = valsY2[17]
valsY2[17] = valsY2[19]
valsY2[19] = save
save = valsY3[17]
valsY3[17] = valsY3[19]
valsY3[19] = save
'''
#plt.plot(valsX[0:10:2], valsY1[0:10:2], "--o", label=r"$\zeta=0.3, \xi=0.1$")
#plt.plot(valsX[1:10:2], valsY1[1:10:2], "--o", label=r"$\zeta=0.5, \xi=0.1$")
#plt.plot(valsX[10:20:2], valsY1[10:20:2], "--o", label=r"$\zeta=0.3, \xi=0.01$")
#plt.plot(valsX[11:20:2], valsY1[11:20:2], "--o", label=r"$\zeta=0.5, \xi=0.01$")
#plt.plot(valsX[20:24:2], valsY1[20:24:2], '--o', label=r'$\zeta=0.3, \xi=0.1$')
#plt.plot(valsX[21:24:2], valsY1[21:24:2], '--o', label=r'$\zeta=0.5, \xi=0.1$')
#plt.plot(valsX[24:30:2], valsY1[24:30:2], '--o', label=r'$\zeta=0.3, \xi=0.01$')
#plt.plot(valsX[25:30:2], valsY1[25:30:2], '--o', label=r'$\zeta=0.5, \xi=0.01$')
#plt.xscale("log")

#plt.plot(valsX[0:10:2], valsY2[0:10:2], '--o', label=r'$\zeta=0.3, \xi=0.1$')
#plt.plot(valsX[1:10:2], valsY2[1:10:2], '--o', label=r'$\zeta=0.5, \xi=0.1$')
#plt.plot(valsX[10:20:2], valsY2[10:20:2], '--o', label=r'$\zeta=0.3, \xi=0.01$')
#plt.plot(valsX[11:20:2], valsY2[11:20:2], '--o', label=r'$\zeta=0.5, \xi=0.01$')
#plt.plot(valsX[20:24:2], valsY2[20:24:2], '--o', label=r'$\zeta=0.3, \xi=0.1$')
#plt.plot(valsX[21:24:2], valsY2[21:24:2], '--o', label=r'$\zeta=0.5, \xi=0.1$')
#plt.plot(valsX[24:30:2], valsY2[24:30:2], '--o', label=r'$\zeta=0.3, \xi=0.01$')
#plt.plot(valsX[25:30:2], valsY2[25:30:2], '--o', label=r'$\zeta=0.5, \xi=0.01$')
#plt.xscale('log')

#plt.plot(valsX[0:10:2], valsY3[0:10:2], '--o', label=r'$\zeta=0.3, \xi=0.1$')
#plt.plot(valsX[1:10:2], valsY3[1:10:2], '--o', label=r'$\zeta=0.5, \xi=0.1$')
#plt.plot(valsX[10:20:2], valsY3[10:20:2], '--o', label=r'$\zeta=0.3, \xi=0.01$')
#plt.plot(valsX[11:20:2], valsY3[11:20:2], '--o', label=r'$\zeta=0.5, \xi=0.01$')
#plt.plot(valsX[20:24:2], valsY3[20:24:2], '--o', label=r'$\zeta=0.3, \xi=0.1$')
#plt.plot(valsX[21:24:2], valsY3[21:24:2], '--o', label=r'$\zeta=0.5, \xi=0.1$')
#plt.plot(valsX[24:30:2], valsY3[24:30:2], '--o', label=r'$\zeta=0.3, \xi=0.01$')
#plt.plot(valsX[25:30:2], valsY3[25:30:2], '--o', label=r'$\zeta=0.5, \xi=0.01$')
#plt.xscale('log')


#plt.ylabel('Average velocity', fontsize=18)
#plt.ylabel('Shape AR', fontsize=18)
#plt.ylabel('misaligned area', fontsize=18)
#plt.ylabel('Defect number', fontsize=18)
#plt.xlabel(r'$J_{QS}$', fontsize=18)
#plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
#plt.legend(loc=(0.01, 0.6), ncols=1, fontsize=12, frameon=False)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/ResultsSumesh/Defect_number_gamma006_JQ.png", transparent=True)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/ResultsSumesh/Defect_number_gamma006_JQ.svg", transparent=True)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/ResultsSumesh/Misaligned_area_gamma006_JQ.png", transparent=True)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/ResultsSumesh/Misaligned_area_gamma006_JQ.svg", transparent=True)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/ResultsSumesh/Shape_AR_gamma006_JQ.png", transparent=True)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/ResultsSumesh/Shape_AR_gamma006_JQ.svg", transparent=True)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/ResultsSumesh/Average_V_gamma006_JQ.png", transparent=True)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/ResultsSumesh/Average_V_gamma006_JQ.svg", transparent=True)

'''
plt.ylabel(r'$C_W$', fontsize=18)
plt.xlabel(r'$R$', fontsize=18)
plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
plt.legend(loc=(0.6, 0.5), ncols=1, fontsize=12, frameon=False)
plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/ResultsSumesh/Corr_W_gamma006_act05_fric001_JQ.png", transparent=True)
plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/ResultsSumesh/Corr_W_gamma006_act05_fric001_JQ.svg", transparent=True)

plt.show()
exit (1)
'''
#Sumesh two tensor models----------------------------------------------------------------------------------------------







#MPF holes paper----------------------------------------------------------------------------------------------
# read output file
n_holes = np.zeros(end-start)
lifetime_holes = np.zeros(end-start)
max_area_holes = np.zeros(end-start)
avg_FTLE = np.zeros(end-start)
weibull_sigma_0 = []


max_stress_all = []
vortex_all = []
survival_all_y = []
survival_all_x = []
logsurvival_all_y = []
logsurvival_all_x = []
nematic_1_all = []
nematic_2_all = []

st_all_x = []
st_all_y = []

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
    avg_distance_nematic_defects_1 = np.zeros(10)
    avg_distance_nematic_defects_2 = np.zeros(10)

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

        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_stats.txt","r")
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


        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_histogram.txt","r")
        for line in fileoutput:
            save=line.split()
            hist_step=int(int(float(save[0]) * 50.01) * 0.5)
            area_histogram[hist_step] += (float(save[1]) * pi * 8 * 8) / (lx[traj] * ly[traj]) * 0.5
            circularity_histogram[hist_step] += float(save[2]) * 0.5
            radius_speed_histogram[hist_step] += float(save[3]) * 0.5
            ani_histogram[hist_step] += float(save[4]) * 0.5
        fileoutput.close()

        '''
        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_histogram_circularity.txt","r")
        for line in fileoutput:
            save=line.split()
            circularity_histogram[int(float(save[0]) * 50.01)] += float(save[1])
        fileoutput.close()

        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_histogram_radius_speed_time.txt","r")
        for line in fileoutput:
            save=line.split()
            radius_speed_histogram[int(float(save[0]) * 50.01)] += float(save[1])
        fileoutput.close()

        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_histogram_ani.txt","r")
        for line in fileoutput:
            save=line.split()
            ani_histogram[int(float(save[0]) * 50.01)] += float(save[1])
        fileoutput.close()
        '''


        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_histogram_radius_speed.txt","r")
        for line in fileoutput:
            save=line.split()
            rspeed_histogram[int(float(save[0]))] += float(save[1])
        fileoutput.close()

        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_stress_histogram_tau10.txt","r")
        for line in fileoutput:
            save=line.split()
            avg_stress_histogram[int(float(save[0]))] += float(save[1])
            max_stress_histogram[int(float(save[0]))] += float(save[2])
        fileoutput.close()

        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_velocity_histogram_tau10.txt","r")
        for line in fileoutput:
            save=line.split()
            avg_distance_velocity_defects[int(float(save[0]))] += float(save[1])
        fileoutput.close()




        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_nematic_histogram_tau10.txt","r")
        for line in fileoutput:
            save=line.split()
            avg_distance_nematic_defects_1[int(float(save[0]))] += float(save[1])
            avg_distance_nematic_defects_2[int(float(save[0]))] += float(save[2])
        fileoutput.close()




        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_velocity_stats.txt","r")
        for line in fileoutput:
            save=line.split()
            avg_FTLE[traj - start] += float(save[2]) / jobs[traj]
        fileoutput.close()


        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_strain_histogram_20bins.txt","r")
        cont_line = 0
        for line in fileoutput:
            save=line.split()
            x_st_20[cont_line] = float(save[0])
            strain_hist_20[cont_line] += float(save[1])
            strain_count_20[cont_line] += float(save[2])
            cont_line+=1
        fileoutput.close()

        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_strain_histogram_30bins.txt","r")
        cont_line = 0
        for line in fileoutput:
            save=line.split()
            x_st_30[cont_line] = float(save[0])
            strain_hist_30[cont_line] += float(save[1])
            strain_count_30[cont_line] += float(save[2])
            cont_line+=1
        fileoutput.close()

        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_strain_histogram_40bins.txt","r")
        cont_line = 0
        for line in fileoutput:
            save=line.split()
            x_st_40[cont_line] = float(save[0])
            strain_hist_40[cont_line] += float(save[1])
            strain_count_40[cont_line] += float(save[2])
            cont_line+=1
        fileoutput.close()

        '''
        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/voids_strain_defectsplusone.txt","r")
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
        '''


    if count_jobs > 0:
        area_histogram[:] = area_histogram[:] / count_jobs
        circularity_histogram[:] = circularity_histogram[:] / count_jobs
        radius_speed_histogram[:] = radius_speed_histogram[:] / count_jobs
        rspeed_histogram[:] = rspeed_histogram[:] / count_jobs
        ani_histogram[:] = ani_histogram[:] / count_jobs

        avg_stress_histogram[:] = avg_stress_histogram[:] / count_jobs
        max_stress_histogram[:] = max_stress_histogram[:] / count_jobs

        avg_distance_velocity_defects[:] = avg_distance_velocity_defects[:] / count_jobs
        avg_distance_nematic_defects_1[:] = avg_distance_nematic_defects_1[:] / count_jobs
        avg_distance_nematic_defects_2[:] = avg_distance_nematic_defects_2[:] / count_jobs

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

        '''
        plt.close()
        plt.figure(figsize=(3*40.179/25.4, 3*25.142/25.4))
        plt.plot(x_lt, ani_histogram, '--o', color='firebrick')
        #ax2.plot(x_lt, area_histogram, '--s', color='green')
        #plt.ylabel(r'$A_{max}/L^2$', fontname='Times New Roman', fontsize=18)
        plt.ylabel('Anisotropy', fontname='Times New Roman', fontsize=18)
        plt.xlabel(r'$\tau_{hole}$', fontname='Times New Roman', fontsize=18)
        plt.xticks(fontname='Times New Roman', fontsize=18)
        plt.yticks(fontname='Times New Roman', fontsize=18)
        plt.subplots_adjust(left=0.21, bottom=0.22, right=0.80, top=0.82)
        fig = plt.gcf()
        ax = plt.gca()
        ax2 = ax.twinx()
        ax2.plot(x_lt, area_histogram, '--s', color='green')
        ax2.set_ylabel(r'$A_{max}/L^2$', fontname='Times New Roman', fontsize=18)
        ax.tick_params(axis='y', colors='firebrick')
        ax2.tick_params(axis='y', colors='green')
        ax.yaxis.label.set_color('firebrick')
        ax2.yaxis.label.set_color('green')
        ax2.spines['right'].set_color('green')
        ax2.spines['left'].set_color('firebrick')
        for tick in ax2.yaxis.get_ticklabels():
            tick.set_fontsize(18)
            tick.set_fontname('Times New Roman')
        plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/figOther.png", transparent=True)
        plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/figOther.svg", transparent=True)
        plt.show()
        exit (1)
        '''


    if variable == 3:
        if total_hist_counts>0:
            survival_probability_hist[:] = survival_probability_hist[:] / total_hist_counts

        ax[0,0].plot(avg_stress_histogram, '-o')
        ax[0,1].plot(max_stress_histogram, '-s')
        max_stress_all.append(max_stress_histogram)
        ax[1,0].plot(avg_distance_velocity_defects, '-^', label=nemself[traj])
        vortex_all.append(avg_distance_velocity_defects)


        nematic_1_all.append(avg_distance_nematic_defects_1)
        nematic_2_all.append(avg_distance_nematic_defects_2)


        xx = x_sp
        yy = np.zeros(len(x_sp))
        Ntotal = sum(survival_probability_counts[:])
        final_bin = -1
        if Ntotal > 0:
            for q in range(len(x_sp)):
                yy[q] = (Ntotal - sum(survival_probability_counts[0:q])) / Ntotal
                if yy[q] < 0.065 and final_bin == -1:
                    final_bin = q

        #ax[1,2].plot(xx, yy, '-o')
        st_all_x.append(xx)
        st_all_y.append(yy)
        final_bin=-3
        if final_bin > -2:
            S_sigma = yy[1:final_bin]
            S_sigma[S_sigma <= 1e-10] = 1e-10
            #lnlnS = np.log(-np.log(S_sigma))
            lnlnS = np.log(S_sigma)
            lnsigma = xx[1:final_bin]
            lnsigma = np.array(lnsigma)
            #lnsigma = np.log(lnsigma)
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
        survival_all_y.append(yy)
        survival_all_x.append(xx)
        final_bin = -3
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
            weibull_sigma_0.append(sigma_0)

            #for q in range(len(yy)):
                #if yy[q] > 0:
                    #yy[q] = -log(yy[q])

            #ax[1,1].plot(xx, yy, '-v')
            #ax[1,1].plot(lnsigma, lnlnS, '-v')
            ax[1,1].plot(sigma, minuslnS, '-v')
            logsurvival_all_y.append(minuslnS)
            logsurvival_all_x.append(sigma)


            line_plot = []
            if traj == end - 1:
                for q in range(len(lnsigma)):
                    line_plot.append(np.exp(2 * lnsigma[q] + intercept - 3))
                ax[1,1].plot(sigma, line_plot, '-')
                logsurvival_all_y.append(line_plot)
                logsurvival_all_x.append(sigma)


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

    #ax[1,1].set_yscale('log')
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
    #ax[0,0].plot(nemself[0:5], lifetime_holes[0:5], '--p')
    #ax[0,0].plot(nemself[5:10], lifetime_holes[5:10], '--p')
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

    #ax[0,1].plot(nemself[0:5], max_area_holes[0:5], '--v')
    #ax[0,1].plot(nemself[5:10], max_area_holes[5:10], '--v')
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

    #ax[0,2].plot(nemself[0:5], n_holes[0:5], '--^')
    #ax[0,2].plot(nemself[5:10], n_holes[5:10], '--^')
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
    #plt.show()
    plt.close()

    '''
    plt.figure(figsize=(3*40.179/25.4, 3*25.142/25.4))
    #plt.plot(gamma[4:8], max_area_holes[0:4], '--o', color='firebrick')
    plt.plot(nemself[0:5], max_area_holes[0:5], '--p', color='firebrick')
    plt.plot(nemself[5:10], max_area_holes[5:10], '--h', color='green')
    plt.plot(nemself[10:15], max_area_holes[10:15], '--D', color='royalblue')
    #plt.legend(['0.1', '0.01', '0.001'], prop={'family':'Times New Roman', 'size':'12'}, loc=(0.05, 0.5), ncols=1, frameon=False)
    #plt.text(0.145, 0.185, r'$\\chi$', fontsize=18, fontname="Times New Roman")
    plt.ylabel(r'$A_{max}/L^2$', fontname='Times New Roman', fontsize=18)
    plt.xlabel(r'$\\zeta$', fontname='Times New Roman', fontsize=18)
    #plt.xlabel(r'$\\gamma$', fontname='Times New Roman', fontsize=18)
    plt.xticks(fontname='Times New Roman', fontsize=18)
    plt.yticks(fontname='Times New Roman', fontsize=18)
    plt.subplots_adjust(left=0.21, bottom=0.225, right=0.985, top=0.995)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/fig3a.png", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/fig3a.svg", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/figSM2a.png", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/figSM2a.svg", transparent=True)
    plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig3a.png", transparent=True)
    plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig3a.svg", transparent=True)
    plt.show()
    '''

    '''
    plt.figure(figsize=(3*40.179/25.4, 3*25.142/25.4))
    #plt.plot(gamma[4:8], max_area_holes[0:4], '--o', color='firebrick')
    plt.plot(nemself[0:5], lifetime_holes[0:5], '--p', color='firebrick')
    plt.plot(nemself[5:10], lifetime_holes[5:10], '--h', color='green')
    plt.plot(nemself[10:15], lifetime_holes[10:15], '--D', color='royalblue')
    #plt.legend(['0.1', '0.01', '0.001'], prop={'family':'Times New Roman', 'size':'12'}, loc=(0.05, 0.5), ncols=1, frameon=False)
    #plt.text(0.145, 0.185, r'$\\chi$', fontsize=18, fontname="Times New Roman")
    plt.ylabel(r'$\tau_{hole}/t_{total}$', fontname='Times New Roman', fontsize=18)
    plt.xlabel(r'$\\zeta$', fontname='Times New Roman', fontsize=18)
    #plt.xlabel(r'$\\gamma$', fontname='Times New Roman', fontsize=18)
    plt.xticks(fontname='Times New Roman', fontsize=18)
    plt.yticks(fontname='Times New Roman', fontsize=18)
    plt.subplots_adjust(left=0.21, bottom=0.225, right=0.985, top=0.995)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/fig3a.png", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/fig3a.svg", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/figSM2a.png", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/figSM2a.svg", transparent=True)
    plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig3New.png", transparent=True)
    plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig3New.svg", transparent=True)
    plt.show()
    '''

    '''
    plt.figure(figsize=(3*40.179/25.4, 3*25.142/25.4))
    plt.plot(nemself[0:5], n_holes[0:5], '--p', color='firebrick')
    plt.plot(nemself[5:10], n_holes[5:10], '--h', color='green')
    plt.plot(nemself[10:15], n_holes[10:15], '--D', color='royalblue')
    #plt.plot(gamma[4:8], n_holes[0:4], '--o', color='firebrick')
    #plt.legend(['0.1', '0.01', '0.001'], prop={'family':'Times New Roman', 'size':'12'}, loc=(0.1, 0.6), ncols=1, frameon=False)
    plt.ylabel(r'$N_{holes}$', fontname='Times New Roman', fontsize=18)
    plt.xlabel(r'$\\zeta$', fontname='Times New Roman', fontsize=18)
    #plt.xlabel(r'$\\gamma$', fontname='Times New Roman', fontsize=18)
    plt.xticks(fontname='Times New Roman', fontsize=18)
    plt.yticks(fontname='Times New Roman', fontsize=18)
    plt.ylim(-0.1,2.7)
    plt.subplots_adjust(left=0.21, bottom=0.225, right=0.985, top=0.995)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/fig3b.png", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/fig3b.svg", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/figSM2b.png", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/figSM2b.svg", transparent=True)
    plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig3b.png", transparent=True)
    plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig3b.svg", transparent=True)
    plt.show()
    '''


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
    #plt.show()
    plt.close()


    #plt.figure(figsize=(5.452423529,4.089317647))
    #plt.figure(figsize=(3*40.179/25.4, 3*25.142/25.4))
    #y_values = [weibull_sigma_0[i] for i in range(len(weibull_sigma_0))] 
    #x_values = [sqrt(1/nemself[i]) for i in range(start, end)] 
    #x_values = [nemself[i] for i in range(start, end)] 
    #plt.plot(x_values[0:5], y_values[0:5], '--o')
    #plt.plot(x_values[5:10], y_values[5:10], '--s')

    #y_values = [weibull_sigma_0[i] for i in range(len(weibull_sigma_0))] 
    #x_values = [fric[i] for i in range(start, end)] 
    #plt.plot(x_values, y_values, '--o')
    #plt.xscale('log')

    '''
    ttt = [-i for i in range(9, -1, -1)]
    plt.figure(figsize=(3*40.179/25.4, 3*25.142/25.4))
    plt.plot(ttt, max_stress_all[11][::-1], '--o', color='firebrick')
    plt.plot(ttt, max_stress_all[10][::-1], '--^', color='green')
    plt.plot(ttt, max_stress_all[9][::-1], '--s', color='royalblue')
    plt.plot(ttt, max_stress_all[8][::-1], '--p', color='goldenrod')
    plt.plot(ttt, max_stress_all[7][::-1], '--h', color='peru')
    plt.plot(ttt, max_stress_all[5][::-1], '--8', color='darkviolet')
    plt.plot(ttt, max_stress_all[3][::-1], '--D', color='gray')
    plt.legend(['0.001', '0.0025', '0.005', '0.0075', '0.01', '0.05', '0.1'], prop={'family':'Times New Roman', 'size':'12'}, loc=(0.08, 0.05), ncols=3, frameon=False)
    #plt.text(-9.2, 0.185, r'$\\chi$', fontsize=18, fontname="Times New Roman")
    plt.ylabel(r"$\\dot{\varepsilon}_{hole}/\\dot{\varepsilon}_{max}$", fontname='Times New Roman', fontsize=18)
    plt.xlabel(r'$t-t_{hole}$', fontname='Times New Roman', fontsize=18)
    plt.xticks(fontname='Times New Roman', fontsize=18)
    plt.yticks(fontname='Times New Roman', fontsize=18)
    plt.subplots_adjust(left=0.2, bottom=0.22, right=0.975, top=0.99)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/fig1d.png", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/fig1d.svg", transparent=True)
    plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig1c.png", transparent=True)
    plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig1c.svg", transparent=True)
    plt.show()
    '''
    #plt.text(-9.35, 0.16, r"$\frac{\xi}{\xi_{cell}}$", fontsize=18, fontname="Times New Roman")


    '''
    ttt = [-i for i in range(9, -1, -1)]
    plt.figure(figsize=(3*40.179/25.4, 3*25.142/25.4))
    plt.plot(ttt, vortex_all[11][::-1]/8, '--o', color='firebrick')
    plt.plot(ttt, vortex_all[10][::-1]/8, '--^', color='green')
    plt.plot(ttt, vortex_all[9][::-1]/8, '--s', color='royalblue')
    plt.plot(ttt, vortex_all[8][::-1]/8, '--p', color='goldenrod')
    plt.plot(ttt, vortex_all[7][::-1]/8, '--h', color='peru')
    plt.plot(ttt, vortex_all[5][::-1]/8, '--8', color='purple')
    plt.plot(ttt, vortex_all[3][::-1]/8, '--D', color='gray')
    #plt.legend(['0.001', '0.0025', '0.005', '0.0075', '0.01', '0.05', '0.1'], prop=legend_font, loc=(0.02, 0.05), ncols=3, frameon=False)
    plt.ylabel(r'$r_{min}^{(spiral \\ core)}/R$', fontname='Times New Roman', fontsize=18)
    plt.xlabel(r'$t-t_{hole}$', fontname='Times New Roman', fontsize=18)
    plt.xticks(fontname='Times New Roman', fontsize=18)
    plt.yticks(fontname='Times New Roman', fontsize=18)
    plt.subplots_adjust(left=0.2, bottom=0.22, right=0.975, top=0.99)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/fig1c.png", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/fig1c.svg", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig1d.png", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig1d.svg", transparent=True)
    plt.show()
    '''


    ttt = [-i for i in range(9, -1, -1)]
    plt.figure(figsize=(3*40.179/25.4, 3*25.142/25.4))
    plt.plot(ttt, nematic_2_all[11][::-1]/8, '--o', color='firebrick')
    plt.plot(ttt, nematic_2_all[10][::-1]/8, '--^', color='green')
    plt.plot(ttt, nematic_2_all[9][::-1]/8, '--s', color='royalblue')
    plt.plot(ttt, nematic_2_all[8][::-1]/8, '--p', color='goldenrod')
    plt.plot(ttt, nematic_2_all[7][::-1]/8, '--h', color='peru')
    plt.plot(ttt, nematic_2_all[5][::-1]/8, '--8', color='purple')
    plt.plot(ttt, nematic_2_all[3][::-1]/8, '--D', color='gray')
    plt.legend(['0.001', '0.0025', '0.005', '0.0075', '0.01', '0.05', '0.1'], prop={'family':'Times New Roman', 'size':'12'}, loc=(0.08, 0.05), ncols=3, frameon=False)
    plt.text(-9.2, 1.2, r'$\chi$', fontsize=18, fontname="Times New Roman")
    plt.ylabel(r'$r_{min}^{(2^{nd} defect)}/R$', fontname='Times New Roman', fontsize=18)
    plt.xlabel(r'$t-t_{hole}$', fontname='Times New Roman', fontsize=18)
    plt.xticks(fontname='Times New Roman', fontsize=18)
    plt.yticks(fontname='Times New Roman', fontsize=18)
    plt.subplots_adjust(left=0.2, bottom=0.22, right=0.975, top=0.99)
    plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/figSM3b.png", transparent=True)
    plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/figSM3b.svg", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig1d.png", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig1d.svg", transparent=True)
    plt.show()



    '''
    plt.figure(figsize=(5.452423529,4.089317647))
    #plt.figure(figsize=(3*40.179/25.4, 3*25.142/25.4))
    #plt.plot(survival_all_x[2], survival_all_y[2], ':o', color='firebrick', markerfacecolor='none')
    #plt.plot(survival_all_x[3], survival_all_y[3], ':^', color='green', markerfacecolor='none')
    #plt.plot(survival_all_x[4], survival_all_y[4], ':s', color='royalblue', markerfacecolor='none')
    #plt.plot(survival_all_x[5], survival_all_y[5], '--p', color='goldenrod', label='0.1')
    #plt.plot(survival_all_x[6], survival_all_y[6], '--h', color='peru', label='0.2')
    #plt.plot(survival_all_x[7], survival_all_y[7], '--o', color='firebrick', label='0.3')
    #plt.plot(survival_all_x[8], survival_all_y[8], '--^', color='green', label='0.4')
    #plt.plot(survival_all_x[9], survival_all_y[9], '--s', color='royalblue', label='0.5')
    plt.plot(logsurvival_all_x[2], logsurvival_all_y[2], ':o', color='firebrick', markerfacecolor='none')
    plt.plot(logsurvival_all_x[3], logsurvival_all_y[3], ':^', color='green', markerfacecolor='none')
    plt.plot(logsurvival_all_x[4], logsurvival_all_y[4], ':s', color='royalblue', markerfacecolor='none')
    plt.plot(logsurvival_all_x[5], logsurvival_all_y[5], '--p', color='goldenrod', label='0.1')
    plt.plot(logsurvival_all_x[6], logsurvival_all_y[6], '--h', color='peru', label='0.2')
    plt.plot(logsurvival_all_x[7], logsurvival_all_y[7], '--o', color='firebrick', label='0.3')
    plt.plot(logsurvival_all_x[8], logsurvival_all_y[8], '--^', color='green', label='0.4')
    plt.plot(logsurvival_all_x[9], logsurvival_all_y[9], '--s', color='royalblue', label='0.5')
    plt.plot(logsurvival_all_x[10], logsurvival_all_y[10], '-', color='black')
    plt.plot([0.00213482, 0.00248559], [0.122669, 0.122669], 'k-')  # horizontal
    plt.plot([0.00248559, 0.00248559], [0.122669, 0.184225], 'k-')  # vertical
    plt.text(0.00213482 + (0.00248559 - 0.00213482), 0.122669 -  (0.184225 - 0.122669), "2", fontsize=18, fontname='Times New Roman', ha='left', va='bottom')
    #plt.plot(logsurvival_all_x[1], logsurvival_all_y[1], '--o', color='firebrick')
    #plt.plot(logsurvival_all_x[0], logsurvival_all_y[0], '--o', color='firebrick')
    plt.text(0.00034, 2.20, r'$\\zeta$', fontsize=18, fontname="Times New Roman")
    plt.legend(prop={'family':'Times New Roman', 'size':'12'}, loc=(0.02, 0.65), ncols=2, frameon=False)
    plt.ylabel(r'$-log(S(\\dot{\varepsilon}))$', fontname='Times New Roman', fontsize=18)
    plt.xlabel(r'$\\dot{\varepsilon}$', fontname='Times New Roman', fontsize=18)
    #plt.ylabel(r'$S(\\dot{\varepsilon})$', fontname='Times New Roman', fontsize=18)
    #plt.xlabel(r'$\\dot{\varepsilon}$', fontname='Times New Roman', fontsize=18)
    plt.xticks(fontname='Times New Roman', fontsize=18)
    plt.yticks(fontname='Times New Roman',  fontsize=18)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([1e-3, 5])
    plt.subplots_adjust(left=0.18, bottom=0.2, right=0.955, top=0.97)
    '''
    #plt.text(0.00016, 0.725, r'$\zeta$', fontsize=18, fontname="Times New Roman")
    #plt.text(0.00007, 0.25, r'$empty: \xi/\xi_{cell}=0.01$', fontsize=14, fontname="Times New Roman")
    #plt.text(0.00007, 0.1, r'$full: \xi/\xi_{cell}=0.001$', fontsize=14, fontname="Times New Roman")
    #plt.legend(prop={'family':'Times New Roman', 'size':'12'}, loc=(0.01, 0.325), ncols=2, frameon=False)
    '''
    fig = plt.gcf()
    ax = plt.gca()
    inset_ax = fig.add_axes([0.675, 0.285, 0.28, 0.24])
    inset_ax.plot(survival_all_x[2], survival_all_y[2], ':o', color='firebrick', ms=3, markerfacecolor='none')
    inset_ax.plot(survival_all_x[3], survival_all_y[3], ':^', color='green', ms=3, markerfacecolor='none')
    inset_ax.plot(survival_all_x[4], survival_all_y[4], ':s', color='royalblue', ms=3, markerfacecolor='none')
    inset_ax.plot(survival_all_x[5], survival_all_y[5], '--p', color='goldenrod', ms=3)
    inset_ax.plot(survival_all_x[6], survival_all_y[6], '--h', color='peru', ms=3)
    inset_ax.plot(survival_all_x[7], survival_all_y[7], '--o', color='firebrick', ms=3)
    inset_ax.plot(survival_all_x[8], survival_all_y[8], '--^', color='green', ms=3)
    inset_ax.plot(survival_all_x[9], survival_all_y[9], '--s', color='royalblue', ms=3)
    inset_ax.set_ylabel(r'$S(\\dot{\varepsilon})$', fontname='Times New Roman', fontsize=12)
    inset_ax.set_xscale('log')
    for label in inset_ax.get_xticklabels():
        label.set_fontproperties('Times New Roman')
        label.set_fontsize(12)
    for label in inset_ax.get_yticklabels():
        label.set_fontproperties('Times New Roman')
        label.set_fontsize(12)

    plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/fig2.png", transparent=True)
    plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/fig2.svg", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig2.png", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig2.svg", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig2_inset.png", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig2_inset.svg", transparent=True)
    plt.show()
    '''


    #plt.text(0.0024, 0.015, r'$\zeta$', fontsize=18, fontname="Times New Roman")
    #plt.text(0.000175, 2.20, r'$full: \xi/\xi_{cell}=0.01$', fontsize=14, fontname="Times New Roman")
    #plt.text(0.000175, 0.9, r'$empty: \xi/\xi_{cell}=0.001$', fontsize=14, fontname="Times New Roman")
    #plt.legend(prop={'family':'Times New Roman', 'size':'12'}, loc=(0.55, 0.05), ncols=2, frameon=False)
    #plt.text(0.0039, 0.95, r'$\zeta$', fontsize=18, fontname="Times New Roman")
    #plt.text(0.00295, 0.5, r'$empty: \xi/\xi_{cell}=0.01$', fontsize=14, fontname="Times New Roman")
    #plt.text(0.00305, 0.375, r'$full: \xi/\xi_{cell}=0.001$', fontsize=14, fontname="Times New Roman")
    #plt.legend(prop={'family':'Times New Roman', 'size':'12'}, loc=(0.52, 0.55), ncols=2, frameon=False)

    #plt.figure(figsize=(3*40.179/25.4, 3*25.142/25.4))
    '''
    plt.figure(figsize=(2*86/25.4, 2*38.3/25.4))
    plt.plot(st_all_x[0], st_all_y[0], ':v', color='firebrick', label='0.3')
    plt.plot(st_all_x[1], st_all_y[1], ':d', color='green', label='0.4')
    plt.plot(st_all_x[2], st_all_y[2], ':H', color='royalblue', label='0.5')
    plt.plot(st_all_x[3], st_all_y[3], '--v', color='firebrick', markerfacecolor='none')
    plt.plot(st_all_x[4], st_all_y[4], '--d', color='green', markerfacecolor='none')
    plt.plot(st_all_x[5], st_all_y[5], '--H', color='royalblue', markerfacecolor='none')
    #plt.plot(st_all_x[6], st_all_y[6], ':p', color='firebrick', label='0.3')
    #plt.plot(st_all_x[7], st_all_y[7], ':h', color='green', label='0.4')
    #plt.plot(st_all_x[8], st_all_y[8], ':D', color='royalblue', label='0.5')
    plt.text(0.33, 0.875, r'$\\zeta$', fontsize=18, fontname="Times New Roman")
    plt.legend(prop={'family':'Times New Roman', 'size':'12'}, loc=(0.4, 0.8), ncols=3, frameon=False)
    plt.ylabel(r'$S(\tau_{hole})$', fontname='Times New Roman', fontsize=18)
    plt.xlabel(r'$\tau_{hole}/t_{total}$', fontname='Times New Roman', fontsize=18)
    plt.xticks(fontname='Times New Roman', fontsize=18)
    plt.yticks(fontname='Times New Roman', fontsize=18)
    plt.subplots_adjust(left=0.15, bottom=0.22, right=0.965, top=0.99)
    #plt.subplots_adjust(left=0.18, bottom=0.22, right=0.955, top=0.99)
    plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/fig3c.png", transparent=True)
    plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/fig3c.svg", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig3c.png", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig3c.svg", transparent=True)
    plt.show()
    '''

    exit(1)

plt.close()

if variable == 4:

    fig = plt.figure(figsize=(5.452423529,4.089317647), constrained_layout=True)
    PD_activity = [0.1, 0.2, 0.3, 0.4, 0.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.1, 0.2, 0.3, 0.4, 0.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.1, 0.2, 0.3, 0.4, 0.5]
    PD_friction = [0.1, 0.1, 0.1, 0.1, 0.1, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.75, 0.5, 0.25, 0.1, 0.075, 0.05, 0.025, 0.01, 0.0075, 0.005, 0.0025, 0.001, 0.05, 0.05, 0.05, 0.05, 0.05, 0.025, 0.025, 0.025, 0.025, 0.025, 0.0075, 0.0075, 0.0075, 0.0075, 0.0075, 0.005, 0.005, 0.005, 0.005, 0.005, 0.0025, 0.0025, 0.0025, 0.0025, 0.0025]
    PD_phase = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 2, 0, 0, 1, 1, 2, 0, 1, 1, 2, 2]

    #colors = cm.Set1(np.linspace(0, 1, 9))
    colors = ['firebrick', 'royalblue', 'green']
    markers = ['o', 's', '^']
    labnames = ['No holes', 'Short-lived', 'Long-lived']

    for i in range(len(PD_friction)):
        if i == 0:
            plt.scatter(PD_activity[i], PD_friction[i], color=colors[PD_phase[i]], marker=markers[PD_phase[i]], label=labnames[0])
        elif i == 7:
            plt.scatter(PD_activity[i], PD_friction[i], color=colors[PD_phase[i]], marker=markers[PD_phase[i]], label=labnames[1])
        elif i == 13:
            plt.scatter(PD_activity[i], PD_friction[i], color=colors[PD_phase[i]], marker=markers[PD_phase[i]], label=labnames[2])
        else:
            plt.scatter(PD_activity[i], PD_friction[i], color=colors[PD_phase[i]], marker=markers[PD_phase[i]])


    #plt.ylabel(r'$\chi$', fontsize=18, fontname='Times New Roman')
    plt.ylabel(r'$\xi/\xi_{cell}$', fontsize=18, fontname='Times New Roman')
    plt.xlabel(r'$\zeta$', fontsize=18, fontname='Times New Roman')
    plt.xticks(fontsize=18, fontname='Times New Roman')
    plt.yticks(fontsize=18, fontname='Times New Roman')
    plt.xticks(np.arange(0.1,0.6,0.1))
    #plt.legend(prop={'family':'Times New Roman', 'size':'12'}, loc=(0.025, 0.72), ncols=1, frameon=True, facecolor='white', edgecolor='black', framealpha=1.0)
    plt.yscale('log')
    plt.ylim([0.0007, 0.15])
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Results7/Phase_diagram_1.png", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Results7/Phase_diagram_1.svg", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/fig3d.png", transparent=True)
    #plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/fig3d.svg", transparent=True)
    plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig3d.png", transparent=True)
    plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Slides/Presentation/fig3d.svg", transparent=True)
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
#MPF holes paper----------------------------------------------------------------------------------------------



'''
#stress_time plot
plt.figure(figsize=(5.452423529,4.089317647))
avg_vel=[]
for job in range(1, 6, 1):
    #fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/stress_time.txt","r")
    if job==4:
        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts50/Job_4/stress_time.txt","r")
    else:
        fileoutput=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(job)+"/stress_time.txt","r")
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
filedata=open("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/params","r")
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
        #fname = "/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts72/Job_"+str(j)+"/velocity_correlations.txt"
        fname = "/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts72/Job_"+str(j)+"/nematic_correlations.txt"
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
        #fname = "/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts72/Job_"+str(j)+"/velocity_correlations.txt"
        fname = "/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts72/Job_"+str(j)+"/nematic_correlations.txt"
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
        #fname = "/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(j)+"/velocity_correlations.txt"
        fname = "/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/Job_"+str(j)+"/nematic_correlations.txt"
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
            #fname = "/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts71/Job_"+str(j)+"/velocity_correlations.txt"
            fname = "/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts71/Job_"+str(j)+"/nematic_correlations.txt"
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
            #fname = "/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts71/Job_"+str(j)+"/velocity_correlations.txt"
            fname = "/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts71/Job_"+str(j)+"/nematic_correlations.txt"
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
plt.ylabel(r'$C_\mathbf{v}$', fontsize=18)
#plt.ylabel(r'$C_\mathbf{Q}$', fontsize=18, fontname='Times New Roman')
plt.xlabel('R', fontsize=18, fontname='Times New Roman')
plt.xticks(fontname='Times New Roman', fontsize=18)
plt.yticks(fontname='Times New Roman', fontsize=18)
#plt.legend([r'$\omega$=0.01',r'$\omega$=0.1',r'$\omega$=0.4'], fontsize=14, loc=(0.,0.65), frameon=False)
plt.legend([r'$\xi_{cell}=0$', r'$\chi=1$',r'$\chi=0.1$',r'$\chi=0.01$', r'$\chi=0.001$',r'$\chi=0.0001$'], loc=(0.6,0.37), prop={'family':'Times New Roman', 'size':'12'}, frameon=False)
plt.subplots_adjust(left=0.235, bottom=0.235, right=0.95, top=0.95)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/C_Q_soft.png", transparent=True)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/C_Q_soft.svg", transparent=True)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/scripts"+str(scripts)+"/mean_velocity_1.png", transparent=True)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/figSM1b.png", transparent=True)
#plt.savefig("/home/diogo/Phase_Field/RheoCell/Work/Analysis/Paper/figSM1b.svg", transparent=True)
plt.show()

exit (1)


'''
#other details for plotting (from older scripts)
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
	plt.xlabel(r'$\\gamma$', fontsize=18, fontname='Times New Roman')
elif variable==2:
	#plt.xlabel(r'$\rho$', fontsize=18, fontname='Times New Roman')
	plt.xlabel(r'$\\gamma$', fontsize=18, fontname='Times New Roman')
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
