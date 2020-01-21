import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.ticker
import os
import time
from numpy.polynomial import polynomial as P

#------------------------------------------------
#--------------functions--------------------------
#------------------------------------------------

#this function test if the argument is a number
def isnumber(value):
    try:
         float(value)
    except ValueError:
         return False
    return True

#this function test if the argument is an integer number
def isint(value):
    try:
         int(value)
    except ValueError:
         return False
    return True

#this function return y(x) for a polynomial with the coefficients and the x in the arguments
def p_function(coefficients,x):
	sum_ = 0
	for i in range(0,len(coefficients)):
		sum_ += coefficients[i]*(x**i)

	return sum_

def separate_curves(H, M, reverse = False):
	downs = []
	ups = []

	aux_d = []
	aux_u = []

	down = True

	max_m = 0
	for i in range(0,len(H)-1):	
		
		if(down):
			aux_d.append([H[i], M[i]])
			if(H[i] < H[i+1]):
				if(M[i] > max_m):
					max_m = M[i]
				aux_u.append([H[i], M[i]])
				down = False		
			
		else:
			aux_u.append([H[i], M[i]])
			if(H[i] > H[i+1] or i == len(H)-2):	#o or eh para pegar o ultimo par de curvas		
				down = True
				if(M[i] > max_m):
					max_m = M[i]
				if(i == len(H)-2): 			#para pegar ultimo termo
					aux_u.append([H[i+1], M[i+1]])
				downs.append(aux_d)
				ups.append(aux_u)
				aux_d = []
				aux_u = []	

	return downs, ups, max_m

def get_hys(fst_d, last_d, fst_u, last_u):
	x = []
	y = []
	for i in range(0,len(last_d)):					
		x.append(last_d[i][0])
		y.append(last_d[i][1])			

	x_d = []
	y_d = []

	for i in range(0,len(fst_d)):
		if fst_d[i][0] > max(x):
			x_d.append(fst_d[i][0])
			y_d.append(fst_d[i][1])

	x_d = x_d + x
	y_d = y_d + y

	x = [] #H
	y = [] #M		
	for i in range(0, len(last_u)):			
		x.append(last_u[i][0])
		y.append(last_u[i][1])		

	x_u = x
	y_u = y

	for i in range(0,len(fst_u)):
		if fst_u[i][0] > max(x_u):
			x_u.append(fst_u[i][0])
			y_u.append(fst_u[i][1])

	return x_d, y_d, x_u, y_u

def max_der(x, y, c):
	der = [p_function(c,item) for item in x]
	max_d = 0

	for i in range(0,len(x)):
		if(max_d < der[i]):
			max_d = der[i]
			Hsw = x[i]		
			M_Hsw = y[i]

	der = np.asarray(der)

	return der, max_d, Hsw, M_Hsw

def extend_curves(curves, x_hys, y_hys):

	hys = [[x,y] for x,y in zip(x_hys,y_hys)]

	new_curves = []

	for curve in curves:
		ex = filter(lambda item: item[0] < min(curve)[0],hys)	
		ex.reverse()	
		new_curves.append(ex+curve)

	return new_curves


class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_orderOfMagnitude(self, nothing):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin, vmax):
        self.format = self.fformat
        if self._useMathText:
            self.format = '$%s$' % matplotlib.ticker._mathdefault(self.format)

#-----------------------------------------------
#------------printing the header---------------
#-----------------------------------------------
print(' _________________________________________________________________________\n')
print(' This program generates delta-M_R plot(s) from recoil loop(s), i,e, closed FORCs,')
print(' each of them started at saturation, and/or a saturated major loop.')
print(' An input data file (two columns, magnetization, M, versus magnetic field,')
print(' H) is required.\n')


print(' dMr v5.0 a stable beta version writen in Python by Luana Lazzarotto Bianchi')
print(' (luhbianchi06@gmail.com) under the supervision of Dr. Julian Geshev,')
print(' Laboratory of Magnetism - LAM (UFRGS, Brazil). For updates, please check')
print(' http://www.if.ufrgs.br/pes/lam/dMr.html')
print(' ________________________________________________________________________\n')


#printing the header of the option
# 1 for recoils/major loops and the shifts
# 2 for recoils loops and/or major loops. In this case the shifts will be caculated

print(' \tIf the data file contains recoil loop(s) only (in this case, you must ')
print(' provide the shift values along the H and M axes, if any), you must  ENTER 1.')
print(' \tAlternatively, if the data file contains one or a sequence of recoil ')
print(' loop(s) together with the respective major hysteresis loop, please ENTER 2.') 
print(' In this case, an image will pop-up showing the major loop, the coercive') 
print(' and swithching fileds and the corresponding shift values.') 
print(' \tThe program suggests the M- and H-shifts obtained from the') 
print(' coercive fields; however, you can enter other shift values.\n') 


print(' ENTER:')
print("\t'1' (recoil loop(s) & the shifts' of the major loop) or")

#------------------------------------------------------------------------------------
option = raw_input("\t'2' (recoil loop(s) & major loop): 1 ") or 1


#testing if the enter is a number
if(isint(option)):
	option = int(option)

#testing if the enter is 1 or 2
if(option != 1 and option != 2):
	option = raw_input("\n Invalid value!\n Please, enter 1 or 2: ") or 1	

	if(isint(option)):
		option = int(option)

	while (option != 2 and option != 1):					
		option = raw_input("\n Invalid value!\n Please, enter 1 or 2: ") or 1	
		if(isint(option)):
			option = int(option)

print('\n')

#------------------------------------------------------------------------------------

FileName = raw_input(" Enter the name of the input data file (.dat or .VHD): rec.dat " ) or ('rec.dat')

name_vector = FileName.split('.')

if(len(name_vector) != 2):
	print(' File without extension')
	time.sleep(5)
	os.sys.exit()

if (name_vector[1] == 'VHD'):
	try:
		InputFile = open(FileName).read().splitlines()
	except IOError:
		print(' File not found. Please enter the correct file')
		time.sleep(5)
		os.sys.exit()

	#select the columns with H and M
	for i in range(0,len(InputFile)-1):
		line = InputFile[i]
		if ': Applied Field For Plot' in line: 
			H_column = int(filter(str.isdigit, line))
		if ': Signal X direction' in line: 
			Mx_column = int(filter(str.isdigit, line))
		if '@@End of Header.' in line: 
			begin_table = i + 5
		if '@@END Data.' in line: 
			end_table = i
			break
 
	name_file = FileName.replace(".VHD", ".dat")	

	OutputFile = open(name_file, 'w')	

	M = np.array([])
	H = np.array([])

	#para pegar as colunas com os dados do arquivo VHD
	for i in range(begin_table, end_table):
		columns = InputFile[i].split()
		M = np.append(M, float(columns[Mx_column]))
		H = np.append(H, float(columns[H_column]))
		OutputFile.write("%s %s\n" % (columns[H_column], columns[Mx_column]))			
		
	OutputFile.close()

elif name_vector[1] == 'dat':
	try:
		H, M = np.loadtxt(FileName, unpack=True)
	except IOError:
		print(' File not found. Please enter the correct file name.')
		time.sleep(10)
		os.sys.exit()
	except ValueError:
		print(' There are non-numeric values in the input data file; these will be ignored.')

		arq = open(FileName,"r")

		M = np.array([])
		H = np.array([])

		for linha in arq:
			valores = linha.split()			
			
			if(isnumber(valores[1]) and isnumber(valores[0])):
				H = np.append(H, float(valores[0]))
				M = np.append(M, float(valores[1]))
							
		arq.close()	
else:
	print(' Wrong file extension.')
	time.sleep(5)
	os.sys.exit()


print('\n')

#------------------------------------------------------------------------------------
try:		
	massa = float(raw_input(" Enter the mass (in g) of the sample in order to have M in emu/g \n (or simply press ENTER for M/Mmax): ") or 0)
except ValueError:
	massa = raw_input("\n Invalid value.\n Enter a NUMBER for the mass (in g) of the sample in order to have M in \nemu/g, or simply press ENTER for M/Mmax:") or 0
	while (not isnumber(massa)):
		massa = raw_input("\n Invalid value.\n Enter a NUMBER for the mass (in g) of the sample in order to have M in \nemu/g, or simply press ENTER for M/Mmax:") or 0

	massa = float(massa)
#------------------------------------------------------------------------------------
print('\n')

if (option == 1):
	#------------------------------------------------------------------------------------
	try:
		Hshift = float(raw_input(" ENTER the value of the shift along the horizontal axis: 0 ") or 0) 
	except ValueError:
		Hshift = raw_input("\n Invalid value.\n Please, enter a NUMBER for the shift along the horizontal axis: 0 ") or 0
		while (not isnumber(Hshift)):		
			Hshift = raw_input("\n Invalid value.\n Please, enter a NUMBER for the shift along the horizontal axis: 0 ") or 0

		Hshift = float(Hshift)
		
	#-----------------------------------------------------------------------------------
	#------------------------------------------------------------------------------------
	try:
		Mshift = float(raw_input(" ENTER the value of the shift along the vertical axis:   0 ") or 0) 
	except ValueError:
		Mshift = raw_input("\n Invalid value.\n Please, enter a NUMBER for the shift along the vertical axis: 0 ") or 0
		while (not isnumber(Mshift)):		
			Mshift = raw_input("\n Invalid value.\n Please, enter a NUMBER for the shift along the vertical axis: 0 ") or 0

		Mshift = float(Mshift)
	print('\n')	
#------------------------------------------------------------------------------------


inv = False

if(H[0] < H[1]):	
	inv = True
	H = -H
	M = -M

#chama a funcao que separa curvas
downs, ups, Mmax = separate_curves(H,M)

#get the hysteresis 
x_d, y_d, x_u, y_u = get_hys(downs[0],downs[-1],ups[0],ups[-1])

#interpola as curvas obtidas
hys_curve_up = interp1d(x_u,y_u, kind='linear')
hys_curve_d = interp1d(x_d,y_d, kind='linear')

#calculate the major area
x = np.arange(max(min(x_u),min(x_d)), min(max(x_u),max(x_d)), 0.1)

y_2 = hys_curve_up(x)
y_1 = hys_curve_d(x)

y = y_1 - y_2

area_major = np.trapz(y,x)

if (option == 2):
	#find the Heb

	dx = 0.1

	x_inter_d = np.arange(min(x_d), max(x_d), dx)
	x_inter_up = np.arange(min(x_u), max(x_u), dx)

	y_inter_d = hys_curve_d(x_inter_d)
	y_inter_up = hys_curve_up(x_inter_up)

	#-------Campo coercivo-----------

	for i in range(1,len(x_inter_d)):
		
		if(y_inter_d[i-1] * y_inter_d[i] < 0 ):		
			Hc_n =  x_inter_d[i]
			break

	for i in range(1,len(x_inter_up)):
		if(y_inter_up[i-1] * y_inter_up[i] < 0 ):
			Hc =  x_inter_up[i]
			break

	
	#-------sw field-----------------------------------------	

	#-----find limits to derivate-------

	#-----------ascending

	for i in range(0, len(x_u)):
		if(y_u[i] > max(M)*0.6):
			limit_uu = x_u[i];
			break

	for i in range(0, len(x_u)):
		if(y_u[i] > min(M)*0.6):
			limit_ud = x_u[i];
			break

	#----------descending
	
	for i in range(0, len(x_d)):		
		if(y_d[i] < max(M)*0.6):			
			limit_du = x_d[i];
			break

	for i in range(0, len(x_d)):
		if(y_d[i] < min(M)*0.6):			
			limit_dd = x_d[i];
			break	

	x_inter_d = np.arange(limit_dd, limit_du, dx)
	y_inter_d = hys_curve_d(x_inter_d)

	x_inter_up = np.arange(limit_ud, limit_uu, dx)
	y_inter_up = hys_curve_up(x_inter_up)

	de = np.gradient(y_inter_up, x_inter_up)
	de2 = np.gradient(y_inter_d, x_inter_d)

	grau = 9

	c1 = P.polyfit(x_inter_up,de,grau,full=True)[0]
	c2 = P.polyfit(x_inter_d,de2,grau,full=True)[0]

	#-------------------------------------------------------
	
	#--------find max der

	pos1 = 0

	de_p ,max1, Hsw, M_Hsw = max_der(x_inter_up, y_inter_up, c1)
	de2_p ,max2, Hsw_n, M_Hsw_n = max_der(x_inter_d, y_inter_d, c2)


	Heb = (Hsw + Hsw_n)/2

	#------- plot Hsw and Hc---------------------------------
	
	if(inv):
		plt.plot(-np.asarray(x_u),-np.asarray(y_u), color='blue')
		plt.plot(-np.asarray(x_d),-np.asarray(y_d), color='blue')
		plt.plot(-x_inter_up, ((de_p-min(de_p))*(max(M)*0.5)/max(de_p)), color='red', )
		plt.plot(-x_inter_d, ((de2_p-min(de2_p))*(max(M)*0.5)/max(de2_p)), color='red')
	else:
		plt.plot(x_u,y_u, color='blue')
		plt.plot(x_d,y_d, color='blue')
		plt.plot(x_inter_up, (de_p-min(de_p))*(max(M)*0.5)/max(de_p), color='red', )
		plt.plot(x_inter_d, (de2_p-min(de2_p))*(max(M)*0.5)/max(de2_p), color='red')

	x = np.arange(limit_dd, limit_uu, dx)
	y_inter_up = hys_curve_up(x)
	y_inter_d = hys_curve_d(x)


	for i in range(0,len(x)):
		if(int(x[i]) == int(Heb)):
			pos1 = i

		if(int(x[i]) == int(Heb)):
			pos2 = i

	Mr1 = y_inter_up[pos1]
	Mr2 = y_inter_d[pos2]

	vertical_center = (Mr1 + Mr2)/2

	if inv:
		label1 = '''$H_{0}$ = {1:.1f}\t\t$H_{2}$ = {3:.1f}'''.format('{sw1}',-Hsw_n, '{c1}', -Hc_n)
		label2 = '''$H_{0}$ = {1:.1f}\t\t$H_{2}$ = {3:.1f}'''.format('{sw2}',-Hsw, '{c2}', -Hc)
		label3 = '''$H_{0}(H_{2})$ = {1:.1f}\t$H_{0}(H_c)$ = {3:.1f}'''.format('{eb}',-Heb,'{sw}',-(Hc_n+Hc)/2)
		
		plt.scatter(-Hsw_n,-M_Hsw_n, marker = 'x', color="black", s=500, label = label1)
		plt.scatter(-Hsw,-M_Hsw, marker = 'x', color="m", s=500, label = label2)
		plt.scatter(-Heb,-vertical_center, marker = 'x', color="c", s=500, label = label3)
		plt.scatter(-Hc_n,0, marker = 'o', color="c", s=50)
		plt.scatter(-Hc,0, marker = 'o', color="c", s=50)
	else:
		label1 = '''$H_{0}$ = {1:.1f}\t\t$H_{2}$ = {3:.1f}'''.format('{sw1}',Hsw_n, '{c1}', Hc_n)
		label2 = '''$H_{0}$ = {1:.1f}\t\t$H_{2}$ = {3:.1f}'''.format('{sw2}',Hsw, '{c2}', Hc)
		label3 = '''$H_{0}(H_{2})$ = {1:.1f}\t$H_{0}(H_c)$ = {3:.1f}'''.format('{eb}',Heb,'{sw}',(Hc_n+Hc)/2)
		
		plt.scatter(Hsw_n,M_Hsw_n, marker = 'x', color="black", s=500, label = label1)
		plt.scatter(Hsw,M_Hsw, marker = 'x', color="m", s=500, label = label2)
		plt.scatter(Heb,vertical_center, marker = 'x', color="c", s=500, label = label3)
		plt.scatter(Hc_n,0, marker = 'o', color="c", s=50)
		plt.scatter(Hc,0, marker = 'o', color="c", s=50)


	plt.axhline(color = 'black')
	plt.axvline(color = 'black')

	plt.xlabel('$H$ (Oe)', fontsize=25)
	plt.ylabel('$M$', fontsize=25)
	plt.xticks(fontsize=18, rotation=0)
	plt.yticks(fontsize=18, rotation=0)
	plt.legend(loc='upper left')

	plt.show()

	if(inv):
		print(''' H_sw1 = {0:.1f}\t\tH_c1 = {1:.1f}'''.format(-Hsw_n, -Hc_n))
		print(''' H_sw2 = {0:.1f}\t\tH_c2 = {1:.1f}'''.format(-Hsw, -Hc))
		print(''' H_eb(H_sw) = {0:.1f}\tH_eb(H_c) = {1:.1f}'''.format(-Heb,-(Hc_n+Hc)/2))
	else:
		print(''' H_sw1 = {0:.1f}\t\tH_c1 = {1:.1f}'''.format(Hsw_n, Hc_n))
		print(''' H_sw2 = {0:.1f}\t\tH_c2 = {1:.1f}'''.format(Hsw, Hc))
		print(''' H_eb(H_sw) = {0:.1f}\tH_eb(H_c) = {1:.1f}'''.format(Heb,(Hc_n+Hc)/2))


	if(inv):
		aux_eb = -(Hc_n+Hc)/2
		str_Heb = str('''{:.1f}'''.format(-(Hc_n+Hc)/2))
	else:
		aux_eb = (Hc_n+Hc)/2
		str_Heb = str('''{:.1f}'''.format((Hc_n+Hc)/2))

	#--------------------------------------------------------------------------
	print('\n')
	try:
		Heb = float(raw_input(" ENTER the value of the shift along the horizontal axis: "+str_Heb+" ") or aux_eb) 
	except ValueError:
		Heb = raw_input("\n Invalid value.\n Please, enter a NUMBER for the shift along the horizontal axis:"+str_Heb+" ") or aux_eb
		while (not isnumber(Heb)):		
			Heb = raw_input("\n Invalid value.\n Please, enter a NUMBER for the shift along the horizontal axis:"+str_Heb+" ") or aux_eb

		Heb = float(Heb)


	
	try:
		vertical_center = float(raw_input(" ENTER the value of the shift along the vertical axis: 0 ") or 0) 
	except ValueError:
		vertical_center = raw_input("\n Invalid value.\n Please, enter a NUMBER for the shift along the vertical axis: 0 ") or 0
		while (not isnumber(vertical_center)):		
			vertical_center = raw_input("\n Invalid value.\n Please, enter a NUMBER for the shift along the vertical axis: 0 ") or 0

		vertical_center = float(vertical_center)

	if inv:
		vertical_center = -vertical_center
		Heb = -Heb


	#--------------------------------------------------------------------------

	#------------------------------------------------------------------------------------
	#----------------------------------------------------------------------------------'''

	file = open("centered_loop(s).dat", "w")

	H = [x-Heb for x in H]
	M = [y-vertical_center for y in M]

	if inv:
		H = [-x for x in H]
		M = [-y for y in M]

	for i in range(0,len(H)):		
		file.write("%s\t%s \n" % (H[i], M[i]))

	file.close()

	Mshift = vertical_center
	Hshift = Heb

	file = open("major-loop_parameters.dat", "w")

	if(inv):
		file.write('''H_sw1 = {0:.2f},\tH_sw2 = {1:.2f},\tH_eb(H_sw) = {2:.2f}
Mr_1 = {3:.2e},\tMr_2 = {4:.2e},\tMshift = {5:.2e}

H_c1 = {6:.2f},\tH_c2 = {7:.2f},\tH_eb(H_c) = {8:.2f}'''.format(-Hsw_n,-Hsw,-Hshift, -Mr2, -Mr1, -Mshift, -Hc_n, -Hc,-(Hc_n+Hc)/2))

	else:
		file.write('''H_sw1 = {0:.2f},\tH_sw2 = {1:.2f},\tH_eb(H_sw) = {2:.2f}
Mr_1 = {3:.2e},\tMr_2 = {4:.2e},\tMshift = {5:.2e}

H_c1 = {6:.2f},\tH_c2 = {7:.2f},\tH_eb(H_c) = {8:.2f}'''.format(Hsw_n,Hsw,Hshift, Mr2, Mr1, Mshift, Hc_n, Hc,(Hc_n+Hc)/2 ))
	
	file.write('''\n\nmajor's area = {0:.2e}'''.format(area_major))

	file.close()

#if we dont have a major loop
if (option == 1):
	file = open("centered_loop(s).dat", "w")
	file2 = open("major-loop_parameters.dat", "w")

	file2.write('''H_shift = {:.2f}
M_shift = {:.2e}'''.format(Hshift, Mshift))

	if(inv):
		Hshift = -Hshift
		Mshift = -Mshift	

	H = [x-Hshift for x in H]
	M = [y-Mshift for y in M]
			
	for i in range(0,len(H)):
		if(inv):
			file.write("%s\t%s \n" % (-H[i], -M[i]))
					
		else:
			file.write("%s\t%s \n" % (H[i], M[i]))
		
	file.close()
	file2.close()


#centralize the major loop
x_d = [x-Hshift for x in x_d]
y_d = [y-Mshift for y in y_d]

x_u = [x-Hshift for x in x_u]
y_u = [y-Mshift for y in y_u]

H, M = np.loadtxt("centered_loop(s).dat", unpack=True)

down = True
inv = False

if(H[0] < H[1]):
	inv = True
	H = -H
	M = -M
	
downs, ups, Mmax = separate_curves(H,M)

reverse_fields = [up[0] for up in ups]


if(option == 1 and len(downs) == 1):
	aux = np.arange(-max(H), min(x_d)-1, 1)
	aux = aux[::-1]

	for x in aux:
		x_d.append(x)
		y_d.append(0)


#extend the up curves
ups = extend_curves(ups, x_d, y_d)



#rotate the major dsc, using that the curve is centralized
x_d_sym = [-x for x in x_d]
y_d_sym = [-y for y in y_d]


Hmax_d = max(x_d)

#interpolate the major dsc, asc and the rotate major
hys_curve_d = interp1d(x_d,y_d, kind='linear')
hys_curve_d_sym = interp1d(x_d_sym,y_d_sym, kind='linear')
hys_curve_up = interp1d(x_u,y_u, kind='linear')


file = open("centered_loop(s).dat", 'w')
file.write("H\tM\n")

for x,y in zip(x_d, y_d):
	file.write("%.2f\t%.2e\n" % (x,y))

for x,y in zip(x_u, y_u):
	file.write("%.2f\t%.2e\n" % (x,y))


file.close()

#-----------------------------------------------------------------

up_curves = []
Hs = []

limit_d = max([up[0] for up in ups])[0] #Hrec field max



for up in ups:

	x = [point[0] for point in up]
	y = [point[1] for point in up]

	f = interp1d(x,y, kind='linear')

	Hs.append(int(max(x)))

	if(Hmax_d < Hs[-1]):
		Hs[-1] = Hmax_d	

	up_curves.append(f) 

if(option == 2):
	Hc_p = Hc - Hc_n #Hc to the step
elif(option == 1):
	Hc_p = Hshift - limit_d

p = Hc_p/50

#if p  > max(step) in the data, p = max(step)
if p > max([abs(x_d[i] - x_d[i-1]) for i in range(1,len(x_d))]):
	p = max([abs(x_d[i] - x_d[i-1]) for i in range(1,len(x_d))])

p = round(p,1)

if p < 0.5:
	p = 0.5

aux_p = p

print('\n')

try:	
	p = float(raw_input(' Enter the field step: '+ str(p) + ' ') or p)
except ValueError:
	p = raw_input('\n Invalid value!\n Press ENTER or type an INTEGER NUMBER for the field step: '+ str(aux_p)+ ' ') or aux_p
	while not isnumber(p):
		p = raw_input('\n Invalid value!\n Press ENTER or type an INTEGER NUMBER for the field step: '+ str(aux_p)+ ' ') or aux_p

	p = float(p)
	#p = round(p,1)

recoil_curves = []
fields = []
for i in range(0,len(reverse_fields)):	
	g = up_curves[i]

	xnew = np.arange(int(limit_d)+1, 1, p)
	ynew = g(xnew)

	#create the M_sym
	recoil_curves.append([[-x,-y] for x, y in zip(xnew, ynew)]) 
	recoil_curves[i].reverse()

	fields.append([-x for x in xnew])
	fields[-1].reverse()


OutputFile1 = open('dMrs.dat', 'w')
info_file = open('major-loop_parameters.dat', 'a')
area_file = open('area.dat', 'w')
info_file.write("\n\n")

area_file.write("Hrec\tHofMax\tdMr_max\tHofMin\tdMr_min\tApos\tAneg\n")
OutputFile1.write("H\tdMr\tHrec\tMasc\tMdsc\n")


#------------------ORGANIZACAO DE LABEL DOS GRAFICOS-----------------
#--------------------------------------------------------------------
label_M = " (Emu/g)"
if(massa == 0):
	label_M = "/M_{max}"	
	massa = Mmax
#--------------------------------------------------------------------
#--------------------------------------------------------------------

#----------------------------------------------#
#----------------------------------------------#
#----------------------------------------------#
#--------------CALCULATE DELTA M_R-------------#
#----------------------------------------------#
#----------------------------------------------#
#----------------------------------------------#
b = np.array([]).reshape(0,3) #Ha, Hb, M''
max_delta = []
curves = []
for i in range(0,len(fields)):
	xnew = fields[i];

	Hr = reverse_fields[i][0] #campo de reversao
	Hr_p = Hr + Hshift
	if inv:
		Hr_p = -Hr_p

	up_curve = up_curves[i]

	x_del = []
	delta = []
	mr = []
	mhys = []

	if(Hr > 0):
		xnew = np.arange(Hr, Hs[i], p)

		y_hys_d = hys_curve_d(xnew)			#M da descida da histerese
		y_up = up_curve(xnew)				#M da subida do forc para Hr

		j = 0
		for j in range(0, len(xnew)):

			Mr = y_up[j]
			Mhys = y_hys_d[j]

			mr.append(Mr)
			mhys.append(Mhys)

			delta_m = -Mhys + Mr
			
			x = xnew[j] + Hshift
			
			if(inv):
				delta_m = -delta_m
				x = -x

			x_del.append(x)			
			delta.append(delta_m)	

			OutputFile1.write("%s\t%s\t%s\t%s\t%s\n" % (x,(delta_m/massa), Hr_p, Mr, Mhys))
			b = np.append(b, [[x, y_up[j], delta_m]], axis=0)	
	else:
		if(len(xnew) == 0 or len(recoil_curves[i]) < 2):
			continue
		if(max(xnew) > Hs[i]):
			xnew = np.arange(min(xnew),Hs[i],p)
			if(len(xnew) == 0):
				continue

		y_hys_d = hys_curve_d(xnew)			#M da descida da histerese
		y_hys_d_sym = hys_curve_d_sym(xnew)	#M da descida espelhada da histerese
		y_up = up_curve(xnew)				#M da subida do forc para Hr

		for j in range(0, len(xnew)):

			Mr = (y_up[j] + recoil_curves[i][j][1])			
			Mhys = (y_hys_d_sym[j] + y_hys_d[j])

			mr.append(Mr)
			mhys.append(Mhys)

			delta_m = -Mhys + Mr
			
			x = xnew[j] + Hshift

			if(inv):
				delta_m = -delta_m
				x = -x
			x_del.append(x)
			
			delta.append(delta_m)	
			OutputFile1.write("%s\t%s\t%s\t%s\t%s\n" % (x,(delta_m/massa), Hr_p, Mr, Mhys))
			b = np.append(b, [[x, y_up[j], delta_m]], axis=0)

	OutputFile1.write("\n\n")

	info_file.write("H_rec = %.2f\n" % (Hr_p))		

	H_max_dmr = [x_del[k] for k in range(0, len(delta)) if delta[k] == max(delta)][0]
	H_min_dmr = [x_del[k] for k in range(0, len(delta)) if delta[k] == min(delta)][0]

	max_delta.append(max(delta))
	max_delta.append(abs(min(delta)))

	if(inv):
		delta.reverse()
		x_del.reverse()

	posPart = np.maximum(delta, 0) 
	negPart = -np.minimum(delta, 0) 

	posArea = np.trapz(posPart, x_del)
	negArea = np.trapz(negPart, x_del)

	info_file.write("H(max_dMr), max_dMr = (%.2f,%.2e) \n" % (H_max_dmr, max(delta)))
	info_file.write("H(min_dMr), min_dMr = (%.2f,%.2e) \n" % (H_min_dmr, min(delta)))
	info_file.write("Neg area = %.2e; Pos area = %.2e \n\n" % (negArea, posArea))

	area_file.write("%.2f\t%.2f\t%.2e\t%.2f\t%.2e\t%.2e\t%.2e\n" % (Hr_p, H_max_dmr, max(delta), H_min_dmr, min(delta), posArea, negArea))

	plt.plot(x_del, delta)
	curves.append([x_del, mr, mhys,Hr_p])

if(inv):
	plt.plot(-(np.asarray(x_u)+ Hshift),-(np.asarray(y_u) + Mshift), color='blue')
	plt.plot(-(np.asarray(x_d) + Hshift),-(np.asarray(y_d)  + Mshift), color='blue')	
else:
	plt.plot((np.asarray(x_u)+ Hshift),np.asarray(y_u)  + Mshift, color='blue')
	plt.plot((np.asarray(x_d) + Hshift),np.asarray(y_d)  + Mshift, color='blue')	
		
OutputFile1.close()
info_file.close()
area_file.close()

plt.axhline(color = 'black')
plt.axvline(color = 'black')

plt.xlabel('$H$ (Oe)', fontsize=25)
plt.ylabel('$\\delta M_R'+label_M+'$', fontsize=25)
plt.xticks(fontsize=18, rotation=0)
plt.yticks(fontsize=18, rotation=0)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.show()


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#--------------dHr------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

file = open('dHr.dat','w')
file.write("Masc\tH(Masc)\tH(Mdsc)\tdHr\tHrec\n")
file.write("\n\n" )	
for curve in curves:
	mr_inter = interp1d(curve[1],curve[0], kind='linear')
	mhys_inter = interp1d(curve[2],curve[0], kind='linear')

	l_d = max(min(curve[1]),min(curve[2]))+1
	l_u = min(max(curve[1]),max(curve[2]))-1

	p = (l_u - l_d)/1000

	x = np.arange(l_d,l_u ,p)

	mr = mr_inter(x)
	mhys = mhys_inter(x)

	for j in range(0, len(x)):
		file.write("%s\t%s\t%s\t%s\t%s\n" % (x[j], mr[j], mhys[j], mr[j] - mhys[j], curve[3]))

file.close()


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#--------------dMr color plot-------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

if(len(reverse_fields)>15):
	xi = np.linspace(b[:,0].min(), b[:,0].max(), 1000)
	yi = np.linspace(b[:,1].min(), b[:,1].max(), 1000) 

	xi,yi = np.meshgrid(xi,yi)
	zi = griddata((b[:,0], b[:,1]), b[:,2], (xi, yi), method='linear')

	x_u = np.asarray(x_u) + Hshift
	x_d = np.asarray(x_d) + Hshift
	y_d = np.asarray(y_d) + Mshift		
	y_u = np.asarray(y_u) + Mshift


	if(inv):
		plt.plot(-(x_u),-(y_u), color='black')
		plt.plot(-(x_d),-(y_d), color='black')	
	else:
		plt.plot(x_u, y_u , color='black')
		plt.plot(x_d, y_d , color='black')	

	
	label_y = '$M$ (emu)'
	label_x = '$H$ (Oe)'

	plt.imshow(zi, extent=(b[:,0].min(), b[:,0].max(), b[:,1].min(), b[:,1].max()), origin='lower', aspect='auto')
	plt.xlabel(label_x,  fontsize=25)
	plt.ylabel(label_y,  fontsize=25)
	plt.xticks(fontsize=18, rotation=0)
	plt.yticks(fontsize=18, rotation=0)
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.xticks(fontsize=18, rotation=0)
	plt.yticks(fontsize=18, rotation=0)

	exp = 0
	ajuste = 10**exp
	while max(max_delta)/ajuste < 1:
		exp -= 1
		ajuste = 10**exp



	cbar = plt.colorbar(format=OOMFormatter(exp, mathText=False))
	cbar.set_label('$\\delta M_R$',  fontsize=25)
	cbar.ax.tick_params(labelsize=18)


	y_2u = []
	y_2d = []
	for i in range(0, len(x_u)):
		y_2u.append(min(y_u))
	for i in range(0, len(x_d)):
		y_2d.append(max(y_d))

	plt.fill_between(x_u, y_u, y_2u, facecolor='white', interpolate=True, linewidth = "0")
	plt.fill_between(x_d, y_d, y_2d, facecolor='white', interpolate=True, linewidth = "0")

	#plt.clim(-1.5*ajuste,1*ajuste)

	cbar.formatter.set_scientific(True)
	cbar.update_ticks()
	plt.show()
