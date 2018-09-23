import pandas as pd
import numpy as np
from bokeh.plotting import show, output_file, figure
import astropy.constants as c
import astropy.units as u
import os
from elements import ELEMENTS as ele 
import collections
from scipy.spatial import cKDTree
import subprocess 
import shutil
pd.options.mode.chained_assignment = None
#top level albedo code file
init = os.getcwd()
path_to_alb_code = os.path.join(os.getcwd(),'AlbedoCode')
import scipy.signal as scisig
def run_albedo(do_clouds, three_chem_files, extinction_file, albedo_file, gravity, name, output_path, scale_ext=1):
	"""
	Top level code to take CH's output, convert to necessary input files for abledo code, 
	and run the albedo code. 

	Parameters
	----------
	do_clouds : bool 
		Turn on and off clouds
	three_chem_files : list of str 
		CH should provide you with three chemistry files. Put them here in list form in any order.
	extinction_file : str 
		CH should also provide you with an extinction file from her cloud code. Reference that here. 
		If you want to run a cloud free example, then just specify: None
	albedo_file : str 
		CH should also provide you with an albedo file to describe the scattering. If you want this to 
		be zero then you should also specify: None 
	gravity : float 
		Gravity in units of m/s2
	output_path : str 
		Where you want things to be saved to.  
	scale_ext : float 
		Default=1, this scaling factor is purely for learning purposes. It scales the total extinction by: 
		new_extinction = scale_ext * old_extinction
	"""
	run = AlbedoCode(three_chem_files, extinction_file, albedo_file, gravity, output_path, scale_ext=scale_ext)
	#makes first input file
	run.make_input_pt()
	#make ind file for indexing opacities 
	run.calc_ind_file(pt_file=os.path.join(path_to_alb_code, 'run', 'input.pt'), 
		ind_file = os.path.join(path_to_alb_code, 'run','input.ind'))

	#make input.cld 
	if do_clouds:
		try:
			os.path.exists(extinction_file)
		except: 
			raise Exception("No Extinction_file Found, Running Cloud Free Case")
			
		run.make_input_cld()
		do_clouds = str(int(do_clouds))
		#run if cloud free is specified 
	else: 
		do_clouds = str(int(do_clouds))

	inputs = os.path.join(path_to_alb_code, 'run', 'int.params')
	input_bak = os.path.join(path_to_alb_code, 'run', 'int.params.bak')
	with open(input_bak, 'r') as input_file, open(inputs, 'w') as out_file:
		for line in input_file:
			if line.find('DO_CLOUD=') != -1:
				out_file.write('DO_CLOUD= '+str(int(do_clouds))+'\n')
			else:
				out_file.write(line)

	#now you are ready to build the process 
	#let's go to the albedo code 
	os.chdir(os.path.join(path_to_alb_code,'run'))
	subprocess.run(['./geom'],stdout=subprocess.DEVNULL)
	os.chdir(init)

	for i in os.listdir(os.path.join(path_to_alb_code,'run')):
		if i.find('00_Spectral_Albedos') == -1: 
			continue
		phase = i[i.rfind('_')+1:i.find('.')]
		alb = pd.read_csv(os.path.join(path_to_alb_code,'run',i), delim_whitespace=True,skiprows=3,usecols=[1,3])
		alb['GEOMALB'] = scisig.medfilt(alb['GEOMALB'],kernel_size=3)
		alb.to_csv(os.path.join(output_path,name+'_'+phase+'.dat'))
	return

class AlbedoCode():
	def __init__(self, three_chem_files, extinction_file, albedo_file, gravity, output_path, scale_ext=1):
		self.three_chem_files = three_chem_files
		self.extinction_file = extinction_file
		self.albedo_file = albedo_file
		self.gravity = gravity
		self.output_path = output_path
		self.scale_ext = scale_ext

	def make_input_pt(self):
		"""
		Functionality to make input.pt 
		1)grabs necessary chemistry from CH three fiels 
		2) converts/interpolates to what we need 
		"""
		#GET THE NEW CHEMISTRY FILES FROM CH  
		self.new1 = pd.read_csv(self.three_chem_files[0],skiprows=3,delim_whitespace=True) 
		new2 = pd.read_csv(self.three_chem_files[1],skiprows=3,delim_whitespace=True) 
		new3 = pd.read_csv(self.three_chem_files[2],skiprows=3,delim_whitespace=True) 
		self.standard_pt = pd.read_csv(os.path.join(path_to_alb_code,'run','standard.pt'), delim_whitespace=True, 
		                         skiprows=[1])

		#FIND THE MOLECULES WE NEED IN HER FILES
		int1 = set(list(self.standard_pt.keys())) & set(list(self.new1.keys()))
		int2 = set(list(self.standard_pt.keys())) & set(list(new2.keys()))
		int3 = set(list(self.standard_pt.keys())) & set(list(new3.keys()))

		#GET THE NEW TEMPERATURE PROFILE FROM HER
		newCH = np.interp(np.log10(self.standard_pt['P']), np.log10(self.new1['p'].values*1e-6), self.new1['T'])
		self.standard_pt['T'] = newCH

		#GET P, K, T IN ORDER TO CALCULATE THE NUMBER DENSITY 
		#WE WILL CONVERT CM^-3 TO MIXING RATIO
		p = (self.standard_pt['P'].values*u.bar).to(u.dyn/u.cm/u.cm).value #dyn/cm2
		self.k = (c.k_B).to(u.erg/u.K).value
		t = self.standard_pt['T'].values
		num_dens = p/self.k/t
		#START WITH CH'S CHEM FILE 1
		for i in int1:
		    if i == 'T':continue
		    #print(i)
		    newCH = np.interp(np.log10(self.standard_pt['P']), np.log10(self.new1['p'].values*1e-6), self.new1[i])
		    self.standard_pt[i] = newCH / num_dens #convert to mixing ratio

		#REPEAT FOR CH'S CHEM FILE 2
		for i in int2:
		    if i == 'T':continue
		    #print(i)
		    newCH = np.interp(np.log10(self.standard_pt['P']), np.log10(new2['p'].values*1e-6), new2[i])
		    self.standard_pt[i] = newCH / num_dens #convert to mixing ratio

		#REPEAT FOR CH'S CHEM FILE 3
		for i in int3:
		    if i == 'T':continue
		    #print(i)
		    newCH = np.interp(np.log10(self.standard_pt['P']), np.log10(new3['p'].values*1e-6), new3[i])
		    self.standard_pt[i] = newCH / num_dens #convert to mixing ratio

		#GET RID OF WEIRD INTERPOLATED BOUNDARIES
		self.standard_pt = self.standard_pt.loc[(self.standard_pt['P']> np.min(self.new1['p'].values*1e-6)) 
		                              & (self.standard_pt['P']<np.max(self.new1['p'].values*1e-6))]

		#PUT EVERYTHING IN OUR FORMATED PT FILE
		self.standard_pt.to_csv(os.path.join(path_to_alb_code, 'run','input.pt')
			,sep = ' ',float_format='%.4E',index=False)

		#LASTLY, STICK IN THE GRAVITY OF THE PLANET 
		f = open(os.path.join(path_to_alb_code, 'run','input.pt'), "r")
		contents = f.readlines()
		f.close()
		contents.insert(1, str(self.gravity)+' \n')
		f = open(os.path.join(path_to_alb_code, 'run','input.pt'), "w")
		contents = "".join(contents)
		f.write(contents)
		f.close()

		#edit program parmaeters to make sure layers are consistent
		inputs = os.path.join(path_to_alb_code,'prog_params')
		input_bak = os.path.join(path_to_alb_code,'prog_params.bak')
		with open(input_bak, 'r') as input_file, open(inputs, 'w') as out_file:
			for line in input_file:
				if line.find('(nlayer=60,nlevel=61)') != -1:
					out_file.write('      PARAMETER (nlayer={0},nlevel={1})\n'.format((self.standard_pt.shape[0]-1)  ,(self.standard_pt.shape[0]) ))
				else:
					out_file.write(line)


	def calc_ind_file(self, pt_file='input.pt', ind_file = 'input.ind'):
		'''
		This function contains the functionality for creating the index file for computing the midpoint 
		values of temperature and pressure. It also grabs the nearest p,t pair in the opacity files. 
		This is used to grab the right file from the freedman opacity grid

		Curerntly uses two hard coded files PTgrid1060.txt and PTgrid736.txt 
		Code depends on these two files being structured as col1: file number, col2: pressure, col3: temperature.

		Parameters
		----------
		pt_file : filename
			File name with first three columns [level, pressure level, temperature level]

		ind_file : filename
			(Optional) output file name to use as input for the albedo code
		'''

		#opcity files
		grid1060 = pd.read_csv(os.path.join(path_to_alb_code, 'run', 'PTgrid1060.txt'), 
			delim_whitespace=True, skiprows = 1, header=None, names=['fn', 'p', 't'], usecols=[0,1,2])
		grid736 =  pd.read_csv(os.path.join(path_to_alb_code, 'run', 'PTgrid736.txt'), 
			delim_whitespace=True, skiprows = 1, header=None, names=['fn', 'p', 't'], usecols=[0,1,2])
		
		#pt file from model
		pt =  pd.read_csv(pt_file, delim_whitespace=True, skiprows = 2, header=None, names=['lvl', 'p_lvl', 't_lvl'], usecols=[0,1,2])

		#evaluate t and p at midpoints of t,p grid
		t_lyr = 0.5*(np.array(pt['t_lvl'][0:-1]) + np.array(pt['t_lvl'][1:]))
		p_lyr = np.sqrt(np.array(pt['p_lvl'][0:-1]) * np.array(pt['p_lvl'][1:]))

		#create spatial tree to find best pt point
		tree1060 = cKDTree(grid1060[['p','t']])
		tree736 = cKDTree(grid736[['p','t']])

		#find nearest neighbors for each t,p pair in the mid point layers 
		indout = pd.DataFrame({'p':[], 't':[], 'x736f':[],'x736p':[],'x736t':[], 'x1060f':[],'x1060p':[],'x1060t':[]})

		for t, p in zip(t_lyr,p_lyr):
			dists1060, indexes1060 = tree1060.query(np.array([p,t]), k=1)
			dists736, indexes736 = tree736.query(np.array([p,t]), k=1)

			x736 = list(grid736.loc[indexes736])
			x1060 = list(grid1060.loc[indexes1060])

			indout = indout.append({'p':'{:.6e}'.format(p), 't':'{:.2f}'.format(t), 'x736f':'{:.0f}'.format(x736[0]),'x736p':'{:.6e}'.format(x736[1]),'x736t':'{:.2f}'.format(x736[2]), 'x1060f':'{:.0f}'.format(x1060[0]),'x1060p':'{:.6e}'.format(x1060[1]),'x1060t':'{:.2f}'.format(x1060[2])},ignore_index=True)
		indout.to_csv(ind_file,columns=['p','t','x736f','x736p','x736t','x1060f','x1060p','x1060t'],header=False, index=False, sep='\t')

	def make_input_cld(self):
		"""
		This function takes the extinction and albedo outputs from CH and converts them into 
		total optical depth per layer
		"""
		#GET CH'S ALBEDO AND EXTINCTION FILE, CLEAN UP, PUT INTO NICE DATAFRAMES
		skip = [0,1]
		skip += [i for i in range(3,322)]

		#FIRST EXTINCTION
		CH_wave = list(pd.read_csv(self.extinction_file, skiprows=skip,header=None,delim_whitespace=True).values[0])
		CH_ext = pd.read_csv(self.extinction_file,skiprows=4,delim_whitespace=True,header=None,names=['p']+CH_wave)
		CH_ext = CH_ext.loc[(CH_ext['p']> np.min(self.standard_pt['P'].values*1e6)) 
		                              & (CH_ext['p']<np.max(self.standard_pt['P'].values*1e6))]

		CH_ext_dropz= CH_ext.loc[CH_ext[0.216]!=0] #dropping zeros so I can do log10 interpolation
		CH_ext_dropz= CH_ext_dropz.drop_duplicates('p') #getting rid of duplicates in CH's grid

		#THEN ALBEDO
		CH_alb = pd.read_csv(self.albedo_file,skiprows=4,delim_whitespace=True,header=None,names=['p']+CH_wave)
		CH_alb = CH_alb.loc[(CH_alb['p']> np.min(self.standard_pt['P'].values*1e6)) 
		                              & (CH_alb['p']<np.max(self.standard_pt['P'].values*1e6))]
		CH_alb= CH_alb.drop_duplicates('p') #getting rid of duplicates in CH's grid, no need for zeros


		#then convert from cm2/g

		#FIRST CALCULATE TOTAL MMW OF ATMOSPHERE FROM CHEMISTRY FILES AS FUNCTION OF ALTITUDE 
		self.standard_pt.rename(columns={'TIO':'TiO'}, inplace=True)
		gas_weights = pd.DataFrame(get_weights(list(self.standard_pt.keys())[4:]),index=[0])
		mixingratios = self.standard_pt[sorted(list(self.standard_pt.keys())[4:])]
		mu_atm = pd.Series(np.dot(gas_weights.values.as_matrix(),mixingratios.values.as_matrix().transpose())[0], index=self.standard_pt['P'])
		#CONVERT TO GRAMS BC CH'S FILE IS CM2/G
		m_u = c.u.to(u.g) #c is astropy.constants.u for atomic mass constant
		mu_atm = mu_atm*m_u.value #grams

		#we define cloud parameters at midpoint of the pt grid, so this is slightly different grid from above for input.pt
		#with one less grid point.. so redefining everything here with _cld to signify the different grid
		self.standard_cld = pd.read_csv(os.path.join(path_to_alb_code,'run','log.eddy'), delim_whitespace=True)
		self.standard_cld = self.standard_cld.loc[(self.standard_cld['P(bar)']> np.min(self.standard_pt['P'].values)) 
		                              & (self.standard_cld['P(bar)']<np.max(self.standard_pt['P'].values))]

		#GET CH'S EXTINCTION IN CM2/G
		new_ext = pd.DataFrame(columns=CH_ext.keys()[1:], index = self.standard_cld['P(bar)'].values)
		#GET CH'S PRESSURE IN DYN/CM2
		p_cld = (self.standard_cld['P(bar)'].values*u.bar).to(u.dyn/u.cm/u.cm).value #dyn/cm2
		#GET CH'S TEMP IN K
		t_cld = np.interp(np.log10(self.standard_cld['P(bar)']), np.log10(self.new1['p'].values*1e-6), self.new1['T'])
		#GET OUR CALCULATED MMW IN GRAMS
		mu_atm_cld = 10**np.interp(np.log10(self.standard_cld['P(bar)']), np.log10(mu_atm.index), np.log10(mu_atm.values))

		#NOW GET THE ALTITUDE SO WE CAN CALCULATE DZ
		z = np.interp(np.log10(self.standard_pt['P']), np.log10(self.new1['p'].values*1e-6), self.new1['z'])

		#GET LENGTH OF EACH LAYER IN CM
		dz = np.abs(np.diff(z)) #cm 

		#interpolate extinction to PT and convert to optical depth
		for i in new_ext.keys():
		    #interpolate
		    ext_new_grid = 10**np.interp(np.log10(self.standard_cld['P(bar)']), np.log10(CH_ext_dropz['p'].values*1e-6), np.log10(CH_ext_dropz[i]))

		    #convert to from extinction in cm2/g to optical detph
		    rho = mu_atm_cld * p_cld /t_cld /self.k #all cgs units
		    new_ext[i] = ext_new_grid * rho * dz#cm2/g * g/cm3 * cm
		    
		#make sure pressures with zeros are zerod out again 
		zeros_start = CH_ext.loc[CH_ext[0.216]==0]['p'].iloc[0]*1e-6
		new_ext.loc[new_ext.index>=zeros_start] = 0 

		#FINALLY GET SCATTERING ALBEDO
		new_alb = pd.DataFrame(columns=CH_ext.keys()[1:], index = self.standard_cld['P(bar)'].values)
		print(CH_ext.keys()[1:])
		for i in new_alb.keys():
			alb_new_grid = 10**np.interp(np.log10(self.standard_cld['P(bar)']), np.log10(CH_alb['p'].values*1e-6), np.log10(CH_alb[i]))
			new_alb[i] = alb_new_grid

		#Finally, interpolate wavelength to create input.cld
		#now all that's left is to reinterp wavelength and build final cld output 
		#starting with new_ext and new_alb 

		Lvl = np.concatenate([[i]*196 for i in range(1,51)])
		Wv = list(range(1,197))*50
		#g0=0 is isotropic scattering
		final_cld = pd.DataFrame(columns=[ 'Lvl-', 'Wv', 'Opd',  'g0',  'w0' ,  'sigma'],
		                        index = range(len(Lvl))) #dont need sigma so adding zeros
		final_cld['Lvl-'] = Lvl
		final_cld['Wv'] = Wv
		final_cld['sigma'] = np.array(Wv)*0 #cross sections only used in transit code
		#isotropic for now
		final_cld['g0'] = np.array(Wv)*0
		wave_196 = sorted(np.concatenate(pd.read_csv(os.path.join(path_to_alb_code,'run','wave_EGP.dat'), delim_whitespace=True, header=None,usecols=[1]).values))

		#lets interpolate each level 
		for i in range(50): 
		    opd = 10**np.interp(np.log10(wave_196), np.log10(np.array([float(ii) for ii in new_ext.keys()])), np.log10(new_ext.iloc[i].values))
		    alb = 10**np.interp(np.log10(wave_196), np.log10(np.array([float(ii) for ii in new_alb.keys()])), np.log10(new_alb.iloc[i].values))
		    final_cld['Opd'].iloc[i*196:(i+1)*196] = opd*self.scale_ext
		    final_cld['w0'].iloc[i*196:(i+1)*196] = alb
		final_cld=final_cld.fillna(0)

		#FINALLY build cld file
		final_cld.to_csv(os.path.join(path_to_alb_code, 'run', 'input.cld'),sep = ' ',float_format='%.4E',index=False)



def get_weights(molecule):
    """Author:Natasha Batalha
    Automatically gets mean molecular weights of any molecule. Requires that 
    user inputs case sensitive molecules i.e. TiO instead of TIO. 

    Parameters
    ----------
    molecule : str or list
        any molecule string e.g. "H2O", "CH4" etc "TiO" or ["H2O", 'CH4']

    Returns
    -------
    dict 
        molecule as keys with molecular weights as value
    """
    weights = {}
    if not isinstance(molecule,list):
        molecule = [molecule]

    for i in molecule:
        if (i.find('-')!=-1) or (i.find('+')!=-1): 
            add = i[np.max([i.find('+'), i.find('-')])]
            i = i[0:np.max([i.find('+'), i.find('-')])]
        else:
            add = ''
        molecule_list = []
        for j in range(0,len(i)):
            try:
                molecule_list += [float(i[j])]
            except: 
                if i[j].isupper(): molecule_list += [i[j]] 
                elif i[j].islower(): molecule_list[j-1] =  molecule_list[j-1] + i[j]
        totmass=0
        for j in range(0,len(molecule_list)): 

            if isinstance(molecule_list[j],str):
                elem = ele[molecule_list[j]]
                try:
                    num = float(molecule_list[j+1])
                except: 
                    num = 1 
                totmass += elem.mass * num
        weights[i+add] = totmass
    return weights

def jupiter(output_path, do_clouds=False, fsed=None , name=None):
	"""
	This runs a classic Jupiter example with different specifications for clouds. 

	Parameters
	----------
	output_path : str 
		path to output directory 
	do_clouds : bool 
		Turns on and off clouds (Default=False)
	fsed : int 
		Sedimentation efficiency. Options are [6, 3, 1]. Only used if do_clouds=True
	name : str 
		If you want to name the file something special: default = jupiter_nocloud, jupiter_cloud_6fsed

	"""
	#opacity index file
	shutil.copy2(os.path.join(path_to_alb_code, 'inputs','jupiter','jupiter.ind'), 
		os.path.join(path_to_alb_code, 'run','input.ind'))
	#pt chem file
	shutil.copy2(os.path.join(path_to_alb_code, 'inputs','jupiter','jupiter.pt'), 
		os.path.join(path_to_alb_code, 'run','input.pt'))
	#make input.cld 
	if do_clouds:
	#pt chem file
		if fsed==None:
			f = int(6)
			print('Selecting Default Sedimentation efficiency= 6.0')
		elif int(fsed) in [6,3,1]:
			f = int(fsed)
		else: 
			raise Exception ('Only fsed=6,3,1 available on git, Email natasha.e.batalha@gmail for other profiles')

		if name==None: name = 'jupiter_cloud_'+str(f)+'fsed'
		shutil.copy2(os.path.join(path_to_alb_code, 'inputs','jupiter','f'+str(f)+'.cld'), 
			os.path.join(path_to_alb_code, 'run','input.cld'))

		do_clouds = str(int(do_clouds))

	else: 
		do_clouds = str(int(do_clouds))
		if name==None: name = 'jupiter_nocloud'


	inputs = os.path.join(path_to_alb_code, 'run', 'int.params')
	input_bak = os.path.join(path_to_alb_code, 'run', 'int.params.bak')
	with open(input_bak, 'r') as input_file, open(inputs, 'w') as out_file:
		for line in input_file:
			if line.find('DO_CLOUD=') != -1:
				out_file.write('DO_CLOUD= '+str(int(do_clouds))+'\n')
			else:
				out_file.write(line)

	#now you are ready to build the process 
	#let's go to the albedo code 
	os.chdir(os.path.join(path_to_alb_code,'run'))
	subprocess.run(['./geom'],stdout=subprocess.DEVNULL)
	os.chdir(init)

	for i in os.listdir(os.path.join(path_to_alb_code,'run')):
		if i.find('00_Spectral_Albedos') == -1: 
			continue
		phase = i[i.rfind('_')+1:i.find('.')]
		alb = pd.read_csv(os.path.join(path_to_alb_code,'run',i), delim_whitespace=True,skiprows=3,usecols=[1,3])
		alb['GEOMALB'] = scisig.medfilt(alb['GEOMALB'],kernel_size=3)
		alb.to_csv(os.path.join(output_path,name+'_'+phase+'.dat'))

