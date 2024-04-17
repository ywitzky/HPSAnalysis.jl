# Apesar de usar o NVE a energia nao esta sendo conservada, provavelmente pq o momento total nao e zerado pelo integrador do NVE. Pesquisar como utilzar o comando de zerar momentos.

import itertools
import math
import numpy as np

import gsd.hoomd
#import hoomd
#import hoomd.md
#import hoomd.pair_plugin
#import hoomd.plugin_insertion

import numpy

import matplotlib.pyplot as plt
import time

#Aminoacids = np.loadtxt("Aminoacid_Parameters_Urry.dat", skiprows=1, dtype = "str")
charge_conversion_factor = 0.02682831841 * np.sqrt(80.0)
bondLength = 3.8


SimName = "ALA40CYS40"
volumeFraction = 0.001
chainDistanceXY = 10.0

eps = 0.2*4.18
kappa = 10.0

dih_epsA = -2.7
dih_epsC = -0.15

kb = 0.00831446262 ### KJ/mol
T = 300.0#K
kT = kb*T
beta = 1.0/kT

rcut = 30.0

Nsteps = 10**7
Nout = 10**4
Ninsert = 10
timestep = 0.025 


#idpStructure = np.loadtxt("ALA40CYS40.dat", delimiter = ",", dtype = "str")
snapshot = gsd.hoomd.Snapshot()

#print(idpStructure)

Nchains = 4
Nbeads = 10#len(idpStructure)
Nparticles = Nchains*Nbeads
nlattice = int(np.ceil(np.sqrt(Nchains)))

chainLength = bondLength * (Nbeads-1)
chainVolume = Nbeads * 0.55**3 * np.pi / 6.0
totalVolume = Nchains * chainVolume
L = (totalVolume / volumeFraction)**(1/3.) * 8
L = 5*L

Lx = 150.0
Ly = 150.0
Lz = 2800.0
print(L)
indexBead = 0

#Create chain

ChainMasses = np.zeros(Nbeads)
ChainCharges = np.zeros(Nbeads)
ChainTypes = np.zeros(Nbeads)
ChainPositions = np.zeros((Nbeads,3))
ChainBonds = np.zeros((Nbeads-1,2), dtype=int)
ChainAngles = np.zeros((Nbeads-2,3), dtype=int)
ChainDihedrals = np.zeros((Nbeads-3,4), dtype=int)

ChainDihedralTypes = np.zeros((Nbeads-3), dtype=int)


TypesList = []
TypesLookup = []


for i in range(Nbeads):
	currentAminoacid = idpStructure[i]
	if currentAminoacid not in TypesLookup:
		TypesLookup.append(currentAminoacid)
		TypesList.append(len(TypesList))

TypesLookupArray = np.array(TypesLookup)

for i in range(Nbeads):
	currentAminoacid = idpStructure[i]
	
	#if currentAminoacid not in TypesLookup:
	#	TypesLookup.append(currentAminoacid)
	#	TypesList.append(len(TypesList))
		
	currentType = np.where(Aminoacids == currentAminoacid)	
	Params = Aminoacids[currentType[0][0]]
	
	ChainMasses[i] = float(Params[1])
	ChainCharges[i] = float(Params[2])/charge_conversion_factor
	TypeIndex = np.where(TypesLookupArray == currentAminoacid)
	ChainTypes[i] = TypeIndex[0][0]
	
	ChainPositions[i,:] = [-float(nlattice)/2.0*chainDistanceXY, -float(nlattice)/2.0*chainDistanceXY, bondLength*i - chainLength/2.0]
	
	#print(Aminoacids[currentType[0][0]],TypesList,TypesLookupArray,TypeIndex[0][0],[0.0, 0.0, bondLength*i])
	
for i in range(Nbeads-1):
	ChainBonds[i,:] = [i,i+1]
	
for i in range(Nbeads-2):
	ChainAngles[i,:] = [i,i+1,i+2]
	
for i in range(Nbeads-3):
	ChainDihedrals[i,:] = [i,i+1,i+2,i+3]

dihTypesList = []
dihTypesLookup = []

#First dihedral
dihResidues = ["O",idpStructure[0],idpStructure[3],idpStructure[4]]
dihTypeName = str(dihResidues[0])+"-"+str(dihResidues[1])+"-"+str(dihResidues[2])+"-"+str(dihResidues[3])
ChainDihedralTypes[0] = 0
dihTypesLookup.append(dihTypeName)
dihTypesList.append(len(dihTypesList))


print(dihResidues,dihTypeName)
for i in range(1,Nbeads-4):
	dihResidues = [idpStructure[i-1],idpStructure[i],idpStructure[i+3],idpStructure[i+4]]
	dihTypeName = str(dihResidues[0])+"-"+str(dihResidues[1])+"-"+str(dihResidues[2])+"-"+str(dihResidues[3])
	
	if dihTypeName not in dihTypesLookup:
		dihTypesLookup.append(dihTypeName)
		dihTypesList.append(len(dihTypesList))


#Last Dihedral
dihResidues = [idpStructure[Nbeads-5],idpStructure[Nbeads-4],idpStructure[Nbeads-2],"O"]
dihTypeName = str(dihResidues[0])+"-"+str(dihResidues[1])+"-"+str(dihResidues[2])+"-"+str(dihResidues[3])
dihTypesLookup.append(dihTypeName)
dihTypesList.append(len(dihTypesList))
print(dihTypesLookup,dihTypesList,ChainDihedralTypes)
ChainDihedralTypes[-1] = dihTypesList[-1]

dihTypesLookupArray = np.array(dihTypesLookup)

DihedralEps = np.zeros(len(dihTypesLookup))
for i in range(len(dihTypesLookup)):
	epsCumulant = 0.0
	dihName = dihTypesLookup[i]
	resName = ""
	print(dihName)
	
	for j in range(len(dihName)):
		if dihName[j] != "-":
			resName += dihName[j]
		else:
			print(resName)
			if resName == "A":
				epsCumulant += dih_epsA
			if resName == "C":
				epsCumulant += dih_epsC			
			resName = ""
	
	print(resName)
	if resName == "A":
		epsCumulant += dih_epsA
	if resName == "C":
		epsCumulant += dih_epsC
		
	DihedralEps[i] = epsCumulant/4.0

print(DihedralEps)
for i in range(1,Nbeads-4):
	dihResidues = [idpStructure[i-1],idpStructure[i],idpStructure[i+3],idpStructure[i+4]]
	dihTypeName = str(dihResidues[0])+"-"+str(dihResidues[1])+"-"+str(dihResidues[2])+"-"+str(dihResidues[3])

	dihTypeIndex = np.where(dihTypesLookupArray == dihTypeName)	
	print(dihTypeIndex,dihTypeName,dihTypesLookupArray)
	ChainDihedralTypes[i] = dihTypeIndex[0][0]


print(dihTypesLookup,dihTypesList,ChainDihedralTypes)

#Replicate Chain

InputPositions = np.zeros((Nparticles,3))
InputMasses = np.zeros(Nparticles)
InputCharges = np.zeros(Nparticles)
InputTypes = np.zeros(Nparticles)
InputBonds = np.zeros((Nchains*(Nbeads-1),2), dtype=int)
InputAngles = np.zeros((Nchains*(Nbeads-2),3), dtype=int)
InputDihedrals = np.zeros((Nchains*(Nbeads-3),4), dtype=int)
InputDihedralIDs = np.zeros((Nchains*(Nbeads-3)), dtype=int)

for i in range(Nchains):
	for j in range(Nbeads):
	
		InputPositions[Nbeads*i+j,0] = ChainPositions[j,0] + (i%nlattice) * chainDistanceXY + np.random.uniform(low=-0.01, high=0.01)
		InputPositions[Nbeads*i+j,1] = ChainPositions[j,1] + (i/nlattice) * chainDistanceXY + np.random.uniform(low=-0.01, high=0.01)
		InputPositions[Nbeads*i+j,2] = ChainPositions[j,2] + np.random.uniform(low=-0.01, high=0.01)
		
		InputMasses[Nbeads*i+j] = ChainMasses[j]
		InputCharges[Nbeads*i+j] = ChainCharges[j]
		InputTypes[Nbeads*i+j] = ChainTypes[j]
		
for i in range(Nchains):
	for j in range(Nbeads-1):
		InputBonds[(Nbeads-1)*i+j,0] =  ChainBonds[j,0] + i*(Nbeads)
		InputBonds[(Nbeads-1)*i+j,1] =  ChainBonds[j,1] + i*(Nbeads)
		
for i in range(Nchains):
	for j in range(Nbeads-2):
		InputAngles[(Nbeads-2)*i+j,0] =  ChainAngles[j,0] + i*(Nbeads)
		InputAngles[(Nbeads-2)*i+j,1] =  ChainAngles[j,1] + i*(Nbeads)
		InputAngles[(Nbeads-2)*i+j,2] =  ChainAngles[j,2] + i*(Nbeads)

for i in range(Nchains):
	for j in range(Nbeads-3):
		InputDihedrals[(Nbeads-3)*i+j,0] =  ChainDihedrals[j,0] + i*(Nbeads)
		InputDihedrals[(Nbeads-3)*i+j,1] =  ChainDihedrals[j,1] + i*(Nbeads)
		InputDihedrals[(Nbeads-3)*i+j,2] =  ChainDihedrals[j,2] + i*(Nbeads)
		InputDihedrals[(Nbeads-3)*i+j,3] =  ChainDihedrals[j,3] + i*(Nbeads)
		InputDihedralIDs[(Nbeads-3)*i+j] =  ChainDihedralTypes[j]

		#print(Nparticles,InputDihedrals[(Nbeads-3)*i+j,:])			
print(TypesLookup)
	
#Create particles
snapshot.particles.N = Nparticles
snapshot.particles.position = InputPositions
snapshot.particles.types = TypesLookup
snapshot.particles.typeid = InputTypes

#print(snapshot.particles.N)
A = [1.0] * Nbeads
for i in range(Nchains-1):
	a = numpy.random.uniform(0.0,1.0)
	if a<10.5:
		A += [1.0] * Nbeads
	else:
		A += [1.0] * Nbeads
		
		
snapshot.particles.diameter	= A
snapshot.configuration.box = [L, L, L, 0, 0, 0]


snapshot.particles.mass = InputMasses
snapshot.particles.charge = InputCharges

# Connect particles with bonds.
snapshot.bonds.N = (Nbeads-1) * Nchains
snapshot.bonds.types = ['O-O']
snapshot.bonds.typeid = np.array([0] * snapshot.bonds.N)
snapshot.bonds.group = InputBonds

# Create Angles
snapshot.angles.N = (Nbeads-2) * Nchains
snapshot.angles.types = ['O-O-O']
snapshot.angles.typeid = np.array([0] * snapshot.angles.N)
snapshot.angles.group = InputAngles

# Create Angles
snapshot.dihedrals.N = (Nbeads-3) * Nchains
snapshot.dihedrals.types = dihTypesLookup#['O-O-O-O']
snapshot.dihedrals.typeid = InputDihedralIDs#np.array([0] * snapshot.dihedrals.N)
snapshot.dihedrals.group = InputDihedrals

#print

with gsd.hoomd.open(name=SimName + "_StartConfiguration.gsd", mode='wb') as f:
    f.append(snapshot)
    
    
# Apply the harmonic potential on the bonds.
harmonic = hoomd.md.bond.Harmonic()
harmonic.params['O-O'] = dict(k=10000, r0=bondLength)


# Apply the Angle potential on the bonds.
N_angle_points = 500
angle_gamma = 0.1
angle_eps = 4.3
angle_phi0a = 1.60
angle_phi0b = 2.27
angle_ka = 106.4
angle_kb = 26.3
angle_x = np.linspace(-0*np.pi,np.pi,N_angle_points)
angle_dx = angle_x[1] - angle_x[0]

def GetAngles(x,ka,kb,phi0a,phi0b,eps,gamma):
	return - (1.0/gamma) * np.log (  np.exp( - gamma*(ka * (x - phi0a)**2 + eps) ) + np.exp( - gamma * kb * (x - phi0b)**2)  )

y_anglePotential = GetAngles(angle_x,angle_ka,angle_kb,angle_phi0a,angle_phi0b,angle_eps,angle_gamma)
y_anglePotential *= 4.18
angle_tau = -np.gradient(y_anglePotential,angle_dx)

anglePotential = hoomd.md.angle.Table(N_angle_points)
anglePotential.params['O-O-O'] = dict(U=y_anglePotential, tau=angle_tau)

#anglePotential = hoomd.md.angle.Harmonic()
#anglePotential.params['O-O-O'] = dict(k=3000000.0, t0=1.7)

#Apply the Dihedral Potential on the beads
N_dih_points = 500
dih_eps = -2.70
dih_ka1 = 11.4
dih_ka2 = 0.15
dih_kb1 = 1.80
dih_kb2 = 0.65
dih_phia1 = 0.90
dih_phia2 = 1.02
dih_phib1 = -1.55
dih_phib2 = -2.50
dih_e0 = 0.27
dih_e1 = 0.14
dih_e2 = 0.40

dih_x = np.linspace(-np.pi,np.pi,N_dih_points)
dih_dx = dih_x[1] - dih_x[0]

def Table1(x,k,phi0):
	return np.exp( - k * (x - phi0)**2 )
def Table2(x,k,phi0,e1):
	return np.exp( - k * (x - phi0)**2 + e1 ) + np.exp( - k * (x - phi0 - 2*np.pi)**2 + e1 )
def Table3(x,ka,kb,phi0a,phi0b,e0,e2):
	return np.exp( - ka * (x - phi0a)**4 + e0 ) + np.exp( - ka * (x - phi0a + 2*np.pi)**4 + e0 ) + np.exp( - kb * (x - phi0b)**4 + e2 ) + np.exp( - kb * (x - phi0b - 2*np.pi)**4 + e2 )
def GetOPLS(x,k1,k2,k3,k4):
	return 0.5 * ( k1*(1+np.cos(x)) + k2*(1-np.cos(2*x)) + k3*(1+np.cos(3*x)) + k4*(1-np.cos(4*x)) )
def GetTauOPLS(x,k1,k2,k3,k4):
	return -0.5 * ( k1*(-np.sin(x)) + k2*(+2.0*np.sin(2*x)) + k3*(-3.0*np.sin(3*x)) + k4*(+4.0*np.sin(4*x)) )
y1 = Table1(dih_x,dih_ka1,dih_phia1)
y2 = Table2(dih_x,dih_kb1,dih_phib1,dih_e1)
y3 = Table3(dih_x,dih_ka2,dih_kb2,dih_phia2,dih_phib2,dih_e0,dih_e2)
#y1 *= np.exp(4.18)
#y2 *= np.exp(4.18)
#y3 *= np.exp(4.18)
y_dihedralPotential = (y1*np.exp(-dih_eps) + y2*np.exp(dih_eps) + y3)
y_dihedralPotential = -np.log(y_dihedralPotential)
#y_dihedralPotential *= 4.18
#y_dihedralPotential = GetOPLS(dih_x,1.0,1.0,1.0,1.0)
tau_dih = -np.gradient(y_dihedralPotential,dih_dx)
tau1 = np.gradient(y1,dih_dx)
tau2 = np.gradient(y2,dih_dx)
tau3 = np.gradient(y3,dih_dx)
 
pott = open("Output_Potentials2.dat","w")
for i in range(N_angle_points):
	pott.write("%f %f %f %f %f %f\n"%(angle_x[i],y_anglePotential[i],angle_tau[i],dih_x[i],y_dihedralPotential[i],tau_dih[i]))
pott.close()	
	
#dihedralPotential = hoomd.md.dihedral.Table(N_dih_points)
#dihedralPotential.params['O-O-O-O'] = dict(U=y_dihedralPotential, tau=tau_dih)

dihedralPotential = hoomd.plugin_insertion.dihedral.DihedralSS(width = N_dih_points, u1 = y1, u2 = y2, u3 = y3, tau1 = tau1, tau2 = tau2, tau3 = tau3, prefactor=4.18)
for i in range(len(dihTypesLookup)):
	dihedralPotential.params[dihTypesLookup[i]] = dict(eps = DihedralEps[i])
	print("dihedralPotential.params["+dihTypesLookup[i]+"] = dict(eps = " + str(DihedralEps[i])+ "))")

#dihedralPotential = hoomd.md.dihedral.OPLS()
#dihedralPotential.params['O-O-O-O'] = dict(k1=1.0, k2=1.0, k3=1.0, k4=1.0)

#print(len(y_dihedralPotential),len(tau_dih),len(y_anglePotential),len(angle_tau))
#print((y_dihedralPotential),(tau_dih),(y_anglePotential),(angle_tau))


# Create Simulation
gpu = hoomd.device.CPU()
sim = hoomd.Simulation(device=gpu, seed=18657)
sim.create_state_from_gsd(filename=SimName + "_StartConfiguration.gsd")








#MD integrators

# dt
integrator = hoomd.md.Integrator(dt=timestep) 

integrator.forces.append(harmonic)
integrator.forces.append(anglePotential)
integrator.forces.append(dihedralPotential)


# neighbor list
cell = hoomd.md.nlist.Cell(buffer=0.4)   

# Interactions
hps = hoomd.plugin_insertion.pair.ExamplePair(nlist=cell)
#print(len(TypesList))
for i in range(len(TypesList)):
	for j in range(i,len(TypesList)):
	
		PairTypes = [snapshot.particles.types[i], snapshot.particles.types[j]]
		A = [np.where(Aminoacids == PairTypes[0])[0][0], np.where(Aminoacids == PairTypes[1])[0][0]]
		print(i,j,PairTypes,A)
		Params_i = Aminoacids[A[0]]
		Params_j = Aminoacids[A[1]]
		#print(i,j,PairTypes,A,Params_i,Params_j)
	
		hps.params[(snapshot.particles.types[i], snapshot.particles.types[j])] = dict(epsilon=eps, sigma=0.5*(float(Params_i[3]) + float(Params_j[3])), kappa=kappa, lamb=0.5*(float(Params_i[4]) + float(Params_j[4])))
		hps.r_cut[(snapshot.particles.types[i], snapshot.particles.types[j])] = rcut# 2*(rcut_max + r_buff) < L
	

integrator.forces.append(hps)


# NVE
#nvt = hoomd.md.methods.NVT(kT=kT, filter=hoomd.filter.All(), tau=1.0)
nvt = hoomd.md.methods.Langevin(filter=hoomd.filter.All(), kT=kT, alpha=0.01)
integrator.methods.append(nvt)

# Assign integrator to simulation
sim.operations.integrator = integrator

# Set initial velocities
sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=kT)


# Computes
thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
sim.operations.computes.append(thermodynamic_properties)
sim.run(0)

class DiameterComposition():

    def __init__(self, sim):
        self.sim = sim
        
    @property
    def composition(self):
    	snapshot = self.sim.state.get_snapshot()
    	return numpy.mean(snapshot.particles.diameter)


#Logger
logger = hoomd.logging.Logger()
logger.add(thermodynamic_properties)

diametercomp = DiameterComposition(sim)
logger[('DiameterComposition', 'composition')] = (diametercomp, 'composition', 'scalar')

print("Before setting updater\n")


molecules = hoomd.plugin_insertion.molecule.Molecules(hoomd.filter.All())
##forcecompute = hoomd._hoomd.ForceCompute(hoomd._hoomd.SystemDefinition)
insert_chain = hoomd.plugin_insertion.update.SetDiametersToOne(hoomd.trigger.Periodic(Ninsert,5), molecules, hps)

insert_chain.beta = beta
insert_chain.mu  =  68.0/(1.0*Nbeads*beta)#1.201*46.8885501668/(1.0*Nbeads)

for i in range(len(TypesList)):
	for j in range(i,len(TypesList)):
	
		PairTypes = [snapshot.particles.types[i], snapshot.particles.types[j]]
		A = [np.where(Aminoacids == PairTypes[0])[0][0], np.where(Aminoacids == PairTypes[1])[0][0]]
		Params_i = Aminoacids[A[0]]
		Params_j = Aminoacids[A[1]]
		print(i,j,PairTypes,A,Params_i,Params_j)
	
		insert_chain.params[(snapshot.particles.types[i], snapshot.particles.types[j])] = dict(epsilon=eps, sigma=0.5*(float(Params_i[3]) + float(Params_j[3])), kappa=kappa, lamb=0.5*(float(Params_i[4]) + float(Params_j[4])))
		insert_chain.r_cut[(snapshot.particles.types[i], snapshot.particles.types[j])] = rcut# 2*(rcut_max + r_buff) < L


#insert_chain.params[('R', 'R')] = dict(epsilon=eps, sigma=sigma, kappa=kappa, lamb=l)
#insert_chain.r_cut[('R', 'R')] = rcut# 2*(rcut_max + r_buff) < L

#insert_chain.params[('R', 'V')] = dict(epsilon=eps, sigma=sigma, kappa=kappa, lamb=l)
#insert_chain.r_cut[('R', 'V')] = rcut# 2*(rcut_max + r_buff) < L

#insert_chain.params[('V', 'V')] = dict(epsilon=eps, sigma=sigma, kappa=kappa, lamb=l)
#insert_chain.r_cut[('V', 'V')] = rcut# 2*(rcut_max + r_buff) < L

insert_chain.bond_params = [10000,bondLength]

#print(insert_chain.params[('R', 'R')])
#print(insert_chain.params[('R', 'V')])
#print(insert_chain.params[('V', 'V')])

teleport_chain = hoomd.plugin_insertion.update.ChainTeleporting(hoomd.trigger.Periodic(Ninsert,0), molecules)
teleport_chain.beta = beta

sim.operations.computes.append(molecules)
#sim.operations += insert_chain
#sim.operations += teleport_chain


print("After setting updater\n")

#table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(period=5000),logger=logger)
#sim.operations.writers.append(table)


# Write gsd file
gsd_writer = hoomd.write.GSD(trigger=hoomd.trigger.Periodic(Nout),filename='log_'+ SimName + '.gsd',filter=hoomd.filter.All(),mode='wb',dynamic=['attribute','property'])
sim.operations.writers.append(gsd_writer)
gsd_writer.log = logger

#print(integrator.forces[2].params.param_dict)
print("Before running simulation\n")
print("#Type timestep molecule Wratio DrawNumber Accept\n")
sim.run(Nsteps)
print("After running simulation\n")


# Write outputs
thermo = open("thermodynamics_"+SimName+".dat","w")
thermo.write("#t diameter pe ke p T\n")

traj = gsd.hoomd.open('log_'+ SimName + '.gsd', 'rb')
print(len(traj))


traj = gsd.hoomd.open('log_'+ SimName + '.gsd', 'rb')
print(len(traj))


for frame in traj:
	thermo.write("%f %f %f %f %f %f\n"%(frame.configuration.step,frame.log['DiameterComposition/composition'],frame.log['md/compute/ThermodynamicQuantities/potential_energy'][0],frame.log['md/compute/ThermodynamicQuantities/kinetic_energy'][0],frame.log['md/compute/ThermodynamicQuantities/pressure'][0],frame.log['md/compute/ThermodynamicQuantities/kinetic_temperature'][0]))


thermo.close()

print("END\n")






"""
SimName = "FR_TestNewHPS"

Aminoacids = np.loadtxt("Aminoacid_Parameters.dat", skiprows=1, dtype = "str")

Nbeads = 30
Nchains = 529
Concentration = 0.001 # 1000 micro mol per liter
Prefix = str(Nchains) + "Chains" + str(Nbeads) + "Beads_"

InputPositions = numpy.loadtxt(Prefix + "Positions.dat")
InputBonds     = numpy.loadtxt(Prefix + "Bonds.dat")
InputAngles    = numpy.loadtxt(Prefix + "Angles.dat")
InputTypes     = numpy.loadtxt(Prefix + "ParticleTypes.dat",dtype='int')
InputBondTypes = numpy.loadtxt(Prefix + "BondTypes.dat")

charge_conversion_factor = 0.02682831841 * np.sqrt(80.0)

eps = 0.2*4.18
kappa = 10.0
bond_length = 3.8

kb = 0.00831446262
T = 300.0

kT = kb*T
beta = 1.0/kT


rcut = 30.0

Nsteps = 10**8
Nout = 100000
Ninsert = 10
timestep = 0.05

N_particles = len(InputTypes)
L = 433.377 #10000.0/6.022 * (1.0*Nchains/Concentration)**(1.0/3.0)

snapshot = gsd.hoomd.Snapshot()
	

#Create particles
snapshot.particles.N = N_particles
snapshot.particles.position = InputPositions
snapshot.particles.types = ['F','R']
snapshot.particles.typeid = InputTypes
#snapshot.particles.mass = [mass] * N_particles
#snapshot.particles.charge = [q] * N_particles

A = [0.0] * Nbeads
for i in range(Nchains-1):
	a = numpy.random.uniform(0.0,1.0)
	if a<1.5:
		A += [0.0] * Nbeads
	else:
		A += [1.0] * Nbeads
		

print(len(A),N_particles)
		
snapshot.particles.diameter	= A#[1.0] * N_particles
#snapshot.particles.body = [0] * N_particles
snapshot.configuration.box = [L, L, L, 0, 0, 0]
#snapshot.configuration.box = [2*L, 2*L, 2*L, L, L, L]

InputMasses = np.zeros(N_particles)
InputCharges = np.zeros(N_particles)

print(L)
for i in range(N_particles):
	AminoacidType = snapshot.particles.types[InputTypes[i]]
	A = np.where(Aminoacids == AminoacidType)
	Params = Aminoacids[A[0][0]]
	InputMasses[i] = float(Params[1])
	InputCharges[i] = float(Params[2])/charge_conversion_factor
	#print(i,AminoacidType,Params)

#print(InputMasses,InputCharges)
snapshot.particles.mass = InputMasses
snapshot.particles.charge = InputCharges

# Connect particles with bonds.
snapshot.bonds.N = len(InputBonds)
snapshot.bonds.types = ['R-R','R-V','V-V']
snapshot.bonds.typeid = InputBondTypes
snapshot.bonds.group = InputBonds

# Connect particles with angles.
#snapshot.angles.N = len(InputAngles)
#snapshot.angles.types = ['A-A-A']
#snapshot.angles.typeid = [0] * len(InputAngles)
#snapshot.angles.group = InputAngles




with gsd.hoomd.open(name="4Chains" + str(Nbeads) + "Beads_Insertion.gsd", mode='wb') as f:
    f.append(snapshot)

# Apply the harmonic potential on the bonds.
harmonic = hoomd.md.bond.Harmonic()
harmonic.params['R-R'] = dict(k=10000, r0=bond_length)
harmonic.params['R-V'] = dict(k=10000, r0=bond_length)
harmonic.params['V-V'] = dict(k=10000, r0=bond_length)


# Apply the harmonic potential on the angles.
#angles = hoomd.md.angle.Harmonic()
#angles.params['A-A-A'] = dict(k=3.0, t0=numpy.pi)
    
    
# Create Simulation
gpu = hoomd.device.GPU()
sim = hoomd.Simulation(device=gpu, seed=18657)
sim.create_state_from_gsd(filename="4Chains" + str(Nbeads) + "Beads_Insertion.gsd")


#MD integrators

# dt
integrator = hoomd.md.Integrator(dt=timestep) 

integrator.forces.append(harmonic)
#integrator.forces.append(angles)


# neighbor list
cell = hoomd.md.nlist.Cell(buffer=0.4)   

# Interactions
lj = hoomd.plugin_insertion.pair.ExamplePair(nlist=cell)

for i in range(2):
	for j in range(i,2):
	
		PairTypes = [snapshot.particles.types[i], snapshot.particles.types[j]]
		A = [np.where(Aminoacids == PairTypes[0])[0][0], np.where(Aminoacids == PairTypes[1])[0][0]]
		Params_i = Aminoacids[A[i]]
		Params_j = Aminoacids[A[j]]
		#print(i,j,PairTypes,A,Params_i,Params_j)
	
		lj.params[(snapshot.particles.types[i], snapshot.particles.types[j])] = dict(epsilon=eps, sigma=0.5*(float(Params_i[3]) + float(Params_j[3])), kappa=kappa, lamb=0.5*(float(Params_i[4]) + float(Params_j[4])))
		lj.r_cut[(snapshot.particles.types[i], snapshot.particles.types[j])] = rcut# 2*(rcut_max + r_buff) < L
	

integrator.forces.append(lj)


# NVE
#nvt = hoomd.md.methods.NVT(kT=kT, filter=hoomd.filter.All(), tau=1.0)
nvt = hoomd.md.methods.Langevin(filter=hoomd.filter.All(), kT=kT, alpha=0.01)
integrator.methods.append(nvt)

# Assign integrator to simulation
sim.operations.integrator = integrator

# Set initial velocities
sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=kT)


# Computes
thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
sim.operations.computes.append(thermodynamic_properties)
sim.run(0)

class DiameterComposition():

    def __init__(self, sim):
        self.sim = sim
        
    @property
    def composition(self):
    	snapshot = self.sim.state.get_snapshot()
    	return numpy.mean(snapshot.particles.diameter)


#Logger
logger = hoomd.logging.Logger()
logger.add(thermodynamic_properties)

diametercomp = DiameterComposition(sim)
logger[('DiameterComposition', 'composition')] = (diametercomp, 'composition', 'scalar')

print("Before setting updater\n")

molecules = hoomd.plugin_insertion.molecule.Molecules(hoomd.filter.All())
##forcecompute = hoomd._hoomd.ForceCompute(hoomd._hoomd.SystemDefinition)
insert_chain = hoomd.plugin_insertion.update.SetDiametersToOne(hoomd.trigger.Periodic(Ninsert,5), molecules, lj)

insert_chain.beta = beta
insert_chain.mu  =  1.62*(46.8885501668-29.204276179999997)/(1.0*Nbeads*beta)#1.201*46.8885501668/(1.0*Nbeads)

for i in range(2):
	for j in range(i,2):
	
		PairTypes = [snapshot.particles.types[i], snapshot.particles.types[j]]
		A = [np.where(Aminoacids == PairTypes[0])[0][0], np.where(Aminoacids == PairTypes[1])[0][0]]
		Params_i = Aminoacids[A[i]]
		Params_j = Aminoacids[A[j]]
		print(i,j,PairTypes,A,Params_i,Params_j)
	
		insert_chain.params[(snapshot.particles.types[i], snapshot.particles.types[j])] = dict(epsilon=eps, sigma=0.5*(float(Params_i[3]) + float(Params_j[3])), kappa=kappa, lamb=0.5*(float(Params_i[4]) + float(Params_j[4])))
		insert_chain.r_cut[(snapshot.particles.types[i], snapshot.particles.types[j])] = rcut# 2*(rcut_max + r_buff) < L


#insert_chain.params[('R', 'R')] = dict(epsilon=eps, sigma=sigma, kappa=kappa, lamb=l)
#insert_chain.r_cut[('R', 'R')] = rcut# 2*(rcut_max + r_buff) < L

#insert_chain.params[('R', 'V')] = dict(epsilon=eps, sigma=sigma, kappa=kappa, lamb=l)
#insert_chain.r_cut[('R', 'V')] = rcut# 2*(rcut_max + r_buff) < L

#insert_chain.params[('V', 'V')] = dict(epsilon=eps, sigma=sigma, kappa=kappa, lamb=l)
#insert_chain.r_cut[('V', 'V')] = rcut# 2*(rcut_max + r_buff) < L

insert_chain.bond_params = [10000,bond_length]

#print(insert_chain.params[('R', 'R')])
#print(insert_chain.params[('R', 'V')])
#print(insert_chain.params[('V', 'V')])

teleport_chain = hoomd.plugin_insertion.update.ChainTeleporting(hoomd.trigger.Periodic(Ninsert,0), molecules)
teleport_chain.beta = beta

sim.operations.computes.append(molecules)
sim.operations += insert_chain
sim.operations += teleport_chain

print("After setting updater\n")

# Write gsd file
gsd_writer = hoomd.write.GSD(trigger=hoomd.trigger.Periodic(Nout),filename='log_'+ SimName + '.gsd',filter=hoomd.filter.All(),mode='wb',dynamic=['attribute','property'])
sim.operations.writers.append(gsd_writer)
gsd_writer.log = logger

#print(integrator.forces[2].params.param_dict)
print("Before running simulation\n")
print("#Type timestep molecule Wratio DrawNumber Accept\n")
sim.run(Nsteps)
print("After running simulation\n")


# Write outputs
thermo = open("thermodynamics_"+SimName+".dat","w")
thermo.write("#t diameter pe ke p T\n")

traj = gsd.hoomd.open('log_'+ SimName + '.gsd', 'rb')
print(len(traj))


traj = gsd.hoomd.open('log_'+ SimName + '.gsd', 'rb')
print(len(traj))


for frame in traj:
	thermo.write("%f %f %f %f %f %f\n"%(frame.configuration.step,frame.log['DiameterComposition/composition'],frame.log['md/compute/ThermodynamicQuantities/potential_energy'][0],frame.log['md/compute/ThermodynamicQuantities/kinetic_energy'][0],frame.log['md/compute/ThermodynamicQuantities/pressure'][0],frame.log['md/compute/ThermodynamicQuantities/kinetic_temperature'][0]))


thermo.close()

print("END\n")
"""

print("Done\n")
