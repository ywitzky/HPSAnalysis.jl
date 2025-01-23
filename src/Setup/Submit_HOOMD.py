import time

import os

import numpy as np

import gsd.hoomd
import hoomd
import hoomd.md
import h5py

import sys
from PythonFuncs import *
from hoomd import ashbaugh_plugin
#from hoomd import plugin_HPS_SS
#import ashbaugh_plugin as aplugin

PWD = os.getcwd()

def run(FolderPath):#, GPUNUM):
    ### Read Input Data
    ### All inputs are in lammps units, have to convert to 
    Params = readParam(f"{FolderPath}/HOOMD_Setup/Params.txt")

    Seqs, NBeads, NChains, InputBonds, InputAngles, InputDihedrals = readSequences(f"{FolderPath}/HOOMD_Setup/Sequences.txt")

    InputPositions, InputTypes, InputCharges, InputMasses, _, Diameter, InputImage = readParticleData(f"{FolderPath}/HOOMD_Setup/Particles.txt", NBeads, Seqs)

    if Params["UseAngles"]:
        dihedral_eps, dihedral_dict, dihedral_list, dihedral_IDs, dihedral_AllIDs = readDihedrals(f"{FolderPath}/HOOMD_Setup/DihedralMap.txt", Seqs, InputTypes)


    IDS, IDToResName, IDToCharge, IDToMass, IDToSigma, IDToLambda = readDictionaries(f"{FolderPath}/HOOMD_Setup/Dictionaries.txt")


    ### constants 
    bondLength = 0.38


    kb = 0.00831446262
    kT = kb*Params["Temp"]

    ### Prepare HOOMD 
    
    tmp = []
    for i in IDS: 
        tmp.append(str(IDToResName[i]))
    Types = np.array(tmp)


    if Params["Device"]=="GPU":
        device = hoomd.device.GPU()
    else:
        device = hoomd.device.CPU()

    sim = hoomd.Simulation(device=device, seed=Params["Seed"])

    if Params["Create_Start_Config"]:
        snapshot = gsd.hoomd.Frame()  
        #if Params["Use_Minimised_GSD"]:
            ### Specify particles
        snapshot.particles.N = NBeads
        snapshot.particles.position = InputPositions.astype(np.float32) ### convert to nm
        snapshot.particles.types = tuple(Types)
        snapshot.particles.typeid = InputTypes.astype(np.int32)
        snapshot.particles.image = InputImage

        snapshot.configuration.box = [Params["Lx"], Params["Ly"], Params["Lz"], 0, 0, 0] #4:6 are tilt

        snapshot.particles.mass = InputMasses.astype(np.float32)
        snapshot.particles.charge = InputCharges.astype(np.float32)

        # Connect particles with bonds.
        snapshot.bonds.N = NBeads-NChains
        snapshot.bonds.types = ['O-O']
        snapshot.bonds.typeid = np.zeros(snapshot.bonds.N, dtype=np.uint32)
        snapshot.bonds.group = InputBonds#.astype(np.uint32)

        ## Create Angles
        snapshot.angles.N = NBeads-2*NChains
        snapshot.angles.types = ['O-O-O']
        snapshot.angles.typeid = np.zeros( snapshot.angles.N, dtype=int)
        snapshot.angles.group = InputAngles

        # Create Angles
        snapshot.dihedrals.N =  NBeads-3*NChains 
        snapshot.dihedrals.types = list(dihedral_list)
        snapshot.dihedrals.typeid =  dihedral_AllIDs
        snapshot.dihedrals.group = InputDihedrals

        with gsd.hoomd.open(name=FolderPath+Params["Simname"] + "_StartConfiguration.gsd", mode='w') as f:
            f.append(snapshot)

        sim.create_state_from_gsd(filename=FolderPath+Params['Simname']+"_StartConfiguration.gsd")
    else:
        print(f"Creat start from: {FolderPath+Params['Simname']}gsd")
        sim.create_state_from_gsd(filename=FolderPath+Params["Simname"]+".gsd")


    forces = []

    ### Harmonic bonds
    harmonic = hoomd.md.bond.Harmonic()
    harmonic.params['O-O'] = dict(k=8033, r0=bondLength) ###calvados2: k=8033kJ/mol/nm^2 k=1000kJ/nm^2 = 10KJ/AA^2
    forces.append(harmonic)

    if Params["UseAngles"]:
        cell2 = hoomd.md.nlist.Tree(buffer=0.4,default_r_cut=2.0, exclusions=['bond', 'angle', 'dihedral', '1-3', '1-4'])
        cell  = hoomd.md.nlist.Tree(buffer=0.4,default_r_cut=0.0, exclusions=['bond', 'angle', 'dihedral', '1-3', '1-4'] )#, mesh=hoomd.mesh.Mesh(), default_r_cut=3.5)
    else:
        cell2 = hoomd.md.nlist.Tree(buffer=0.4,default_r_cut=2.0, exclusions=['bond'])
        cell  = hoomd.md.nlist.Tree(buffer=0.4,default_r_cut=0.0, exclusions=['bond'] )#, mesh=hoomd.mesh.Mesh(), default_r_cut=3.5)
        #cell  = hoomd.md.nlist.Cell(buffer=0.4,default_r_cut=0.0, exclusions=['bond'] )#, mesh=hoomd.mesh.Mesh(), default_r_cut=3.5)


    if Params["UseCharge"]:
        # # electrostatics forces
        yukawa = hoomd.md.pair.Yukawa(nlist=cell)
        for i in IDToResName.keys():
            res_i = IDToResName[i]
            for j in IDToResName.keys():
                res_j = IDToResName[j]
                yukawa.params[(res_i,res_j)] = dict(epsilon=IDToCharge[i]*IDToCharge[j]*Params["epsilon_r"], kappa=Params["kappa"]) ### (KJ , 1/nm)d

                if ((IDToCharge[i]==0) or (IDToCharge[j]==0)):
                    yukawa.r_cut[(res_i,res_j)] = 0.0
                else:
                    yukawa.r_cut[(res_i,res_j)] = 3.5
        forces.append(yukawa)
    
    # # nonbonded: ashbaugh-hatch potential
    ash = ashbaugh_plugin.pair.AshbaughPair(nlist=cell2, default_r_cut=2.0, default_r_on=0.0)
    for i in IDToResName.keys():
        res_i = IDToResName[i]
        for j in IDToResName.keys():
            res_j = IDToResName[j]
            ash.params[(res_i, res_j)] = {"epsilon":0.8368, "sigma":round((IDToSigma[i]+IDToSigma[j])/20.0,5), "lam":(IDToLambda[i]+IDToLambda[j])/2.0}#, "lam":(IDToLambda[i]+IDToLambda[j])/2.0} ### convert to nm as well
            ash.r_cut[(res_i,res_j)] = 2.0
    forces.append(ash)


    if Params["UseAngles"]:
        # ANGLE POTENTIAL 
        anglePotential = getAnglePotential()
        forces.append(anglePotential)

        ### DIHEDRAL
        dihedralPotential = getDihedralPotential(dihedral_list, dihedral_IDs, dihedral_eps)
        forces.append(dihedralPotential)


    if Params["Alt_GSD_Start"]!="-":  ### not needed?
        sim = hoomd.Simulation(device=device, seed=Params["Seed"])
        sim.create_state_from_gsd(filename=FolderPath+Params["Alt_GSD_Start"])


    ### logging quantitites
    thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    sim.operations.computes.append(thermodynamic_properties)

    logger = hoomd.logging.Logger(categories=['scalar', 'string'])
    #logger.add(thermodynamic_properties)

    gsd_writer = hoomd.write.GSD(trigger=hoomd.trigger.Periodic(Params["NOut"]),filename=FolderPath+Params["Trajectory"] ,filter=hoomd.filter.All(),mode='wb',dynamic=['particles/position', 'particles/image'])
    sim.operations.writers.append(gsd_writer)
    gsd_writer.log = logger

    #logger = hoomd.logging.Logger(categories=['scalar', 'string']) #'sequence'
    logger.add(sim, quantities=['timestep', 'walltime', 'tps'])
    table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(100000), logger=logger)
    sim.operations.writers.append(table)
        
    ### track presure
    logger_pres = hoomd.logging.Logger(categories=['scalar', 'sequence'])
    logger_pres.add(sim, quantities=["timestep"])
    logger_pres.add(thermodynamic_properties, quantities=['kinetic_temperature','kinetic_energy', 'potential_energy','pressure', 'pressure_tensor'])
    hdf5_writer = hoomd.write.HDF5Log(
        trigger=hoomd.trigger.Periodic(100), filename=FolderPath+'pressure.h5', mode='w', logger=logger_pres
    )
    sim.operations.writers.append(hdf5_writer)


    ### apply langevin
    integrator = hoomd.md.Integrator(dt=Params["dt"]) 
    nvt = hoomd.md.methods.Langevin(filter=hoomd.filter.All(), kT=kT) ### ps^-1
    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=kT/100.0)
    integrator.methods = [nvt]
    integrator.dt = Params["dt"]


    for i in IDToResName.keys():
        name = IDToResName[i]
        nvt.gamma[name] = IDToMass[i]*10.0**-5
        nvt.gamma_r[name] = (0.0, 0.0, 0.0)
    sim.operations.integrator=integrator

	#

    sim.operations.integrator.forces=forces

    print("Before simulation\n")

    ### optimise cell list buffer
    now = sim.timestep
    if now is None:
        now = 0
    hoomd.md.tune.NeighborListBuffer(trigger=hoomd.trigger.Before(now+20000), nlist=cell , maximum_buffer=1.0, solver=hoomd.tune.GradientDescent())
    hoomd.md.tune.NeighborListBuffer(trigger=hoomd.trigger.Before(now+20000), nlist=cell2, maximum_buffer=1.0, solver=hoomd.tune.GradientDescent())

    ### start simulation
    sim.run(Params["NSteps"])
    
    print(f"TPS: {sim.tps:0.5g}")




if __name__ == '__main__':
    if len(sys.argv)<2:
        print("Need folder of input parameters as second argument.")
    else:
        print(sys.argv[1])#, sys.argv[2])
        run(sys.argv[1])#, int(sys.argv[2]))

 
