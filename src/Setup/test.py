from PythonFuncs import *

Seqs, NBeads, NChains, InputBonds, InputAngles, InputDihedrals = readSequences("/localscratch/HPS_DATA/HPS-Alpha/HOOMD/RS41_PHOS6/320K/HOOMD_Setup/Sequences.txt")
#print(InputAngles[300:400,:])

InputPositions, InputTypes, InputCharges, InputMasses, Types, Diameter = readParticleData("/localscratch/HPS_DATA/HPS-Alpha/HOOMD/RS41_PHOS6/320K/HOOMD_Setup/Particles.txt", NBeads, Seqs)

dihedrals, dihedral_dict, dihedral_list, dihedral_IDs = readDihedrals("/localscratch/HPS_DATA/HPS-Alpha/HOOMD/RS41_PHOS6/320K/HOOMD_Setup/DihedralMap.txt", Seqs, InputTypes)

Params = readParam("/localscratch/HPS_DATA/HPS-Alpha/HOOMD/RS41_PHOS6/320K/HOOMD_Setup/Params.txt")

IDToResName, IDToCharge, IDToMass, IDToSigma, IDToLambda = readDictionaries("/localscratch/HPS_DATA/HPS-Alpha/HOOMD/RS41_PHOS6/320K/HOOMD_Setup/Dictionaries.txt")

print(IDToResName)
print(IDToCharge)
print(IDToMass)
print(IDToSigma)
print(IDToLambda)