
import numpy as np
import hoomd
import hoomd.md
import copy

from hoomd import ashbaugh_plugin

def determineDihedrals( Sequences,IDs, dihedral_dict, MixingRule="1-1001-1"):
    IDs = IDs.astype(int)+1
    atomid=-1
    id_list = []
    for (SeqId, Seq) in enumerate(Sequences):
        seq_keys = []
        for (ResId,Res) in enumerate(Seq[:]):
            atomid +=1
            if (ResId>(len(Seq)-4) ):
                continue
            else:
                if MixingRule=="1-1001-1":
                    Res_min = IDs[atomid-1] if (ResId-1)>=0 else 0
                    Res_max = IDs[atomid+4] if (ResId)<(len(Seq)-4) else -1
                    key = (np.sort([Res_min,IDs[atomid],IDs[atomid+3], Res_max]))
                    key_str = f"{key[0]}-{key[1]}-{key[2]}-{key[3]}"
                elif MixingRule=="1001":
                    key = (np.sort([IDs[atomid],IDs[atomid+3]]))
                    key_str = f"{key[0]}-{key[1]}"
                elif MixingRule=="0110":
                    key = (np.sort([IDs[atomid+1],IDs[atomid+2]]))
                    key_str = f"{key[0]}-{key[1]}"
                id_list.append(dihedral_dict[key_str])
    return id_list
    
def readSequences(fileName):
    file  = open(fileName,'r')
    NBeads=0
    tmp = ""
    Seqs = []
    for line in file.readlines():
        Seqs.append(line.strip())
        NBeads+=len(line.strip())
        tmp+=line.strip()

    NChains = len(Seqs)

    #InputBonds = np.zeros((NBeads-NChains,2), dtype=int)
    InputBonds=[]
    InputAngles = np.zeros((NBeads-2*NChains,3), dtype=int)
    InputDihedrals = np.zeros((NBeads-3*NChains,4), dtype=int)

    cb =0
    ca =0
    cd =0
    cnt = 0
    for Seq in Seqs:
        for (i,aa) in enumerate(Seq[:-1]):
            #InputBonds[cb+i,:] = [cnt+i,cnt+i+1]
            InputBonds.append([cnt+i, cnt+i+1])

        for (i,aa) in enumerate(Seq[:-2]):
            InputAngles[ca+i,:] = [cnt+i,cnt+i+1,cnt+i+2]

        for (i,aa) in enumerate(Seq[:-3]):
            InputDihedrals[cd+i,:] = [cnt+i,cnt+i+1,cnt+i+2, cnt+i+3]
        cb += len(Seq)-1
        ca += len(Seq)-2
        cd += len(Seq)-3
        cnt += len(Seq)

    return Seqs, NBeads, NChains, InputBonds, InputAngles, InputDihedrals

def readParticleData(fileName, Nparticles, Sequences):
    InputPositions = np.zeros((Nparticles,3), dtype=np.float32)
    InputImage = np.zeros((Nparticles,3), dtype=np.int32)

    i, InputTypes, x,y,z, InputCharges,InputMasses, Diameter, ix, iy,iz = np.genfromtxt(fileName, delimiter=",", comments="#", unpack=True)
    InputPositions[:,0] = x/10.0
    InputPositions[:,1] = y/10.0
    InputPositions[:,2] = z/10.0
    InputTypes -= 1 #convert from array start 1 to 0

    InputImage[:,0]=ix
    InputImage[:,1]=iy
    InputImage[:,2]=iz

    TypeToID=dict()
    cnt=0
    for Seq in Sequences:
        for aa in Seq:
            TypeToID[aa] = InputTypes[cnt]
        cnt+=1
    IDToType= {v: k for k, v in TypeToID.items()}
    Types =np.array([IDToType[i] for i in sorted(IDToType.keys())])

    return InputPositions, InputTypes, InputCharges, InputMasses, Types, Diameter, InputImage

def readDictionaries(filename):
    ID,resname, q, mass, sigma, lamda_val = np.genfromtxt(filename, delimiter=", ", comments="//", unpack=True, dtype=str)#int, converters={0: lambda x: int(x), 1: lambda x: str(x).strip(), 2: lambda x: float(x), 3: lambda x: float(x), 4: lambda x: float(x), 5: lambda x: float(x) })
    ID=ID.astype(int)-1
    resname = np.array(resname)
    IDToResName = dict(zip(ID, resname.astype(str)))
    IDToCharge  = dict(zip(ID, q.astype(float)))
    IDToMass    = dict(zip(ID, mass.astype(float)))
    IDToSigma   = dict(zip(ID, sigma.astype(float)))
    IDToLambda  = dict(zip(ID, lamda_val.astype(float)))
    return ID, IDToResName, IDToCharge, IDToMass, IDToSigma, IDToLambda

def readDihedrals(fileName, Sequences, IDs):
    file  = open(fileName,'r')
    line = file.readline()
    dihedrals = np.array(line.strip("# [").strip("] \n").split(", "), dtype=float)
    dihedral_dict= dict()
    for line in file.readlines():
        data = np.array(line.split(), dtype=int)
        dihedral_dict[f"{data[0]}-{data[1]}-{data[2]}-{data[3]}"] = int(data[4])-1
    
    IDs = np.array(IDs, dtype=int)

    dihedral_list= np.array(sorted(dihedral_dict, key=dihedral_dict.get)) # dihTypesLookup
    dihedral_IDs = np.array(sorted(dihedral_dict.values()), dtype=int)

    id_list = determineDihedrals(Sequences,IDs, dihedral_dict)

    return dihedrals, dihedral_dict, dihedral_list, dihedral_IDs, id_list

def parseKeywords(keyword, value):
    if keyword in ["Simname", "Trajectory","Alt_GSD_Start", "Device"]:
        return value
    if keyword in ["Seed", "NSteps", "NOut"]:
        return int(value)
    if keyword in ["dt", "Lx", "Ly", "Lz", "Temp","epsilon_r", "kappa"]:
        return float(value)
    if keyword in ["Minimise", "UseAngles", "UseCharge", "Create_Start_Config"]:
        return value=="True" or value=="true"
    return value

def readParam(filename):
    ParamDict = dict()
    ParamDict["Create_Start_Config"] = True
    file = open(filename, "r")
    for line in file.readlines():
        key, value = line.split()
        key = key.strip(":")
        ParamDict[key] = parseKeywords(key, value)

    if ParamDict["Alt_GSD_Start"]=="-" and ParamDict["Minimise"]:
        ParamDict["Use_Minimised_GSD"]=True
    else:
        ParamDict["Use_Minimised_GSD"]=False

    return ParamDict
    
def getAnglePotential():
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
    return anglePotential

def getDihedralPotential(dihedral_list, dihedral_IDs, dihedral_eps):
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

    def Dihedral_Table1(x,k,phi0):
        return np.exp( - k * (x - phi0)**2 )
    def Dihedral_Table2(x,k,phi0,e1):
        return np.exp( - k * (x - phi0)**2 + e1 ) + np.exp( - k * (x - phi0 - 2*np.pi)**2 + e1 )
    def Dihedral_Table3(x,ka,kb,phi0a,phi0b,e0,e2):
        return np.exp( - ka * (x - phi0a)**4 + e0 ) + np.exp( - ka * (x - phi0a + 2*np.pi)**4 + e0 ) + np.exp( - kb * (x - phi0b)**4 + e2 ) + np.exp( - kb * (x - phi0b - 2*np.pi)**4 + e2 )
    def GetOPLS(x,k1,k2,k3,k4):
        return 0.5 * ( k1*(1+np.cos(x)) + k2*(1-np.cos(2*x)) + k3*(1+np.cos(3*x)) + k4*(1-np.cos(4*x)) )
    def GetTauOPLS(x,k1,k2,k3,k4):
        return -0.5 * ( k1*(-np.sin(x)) + k2*(+2.0*np.sin(2*x)) + k3*(-3.0*np.sin(3*x)) + k4*(+4.0*np.sin(4*x)) )

    y1 = Dihedral_Table1(dih_x,dih_ka1,dih_phia1)
    y2 = Dihedral_Table2(dih_x,dih_kb1,dih_phib1,dih_e1)
    y3 = Dihedral_Table3(dih_x,dih_ka2,dih_kb2,dih_phia2,dih_phib2,dih_e0,dih_e2)
    
    y_dihedralPotential = (y1*np.exp(-dih_eps) + y2*np.exp(dih_eps) + y3)
    y_dihedralPotential = -np.log(y_dihedralPotential)
    
    tau1 = np.gradient(y1,dih_dx)
    tau2 = np.gradient(y2,dih_dx)
    tau3 = np.gradient(y3,dih_dx)
 

    dihedralPotential = hoomd.plugin_HPS_SS.dihedral.DihedralSS(width = N_dih_points, u1 = y1, u2 = y2, u3 = y3, tau1 = tau1, tau2 = tau2, tau3 = tau3, prefactor=4.18) ### prefactor converts to kJ to kcal

    for key, id in zip(dihedral_list, dihedral_IDs):
        dihedralPotential.params[key] = dict(eps = dihedral_eps[id])

    return dihedralPotential

def minimiseSystem(sim,cell, forcesWithoutLJ, IDToResName, IDToLambda, IDToSigma, kT, Params):
    forces = copy.deepcopy(forcesWithoutLJ)
    # # nonbonded: ashbaugh-hatch potential
    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=kT)
    fire = hoomd.md.minimize.FIRE(dt=0.001, force_tol=1e-3,angmom_tol=1e-3,energy_tol=1e-3)
    fire.methods.append(hoomd.md.methods.ConstantVolume(hoomd.filter.All()))
    sim.operations.integrator = fire

    if True:
        ash =hoomd.md.pair.LJ(nlist=cell, default_r_on=0.0, mode="shift", tail_correction=False)
        for i in IDToResName.keys():
            res_i = IDToResName[i]
            for j in IDToResName.keys():
                res_j = IDToResName[j]
                ash.params[(res_i, res_j)] = {"epsilon":0.8368, "sigma":0.1}
                ash.r_cut[(res_i,res_j)] = 2.0

        sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=kT)

        
        print("first minimisation")
        #dt = 0.005

        f_tmp = copy.deepcopy(forcesWithoutLJ)
        f_tmp.append(ash)
        sim.operations.integrator.forces= f_tmp
        cnt=0
        while not(fire.converged)and cnt<100:
            if (cnt%10)==0: print(cnt)
            sim.run(100)
            cnt+=1
        print("done")




        ash =hoomd.md.pair.LJ(nlist=cell, default_r_on=0.0, mode="shift", tail_correction=False)
        for i in IDToResName.keys():
            res_i = IDToResName[i]
            for j in IDToResName.keys():
                res_j = IDToResName[j]
                ash.params[(res_i, res_j)] = {"epsilon":0.8368, "sigma":0.2}
                ash.r_cut[(res_i,res_j)] = 2.0

        sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=kT)

        
        print("first minimisation")
        #dt = 0.005
        fire = hoomd.md.minimize.FIRE(dt=0.001, force_tol=1e-2,angmom_tol=1e-2,energy_tol=1e-4)
        fire.methods.append(hoomd.md.methods.ConstantVolume(hoomd.filter.All()))
        sim.operations.integrator = fire
        f_tmp = copy.deepcopy(forcesWithoutLJ)
        f_tmp.append(ash)
        sim.operations.integrator.forces= f_tmp
        cnt=0
        while not(fire.converged) and cnt<100:
            if (cnt%10)==0: print(cnt)
            sim.run(100)
            cnt+=1
        print("done")


    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=kT/100.0)

    ash =hoomd.md.pair.LJ(nlist=cell, default_r_on=0.0, mode="shift", tail_correction=False)
    for i in IDToResName.keys():
        res_i = IDToResName[i]
        for j in IDToResName.keys():
            res_j = IDToResName[j]
            ash.params[(res_i, res_j)] = {"epsilon":0.8, "sigma":0.40}
            ash.r_cut[(res_i,res_j)] = 2.0


    print("second minimisation")
    f_tmp = copy.deepcopy(forcesWithoutLJ)
    f_tmp.append(ash)
    sim.operations.integrator.forces = f_tmp
    cnt=0
    #fire.reset()
    while not(fire.converged) and cnt<100:
        if (cnt%10)==0: print(cnt)
        sim.run(100)
        cnt+=1
    print("done")

    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=kT/1.0)

    ash = ashbaugh_plugin.pair.AshbaughPair(nlist=cell, default_r_cut=2.0, default_r_on=0.0)
    for i in IDToResName.keys():
        res_i = IDToResName[i]
        for j in IDToResName.keys():
            res_j = IDToResName[j]
            ash.params[(res_i, res_j)] = {"epsilon":0.8368, "sigma":round((IDToSigma[i]+IDToSigma[j])/20.0,5), "lam":(IDToLambda[i]+IDToLambda[j])/2.0} ### convert to nm as well
            ash.r_cut[(res_i,res_j)] = 2.0


    print("third minimisation")
    f_tmp = copy.deepcopy(forcesWithoutLJ)
    f_tmp.append(ash)
    sim.operations.integrator.forces = f_tmp
    cnt=0
    fire.reset()
    while not(fire.converged) and cnt<100:
        if (cnt%10==0):
            print(cnt)
        sim.run(100)
        cnt+=1
    print("done")

    if False:
        print("md with low time step")

        integrator = hoomd.md.Integrator(dt=Params["dt"]) 
        sim.operations.integrator = integrator
        nvt = hoomd.md.methods.Langevin(filter=hoomd.filter.All(), kT=kT/100) ### ps^-1
        integrator.methods = [nvt]
        sim.operations.integrator.forces = f_tmp

        
        integrator.dt = Params["dt"]/10.0
        sim.run(100000)
        print("md: 1/5 done")

        integrator.dt = Params["dt"]/5.0
        sim.run(200000)
        print("md: 2/5 done")


        integrator.dt = Params["dt"]/2.0
        sim.run(300000)
        print("md: 3/5 done")


        integrator.dt = Params["dt"]/1.5
        sim.run(400000)
        print("md: 4/5 done")

        integrator.dt = Params["dt"]/1.2
        sim.run(500000)
        print("md: fully done")

        integrator.dt = Params["dt"]/1.0
        sim.run(1000000)
        print("md: fully done2")

    return sim









