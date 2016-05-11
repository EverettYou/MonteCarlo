import itertools
import operator
import math
import numpy
'''lattice site object, for construction of lattice'''
class site(object):
    # site constructor given (d,x,f)
    def __init__(self, d, x, f):
        self.d = d # depth
        self.f = f # feature space coordinate
        self.x = x # real space coordinate
        self.UV = [] # UV sites
        self.IR = [] # IR sites
        self.label = None
        self.index = None
    def __repr__(self):
        return '<%d,%d,%d>'%(self.d,self.f,self.x)
    def sortkey(self):
        return (self.label,-self.d,self.f,self.x)
'''lattice system'''
class lattice(object):
    # lattice constructor given UV system size L
    def __init__(self, L, branches=2):
        self.L = L # keep system size L
        self.branches = branches # set branches at each RG step
        # initial layer of sites
        layer0 = [site(0,x,0) for x in range(L)]
        self.ftop = dict() # initialize feature top register
        # build RG tree
        self.depth = 0
        self.sites = self.RGtree(layer0)
        self.Nsite = len(self.sites)
        # on return self.depth has been set to the depth of the tree
        # label sites: A - A sublattice, B - B sublattice, C - boundary
        self.NA, self.NB, self.NC = 0,0,0
        for a in self.sites:
            if a.d == 0: # the most UV layer is the boundary
                a.label = 'C'
                self.NC += 1
            elif a.d%2 == self.depth%2:
                # A sublattice must contain the most IR layer
                a.label = 'A'
                self.NA += 1
            else: # B sublattice contains the rest
                a.label = 'B'
                self.NB += 1
        # sort sites and make index
        self.sites.sort(key = operator.methodcaller('sortkey'))
        self.site = dict()
        for i, a in zip(range(self.Nsite), self.sites):
            a.index = i
            self.site[(a.d,a.f,a.x)] = a
    def __iter__(self):
        return iter(self.sites)
    def __getitem__(self, key):
        return self.site[key]
    # recursive build RG tree
    def RGtree(self, this_layer, d=1):
        # if this_layer has exausted
        if len(this_layer) <= 1:
            self.depth = max(self.depth, d-1)
            return this_layer
        # get the top index of features in this depth
        try:
            f = self.ftop[d]
        except: # if failed, set new feature top
            self.ftop[d] = 0
            f = 0
        # prepare next layers
        next_layers = []
        x = 0 # initialized next layer real space
        df = 0 # number of new features
        # partition this_layer by into groups of length = branches
        for group in itertools.zip_longest(*[iter(this_layer)]*self.branches):
            # partition UV sites into UV groups
            grp_UV = [a for a in group if a is not None]
            # construct IR sites of the same size as UV in the group
            grp_IR = [site(d,x,f+f1) for f1 in range(len(grp_UV))]
            #print('%d,%d,%d'%(d,x,f),grp_UV,grp_IR)
            # link UV and IR sites
            for a, b in itertools.product(grp_UV, grp_IR):
                a.IR.append(b)
                b.UV.append(a)
            # add IR sites to the next layers
            next_layers.append(grp_IR)
            x += 1
            df = max(df, len(grp_UV))
        # ftop increase by the number of features developed in this step
        self.ftop[d] += df
        sites = this_layer
        for layer_raw in itertools.zip_longest(*next_layers):
            layer = [a for a in layer_raw if a is not None]
            sites.extend(self.RGtree(layer, d+1))
        return sites
'''symmetric group system'''
class group(object):
    # construct symmetric group S_n given n
    def __init__(self, n):
        self.n = n
        self.elements = list(itertools.permutations(range(n)))
        self.order = len(self.elements) # group order
        self.index = {g:i for g, i in zip(self.elements, range(self.order))}
        self.build_chi() # build chi table
    def __iter__(self):
        return iter(self.elements)
    # group multiplication
    def multiply(self, g1, g2):
        return tuple(g1[i] for i in g2)
    # group inversion
    def inverse(self, g):
        return tuple(sorted(range(self.n), key=lambda k: g[k]))
    # permutation character (number of cycles)
    def character(self, g):
        # setup a set of objects
        objs = set(range(self.n))
        cycles = 0 # initially no cycles
        while objs: # if objects not empty
            obj = objs.pop() # pop out an object
            obj = g[obj] # get its image object under g action
            while obj in objs: # if obj is still in objs
                objs.remove(obj) # remove obj
                obj = g[obj] # and get the next image
            # while loop terminates if a full cycles has been gone through
            cycles += 1 # cycles +1
        return cycles
    # build chi table: chi(g1, g2) = ln Tr g1*g2^(-1)
    def build_chi(self):
        g1s = self.elements
        # prepare the inversions
        g2s = [self.inverse(g) for g in self.elements]
        # prepare the characters
        chis = [self.character(g) for g in self.elements]
        # construct chi table
        self.chi = [[chis[self.index[self.multiply(g1,g2)]] for g1 in g1s] for g2 in g2s]
'''statistical mechanics model system
Generic model:
* onsite dof: g in range(model.dof)
* Hamiltonian H = - sum_{ij} J_{ij} chi(g_i, g_j)'''
import MC # FORTRAN extension: Monte Carlo kernel
class model(object):
    def __init__(self):
        MC.core.init() # kernel initialization (allocation)
    # return current configuration
    # example: model.config(), model.config([0,1,2])
    def config(self, ind=None):
        if ind is None:
            return MC.core.config-1
        else:
            return MC.core.config[ind]-1
    # run MC for steps, under mode 0 or 1
    def run(self, steps=1, mode=0):
        if mode==0: # MC without update physical observibles
            MC.core.run(steps, 0)
        elif mode==1: # MC with physical observibles updated
            MC.core.run(steps, 1)
        else: # unrecognised mode
            raise ValueError('The mode parameter of model.run only takes 0 or 1.')
        return self
    # return energy
    def energy(self):
        if numpy.isnan(MC.core.energy):
            MC.core.get_energy()
        return numpy.asscalar(MC.core.energy)
    # return histogram of bulk spin states (boundary not counted)
    def hist(self):
        if MC.core.hist[0]<0:
            MC.core.get_hist()
        return MC.core.hist
    # take measurement for steps, monitoring specified spins
    def measure(self, steps=1, monitor=None):
        if monitor is not None:
            MC.physics.monitor = numpy.array(monitor)+1
        MC.physics.measure(steps)
        result = {'energy1': numpy.asscalar(MC.physics.energy1), # energy 1st moment
                  'energy2': numpy.asscalar(MC.physics.energy2), # energy 2nd moment
                  'magnet1': MC.physics.magn1, # bulk magnetization 1st moment
                  'magnet2': MC.physics.magn2, # bulk magnetization 2nd moment 
                  'spins':MC.physics.spins,'steps': steps} # monitored spin states
        return result
# inheritances:
class hypertree_model(model):
    def __init__(self, L=0, n=0, Ks=(1.,1.), branches=2):
        # setup the lattice
        latt = lattice(L, branches=branches)
        self.latt = latt
        # setup the permutation group
        grp = group(n)
        self.grp = grp
        dof = grp.order
        # passing data to FORTRAN MC core
        MC.core.dof = dof
        MC.core.nsite = latt.Nsite
        MC.core.na, MC.core.nb, MC.core.nc = latt.NA, latt.NB, latt.NC
        # initialize configuration
        MC.core.config = numpy.array([1]*latt.Nsite)
        # chi is adjusted by dof
        MC.core.chi = numpy.array(grp.chi)-dof
        # build Hamiltonian
        H = []
        for a in latt:
            i = a.index # take site a index
            bs = a.IR # get its IR sites
            for b, K in zip(bs, Ks):
                j = b.index # take site b index
                H.append((i,j,K))
                H.append((j,i,K))
        ilst, jlst, Klst = tuple(zip(*sorted(H)))
        MC.core.irng = self.make_irng(ilst)+1
        MC.core.jlst = numpy.array(jlst)+1
        MC.core.klst = numpy.array(Klst)
        model.__init__(self)
    # convert ilst to irng
    def make_irng(self, ilst):
        imax = 0
        irng = [0]
        for k, i in zip(range(len(ilst)),ilst):
            while i > imax:
                irng.append(k)
                imax += 1
        irng.append(k+1)
        return numpy.array(irng)