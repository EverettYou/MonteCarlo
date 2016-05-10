import itertools
import operator
import math
import numpy
'''lattice site object, for construction of lattice'''
class site(object):
    # site constructor given (d,x,f)
    def __init__(self, d, x, f):
        self.d = d # depth
        self.x = x # real space coordinate
        self.f = f # feature space coordinate
        self.UV = [] # UV sites
        self.IR = [] # IR sites
        self.index = None
    def __repr__(self):
        return '<%d,%d,%d>'%(self.d,self.x,self.f)
    def sortkey(self):
        return (self.d,self.f,self.x)
'''lattice system'''
class lattice(object):
    # lattice constructor given UV system size L
    def __init__(self, L, branches=2):
        self.branches = branches # set branches at each RG step
        # initial layer of sites
        layer0 = [site(0,x,0) for x in range(L)]
        self.ftop = dict() # initialize feature top register
        # build RG tree
        self.sites = sorted(self.RGtree(layer0), key=operator.methodcaller('sortkey'))
        # index every site and partition to even/odd sublattices
        self.Nsite = len(self.sites)
        # depth of the last site is the depth of the lattice
        self.depth = self.sites[-1].d
        self.sublattice = {'boundary':[], 'even':[], 'odd':[]}
        for i, a in zip(range(self.Nsite), self.sites):
            a.index = i
            if a.d == 0:
                key = 'boundary'
            elif a.d%2 == 0:
                key = 'even'
            else:
                key = 'odd'
            self.sublattice[key].append(a)
    # recursive build RG tree
    def RGtree(self, this_layer, d=1):
        # if this_layer has exausted
        if len(this_layer) <= 1:
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
* Hamiltonian H = - sum_{ij} J_{ij} chi(g_i, g_j)
    H is sparse and uses a row compressed storage scheme'''
class model(object):
    def __init__(self):
        self.dof = 0 # onsite degrees of freedom
        self.Nsite = 0 # number of sites
        self.config = [] # configuration
        self.chi = numpy.array([[0]]) # chi table, as numpy array
        self.H = [[]] # Hamiltonian in row compressed storage
        self.route = [] # update route
        self.boundary = [] # boundary sites
    # one step of MC update
    def MCstep(self):
        # each MC step goes over the update route once
        for i in self.route:
            h = sum(J*self.chi[self.config[j]] for j, J in self.H[i])
            w = numpy.exp(h)
            self.config[i] = numpy.random.choice(self.dof,p=w/sum(w))
    # calculate energy of the system by config
    def energy(self):
        energy = 0
        for i in range(self.Nsite):
            chii = self.chi[self.config[i]]
            energy -= sum(J*chii[self.config[j]] for j, J in self.H[i])
        return energy
class hypertree_model(model):
    def __init__(self, L=0, dof=0, Js=(1.,1.), branches=2):
        # setup the lattice
        latt = lattice(L, branches=branches)
        # setup the permutation group
        G = group(dof)
        # initialization
        self.dof = dof
        self.Nsite = latt.Nsite
        self.config = [0]*latt.Nsite
        # chi is adjusted by dof
        self.chi = numpy.array(G.chi)-self.dof
        # build Hamiltonian
        self.H = [[] for i in range(latt.Nsite)]
        for i in range(latt.Nsite):
            a = latt.sites[i] # take site i
            bs = a.IR # get its IR sites
            for b, J in zip(bs, Js):
                j = b.index
                self.H[i].append((j,J))
                self.H[j].append((i,J))
        # assign update route (use the sublattice structure)
        self.route = []
        for key in ('odd','even'):
            for a in latt.sublattice[key]:
                self.route.append(a.index)
        # assign boundary sites
        self.boundary = [a.index for a in latt.sublattice['boundary']]