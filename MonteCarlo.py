import itertools
import operator
import math
import numpy
''' ----- Permutation Group System ----- 
attributes: n, elements, order, index, chi'''
class Group(object):
    # construct symmetric group S_n given n
    def __init__(self, n):
        self.n = n
        self.elements = list(itertools.permutations(range(n)))
        self.dof = len(self.elements) # group order
        self.index = {g:i for i,g in enumerate(self.elements)}
        self.build_chi() # build chi table
        # after building chi table, self.chi is constructed
    def __repr__(self):
        return '<%s.Group S_%d (%d)>'%(__name__,self.n,self.dof)
    def __iter__(self):
        return iter(self.elements)
    # group multiplication
    def multiply(self, g1, g2):
        return tuple(g1[i] for i in g2)
    # group inversion
    def inverse(self, g):
        return tuple(sorted(range(self.n), key=lambda k: g[k]))
    # permutation character (number of cycles)
    # the range: character =  1,2,...,n
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
        # prepare the characters adjusted by n, such that chi <= 0
        # the reason to adjust by n is to avoid overflow of partition weight at low temperature
        chis = [self.character(g)-self.n for g in self.elements]
        # construct chi table
        self.chi = [[chis[self.index[self.multiply(g1,g2)]] for g1 in g1s] for g2 in g2s]
''' ----- Lattice System ----- 
attributes: nsite, na, nb, nc, sites, coordmap,
            (Action related) irng, jlst, klst, nlst'''
class Lattice(object):
    # initialize an empty lattice
    def __init__(self):
        # get lattice sites (container)
        self.sites = self.get_sites()
        # --- counting sites ---
        self.nsite = len(self.sites) # put the number of sites
        # make sublattice partition, and count sublattice sites
        self.na, self.nb, self.nc = 0,0,0
        for site in self.sites:
            label = site.get_label(self)
            # label sites: A - A sublattice, B - B sublattice, C - boundary
            if label == 'A':
                self.na += 1
            elif label == 'B':
                self.nb += 1
            elif label == 'C':
                self.nc += 1
        assert(self.na+self.nb+self.nc == self.nsite),'nsite=%d is not equal to the toal of na,nb,nc=%d,%d,%d.'%(self.nsite,self.na,self.nb,self.nc)
        # --- sorting and indexing ---
        # sort by site.sortkey (Site object must provide sortkey method)
        self.sites.sort(key = operator.methodcaller('sortkey'))
        # index the sites and setup coordinate map
        self.coordmap = {} # coordinate map: coordinate tuple => site object
        for index, site in enumerate(self.sites):
            site.index = index
            self.coordmap[site.coordinate] = site
        # --- build Action ---
        # get Action terms
        S = self.Action()
        if S: # if not null
            ilst, jlst, klst = tuple(zip(*sorted(S)))
        else: # for null Action
            ilst, jlst, klst = [], [], []
        # length of the list, by zip construction all lists are of equal len
        nlst = len(ilst)
        # now ilst, jlst, klst separately hold the sequence of (i,j,K) for all terms
        # ilst has been sorted in assending order, convert it to irng
        imax = 0
        irng = [0]
        p = -1
        for p, i in enumerate(ilst):
            while i > imax:
                irng.append(p)
                imax += 1
        # pad irng by p+1 until len(irng) = nsite+1
        irng.extend([p+1]*(self.nsite-imax))
        # set irng,jlst,klst,nlst as numpy arrays
        # which will be passed to FORTRAN MC module in the Model class
        self.irng = numpy.array(irng)+1
        self.jlst = numpy.array(jlst)+1
        self.klst = numpy.array(klst)
        self.nlst = nlst
        # Note: +1 is needed to convert irng and jlst to FORTRAN index
        #       missing +1 will lead to bus error in the runtime
    def __repr__(self):
        return '<%s.Lattice %d:[%d|%d|%d]>'%(__name__,self.nsite,self.na,self.nb,self.nc)
    # for site in Lattice will iterate over all sites
    def __iter__(self):
        return iter(self.sites)
    # Lattice(coord) will return the site specified by the coordinate tuple
    def __getitem__(self, coord):
        return self.coordmap[coord]
    # Note: the inheritant class must redefine the get_sites and build_Action method!
    # build an empty lattice
    def get_sites(self):
        return []
    # return an empty Action
    def Action(self):
        return []
# Inheritant class: HypertreeLattice
class HypertreeLattice(Lattice):
    # Site class specific to HypertreeLattice
    class Site(object):
        # site constructor given coord=(d,f,x)
        def __init__(self, d,f,x):
            self.d = d # depth
            self.f = f # feature space coordinate
            self.x = x # real space coordinate
            self.coordinate = (d,f,x) # coordinate tuple
            self.UV = [] # UV sites
            self.IR = [] # IR sites
            self.label = None
            self.index = None
        def __repr__(self):
            if self.index is None:
                return '<%d,%d,%d>'%(self.d,self.f,self.x)
            else:
                return '%d:<%d,%d,%d>'%(self.index,self.d,self.f,self.x)
        def get_label(self, lattice):
            if self.d == 0: # the ultimate UV layer -> boundary
                self.label = 'C'
            elif self.d%2 == lattice.depth%2: # A sublattice must contain the ultimate IR layer
                self.label = 'A'
            else: # B sublattice contains the rest
                self.label = 'B'
            return self.label
        def sortkey(self):
            return (self.label,-self.d,self.f,self.x)
    # lattice constructor
    def __init__(self, L=0, Ks=[], branches=2):
        # broadcast parameters
        self.L = L
        self.Ks = Ks
        self.branches = branches
        # call Lattice to initialize
        super().__init__()
    # build hypertree lattice
    def get_sites(self):
        # build hypertree
        self.ftop = {} # initialize feature top register
        self.depth = 0
        # on return, self.depth will be set to the depth of the tree
        return self.hypertree([self.Site(0,0,x) for x in range(self.L)])
    # build Action
    def Action(self):
        S = [] # container for Action terms
        # each Action term is a tuple (i,j,K)
        for site in self.sites:
            i = site.index # take site index
            IRsites = site.IR # get its IR sites
            for IRsite, K in zip(IRsites, self.Ks):
                j = IRsite.index # take IR site index
                S.append((i,j,K))
                S.append((j,i,K))
        return S
    # recursive build hypertree
    def hypertree(self, this_layer, d=1):
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
            grp_UV = [site for site in group if site is not None]
            # construct IR sites of the same size as UV in the group
            grp_IR = [self.Site(d,f+f1,x) for f1 in range(len(grp_UV))]
            # link UV and IR sites
            for site_UV, site_IR in itertools.product(grp_UV, grp_IR):
                site_UV.IR.append(site_IR)
                site_IR.UV.append(site_UV)
            # add IR sites to the next layers
            next_layers.append(grp_IR)
            x += 1
            df = max(df, len(grp_UV))
        # ftop increase by the number of features developed in this step
        self.ftop[d] += df
        sites = this_layer
        for layer_raw in itertools.zip_longest(*next_layers):
            layer = [site for site in layer_raw if site is not None]
            sites.extend(self.hypertree(layer, d+1))
        return sites
# Inheritant class: SquareLattice
class SquareLattice(Lattice):
    # Site class specific to SquareLattice
    class Site(object):
        # site constructor given (x, y)
        def __init__(self, x, y):
            self.x = x
            self.y = y
            self.coordinate = (x,y)
            self.neighbors = []
            self.label = None
            self.index = None
        def __repr__(self):
            if self.index is None:
                return '<%d,%d>'%(self.x,self.y)
            else:
                return '%d<%d,%d>'%(self.index,self.x,self.y)
        def get_label(self, lattice):
            if (self.x + self.y)%2: # odd lattice
                self.label = 'B'
            else: # even lattice
                self.label = 'A'
            return self.label
        def sortkey(self):
            return (self.label,self.x,self.y)
    # lattice constructor
    def __init__(self, L=0, Ks=[]):
        # broadcast parameters
        if L%2 != 0:
            raise ValueError('L must be even to ensure that the lattice bipartite.')
        self.L = L
        self.Ks = Ks
        # call Lattice to initialize
        super().__init__()
    # build square lattice
    def get_sites(self):
        return [self.Site(x,y) for x in range(self.L) for y in range(self.L)]
    # build Action
    def Action(self):
        S = [] # container for Action terms
        # each Action term is a tuple (i,j,K)
        for site in self.sites:
            i = site.index # take site index
            # collect neighbors
            site.neighbors.append([self[((site.x+dx)%self.L, (site.y+dy)%self.L)]
                                   for (dx,dy) in ((-1,0),(0,-1),(0,1),(1,0))])
            # add adjencent terms to Action
            for jsites, K in zip(site.neighbors,self.Ks):
                S.extend([(i,jst.index,K) for jst in jsites])
        return S
''' ----- Model System -----
attributes:
    system: dof, nsite, na, nb, nc, nlst, chi, irng, jlst, klst
    state: action, beta, config, hist
    data: nspin, monitor, energy1, energy2, magnet1, magnet2, spins
          nser, eser, mser
methods: run(), measure()
Generic statistical mechanics model:
* onsite dof: g in range(model.dof)
* Action S = - sum_{ij} K_{ij} chi(g(i), g(j))'''
import MC # FORTRAN extension: Monte Carlo kernel
class Model(object):
    # specify system parameters
    system_parameters = ('dof','nsite','na','nb','nc','nlst','chi','irng','jlst','klst')
    # Model constructor given system parameter
    # system must be a dict containing keys specified in Model._system_parameters
    def __init__(self, system, state={}, data={}):
        # passing system parameters to MC.core
        if not all(key in system for key in Model.system_parameters):
            raise ValueError('The dictionary %s passed in as system does not contains all the keys required in %s.'%(repr(system), repr(Model.system_parameters)))
        self.system = system
        # kernel initialization
        MC.core.init() # private workspace allocated
        # initialize state parameters
        if state == {}:
            self.state = {'action':'unknown','beta':'default','config':'FM','hist':'unknown'}
        else:
            self.state = state
        self.data = data
    def __getattr__(self, attrname):
        if attrname == 'config':
            if MC.core.config is None:
                return None
            else:
                return MC.core.config-1 # shift back to python index convesion
        elif attrname == 'action':
            return numpy.asscalar(MC.core.action)
        elif attrname == 'hist':
            return MC.core.hist
        elif attrname == 'beta':
            return numpy.asscalar(MC.core.beta)
        elif attrname in Model.system_parameters:
            return getattr(MC.core, attrname)
        elif attrname == 'state':
            return {key:getattr(self,key) for key in ('action','beta','config','hist')}
        elif attrname == 'system':
            return {key:getattr(self,key) for key in Model.system_parameters}
        else:
            raise AttributeError("'%s' object has no attribute '%s'"%(self.__class__.__name__, attrname))
    def __setattr__(self, attrname, attrval):
        if attrname == 'config':
            if isinstance(attrval, numpy.ndarray):
                MC.core.config = attrval+1 # shift to FORTRAN index convension
            elif isinstance(attrval, (list, tuple)):
                MC.core.config = numpy.array(attrval)+1 # shift to FORTRAN index convension
            elif attrval in {'FM','uniform'}:
                MC.core.config = numpy.full(self.nsite,1,dtype=numpy.int_)
            elif attrval in {'PM','random'}:
                MC.core.config = numpy.random.randint(self.dof, size=self.nsite)+1
            else:
                raise ValueError('Illigal value %s for config.'%repr(attrval))
        elif attrname == 'action':
            if attrval in {'unknown'}:
                MC.core.get_action()
            else:
                MC.core.action = attrval
        elif attrname == 'hist':
            if attrval in {'unknown'}:
                MC.core.hist = numpy.empty(self.dof)
                MC.core.get_hist()
            else:
                MC.core.hist = attrval
        elif attrname == 'beta':
            if attrval in {'default'}:
                MC.core.beta = 1. # to initialize
            else:
                MC.core.set_beta(attrval)
        elif attrname in Model.system_parameters:
            setattr(MC.core, attrname, attrval)
        elif attrname == 'state':
            # order matters! config must be set before beta, action and hist
            # otherwise call to get_action(), get_hist() will bus error
            for key in ('config','beta','action','hist'):
                if key in attrval:
                    setattr(self, key, attrval[key])
        elif attrname == 'system':
            for key in attrval:
                setattr(self, key, attrval[key])
        else:
            self.__dict__[attrname] = attrval
    # run MC for steps, under mode 0 or 1
    def run(self, steps=1):
        MC.core.run(steps)
        return self
    # take measurement for steps, monitoring specified spins
    def measure(self, steps=1, monitor=None):
        if monitor is None:
            MC.data.nspin = 0
            MC.data.monitor = numpy.array([],dtype=numpy.int_)
        else:
            if all(0<= i < self.lattice.nsite for i in monitor):
                MC.data.nspin = len(monitor)
                MC.data.monitor = numpy.array(monitor)+1
            else:
                raise ValueError('The monitor %s of Model.measure out of lattice range.'%repr(monitor))
        MC.data.measure(steps)
        self.data = {
            'energy1': numpy.asscalar(MC.data.energy1), # energy 1st moment
            'energy2': numpy.asscalar(MC.data.energy2), # energy 2nd moment
            'magnet1': MC.data.magnet1, # bulk magnetization 1st moment
            'magnet2': MC.data.magnet2, # bulk magnetization 2nd moment 
            'spins':MC.data.spins,'steps': steps # monitored spin states
        }
        return self.data
# Inheritant class: LatticeModel - construct model from Lattice and Group
# LatticeModel has two more attributes: lattice and group
class LatticeModel(Model):
    def __init__(self, lattice=Lattice(), group=Group(1)):
        # save lattice and group as attribute
        self.lattice = lattice
        self.group = group
        # construct system parameters
        para = {**vars(lattice),**vars(group)}
        system = {name:para[name] for name in Model.system_parameters}
        # call Model method to initialize
        super().__init__(system)