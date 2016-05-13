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
        self.index = {g:i for g, i in zip(self.elements, range(self.dof))}
        self.build_chi() # build chi table
        # after building chi table, self.chi is constructed
    def __repr__(self):
        return repr(self.elements)
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
            (Hamiltonian related) irng, jlst, klst, nlst'''
class Lattice(object):
    # initialize an empty lattice
    def __init__(self):
        self.empty()
        self.build_Hamiltonian()
    def empty(self):
        self.nsite = 0 # number of lattice sites
        self.na, self.nb, self.nc = 0,0,0 # sublattice sites partition
        self.sites = [] # lattice sites (container)
        self.coordmap = {} # coordinate map: coordinate tuple => site object
    def __repr__(self):
        return repr(self.sites)
    # for site in Lattice will iterate over all sites
    def __iter__(self):
        return iter(self.sites)
    # Lattice(coord) will return the site specified by the coordinate tuple
    def __getitem__(self, coord):
        return self.coordmap[coord]
    # build_Hamiltonian returns the tuple (irng, jlst, klst, nlst) as numpy array
    # which will be passed to FORTRAN MC extension core in the Model class
    # From these data, Hamiltonian is constructed as (in a row compressed fassion)
    #     slice(i) := irng[i]:irng[i+1]-1
    #     H = - sum_i [chi(g(i),g(j)) for j in jlst[slice(i)]] dot klst[slice(i)]
    #     nlst specifies the length of jlst and klst
    # Here we provide a phantom version for null Hamiltonian
    # Note: the inheritant class must redefine the build_Hamiltonian method!
    def build_Hamiltonian(self):
        self.irng = numpy.array([1]*(self.nsite+1))
        self.jlst = numpy.array([],dtype=numpy.int_)
        self.klst = numpy.array([],dtype=numpy.float_)
        self.nlst = 0
# Inheritant class: HypertreeLattice
class HypertreeLattice(Lattice):
    # Site class specific to HypertreeLattice
    class Site(object):
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
    # lattice constructor given UV system size L
    def __init__(self, L=0, Ks=[], branches=2):
        Lattice.empty(self) # setup an empty lattice
        self.L = L # keep system size L
        self.branches = branches # set branches at each RG step
        # UV layer of sites
        layer0 = [self.Site(0,x,0) for x in range(L)]
        self.ftop = {} # initialize feature top register
        # build hypertree
        self.depth = 0
        self.sites = self.hypertree(layer0)
        self.nsite = len(self.sites)
        # on return self.depth has been set to the depth of the tree
        # label sites: A - A sublattice, B - B sublattice, C - boundary
        for site in self.sites:
            if site.d == 0: # the ultimate UV layer is the boundary
                site.label = 'C'
                self.nc += 1
            elif site.d%2 == self.depth%2:
                # A sublattice must contain the ultimate IR layer
                site.label = 'A'
                self.na += 1
            else: # B sublattice contains the rest
                site.label = 'B'
                self.nb += 1
        # sort sites and setup coordinate mapping
        self.sites.sort(key = operator.methodcaller('sortkey'))
        for index, site in zip(range(self.nsite), self.sites):
            site.index = index
            self.coordmap[(site.d,site.f,site.x)] = site
        # build Hamiltonian, after building,
        # self.irng,.jlst,.klst,.nlst will be set
        self.build_Hamiltonian(Ks)
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
            grp_IR = [self.Site(d,x,f+f1) for f1 in range(len(grp_UV))]
            #print('%d,%d,%d'%(d,x,f),grp_UV,grp_IR)
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
    # build Hamiltonian
    def build_Hamiltonian(self, Ks):
        # energy parameters are passed in from Ks (as list or array)
        H = [] # container for Hamiltonian terms
        # each Hamiltonian term is a tuple (i,j,K)
        for site in self.sites:
            i = site.index # take site index
            IRsites = site.IR # get its IR sites
            for IRsite, K in zip(IRsites, Ks):
                j = IRsite.index # take IR site index
                H.append((i,j,K))
                H.append((j,i,K))
        if H:
            ilst, jlst, klst = tuple(zip(*sorted(H)))
        else:
            ilst, jlst, klst = [], [], []
        nlst = len(ilst)
        # now ilst, jlst, klst separately hold (i,j,K) for all terms
        # ilst has been sorted in assending order, convert to irng
        imax = 0
        irng = [0]
        p = -1
        for p, i in zip(range(nlst),ilst):
            while i > imax:
                irng.append(p)
                imax += 1
        irng.extend([p+1]*(self.nsite-imax))
        # return tuple (irng, jlst, klst)
        # Note: +1 is needed to convert irng and jlst to FORTRAN index
        #       missing +1 will lead to bus error in the runtime
        self.irng = numpy.array(irng)+1
        self.jlst = numpy.array(jlst)+1
        self.klst = numpy.array(klst)
        self.nlst = nlst
''' ----- Model System -----
attributes:
    system: dof, nsite, na, nb, nc, nlst, chi, irng, jlst, klst
    state: energy, config, hist
    data: nspin, monitor, energy1, energy2, magnet1, magnet2, spins
methods: run(), measure()
Generic statistical mechanics model:
* onsite dof: g in range(model.dof)
* Hamiltonian H = - sum_{ij} K_{ij} chi(g(i), g(j))'''
import MC # FORTRAN extension: Monte Carlo kernel
class Model(object):
    # specify system parameters
    system_parameters = ('dof','nsite','na','nb','nc','nlst','chi','irng','jlst','klst')
    # Model constructor given system parameter
    # system must be a dict containing keys specified in Model._system_parameters
    def __init__(self, system, state={}, data={}):
        # passing system parameters to MC.core
        self.system = system
        # kernel initialization
        MC.core.init() # private workspace allocated
        # initialize state parameters
        if state == {}:
            self.state = {'config':'FM','energy':'unknown','hist':'unknown'}
        else:
            self.state = state
        self.data = data
    def __getattr__(self, attrname):
        if attrname == 'config':
            if MC.core.config is None:
                return None
            else:
                return MC.core.config-1 # shift back to python index convesion
        elif attrname == 'energy':
            if numpy.isnan(MC.core.energy):
                MC.core.get_energy()
            return numpy.asscalar(MC.core.energy)
        elif attrname == 'hist':
            if MC.core.hist[0] < 0:
                MC.core.get_hist()
            return MC.core.hist
        elif attrname in Model.system_parameters:
            return getattr(MC.core, attrname)
        elif attrname == 'state':
            return {key:getattr(self,key) for key in ('config','energy','hist')}
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
        elif attrname == 'energy':
            if attrval in {'unknown', 'nan'}:
                MC.core.energy = numpy.nan
            else:
                MC.core.energy = attrval
        elif attrname == 'hist':
            if attrval in {'unknown'}:
                MC.core.hist = numpy.empty(self.dof)
                MC.core.hist[0] = -1
            else:
                MC.core.hist = attrval
        elif attrname in Model.system_parameters:
            setattr(MC.core, attrname, attrval)
        elif attrname in {'state', 'system'}:
            for key in attrval:
                setattr(self, key, attrval[key])
        else:
            self.__dict__[attrname] = attrval
    # run MC for steps, under mode 0 or 1
    def run(self, steps=1, mode=0):
        if mode==0: # MC without update physical observibles
            MC.core.run(steps, 0)
        elif mode==1: # MC with physical observibles updated
            MC.core.run(steps, 1)
        else: # unrecognised mode
            raise ValueError('The mode %s of Model.run should be 0 or 1.'%repr(mode))
        return self
    # take measurement for steps, monitoring specified spins
    def measure(self, steps=1, monitor=None):
        if monitor is None:
            MC.physics.nspin = 0
            MC.physics.monitor = numpy.array([],dtype=numpy.int_)
        else:
            if all(0<= i < self.lattice.nsite for i in monitor):
                MC.physics.nspin = len(monitor)
                MC.physics.monitor = numpy.array(monitor)+1
            else:
                raise ValueError('The monitor %s of Model.measure out of lattice range.'%repr(monitor))
        MC.physics.measure(steps)
        self.data = {
            'energy1': numpy.asscalar(MC.physics.energy1), # energy 1st moment
            'energy2': numpy.asscalar(MC.physics.energy2), # energy 2nd moment
            'magnet1': MC.physics.magnet1, # bulk magnetization 1st moment
            'magnet2': MC.physics.magnet2, # bulk magnetization 2nd moment 
            'spins':MC.physics.spins,'steps': steps # monitored spin states
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
        para = {key:val for key,val in itertools.chain(lattice.__dict__.items(),group.__dict__.items())}
        system = {name:para[name] for name in Model.system_parameters}
        # call Model method to initialize
        Model.__init__(self, system)