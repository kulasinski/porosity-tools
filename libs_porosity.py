import getopt, os, sys, math
import numpy as np
from scipy.optimize import minimize
from matplotlib import pyplot as plt

def print_help():
	print "   This code analyzes pore structure of an MD structure in Gromacs format."
	print "   Usage:"
	print "        -f <input.gro>  # input structure file, in .gro format"
	print "           <grid.npy>   # or previously calculated grid in .npy format"
	print "        -o <out>        # generic name for output files"
	print "        -d <output>     # output directory"
	print "        -r <float>      # resolution of the grid elements in nm, default: 0.01 nm"
	print "        -m <float>      # radius of the probe molecule in nm,"
	print "                          default is 0.14 nm (H2O radius)."
	print "                          If radius is 0 Van der Waals surface is probed."
	print "        -g <gro|npy|no> # save the grid to .gro file, numpy .npy file,"
	print "                          or don't (default)."
	print "        -p <y/n>        # if your structure is fully periodic, default is yes."
	print "        -c <y/n>        # calculate chord and diameter distribution"
	print "        -t <x|y|z>      # calculate tortuosity in given direction"
def process_input(arguments):
	input_structure = 'conf.gro'
	output = 'out'
	odir = 'output'
	spacing = 0.01 # nm
	radius = 0.14 # H2O molecule
	ifsavegrid = 'no'
	ifperiodic = True
	ifdistribution = False
	tdir = ''

	if len(arguments)<2:
		print_help()
		sys.exit(0)
	try:
		opts, args = getopt.getopt(arguments[1:], "f:o:d:r:m:g:p:c:t:no_usable")
		# print opts
	except getopt.GetoptError:
		print "   Inappropriate use."
		print_help()
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print_help()
			sys.exit(0)
		elif opt in "-f":
			input_structure = arg
		elif opt in "-o":
			output = arg
		elif opt in "-d":
			odir = arg
		elif opt in "-r":
			spacing = float(arg)
			if spacing == 0:
				print "   ERROR: Spacing must be greater than 0! Exiting."
				sys.exit(0)
		elif opt in "-m":
			radius = float(arg)
			if radius < 0:
				print "   ERROR: Radius must be larger or equal zero! Exiting."
				sys.exit(0)
		elif opt in "-g":
			if arg in ('gro','.gro'):
				ifsavegrid = 'gro'
			elif arg in ('npy','.npy'):
				ifsavegrid = 'npy'
			if arg in ('No','no'):
				ifsavegrid = 'no'
		elif opt in "-p":
			if arg in ('n','N','No','no'):
				ifperiodic = False
		elif opt in "-c":
			if arg in ('y','Y','yes','Yes'):
				ifdistribution = True
		elif opt in "-t":
			if arg in ('x','X','y','Y','z','Z'):
				tdir = arg.lower()
	return input_structure, output, odir, spacing, radius, ifsavegrid, ifperiodic, ifdistribution, tdir

def read_structure(filename):
	print '   Reading',filename,'...\r',
	fin = open(filename,'r') # in .gro format without velocities
	ll = fin.readlines()
	fin.close()

	# add here later atom type
	N = int(ll[1])
	X = float(ll[-1].split()[0])
	Y = float(ll[-1].split()[1])
	Z = float(ll[-1].split()[2])

	# read coordinates and wrap them up to a untilted shape
	coords = [None]*3*N
	atomtype = [None]*N
	for i in range(0,N):
		lll = ll[i+2]
		atomtype[i] = ll[i+2][9:15].strip()[0]
		# dividing modulo is to wrap up the coords
		coords[3*i]   = float(lll[20:28])%X
		coords[3*i+1] = float(lll[28:36])%Y
		coords[3*i+2] = float(lll[36:44])%Z

	print '   Number of atoms:', N
	print '   System size: %.5f x %.5f x %.5f' % (X, Y, Z)

	return atomtype, coords, X, Y, Z

def distance(x,y,z,X,Y,Z):
	return math.sqrt((x-X)*(x-X)+(y-Y)*(y-Y)+(z-Z)*(z-Z))
def remove_grid_points_pbc(G, atomtype, coords, X,Y,Z, NX,NY,NZ, spacing, radius):
	mindistdict = {'H':0.12,'C':0.17,'N':0.155,'O':0.152,'F':0.147,'P':0.18,'S':0.18,'Cl':0.175}
	probemol_radius = radius # r=0.14 nm for H2O
	overlap = 0.2
	for a in range(0,len(coords)/3):
		sys.stdout.write('   Calculating porosity grid:  %d%%\r' % int(100*float(a)/(len(coords)/3)+1))
		sys.stdout.flush()
		ax = coords[3*a]
		ay = coords[3*a+1]
		az = coords[3*a+2]
		mindist=mindistdict[atomtype[a]]
		dist = (mindist+probemol_radius)*(1-overlap)
		for i in range(int((ax-dist)/spacing)-1,int((ax+dist)/spacing)+2):
			for j in range(int((ay-dist)/spacing)-1,int((ay+dist)/spacing)+2):
				for k in range(int((az-dist)/spacing)-1,int((az+dist)/spacing)+2):
					x = i*spacing
					y = j*spacing
					z = k*spacing
					i_= int(x/spacing+0.001)%NX
					j_= int(y/spacing+0.001)%NY
					k_= int(z/spacing+0.001)%NZ
					if G[i_,j_,k_]>0:
						if distance(ax,ay,az,x,y,z)<=dist:
							G[i_,j_,k_]=0
	return G
def remove_grid_points(G, atomtype, coords, X,Y,Z, NX,NY,NZ, spacing, radius):
	mindistdict = {'H':0.12,'C':0.17,'N':0.155,'O':0.152,'F':0.147,'P':0.18,'S':0.18,'Cl':0.175}
	probemol_radius = radius # r=0.14 nm for H2O
	for a in range(0,len(coords)/3):
		sys.stdout.write('   Calculating porosity grid:  %d%%\r' % int(100*float(a)/(len(coords)/3)+1))
		sys.stdout.flush()
		ax = coords[3*a]
		ay = coords[3*a+1]
		az = coords[3*a+2]
		mindist=mindistdict[atomtype[a]]
		dist = mindist+probemol_radius
		for i in range(max(int((ax-dist)/spacing)-1,0),min(int((ax+dist)/spacing)+2,NX)):
			for j in range(max(int((ay-dist)/spacing)-1,0),min(int((ay+dist)/spacing)+2,NY)):
				for k in range(max(int((az-dist)/spacing)-1,0),min(int((az+dist)/spacing)+2,NZ)):
					if G[i,j,k]>0:
						x = i*spacing
						y = j*spacing
						z = k*spacing
						if distance(ax,ay,az,x,y,z)<=dist:
							G[i,j,k]=0
	return G
def save_grid2gro(filename,G,nx,ny,nz,d):
	n = 0 # current atom number
	N = G.sum() # all atoms
	x = 0 # atom coordinates
	y = 0
	z = 0

	fout = open(filename, 'w')
	fout.write('GRID: atoms are porosity network elements\n')
	fout.write(' %d\n' % N)

	for k in range(0,nz):
		z = k*d
		for j in range(0,ny):
			y = j*d
			for i in range(0,nx):
				x = i*d
				if G[i,j,k] != 0:
					fout.write('    1RESX     C%s   %.3f   %.3f   %.3f\n' % (str(n+1).rjust(5)[-5:], x,y,z))
					n = n + 1
	fout.write('  %f  %f  %f\n' % (nx*d, ny*d, nz*d))
	fout.close()

def get_surf(G,NX,NY,NZ,d):
	grad = np.gradient(G)
	gx,gy,gz = np.abs(grad[0]),np.abs(grad[1]),np.abs(grad[2])
	surf = 0
	for i in range(NX):
		for j in range(NY):
			for k in range(NZ):
				if gx[i,j,k]>0 or gy[i,j,k]>0 or gz[i,j,k]>0:
					surf += 1
	return surf*d*d

def split_chord(ch):
	N = ch.shape[0]
	grad  = ch[1:N]-ch[0:N-1]
	start = np.where(grad==-1)[0]
	end   = np.where(grad== 1)[0]
	# print 'chord',ch
	l_list = []
	for i in range(0,start.shape[0]):
		temp = np.where(end>start[i])[0]
		if temp.shape[0]>0:
			e = temp[0]
			l = end[e]-start[i]
			l_list.append(l)
	# print 'list',l_list
	return l_list
def get_chord_vec_full(N,grid,spacing,ch):
	ch = int(ch/spacing)
	NX,NY,NZ = grid.shape[0],grid.shape[1],grid.shape[2]
	chords = np.zeros([N,max(NX,NY,NZ)])
	for i in range(0,N):
		# toss chord direction
		dir = np.random.randint(3)
		if dir == 0:
			Y = np.random.randint(NY)
			Z = np.random.randint(NZ)
			chords[i,0:NX] = grid[:,Y,Z]
		elif dir == 1:
			X = np.random.randint(NX)
			Y = np.random.randint(NY-ch)
			Z = np.random.randint(NZ)
			chords[i,0:NY] = grid[X,:,Z]
		else:
			X = np.random.randint(NX)
			Y = np.random.randint(NY)
			Z = np.random.randint(NZ-ch)
			chords[i,0:NZ] = grid[X,Y,:]
	return chords
def get_pore_distr(chords,spacing):
	N = chords.shape[0]
	radii = []
	for i in range(0,N):
		radius = split_chord(chords[i,:])
		radii  = radii + radius
		sys.stdout.write('       Evaluating chords  %d%%\r' % int(float(i+1)/N*100))
		sys.stdout.flush()
	edges = np.arange(-0.5,chords.shape[1]+0.5,1)
	p, edges = np.histogram(radii,edges,normed=True)
	centers = np.arange(0,chords.shape[1])*spacing
	# trim zeros at the end
	p = np.trim_zeros(p,'b')
	centers = centers[0:p.shape[0]]
	return p,centers
def pd2pc(c, *pd):
	N = len(pd)
	pc = np.zeros(N)
	for i in range(1,N): # skip 0
		# pc[0:i+1] = pc[0:i+1] + pd[i]*c[0:i+1]/(c[i]*c[i])# p(c) = p(di)*c/c_i^2
		pc[0:i+1] = pc[0:i+1] + pd[i]*c[0:i+1]/(c[i]*c[i])# p(c) = p(di)*c/c_i^2
		# add penalty
		# if pd[i] < 0:
		# 	pc[0:i+1] = pc[0:i+1]+pd[i]*pd[i]*10000 # square the negative pd(i) and increase enormously the p(c)
	return pc
def optimize_pd(par0, xdata, ydata):
	N = par0.shape[0]
	# par, parcov = optimization.curve_fit(pd2pc, xdata, ydata, p0=par0)
	err = lambda p: np.mean(((pd2pc(xdata,*p)-ydata)/1)**2)
	bounds = []
	for i in range(0,N):
		bounds.append((-1,None))
	par = minimize(err, par0, bounds=bounds, method="L-BFGS-B").x
	return np.asarray(par)
def transform(c,pc_copy):
	pc = np.copy(pc_copy)
	# given is range of chords, c, and its probability, pd = p(c)
	N = c.shape[0]
	pd = np.zeros(N) # probability of d
	if False: # direct calculation of p(d)
		for i in range(0,N):
			c_curr, pc_curr = c[N-1-i], pc[N-1-i]
			if c_curr > 0.0:
				pd_curr   = pc_curr*c_curr # since p(c) ~ p(d_i)*c / d_i^2 and hence p(c=d_i) ~ p(d_i)/d_i
				delta_pc  = c[0:N-i]*pd_curr/(c_curr*c_curr)
				pc[0:N-i] = pc[0:N-i] - delta_pc # reduced p(c), original
				pd[N-1-i] = pd_curr
		pd = smooth(pd,10)
		# pd = shift2avoid0s(pd)
	if True: # least squares
		pd = optimize_pd(np.copy(pc_copy), c, pc)
		# c, pd = rebin(c,pd,2)
	pd = pd/pd.sum() # normalize
	return c, pd

def savefig(V,name):
    # figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    imgplot = plt.imshow(V) # ,interpolation='none')
    plt.colorbar()
    plt.savefig('%s.png' % name)
    plt.close()
def guessU(U,Nx):
    for i in range(0,Nx):
        U[i,:,:] = 1-1.0*i/Nx
    U[ 0,:,:] = 1
    U[-1,:,:] = 0
    return U
def equilibrateU(U,grid,Nx,Ny,Nz, plane,odir,output,spacing):
	x,y,z = np.arange(0,Nx), np.arange(0,Ny), np.arange(0,Nz)
	ftort = open(odir+'/'+output+'_tortuosity-%splane.dat' % plane,'w')
	xf,yf,zf = np.mod(x+1,Nx), np.mod(y+1,Ny), np.mod(z+1,Nz) # array shifted forward
	xb,yb,zb = np.mod(x-1,Nx), np.mod(y-1,Ny), np.mod(z-1,Nz) # array shifted backward
	neighbors = grid[xf,:,:]+grid[xb,:,:]+grid[:,yf,:]+grid[:,yb,:]+grid[:,:,zf]+grid[:,:,zb]  # number of non-zero neighbors
	# make zero-neighors disappear
	arr = np.where(neighbors == 0)
	neighbors[arr] = 9999
	tortvec = []
	# loop over equilibrating U and stop when reached specified accuracy
	old_t = 0.0
	U[0,:,:] = 1 # just in case
	t_prev, t_prev_prev, fluct = 0.,0.,1.
	for i in range(0,10*Nx*Ny*Nz): # Nx*Ny*Nz is the limit
		Ug = U*grid # by this multiplication ignore the empty sites, i.e. where grid==0
		U[1:-1,:,:] = ( Ug[0:-2,:,:] + Ug[2:Nx,:,:] + Ug[1:-1,yb,:] + Ug[1:-1,yf,:] + Ug[1:-1,:,zb] + Ug[1:-1,:,zf] )/neighbors[1:-1,:,:]
		I = Ug[-2,:,:].sum()
		Rideal = float(Nx-1)/(Ny*Nz)
		Reff = 1/I # R = U/I
		tortuosity = Reff/Rideal
		sys.stdout.write('  (%d%%) ~conductivity^-1: %.4e (fluct %.2f%%)  \r ' % (int(10*i/(Nx*Ny*Nz)),tortuosity,100*fluct))
		sys.stdout.flush()
		if i%500 == 0: # save progress every N steps
			fluct = np.std(np.array([tortuosity, t_prev, t_prev_prev]))/tortuosity
			if fluct < 0.01: #0.0001: # less than 0.01%b
				break
			t_prev_prev = t_prev
			t_prev = tortuosity
			savefig(U[1:,:,:].sum(axis=2)/Nz,odir+'/'+output+'_U-%splane' % (plane))
			ftort.write('%d    %.4f\n' % (i,tortuosity))
			tortvec.append(tortuosity)
			plt.plot(np.linspace(0,len(tortvec)-1,len(tortvec)),np.asarray(tortvec),linewidth=4.0,alpha=0.5)
			plt.savefig(odir+'/'+output+'_tortuosity-%splane.png' % plane)
			plt.close()
			Ug = U*grid
			plt.plot(np.arange(Nx)*spacing,Ug.sum(axis=2).sum(axis=1)/grid.sum(axis=2).sum(axis=1),'g-',linewidth=4.0,alpha=0.5)
			plt.savefig(odir+'/'+output+'_Uprofile-%splane.png' % plane)
			plt.close()
	ftort.close()
	return U
def drawI(U,grid, Nx,Ny,Nz, epsilon):
    y, z  = np.arange(0,Ny), np.arange(0,Nz)
    yf,zf = np.mod(y+1,Ny),  np.mod(z+1,Nz) # array shifted forward
    yb,zb = np.mod(y-1,Ny),  np.mod(z-1,Nz) # array shifted backward
    # Ug = U*grid

    Imap,Iim,Iip,Ijm,Ijp,Ikm,Ikp = np.zeros([Nx,Ny,Nz]),np.zeros([Nx,Ny,Nz]),np.zeros([Nx,Ny,Nz]),np.zeros([Nx,Ny,Nz]),np.zeros([Nx,Ny,Nz]),np.zeros([Nx,Ny,Nz]),np.zeros([Nx,Ny,Nz])
    Iim[1:Nx,:,:] = U[1:Nx,:,:] - U[0:-1,:,:]
    Iip[0:-1,:,:] = U[0:-1,:,:] - U[1:Nx,:,:]
    Ijm[:,:,:]    = U[:,:,:] - U[:,yb,:]
    Ijp[:,:,:]    = U[:,:,:] - U[:,yf,:]
    Ikm[:,:,:]    = U[:,:,:] - U[:,:,zb]
    Ikp[:,:,:]    = U[:,:,:] - U[:,:,zf]

    Imap = np.abs(Iim)+np.abs(Iip)+np.abs(Ijm)+np.abs(Ijp)+np.abs(Ikm)+np.abs(Ikp)
    Imap[1:-1,:,:] = Imap[1:-1,:,:]/2
    Imap = Imap*grid
    Imap[np.where(Imap < epsilon)] = 0
    eff_porosity = float(np.where(Imap > 0)[0].shape[0])/(Nx*Ny*Nz)
    print '    Effective porosity: %.2f%%                         ' % (100*eff_porosity)
    return eff_porosity
