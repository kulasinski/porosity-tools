from libs_porosity import *

if __name__ == "__main__":
    input_structure, output, odir, spacing, radius, ifsavegrid, ifperiodic, ifdistribution, tdir \
                                                        = process_input(sys.argv)
    if not os.path.exists(odir+'/'):
        os.system('mkdir %s' % odir)
    flog = open(odir+'/'+output+'.log','w')
    flog.write('Input structure: %s\n' % input_structure)
    if input_structure.endswith('.gro'):
        # read the structure coordinates
        atomtype, coords, X, Y, Z = read_structure(input_structure)
        flog.write('Number of atoms: %d\n' % len(atomtype))
        # determine the grid elements
        NX = int(X/spacing)
        NY = int(Y/spacing)
        NZ = int(Z/spacing)
        flog.write('Radius of the probe molecule: %f nm\n' % radius)
        flog.write('Grid spacing: %f nm\n' % spacing)
        flog.write('Grid dimensions: %d x %d x %d\n' % (NX,NY,NZ))
        flog.write('Number of grid elements: %.1g (%d)\n' % (NX*NY*NZ,NX*NY*NZ))
        # initialize grid, 1 means the space is free
        G = np.ones([NX,NY,NZ])
        # remove the grid points where atoms are
        if ifperiodic:
            G = remove_grid_points_pbc(G, atomtype, coords, X,Y,Z, NX,NY,NZ, spacing, radius)
        else:
            G = remove_grid_points(G, atomtype, coords, X,Y,Z, NX,NY,NZ, spacing, radius)
    elif input_structure.endswith('.npy'):
        # load grid from file
        G = np.load(input_structure)
        NX,NY,NZ = G.shape
    print '   Grid dimensions: %d x %d x %d' % (NX,NY,NZ)
    print '   Number of grid elements: %.1g (%d)' % (NX*NY*NZ,NX*NY*NZ)
    Gsum = G.sum()
    vol = Gsum*spacing*spacing*spacing
    print '   Remaining number of gridpoints: %.1g (%d)' % (Gsum,Gsum)
    print '   Porosity:',round(float(Gsum)/(NX*NY*NZ)*100,1),'%'
    print '   Pore volume:',vol,'nm^3'
    flog.write('Remaining number of grid elements: %.1g (%d)\n' % (Gsum,Gsum))
    flog.write('Porosity: %.1f%%\n' % round(float(Gsum)/(NX*NY*NZ)*100,1))
    flog.write('Pore volume: %f nm^3\n' % float(vol))
    if ifsavegrid=='npy':
        np.save(odir+'/'+output+'_grid',G)
    elif ifsavegrid=='gro':
        save_grid2gro(odir+'/'+output+'_grid.gro',G, NX,NY,NZ, spacing)
    totalsurf=get_surf(G,NX,NY,NZ,spacing)
    print '   Total surface area: %.2f nm^2' % totalsurf
    print '   Surface to volume ratio: %.2f nm^-1' % float(totalsurf/vol)
    flog.write('Total surface area: %.2f nm^2\n' % totalsurf)
    flog.write('Surface to volume ratio: %.2f nm^-1\n' % float(totalsurf/vol))
    flog.close()
    print ''

    if ifdistribution:
        print '   Calculating pore size distribution'
        ch_max = (min(NX,NY,NZ)-1)*spacing
        print '     Obtaining chords'
        chords = get_chord_vec_full(NX*NY*NZ,G,spacing,ch_max)
        print '     Calculating p(c)'
        p,ch = get_pore_distr(chords,spacing)
        print '     Transformation p(c)->Q(D)'
        d,Q = transform(ch,p)
        fname = odir+'/'+output+'_pore_distribution.csv'
        print '     Saving to',fname
        np.savetxt(fname, np.vstack((ch,p/spacing,d,Q/spacing)).T,delimiter=',', fmt='%5e',\
            header='chord length c (nm), probability p(c), pore diameter D (nm), probability Q(D)')
        plt.plot(ch,p/spacing,'b--',label='p(c)')
        plt.plot(d,Q/spacing,'k-',label='Q(D)')
        plt.xlabel('ch | D (nm)')
        plt.ylabel('p(c) | Q(D)')
        plt.legend()
        plt.title('Probability distribution function')
        plt.savefig(odir+'/'+output+'_pore_distribution.png',dpi=150)
        print ''

    if tdir in ('x','y','z'):
        if tdir == 'x':
            plane = 'yz'
        elif tdir == 'y':
            plane = 'zx'
        elif tdir == 'z':
            plane = 'xy'
        print '   Calculating tortuosity in %s direction (%s plane)' % (tdir,plane)
        if plane == 'xy':
            pass
        elif plane == 'yz':
            grid = np.swapaxes(G,0,2) # (x,y,z)->(z,y,x)
            grid = np.swapaxes(G,0,1) # (z,y,x)->(y,z,x)
        elif plane == 'zx':
            grid = np.swapaxes(G,0,2) # (x,y,z)->(z,y,x)
            grid = np.swapaxes(G,1,2) # (z,y,x)->(z,x,y)
        NX,NY,NZ = grid.shape
        epsilon = 1.0/NX/1000 # negligible current, used to estimate effective porosity. 1000x smaller than average U drop per element
        U = np.zeros([NX,NY,NZ])
        savefig(grid[1:,:,:].sum(axis=2)/NZ,odir+'/'+output+'_ave_porosity-%splane' % plane)

        U = guessU(U,NX)
        U = equilibrateU(U,grid,NX,NY,NZ, plane,odir,output,spacing)
        np.save(odir+'/'+output+'_U-matrix',U)
        p = drawI(U, grid, NX,NY,NZ, epsilon)

        # calculate currents
        I = (U[-2,:,:]*grid[-2,:,:]).sum()
        Rideal = float(NX-1)/(NY*NZ)
        Reff = 1/I # R = U/I
        tortuosity = Reff/Rideal
        savefig(U.sum(axis=2)/NZ,odir+'/'+output+'_U-%splane' % plane)
        if False:
            print 'Grid:\n',grid
            print 'U:\n',np.array_str(U*grid,precision=2) # zeros represent also the empty sites
            print '  occupancy:',float(grid.sum())/(NX*NY*NZ)
            print '  Isum: %.4f' % I
            print '  Rideal: %.4f' % Rideal
            print '  Reff: %.4f' % Reff
        t_real = p*tortuosity
        print '     Conductivity^-1: %.4f   ' % tortuosity
        print '     Tortuosity: %.4f     ' % t_real
        print ''

    print "   ==========================================="
    print "   |*| Thank you for using POROSITY TOOLS |*|"
    print "   |*| If have any concerns or questions, |*|"
    print "   |*| e-mail me: kulasinski AT gmail.com |*|"
    print "   ==========================================="
    print "                  Please cite:                  "
    print "   Kulasinski K. and Guyer R., Quantification of"
    print "   nanopore networks: application to amorphous  "
    print "   polymers, J. Comput. Chem., in review        "
    print "   ===========================================  "
