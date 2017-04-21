import MDAnalysis
import numpy.linalg as nl
import numpy
from MDAnalysis.analysis.leaflet import LeafletFinder


#calculation of lipid tilt angle vs bilayer normal with MDAnalysis (https://pythonhosted.org/MDAnalysis/index.html) by Yoav Atsmon Raz. Any questions are welcomed to at my github repo.  

#read the trajectory
u = MDAnalysis.Universe('lbpa_sn3_md.1.pdb', 'mdt.trr')
out_hg = open('tilt_angle_hg.xvg', 'w')

#obtain center of mass of the computed bilayer to allow seperate calculation for the upper and lower leaflets. In this example I got the value 55 which is used in the atom selection lines accordingly. You should run the script the first time  just for this value. Alternatively you can get it by going into python prompt and inputting these lines after you imported MDAnalysis.

L = LeafletFinder('lbpa_sn3_md.1.pdb', 'name P65')
leaflet0 = L.groups(0)
leaflet1 = L.groups(1)
atoms0=leaflet0.residues.atoms
atoms1=leaflet1.residues.atoms
print >> out_hg, '# angle[deg]'

cen_atoms=u.select_atoms('resname LBPA')
cen=cen_atoms.center_of_mass()
#atom selection
sel_a1a = u.select_atoms('prop z > 55 and resname LBPA and (name P65)')
sel_a2a = u.select_atoms('prop z > 55 and resname LBPA and (name O67)')
sel_a3a = u.select_atoms('prop z > 55 and resname LBPA and (name O66)')

sel_a1b = u.select_atoms('prop z < 55 and resname LBPA and (name P65)')
sel_a2b = u.select_atoms('prop z < 55 and resname LBPA and (name O67)')
sel_a3b = u.select_atoms('prop z < 55 and resname LBPA and (name O66)')

#input normals to calculate the dot products from
z1_ax = numpy.array([0, 0, 1])
z2_ax = numpy.array([0, 0, -1])


#calculate for each frame in the trajectory
for ts in u.trajectory:
#upper leaflet calculation
    pos_a1a = sel_a1a.positions
    pos_a2a = sel_a2a.positions
    pos_a3a = sel_a3a.positions
    r1a = sel_a2a.positions - sel_a1a.positions
    r2a = sel_a3a.positions - sel_a1a.positions
    r_comba = r1a + r2a
    r_combaz = numpy.inner(r_comba,z1_ax)
    factor_ofa = numpy.power(r_comba,2)
    factor_of2a = numpy.sqrt(numpy.sum(factor_ofa,axis=1))
    cos_ofa = numpy.divide(r_combaz,factor_of2a)
    a = numpy.arccos(cos_ofa)
    fin = numpy.rad2deg(a) 
    fin2 = numpy.mean(fin)
    fin2_std = numpy.std(fin)
#lower leaflet calculation
    pos_a1b = sel_a1b.positions
    pos_a2b = sel_a2b.positions
    pos_a3b = sel_a3b.positions
    r1b = sel_a2b.positions - sel_a1b.positions
    r2b = sel_a3b.positions - sel_a1b.positions
    r_combb = r1b + r2b
    r_combbz =  numpy.inner(r_combb,z2_ax)
    factor_ofb = numpy.power(r_combb,2)
    factor_of2b = numpy.sqrt(numpy.sum(factor_ofb,axis=1))
    cos_ofb = numpy.divide(r_combbz,factor_of2b)
    b = numpy.arccos(cos_ofb)
    finb = numpy.rad2deg(b)
    fin2b = numpy.mean(finb)
    fin2b_std = numpy.std(finb)
    fin_avg = (fin2 + fin2b)/2
#output of results
    print >> out_hg, ts.frame,
    print >> out_hg, '%f %f %f %f' % (fin2, fin2_std ,fin2b, fin2b_std),
    print >> out_hg, '\n',
