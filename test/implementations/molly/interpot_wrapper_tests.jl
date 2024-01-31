
### set up small nonorthogonal TiAl system 
atom1 = AtomsBase.Atom(:Ti, (@SVector [1.47692848, 1.36039173,1.98890986]) * u"Å")
atom2 = AtomsBase.Atom(:Al, (@SVector [-0.07703613, 0.04783871, -0.03196409]) * u"Å")
atoms = [atom1, atom2]
box = [[2.8186357315284174, 0.0, 0.0], 
       [0.0, 2.8186357315284174, 0.0], 
       [0.0, 0.1959273402867757, 4.064092262967999]]u"Å"
bcs = [Periodic(), Periodic(), Periodic()]
small_nonorth_sys = FlexibleSystem(atoms, box, bcs);


### set up ACE potential
ace_tial = ACE(species           = [:Ti, :Al],
          body_order        = 3,
          polynomial_degree = 6,
          wL                = 1.5,
          csp               = 1.0,
          r0                = 2.9,
          rcutoff           = 5.5 )

coeffs_tial = [0.025591381519026762, 0.03527125788791906, 0.030612348271342734, 0.06646319273427727, -0.007326117584533097, 0.0021476223879284516, 0.007005564677788114, 0.007769244005502122, 0.0033019143419533176, -0.024266397752728517, -0.03511058549188825, -0.004486106026460792, 0.07649832144216986, -0.017295005584346018, 0.0348519410800185, -0.0026935081045876344, -0.0036996351104090115, 0.04250332320742131, -0.06611126598243479, 0.07452744399669442, -0.08807382022058645, 0.006553101218837121, -0.02825330387435087, -0.005070437887508557, 0.017488241826946662, -0.041461388491636234, -0.050152804966179194, 0.014554551186620662, 0.005494466857846328, 0.03395869840669037, -0.12004390275966798, 0.07758118243125994, 0.024624168020804672, 0.0006581992277555695, -0.002196641935532242, 0.03231551745953874, 0.0005431753297032715, 0.009602374511533056, 0.028907266845791348, 0.03557855646347803, 0.000832998634326787, 0.019238505326450918, -0.007863457928993406, 0.03497657242548427, -0.058485491203844206, -0.025527625067137013, -0.003851837725125408, 0.019472633328804008, -0.04975455754968226, 0.008243807089528446, 0.020612783411412677, -0.07411984524326856, 0.007005564677788615, 0.0077692440055020795, 0.003301914341951156, -0.024266397752728874, -0.03511058549188732, -0.004486106026461139, 0.004050167320520888, 0.011275083723136878, -0.009533633282696359, -0.016652089366136488, 0.005947187981081792, -0.0086386178798077, 0.027556838876613178, -0.01794755394550558, 0.0328518497817209, -0.0444008944069233, -0.04142464521909121, 0.014939466653767677, -0.0013061815492572404, -0.008399904141687925, -0.013070571180237286, 0.07022679858374972, 0.03655463426663164, -0.02425878114877371, -0.013089322405632224, 0.007663504514768707, -0.0006932536563853398, -0.015392489165582057, -0.005333834581033578, 0.000966860983206308, -0.06259246382383571, 0.04896321372972445, -0.012976299346956766, 0.00575471543255263, -0.010710826925328487, 0.009130893987440367, 0.025455356200891677, 0.03737467186835743, -0.061410072131176816, 0.05891873535070835, 0.02179281899408886, 0.02640823532251172, 0.0002904232787473565, -0.028649695323579742, -0.018788426151163752, -0.004911376520223526, 0.034242726688142995, -0.008960425717451335, -0.04627434272845332, 0.05321984323617818, 0.013802856612787245, 0.023560957961937884]

lb_tial = LBasisPotential(coeffs_tial,[0.0],ace_tial)

ip_forces = InteratomicPotentials.force(small_nonorth_sys, lb_tial);
ip_pe = InteratomicPotentials.potential_energy(small_nonorth_sys, lb_tial);

# set up wrapper
inter_ace_tial = InteratomicPotentialInter(lb_tial,u"eV", u"Å")
inter_ace_unitless = InteratomicPotentialInter(lb_tial,NoUnits, NoUnits)

# TODO: Add tests to ensure proper errors are raised if wrapper not constructed properly

# Check that the wrapper struct is set up properly
@test typeof(inter_ace_tial.potential) <: AbstractPotential
@test dimension(inter_ace_tial.energy_units) == dimension(u"J")
@test dimension(inter_ace_tial.length_units) == dimension(u"m")
@test dimension(inter_ace_unitless.energy_units) == dimension(inter_ace_unitless.length_units) == NoDims

# Check direct Molly.forces and Molly.potential_energy call
molly_naive_forces = Molly.forces(inter_ace_tial, small_nonorth_sys)
molly_naive_pe     = Molly.potential_energy(inter_ace_tial,small_nonorth_sys)
@test ip_forces*u"eV/Å" == molly_naive_forces
@test ip_pe*u"eV" == molly_naive_pe

molly_naive_f_unitless  = Molly.forces(inter_ace_unitless, small_nonorth_sys)
molly_naive_pe_unitless = Molly.potential_energy(inter_ace_unitless,small_nonorth_sys)

@test molly_naive_f_unitless == ip_forces
@test molly_naive_pe_unitless == ip_pe

#set up Molly.System 
general_inters_tial = (inter_ace_tial,)

m_sys_tial = System(Molly.System(small_nonorth_sys,
                                 u"eV",
                                 u"eV/Å");
                    general_inters=general_inters_tial, 
                    loggers=(force=ForceLogger(typeof(1.0u"eV/Å"), 1),
                             energy=PotentialEnergyLogger(typeof(1.0u"eV"),1)
                    ));

simulator = VelocityVerlet(
    dt=0.001u"ps",
    coupling=AndersenThermostat(80u"K", 1.0u"ps"),
)

# test starting configuation force and energy
simulate!(m_sys_tial,simulator,0);
f0_tial = m_sys_tial.loggers.force.history[1]
e0_tial = m_sys_tial.loggers.energy.history[1]

@test all(Vector.(f0_tial) .≈ ip_forces * u"eV/Å") # slight differences in forces
@test e0_tial ≈ ip_pe * u"eV"


# run simulation, check forces and energies 
simulate!(m_sys_tial,simulator,100);

f100_tial = m_sys_tial.loggers.force.history[end]
e100_tial = m_sys_tial.loggers.energy.history[end]

ip_forces_100 = InteratomicPotentials.force(m_sys_tial,lb_tial)
ip_pe_100 = InteratomicPotentials.potential_energy(m_sys_tial,lb_tial)

@test all(Vector.(f100_tial) .≈ ip_forces_100 * u"eV/Å")
@test e100_tial ≈ ip_pe_100 * u"eV"


# TODO: test a unitless molly system as well
# TODO: integration test using AtomsIO and staticAtoms

