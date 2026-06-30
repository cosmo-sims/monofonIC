[setup]
GridRes = 8
BoxLength = 100
zstart = 49
LPTorder = 1
DoBaryons = no
ParticleLoad = sc

[cosmology]
Omega_m = 0.3099
Omega_b = 0.0488911
Omega_L = 0.6900210718
H0 = 67.742
Tcmb = 2.7255
n_s = 0.96822
A_s = 2.1064e-9
k_p = 0.05
transfer = eisenstein
ztarget = 2.5
BackscalingRadiation = @BACKSCALING_RADIATION@
BackscalingNeutrinos = @BACKSCALING_NEUTRINOS@
m_nu1 = @M_NU1@
m_nu2 = @M_NU2@
m_nu3 = @M_NU3@

[random]
generator = NGENIC
seed = 12345

[execution]
NumThreads = 1
WriteBackscalingDiagnostics = true

[output]
format = generic
filename = @TEST_NAME@.hdf5
