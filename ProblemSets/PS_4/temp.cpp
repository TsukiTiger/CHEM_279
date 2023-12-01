int CNDO::updateG(){
    int num_atoms = mol.getnum_atoms ();
    arma::vec P_t = arma::zeros (num_atoms);
    size_t k_AO = 0;
    for (size_t k=0; k < num_atoms; k++)
    {
        for (auto ao : mol.mAtoms[k].mAOs)
        {
            P_t(k) += Pa(k_AO, k_AO) + Pb(k_AO, k_AO);
            k_AO++;
        }
    }
    P_t.print("P_t");
    k_AO = 0;
    for (size_t k = 0; k < num_atoms; k++)
    {
        double gammaAA = gamma(k, k);
        for (auto ao_A : mol.mAtoms[k].mAOs)
        {
            Ga(k_AO, k_AO) = (P_t(k) - Pa(k_AO, k_AO));
            Gb(k_AO, k_AO) = (P_t(k) - Pb(k_AO, k_AO));
            size_t j_AO = 0;
            for (size_t j = 0; j<num_atoms; j++)
            {
                double gammaAB = gamma(k, j);
                if (k!= j)
                {
                    Ga(k_AO, k_AO) +=P_t(j)*gammaAB;
                    Gb(k_AO, k_AO) +=P_t(j)*gammaAB;
                }
                for (auto ao : mol.mAtoms[j].mAOs)
                {
                    if (k_AO != j_AO)
                    {
                        Ga(k_AO, j_AO) = -gammaAB *Pa(k_AO, j_AO);
                        Gb(k_AO, j_AO) = -gammaAB *Pb(k_AO, j_AO);
                    }
                    j_AO++;
                }
            }
            k_AO++;
        }
    }
    Ga.print("Ga");
    Gb.print("Gb");
    return 0;
}


int CNDO: : run ()
arma: imat Fa = H_core + Ga;
arma: :mat Fb = H_core + Gb;
arma::mat Pa_old, Pb_old;
size_t k = 0;
for
(; k < max_iter; k++)
cout « "Iteration:"<< k << endl;
// Pa. print ("Pa");
/1 Ga.print("Ga");
Fa-print ("Fa") ;
Fb. print ("Fb");
Pa_old = Pa;
Pb_old = Pb;
arma: :eig_sym (Ea, Ca, Fa);
arma::eig_sym (Eb, Cb, Fb):
cout « "after solving eigen equation: " « k « endl;
Ca. print ("Ca");
Cb.print ("Cb" );
cout " «p « " < q << endl;
Pa = Ca. cols (0, p - 1) *
     Ca. cols (0, p - 1). t/);
if
(q > 0)
Pb = Cb. cols (0, q - 1) * Cb. cols (0, q - 1). t();
else
Pb. zeros ();
if (arma: :approx_equal (Pa, Pa_old,
"absdiff", tol) &&
arma: :approx_equal(Pb, Pb_old, "absdiff", tol))
break;
Pa. print ("Pa_new");
Pb. print ("Pb_new");
updateGi
        Fa = H_core + Ga;
Fb = H_core + Gb;
if
(k == max_iter)
{
cout << "Error: the job could not be finished in " « max_iter « "iterations. In";
return 1;
Ea. print ("Ea");
Eb. print ("Eb");
Ca. print ("Ca");
Cb. print ("Cb");
return 0;}

double CNDO: :getEnergy ()
{
        arma::mat Ptotal = Pa + Pb;

Ee = arma:: dot (Pa, Ga) / 2. + arma:: dot (Pb, Gb) / 2.;

Ee += arma: :dot (Ptotal, H_core);
Etotal = Ee + Ec;


cout « "Nuclear Repulsion Energy

cout « "Electron Energy is " « Ee « " ev." « endl;
return Etotal;}

CNDO:init()
{
    int
    num_atoms = mol-getnum_atoms);
//Do H_core part
    size_t k_A0 = 0;
    for (size_t k =
            0; k < num_atoms; k++)
    {
        Atom
                &A_AO = mol. mAtoms [k];
        double ZA =
                double (A_AO. VAN);
        double
                gammaAA = gamma (k, k);
        CNDO_para A_para = CNDO_para_map LA_AO. namel;
        for (auto ao_A : A_AO. mAOs)
        {
            H_core(k_A0, k_AO) = -A_para. IA[ao_A.get_lable()] - (ZA - 0.5) * gammaAA;
// H_core. print ("H_core");
            size_t j_A0 = 0;
            for (size_t j = 0; j < num_atoms; j++)
            {
                Atom &B_AO = mol. mAtoms Ij;
                CNDO_para B_para = CNDO_para_map [B_AO. namel ;
                double averagebeta = -(A_para. beta + B_para. beta) / 2.;
                if (k!= j)
                    H_core (k_AO, k_A0) == double(B_AO.VAN) * gamma (k, j);
// cout<< k_A0 < " " « H_core(k_A0,
                /I cout<< k_A0 < "" << H_core(k_A0, k_A0)« endl;
                for (auto aos: B_AO. mAOs)
                if (kAO != j_A0)
                    A core KAO, JAO) = averagebeta * SIKAO, JA);
                j_A0++;
            }
            k_A0++;
            if
                    (k_AO != dim)
            {
                cout « "warn! the number of AOs is wrong." « endl; return 1;
            }
            H_core-print ("H_core");

            return 0;}