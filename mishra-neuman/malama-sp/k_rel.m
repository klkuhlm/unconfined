function u = k_rel(z)

ak = 0.5;
b2 = 20.0;
phi_k =0.10;
phi_a = 0.020;
u = exp(-ak*(phi_k - phi_a + z - b2));