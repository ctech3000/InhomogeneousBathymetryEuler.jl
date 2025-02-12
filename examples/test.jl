#= b(x) = eval_bath(bath,x)
b_prime(x) = eval_bath(bath,x,1)
b_pprime(x) = eval_bath(bath,x,2)
u_0(x,z) = cos(x-domain.x_R)*cosh(z-b(x))
u_0_dx(x,z) = -sin(x-domain.x_R)*cosh(z-b(x)) - cos(x-domain.x_R)*sinh(z-b(x))*b_prime(x)
u_0_dz(x,z) = cos(x-domain.x_R)*sinh(z-b(x))
f_0(x,z) = -cos(x-domain.x_R)*cosh(z-b(x)) + sin(x-domain.x_R)*sinh(z-b(x))*b_prime(x) + sin(x-domain.x_R)*sinh(z-b(x))*b_prime(x) - cos(x-domain.x_R)*sinh(z-b(x))*b_pprime(x) + cos(x-domain.x_R)*cosh(z-b(x))*b_prime(x)^2 + cos(x-domain.x_R)*sinh(z-b(x))
 =#